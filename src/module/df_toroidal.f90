!=====================================================================!
!                                                                     !
!   module df_toroidal  [Nov 2026]  keyword density_fitting           !
!                                                                     !
!   Density fitting (Coulomb metric) for the toroidal SAO RHF and     !
!   k-MP2. The 4-index integrals are never formed:                    !
!                                                                     !
!     (mu nu|la si) ~ sum_Q B^Q_{mu nu} B^Q_{la si}                   !
!     B^Q_{mu nu}   = sum_P (mu nu|P) [V^-1/2]_{PQ} , V_PQ = (P|Q)    !
!                                                                     !
!   Auxiliary basis: file AuxBasis (format of the file Basis), made   !
!   by tools/make_auxbasis.py with PySCF's rule (JK-fit basis of the  !
!   orbital basis, e.g. def2-svp-jkfit for sto-3g). Toroidal          !
!   Cartesian functions, normalized like the AOs.                     !
!                                                                     !
!   3- and 2-index integrals: the quartet kernel with a dummy s       !
!   function of exponent 0 (the constant 1) next to each auxiliary    !
!   function:  (mu nu|P) = (mu nu|P 1) , (P|Q) = (P 1|Q 1).           !
!                                                                     !
!   Translation symmetry: V is block-circulant over the cells, so     !
!   V^-1/2 is made per k from V(k) = sum_R (P0|Q R) e^{ikR}           !
!   (eigenvalues below 1e-10 dropped, as PySCF's GDF) and brought     !
!   back to real space. Only the rows mu in the unit cell are kept:   !
!     df_B(mu0, nu, Q R)  (real, nfuc x nBas x naux)                  !
!   and  B^{Q R'}_{la R, si} = df_B(la0, si-R, Q (R'-R)).             !
!                                                                     !
!   df_JK  : J and K (rows mu <= nfuc) from df_B, every SCF iteration !
!   df_kMO : (ov|ov) for k-MP2 via L^Q_{pq}(k1,k2)                    !
!            (pk1 qk2|rk3 sk4) = 1/N sum_Q L^Q_pq(k1,k2) L^Q_rs(k3,k4)!
!=====================================================================!

module df_toroidal

      use classification_ERI, only: ERI_function
      use iso_c_binding, only: c_int
      implicit none

      private
      public :: df_init, df_JK, df_kMO, df_ready, df_build_chol
      public :: df_sym                        ! [Nov 2026] density_fitting + symmetry (HF)

      logical, save :: df_ready = .false.
      integer, save :: df_nau = 0, df_naux = 0, df_ncell = 0
      integer, save :: nb_saved = 1                 ! [Nov 2026] BLAS thread count saved by blas_serial
      integer, allocatable, save :: df_coord(:,:)
      double precision, allocatable, save :: df_B(:,:,:)
      ! [Nov 2026] k-space B for the exchange (HF only): k2 over a half set df_kh, conj for the rest
      complex(8), allocatable, save :: df_Bh(:,:,:,:,:)
      integer, save :: df_nkh = 0
      integer, allocatable, save :: df_kh(:), df_khi(:)
      complex(8), allocatable, save :: df_phase(:,:)
      ! [Nov 2026] density_fitting + symmetry (df_init_sym / df_JK_sym)
      logical, save :: df_sym = .false.
      integer, save :: dfs_npiv = 0, dfs_nirr = 0
      complex(8), allocatable, save :: dfs_big(:,:)          ! Bh(a,q,b,k1) for the irreducible k2
      double precision, allocatable, save :: dfs_B0(:,:)     ! B(mu0, nu, Q; k2 = 0), real
      integer, allocatable, save :: dfs_kirr(:), dfs_krep(:), dfs_kg(:), dfs_kcj(:), dfs_kslot(:), dfs_krot(:,:)
      integer, allocatable, save :: dfs_aqp(:,:), dfs_aqc(:,:), dfs_aop(:,:), dfs_aoc(:,:)
      integer(1), allocatable, save :: dfs_aqs(:,:), dfs_aos(:,:)
      double precision, parameter :: metric_tol_aux = 1.0d-10   ! [Nov 2026] pivoted Cholesky of (P|Q), Gaussian aux basis
      double precision, parameter :: metric_tol_cd  = 1.0d-12   ! [Nov 2026] same, Cholesky fitting set
      double precision, parameter :: df_schwarz    = 1.0d-13
      double precision, parameter :: etb_beta      = 2.0d0     ! PySCF ETB_BETA
      ! inputs kept for a rebuild with the MP2-fit basis (file AuxBasis_MP2)
      integer, save :: sv_ndp = 0, sv_iout = 6
      type(ERI_function), allocatable, save :: sv_AO(:)
      double precision, allocatable, save :: sv_geom(:,:)
      integer, allocatable, save :: sv_charge(:)
      integer, allocatable, save :: sv_aoat(:)      ! [Nov 2026] atom of each AO (for the fitting-error report)

contains

      !-----------------------------------------------------------------!
      ! cell arithmetic on the nx*ny*nz grid (idx = iz nx ny + iy nx + ix + 1)
      !-----------------------------------------------------------------!
      integer function cell_of(c)
        use unitcell_module, only: nx, ny, nz
        integer, intent(in) :: c(3)
        integer :: n1, n2, n3
        n1 = max(nx,1) ; n2 = max(ny,1) ; n3 = max(nz,1)
        cell_of = modulo(c(3),n3)*n1*n2 + modulo(c(2),n2)*n1 + modulo(c(1),n1) + 1
      end function cell_of

      integer function csub(r1, r2)            ! cell of R1 - R2
        integer, intent(in) :: r1, r2
        csub = cell_of(df_coord(:,r1) - df_coord(:,r2))
      end function csub

      integer function cadd(r1, r2)            ! cell of R1 + R2
        integer, intent(in) :: r1, r2
        cadd = cell_of(df_coord(:,r1) + df_coord(:,r2))
      end function cadd

      ! half set of k-points: one of each pair (k,-k); df_khi(k) = position of k in the set, 0 if -k is in it
      subroutine df_khalf(ncell)
        integer, intent(in) :: ncell
        integer :: k
        if (allocated(df_kh)) deallocate(df_kh, df_khi)
        allocate(df_kh(ncell), df_khi(ncell))
        df_nkh = 0 ; df_khi = 0
        do k = 1, ncell
          if (csub(1,k) >= k) then
            df_nkh = df_nkh + 1 ; df_kh(df_nkh) = k ; df_khi(k) = df_nkh
          end if
        end do
      end subroutine df_khalf

      !-----------------------------------------------------------------!
      ! read the shells of one atom from the file fname (format of Basis) !
      !-----------------------------------------------------------------!
      subroutine read_shells_file(at, fname)
        use atom_basis
        type(atom), intent(inout) :: at
        character(*), intent(in)  :: fname
        character(len=1000) :: line, key
        character(len=10)   :: tag1, tag2
        integer :: ios, ne, nc, i, j, u
        logical :: inside, found
        double precision, allocatable :: ex(:), co(:,:)

        write(tag1,'(A,I0)') "A ", at%charge
        write(tag2,'(A,I0)') "a ", at%charge
        at%num_s_function = 0 ; at%num_exponent_s = 0
        at%num_p_function = 0 ; at%num_exponent_p = 0
        at%num_d_function = 0 ; at%num_exponent_d = 0
        at%num_f_function = 0 ; at%num_exponent_f = 0
        at%num_g_function = 0 ; at%num_exponent_g = 0

        open(newunit=u, file=fname, status='old', action='read', iostat=ios)
        if (ios /= 0) then
          write(*,'(3a)') ' keyword density_fitting: cannot open the file ', fname, &
                          ' (make it with tools/make_auxbasis.py)'
          error stop
        end if
        inside = .false. ; found = .false.
        do
          read(u,'(A)',iostat=ios) line
          if (ios /= 0) exit
          if (.not. inside) then
            if (trim(adjustl(line)) == trim(tag1) .or. trim(adjustl(line)) == trim(tag2)) then
              inside = .true. ; found = .true.
            end if
            cycle
          end if
          key = adjustl(line)
          if (key(1:1) == 'A' .or. key(1:1) == 'a') exit
          if (key(1:1) /= '$') cycle
          if (index(key, '-TYPE FUNCTIONS') == 0) cycle
          read(u,*) ne , nc
          allocate(ex(ne), co(ne,nc))
          do i = 1 , ne
            read(u,*) ex(i) , (co(i,j), j=1,nc)
          end do
          select case (key(3:3))
          case ('S','s')
            at%num_exponent_s = ne ; at%num_s_function = nc ; at%exponent_s = ex ; at%coefficient_s = co
          case ('P','p')
            at%num_exponent_p = ne ; at%num_p_function = nc ; at%exponent_p = ex ; at%coefficient_p = co
          case ('D','d')
            at%num_exponent_d = ne ; at%num_d_function = nc ; at%exponent_d = ex ; at%coefficient_d = co
          case ('F','f')
            at%num_exponent_f = ne ; at%num_f_function = nc ; at%exponent_f = ex ; at%coefficient_f = co
          case ('G','g')
            at%num_exponent_g = ne ; at%num_g_function = nc ; at%exponent_g = ex ; at%coefficient_g = co
          case default
            write(*,'(2a)') ' keyword density_fitting: shell type not supported in AuxBasis: ', trim(key)
            error stop
          end select
          deallocate(ex, co)
        end do
        close(u)
        if (.not. allocated(at%exponent_s)) allocate(at%exponent_s(0), at%coefficient_s(0,0))
        if (.not. allocated(at%exponent_p)) allocate(at%exponent_p(0), at%coefficient_p(0,0))
        if (.not. allocated(at%exponent_d)) allocate(at%exponent_d(0), at%coefficient_d(0,0))
        if (.not. allocated(at%exponent_f)) allocate(at%exponent_f(0), at%coefficient_f(0,0))
        if (.not. allocated(at%exponent_g)) allocate(at%exponent_g(0), at%coefficient_g(0,0))
        if (.not. found) then
          write(*,'(a,i0,2a)') ' keyword density_fitting: no auxiliary basis for atom "A ', at%charge, &
                               '" in the file ', fname
          error stop
        end if
      end subroutine read_shells_file

      !-----------------------------------------------------------------!
      ! Even-tempered auxiliary set of one atom from its orbital basis,  !
      ! PySCF's rule for a basis without a fitted set                    !
      ! (pyscf.df.addons._aug_etb_element, default version-2.6 branch):  !
      !   per orbital l: emin_l, emax_l over primitives with             !
      !                  max_contraction |c| > 1e-3                      !
      !   l_max  = min(l_max, max_shells(Z)) , max_shells = 1 (H-Be),    !
      !            2 (B-Ca), 3 (Sc-La), 4 (Ce-) ; 0 for ghost atoms      !
      !   pair (li,lj): emax = 2 sqrt(emax_li emax_lj) , emin likewise   !
      !   L = 0..2 l_max: emax_L = max, emin_L = min over li+lj = L      !
      !   n_L = ceil( ln((emax_L+emin_L)/emin_L) / ln beta )             !
      !   exponents emin_L beta^i , i = n_L-1 .. 0 (uncontracted)        !
      !-----------------------------------------------------------------!
      subroutine etb_atom(orb, aux, beta)
        use atom_basis
        type(atom), intent(in)    :: orb
        type(atom), intent(inout) :: aux
        double precision, intent(in) :: beta
        double precision :: emin_l(0:4), emax_l(0:4), eL_min, eL_max, e2max, e2min
        double precision, allocatable :: ex(:), co(:,:)
        integer :: l, lmax, lmax_aux, ms, li, lj, LL, n, i, z

        emin_l = 1.d99 ; emax_l = 0.d0 ; lmax = -1
        call range_l(0, orb%exponent_s, orb%coefficient_s)
        call range_l(1, orb%exponent_p, orb%coefficient_p)
        call range_l(2, orb%exponent_d, orb%coefficient_d)
        call range_l(3, orb%exponent_f, orb%coefficient_f)
        call range_l(4, orb%exponent_g, orb%coefficient_g)
        if (lmax < 0) error stop 'density_fitting: an atom has no orbital basis functions'
        z = orb%charge
        if (z <= 0) then
          ms = 0
        else if (z <= 4) then
          ms = 1
        else if (z <= 20) then
          ms = 2
        else if (z <= 57) then
          ms = 3
        else
          ms = 4
        end if
        lmax = min(lmax, ms) ; lmax_aux = 2*lmax
        if (lmax_aux > 4) error stop 'density_fitting: the generated auxiliary set needs l > 4 (only s..g supported)'

        aux%num_s_function = 0 ; aux%num_exponent_s = 0
        aux%num_p_function = 0 ; aux%num_exponent_p = 0
        aux%num_d_function = 0 ; aux%num_exponent_d = 0
        aux%num_f_function = 0 ; aux%num_exponent_f = 0
        aux%num_g_function = 0 ; aux%num_exponent_g = 0
        allocate(aux%exponent_s(0), aux%coefficient_s(0,0), aux%exponent_p(0), aux%coefficient_p(0,0), &
                 aux%exponent_d(0), aux%coefficient_d(0,0), aux%exponent_f(0), aux%coefficient_f(0,0), &
                 aux%exponent_g(0), aux%coefficient_g(0,0))
        do LL = 0, lmax_aux
          eL_max = -1.d0 ; eL_min = huge(1.d0)
          do li = 0, lmax
            lj = LL - li
            if (lj < 0 .or. lj > lmax) cycle
            e2max = 2.d0*dsqrt(emax_l(li)*emax_l(lj))
            e2min = 2.d0*dsqrt(emin_l(li)*emin_l(lj))
            eL_max = max(eL_max, e2max) ; eL_min = min(eL_min, e2min)
          end do
          n = ceiling( log((eL_max+eL_min)/eL_min) / log(beta) )
          if (n <= 0) cycle
          allocate(ex(n), co(n,n))
          co = 0.d0
          do i = 1, n
            ex(i) = eL_min*beta**(n-i) ; co(i,i) = 1.d0
          end do
          select case (LL)
          case (0) ; aux%num_exponent_s = n ; aux%num_s_function = n ; aux%exponent_s = ex ; aux%coefficient_s = co
          case (1) ; aux%num_exponent_p = n ; aux%num_p_function = n ; aux%exponent_p = ex ; aux%coefficient_p = co
          case (2) ; aux%num_exponent_d = n ; aux%num_d_function = n ; aux%exponent_d = ex ; aux%coefficient_d = co
          case (3) ; aux%num_exponent_f = n ; aux%num_f_function = n ; aux%exponent_f = ex ; aux%coefficient_f = co
          case (4) ; aux%num_exponent_g = n ; aux%num_g_function = n ; aux%exponent_g = ex ; aux%coefficient_g = co
          end select
          deallocate(ex, co)
        end do

      contains

        subroutine range_l(l, e, c)
          integer, intent(in) :: l
          double precision, intent(in) :: e(:), c(:,:)
          integer :: k
          do k = 1, size(e)
            if (size(c,2) == 0) cycle
            if (maxval(abs(c(k,:))) > 1.d-3) then
              emin_l(l) = min(emin_l(l), e(k)) ; emax_l(l) = max(emax_l(l), e(k))
              lmax = max(lmax, l)
            end if
          end do
        end subroutine range_l

      end subroutine etb_atom

      ! the generated set, once per element

      subroutine print_etb(n_atoms, aat, iout)
        use atom_basis
        integer, intent(in) :: n_atoms, iout
        type(atom), intent(in) :: aat(n_atoms)
        integer :: i, j
        logical :: seen
        character(len=1), parameter :: lname(0:4) = (/ 'S','P','D','F','G' /)
        write(iout,'(a)') ' auxiliary basis: even-tempered, generated from the orbital basis'
        write(iout,'(a)') ' (PySCF rule for a basis without a fitted set, beta = 2.0;'
        write(iout,'(a)') '  a file AuxBasis in the run folder is used instead when present)'
        do i = 1, n_atoms
          seen = .false.
          do j = 1, i-1
            if (aat(j)%charge == aat(i)%charge) seen = .true.
          end do
          if (seen) cycle
          write(iout,'(a,i4,a)', advance='no') '   Z =', aat(i)%charge, ' :'
          call one(0, aat(i)%exponent_s) ; call one(1, aat(i)%exponent_p) ; call one(2, aat(i)%exponent_d)
          call one(3, aat(i)%exponent_f) ; call one(4, aat(i)%exponent_g)
          write(iout,'(a)') ''
        end do
      contains
        subroutine one(l, e)
          integer, intent(in) :: l
          double precision, intent(in) :: e(:)
          if (size(e) == 0) return
          write(iout,'(2x,i0,a,a,es10.3,a,es10.3,a)', advance='no') size(e), lname(l), ' (', minval(e), ' -', maxval(e), ')'
        end subroutine one
      end subroutine print_etb

      !-----------------------------------------------------------------!
      ! shells of a function list (same rule as the ERI driver)          !
      !-----------------------------------------------------------------!
      subroutine make_shells(n, F, nsh, first, ncp)
        use clifford_tools, only : orbital_to_lxyz
        integer, intent(in) :: n
        type(ERI_function), intent(in) :: F(n)
        integer, intent(out) :: nsh, first(n), ncp(n)
        integer :: i, j, lt, nc, lxyz(3)
        logical :: same
        nsh = 0 ; i = 1
        do while (i <= n)
          call orbital_to_lxyz(F(i)%orbital, lxyz)
          lt = sum(lxyz) ; nc = (lt+1)*(lt+2)/2
          same = (i + nc - 1 <= n)
          if (same) then
            do j = i+1 , i+nc-1
              call orbital_to_lxyz(F(j)%orbital, lxyz)
              same = same .and. sum(lxyz) == lt .and. F(j)%x == F(i)%x .and. F(j)%y == F(i)%y &
                     .and. F(j)%z == F(i)%z .and. size(F(j)%exponent) == size(F(i)%exponent)
              if (same) same = all(F(j)%exponent == F(i)%exponent)
            end do
          end if
          if (.not. same) nc = 1
          nsh = nsh + 1 ; first(nsh) = i ; ncp(nsh) = nc
          i = i + nc
        end do
      end subroutine make_shells

      !-----------------------------------------------------------------!
      ! df_init: auxiliary basis, (P|Q), V^-1/2, (mu0 nu|P), df_B        !
      ! (after direct_init: AO shells and Schwarz bounds)                !
      !-----------------------------------------------------------------!
      subroutine df_init(ndp, n_atoms, atoms, geometry, nBas, AO, iout, fname_in)

        use atom_basis
        use unitcell_module, only: nfuc, nx, ny, nz, N_cell
        use torus_init     , only: ax, ay, az
        use keywords       , only: c_df_etb, df_cd_tol, s_eig_min   ! [Nov 2026]
        use keywords       , only: c_symmetry, c_MP2, c_ccsdt, c_fast_ccsdt   ! [Nov 2026] density_fitting + symmetry
        use symmetry_clifford, only: sym_ready, sym_nop
        use clifford_tools , only: orbital_to_lxyz, setup_t_tables, bessel_nmax_estimate, ps_on
        use eri_direct     , only: d_nsh, d_nshuc, d_sh_first, d_sh_n, d_Q, d_Qsh
        use constants_module, only: pi
        use omp_lib

        integer           , intent(in) :: ndp, n_atoms, nBas, iout
        type(atom)        , intent(in) :: atoms(n_atoms)
        double precision  , intent(in) :: geometry(n_atoms,3)
        type(ERI_function), intent(in) :: AO(nBas)
        character(*), optional, intent(in) :: fname_in
        character(len=64) :: fname
        logical :: use_file

        type(atom), allocatable :: aat(:)
        type(ERI_function), allocatable :: AUX(:)
        type(ERI_function) :: dumP(1), dumQ(1)
        integer, allocatable :: afirst(:), ancp(:)
        integer :: i, naux, nash, nashuc, ncell, nau, ix, iy, iz, idx, k, r, rp, rq, SP, SQ, SA, SB, it, nt
        integer :: a, b, p, q, lxyz(3), lmax, ncap, ndrop, ndrop_tot, info, lwork, nb
        double precision :: emax, amin, a_dir(3), theta, wmin, w0, t_V, t_3, t_B, invN
        double precision, allocatable :: Vm(:,:), Qa(:), T3(:,:,:), blk(:,:,:,:), w(:), rwork(:)
        double precision, allocatable :: Cph(:,:), Sph(:,:), Ar(:,:), Ai(:,:)
        complex(8), allocatable :: phase(:,:), Mk(:,:), work(:), fk(:,:,:), Zk(:,:), Wk(:,:,:), Wt(:,:,:,:)
        integer :: nu0, m1, kh
        logical :: do_bh
        integer, allocatable :: tq(:,:)
        logical :: cdmode
        double precision :: cd_tol_eff, mtol
        double precision, allocatable :: cd_Vm(:,:), cd_T3(:,:,:)

        fname = 'AuxBasis'
        if (present(fname_in)) then
          fname = fname_in
        else
          sv_ndp = ndp ; sv_iout = iout
          if (allocated(sv_AO)) deallocate(sv_AO, sv_geom, sv_charge)
          allocate(sv_AO(nBas), sv_geom(n_atoms,3), sv_charge(n_atoms))
          sv_AO = AO ; sv_geom = geometry ; sv_charge = atoms(:)%charge
          if (allocated(sv_aoat)) deallocate(sv_aoat)
          allocate(sv_aoat(nBas))
          a = 0
          do i = 1, n_atoms
            do p = 1, atoms(i)%num_s_function + 3*atoms(i)%num_p_function + 6*atoms(i)%num_d_function &
                      + 10*atoms(i)%num_f_function + 15*atoms(i)%num_g_function
              a = a + 1 ; if (a <= nBas) sv_aoat(a) = i
            end do
          end do
          if (a /= nBas) sv_aoat = 0                    ! unknown layout: report without element labels
        end if
        ncell = N_cell ; df_ncell = ncell
        if (allocated(df_coord)) deallocate(df_coord)
        allocate(df_coord(3,ncell))
        do iz = 0, max(nz,1)-1
          do iy = 0, max(ny,1)-1
            do ix = 0, max(nx,1)-1
              idx = iz*max(nx,1)*max(ny,1) + iy*max(nx,1) + ix + 1
              df_coord(:,idx) = (/ ix, iy, iz /)
            end do
          end do
        end do
        allocate(phase(ncell,ncell))
        do k = 1, ncell
          do r = 1, ncell
            theta = 2.0d0*pi*( dble(df_coord(1,k)*df_coord(1,r))/dble(max(nx,1)) &
                             + dble(df_coord(2,k)*df_coord(2,r))/dble(max(ny,1)) &
                             + dble(df_coord(3,k)*df_coord(3,r))/dble(max(nz,1)) )
            phase(k,r) = exp(cmplx(0.d0,theta,kind=8))
          end do
        end do
        invN = 1.d0/dble(ncell)
        if (allocated(df_phase)) deallocate(df_phase)
        allocate(df_phase(ncell,ncell)) ; df_phase = phase

        ! ---------------- auxiliary functions ---------------- !
        ! [Nov 2026] default: Cholesky fitting set (orbital products chosen by a
        ! pivoted Cholesky, tolerance cd_tol_eff). A file AuxBasis in the run folder,
        ! the MP2 file AuxBasis_MP2, or keyword df_etb select a Gaussian aux basis.
        inquire(file=trim(fname), exist=use_file)
        cdmode = .not. (use_file .or. present(fname_in) .or. c_df_etb)
        cd_tol_eff = min(df_cd_tol, 1.d-4*s_eig_min**2)
        ! [Nov 2026] keywords density_fitting + symmetry (HF): symmetric Cholesky fitting set,
        ! B(k1,k2) for the irreducible k2 only (see df_init_sym)
        df_sym = .false.
        if (cdmode .and. .not. present(fname_in) .and. c_symmetry .and. sym_ready) then
          if (sym_nop > 1) then
            if (c_MP2 .or. c_ccsdt .or. c_fast_ccsdt) then
              write(iout,'(a)') ' density_fitting + symmetry: MP2/CC need the full B -> fitting WITHOUT symmetry'
            else
              call df_init_sym(ndp, nBas, AO, ncell, phase, cd_tol_eff, iout)
              deallocate(phase)
              return
            end if
          end if
        end if
        allocate(aat(n_atoms))
        if (cdmode) then
          call cd_build(ndp, nBas, AO, ncell, cd_tol_eff, iout, nau, cd_Vm, cd_T3)
          naux = nau*ncell ; nash = 0 ; nashuc = 0 ; use_file = .true.
          fname = 'Cholesky fitting set (generated)'
          allocate(AUX(0), afirst(0), ancp(0))
        else
        naux = 0
        inquire(file=trim(fname), exist=use_file)
        if (.not. present(fname_in) .and. .not. use_file) fname = 'even-tempered (generated)'
        do i = 1, n_atoms
          aat(i)%charge = atoms(i)%charge
          if (use_file) then
            call read_shells_file(aat(i), trim(fname))
          else
            call etb_atom(atoms(i), aat(i), etb_beta)   ! PySCF aug_etb rule from the orbital basis
          end if
          naux = naux + aat(i)%num_s_function + 3*aat(i)%num_p_function + 6*aat(i)%num_d_function &
                      + 10*aat(i)%num_f_function + 15*aat(i)%num_g_function
        end do
        allocate(AUX(naux))
        call classification(n_atoms, naux, geometry, aat, AUX)
        select case (ndp)
        case (1) ; call normalize_AO_toroidal_1D(naux, AUX, .false.)
        case (2) ; call normalize_AO_toroidal_2D(naux, AUX, .false.)
        case (3) ; call normalize_AO_toroidal_3D(naux, AUX, .false.)
        end select
        allocate(afirst(naux), ancp(naux))
        call make_shells(naux, AUX, nash, afirst, ancp)
        nau = naux/ncell ; nashuc = nash/ncell
        if (nau*ncell /= naux .or. nashuc*ncell /= nash) &
          error stop 'keyword density_fitting: the auxiliary functions must be the same in every cell'
        do SP = 1, nash - nashuc
          if (ancp(SP+nashuc) /= ancp(SP) .or. afirst(SP+nashuc) /= afirst(SP) + nau) &
            error stop 'keyword density_fitting: the auxiliary functions must be the same in every cell'
        end do
        end if
        df_nau = nau ; df_naux = naux

        ! t-tables large enough for the auxiliary exponents and l
        lmax = 0 ; emax = 0.d0
        do i = 1, nBas
          call orbital_to_lxyz(AO(i)%orbital, lxyz)
          lmax = max(lmax, maxval(lxyz(1:ndp))) ; emax = max(emax, 2.d0*maxval(AO(i)%exponent))
        end do
        do i = 1, size(AUX)
          call orbital_to_lxyz(AUX(i)%orbital, lxyz)
          lmax = max(lmax, maxval(lxyz(1:ndp))) ; emax = max(emax, maxval(AUX(i)%exponent))
        end do
        a_dir = (/ ax, ay, az /) ; amin = minval(a_dir(1:ndp))
        ncap = bessel_nmax_estimate(2.d0*emax/(amin*amin)) + 4*lmax + 2
        call setup_t_tables(ncap, ndp)
        ps_on = .false.

        write(iout,'(a)') ''
        write(iout,'(a)') '------------------------------------------------------------'
        write(iout,'(2a)') ' Density fitting (keyword density_fitting), file ', trim(fname)
        write(iout,'(a)') '------------------------------------------------------------'
        if (.not. use_file .and. .not. cdmode) call print_etb(n_atoms, aat, iout)
        if (cdmode) then
          write(iout,'(a)')        ' auxiliary set: orbital products chosen by pivoted Cholesky (Cholesky fitting)'
          write(iout,'(a,es10.2)') '   requested tolerance (1e-8, df_cd_6: 1e-6, df_cd_10: 1e-10) : ', df_cd_tol
          write(iout,'(a,es10.2)') '   smallest overlap eigenvalue lambda_min(S)               : ', s_eig_min
          write(iout,'(a,es10.2)') '   tolerance used = min(requested, 1e-4 lambda_min(S)^2)    : ', cd_tol_eff
        end if
        write(iout,'(a,i10)') ' auxiliary functions per cell     : ', nau
        write(iout,'(a,i10)') ' auxiliary functions (total)      : ', naux
        flush(iout)

        ! ---------------- metric (P0|Q R) ---------------- !
        w0 = omp_get_wtime()
        allocate(Vm(nau,naux))
        Vm = 0.d0
        !$omp parallel do default(shared) private(SP,SQ,blk,dumP,dumQ,p,q) schedule(dynamic,1) collapse(2)
        do SP = 1, nashuc
          do SQ = 1, nash
            call dummy_at(AUX(afirst(SP)), dumP(1))
            call dummy_at(AUX(afirst(SQ)), dumQ(1))
            allocate(blk(ancp(SP),1,ancp(SQ),1))
            call ERI_shell_block_toroidal(ndp, ancp(SP), AUX(afirst(SP):), 1, dumP, &
                                               ancp(SQ), AUX(afirst(SQ):), 1, dumQ, blk)
            do q = 1, ancp(SQ)
              do p = 1, ancp(SP)
                Vm(afirst(SP)+p-1, afirst(SQ)+q-1) = blk(p,1,q,1)
              end do
            end do
            deallocate(blk)
          end do
        end do
        !$omp end parallel do
        if (cdmode) Vm = cd_Vm
        allocate(Qa(naux))
        do r = 1, ncell
          do p = 1, nau
            Qa((r-1)*nau+p) = dsqrt(dabs(Vm(p,p)))
          end do
        end do

        ! V^-1/2 per k, back to real space: mf(P,Q,R) = [V^-1/2](P0, Q R)
        allocate(fk(nau,nau,ncell), Mk(nau,nau), w(nau), rwork(3*nau), work(1))
        ndrop_tot = 0 ; wmin = huge(1.d0)
        mtol = metric_tol_aux ; if (cdmode) mtol = metric_tol_cd
        call zheev('V','U', nau, Mk, nau, w, work, -1, rwork, info)
        lwork = max(1, int(dble(work(1)))) ; deallocate(work) ; allocate(work(lwork))
        do k = 1, ncell
          Mk = (0.d0,0.d0)
          do r = 1, ncell
            Mk = Mk + Vm(:, (r-1)*nau+1:r*nau)*phase(k,r)
          end do
          Mk = 0.5d0*(Mk + conjg(transpose(Mk)))
          ! [Nov 2026] pivoted Cholesky of V(k) (was: eigenvalues, drop < 1e-10):
          ! redundant aux functions are skipped, no division by tiny eigenvalues
          if (csub(1,k) < k) then                      ! V(-k) = conj V(k): keep W(-k) = conj W(k)
            fk(:,:,k) = conjg(fk(:,:,csub(1,k)))
            cycle
          end if
          call chol_f(nau, Mk, mtol, fk(:,:,k), ndrop, wmin)
          ndrop_tot = max(ndrop_tot, ndrop)
        end do
        deallocate(Mk, w, rwork, work, Vm)
        t_V = omp_get_wtime() - w0
        write(iout,'(a,es12.4)') ' metric (P|Q)(k): pivoted Cholesky, tolerance ', mtol
        write(iout,'(a,es12.4)') ' smallest pivot kept               : ', wmin
        write(iout,'(a,i10)')    ' aux functions skipped (max over k): ', ndrop_tot
        flush(iout)

        ! ---------------- (mu0 nu | P) ---------------- !
        w0 = omp_get_wtime()
        allocate(T3(nfuc,nBas,naux))
        T3 = 0.d0
        allocate(tq(2,d_nshuc*d_nsh)) ; nt = 0
        do SA = 1, d_nshuc
          do SB = 1, d_nsh
            nt = nt + 1 ; tq(:,nt) = (/ SA, SB /)
          end do
        end do
        !$omp parallel do default(shared) private(it,SA,SB,SP,blk,dumP,a,b,p) schedule(dynamic,1)
        do it = 1, nt
          SA = tq(1,it) ; SB = tq(2,it)
          do SP = 1, nash
            if (d_Qsh(SA,SB)*maxval(Qa(afirst(SP):afirst(SP)+ancp(SP)-1)) < df_schwarz) cycle
            call dummy_at(AUX(afirst(SP)), dumP(1))
            allocate(blk(d_sh_n(SA),d_sh_n(SB),ancp(SP),1))
            call ERI_shell_block_toroidal(ndp, d_sh_n(SA), AO(d_sh_first(SA):), d_sh_n(SB), AO(d_sh_first(SB):), &
                                               ancp(SP), AUX(afirst(SP):), 1, dumP, blk)
            do p = 1, ancp(SP)
              do b = 1, d_sh_n(SB)
                do a = 1, d_sh_n(SA)
                  T3(d_sh_first(SA)+a-1, d_sh_first(SB)+b-1, afirst(SP)+p-1) = blk(a,b,p,1)
                end do
              end do
            end do
            deallocate(blk)
          end do
        end do
        !$omp end parallel do
        if (cdmode) T3 = cd_T3
        deallocate(tq)
        t_3 = omp_get_wtime() - w0

        ! ---------------- df_B = T3 V^-1/2 , in k-space ---------------- !
        ! [Nov 2026] B(R') = sum_RP T3(RP) mf(R'-RP) is a convolution over cells:
        !   W(k) = [sum_RP T3(RP) e^{ik RP}] f(k) ,  df_B(R') = 1/N Re sum_k W(k) e^{-ik R'}
        ! done for one column nu0 (all its cell images) at a time.
        ! For the SCF exchange it also keeps  Bh(mu,Q,nu; k1,k2) = sum_{R1,R2} df_B(mu, nu R1, Q R2) e^{i k1 R1} e^{i k2 R2}
        ! (k2 over half of the k-points; the other half is the complex conjugate).
        w0 = omp_get_wtime()
        if (allocated(df_B)) deallocate(df_B)
        allocate(df_B(nfuc,nBas,naux))
        do_bh = .not. present(fname_in)
        call df_khalf(ncell)
        if (allocated(df_Bh)) deallocate(df_Bh)
        if (do_bh) allocate(df_Bh(nfuc,nau,nfuc,ncell,df_nkh))
        allocate(Cph(ncell,ncell), Sph(ncell,ncell))
        Cph = dble(transpose(phase)) ; Sph = aimag(transpose(phase))      ! (R, k)
        ! T3(nfuc,nBas,naux) is a (nfuc*nBas*nau, ncell) matrix: aux index = (RP-1)*nau + P.
        ! Only large BLAS calls, outside OpenMP regions (no nested threading).
        m1 = nfuc*nBas*nau ; nb = nfuc*nBas
        allocate(Ar(m1,ncell), Ai(m1,ncell))
        call dgemm('N','N', m1, ncell, ncell, 1.d0, T3, m1, Cph, ncell, 0.d0, Ar, m1)
        call dgemm('N','N', m1, ncell, ncell, 1.d0, T3, m1, Sph, ncell, 0.d0, Ai, m1)
        deallocate(T3)
        allocate(Zk(nb,nau), Wk(nb,nau,1))
        if (do_bh) allocate(Wt(nfuc,nau,nfuc,ncell))
        do k = 1, ncell
          Zk = reshape(cmplx(Ar(:,k), Ai(:,k), kind=8), (/ nb, nau /))
          call zgemm('N','N', nb, nau, nau, (1.d0,0.d0), Zk, nb, fk(1,1,k), nau, (0.d0,0.d0), Wk, nb)
          Ar(:,k) = reshape(dble(Wk), (/ m1 /)) ; Ai(:,k) = reshape(aimag(Wk), (/ m1 /))
          if (do_bh) then
            kh = df_khi(k)
            if (kh > 0) then                                 ! Bh(mu,Q,nu0; k1,k) = sum_R1 W(k)(mu, nu0 R1, Q) e^{i k1 R1}
              !$omp parallel do default(shared) private(r,nu0,q) collapse(2)
              do r = 1, ncell
                do nu0 = 1, nfuc
                  do q = 1, nau
                    Wt(:,q,nu0,r) = Wk((r-1)*nfuc*nfuc + (nu0-1)*nfuc + 1 : (r-1)*nfuc*nfuc + nu0*nfuc, q, 1)
                  end do
                end do
              end do
              !$omp end parallel do
              call zgemm('N','T', nfuc*nau*nfuc, ncell, ncell, (1.d0,0.d0), Wt, nfuc*nau*nfuc, phase, ncell, &
                         (0.d0,0.d0), df_Bh(1,1,1,1,kh), nfuc*nau*nfuc)
            end if
          end if
        end do
        deallocate(Zk, Wk)
        if (do_bh) deallocate(Wt)
        ! df_B(R') = 1/N sum_k [Re W(k) cos(kR') + Im W(k) sin(kR')]
        call dgemm('N','T', m1, ncell, ncell, invN, Ar, m1, Cph, ncell, 0.d0, df_B, m1)
        call dgemm('N','T', m1, ncell, ncell, invN, Ai, m1, Sph, ncell, 1.d0, df_B, m1)
        deallocate(Ar, Ai)
        deallocate(fk, Qa, AUX, afirst, ancp, phase, Cph, Sph)
        t_B = omp_get_wtime() - w0

        write(iout,'(a,f10.2,a)') ' wall time (P|Q) and V^-1/2       : ', t_V, ' s'
        write(iout,'(a,f10.2,a)') ' wall time (mu nu|P)              : ', t_3, ' s'
        write(iout,'(a,f10.2,a)') ' wall time B = (mu nu|P) V^-1/2   : ', t_B, ' s'
        write(iout,'(a,f10.3,a)') ' memory of B (kept)               : ', 8.d0*dble(nfuc)*dble(nBas)*dble(naux)/1024.d0**3, ' GB'
        if (do_bh) write(iout,'(a,f10.3,a)') ' memory of B(k1,k2) for exchange  : ', &
                   16.d0*dble(nfuc)*dble(nfuc)*dble(nau)*dble(ncell)*dble(df_nkh)/1024.d0**3, ' GB (freed after HF)'
        write(iout,'(a,f10.3,a)') ' (stored 4-index integrals would be ', 8.d0*dble(nfuc)*dble(nBas)**3/1024.d0**3, ' GB)'
        call df_fit_report(nBas, AO, iout)
        write(iout,'(a)') '------------------------------------------------------------'
        flush(iout)
        df_ready = .true.
      end subroutine df_init


      !-----------------------------------------------------------------!
      ! [Nov 2026] fitting error of each orbital product, Coulomb metric: !
      !   delta(mu0,nu) = (mu0 nu|mu0 nu) - sum_Q B^Q_{mu0 nu}^2  >= 0    !
      ! = self-repulsion of the part of mu*nu the auxiliary set misses.  !
      ! (mu nu|mu nu) = d_Q(mu,nu)^2 from the Schwarz factors.           !
      ! Grouped by (element, l) of mu and nu, same atom / two atoms.     !
      !-----------------------------------------------------------------!
      subroutine df_fit_report(nBas, AO, iout)
        use unitcell_module, only: nfuc
        use eri_direct     , only: d_Q
        integer           , intent(in) :: nBas, iout
        type(ERI_function), intent(in) :: AO(nBas)
        integer, parameter :: maxg = 400
        double precision, parameter :: rel_floor = 1.d-6     ! products with (mu nu|mu nu) below this: no relative error
        character(len=1), parameter :: lname(0:4) = (/ 'S','P','D','F','G' /)
        integer :: mu, nu, ig, ng, key(5), keys(5,maxg), cnt(maxg), t
        double precision :: dmax(maxg), dsum(maxg), rmax(maxg), ex, del, rel, tot, totex, dworst, rworst
        integer :: iworst(2), irworst(2)
        double precision, allocatable :: d(:,:)

        allocate(d(nfuc,nBas))
        !$omp parallel do default(shared) private(nu,mu) collapse(2)
        do nu = 1, nBas
          do mu = 1, nfuc
            d(mu,nu) = d_Q(mu,nu)**2 - dot_product(df_B(mu,nu,:), df_B(mu,nu,:))
          end do
        end do
        !$omp end parallel do

        ng = 0 ; cnt = 0 ; dmax = 0.d0 ; dsum = 0.d0 ; rmax = 0.d0
        tot = 0.d0 ; totex = 0.d0 ; dworst = -1.d0 ; rworst = -1.d0 ; iworst = 0 ; irworst = 0
        do nu = 1, nBas
          do mu = 1, nfuc
            ex = d_Q(mu,nu)**2 ; del = d(mu,nu)
            tot = tot + del ; totex = totex + ex
            rel = 0.d0
            if (ex > rel_floor) rel = del/ex
            if (del > dworst) then
              dworst = del ; iworst = (/ mu, nu /)
            end if
            if (rel > rworst) then
              rworst = rel ; irworst = (/ mu, nu /)
            end if
            call group_key(mu, nu, key)
            ig = 0
            do t = 1, ng
              if (all(keys(:,t) == key)) then
                ig = t ; exit
              end if
            end do
            if (ig == 0) then
              if (ng == maxg) cycle
              ng = ng + 1 ; ig = ng ; keys(:,ig) = key
            end if
            cnt(ig) = cnt(ig) + 1 ; dsum(ig) = dsum(ig) + del
            dmax(ig) = max(dmax(ig), del) ; rmax(ig) = max(rmax(ig), rel)
          end do
        end do

        write(iout,'(a)') ''
        write(iout,'(a)') ' fitting error of the orbital products (Coulomb metric), mu in cell 0, nu anywhere:'
        write(iout,'(a)') '   delta(mu,nu) = (mu nu|mu nu) - sum_Q (B^Q_mu,nu)^2  [Eh] ,  relative = delta/(mu nu|mu nu)'
        write(iout,'(a)') '   (relative error only for products with (mu nu|mu nu) > 1e-6)'
        write(iout,'(a,es11.3,a,es11.3)') '   sum of delta            : ', tot, '   (sum of (mu nu|mu nu): ', totex
        write(iout,'(a,es11.3,a,a)')       '   largest delta           : ', dworst, '   for ', pair_label(iworst(1), iworst(2))
        write(iout,'(a,es11.3,a,a)')       '   largest relative delta  : ', rworst, '   for ', pair_label(irworst(1), irworst(2))
        write(iout,'(a)') '   by type of product (Z = element, l; centres: one atom, or two atoms / two cells):'
        write(iout,'(a)') '      Z  l    Z  l   centres     pairs    max delta    sum delta   max relative'
        call sort_groups()
        do ig = 1, ng
          write(iout,'(3x,i4,2x,a1,2x,i4,2x,a1,3x,a8,i10,3es13.3)') keys(1,ig), lname(keys(2,ig)), keys(3,ig), &
               lname(keys(4,ig)), merge('one     ', 'two     ', keys(5,ig) == 1), cnt(ig), dmax(ig), dsum(ig), rmax(ig)
        end do
        deallocate(d)

      contains

        subroutine group_key(m, n, kk)
          integer, intent(in)  :: m, n
          integer, intent(out) :: kk(5)
          integer :: z1, z2, l1, l2, s
          call ao_info(m, z1, l1) ; call ao_info(n, z2, l2)
          s = 2
          if (AO(m)%x == AO(n)%x .and. AO(m)%y == AO(n)%y .and. AO(m)%z == AO(n)%z) s = 1
          if (z1 > z2 .or. (z1 == z2 .and. l1 > l2)) then      ! (Z,l) pair unordered
            kk = (/ z2, l2, z1, l1, s /)
          else
            kk = (/ z1, l1, z2, l2, s /)
          end if
        end subroutine group_key

        subroutine ao_info(i, z, l)
          integer, intent(in)  :: i
          integer, intent(out) :: z, l
          z = 0
          if (allocated(sv_aoat)) then
            if (size(sv_aoat) == nBas) then
              if (sv_aoat(i) > 0) z = sv_charge(sv_aoat(i))
            end if
          end if
          l = index('spdfg', AO(i)%orbital(1:1)) - 1
          if (l < 0) l = 0
        end subroutine ao_info

        character(len=40) function pair_label(m, n)
          integer, intent(in) :: m, n
          integer :: z1, z2, l1, l2
          if (m == 0) then
            pair_label = '-' ; return
          end if
          call ao_info(m, z1, l1) ; call ao_info(n, z2, l2)
          write(pair_label,'(a,i0,a,i0,a,a,a,i0,a,i0,a,a,a)') 'AO ', m, ' (Z=', z1, ' ', trim(AO(m)%orbital), &
                ') x AO ', n, ' (Z=', z2, ' ', trim(AO(n)%orbital), ')'
        end function pair_label

        subroutine sort_groups()                           ! largest max delta first
          integer :: i, j, ki(5), ci
          double precision :: a1, a2, a3
          do i = 2, ng
            ki = keys(:,i) ; ci = cnt(i) ; a1 = dmax(i) ; a2 = dsum(i) ; a3 = rmax(i)
            j = i - 1
            do while (j >= 1)
              if (dmax(j) >= a1) exit
              keys(:,j+1) = keys(:,j) ; cnt(j+1) = cnt(j) ; dmax(j+1) = dmax(j) ; dsum(j+1) = dsum(j) ; rmax(j+1) = rmax(j)
              j = j - 1
            end do
            keys(:,j+1) = ki ; cnt(j+1) = ci ; dmax(j+1) = a1 ; dsum(j+1) = a2 ; rmax(j+1) = a3
          end do
        end subroutine sort_groups

      end subroutine df_fit_report

      ! s function of exponent 0 (the constant 1) at the centre of F
      subroutine dummy_at(F, D)
        type(ERI_function), intent(in)  :: F
        type(ERI_function), intent(out) :: D
        D%x = F%x ; D%y = F%y ; D%z = F%z
        D%orbital = 's'
        D%exponent = (/ 0.d0 /)
        D%coefficient = (/ 1.d0 /)
      end subroutine dummy_at

      !-----------------------------------------------------------------!
      ! J(mu0,nu) =        sum_{la si} (mu0 nu|la si) P(la,si)           !
      ! K(mu0,si) = -1/2 sum_{nu la} (mu0 nu|la si) P(nu,la)             !
      !-----------------------------------------------------------------!
      subroutine df_JK(nBas, P, J, K)

        use unitcell_module, only: nfuc

        integer         , intent(in)  :: nBas
        double precision, intent(in)  :: P(nBas,nBas)
        double precision, intent(out) :: J(nBas,nBas), K(nBas,nBas)

        integer :: nau, naux, ncell, nb, r, rs, rq, q, iq, rr, rsg, s
        double precision, allocatable :: Ps(:,:), G(:,:), gam(:), wk(:,:), ew(:), rw(:)
        complex(8), allocatable :: pk(:,:), ck(:,:,:), Kt(:,:,:), L(:,:,:), Lw(:,:,:), Kloc(:,:), zw(:), cc(:,:)
        integer, allocatable :: nok(:)
        integer :: ik, kk, k1, k2, mk1, no, nmax, n, i, info

        if (.not. df_ready) error stop 'df_JK: df_init not called'
        if (df_sym) then                              ! [Nov 2026] density_fitting + symmetry
          call df_JK_sym(nBas, P, J, K)
          return
        end if
        nau = df_nau ; naux = df_naux ; ncell = df_ncell ; nb = nfuc*nBas
        J = 0.d0 ; K = 0.d0

        ! ---- J : gamma(Q R') = sum_{la si} B^{Q R'}_{la si} P(la,si) ---- !
        ! B^{Q R'}_{la0+R , si} = df_B(la0, si-R, Q (R'-R))
        allocate(Ps(nb,ncell), G(naux,ncell), gam(naux))
        do r = 1, ncell
          do rs = 1, ncell                                   ! si' = (rs-1)nfuc + s  ->  si = si' + R
            Ps(1+(rs-1)*nfuc*nfuc : rs*nfuc*nfuc, r) = &
              reshape(P((r-1)*nfuc+1:r*nfuc, (cadd(rs,r)-1)*nfuc+1:cadd(rs,r)*nfuc), (/ nfuc*nfuc /))
          end do
        end do
        call dgemm('T','N', naux, ncell, nb, 1.d0, df_B, nb, Ps, nb, 0.d0, G, naux)
        gam = 0.d0
        do r = 1, ncell
          do rq = 1, ncell                                   ! G(Q R'', R) -> gamma(Q, R''+R)
            gam((cadd(rq,r)-1)*nau+1:cadd(rq,r)*nau) = gam((cadd(rq,r)-1)*nau+1:cadd(rq,r)*nau) &
                                                      + G((rq-1)*nau+1:rq*nau, r)
          end do
        end do
        call dgemv('N', nb, naux, 1.d0, df_B, nb, gam, 1, 0.d0, J, 1)
        ! J was filled as a (nfuc,nBas) block in the first nfuc*nBas elements: move it to rows 1..nfuc
        call spread_rows(nfuc, nBas, J)
        deallocate(Ps, G, gam)

        ! ---- K in k-space [Nov 2026] ---- !
        ! K(mu0, si R) = 1/N sum_K Kt(K) e^{-iKR} ,
        ! Kt(K) = -1/(2N) sum_{k1,Q} M P(k1) M^H ,  M = Bh_Q(k1, K-k1) ,  P(k1) = sum_n w_n c_n c_n^H
        if (.not. allocated(df_Bh)) error stop 'df_JK: k-space B not built'
        allocate(pk(nfuc,nfuc), ck(nfuc,nfuc,ncell), wk(nfuc,ncell), nok(ncell), ew(nfuc), &
                 rw(3*nfuc), zw(4*nfuc), Kt(nfuc,nfuc,ncell))
        do ik = 1, ncell
          pk = (0.d0,0.d0)
          do r = 1, ncell
            pk = pk + P(1:nfuc, (r-1)*nfuc+1:r*nfuc)*df_phase(ik,r)
          end do
          pk = 0.5d0*(pk + conjg(transpose(pk)))
          call zheev('V','U', nfuc, pk, nfuc, ew, zw, 4*nfuc, rw, info)
          if (info /= 0) error stop 'df_JK: zheev failed for P(k)'
          n = 0
          do i = 1, nfuc
            if (dabs(ew(i)) > 1.d-10) then
              n = n + 1 ; ck(:,n,ik) = pk(:,i) ; wk(n,ik) = ew(i)
            end if
          end do
          nok(ik) = n
        end do
        nmax = max(1, maxval(nok))
        nb_saved = blas_serial()                    ! [Nov 2026] no BLAS threads inside OpenMP
        !$omp parallel default(shared) private(kk,k1,k2,mk1,no,L,Lw,Kloc,i,cc)
        allocate(L(nfuc,nau,nmax), Lw(nfuc,nau,nmax), Kloc(nfuc,nfuc), cc(nfuc,nmax))
        !$omp do schedule(dynamic,1)
        do kk = 1, ncell
          Kloc = (0.d0,0.d0)
          do k1 = 1, ncell
            no = nok(k1) ; if (no == 0) cycle
            k2 = csub(kk,k1)
            if (df_khi(k2) > 0) then
              call zgemm('N','N', nfuc*nau, no, nfuc, (1.d0,0.d0), df_Bh(1,1,1,k1,df_khi(k2)), nfuc*nau, &
                         ck(1,1,k1), nfuc, (0.d0,0.d0), L, nfuc*nau)
            else                                             ! Bh(k1,k2) = conj Bh(-k1,-k2)
              mk1 = csub(1,k1)
              cc(:,1:no) = conjg(ck(:,1:no,k1))
              call zgemm('N','N', nfuc*nau, no, nfuc, (1.d0,0.d0), df_Bh(1,1,1,mk1,df_khi(csub(1,k2))), nfuc*nau, &
                         cc, nfuc, (0.d0,0.d0), L, nfuc*nau)
              L(:,:,1:no) = conjg(L(:,:,1:no))
            end if
            do i = 1, no
              Lw(:,:,i) = L(:,:,i)*wk(i,k1)
            end do
            call zgemm('N','C', nfuc, nfuc, nau*no, (1.d0,0.d0), Lw, nfuc, L, nfuc, (1.d0,0.d0), Kloc, nfuc)
          end do
          Kt(:,:,kk) = Kloc*(-0.5d0/dble(ncell))
        end do
        !$omp end do
        deallocate(L, Lw, Kloc, cc)
        !$omp end parallel
        call blas_restore(nb_saved)
        do r = 1, ncell
          pk = (0.d0,0.d0)
          do kk = 1, ncell
            pk = pk + Kt(:,:,kk)*conjg(df_phase(kk,r))
          end do
          K(1:nfuc, (r-1)*nfuc+1:r*nfuc) = dble(pk)/dble(ncell)
        end do
        deallocate(pk, ck, wk, nok, ew, rw, zw, Kt)

      end subroutine df_JK

      ! A holds a (nfuc,nBas) array in its first nfuc*nBas elements -> rows 1..nfuc of A(nBas,nBas)
      subroutine spread_rows(nfuc, nBas, A)
        integer, intent(in) :: nfuc, nBas
        double precision, intent(inout) :: A(nBas*nBas)
        double precision, allocatable :: t(:,:)
        integer :: j
        allocate(t(nfuc,nBas))
        t = reshape(A(1:nfuc*nBas), (/ nfuc, nBas /))
        A = 0.d0
        do j = 1, nBas
          A((j-1)*nBas+1:(j-1)*nBas+nfuc) = t(:,j)
        end do
        deallocate(t)
      end subroutine spread_rows

      !-----------------------------------------------------------------!
      ! (ov|ov) in k-space for k-MP2 (same output as                     !
      ! transform_ERI_AO_to_kMO):                                        !
      !   F^Q_ab(k,kap) = sum_{R,RQ} df_B(a0,bR,Q RQ) e^{ikR} e^{-i kap RQ}
      !   L^Q_pq(k1,k2) = sum_ab C1^*(a,p) C2(b,q) F^Q_ab(k2, k2-k1)     !
      !   ERI(p,q,r,s,k1,k2,k3) = 1/N sum_Q L^Q_pq(k1,k2) L^Q_rs(k3,k4)  !
      !   k4 = k1 - k2 + k3                                              !
      !-----------------------------------------------------------------!
      subroutine df_kMO(nBas, coeff_SAO, p_lo, p_hi, q_lo, q_hi, ERI_out, iout)

        use unitcell_module, only: nfuc, nx, ny, nz, N_cell
        use constants_module, only: pi
        use omp_lib

        integer   , intent(in)  :: nBas, iout
        complex(8), intent(in)  :: coeff_SAO(nBas,nBas)
        integer   , intent(in)  :: p_lo, p_hi, q_lo, q_hi
        complex(8), intent(out) :: ERI_out(p_lo:p_hi,q_lo:q_hi,p_lo:p_hi,q_lo:q_hi,N_cell,N_cell,N_cell)

        integer :: nau, naux, ncell, no, nv, kap, k, k1, k2, k3, k4, r, rq, q, it, off1, off2, iq
        double precision :: theta, w0, t_L, t_E
        complex(8), allocatable :: phase(:,:), L(:,:,:,:,:), Lk(:,:,:,:)
        complex(8), parameter :: one=(1.d0,0.d0), zero=(0.d0,0.d0)

        if (.not. df_ready) error stop 'df_kMO: df_init not called'
        if (allocated(df_Bh)) deallocate(df_Bh)          ! exchange is done: free the k-space B
        call df_rebuild_mp2()
        nau = df_nau ; naux = df_naux ; ncell = df_ncell
        no = p_hi - p_lo + 1 ; nv = q_hi - q_lo + 1
        allocate(phase(ncell,ncell))
        do k = 1, ncell
          do r = 1, ncell
            theta = 2.0d0*pi*( dble(df_coord(1,k)*df_coord(1,r))/dble(max(nx,1)) &
                             + dble(df_coord(2,k)*df_coord(2,r))/dble(max(ny,1)) &
                             + dble(df_coord(3,k)*df_coord(3,r))/dble(max(nz,1)) )
            phase(k,r) = exp(cmplx(0.d0,theta,kind=8))
          end do
        end do

        w0 = omp_get_wtime()
        allocate(L(p_lo:p_hi,q_lo:q_hi,nau,ncell,ncell))     ! L(p,q,Q,k1,k2)
        allocate(Lk(p_lo:p_hi,q_lo:q_hi,nau,ncell))
        do kap = 1, ncell
          call df_L_kap(nBas, coeff_SAO, p_lo, p_hi, q_lo, q_hi, kap, phase, Lk)
          do k2 = 1, ncell
            L(:,:,:,csub(k2,kap),k2) = Lk(:,:,:,k2)
          end do
        end do
        deallocate(Lk)
        t_L = omp_get_wtime() - w0

        w0 = omp_get_wtime()
        nb_saved = blas_serial()                    ! [Nov 2026] no BLAS threads inside OpenMP
        !$omp parallel do default(shared) private(it,k1,k2,k3,k4) schedule(dynamic,1)
        do it = 0, ncell**3-1
          k1 = it/(ncell*ncell) + 1
          k2 = mod(it/ncell, ncell) + 1
          k3 = mod(it, ncell) + 1
          k4 = cadd(csub(k1,k2), k3)
          call zgemm('N','T', no*nv, no*nv, nau, cmplx(1.d0/dble(ncell),0.d0,kind=8), L(p_lo,q_lo,1,k1,k2), no*nv, &
                     L(p_lo,q_lo,1,k3,k4), no*nv, zero, ERI_out(p_lo,q_lo,p_lo,q_lo,k1,k2,k3), no*nv)
        end do
        !$omp end parallel do
        call blas_restore(nb_saved)
        t_E = omp_get_wtime() - w0
        deallocate(L, phase)

        write(iout,'(a,f10.2,a,f10.2,a)') '  [density fitting] MP2: L^Q_ia(k1,k2) ', t_L, ' s , (ov|ov) ', t_E, ' s'
        write(iout,'(a,f9.3,a,f9.3,a)') '  [density fitting] MP2 memory: (ov|ov) ', &
             16.d0*dble(no*nv)**2*dble(ncell)**3/1024.d0**3, ' GB , L ', &
             16.d0*dble(no*nv)*dble(nau)*dble(ncell)**2/1024.d0**3, ' GB'
        flush(iout)
      end subroutine df_kMO


      !-----------------------------------------------------------------!
      ! L^Q_pq(k1,k2) for one momentum transfer kap = k2 - k1:           !
      !   Lk(p,q,Q,k2) , k1 = k2 - kap                                   !
      !-----------------------------------------------------------------!
      subroutine df_L_kap(nBas, coeff_SAO, p_lo, p_hi, q_lo, q_hi, kap, phase, Lk)
        use unitcell_module, only: nfuc
        integer   , intent(in)  :: nBas, p_lo, p_hi, q_lo, q_hi, kap
        complex(8), intent(in)  :: coeff_SAO(nBas,nBas), phase(df_ncell,df_ncell)
        complex(8), intent(out) :: Lk(p_lo:p_hi,q_lo:q_hi,df_nau,df_ncell)
        integer :: nau, ncell, no, nv, q, rq, k, r, k1, k2, off1, off2
        complex(8), allocatable :: Bk(:,:,:), F(:,:,:,:), tmp(:,:,:), C2(:,:)
        complex(8), parameter :: one=(1.d0,0.d0), zero=(0.d0,0.d0)
        nau = df_nau ; ncell = df_ncell
        no = p_hi - p_lo + 1 ; nv = q_hi - q_lo + 1
        allocate(Bk(nfuc,nBas,nau), F(nfuc,nfuc,nau,ncell))
        ! Bk(a, b, Q) = sum_RQ df_B(a,b,Q RQ) e^{-i kap RQ}
        !$omp parallel do default(shared) private(q,rq)
        do q = 1, nau
          Bk(:,:,q) = zero
          do rq = 1, ncell
            Bk(:,:,q) = Bk(:,:,q) + df_B(:,:,(rq-1)*nau+q)*conjg(phase(kap,rq))
          end do
        end do
        !$omp end parallel do
        ! F(a, b0, Q, k) = sum_R Bk(a, b0+R, Q) e^{ikR}
        !$omp parallel do default(shared) private(k,q,r) collapse(2)
        do k = 1, ncell
          do q = 1, nau
            F(:,:,q,k) = zero
            do r = 1, ncell
              F(:,:,q,k) = F(:,:,q,k) + Bk(:,(r-1)*nfuc+1:r*nfuc,q)*phase(k,r)
            end do
          end do
        end do
        !$omp end parallel do
        deallocate(Bk)
        nb_saved = blas_serial()                    ! [Nov 2026] no BLAS threads inside OpenMP
        !$omp parallel default(shared) private(k2,k1,off1,off2,tmp,C2,q)
        allocate(tmp(nfuc,q_lo:q_hi,nau), C2(nfuc,q_lo:q_hi))
        !$omp do schedule(dynamic,1)
        do k2 = 1, ncell
          k1 = csub(k2, kap)
          off1 = (k1-1)*nfuc ; off2 = (k2-1)*nfuc
          C2 = coeff_SAO(off2+1:off2+nfuc, off2+q_lo:off2+q_hi)
          do q = 1, nau
            call zgemm('N','N', nfuc, nv, nfuc, one, F(1,1,q,k2), nfuc, C2, nfuc, zero, tmp(1,q_lo,q), nfuc)
          end do
          call zgemm('C','N', no, nv*nau, nfuc, one, coeff_SAO(off1+1,off1+p_lo), nBas, tmp, nfuc, zero, &
                     Lk(p_lo,q_lo,1,k2), no)
        end do
        !$omp end do
        deallocate(tmp, C2)
        !$omp end parallel
        call blas_restore(nb_saved)
        deallocate(F)
      end subroutine df_L_kap

      !-----------------------------------------------------------------!
      ! [Nov 2026] density_fitting + ccsdt: the DF vectors replace the   !
      ! Cholesky vectors of module chol_store (same use in the CC):      !
      !   g(p,q,r,s;k1,k2,k3) = 1/N sum_Q L^Q_pq(k1,k2) L^Q_rs(k3,k4)    !
      !                       = sum_Q X^Q_pq(k1,k2) conj(X^Q_sr(k4,k3))  !
      !   X = L/sqrt(N)  (conj L^Q_sr(k4,k3) = L^Q_rs(k3,k4): real B,    !
      !   symmetric pairs), so LK(K)%x(Q,p,q,k2) = X^Q_pq(k2-K,k2).      !
      ! df_B is freed afterwards (no SCF any more).                      !
      !-----------------------------------------------------------------!
      subroutine df_build_chol(nBas, coeff_SAO, iout)
        use unitcell_module, only: nfuc, nx, ny, nz
        use constants_module, only: pi
        use chol_store, only: chol_init_ext, LK, kdiff, chol_bytes
        use omp_lib
        integer   , intent(in) :: nBas, iout
        complex(8), intent(in) :: coeff_SAO(nBas,nBas)
        integer :: ncell, nau, kap, k2, k1, k, r, q, iK
        double precision :: theta, w0, sc
        complex(8), allocatable :: phase(:,:), Lb(:,:,:,:)

        if (.not. df_ready) error stop 'df_build_chol: df_init not called'
        if (allocated(df_Bh)) deallocate(df_Bh)          ! exchange is done: free the k-space B
        w0 = omp_get_wtime()
        call df_rebuild_mp2()
        ncell = df_ncell ; nau = df_nau
        allocate(phase(ncell,ncell))
        do k = 1, ncell
          do r = 1, ncell
            theta = 2.0d0*pi*( dble(df_coord(1,k)*df_coord(1,r))/dble(max(nx,1)) &
                             + dble(df_coord(2,k)*df_coord(2,r))/dble(max(ny,1)) &
                             + dble(df_coord(3,k)*df_coord(3,r))/dble(max(nz,1)) )
            phase(k,r) = exp(cmplx(0.d0,theta,kind=8))
          end do
        end do
        call chol_init_ext(nfuc, nau)
        sc = 1.d0/dsqrt(dble(ncell))
        allocate(Lb(nfuc,nfuc,nau,ncell))
        do kap = 1, ncell
          call df_L_kap(nBas, coeff_SAO, 1, nfuc, 1, nfuc, kap, phase, Lb)
          do k2 = 1, ncell
            k1 = csub(k2, kap)
            iK = kdiff(k2,k1)
            if (iK /= kap) error stop 'df_build_chol: momentum index mismatch'
            do q = 1, nau
              LK(iK)%x(q,:,:,k2) = sc*Lb(:,:,q,k2)
            end do
          end do
        end do
        deallocate(Lb, phase)
        if (allocated(df_B)) deallocate(df_B)
        df_ready = .false.
        write(iout,'(a,i7,a,f10.3,a,f9.2,a)') '  [density fitting] CC vectors: ', nau, ' per momentum transfer, ', &
             chol_bytes()/1024.d0**3, ' GB , time ', omp_get_wtime()-w0, ' s'
        flush(iout)
      end subroutine df_build_chol

      !-----------------------------------------------------------------!
      ! optional file AuxBasis_MP2 (e.g. the RI/MP2-fit basis, as PySCF's !
      ! molecular DF-MP2): B is rebuilt with it before the MP2 step      !
      !-----------------------------------------------------------------!
      subroutine df_rebuild_mp2()
        use atom_basis
        logical :: ex
        integer :: i
        type(atom), allocatable :: at(:)
        inquire(file='AuxBasis_MP2', exist=ex)
        if (.not. ex) return
        allocate(at(size(sv_charge)))
        do i = 1, size(sv_charge)
          at(i)%charge = sv_charge(i)
        end do
        call df_init(sv_ndp, size(sv_charge), at, sv_geom, size(sv_AO), sv_AO, sv_iout, 'AuxBasis_MP2')
        deallocate(at)
      end subroutine df_rebuild_mp2





      !-----------------------------------------------------------------!
      ! [Nov 2026] pivoted Cholesky of the metric V(k) (Hermitian PSD):      !
      !   V(piv,piv) ~ L L^H , rank m where the largest remaining        !
      !   diagonal (Schur complement) < tol                              !
      !   f(piv(i), j) = [L^-H](i,j) , j <= m ; f = 0 otherwise          !
      ! so B = T3 f fits with the m pivot functions only:  f f^H = V_pp^-1
      !-----------------------------------------------------------------!
      subroutine chol_f(n, V, tol, f, ndrop, dmin)
        integer, intent(in) :: n
        complex(8), intent(inout) :: V(n,n)
        double precision, intent(in) :: tol
        complex(8), intent(out) :: f(n,n)
        integer, intent(out) :: ndrop
        double precision, intent(inout) :: dmin
        integer :: piv(n), m, info, i, j
        double precision :: work(2*n)
        complex(8), allocatable :: L(:,:)
        call zpstrf('L', n, V, n, piv, m, tol, work, info)
        if (info < 0) error stop 'chol_f: zpstrf failed'
        allocate(L(m,m)) ; L = (0.d0,0.d0)
        do j = 1, m
          do i = j, m
            L(i,j) = V(i,j)
          end do
        end do
        dmin = min(dmin, dble(L(m,m))**2)
        call ztrtri('L', 'N', m, L, m, info)            ! L^-1
        if (info /= 0) error stop 'chol_f: ztrtri failed'
        f = (0.d0,0.d0)
        do j = 1, m
          do i = 1, m
            f(piv(i), j) = conjg(L(j,i))                 ! (L^-1)^H
          end do
        end do
        ndrop = n - m
        deallocate(L)
      end subroutine chol_f



      !-----------------------------------------------------------------!
      ! [Nov 2026] Cholesky-based fitting set: the auxiliary "functions" are !
      ! orbital products I = (p0 q) (p in cell 0, q anywhere) chosen by  !
      ! a pivoted Cholesky of the Gamma-point lattice-summed Coulomb     !
      ! matrix  V(I,J) = sum_R (I_-R | J_0)  down to the tolerance tol.  !
      ! The set is replicated in every cell, so the usual DF pipeline    !
      ! (per-k metric, B in k-space) runs unchanged. tol -> 0 gives the  !
      ! exact integrals (all products are then in the fitting set).     !
      ! Returns Vm(P0, Q R) and T3(mu0, nu, Q R) = (mu0 nu | Q_R).       !
      !-----------------------------------------------------------------!
      subroutine cd_build(ndp, nBas, AO, ncell, tol, iout, nau, Vout, Tout)
        use unitcell_module, only: nfuc
        use eri_direct     , only: d_nsh, d_nshuc, d_sh_first, d_sh_n, d_Q, d_Qsh, set_ps
        use clifford_tools , only: ps_on
        use keywords       , only: c_prim_screen
        use omp_lib
        integer, intent(in) :: ndp, nBas, ncell, iout
        type(ERI_function), intent(in) :: AO(nBas)
        double precision, intent(in) :: tol
        integer, intent(out) :: nau
        double precision, allocatable, intent(out) :: Vout(:,:), Tout(:,:,:)
        type colT
          double precision, allocatable :: a(:,:)
        end type colT
        type(colT), allocatable :: col(:)
        integer, allocatable :: sh_of(:), trT(:,:), piv(:), sp(:,:)
        double precision, allocatable :: D(:), L(:,:), c(:), blk(:,:,:,:)
        integer :: ncand, i, r, rn, fn, ci, i0, SA, SB, SAp, SBp, a, b, k, J, ii, p, q, maxr, nsp, it, P_, Q_, nc, ncols
        double precision :: w0, dj, t_d, t_c, t_u, t1
        integer :: nsp_done, nreuse
        logical, allocatable :: done(:,:)

        w0 = omp_get_wtime()
        allocate(sh_of(nBas))
        do SA = 1, d_nsh
          do a = 1, d_sh_n(SA)
            sh_of(d_sh_first(SA)+a-1) = SA
          end do
        end do
        ! trT(fn, r) = function fn translated by -R(r)
        allocate(trT(nBas, ncell))
        do r = 1, ncell
          rn = csub(1, r)
          do fn = 1, nBas
            ci = (fn-1)/nfuc + 1 ; i0 = mod(fn-1, nfuc) + 1
            trT(fn, r) = (cadd(ci, rn)-1)*nfuc + i0
          end do
        end do
        ncand = nfuc*nBas
        allocate(D(ncand)) ; D = 0.d0

        ! ---- diagonal D(ii) = sum_R (I_-R | I_0) ---- !
        allocate(sp(2, d_nshuc*d_nsh)) ; nsp = 0
        do SA = 1, d_nshuc
          do SB = 1, d_nsh
            if (dble(ncell)*d_Qsh(SA,SB)**2 < tol) cycle
            nsp = nsp + 1 ; sp(:,nsp) = (/ SA, SB /)
          end do
        end do
        !$omp parallel do default(shared) private(it,SA,SB,r,SAp,SBp,blk,a,b,ii) schedule(dynamic,1)
        do it = 1, nsp
          SA = sp(1,it) ; SB = sp(2,it)
          allocate(blk(d_sh_n(SA),d_sh_n(SB),d_sh_n(SA),d_sh_n(SB)))
          do r = 1, ncell
            SAp = sh_of(trT(d_sh_first(SA), r)) ; SBp = sh_of(trT(d_sh_first(SB), r))
            call ERI_shell_block_toroidal(ndp, d_sh_n(SA), AO(d_sh_first(SAp):), d_sh_n(SB), AO(d_sh_first(SBp):), &
                                               d_sh_n(SA), AO(d_sh_first(SA):),  d_sh_n(SB), AO(d_sh_first(SB):), blk)
            do b = 1, d_sh_n(SB)
              do a = 1, d_sh_n(SA)
                ii = (d_sh_first(SB)+b-2)*nfuc + d_sh_first(SA)+a-1
                D(ii) = D(ii) + blk(a,b,a,b)
              end do
            end do
          end do
          deallocate(blk)
        end do
        !$omp end parallel do
        deallocate(sp)
        t_d = omp_get_wtime()-w0 ; t_c = 0 ; t_u = 0 ; nsp_done = 0

        ! ---- pivoted Cholesky ---- !
        maxr = min(ncand, 20000)
        allocate(L(ncand, 64), piv(maxr), c(ncand), col(ncand), done(d_nshuc, d_nsh))
        done = .false. ; nreuse = 0
        k = 0 ; ncols = 0
        do
          J = maxloc(D, 1) ; dj = D(J)
          if (dj < tol .or. k == maxr) exit
          t1 = omp_get_wtime()
          if (.not. allocated(col(J)%a)) call shellpair_cols(J)
          t_c = t_c + omp_get_wtime()-t1 ; t1 = omp_get_wtime()
          ! c(ii) = sum_R col_J(p_-R, q_-R)
          !$omp parallel do default(shared) private(ii,p,q,r)
          do ii = 1, ncand
            p = mod(ii-1,nfuc) + 1 ; q = (ii-1)/nfuc + 1
            c(ii) = 0.d0
            do r = 1, ncell
              c(ii) = c(ii) + col(J)%a(trT(p,r), trT(q,r))
            end do
          end do
          !$omp end parallel do
          k = k + 1
          if (k > size(L,2)) call grow_L()
          if (k > 1) call dgemv('N', ncand, k-1, -1.d0, L, ncand, L(J,1), ncand, 1.d0, c, 1)
          L(:,k) = c/dsqrt(dj)
          D = D - L(:,k)**2
          D(J) = 0.d0
          piv(k) = J
          t_u = t_u + omp_get_wtime()-t1
        end do
        nau = k
        write(iout,'(a,i10,a,i10,a,es9.2)') ' Cholesky fitting set: ', nau, ' products per cell (of', ncand, &
                                            ' candidates), largest remaining diagonal ', maxval(D)
        write(iout,'(a,i10,a,i10,a)') ' columns of integrals             : ', ncols, '  (', nsp_done, ' shell pairs)'
        write(iout,'(a,i10)') ' shell blocks copied, not computed: ', nreuse
        write(iout,'(a,3f9.2,a)') ' wall time diagonal/columns/update: ', t_d, t_c, t_u, ' s'

        ! ---- Vm(P0, Q R) = (P_-R | Q_0) , T3(mu0, nu, Q R) = (mu_-R nu_-R | Q_0) ---- !
        allocate(Vout(nau, nau*ncell), Tout(nfuc, nBas, nau*ncell))
        !$omp parallel do default(shared) private(Q_,r,J,P_,ii,p,q,i,fn) collapse(2)
        do r = 1, ncell
          do Q_ = 1, nau
            J = piv(Q_)
            do P_ = 1, nau
              ii = piv(P_) ; p = mod(ii-1,nfuc) + 1 ; q = (ii-1)/nfuc + 1
              Vout(P_, (r-1)*nau+Q_) = col(J)%a(trT(p,r), trT(q,r))
            end do
            do fn = 1, nBas
              do i = 1, nfuc
                Tout(i, fn, (r-1)*nau+Q_) = col(J)%a(trT(i,r), trT(fn,r))
              end do
            end do
          end do
        end do
        !$omp end parallel do
        do ii = 1, ncand
          if (allocated(col(ii)%a)) deallocate(col(ii)%a)
        end do
        deallocate(col, L, D, c, piv, trT, sh_of, done)
        write(iout,'(a,f10.2,a)') ' wall time Cholesky selection     : ', omp_get_wtime()-w0, ' s'

      contains

        subroutine grow_L()
          double precision, allocatable :: t(:,:)
          allocate(t(ncand, 2*size(L,2)))
          t(:,1:size(L,2)) = L
          call move_alloc(t, L)
        end subroutine grow_L

        ! all columns (a b | c d)_{a,b in supercell} for the shell pair of the product J
        ! [Nov 2026] a block (S1 S2|SC SD) whose bra pair, moved to cell 0, was already
        ! a pivot shell pair is copied from that pair's columns (bra-ket symmetry);
        ! the others are computed with the quartet moved so that S1 is in cell 0
        ! (needed for the primitive screening tables of keyword prim_screen)
        subroutine shellpair_cols(Jp)
          integer, intent(in) :: Jp
          integer :: SC, SD, nC, nD, cc, dd, JJ, S1, S2, x, y, r1, T1, T2, TC, TD, KK, fx, fy, nr, r2, U1, U2, UC, UD
          double precision, allocatable :: bk(:,:,:,:)
          SC = sh_of(mod(Jp-1,nfuc)+1) ; SD = sh_of((Jp-1)/nfuc+1)
          nsp_done = nsp_done + 1
          nC = d_sh_n(SC) ; nD = d_sh_n(SD)
          do dd = 1, nD
            do cc = 1, nC
              JJ = (d_sh_first(SD)+dd-2)*nfuc + d_sh_first(SC)+cc-1
              allocate(col(JJ)%a(nBas,nBas)) ; col(JJ)%a = 0.d0
              ncols = ncols + 1
            end do
          end do
          ps_on = c_prim_screen
          nr = 0
          !$omp parallel do default(shared) private(S1,S2,bk,x,y,cc,dd,JJ,r1,T1,T2,TC,TD,KK,fx,fy,r2,U1,U2,UC,UD) &
          !$omp reduction(+:nr) schedule(dynamic,1)
          do S1 = 1, d_nsh
            do S2 = S1, d_nsh
              if (d_Qsh(S1,S2)*d_Qsh(SC,SD) < 1.d-15) cycle
              allocate(bk(d_sh_n(S1),d_sh_n(S2),nC,nD))
              ! move the quartet so that S1 is in cell 0
              r1 = (d_sh_first(S1)-1)/nfuc + 1
              T1 = sh_of(trT(d_sh_first(S1), r1)) ; T2 = sh_of(trT(d_sh_first(S2), r1))
              TC = sh_of(trT(d_sh_first(SC), r1)) ; TD = sh_of(trT(d_sh_first(SD), r1))
              if (done(T1, T2)) then
                ! (S1 S2|SC SD) = (TC TD|T1 T2) = column (x0 y) of pair (T1,T2) at (c-R, d-R)
                do y = 1, d_sh_n(S2)
                  do x = 1, d_sh_n(S1)
                    KK = (d_sh_first(T2)+y-2)*nfuc + d_sh_first(T1)+x-1
                    do dd = 1, nD
                      do cc = 1, nC
                        bk(x,y,cc,dd) = col(KK)%a(d_sh_first(TC)+cc-1, d_sh_first(TD)+dd-1)
                      end do
                    end do
                  end do
                end do
                nr = nr + 1
              else
                ! the same with S2 moved to cell 0: (S1 S2|SC SD) = (UC UD|U1 U2), U1 = S2-R2, U2 = S1-R2
                r2 = (d_sh_first(S2)-1)/nfuc + 1
                U1 = sh_of(trT(d_sh_first(S2), r2)) ; U2 = sh_of(trT(d_sh_first(S1), r2))
                UC = sh_of(trT(d_sh_first(SC), r2)) ; UD = sh_of(trT(d_sh_first(SD), r2))
                if (done(U1, U2)) then
                  do y = 1, d_sh_n(S2)
                    do x = 1, d_sh_n(S1)
                      KK = (d_sh_first(U2)+x-2)*nfuc + d_sh_first(U1)+y-1
                      do dd = 1, nD
                        do cc = 1, nC
                          bk(x,y,cc,dd) = col(KK)%a(d_sh_first(UC)+cc-1, d_sh_first(UD)+dd-1)
                        end do
                      end do
                    end do
                  end do
                  nr = nr + 1
                else
                if (ps_on) call set_ps(T1, T2, TC, TD)
                call ERI_shell_block_toroidal(ndp, d_sh_n(S1), AO(d_sh_first(T1):), d_sh_n(S2), AO(d_sh_first(T2):), &
                                                   nC, AO(d_sh_first(TC):), nD, AO(d_sh_first(TD):), bk)
                end if
              end if
              do dd = 1, nD
                do cc = 1, nC
                  JJ = (d_sh_first(SD)+dd-2)*nfuc + d_sh_first(SC)+cc-1
                  do y = 1, d_sh_n(S2)
                    fy = d_sh_first(S2)+y-1
                    do x = 1, d_sh_n(S1)
                      fx = d_sh_first(S1)+x-1
                      col(JJ)%a(fx, fy) = bk(x,y,cc,dd)
                      col(JJ)%a(fy, fx) = bk(x,y,cc,dd)
                    end do
                  end do
                end do
              end do
              deallocate(bk)
            end do
          end do
          !$omp end parallel do
          ps_on = .false.
          nreuse = nreuse + nr
          done(SC, SD) = .true.
        end subroutine shellpair_cols

      end subroutine cd_build


      !-----------------------------------------------------------------!
      ! [Nov 2026] OpenBLAS (pthreads build, e.g. Homebrew) starts its   !
      ! own threads in every zgemm, also inside our OpenMP regions:      !
      ! threads x threads oversubscription. Switch BLAS to one thread    !
      ! around those regions and restore afterwards.                     !
      !-----------------------------------------------------------------!
      integer function blas_serial()
        interface
          integer(c_int) function openblas_get_num_threads() bind(C, name='openblas_get_num_threads')
            import :: c_int
          end function openblas_get_num_threads
          subroutine openblas_set_num_threads(n) bind(C, name='openblas_set_num_threads')
            import :: c_int
            integer(c_int), value :: n
          end subroutine openblas_set_num_threads
        end interface
        blas_serial = int(openblas_get_num_threads())
        call openblas_set_num_threads(1_c_int)
      end function blas_serial

      subroutine blas_restore(n)
        integer, intent(in) :: n
        interface
          subroutine openblas_set_num_threads(m) bind(C, name='openblas_set_num_threads')
            import :: c_int
            integer(c_int), value :: m
          end subroutine openblas_set_num_threads
        end interface
        call openblas_set_num_threads(int(max(n,1), c_int))
      end subroutine blas_restore


      !=================================================================!
      !                                                                 !
      ! [Nov 2026] density_fitting + symmetry (HF)                      !
      !                                                                 !
      ! Same fitting as the default Cholesky fitting set (cd_build),    !
      ! but the set of products is closed under the point operations   !
      ! (whole orbits are added), and the k-space B is stored only for  !
      ! the irreducible k2 (point operations + time reversal):          !
      !                                                                 !
      !   Bh(a,q,b;k1,k2) = sum_{R1,R2} B(a0, b R1, q R2) e^{i k1 R1}   !
      !                                               e^{i k2 R2}       !
      !   B = T3 V^-1/2 (symmetric V^-1/2 per k2: covariant)            !
      !                                                                 !
      ! An operation g: a0 -> (pa, ca), sign sa (AO), q0 -> (pq, cq),   !
      ! sign sq (aux product) gives                                     !
      !   Bh(pa,pq,pb; Rk1,Rk2) = sa sb sq e^{i Rk1.(cb-ca)}            !
      !                           e^{i Rk2.(cq-ca)} Bh(a,q,b;k1,k2)     !
      !   Bh(k1,k2) = conj Bh(-k1,-k2)                                  !
      ! so every k2 slice is rebuilt from a stored one when K needs it. !
      !                                                                 !
      ! Memory: Bh for the irreducible k2 only                          !
      !   (nfuc * nau * nfuc * N_cell * n_irr complex) + B(k2=0) real.  !
      ! No real-space B (nfuc*nBas*nau*N_cell), no T3 copies, no full   !
      ! columns kept: the columns of the orbit representatives are      !
      ! computed twice (selection, then T3), images by permutation.     !
      !                                                                 !
      ! J = from B(k2 = 0) only, K = one k2 slice at a time.            !
      ! Exact for any density (no symmetric-P requirement).             !
      !=================================================================!
      subroutine df_init_sym(ndp, nBas, AO, ncell, phase, tol, iout)
        use unitcell_module  , only: nfuc
        use eri_direct       , only: d_nsh, d_nshuc, d_sh_first, d_sh_n, d_Q, d_Qsh, set_ps
        use clifford_tools   , only: ps_on
        use keywords         , only: c_prim_screen
        use symmetry_clifford, only: sym_nop, sym_perm, sym_sign, sym_R
        use omp_lib
        integer           , intent(in) :: ndp, nBas, ncell, iout
        type(ERI_function), intent(in) :: AO(nBas)
        complex(8)        , intent(in) :: phase(ncell,ncell)
        double precision  , intent(in) :: tol

        integer :: nG, ncand, ii, p, q, r, g, a, b, k, kk, m, nr, i, j, SA, SB, SC, SD, nC, nD, cc, dd, it
        integer :: nsp, JJ, nsl, ksl, npiv, kL, maxp, nb, ik, nirr, mm, xa_, ya_, xr_, yr_, cshift, s1, ninfo
        integer, allocatable :: fc(:), fl(:), trT(:,:), sh_of(:), invp(:,:), imgr(:,:), cslot(:), sp(:,:)
        integer(1), allocatable :: isg(:,:)
        integer, allocatable :: pcand(:), prep(:), pg(:), pswap(:), pividx(:), mlist(:), mg(:)
        logical, allocatable :: known(:), isrep_sp(:), sel(:)
        double precision, allocatable :: D(:), L(:,:), cst(:,:), cvec(:), blk(:,:,:,:), colw(:,:,:)
        double precision, allocatable :: X(:,:), Vout(:,:,:), Xr(:,:), Xi(:,:), Cph(:,:), Sph(:,:), ew(:), rwork(:)
        complex(8), allocatable :: Vk(:,:), Vih(:,:), Wk(:,:), Wt(:,:,:,:), work(:)
        double precision :: dj, w0, t_d, t_a, t_b, t_c, dmax0, dep_tol, eig_cut, lmin, lmax
        integer :: nkeep, ndep, nrep_sp, nrep_c, lwork, info
        logical :: isnew

        w0 = omp_get_wtime()
        nG = sym_nop ; ncand = nfuc*nBas ; nb = ncand
        write(iout,'(a)') ''
        write(iout,'(a)') '------------------------------------------------------------'
        write(iout,'(a)') ' Density fitting with symmetry (keywords density_fitting + symmetry)'
        write(iout,'(a)') '------------------------------------------------------------'
        write(iout,'(a,i0,a)') ' point operations used: ', nG, ' (with their translations)'
        flush(iout)

        ! ---- cell tables ---- !
        allocate(fc(nBas), fl(nBas), trT(nBas,ncell), sh_of(nBas), invp(nBas,nG))
        do m = 1, nBas
          fc(m) = (m-1)/nfuc + 1 ; fl(m) = mod(m-1,nfuc) + 1
        end do
        do r = 1, ncell
          do m = 1, nBas
            trT(m,r) = (csub(fc(m), r)-1)*nfuc + fl(m)          ! m shifted by -R(r)
          end do
        end do
        do SA = 1, d_nsh
          sh_of(d_sh_first(SA):d_sh_first(SA)+d_sh_n(SA)-1) = SA
        end do
        do g = 1, nG
          do m = 1, nBas
            invp(sym_perm(m,g), g) = m
          end do
        end do

        ! ---- candidate products ii = (q-1)*nfuc + p : images and swap ---- !
        allocate(imgr(ncand,nG), isg(ncand,nG))
        !$omp parallel do default(shared) private(ii,p,q,g,a,b) schedule(static)
        do ii = 1, ncand
          p = mod(ii-1,nfuc) + 1 ; q = (ii-1)/nfuc + 1
          do g = 1, nG
            a = sym_perm(p,g) ; b = sym_perm(q,g)
            imgr(ii,g) = (trT(b, fc(a))-1)*nfuc + fl(a)
            isg(ii,g)  = int(sym_sign(p,g)*sym_sign(q,g), 1)
          end do
        end do
        !$omp end parallel do

        ! ---- diagonal D(ii) = sum_R (I_-R | I_0), shell pairs representative under G ---- !
        t_d = omp_get_wtime()
        allocate(D(ncand), known(ncand), isrep_sp(d_nshuc*d_nsh))
        D = 0.d0 ; known = .false. ; isrep_sp = .false.
        allocate(sp(2, d_nshuc*d_nsh)) ; nsp = 0 ; nrep_sp = 0
        do SA = 1, d_nshuc
          do SB = 1, d_nsh
            if (dble(ncell)*d_Qsh(SA,SB)**2 < tol) cycle
            if (.not. sp_is_rep(SA, SB)) cycle
            nsp = nsp + 1 ; sp(:,nsp) = (/ SA, SB /)
          end do
        end do
        nrep_sp = nsp
        !$omp parallel do default(shared) private(it,SA,SB,r,SC,SD,blk,a,b,ii) schedule(dynamic,1)
        do it = 1, nsp
          SA = sp(1,it) ; SB = sp(2,it)
          allocate(blk(d_sh_n(SA),d_sh_n(SB),d_sh_n(SA),d_sh_n(SB)))
          do r = 1, ncell
            SC = sh_of(trT(d_sh_first(SA), r)) ; SD = sh_of(trT(d_sh_first(SB), r))
            call ERI_shell_block_toroidal(ndp, d_sh_n(SA), AO(d_sh_first(SC):), d_sh_n(SB), AO(d_sh_first(SD):), &
                                               d_sh_n(SA), AO(d_sh_first(SA):), d_sh_n(SB), AO(d_sh_first(SB):), blk)
            do b = 1, d_sh_n(SB)
              do a = 1, d_sh_n(SA)
                ii = (d_sh_first(SB)+b-2)*nfuc + d_sh_first(SA)+a-1
                D(ii) = D(ii) + blk(a,b,a,b)
              end do
            end do
          end do
          do b = 1, d_sh_n(SB)
            do a = 1, d_sh_n(SA)
              known((d_sh_first(SB)+b-2)*nfuc + d_sh_first(SA)+a-1) = .true.
            end do
          end do
          deallocate(blk)
        end do
        !$omp end parallel do
        deallocate(sp)
        ! the other candidates: images (and swaps) of the computed ones; candidates of shell pairs
        ! below the Schwarz bound stay 0
        do ii = 1, ncand
          if (.not. known(ii)) cycle
          if (.not. (sp_rep_of(ii))) cycle
          do g = 1, nG
            r = imgr(ii,g)
            if (.not. known(r)) then
              D(r) = D(ii) ; known(r) = .true.
            end if
            r = swapc(r)
            if (.not. known(r)) then
              D(r) = D(ii) ; known(r) = .true.
            end if
          end do
        end do
        do ii = 1, ncand                                   ! only canonical products are candidates
          if (swapc(ii) < ii) D(ii) = 0.d0
        end do
        t_d = omp_get_wtime() - t_d
        dmax0 = maxval(D)

        ! ---- pivoted Cholesky with whole orbits ---- !
        t_a = omp_get_wtime()
        dep_tol = 1.d-3*tol
        maxp = ncand                                   ! at most every canonical product
        allocate(L(ncand, 64), cvec(ncand), cslot(ncand), cst(ncand, 16), sel(ncand))
        allocate(pcand(maxp), prep(maxp), pg(maxp), pswap(maxp), pividx(ncand), mlist(nG), mg(nG))
        cslot = 0 ; nsl = 0 ; sel = .false. ; pividx = 0
        kL = 0 ; npiv = 0 ; ndep = 0 ; nrep_c = 0
        do
          JJ = maxloc(D, 1) ; dj = D(JJ)
          if (dj < tol) exit
          ! c-vectors of the shell pair of JJ (computed once, kept)
          if (cslot(JJ) == 0) call pair_cvecs(JJ)
          nrep_c = nrep_c + 1
          ! orbit members (canonical), JJ first
          mm = 0
          do g = 1, nG
            r = imgr(JJ,g) ; m = min(r, swapc(r))
            if (sel(m)) cycle
            isnew = .true.
            do i = 1, mm
              if (mlist(i) == m) isnew = .false.
            end do
            if (.not. isnew) cycle
            mm = mm + 1 ; mlist(mm) = m ; mg(mm) = g
          end do
          do i = 1, mm
            m = mlist(i) ; g = mg(i)
            ! c-vector of the member: permutation of the one of JJ
            !$omp parallel do default(shared) private(ii) schedule(static)
            do ii = 1, ncand
              cvec(imgr(ii,g)) = dble(isg(ii,g)*isg(JJ,g))*cst(ii, cslot(JJ))
            end do
            !$omp end parallel do
            if (kL > 0) call dgemv('N', ncand, kL, -1.d0, L, ncand, L(m,1), ncand, 1.d0, cvec, 1)
            npiv = npiv + 1
            pcand(npiv) = m ; prep(npiv) = JJ ; pg(npiv) = g
            pswap(npiv) = merge(1, 0, m /= imgr(JJ,g))
            pividx(m) = npiv ; sel(m) = .true.
            if (D(m) > dep_tol) then
              kL = kL + 1
              if (kL > size(L,2)) call grow2(L)
              L(:,kL) = cvec/dsqrt(D(m))
              D = D - L(:,kL)**2
            else
              ndep = ndep + 1
            end if
            D(m) = 0.d0
          end do
        end do
        t_a = omp_get_wtime() - t_a
        deallocate(L, cvec, cst, cslot, D, known)
        write(iout,'(a,i10)')     ' products in the fitting set (per cell): ', npiv
        write(iout,'(a,i10,a)')   '   of which nearly dependent orbit members: ', ndep, ' (kept for the symmetry)'
        write(iout,'(a,i10)')     ' orbits selected (columns computed)    : ', nrep_c
        write(iout,'(a,i10)')     ' representative shell pairs (diagonal) : ', nrep_sp
        write(iout,'(a,es10.2)')  ' Cholesky tolerance                    : ', tol
        flush(iout)

        ! ---- irreducible k points ---- !
        call df_kirr(ncell, nG, nirr)
        write(iout,'(a,i0,a,i0)') ' k points: ', ncell, '   irreducible (point operations + time reversal): ', nirr
        flush(iout)

        ! ---- aux maps for the operations ---- !
        if (allocated(dfs_aqp)) deallocate(dfs_aqp, dfs_aqc, dfs_aqs, dfs_aop, dfs_aoc, dfs_aos)
        allocate(dfs_aqp(npiv,nG), dfs_aqc(npiv,nG), dfs_aqs(npiv,nG))
        allocate(dfs_aop(nfuc,nG), dfs_aoc(nfuc,nG), dfs_aos(nfuc,nG))
        ninfo = 0
        do g = 1, nG
          do a = 1, nfuc
            dfs_aop(a,g) = fl(sym_perm(a,g)) ; dfs_aoc(a,g) = fc(sym_perm(a,g))
            dfs_aos(a,g) = int(sym_sign(a,g),1)
          end do
          do k = 1, npiv
            ii = pcand(k) ; p = mod(ii-1,nfuc) + 1 ; q = (ii-1)/nfuc + 1
            r = imgr(ii,g) ; m = min(r, swapc(r))
            dfs_aqp(k,g) = pividx(m)
            if (pividx(m) == 0) ninfo = ninfo + 1
            if (m == r) then
              dfs_aqc(k,g) = fc(sym_perm(p,g))
            else
              dfs_aqc(k,g) = fc(sym_perm(q,g))
            end if
            dfs_aqs(k,g) = isg(ii,g)
          end do
        end do
        if (ninfo > 0) error stop 'density_fitting + symmetry: the fitting set is not closed (bug)'

        ! ---- T3(k2) for the irreducible k2, V(R) ---- !
        t_b = omp_get_wtime()
        if (allocated(dfs_big)) deallocate(dfs_big)
        allocate(dfs_big(nb*npiv, nirr)) ; dfs_big = (0.d0,0.d0)
        allocate(Vout(npiv, npiv, ncell), X(nb, ncell), Xr(nb, nirr), Xi(nb, nirr), Cph(ncell,nirr), Sph(ncell,nirr))
        do ik = 1, nirr
          do r = 1, ncell
            Cph(r,ik) = dble(phase(dfs_kirr(ik),r)) ; Sph(r,ik) = aimag(phase(dfs_kirr(ik),r))
          end do
        end do
        ! representative candidates (one column computation per shell pair)
        if (allocated(colw)) deallocate(colw)
        it = 0                                                ! shell pair of the columns in colw
        do k = 1, npiv
          if (prep(k) /= pcand(k) .or. pg(k) /= 1) cycle      ! k is the representative of its orbit
          JJ = pcand(k)
          SC = sh_of(mod(JJ-1,nfuc)+1) ; SD = sh_of((JJ-1)/nfuc+1)
          if (it /= (SD-1)*d_nsh + SC) then
            call pair_cols(JJ)                                ! colw(:,:,component of JJ)
            it = (SD-1)*d_nsh + SC
          end if
          cc = mod(JJ-1,nfuc)+1 - d_sh_first(SC) + 1
          dd = (JJ-1)/nfuc+1 - d_sh_first(SD) + 1
          s1 = (dd-1)*d_sh_n(SC) + cc
          do kk = 1, npiv
            if (prep(kk) /= JJ) cycle
            g = pg(kk) ; ii = pcand(kk)
            p = mod(JJ-1,nfuc) + 1 ; q = (JJ-1)/nfuc + 1
            if (pswap(kk) == 0) then
              cshift = fc(sym_perm(p,g))
            else
              cshift = fc(sym_perm(q,g))
            end if
            ! X(mu0 nu, R) = (mu_-R nu_-R | m_0) = s_x s_y s_J col_J(x,y), x = g^-1(x' + c)
            !$omp parallel do default(shared) private(r,ii,p,q,xr_,yr_,xa_,ya_) schedule(static)
            do r = 1, ncell
              do ii = 1, ncand
                p = mod(ii-1,nfuc) + 1 ; q = (ii-1)/nfuc + 1
                xr_ = trT(p, r) ; yr_ = trT(q, r)
                xa_ = invp((cadd(fc(xr_),cshift)-1)*nfuc + fl(xr_), g)
                ya_ = invp((cadd(fc(yr_),cshift)-1)*nfuc + fl(yr_), g)
                X(ii, r) = dble(sym_sign(xa_,g)*sym_sign(ya_,g)*isg(JJ,g))*colw(xa_, ya_, s1)
              end do
            end do
            !$omp end parallel do
            do kL = 1, npiv
              Vout(kL, kk, :) = X(pcand(kL), :)
            end do
            call dgemm('N','N', nb, nirr, ncell, 1.d0, X, nb, Cph, ncell, 0.d0, Xr, nb)
            call dgemm('N','N', nb, nirr, ncell, 1.d0, X, nb, Sph, ncell, 0.d0, Xi, nb)
            do ik = 1, nirr
              dfs_big((kk-1)*nb+1:kk*nb, ik) = cmplx(Xr(:,ik), Xi(:,ik), kind=8)
            end do
          end do
        end do
        if (allocated(colw)) deallocate(colw)
        deallocate(X, Xr, Xi, Cph, Sph)
        t_b = omp_get_wtime() - t_b

        ! ---- per irreducible k2: V^-1/2 (eigenvalues > eig_cut), W = T3 V^-1/2, Bh slice ---- !
        t_c = omp_get_wtime()
        eig_cut = 1.d-2*tol
        allocate(Vk(npiv,npiv), Vih(npiv,npiv), ew(npiv), rwork(3*npiv), work(1))
        call zheev('V','U', npiv, Vk, npiv, ew, work, -1, rwork, info)
        lwork = max(1, int(dble(work(1)))) ; deallocate(work) ; allocate(work(lwork))
        allocate(Wk(nb,npiv), Wt(nfuc,npiv,nfuc,ncell))
        if (allocated(dfs_B0)) deallocate(dfs_B0)
        allocate(dfs_B0(nb, npiv))
        lmin = huge(1.d0) ; lmax = 0.d0 ; nkeep = npiv
        do ik = 1, nirr
          kk = dfs_kirr(ik)
          Vk = (0.d0,0.d0)
          do r = 1, ncell
            Vk = Vk + Vout(:,:,r)*phase(kk,r)
          end do
          Vk = 0.5d0*(Vk + conjg(transpose(Vk)))
          call zheev('V','U', npiv, Vk, npiv, ew, work, lwork, rwork, info)
          if (info /= 0) error stop 'df_init_sym: zheev of V(k) failed'
          Vih = (0.d0,0.d0)
          do i = 1, npiv
            if (ew(i) > eig_cut) then
              lmin = min(lmin, ew(i))
              do j = 1, npiv
                Vih(:,j) = Vih(:,j) + Vk(:,i)*conjg(Vk(j,i))/dsqrt(ew(i))
              end do
            end if
          end do
          nkeep = min(nkeep, count(ew > eig_cut)) ; lmax = max(lmax, maxval(ew))
          call zgemm('N','N', nb, npiv, npiv, (1.d0,0.d0), dfs_big(1,ik), nb, Vih, npiv, (0.d0,0.d0), Wk, nb)
          if (kk == 1) dfs_B0 = dble(Wk)                     ! k2 = 0: real
          ! Wt(mu,Q,nu0,R1) = W(mu0, nu0 R1, Q) ; Bh(.,k1) = sum_R1 Wt(.,R1) e^{i k1 R1}
          !$omp parallel do default(shared) private(r,b,q) collapse(2)
          do r = 1, ncell
            do b = 1, nfuc
              do q = 1, npiv
                Wt(:,q,b,r) = Wk((r-1)*nfuc*nfuc + (b-1)*nfuc + 1 : (r-1)*nfuc*nfuc + b*nfuc, q)
              end do
            end do
          end do
          !$omp end parallel do
          call zgemm('N','T', nfuc*npiv*nfuc, ncell, ncell, (1.d0,0.d0), Wt, nfuc*npiv*nfuc, phase, ncell, &
                     (0.d0,0.d0), dfs_big(1,ik), nfuc*npiv*nfuc)
        end do
        deallocate(Vk, Vih, ew, rwork, work, Wk, Wt, Vout)
        t_c = omp_get_wtime() - t_c

        dfs_npiv = npiv ; dfs_nirr = nirr ; df_nau = npiv ; df_naux = npiv*ncell
        df_sym = .true. ; df_ready = .true.
        deallocate(fc, fl, trT, sh_of, invp, imgr, isg, isrep_sp, sel, pcand, prep, pg, pswap, pividx, mlist, mg)
        write(iout,'(a,2es12.4)') ' eigenvalues of V(k) kept (min, max)  : ', lmin, lmax
        write(iout,'(a,es10.2,a,i0)') ' eigenvalue cut ', eig_cut, ' ; smallest number kept over k: ', nkeep
        write(iout,'(a,f10.2,a)') ' wall time diagonal                   : ', t_d, ' s'
        write(iout,'(a,f10.2,a)') ' wall time selection (columns)        : ', t_a, ' s'
        write(iout,'(a,f10.2,a)') ' wall time T3(k) (columns again)      : ', t_b, ' s'
        write(iout,'(a,f10.2,a)') ' wall time V^-1/2 and B(k1,k2)        : ', t_c, ' s'
        write(iout,'(a,f10.3,a)') ' memory of B(k1,k2), irreducible k2   : ', &
             16.d0*dble(nb)*dble(npiv)*dble(nirr)/1024.d0**3, ' GB'
        write(iout,'(a,f10.3,a)') '   (all k2 would be                   : ', &
             16.d0*dble(nb)*dble(npiv)*dble(ncell)/1024.d0**3, ' GB)'
        write(iout,'(a,f10.2,a)') ' wall time density fitting setup      : ', omp_get_wtime()-w0, ' s'
        write(iout,'(a)') '------------------------------------------------------------'
        flush(iout)

      contains

        integer function swapc(i0)                      ! (p0, q R) -> (q0, p -R)
          integer, intent(in) :: i0
          integer :: p0, q0
          p0 = mod(i0-1,nfuc) + 1 ; q0 = (i0-1)/nfuc + 1
          swapc = (trT(p0, fc(q0))-1)*nfuc + fl(q0)
        end function swapc

        ! shell-pair key, normalized so that the first shell is in cell 0
        integer function pkey(xa, xb)
          integer, intent(in) :: xa, xb
          pkey = (sh_of(trT(xb, fc(xa)))-1)*d_nshuc + sh_of(trT(xa, fc(xa)))
        end function pkey

        logical function sp_is_rep(S1, S2)
          integer, intent(in) :: S1, S2
          integer :: g0, k0, a0, b0
          k0 = (S2-1)*d_nshuc + S1
          sp_is_rep = .false.
          do g0 = 1, nG
            a0 = sym_perm(d_sh_first(S1), g0) ; b0 = sym_perm(d_sh_first(S2), g0)
            if (pkey(a0,b0) < k0 .or. pkey(b0,a0) < k0) return
          end do
          sp_is_rep = .true.
        end function sp_is_rep

        logical function sp_rep_of(i0)                 ! candidate belongs to a computed (rep) shell pair
          integer, intent(in) :: i0
          sp_rep_of = known(i0)
        end function sp_rep_of

        subroutine grow2(A)
          double precision, allocatable, intent(inout) :: A(:,:)
          double precision, allocatable :: tt(:,:)
          allocate(tt(size(A,1), 2*size(A,2)))
          tt(:,1:size(A,2)) = A
          call move_alloc(tt, A)
        end subroutine grow2

        ! c(ii) = sum_R (ii_-R | J_0) for every component J of the shell pair of J0
        subroutine pair_cvecs(J0)
          integer, intent(in) :: J0
          integer :: S_C, S_D, n_C, n_D, c_, d_, comp, Jc, i_, p_, q_, r_
          call pair_cols(J0)
          S_C = sh_of(mod(J0-1,nfuc)+1) ; S_D = sh_of((J0-1)/nfuc+1)
          n_C = d_sh_n(S_C) ; n_D = d_sh_n(S_D)
          do d_ = 1, n_D
            do c_ = 1, n_C
              comp = (d_-1)*n_C + c_
              Jc = (d_sh_first(S_D)+d_-2)*nfuc + d_sh_first(S_C)+c_-1
              if (cslot(Jc) /= 0) cycle
              nsl = nsl + 1
              if (nsl > size(cst,2)) call grow2(cst)
              cslot(Jc) = nsl
              !$omp parallel do default(shared) private(i_,p_,q_,r_) schedule(static)
              do i_ = 1, ncand
                p_ = mod(i_-1,nfuc) + 1 ; q_ = (i_-1)/nfuc + 1
                cst(i_, nsl) = 0.d0
                do r_ = 1, ncell
                  cst(i_, nsl) = cst(i_, nsl) + colw(trT(p_,r_), trT(q_,r_), comp)
                end do
              end do
              !$omp end parallel do
            end do
          end do
          deallocate(colw) ; allocate(colw(1,1,1))
        end subroutine pair_cvecs

        ! colw(a, b, comp) = (a b | c d) , a, b in the supercell, (c,d) the components of the
        ! shell pair of J0 (c in cell 0); quartets moved so that their first shell is in cell 0
        subroutine pair_cols(J0)
          integer, intent(in) :: J0
          integer :: S_C, S_D, n_C, n_D, S1, S2, x_, y_, c_, d_, r1_, T1, T2, TC, TD, fx, fy, comp
          double precision, allocatable :: bk(:,:,:,:)
          S_C = sh_of(mod(J0-1,nfuc)+1) ; S_D = sh_of((J0-1)/nfuc+1)
          n_C = d_sh_n(S_C) ; n_D = d_sh_n(S_D)
          if (allocated(colw)) deallocate(colw)
          allocate(colw(nBas, nBas, n_C*n_D)) ; colw = 0.d0
          ps_on = c_prim_screen
          !$omp parallel do default(shared) private(S1,S2,bk,x_,y_,c_,d_,r1_,T1,T2,TC,TD,fx,fy,comp) &
          !$omp schedule(dynamic,1)
          do S1 = 1, d_nsh
            do S2 = S1, d_nsh
              if (d_Qsh(S1,S2)*d_Qsh(S_C,S_D) < 1.d-15) cycle
              allocate(bk(d_sh_n(S1),d_sh_n(S2),n_C,n_D))
              r1_ = fc(d_sh_first(S1))
              T1 = sh_of(trT(d_sh_first(S1), r1_)) ; T2 = sh_of(trT(d_sh_first(S2), r1_))
              TC = sh_of(trT(d_sh_first(S_C), r1_)) ; TD = sh_of(trT(d_sh_first(S_D), r1_))
              if (ps_on) call set_ps(T1, T2, TC, TD)
              call ERI_shell_block_toroidal(ndp, d_sh_n(S1), AO(d_sh_first(T1):), d_sh_n(S2), AO(d_sh_first(T2):), &
                                                 n_C, AO(d_sh_first(TC):), n_D, AO(d_sh_first(TD):), bk)
              do d_ = 1, n_D
                do c_ = 1, n_C
                  comp = (d_-1)*n_C + c_
                  do y_ = 1, d_sh_n(S2)
                    fy = d_sh_first(S2)+y_-1
                    do x_ = 1, d_sh_n(S1)
                      fx = d_sh_first(S1)+x_-1
                      colw(fx, fy, comp) = bk(x_,y_,c_,d_)
                      colw(fy, fx, comp) = bk(x_,y_,c_,d_)
                    end do
                  end do
                end do
              end do
              deallocate(bk)
            end do
          end do
          !$omp end parallel do
          ps_on = .false.
        end subroutine pair_cols

      end subroutine df_init_sym

      !-----------------------------------------------------------------!
      ! irreducible k points: dfs_kirr(1:nirr) representatives (Gamma   !
      ! first), and for every k: k = (conj ? -1 : 1) * R_g k_rep         !
      !-----------------------------------------------------------------!
      subroutine df_kirr(ncell, nG, nirr)
        use symmetry_clifford, only: sym_R
        integer, intent(in)  :: ncell, nG
        integer, intent(out) :: nirr
        integer :: k, g, kr, kneg
        character(len=8) :: envv
        if (allocated(dfs_kirr)) deallocate(dfs_kirr, dfs_krep, dfs_kg, dfs_kcj, dfs_kslot, dfs_krot)
        allocate(dfs_kirr(ncell), dfs_krep(ncell), dfs_kg(ncell), dfs_kcj(ncell), dfs_kslot(ncell), dfs_krot(ncell,nG))
        do g = 1, nG
          do k = 1, ncell
            dfs_krot(k,g) = cell_of(matmul(sym_R(:,:,g), df_coord(:,k)))
          end do
        end do
        dfs_krep = 0 ; dfs_kslot = 0 ; nirr = 0
        ! test switch: CIPOP_DFSYM_ALLK=1 stores every k2 (no reconstruction by symmetry)
        call get_environment_variable('CIPOP_DFSYM_ALLK', envv)
        if (trim(envv) == '1') then
          do k = 1, ncell
            nirr = nirr + 1 ; dfs_kirr(nirr) = k ; dfs_kslot(k) = nirr
            dfs_krep(k) = k ; dfs_kg(k) = 1 ; dfs_kcj(k) = 0
          end do
          return
        end if
        do k = 1, ncell
          if (dfs_krep(k) /= 0) cycle
          nirr = nirr + 1 ; dfs_kirr(nirr) = k ; dfs_kslot(k) = nirr
          do g = 1, nG                                   ! the star of k (+ time reversal)
            kr = dfs_krot(k,g)
            if (dfs_krep(kr) == 0) then
              dfs_krep(kr) = k ; dfs_kg(kr) = g ; dfs_kcj(kr) = 0
            end if
          end do
          do g = 1, nG
            kneg = csub(1, dfs_krot(k,g))
            if (dfs_krep(kneg) == 0) then
              dfs_krep(kneg) = k ; dfs_kg(kneg) = g ; dfs_kcj(kneg) = 1
            end if
          end do
        end do
      end subroutine df_kirr

      !-----------------------------------------------------------------!
      ! slice S(a,q,b,k1) = Bh(a,q,b;k1,k2) for any k2, from the stored  !
      ! irreducible slice (see df_init_sym)                              !
      !-----------------------------------------------------------------!
      subroutine dfs_slice(nfuc, nau, ncell, k2, S)
        integer   , intent(in)  :: nfuc, nau, ncell, k2
        complex(8), intent(out) :: S(nfuc, nau, nfuc, ncell)
        complex(8), allocatable :: T(:,:,:,:)
        integer :: kr, g, a, b, q, k1, k1r, pa, pb, pq, ca, cb, cq, ktgt
        double precision :: sgn
        complex(8) :: ph2a
        kr = dfs_krep(k2) ; g = dfs_kg(k2)
        if (kr == k2 .and. dfs_kslot(k2) /= 0) then
          call copy_slot(dfs_kslot(k2))
          return
        end if
        ! target T = Bh(.; R k1', R kr) from the stored slice of kr
        allocate(T(nfuc, nau, nfuc, ncell))
        ktgt = dfs_krot(kr, g)
        !$omp parallel do default(shared) private(k1,k1r,a,b,q,pa,pb,pq,ca,cb,cq,sgn,ph2a) schedule(static)
        do k1 = 1, ncell
          k1r = dfs_krot(k1, g)
          do b = 1, nfuc
            pb = dfs_aop(b,g) ; cb = dfs_aoc(b,g)
            do q = 1, nau
              pq = dfs_aqp(q,g) ; cq = dfs_aqc(q,g)
              do a = 1, nfuc
                pa = dfs_aop(a,g) ; ca = dfs_aoc(a,g)
                sgn = dble(dfs_aos(a,g)*dfs_aos(b,g)*dfs_aqs(q,g))
                ph2a = df_phase(k1r, csub(cb,ca))*df_phase(ktgt, csub(cq,ca))
                T(pa,pq,pb,k1r) = sgn*ph2a*bigval(a,q,b,k1)
              end do
            end do
          end do
        end do
        !$omp end parallel do
        if (dfs_kcj(k2) == 0) then
          S = T
        else                                               ! Bh(k1,k2) = conj Bh(-k1,-k2), -k2 = ktgt
          do k1 = 1, ncell
            S(:,:,:,k1) = conjg(T(:,:,:,csub(1,k1)))
          end do
        end if
        deallocate(T)
      contains
        complex(8) function bigval(a_, q_, b_, k_)
          integer, intent(in) :: a_, q_, b_, k_
          integer(8) :: idx
          idx = int(a_,8) + int(nfuc,8)*(int(q_-1,8) + int(nau,8)*(int(b_-1,8) + int(nfuc,8)*int(k_-1,8)))
          bigval = dfs_big(idx, dfs_kslot(kr))
        end function bigval
        subroutine copy_slot(sl)
          integer, intent(in) :: sl
          integer(8) :: n8
          n8 = int(nfuc,8)*int(nau,8)*int(nfuc,8)*int(ncell,8)
          call zcopy_big(n8, dfs_big(1,sl), S)
        end subroutine copy_slot
      end subroutine dfs_slice

      subroutine zcopy_big(n8, src, dst)
        integer(8), intent(in) :: n8
        complex(8), intent(in) :: src(n8)
        complex(8), intent(out) :: dst(n8)
        dst = src
      end subroutine zcopy_big

      !-----------------------------------------------------------------!
      ! J and K with the symmetric fitting (rows mu <= nfuc)             !
      !-----------------------------------------------------------------!
      subroutine df_JK_sym(nBas, P, J, K)
        use unitcell_module, only: nfuc
        integer         , intent(in)  :: nBas
        double precision, intent(in)  :: P(nBas,nBas)
        double precision, intent(out) :: J(nBas,nBas), K(nBas,nBas)
        integer :: nau, ncell, nb, r, k1, k2, kk, ik, no, nmax, n, i, info
        double precision, allocatable :: gam(:), Jv(:), wk(:,:), ew(:), rw(:)
        complex(8), allocatable :: pk(:,:), ck(:,:,:), Kt(:,:,:), S(:,:,:,:), L(:,:,:), Lw(:,:,:), zw(:)
        integer, allocatable :: nok(:)

        nau = dfs_npiv ; ncell = df_ncell ; nb = nfuc*nBas
        J = 0.d0 ; K = 0.d0

        ! ---- J : gamma_Q = sum_{l0,si} B0(l0,si,Q) P(l0,si) ; J(mu0,nu) = sum_Q B0 gamma ---- !
        allocate(gam(nau), Jv(nb))
        call dgemv('T', nb, nau, 1.d0, dfs_B0, nb, P(1:nfuc,1:nBas), 1, 0.d0, gam, 1)
        call dgemv('N', nb, nau, 1.d0, dfs_B0, nb, gam, 1, 0.d0, Jv, 1)
        J(1:nfuc,1:nBas) = reshape(Jv, (/ nfuc, nBas /))
        deallocate(gam, Jv)

        ! ---- K : Kt(K) = -1/(2N) sum_{k1} M P(k1) M^H , M = Bh(k1, K-k1), one k2 slice at a time ---- !
        allocate(pk(nfuc,nfuc), ck(nfuc,nfuc,ncell), wk(nfuc,ncell), nok(ncell), ew(nfuc), &
                 rw(3*nfuc), zw(4*nfuc), Kt(nfuc,nfuc,ncell))
        do ik = 1, ncell
          pk = (0.d0,0.d0)
          do r = 1, ncell
            pk = pk + P(1:nfuc, (r-1)*nfuc+1:r*nfuc)*df_phase(ik,r)
          end do
          pk = 0.5d0*(pk + conjg(transpose(pk)))
          call zheev('V','U', nfuc, pk, nfuc, ew, zw, 4*nfuc, rw, info)
          if (info /= 0) error stop 'df_JK_sym: zheev failed for P(k)'
          n = 0
          do i = 1, nfuc
            if (dabs(ew(i)) > 1.d-10) then
              n = n + 1 ; ck(:,n,ik) = pk(:,i) ; wk(n,ik) = ew(i)
            end if
          end do
          nok(ik) = n
        end do
        nmax = max(1, maxval(nok))
        Kt = (0.d0,0.d0)
        allocate(S(nfuc, nau, nfuc, ncell))
        do k2 = 1, ncell
          call dfs_slice(nfuc, nau, ncell, k2, S)
          nb_saved = blas_serial()
          !$omp parallel default(shared) private(k1,kk,no,L,Lw,i)
          allocate(L(nfuc,nau,nmax), Lw(nfuc,nau,nmax))
          !$omp do schedule(dynamic,1)
          do k1 = 1, ncell
            no = nok(k1) ; if (no == 0) cycle
            kk = cadd(k1, k2)
            call zgemm('N','N', nfuc*nau, no, nfuc, (1.d0,0.d0), S(1,1,1,k1), nfuc*nau, ck(1,1,k1), nfuc, &
                       (0.d0,0.d0), L, nfuc*nau)
            do i = 1, no
              Lw(:,:,i) = L(:,:,i)*wk(i,k1)
            end do
            call zgemm('N','C', nfuc, nfuc, nau*no, cmplx(-0.5d0/dble(ncell),0.d0,kind=8), Lw, nfuc, L, nfuc, &
                       (1.d0,0.d0), Kt(1,1,kk), nfuc)
          end do
          !$omp end do
          deallocate(L, Lw)
          !$omp end parallel
          call blas_restore(nb_saved)
        end do
        deallocate(S)
        do r = 1, ncell
          pk = (0.d0,0.d0)
          do kk = 1, ncell
            pk = pk + Kt(:,:,kk)*conjg(df_phase(kk,r))
          end do
          K(1:nfuc, (r-1)*nfuc+1:r*nfuc) = dble(pk)/dble(ncell)
        end do
        deallocate(pk, ck, wk, nok, ew, rw, zw, Kt)
      end subroutine df_JK_sym

end module df_toroidal
