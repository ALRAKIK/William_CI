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
      implicit none

      private
      public :: df_init, df_JK, df_kMO, df_ready, df_build_chol

      logical, save :: df_ready = .false.
      integer, save :: df_nau = 0, df_naux = 0, df_ncell = 0
      integer, allocatable, save :: df_coord(:,:)
      double precision, allocatable, save :: df_B(:,:,:)
      ! [Nov 2026] k-space B for the exchange (HF only): k2 over a half set df_kh, conj for the rest
      complex(8), allocatable, save :: df_Bh(:,:,:,:,:)
      integer, save :: df_nkh = 0
      integer, allocatable, save :: df_kh(:), df_khi(:)
      complex(8), allocatable, save :: df_phase(:,:)
      double precision, parameter :: lindep_thresh = 1.0d-10
      double precision, parameter :: df_schwarz    = 1.0d-13
      double precision, parameter :: etb_beta      = 2.0d0     ! PySCF ETB_BETA
      ! inputs kept for a rebuild with the MP2-fit basis (file AuxBasis_MP2)
      integer, save :: sv_ndp = 0, sv_iout = 6
      type(ERI_function), allocatable, save :: sv_AO(:)
      double precision, allocatable, save :: sv_geom(:,:)
      integer, allocatable, save :: sv_charge(:)

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

        fname = 'AuxBasis'
        if (present(fname_in)) then
          fname = fname_in
        else
          sv_ndp = ndp ; sv_iout = iout
          if (allocated(sv_AO)) deallocate(sv_AO, sv_geom, sv_charge)
          allocate(sv_AO(nBas), sv_geom(n_atoms,3), sv_charge(n_atoms))
          sv_AO = AO ; sv_geom = geometry ; sv_charge = atoms(:)%charge
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
        allocate(aat(n_atoms))
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
        df_nau = nau ; df_naux = naux

        ! t-tables large enough for the auxiliary exponents and l
        lmax = 0 ; emax = 0.d0
        do i = 1, nBas
          call orbital_to_lxyz(AO(i)%orbital, lxyz)
          lmax = max(lmax, maxval(lxyz(1:ndp))) ; emax = max(emax, 2.d0*maxval(AO(i)%exponent))
        end do
        do i = 1, naux
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
        if (.not. use_file) call print_etb(n_atoms, aat, iout)
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
        allocate(Qa(naux))
        do r = 1, ncell
          do p = 1, nau
            Qa((r-1)*nau+p) = dsqrt(dabs(Vm(p,p)))
          end do
        end do

        ! V^-1/2 per k, back to real space: mf(P,Q,R) = [V^-1/2](P0, Q R)
        allocate(fk(nau,nau,ncell), Mk(nau,nau), w(nau), rwork(3*nau), work(1))
        ndrop_tot = 0 ; wmin = huge(1.d0)
        call zheev('V','U', nau, Mk, nau, w, work, -1, rwork, info)
        lwork = max(1, int(dble(work(1)))) ; deallocate(work) ; allocate(work(lwork))
        do k = 1, ncell
          Mk = (0.d0,0.d0)
          do r = 1, ncell
            Mk = Mk + Vm(:, (r-1)*nau+1:r*nau)*phase(k,r)
          end do
          Mk = 0.5d0*(Mk + conjg(transpose(Mk)))
          call zheev('V','U', nau, Mk, nau, w, work, lwork, rwork, info)
          if (info /= 0) error stop 'keyword density_fitting: zheev failed for the metric'
          wmin = min(wmin, w(1))
          ndrop = 0
          do p = 1, nau
            if (w(p) > lindep_thresh) then
              Mk(:,p) = Mk(:,p)/dsqrt(dsqrt(w(p)))      ! U diag(w^-1/4): f = (U w^-1/4)(U w^-1/4)^H
            else
              Mk(:,p) = (0.d0,0.d0) ; ndrop = ndrop + 1
            end if
          end do
          ndrop_tot = ndrop_tot + ndrop
          fk(:,:,k) = matmul(Mk, conjg(transpose(Mk)))       ! [V(k)]^-1/2
        end do
        deallocate(Mk, w, rwork, work, Vm)
        t_V = omp_get_wtime() - w0
        write(iout,'(a,es12.4)') ' smallest eigenvalue of (P|Q)(k)  : ', wmin
        write(iout,'(a,i10)')    ' dropped (eigenvalue < 1e-10)     : ', ndrop_tot
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
        write(iout,'(a)') '------------------------------------------------------------'
        flush(iout)
        df_ready = .true.
      end subroutine df_init


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

end module df_toroidal
