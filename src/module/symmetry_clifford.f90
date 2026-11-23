!=====================================================================!
!                                                                     !
!   module symmetry_clifford   [Nov 2026]   keyword symmetry          !
!                                                                     !
!   Step 1 of the geometric symmetry: DETECTION and CHECK only.       !
!   Nothing in the calculation uses it yet (no time saved yet).       !
!                                                                     !
!   Which operations can be used with Clifford (toroidal) Gaussians   !
!   ---------------------------------------------------------------   !
!   A Clifford Gaussian is a product over the Cartesian directions:   !
!     phi(r) = prod_i X_i(r_i - A_i)^l_i  exp(-alpha sum_i G_i)       !
!     X_i(u) = sin(a_i u)/a_i , G_i = 2(1 - cos(a_i u))/a_i^2 (period)!
!     X_i(u) = u              , G_i = u^2                (no period)  !
!   X is odd and G is even in u, so an operation                     !
!     g : r -> R r + t ,  R = signed permutation of x, y, z          !
!   ((R r)_i = s_i r_p(i), s_i = +-1) maps a Clifford Gaussian on a   !
!   Clifford Gaussian, exactly:                                      !
!     phi_{A,l}(g^-1 r) = prod_i s_i^{l_p(i)}  phi_{RA+t, l'}(r),     !
!     l'_i = l_p(i)                                                  !
!   provided a periodic axis is only exchanged with a periodic axis  !
!   of the SAME length (a_i = a_p(i)), and a non-periodic axis with a !
!   non-periodic one. Other crystal operations (3- and 6-fold axes of !
!   hexagonal cells, mirrors not perpendicular to an axis) are NOT    !
!   symmetries of the toroidal integrals and are not searched.        !
!   The 48 signed permutations contain: E, i, mirrors x=0 / x=y,       !
!   C2, C4 about the axes, C2 about face diagonals, C3 about the body !
!   diagonals, S4, S6: the cubic group Oh and all its subgroups.      !
!                                                                     !
!   Detection: every signed permutation R allowed by the cell, every  !
!   translation t that sends atom 1 on a unit-cell atom of the same   !
!   nuclear charge; kept if all unit-cell atoms land on atoms of the  !
!   same charge (periodic components modulo L, others exactly). t is  !
!   reduced modulo the cell length d in the periodic directions, so   !
!   the list is the space group modulo the lattice translations (the  !
!   N_cell translations are already used by SAO / k-points).          !
!                                                                     !
!   AO representation: every operation is a SIGNED PERMUTATION of the !
!   Cartesian AOs (same atom shell, permuted lxyz, sign above).       !
!   Check: max |s_m s_n M(pm,pn) - M(m,n)| for S, T, V and the same   !
!   for the ERIs (all, or random quartets: ~2e7 checks in total).    !
!                                                                     !
!=====================================================================!

module symmetry_clifford

      implicit none
      private

      public :: sym_run, sym_setup, sym_nop, sym_R, sym_t, sym_perm, sym_sign, sym_ok, sym_pg_name
      public :: sym_ready, sym_eri_note
      public :: sym_pairs_on, sym_pair_is_rep, sym_pair_fill   ! [Nov 2026] one-electron matrices

      integer         , save              :: sym_nop = 0       ! space-group operations mod lattice translations
      integer         , save, allocatable :: sym_R(:,:,:)      ! (3,3,nop) signed permutation matrices
      double precision, save, allocatable :: sym_t(:,:)        ! (3,nop) translation (bohr)
      integer         , save, allocatable :: sym_perm(:,:)     ! (nBas,nop) AO image
      integer         , save, allocatable :: sym_sign(:,:)     ! (nBas,nop) +-1
      logical         , save, allocatable :: sym_ok(:)         ! passed the numerical check
      character(len=8), save              :: sym_pg_name = 'C1'
      logical         , save              :: sym_ready = .false.  ! sym_setup done (ops + AO maps)
      integer         , save              :: sym_nbas  = 0
      character(len=200), save            :: sym_eri_note = ''      ! ERI driver report (printed by sym_run)
      ! [Nov 2026] AO pairs (i in the unit cell, j): representative of the orbit under the
      ! point operations x translations x (i<->j). sym_prep(i,j) = (ir-1)*nBas + jr of the
      ! representative (ir <= nfuc), sym_psg = sign s_i s_j of the operation.
      logical         , save              :: sym_pairs_on = .false.
      integer(8)      , save, allocatable :: sym_prep(:,:)
      integer(1)      , save, allocatable :: sym_psg(:,:)
      integer         , save              :: sym_pnfuc = 0

      double precision, parameter :: pos_tol = 1.d-6           ! bohr
      double precision, parameter :: chk_tol = 1.d-8           ! Hartree (integral check)

      ! geometry of the run
      integer :: ndp, ncel(3), nat, nauc
      double precision :: Lp(3), dcell(3)
      logical :: per(3)

contains

      !---------------------------------------------------------------!
      ! driver: detect, build AO maps, check S, T, V, ERI, print      !
      !---------------------------------------------------------------!
      subroutine sym_run(calc_type, n_atoms, geometry, charge, nBas, nfuc_in, AO, S, T, V, &
                         ERI, n1, have_eri, iout)
        use classification_ERI, only: ERI_function
        use unitcell_module   , only: nx, ny, nz
        use omp_lib           , only: omp_get_wtime
        character(len=*)  , intent(in) :: calc_type
        integer           , intent(in) :: n_atoms, nBas, nfuc_in, n1, iout
        double precision  , intent(in) :: geometry(n_atoms,3)
        integer           , intent(in) :: charge(n_atoms)
        type(ERI_function), intent(in) :: AO(nBas)
        double precision  , intent(in) :: S(nBas,nBas), T(nBas,nBas), V(nBas,nBas)
        double precision  , intent(in) :: ERI(n1,nBas,nBas,nBas)
        logical           , intent(in) :: have_eri

        integer :: k, nbad
        double precision :: eS, eT, eV, eE
        double precision :: w0, w1
        integer :: ndp_in

        w0 = omp_get_wtime()
        select case (trim(calc_type))
        case ('Tori1D') ; ndp_in = 1
        case ('Tori2D') ; ndp_in = 2
        case ('Tori3D') ; ndp_in = 3
        case default
          write(iout,'(a)') ' keyword symmetry: only for Tori1D / Tori2D / Tori3D, ignored'
          return
        end select
        ! the operations are normally found before the ERIs (sym_setup in Tori1D/2D/3D);
        ! the AO maps only depend on the AO order, so the unnormalized list gives the same maps
        if (.not. sym_ready .or. sym_nbas /= nBas) call sym_setup(ndp_in, n_atoms, geometry, charge, nBas, AO)

        write(iout,*)
        write(iout,'(a)') ' ============================================================================='
        write(iout,'(a)') '   Geometric symmetry (keyword symmetry)   [detection + check only]'
        write(iout,'(a)') ' ============================================================================='
        write(iout,'(a,i0,a,i0,a)') '   periodic directions: ', ndp, '    atoms per cell: ', nauc
        write(iout,'(a)') '   usable operations: signed permutations of x, y, z (+ translation)'
        write(iout,'(a)') '   check: max |s_m s_n M(Pm,Pn) - M(m,n)| over the AO pairs / quartets'
        write(iout,*)
        write(iout,'(a)') '    #  operation                      type   translation t (bohr)                 '// &
                          '   dS        dT        dV        dERI     ok'
        nbad = 0
        do k = 1, sym_nop
          eS = mat_dev(nBas, S, k) ; eT = mat_dev(nBas, T, k) ; eV = mat_dev(nBas, V, k)
          eE = -1.d0
          if (have_eri) eE = eri_dev(nBas, nfuc_in, n1, ERI, k)
          sym_ok(k) = max(eS, eT, eV, eE) <= chk_tol .and. all(sym_perm(:,k) > 0)
          if (.not. sym_ok(k)) nbad = nbad + 1
          write(iout,'(i5,2x,a28,2x,a5,2x,3f12.6,2x,4es10.2,3x,l1)') k, op_string(sym_R(:,:,k)), &
               op_type(sym_R(:,:,k), sym_t(:,k)), sym_t(:,k), eS, eT, eV, eE, sym_ok(k)
        end do
        sym_pg_name = point_group_name()
        write(iout,*)
        write(iout,'(a,a,a,i0,a)') '   Point group (Clifford-usable part): ', trim(sym_pg_name), &
                                   '   order ', n_point_ops(), ' (distinct R)'
        write(iout,'(a,i0,a,i0)') '   space-group operations mod lattice translations: ', sym_nop, &
                                  '   x lattice translations: ', nx*ny*nz
        if (.not. have_eri) write(iout,'(a)') '   (ERIs not stored: ERI check skipped, dERI = -1)'
        if (len_trim(sym_eri_note) > 0) write(iout,'(a)') '   ERIs: '//trim(sym_eri_note)
        if (nbad == 0) then
          write(iout,'(a,es8.1,a)') '   symmetry check: all operations OK (tolerance ', chk_tol, ' Eh)'
        else
          write(iout,'(a,i0,a)') '   symmetry check: ', nbad, ' operation(s) FAILED the integral check'
          write(iout,'(a)') '   (geometry symmetric only within the position tolerance, or a bug: see dS..dERI)'
        end if
        w1 = omp_get_wtime()
        write(iout,'(a,f10.2,a)') '   wall time for symmetry detection + check: ', w1 - w0, ' s'
        write(iout,'(a)') ' ============================================================================='
        write(iout,*)
        flush(iout)
      end subroutine sym_run

      !---------------------------------------------------------------!
      ! detection of the operations + AO maps (no integrals needed).  !
      ! Called in Tori1D/2D/3D before the ERIs, so that the ERI driver!
      ! can use the operations (keyword symmetry).                    !
      !---------------------------------------------------------------!
      subroutine sym_setup(ndp_in, n_atoms, geometry, charge, nBas, AO)
        use classification_ERI, only: ERI_function
        use torus_init        , only: Lx, Ly, Lz
        use unitcell_module   , only: nx, ny, nz
        integer           , intent(in) :: ndp_in, n_atoms, nBas
        double precision  , intent(in) :: geometry(n_atoms,3)
        integer           , intent(in) :: charge(n_atoms)
        type(ERI_function), intent(in) :: AO(nBas)
        integer :: k
        ndp  = ndp_in
        per  = [ (k <= ndp, k = 1, 3) ]
        Lp   = [Lx, Ly, Lz]
        ncel = [nx, ny, nz]
        dcell = 0.d0
        do k = 1, ndp
          dcell(k) = Lp(k)/dble(ncel(k))
        end do
        nat  = n_atoms
        nauc = n_atoms / (nx*ny*nz)
        call detect(geometry, charge)
        call build_ao_maps(geometry, nBas, AO)
        sym_ready = .true.
        sym_nbas  = nBas
        call build_pairs(nBas)
      end subroutine sym_setup

      !---------------------------------------------------------------!
      ! [Nov 2026] representatives of the AO pairs (one-electron      !
      ! matrices S, T, V: only the representatives are computed, the  !
      ! others are copies with the sign s_i s_j, see sym_pair_fill)   !
      !---------------------------------------------------------------!
      subroutine build_pairs(nBas)
        use unitcell_module, only: nfuc
        integer, intent(in) :: nBas
        integer :: i, j, k, a, b, c(3), ia, jb
        integer(8) :: key, best
        integer :: sg
        sym_pairs_on = .false.
        if (allocated(sym_prep)) deallocate(sym_prep, sym_psg)
        if (sym_nop <= 1 .or. any(sym_perm <= 0)) return
        if (nBas /= nfuc*ncel(1)*ncel(2)*ncel(3)) return
        allocate(sym_prep(nfuc,nBas), sym_psg(nfuc,nBas))
        !$omp parallel do default(shared) private(i,j,k,a,b,c,ia,jb,key,best,sg) schedule(static)
        do j = 1, nBas
          do i = 1, nfuc
            best = huge(1_8) ; sg = 1
            do k = 1, sym_nop
              a = sym_perm(i,k) ; b = sym_perm(j,k)
              ! (a, b) with a brought to cell 0
              call cell_of(a, nfuc, c)
              ia = shift_ao(a, nfuc, c) ; jb = shift_ao(b, nfuc, c)
              key = int(ia-1,8)*int(nBas,8) + int(jb,8)
              if (key < best) then
                best = key ; sg = sym_sign(i,k)*sym_sign(j,k)
              end if
              ! (b, a): the matrices are symmetric
              call cell_of(b, nfuc, c)
              ia = shift_ao(b, nfuc, c) ; jb = shift_ao(a, nfuc, c)
              key = int(ia-1,8)*int(nBas,8) + int(jb,8)
              if (key < best) then
                best = key ; sg = sym_sign(i,k)*sym_sign(j,k)
              end if
            end do
            sym_prep(i,j) = best ; sym_psg(i,j) = int(sg,1)
          end do
        end do
        !$omp end parallel do
        sym_pnfuc = nfuc
        sym_pairs_on = .true.
      end subroutine build_pairs

      logical function sym_pair_is_rep(i, j, nBas)
        integer, intent(in) :: i, j, nBas
        sym_pair_is_rep = (sym_prep(i,j) == int(i-1,8)*int(nBas,8) + int(j,8))
      end function sym_pair_is_rep

      ! M(i,j), i <= nfuc: copy of the representative with the sign
      subroutine sym_pair_fill(nBas, M)
        integer         , intent(in)    :: nBas
        double precision, intent(inout) :: M(nBas,nBas)
        integer :: i, j, ir, jr
        integer(8) :: key
        !$omp parallel do default(shared) private(i,j,ir,jr,key) schedule(static)
        do j = 1, nBas
          do i = 1, sym_pnfuc
            key = sym_prep(i,j)
            ir = int((key-1)/int(nBas,8)) + 1 ; jr = int(mod(key-1, int(nBas,8))) + 1
            if (ir /= i .or. jr /= j) M(i,j) = dble(sym_psg(i,j))*M(ir,jr)
          end do
        end do
        !$omp end parallel do
      end subroutine sym_pair_fill

      !---------------------------------------------------------------!
      ! difference of two positions in direction i (mod L if periodic)!
      !---------------------------------------------------------------!
      double precision function pdiff(x, i, period)
        double precision, intent(in) :: x, period
        integer, intent(in) :: i
        if (per(i)) then
          pdiff = x - period*dnint(x/period)
        else
          pdiff = x
        end if
      end function pdiff

      !---------------------------------------------------------------!
      ! the 48 signed permutations allowed by the cell                !
      !---------------------------------------------------------------!
      subroutine candidate_R(nR, Rl)
        integer, intent(out) :: nR, Rl(3,3,48)
        integer :: perms(3,6), ip, is, i, s(3), p(3)
        logical :: good
        perms = reshape([1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1], [3,6])
        nR = 0
        do ip = 1, 6
          p = perms(:,ip)
          good = .true.
          do i = 1, 3
            if (per(i) .neqv. per(p(i))) good = .false.
            if (per(i) .and. per(p(i))) then
              if (abs(Lp(i) - Lp(p(i))) > 1.d-10*Lp(i) .or. ncel(i) /= ncel(p(i))) good = .false.
            end if
          end do
          if (.not. good) cycle
          do is = 0, 7
            s = [1 - 2*mod(is,2), 1 - 2*mod(is/2,2), 1 - 2*mod(is/4,2)]
            nR = nR + 1
            Rl(:,:,nR) = 0
            do i = 1, 3
              Rl(i, p(i), nR) = s(i)
            end do
          end do
        end do
      end subroutine candidate_R

      !---------------------------------------------------------------!
      ! image of atom a under (R,t): atom index or 0                  !
      !---------------------------------------------------------------!
      integer function image_atom(a, R, t, geometry, charge)
        integer, intent(in) :: a, R(3,3), charge(nat)
        double precision, intent(in) :: t(3), geometry(nat,3)
        double precision :: y(3), dv(3)
        integer :: l, i, c(3), b
        y = matmul(dble(R), geometry(a,:)) + t
        image_atom = 0
        do l = 1, nauc                                   ! unit-cell atom l + lattice vector?
          if (charge(l) /= charge(a)) cycle
          dv = y - geometry(l,:)
          c = 0
          do i = 1, 3
            if (per(i)) then
              c(i) = nint(dv(i)/dcell(i))
              if (abs(dv(i) - c(i)*dcell(i)) > pos_tol) goto 10
              c(i) = modulo(c(i), ncel(i))
            else
              if (abs(dv(i)) > pos_tol) goto 10
            end if
          end do
          b = ((c(3)*ncel(2) + c(2))*ncel(1) + c(1))*nauc + l
          if (b >= 1 .and. b <= nat) then
            image_atom = b
            return
          end if
10        continue
        end do
      end function image_atom

      !---------------------------------------------------------------!
      ! detection of the operations (R,t), t mod d                    !
      !---------------------------------------------------------------!
      subroutine detect(geometry, charge)
        double precision, intent(in) :: geometry(nat,3)
        integer         , intent(in) :: charge(nat)
        integer :: nR, Rl(3,3,48), ir, j, a, i, k, nmax
        double precision :: t(3)
        integer, allocatable :: Rt(:,:,:)
        double precision, allocatable :: tt(:,:)
        logical :: good, dup

        call candidate_R(nR, Rl)
        nmax = nR*max(nauc,1)
        allocate(Rt(3,3,nmax), tt(3,nmax))
        sym_nop = 0
        do ir = 1, nR
          do j = 1, nauc
            if (charge(j) /= charge(1)) cycle
            t = geometry(j,:) - matmul(dble(Rl(:,:,ir)), geometry(1,:))
            do i = 1, 3                                    ! reduce t into [0, d)
              if (per(i)) then
                t(i) = modulo(t(i), dcell(i))
                if (dcell(i) - t(i) < pos_tol) t(i) = 0.d0
              end if
            end do
            good = .true.
            do a = 1, nauc
              if (image_atom(a, Rl(:,:,ir), t, geometry, charge) == 0) then
                good = .false. ; exit
              end if
            end do
            if (.not. good) cycle
            dup = .false.
            do k = 1, sym_nop
              if (all(Rt(:,:,k) == Rl(:,:,ir))) then
                if (all([(abs(pdiff(t(i)-tt(i,k), i, dcell(i))) < pos_tol, i = 1, 3)])) dup = .true.
              end if
            end do
            if (dup) cycle
            sym_nop = sym_nop + 1
            Rt(:,:,sym_nop) = Rl(:,:,ir) ; tt(:,sym_nop) = t
          end do
        end do
        if (allocated(sym_R)) deallocate(sym_R, sym_t, sym_ok)
        allocate(sym_R(3,3,sym_nop), sym_t(3,sym_nop), sym_ok(sym_nop))
        sym_R = Rt(:,:,1:sym_nop) ; sym_t = tt(:,1:sym_nop) ; sym_ok = .false.
        deallocate(Rt, tt)
      end subroutine detect

      !---------------------------------------------------------------!
      ! AO maps: perm(m) = image AO, sign(m) = +-1 (0/-1 if not found)!
      !---------------------------------------------------------------!
      subroutine build_ao_maps(geometry, nBas, AO)
        use classification_ERI, only: ERI_function
        double precision  , intent(in) :: geometry(nat,3)
        integer           , intent(in) :: nBas
        type(ERI_function), intent(in) :: AO(nBas)
        integer, allocatable :: ao_atom(:), first(:), last(:), lx(:,:), occ(:), amap(:)
        integer :: m, a, k, b, n, i, lnew(3), sgn, p(3), s(3), cnt
        integer, allocatable :: chg_dummy(:)
        double precision :: dv(3)

        allocate(ao_atom(nBas), first(nat), last(nat), lx(3,nBas), occ(nBas), amap(nat))
        first = 0 ; last = -1
        do m = 1, nBas                                       ! AO -> atom (AOs are atom-contiguous)
          ao_atom(m) = 0
          do a = 1, nat
            dv = [AO(m)%x, AO(m)%y, AO(m)%z] - geometry(a,:)
            if (maxval(abs(dv)) < 1.d-8) then
              ao_atom(m) = a ; exit
            end if
          end do
          a = ao_atom(m)
          if (a > 0) then
            if (first(a) == 0) first(a) = m
            last(a) = m
          end if
          call lxyz_of(AO(m)%orbital, lx(:,m))
        end do
        do m = 1, nBas                                       ! occurrence of the same lxyz in its atom
          occ(m) = 0
          a = ao_atom(m)
          if (a == 0) cycle
          do n = first(a), m
            if (all(lx(:,n) == lx(:,m))) occ(m) = occ(m) + 1
          end do
        end do

        if (allocated(sym_perm)) deallocate(sym_perm, sym_sign)
        allocate(sym_perm(nBas,sym_nop), sym_sign(nBas,sym_nop))
        allocate(chg_dummy(nat))
        chg_dummy = 0                                        ! charges already checked in detect
        do k = 1, sym_nop
          do a = 1, nat
            amap(a) = image_atom_nocharge(a, sym_R(:,:,k), sym_t(:,k), geometry)
          end do
          p = 1 ; s = 1
          do i = 1, 3
            do b = 1, 3
              if (sym_R(i,b,k) /= 0) then
                p(i) = b ; s(i) = sym_R(i,b,k)
              end if
            end do
          end do
          do m = 1, nBas
            sym_perm(m,k) = 0 ; sym_sign(m,k) = 0
            a = ao_atom(m)
            if (a == 0) cycle
            b = amap(a)
            if (b == 0) cycle
            sgn = 1
            do i = 1, 3
              lnew(i) = lx(p(i), m)
              if (s(i) < 0 .and. mod(lnew(i),2) == 1) sgn = -sgn
            end do
            cnt = 0
            do n = first(b), last(b)
              if (all(lx(:,n) == lnew)) then
                cnt = cnt + 1
                if (cnt == occ(m)) then
                  if (same_shell(AO(m), AO(n))) then
                    sym_perm(m,k) = n ; sym_sign(m,k) = sgn
                  end if
                  exit
                end if
              end if
            end do
          end do
        end do
        deallocate(ao_atom, first, last, lx, occ, amap, chg_dummy)
      end subroutine build_ao_maps

      integer function image_atom_nocharge(a, R, t, geometry)
        integer, intent(in) :: a, R(3,3)
        double precision, intent(in) :: t(3), geometry(nat,3)
        double precision :: y(3), dv(3)
        integer :: l, i, c(3), b
        y = matmul(dble(R), geometry(a,:)) + t
        image_atom_nocharge = 0
        do l = 1, nauc
          dv = y - geometry(l,:)
          c = 0
          do i = 1, 3
            if (per(i)) then
              c(i) = nint(dv(i)/dcell(i))
              if (abs(dv(i) - c(i)*dcell(i)) > pos_tol) goto 20
              c(i) = modulo(c(i), ncel(i))
            else
              if (abs(dv(i)) > pos_tol) goto 20
            end if
          end do
          b = ((c(3)*ncel(2) + c(2))*ncel(1) + c(1))*nauc + l
          if (b >= 1 .and. b <= nat) then
            image_atom_nocharge = b
            return
          end if
20        continue
        end do
      end function image_atom_nocharge

      logical function same_shell(A, B)
        use classification_ERI, only: ERI_function
        type(ERI_function), intent(in) :: A, B
        same_shell = .false.
        if (size(A%exponent) /= size(B%exponent)) return
        if (any(abs(A%exponent - B%exponent) > 1.d-12*abs(A%exponent))) return
        if (any(abs(A%coefficient - B%coefficient) > 1.d-12*max(1.d0,abs(A%coefficient)))) return
        same_shell = .true.
      end function same_shell

      subroutine lxyz_of(label, l)
        character(len=*), intent(in) :: label
        integer, intent(out) :: l(3)
        integer :: c
        l = 0
        do c = 2, len_trim(label)
          select case (label(c:c))
          case ('x','X') ; l(1) = l(1) + 1
          case ('y','Y') ; l(2) = l(2) + 1
          case ('z','Z') ; l(3) = l(3) + 1
          end select
        end do
      end subroutine lxyz_of

      !---------------------------------------------------------------!
      ! AO index shifted by whole cells so that it lies in cell 0      !
      ! (AO index = (cell-1)*nfuc + local, x fastest)                  !
      !---------------------------------------------------------------!
      subroutine cell_of(m, nfuc, c)
        integer, intent(in) :: m, nfuc
        integer, intent(out) :: c(3)
        integer :: lin
        lin = (m-1)/nfuc
        c(1) = mod(lin, ncel(1)) ; lin = lin/ncel(1)
        c(2) = mod(lin, ncel(2)) ; c(3) = lin/ncel(2)
      end subroutine cell_of

      integer function shift_ao(m, nfuc, dc)
        integer, intent(in) :: m, nfuc, dc(3)
        integer :: c(3), loc
        call cell_of(m, nfuc, c)
        loc = mod(m-1, nfuc) + 1
        c = modulo(c - dc, ncel)
        shift_ao = ((c(3)*ncel(2) + c(2))*ncel(1) + c(1))*nfuc + loc
      end function shift_ao

      !---------------------------------------------------------------!
      ! max |s_m s_n M(Pm,Pn) - M(m,n)|                               !
      !---------------------------------------------------------------!
      double precision function mat_dev(nBas, M, k)
        integer, intent(in) :: nBas, k
        double precision, intent(in) :: M(nBas,nBas)
        integer :: i, j
        double precision :: d
        d = 0.d0
        if (any(sym_perm(:,k) == 0)) then
          mat_dev = huge(1.d0) ; return
        end if
        !$omp parallel do default(shared) private(i,j) reduction(max:d)
        do j = 1, nBas
          do i = 1, nBas
            d = max(d, abs(sym_sign(i,k)*sym_sign(j,k)*M(sym_perm(i,k),sym_perm(j,k)) - M(i,j)))
          end do
        end do
        !$omp end parallel do
        mat_dev = d
      end function mat_dev

      !---------------------------------------------------------------!
      ! ERI check. ERI(n1,nBas,nBas,nBas): n1 = nBas (all indices) or  !
      ! n1 = nfuc (SAO: first index in cell 0; the image quartet is    !
      ! translated back so that its first index is in cell 0).         !
      ! All quartets if n1*nBas^3 <= max(2e5, 2e7/nop), else as many   !
      ! random ones.                                                   !
      !---------------------------------------------------------------!
      double precision function eri_dev(nBas, nfuc, n1, ERI, k)
        integer, intent(in) :: nBas, nfuc, n1, k
        double precision, intent(in) :: ERI(n1,nBas,nBas,nBas)
        integer :: p, q, r, s, pp(4), it
        integer(8) :: ntot, nsamp
        double precision :: d, val
        double precision, allocatable :: rnd(:,:)
        if (any(sym_perm(:,k) == 0)) then
          eri_dev = huge(1.d0) ; return
        end if
        d = 0.d0
        ntot  = int(n1,8)*int(nBas,8)**3
        nsamp = max(200000_8, 20000000_8/int(max(sym_nop,1),8))    ! <= ~2e7 quartets in total
        if (ntot <= nsamp) then
          !$omp parallel do default(shared) private(p,q,r,s,pp,val) reduction(max:d) collapse(2)
          do s = 1, nBas
            do r = 1, nBas
              do q = 1, nBas
                do p = 1, n1
                  call image_quartet(p, q, r, s, pp, val)
                  d = max(d, abs(val*ERI(pp(1),pp(2),pp(3),pp(4)) - ERI(p,q,r,s)))
                end do
              end do
            end do
          end do
          !$omp end parallel do
        else
          allocate(rnd(4,nsamp))
          call random_number(rnd)
          !$omp parallel do default(shared) private(it,p,q,r,s,pp,val) reduction(max:d)
          do it = 1, int(nsamp)
            p = min(n1,   1 + int(rnd(1,it)*n1))
            q = min(nBas, 1 + int(rnd(2,it)*nBas))
            r = min(nBas, 1 + int(rnd(3,it)*nBas))
            s = min(nBas, 1 + int(rnd(4,it)*nBas))
            call image_quartet(p, q, r, s, pp, val)
            d = max(d, abs(val*ERI(pp(1),pp(2),pp(3),pp(4)) - ERI(p,q,r,s)))
          end do
          !$omp end parallel do
          deallocate(rnd)
        end if
        eri_dev = d
      contains
        subroutine image_quartet(p, q, r, s, pp, sg)
          integer, intent(in) :: p, q, r, s
          integer, intent(out) :: pp(4)
          double precision, intent(out) :: sg
          integer :: c(3), i
          pp = [sym_perm(p,k), sym_perm(q,k), sym_perm(r,k), sym_perm(s,k)]
          sg = dble(sym_sign(p,k)*sym_sign(q,k)*sym_sign(r,k)*sym_sign(s,k))
          if (n1 < nBas) then                               ! bring the first index to cell 0
            call cell_of(pp(1), nfuc, c)
            do i = 1, 4
              pp(i) = shift_ao(pp(i), nfuc, c)
            end do
          end if
        end subroutine image_quartet
      end function eri_dev

      !---------------------------------------------------------------!
      ! printing helpers                                               !
      !---------------------------------------------------------------!
      character(len=28) function op_string(R)
        integer, intent(in) :: R(3,3)
        character(len=1), parameter :: ax(3) = ['x','y','z']
        character(len=9) :: c(3)
        integer :: i, j
        do i = 1, 3
          do j = 1, 3
            if (R(i,j) == 1)  c(i) = ' '//ax(j)
            if (R(i,j) == -1) c(i) = '-'//ax(j)
          end do
        end do
        op_string = '('//trim(c(1))//','//trim(c(2))//','//trim(c(3))//')'
      end function op_string

      character(len=5) function op_type(R, t)
        integer, intent(in) :: R(3,3)
        double precision, intent(in) :: t(3)
        integer :: tr, dt
        tr = R(1,1) + R(2,2) + R(3,3)
        dt = det3(R)
        op_type = '?'
        if (dt == 1) then
          select case (tr)
          case (3)  ; op_type = 'E'
          case (-1) ; op_type = 'C2'
          case (0)  ; op_type = 'C3'
          case (1)  ; op_type = 'C4'
          end select
        else
          select case (tr)
          case (-3) ; op_type = 'i'
          case (1)  ; op_type = 'm'
          case (-1) ; op_type = 'S4'
          case (0)  ; op_type = 'S6'
          end select
        end if
        if (any(abs(t) > pos_tol) .and. op_type /= 'i') op_type = trim(op_type)//'+t'
        if (any(abs(t) > pos_tol) .and. op_type == 'i') op_type = 'i+t'
      end function op_type

      integer function det3(R)
        integer, intent(in) :: R(3,3)
        det3 = R(1,1)*(R(2,2)*R(3,3)-R(2,3)*R(3,2)) - R(1,2)*(R(2,1)*R(3,3)-R(2,3)*R(3,1)) &
             + R(1,3)*(R(2,1)*R(3,2)-R(2,2)*R(3,1))
      end function det3

      integer function n_point_ops()
        integer :: k, j
        logical :: new
        n_point_ops = 0
        do k = 1, sym_nop
          new = .true.
          do j = 1, k-1
            if (all(sym_R(:,:,j) == sym_R(:,:,k))) new = .false.
          end do
          if (new) n_point_ops = n_point_ops + 1
        end do
      end function n_point_ops

      ! Schoenflies name of the point group (subgroups of Oh made of signed permutations)
      character(len=8) function point_group_name()
        integer :: k, j, tr, dt, h, ni, nm, nc2, nc3, nc4, ns4, ns6
        logical :: new
        h = 0 ; ni = 0 ; nm = 0 ; nc2 = 0 ; nc3 = 0 ; nc4 = 0 ; ns4 = 0 ; ns6 = 0
        do k = 1, sym_nop
          new = .true.
          do j = 1, k-1
            if (all(sym_R(:,:,j) == sym_R(:,:,k))) new = .false.
          end do
          if (.not. new) cycle
          h = h + 1
          tr = sym_R(1,1,k) + sym_R(2,2,k) + sym_R(3,3,k)
          dt = det3(sym_R(:,:,k))
          if (dt == 1) then
            if (tr == -1) nc2 = nc2 + 1
            if (tr == 0)  nc3 = nc3 + 1
            if (tr == 1)  nc4 = nc4 + 1
          else
            if (tr == -3) ni = ni + 1
            if (tr == 1)  nm = nm + 1
            if (tr == -1) ns4 = ns4 + 1
            if (tr == 0)  ns6 = ns6 + 1
          end if
        end do
        point_group_name = '?'
        select case (h)
        case (48) ; point_group_name = 'Oh'
        case (24)
          if (ni == 1) then ; point_group_name = 'Th'
          else if (nm > 0) then ; point_group_name = 'Td'
          else ; point_group_name = 'O' ; end if
        case (16) ; point_group_name = 'D4h'
        case (12)
          if (ni == 1) then ; point_group_name = 'D3d'
          else ; point_group_name = 'T' ; end if
        case (8)
          if (nc4 > 0 .and. ni == 1) then ; point_group_name = 'C4h'
          else if (nc4 > 0 .and. nm == 4) then ; point_group_name = 'C4v'
          else if (nc4 > 0) then ; point_group_name = 'D4'
          else if (ns4 > 0) then ; point_group_name = 'D2d'
          else ; point_group_name = 'D2h' ; end if
        case (6)
          if (ni == 1) then ; point_group_name = 'S6'
          else if (nm > 0) then ; point_group_name = 'C3v'
          else ; point_group_name = 'D3' ; end if
        case (4)
          if (nc4 > 0) then ; point_group_name = 'C4'
          else if (ns4 > 0) then ; point_group_name = 'S4'
          else if (ni == 1) then ; point_group_name = 'C2h'
          else if (nm == 2) then ; point_group_name = 'C2v'
          else ; point_group_name = 'D2' ; end if
        case (3) ; point_group_name = 'C3'
        case (2)
          if (ni == 1) then ; point_group_name = 'Ci'
          else if (nm == 1) then ; point_group_name = 'Cs'
          else ; point_group_name = 'C2' ; end if
        case (1) ; point_group_name = 'C1'
        end select
      end function point_group_name

end module symmetry_clifford
