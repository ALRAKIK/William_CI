!=====================================================================!
!   [Nov 2026] keyword fast_quad : start-up self-check                !
!                                                                     !
!   fast_quad uses 32/40 mapped Gauss-Legendre t-points instead of    !
!   the 64/128 point tables. That is accurate for soft basis sets     !
!   (STO-3G, H chains: ~1e-11) but NOT for very tight primitives:     !
!   cc-pVTZ carbon (exponents 8236, 1235) gave 1.7e-4 Eh in HF.       !
!                                                                     !
!   Here a sample of nuclear attraction elements and a sample of ERI  !
!   shell quartets are computed with both rules (V is the most        !
!   sensitive: point nucleus = infinite exponent).                    !
!   [Nov 2026] V sample (was: the whole matrix twice, too expensive   !
!   for large cells, e.g. 1728 atoms): every diagonal V(i,i) of the   !
!   unit cell (always includes the tightest functions, at most 2000), !
!   200 random pairs inside the unit cell and 200 random pairs        !
!   (i in the unit cell, j anywhere). ERI sample:                     !
!   every on-site quartet (SI SI|SI SI) of the unit cell (this always !
!   includes the tightest shells) plus up to 300 pseudo-random        !
!   quartets. If the largest difference exceeds tol, fast_quad is     !
!   switched off for the whole run (nuclear attraction and ERIs) and  !
!   the usual 64/128 point tables are used. A note is printed in      !
!   both cases. Nothing to set by hand.                               !
!                                                                     !
!   Call: after normalisation of the AOs and Torus_def, before the    !
!   nuclear attraction (Tori1D / Tori2D / Tori3D).                    !
!=====================================================================!
subroutine fast_quad_self_check(ndp, nBas, AO, n_atoms, geometry, atoms)

      use files,              only : outfile
      use unitcell_module,    only : nfuc
      use classification_ERI
      use atom_basis
      use torus_init,         only : ax, ay, az
      use keywords,           only : c_nquad_128
      use precomputed_bessel, only : fast_quad, nquad_point
      use clifford_tools,     only : orbital_to_lxyz, init_clifford, setup_t_tables, &
                                     bessel_nmax_estimate

      implicit none

      integer           , intent(in) :: ndp, nBas
      type(ERI_function), intent(in) :: AO(nBas)
      integer           , intent(in) :: n_atoms
      double precision  , intent(in) :: geometry(n_atoms,3)
      type(atom)        , intent(in) :: atoms(n_atoms)

      double precision, parameter :: tol = 1.d-10
      integer,          parameter :: nrand = 300
      integer,          allocatable :: sh_first(:), sh_n(:)
      double precision, allocatable :: b1(:,:,:,:), b2(:,:,:,:)
      integer          :: i, j, nsh, nshuc, lxyz(3), lt, ncomp, q(4), ns, it, lmax, ncap
      integer          :: nq_fast, nq_full
      integer(8)       :: seed
      logical          :: same, ok
      double precision :: dmax, emax, amin, a_dir(3)
      integer, allocatable :: qs(:,:)
      double precision, allocatable :: buf(:)
      integer          :: ntot, off, nb
      double precision, allocatable :: V1(:), V2(:)
      double precision :: dV
      integer, parameter :: nv_rand = 200
      integer          :: npv, ip
      integer, allocatable :: pv(:,:), lxv(:,:)

      if (.not. fast_quad) return

      nq_fast = nquad_point
      nq_full = merge(128, 64, c_nquad_128)

      ! ---- same table size as the ERI driver ---- !
      call init_clifford()
      lmax = 0 ; emax = 0.d0
      do i = 1 , nBas
        call orbital_to_lxyz(AO(i)%orbital, lxyz)
        lmax = max(lmax, maxval(lxyz(1:ndp)))
        emax = max(emax, maxval(AO(i)%exponent))
      end do
      a_dir = (/ ax, ay, az /)
      amin  = minval(a_dir(1:ndp))
      ncap = bessel_nmax_estimate(4.d0*emax/(amin*amin)) + 2*lmax + 1

      ! ---- shells (same rule as the ERI driver) ---- !
      allocate(sh_first(nBas), sh_n(nBas))
      nsh = 0 ; i = 1
      do while (i <= nBas)
        call orbital_to_lxyz(AO(i)%orbital, lxyz)
        lt = sum(lxyz) ; ncomp = (lt+1)*(lt+2)/2
        same = (i + ncomp - 1 <= nBas)
        if (same) then
          do j = i+1 , i+ncomp-1
            call orbital_to_lxyz(AO(j)%orbital, lxyz)
            same = same .and. sum(lxyz) == lt .and. AO(j)%x == AO(i)%x .and. AO(j)%y == AO(i)%y &
                   .and. AO(j)%z == AO(i)%z .and. size(AO(j)%exponent) == size(AO(i)%exponent)
            if (same) same = all(AO(j)%exponent == AO(i)%exponent)
          end do
        end if
        if (.not. same) ncomp = 1
        nsh = nsh + 1
        sh_first(nsh) = i ; sh_n(nsh) = ncomp
        i = i + ncomp
      end do
      nshuc = 0
      do i = 1 , nsh
        if (sh_first(i) <= nfuc) nshuc = nshuc + 1
      end do
      nshuc = max(nshuc, 1)

      ! ---- sample of quartets ---- !
      ns = nshuc + nrand
      allocate(qs(4,ns))
      seed = 12345_8
      do it = 1 , ns
        if (it <= nshuc) then
          qs(:,it) = it
        else
          qs(1,it) = 1 + int(mod(next(), int(nshuc,8)))
          qs(2,it) = 1 + int(mod(next(), int(nsh,8)))
          qs(3,it) = 1 + int(mod(next(), int(nsh,8)))
          qs(4,it) = 1 + int(mod(next(), int(nsh,8)))
        end if
      end do

      ! ---- all sample quartets with the fast rule, then with the full tables ---- !
      ! (tables rebuilt only twice; the fast-rule blocks are kept in buf)
      ntot = 0
      do it = 1 , ns
        q = qs(:,it)
        ntot = ntot + sh_n(q(1))*sh_n(q(2))*sh_n(q(3))*sh_n(q(4))
      end do
      allocate(buf(ntot))

      ! ---- sample of nuclear attraction elements (i <= nfuc) ---- !
      npv = min(nfuc, 2000) + 2*nv_rand
      allocate(pv(2,npv), V1(npv), V2(npv), lxv(3,nBas))
      do i = 1 , nBas
        call orbital_to_lxyz(AO(i)%orbital, lxv(:,i))
      end do
      do ip = 1 , min(nfuc, 2000)
        pv(:,ip) = ip                                                    ! diagonal V(i,i)
      end do
      do ip = min(nfuc, 2000) + 1 , npv
        pv(1,ip) = 1 + int(mod(next(), int(nfuc,8)))
        if (ip <= min(nfuc, 2000) + nv_rand) then
          pv(2,ip) = 1 + int(mod(next(), int(nfuc,8)))                  ! same unit cell
        else
          pv(2,ip) = 1 + int(mod(next(), int(nBas,8)))                  ! anywhere
        end if
      end do

      fast_quad = .true.  ; nquad_point = nq_fast
      call setup_t_tables(ncap, ndp)
      call nuclear_attraction(V1)
      off = 0
      do it = 1 , ns
        q = qs(:,it)
        allocate(b1(sh_n(q(1)),sh_n(q(2)),sh_n(q(3)),sh_n(q(4))))
        call ERI_shell_block_toroidal(ndp,sh_n(q(1)),AO(sh_first(q(1)):),sh_n(q(2)),AO(sh_first(q(2)):), &
                                      sh_n(q(3)),AO(sh_first(q(3)):),sh_n(q(4)),AO(sh_first(q(4)):),b1)
        nb = size(b1)
        buf(off+1:off+nb) = reshape(b1, [nb])
        off = off + nb
        deallocate(b1)
      end do

      fast_quad = .false. ; nquad_point = nq_full
      call setup_t_tables(ncap, ndp)
      call nuclear_attraction(V2)
      dV = maxval(abs(V1 - V2))
      deallocate(V1, V2, pv, lxv)
      dmax = 0.d0 ; off = 0
      do it = 1 , ns
        q = qs(:,it)
        allocate(b2(sh_n(q(1)),sh_n(q(2)),sh_n(q(3)),sh_n(q(4))))
        call ERI_shell_block_toroidal(ndp,sh_n(q(1)),AO(sh_first(q(1)):),sh_n(q(2)),AO(sh_first(q(2)):), &
                                      sh_n(q(3)),AO(sh_first(q(3)):),sh_n(q(4)),AO(sh_first(q(4)):),b2)
        nb = size(b2)
        dmax = max(dmax, maxval(abs(buf(off+1:off+nb) - reshape(b2, [nb]))))
        off = off + nb
        deallocate(b2)
      end do
      deallocate(buf)

      ! ---- decision: leave the right tables in place for the whole run ---- !
      ok = (dmax <= tol) .and. (dV <= tol)
      if (ok) then
        fast_quad = .true.  ; nquad_point = nq_fast
      else
        fast_quad = .false. ; nquad_point = nq_full
      end if
      call setup_t_tables(ncap, ndp)

      write(outfile,'(a,i0,a,es9.2,a,es9.2,a)') ' fast_quad self-check: ', ns, ' ERI shell quartets, largest difference ', &
            dmax, ' ; nuclear attraction (sample) ', dV, merge(' -> fast_quad ON                                ', &
                                                     ' -> above 1e-10: fast_quad OFF, full tables used', ok)
      write(*,'(a,i0,a,es9.2,a,es9.2,a)') ' fast_quad self-check: ', ns, ' ERI shell quartets, largest difference ', &
            dmax, ' ; nuclear attraction (sample) ', dV, merge(' -> fast_quad ON                                ', &
                                                     ' -> above 1e-10: fast_quad OFF, full tables used', ok)
      flush(outfile)
      deallocate(sh_first, sh_n, qs)

contains

      ! [Nov 2026] only the sampled elements (pv), with the t-tables in place (same
      ! table size as the nuclear attraction matrix routines: same ncap rule)
      subroutine nuclear_attraction(Vm)
        double precision, intent(out) :: Vm(npv)
        double precision :: r1(3), r2(3)
        integer :: k, a, b
        !$omp parallel do default(shared) private(k,a,b,r1,r2) schedule(dynamic,4)
        do k = 1 , npv
          a = pv(1,k) ; b = pv(2,k)
          r1 = (/ AO(a)%x , AO(a)%y , AO(a)%z /)
          r2 = (/ AO(b)%x , AO(b)%y , AO(b)%z /)
          select case (ndp)
            case (1) ; call nuclear_attraction_integral_toroidal_1D(lxv(:,a),lxv(:,b),n_atoms,geometry,atoms,r1,r2,AO(a),AO(b),Vm(k))
            case (2) ; call nuclear_attraction_integral_toroidal_2D(lxv(:,a),lxv(:,b),n_atoms,geometry,atoms,r1,r2,AO(a),AO(b),Vm(k))
            case (3) ; call nuclear_attraction_integral_toroidal_3D(lxv(:,a),lxv(:,b),n_atoms,geometry,atoms,r1,r2,AO(a),AO(b),Vm(k))
          end select
        end do
        !$omp end parallel do
      end subroutine nuclear_attraction

      integer(8) function next()
        seed = mod(seed*1103515245_8 + 12345_8, 2147483648_8)     ! 0 <= seed < 2^31, no overflow
        next = seed/16
      end function next

end subroutine fast_quad_self_check
