!=====================================================================!
!                                                                     !
!       Toroidal 1D/2D/3D : TWO-ELECTRON INTEGRALS (driver)           !
!                                                                     !
!   ERI_integral_toroidal (ndp, nBas, AO, ERI)                        !
!     ndp = 1 (Tori1D), 2 (Tori2D), 3 (Tori3D)                        !
!                                                                     !
!   1) t-tables (setup_t_tables)                                      !
!   2) shells = consecutive Cartesian components with the same       !
!      centre and exponents (as built by classification_tor)          !
!   3) Schwarz bounds from the diagonal shell-pair blocks             !
!   4) loop over shell quartets (bra shell SI in the unit cell),      !
!      ERI_shell_block_toroidal, keep i<=nfuc, i<=j, k<=l, i<=k and   !
!      store with the 8-fold symmetry (chemist or physicist)          !
!   5) normal path: rows i > nfuc by translation, in place            !
!                                                                     !
!   [Nov 2026] symmetry-unique quartets: a shell quartet is computed  !
!   only if it is the smallest member of its orbit under the 8        !
!   permutations x the lattice translations (cells x fastest, then y, !
!   then z). Each value is scattered to all 8 permutations, shifted   !
!   so that the first index lies in the unit cell. ~3.5x fewer kernel !
!   calls. Falls back to the original loop if the shells are not the  !
!   same in every cell.                                               !
!                                                                     !
!   AO must be the SAME normalized list used for S, T and V.          !
!=====================================================================!

subroutine ERI_integral_toroidal(ndp,nBas,AO,two_electron_integrals)

      use files
      use keywords
      use omp_lib
      use torus_init
      use classification_ERI
      use unitcell_module
      use debug_tools
      use clifford_tools, only : init_clifford, orbital_to_lxyz, setup_t_tables, bessel_nmax_estimate, &
                                 ps_on, ps_noscreen, ps_qp, ps_bra, ps_ket

      implicit none

      ! - input - !

      integer           ,intent(in)  :: ndp                 ! periodic directions 1..ndp
      integer           ,intent(in)  :: nBas
      type(ERI_function),intent(in)  :: AO(nBas)          ! normalized AOs (same as S, T, V)

      ! - output - !

      double precision  ,intent(out) :: two_electron_integrals(ERI_size,nBas,nBas,nBas)

      ! - local - !

      double precision, parameter    :: schwarz_thresh = 1.0d-13
      double precision,allocatable   :: Q_schwarz(:,:), Qsh(:,:), blk(:,:,:,:)
      integer, allocatable           :: sh_first(:), sh_n(:), ip_I(:), ip_J(:)
      integer                        :: i, j, k, l, a, b, c, d, nsh, lxyz(3), lt, ncomp, lmax, ncap, npair, ip
      integer                        :: SI, SJ, SK, SL
      double precision               :: emax, amin, a_dir(3), val
      logical                        :: same
      logical                        :: sym_ok
      integer                        :: nshuc, q4(4), qc(4), pp(4), ii, iperm, ncl(3)
      double precision               :: Qmax
      integer                        :: mp, p1, p2, ab(4)
      type(ERI_function), allocatable :: A1(:), B1(:)

      !-----------------------------------------------------------------!


      ! ---- t-table: largest Bessel order needed by any pair ---- !
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
      call setup_t_tables(ncap, ndp)

      ! ---- shells: consecutive Cartesian components, same centre and exponents ---- !
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
        if (.not. same) ncomp = 1          ! fall back: single function treated as its own "shell"
        nsh = nsh + 1
        sh_first(nsh) = i ; sh_n(nsh) = ncomp
        i = i + ncomp
      end do

      ! ---- Schwarz: diagonal shell-pair blocks ---- !
      allocate(Q_schwarz(nBas,nBas), Qsh(nsh,nsh))
      Q_schwarz = 0.d0
      !$omp parallel do default(shared) private(SI,SJ,blk,a,b,i,j) schedule(dynamic,1)
      do SI = 1 , nsh
        do SJ = SI , nsh
          allocate(blk(sh_n(SI),sh_n(SJ),sh_n(SI),sh_n(SJ)))
          call ERI_shell_block_toroidal(ndp,sh_n(SI),AO(sh_first(SI):),sh_n(SJ),AO(sh_first(SJ):), &
                                     sh_n(SI),AO(sh_first(SI):),sh_n(SJ),AO(sh_first(SJ):),blk)
          do a = 1 , sh_n(SI)
            do b = 1 , sh_n(SJ)
              i = sh_first(SI)+a-1 ; j = sh_first(SJ)+b-1
              Q_schwarz(i,j) = dsqrt(dabs(blk(a,b,a,b)))
              Q_schwarz(j,i) = Q_schwarz(i,j)
            end do
          end do
          deallocate(blk)
        end do
      end do
      !$omp end parallel do
      do SI = 1 , nsh
        do SJ = 1 , nsh
          Qsh(SI,SJ) = maxval(Q_schwarz(sh_first(SI):sh_first(SI)+sh_n(SI)-1, sh_first(SJ):sh_first(SJ)+sh_n(SJ)-1))
        end do
      end do

      ! ---- symmetry-unique shell quartets (translation x 8-fold) ---- !
      ncl = (/ max(nx,1), max(ny,1), max(nz,1) /)
      nshuc = 0
      do SI = 1 , nsh
        if (sh_first(SI) <= nfuc) nshuc = nshuc + 1
      end do
      sym_ok = (nsh == nshuc*ncl(1)*ncl(2)*ncl(3)) .and. (nBas == nfuc*ncl(1)*ncl(2)*ncl(3))
      if (sym_ok) then
        do SI = 1 , nsh - nshuc               ! every cell repeats the shells of the unit cell
          if (sh_n(SI+nshuc) /= sh_n(SI) .or. sh_first(SI+nshuc) /= sh_first(SI) + nfuc) sym_ok = .false.
        end do
      end if

      if (sym_ok) then

        Qmax = maxval(Qsh)

        ! ---- keyword prim_screen: primitive-pair Schwarz table ---- !
        ! built with the kernel's own prefactor screen switched off, so  !
        ! that small diagonal integrals are not set to zero.             !
        ps_on = .false.
        if (c_prim_screen) then
          mp = 0
          do SI = 1 , nsh
            mp = max(mp, size(AO(sh_first(SI))%exponent))
          end do
          allocate(ps_qp(mp,mp,nshuc,nsh)) ; ps_qp = 0.d0
          ps_noscreen = .true.
          !$omp parallel do default(shared) private(SI,SJ,p1,p2,A1,B1,blk,a,b) schedule(dynamic,1) collapse(2)
          do SI = 1 , nshuc
            do SJ = 1 , nsh
              allocate(A1(sh_n(SI)), B1(sh_n(SJ)), blk(sh_n(SI),sh_n(SJ),sh_n(SI),sh_n(SJ)))
              do p1 = 1 , size(AO(sh_first(SI))%exponent)
                do p2 = 1 , size(AO(sh_first(SJ))%exponent)
                  do a = 1 , sh_n(SI)
                    A1(a) = AO(sh_first(SI)+a-1)
                    A1(a)%exponent    = (/ AO(sh_first(SI)+a-1)%exponent(p1) /)
                    A1(a)%coefficient = (/ AO(sh_first(SI)+a-1)%coefficient(p1) /)
                  end do
                  do b = 1 , sh_n(SJ)
                    B1(b) = AO(sh_first(SJ)+b-1)
                    B1(b)%exponent    = (/ AO(sh_first(SJ)+b-1)%exponent(p2) /)
                    B1(b)%coefficient = (/ AO(sh_first(SJ)+b-1)%coefficient(p2) /)
                  end do
                  call ERI_shell_block_toroidal(ndp,sh_n(SI),A1,sh_n(SJ),B1,sh_n(SI),A1,sh_n(SJ),B1,blk)
                  do a = 1 , sh_n(SI)
                    do b = 1 , sh_n(SJ)
                      ps_qp(p1,p2,SI,SJ) = max(ps_qp(p1,p2,SI,SJ), dsqrt(dabs(blk(a,b,a,b))))
                    end do
                  end do
                end do
              end do
              deallocate(A1, B1, blk)
            end do
          end do
          !$omp end parallel do
          ps_noscreen = .false.
          ps_on = .true.
        end if
        two_electron_integrals = 0.d0

        !$omp parallel do default(shared) private(SI,SJ,SK,SL,blk,a,b,c,d,i,j,k,l,val,q4,qc,pp,ii,iperm,ab) &
        !$omp schedule(dynamic,1) collapse(2)
        do SI = 1 , nshuc
          do SJ = 1 , nsh
            if (mod(SJ-1,nshuc)+1 < SI) cycle                         ! canonical: SI is the smallest local index
            if (Qsh(SI,SJ)*Qmax < schwarz_thresh) cycle
            do SK = 1 , nsh
              if (mod(SK-1,nshuc)+1 < SI) cycle
              do SL = 1 , nsh
                if (mod(SL-1,nshuc)+1 < SI) cycle
                if (Qsh(SI,SJ)*Qsh(SK,SL) < schwarz_thresh) cycle
                q4 = (/ SI, SJ, SK, SL /)
                call canonical_quartet_toroidal(q4, nsh, nshuc, ncl, qc)
                if (any(qc /= q4)) cycle
                allocate(blk(sh_n(SI),sh_n(SJ),sh_n(SK),sh_n(SL)))
                if (ps_on) then                                   ! primitive tables of the bra and ket pairs
                  ps_bra => ps_qp(:,:,SI,SJ)
                  ab = (/ SK, SL, SK, SK /)
                  call translate_to_unit_cell(ab, nshuc, ncl)     ! (SK,SL) shifted so that SK is in the unit cell
                  ps_ket => ps_qp(:,:,ab(1),ab(2))
                end if
                call ERI_shell_block_toroidal(ndp,sh_n(SI),AO(sh_first(SI):),sh_n(SJ),AO(sh_first(SJ):), &
                                           sh_n(SK),AO(sh_first(SK):),sh_n(SL),AO(sh_first(SL):),blk)
                do d = 1 , sh_n(SL)
                 do c = 1 , sh_n(SK)
                  do b = 1 , sh_n(SJ)
                   do a = 1 , sh_n(SI)
                    i = sh_first(SI)+a-1 ; j = sh_first(SJ)+b-1
                    k = sh_first(SK)+c-1 ; l = sh_first(SL)+d-1
                    if (Q_schwarz(i,j)*Q_schwarz(k,l) < schwarz_thresh) cycle   ! same screen as before
                    val = blk(a,b,c,d)
                    do iperm = 1 , 8
                      call permute_quartet(iperm, (/ i, j, k, l /), pp)
                      call translate_to_unit_cell(pp, nfuc, ncl)     ! (ij|kl) chemist, pp(1) in the unit cell
                      if (c_physicist) then
                        two_electron_integrals(pp(1),pp(3),pp(2),pp(4)) = val
                      else
                        two_electron_integrals(pp(1),pp(2),pp(3),pp(4)) = val
                      end if
                    end do
                   end do
                  end do
                 end do
                end do
                deallocate(blk)
              end do
            end do
          end do
        end do
        !$omp end parallel do
        ps_on = .false.
        if (allocated(ps_qp)) deallocate(ps_qp)

      else

      ! ---- bra shell pairs (SI <= SJ, SI touching the unit cell) ---- !
      npair = 0
      do SI = 1 , nsh
        if (sh_first(SI) > nfuc) cycle
        do SJ = SI , nsh
          npair = npair + 1
        end do
      end do
      allocate(ip_I(npair), ip_J(npair))
      npair = 0
      do SI = 1 , nsh
        if (sh_first(SI) > nfuc) cycle
        do SJ = SI , nsh
          npair = npair + 1
          ip_I(npair) = SI ; ip_J(npair) = SJ
        end do
      end do

      ! rows i <= nfuc are computed here, directly in the output array.
      ! SAO : the output IS (nfuc,nBas,nBas,nBas) ; normal : the rows
      ! i > nfuc are filled afterwards by translation (in place).
      two_electron_integrals = 0.d0

      !$omp parallel do default(shared) private(ip,SI,SJ,SK,SL,blk,a,b,c,d,i,j,k,l,val) schedule(dynamic,1)
      do ip = 1 , npair
        SI = ip_I(ip) ; SJ = ip_J(ip)
        do SK = 1 , nsh
          if (sh_first(SI) > sh_first(SK) + sh_n(SK) - 1) cycle        ! need some i <= k
          do SL = SK , nsh
            if (Qsh(SI,SJ)*Qsh(SK,SL) < schwarz_thresh) cycle
            allocate(blk(sh_n(SI),sh_n(SJ),sh_n(SK),sh_n(SL)))
            call ERI_shell_block_toroidal(ndp,sh_n(SI),AO(sh_first(SI):),sh_n(SJ),AO(sh_first(SJ):), &
                                       sh_n(SK),AO(sh_first(SK):),sh_n(SL),AO(sh_first(SL):),blk)
            do d = 1 , sh_n(SL)
             do c = 1 , sh_n(SK)
              do b = 1 , sh_n(SJ)
               do a = 1 , sh_n(SI)
                i = sh_first(SI)+a-1 ; j = sh_first(SJ)+b-1
                k = sh_first(SK)+c-1 ; l = sh_first(SL)+d-1
                ! exactly the quartets of the original driver
                if (i > nfuc .or. j < i .or. l < k .or. i > k) cycle
                if (Q_schwarz(i,j)*Q_schwarz(k,l) < schwarz_thresh) cycle
                val = blk(a,b,c,d)
                if (c_physicist) then
                  if (i <= nfuc) two_electron_integrals(i,k,j,l) = val
                  if (j <= nfuc) two_electron_integrals(j,k,i,l) = val
                  if (i <= nfuc) two_electron_integrals(i,l,j,k) = val
                  if (j <= nfuc) two_electron_integrals(j,l,i,k) = val
                  if (k <= nfuc) two_electron_integrals(k,i,l,j) = val
                  if (l <= nfuc) two_electron_integrals(l,i,k,j) = val
                  if (k <= nfuc) two_electron_integrals(k,j,l,i) = val
                  if (l <= nfuc) two_electron_integrals(l,j,k,i) = val
                else
                  if (i <= nfuc) two_electron_integrals(i,j,k,l) = val
                  if (j <= nfuc) two_electron_integrals(j,i,k,l) = val
                  if (i <= nfuc) two_electron_integrals(i,j,l,k) = val
                  if (j <= nfuc) two_electron_integrals(j,i,l,k) = val
                  if (k <= nfuc) two_electron_integrals(k,l,i,j) = val
                  if (l <= nfuc) two_electron_integrals(l,k,i,j) = val
                  if (k <= nfuc) two_electron_integrals(k,l,j,i) = val
                  if (l <= nfuc) two_electron_integrals(l,k,j,i) = val
                end if
               end do
              end do
             end do
            end do
            deallocate(blk)
          end do
        end do
      end do
      !$omp end parallel do

      deallocate(ip_I, ip_J)

      end if   ! sym_ok

      if (.not. c_SAO) then
        call symmetry_of_integrals_ERI(nBas,nfuc,two_electron_integrals)        ! rows i > nfuc, in place
        if (c_debug) call check_eri_symmetry(nBas, two_electron_integrals, outfile)
      end if

      deallocate(sh_first, sh_n, Q_schwarz, Qsh)

end subroutine ERI_integral_toroidal

!---------------------------------------------------------------------!
! 8 index permutations of a chemist quartet (ij|kl)                   !
!---------------------------------------------------------------------!
subroutine permute_quartet(iperm, q, p)

      implicit none
      integer, intent(in)  :: iperm, q(4)
      integer, intent(out) :: p(4)

      select case (iperm)
      case (1) ; p = (/ q(1), q(2), q(3), q(4) /)
      case (2) ; p = (/ q(2), q(1), q(3), q(4) /)
      case (3) ; p = (/ q(1), q(2), q(4), q(3) /)
      case (4) ; p = (/ q(2), q(1), q(4), q(3) /)
      case (5) ; p = (/ q(3), q(4), q(1), q(2) /)
      case (6) ; p = (/ q(4), q(3), q(1), q(2) /)
      case (7) ; p = (/ q(3), q(4), q(2), q(1) /)
      case default ; p = (/ q(4), q(3), q(2), q(1) /)
      end select

end subroutine permute_quartet

!---------------------------------------------------------------------!
! Shift the four indices (functions or shells, nblk per cell)     !
! by the same lattice vector so that p(1) is in the unit cell.   !
! Cells are numbered with x fastest, then y, then z (super_molecule). !
!---------------------------------------------------------------------!
subroutine translate_to_unit_cell(p, nblk, ncl)

      implicit none
      integer, intent(inout) :: p(4)
      integer, intent(in)    :: nblk, ncl(3)
      integer :: ii, c, loc, cx, cy, cz, tx, ty, tz

      c  = (p(1)-1)/nblk
      tx = mod(c, ncl(1)) ; ty = mod(c/ncl(1), ncl(2)) ; tz = c/(ncl(1)*ncl(2))
      do ii = 1 , 4
        c   = (p(ii)-1)/nblk
        loc = mod(p(ii)-1, nblk)
        cx  = modulo(mod(c, ncl(1))         - tx, ncl(1))
        cy  = modulo(mod(c/ncl(1), ncl(2))  - ty, ncl(2))
        cz  = modulo(c/(ncl(1)*ncl(2))      - tz, ncl(3))
        p(ii) = (cx + ncl(1)*(cy + ncl(2)*cz))*nblk + loc + 1
      end do

end subroutine translate_to_unit_cell

!---------------------------------------------------------------------!
! Smallest member (lexicographic) of the orbit of a shell quartet     !
! under the 8 permutations x lattice translations.                    !
!---------------------------------------------------------------------!
subroutine canonical_quartet_toroidal(q, nsh, nshuc, ncl, best)

      implicit none
      integer, intent(in)  :: q(4), nsh, nshuc, ncl(3)
      integer, intent(out) :: best(4)
      integer :: p(4), iperm, ii

      best = huge(1)
      do iperm = 1 , 8
        call permute_quartet(iperm, q, p)
        call translate_to_unit_cell(p, nshuc, ncl)
        do ii = 1 , 4
          if (p(ii) < best(ii)) then
            best = p ; exit
          end if
          if (p(ii) > best(ii)) exit
        end do
      end do

end subroutine canonical_quartet_toroidal
