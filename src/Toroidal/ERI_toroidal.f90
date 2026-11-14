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
      use clifford_tools, only : init_clifford, orbital_to_lxyz, setup_t_tables, bessel_nmax_estimate

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
      double precision               :: emax, amin, a_dir(3), start_time, end_time, val
      logical                        :: same

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

      start_time = omp_get_wtime()

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

      end_time = omp_get_wtime()
      write(outfile,'(A65,5X,F12.3,a)') '2 el-integrals wall time = ', end_time - start_time, " s"
      write(outfile,'(6x,a,i8,a,i8)') 'shells: ', nsh, '   bra shell pairs: ', npair
      flush(outfile)

      if (.not. c_SAO) then
        call symmetry_of_integrals_ERI(nBas,nfuc,two_electron_integrals)        ! rows i > nfuc, in place
        if (c_debug) call check_eri_symmetry(nBas, two_electron_integrals, outfile)
      end if

      deallocate(sh_first, sh_n, Q_schwarz, Qsh, ip_I, ip_J)

end subroutine ERI_integral_toroidal
