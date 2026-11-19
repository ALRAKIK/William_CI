!=====================================================================!
!   [Nov 2026] keyword rys_quad : start-up self-check                 !
!                                                                     !
!   Computes a sample of shell quartets twice, with the Rys-type      !
!   nodes and with the global t tables, and switches rys_quad off     !
!   for the whole run when the largest difference exceeds 1e-11.      !
!   Sample: every on-site quartet (SI SI|SI SI) of the unit cell plus  !
!   up to 300 pseudo-random quartets (SI in the unit cell).           !
!   Nothing to set by hand: if the check fails, the run continues     !
!   with the usual quadrature and a note is printed.                  !
!=====================================================================!
subroutine rys_self_check(ndp, nBas, AO)

      use files,             only : outfile
      use unitcell_module,   only : nfuc
      use classification_ERI
      use clifford_tools,    only : rys_on, orbital_to_lxyz

      implicit none

      integer           , intent(in) :: ndp, nBas
      type(ERI_function), intent(in) :: AO(nBas)

      double precision, parameter :: tol = 1.d-11
      integer,          parameter :: nrand = 300
      integer,          allocatable :: sh_first(:), sh_n(:)
      double precision, allocatable :: b1(:,:,:,:), b2(:,:,:,:)
      integer          :: i, j, nsh, nshuc, lxyz(3), lt, ncomp, q(4), ns, it
      integer(8)       :: seed
      logical          :: same
      double precision :: dmax

      if (.not. rys_on) return

      ! shells (same rule as the ERI driver)
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

      dmax = 0.d0 ; ns = 0 ; seed = 12345_8
      do it = 1 , nshuc + nrand
        if (it <= nshuc) then
          q = it
        else
          q(1) = 1 + int(mod(next(), int(nshuc,8)))
          q(2) = 1 + int(mod(next(), int(nsh,8)))
          q(3) = 1 + int(mod(next(), int(nsh,8)))
          q(4) = 1 + int(mod(next(), int(nsh,8)))
        end if
        allocate(b1(sh_n(q(1)),sh_n(q(2)),sh_n(q(3)),sh_n(q(4))), b2(sh_n(q(1)),sh_n(q(2)),sh_n(q(3)),sh_n(q(4))))
        rys_on = .true.
        call ERI_shell_block_toroidal(ndp,sh_n(q(1)),AO(sh_first(q(1)):),sh_n(q(2)),AO(sh_first(q(2)):), &
                                      sh_n(q(3)),AO(sh_first(q(3)):),sh_n(q(4)),AO(sh_first(q(4)):),b1)
        rys_on = .false.
        call ERI_shell_block_toroidal(ndp,sh_n(q(1)),AO(sh_first(q(1)):),sh_n(q(2)),AO(sh_first(q(2)):), &
                                      sh_n(q(3)),AO(sh_first(q(3)):),sh_n(q(4)),AO(sh_first(q(4)):),b2)
        dmax = max(dmax, maxval(abs(b1 - b2)))
        ns = ns + 1
        deallocate(b1, b2)
      end do

      rys_on = (dmax <= tol)
      write(outfile,'(a,i0,a,es9.2,a)') ' rys_quad self-check: ', ns, ' shell quartets, largest difference ', dmax, &
                                        merge(' -> Rys-type nodes ON          ', ' -> above 1e-11: rys_quad OFF  ', rys_on)
      write(*,'(a,i0,a,es9.2,a)') ' rys_quad self-check: ', ns, ' shell quartets, largest difference ', dmax, &
                                  merge(' -> Rys-type nodes ON          ', ' -> above 1e-11: rys_quad OFF  ', rys_on)
      flush(outfile)
      deallocate(sh_first, sh_n)

contains

      integer(8) function next()
        seed = mod(seed*1103515245_8 + 12345_8, 2147483648_8)     ! 0 <= seed < 2^31, no overflow
        next = seed/16
      end function next

end subroutine rys_self_check
