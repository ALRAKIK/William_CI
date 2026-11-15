!---------------------------------------------------------------------!
! [Nov 2026] keyword fast_ccsdt_real                                  !
! gr(i) = real(gc(i)) and imax = max |aimag(gc(i))|, rmax = max |gc|  !
! (threaded, no temporary of the size of the integrals)               !
!---------------------------------------------------------------------!
subroutine eri_to_real(n, gc, gr, imax, rmax)

      use constants_module, only: dp, dpc

      implicit none

      integer(8)  , intent(in)  :: n
      complex(dpc), intent(in)  :: gc(n)
      real(dp)    , intent(out) :: gr(n)
      real(dp)    , intent(out) :: imax, rmax
      integer(8) :: i

      imax = 0.0_dp ; rmax = 0.0_dp
      !$omp parallel do reduction(max:imax,rmax)
      do i = 1, n
        gr(i) = real(gc(i), dp)
        imax  = max(imax, abs(aimag(gc(i))))
        rmax  = max(rmax, abs(gr(i)))
      end do
      !$omp end parallel do

end subroutine eri_to_real
