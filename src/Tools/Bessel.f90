subroutine bessel_I_scaled_backward(Nmax, x, I_scaled)

      use bessel_functions

      implicit none
      integer, intent(in)           :: Nmax
      double precision, intent(in)  :: x
      double precision, intent(out) :: I_scaled(0:Nmax)
      integer                       :: n
      double precision              :: scale_factor , I_0_normalized

      ! ============================================================
      ! SPECIAL CASE: x = 0
      ! I_0(0) = 1, I_n(0) = 0 for n > 0
      ! Scaled versions: exp(-0) * I_n(0) = I_n(0)
      ! ============================================================
      if (x < 1.d-10) then
        I_scaled(0) = 1.d0
        do n = 1, Nmax
          I_scaled(n) = 0.d0
        end do
        return
      end if      
      ! ============================================================
      ! NORMAL CASE: x > 0
      ! ============================================================
      
      ! Initialize
      I_scaled(Nmax)    = iv_scaled(Nmax, x)
      I_scaled(Nmax-1)  = iv_scaled(Nmax-1, x)
      I_0_normalized    = iv_scaled(0, x)
      
      ! Backward recursion
      do n = Nmax-1, 1, -1
        I_scaled(n-1) = I_scaled(n+1) + (2.d0 * dble(n) / x) * I_scaled(n)
      end do
      
      ! Normalize
      scale_factor = I_0_normalized / I_scaled(0)
      
      ! Apply normalization
      I_scaled(0:Nmax) = I_scaled(0:Nmax) * scale_factor  ! Array notation

      
end subroutine bessel_I_scaled_backward