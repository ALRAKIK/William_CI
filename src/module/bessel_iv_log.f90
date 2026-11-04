!------------------------------------------------------------------------------
! This file is part of CUSF - Fortran port
! 
! CUSF - Math library for Special Functions
! Fortran translation of iv_log.tpp
!
! CUSF is free software; you can redistribute it and/or modify it under the 
! terms of the standard MIT license.
!------------------------------------------------------------------------------

module bessel_functions
      
      use u_funcs
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

      implicit none
      private
    
      ! Public interface
      public                      :: iv_log , iv_scaled , iv_log_scaled
    
      ! Constants
      integer, parameter          :: N_TERMS = 41
      double precision, parameter :: LOG_INV_SQRT_2PI = -0.91893853320467267d0
      double precision, parameter :: LOG_2 = 0.69314718055994529d0
      double precision, parameter :: HUGE_VAL = huge(1.0d0)
    
      contains

      !---------------------------------------------------------------------------
      ! iv_log(v, x)
      !
      ! Return log of the modified Bessel function of the first kind of order v
      ! for real parameters by selecting appropriate numerical method.
      !
      ! Input:
      !   v - Order of the Bessel function
      !   x - Argument of the Bessel function
      ! Output:
      !   log(I_v(x))
      !---------------------------------------------------------------------------

      pure elemental function iv_log(v_in, x) result(log_Iv)

        double precision, intent(in) :: v_in, x
        double precision             :: log_Iv
        
        ! Local variables

        double precision             :: v
        double precision             :: log_x, log_v
        double precision             :: terms(0:N_TERMS-1)
        double precision             :: v_power(0:13)
        double precision             :: t_power(0:N_TERMS-1)
        double precision             :: mu, curr_term, inner_sum
        double precision             :: x_prime, x_prime_2, sqrt_1_plus_x2, t, eta
        double precision             :: log_x2_over_4, x2_over_4, inv_x2_over_4, v_inv_x2_over_4
        double precision             :: sum_terms

        double precision,parameter   :: pi = 3.14159265358979323846D00

        integer                      :: k, c, i
        integer                      :: peak_k
        
        double precision             :: u, asinh_u, t1, v_eta, s

        
        ! Input validation
        if (x < 0.0d0 .or. (v_in < 0.0d0 .and. v_in /= floor(v_in))) then
            log_Iv = ieee_value(0.d0, ieee_quiet_nan)
            return
        else if (x == 0.0d0) then
            if (v_in == 0.0d0) then
                log_Iv = 0.0d0
            else
                log_Iv = -HUGE_VAL
            end if
            return
        end if
        
        ! Handle negative order (for integer v)
        v = v_in
        if (v < 0.0d0) v = -v  ! I_{-n} = I_{n} for integer n
        
        log_Iv = 0.0d0
        log_x = log(x)
        log_v = log(v)
        
        if (v == 0.0d0) log_v = -1000.0d0
        
        ! Select numerical method based on v and x
        
        ! Method 1: μ order 3
        
        if ( ( (x > 1.4d3) .and. (v < 3.05d0) ) .or. &                     ! orig
           ( ((0.6229d0*log_x - 3.2318d0) > log_v) .and. (v > 3.1d0) &     ! orig
            .and. (v**8 / (384d0*x**4) < 1d-16) ) ) then                   ! new
            
            mu = 4.0d0 * v * v
            curr_term = 1.0d0
            inner_sum = curr_term
            
            k = 1
            c = 1

            do while (k <= 5)
                curr_term = curr_term * (-(mu - dble(k*k)) / (dble(c) * 8.0d0 * x))
                inner_sum = inner_sum + curr_term
                k = k + 2
                c = c + 1
            end do

            
            log_Iv = LOG_INV_SQRT_2PI - 0.5d0*log_x + log(abs(inner_sum))


        ! Method 2: μ order 20
        else if ( ((x > 30.0d0) .and. (v < 15.3919d0)) .or. &
          ((0.5113d0*log_x + 0.7939d0 > log_v) .and. (x > 59.6925d0) &
           .and. (v < dsqrt(3.0d0*x)) ) ) then
            
            mu = 4.0d0 * v * v
            curr_term = 1.0d0
            inner_sum = curr_term
            
            k = 1
            c = 1
            do while (k <= 39)
                curr_term = curr_term * (-(mu - dble(k*k)) / (dble(c) * 8.0d0 * x))
                inner_sum = inner_sum + curr_term
                k = k + 2
                c = c + 1
            end do
            
            log_Iv =  LOG_INV_SQRT_2PI - 0.5d0*log_x + log(abs(inner_sum))
        
        ! Method 3: uk order 4
        else if ((x > 274.2377d0 .and. v > 0.3d0) .or. (v > 163.6993d0)) then
            
            !----- stable  v*(eta - x')  --------------------------------
            u = v / x
            if (u < 0.5d0) then
                asinh_u = asinh(u)
                t1      = u / (sqrt(1.0d0 + u*u) + 1.0d0)
                v_eta   = v * (t1 - asinh_u)
                x_prime   = 1.0d0 / u
                x_prime_2 = 1.0d0 / (u*u)
            else
                x_prime   = 1.0d0 / u
                x_prime_2 = x_prime * x_prime
                s         = sqrt(1.0d0 + x_prime_2)
                v_eta     = v * (s - x_prime + log(x_prime/(1.0d0 + s)))
            end if
            sqrt_1_plus_x2 = sqrt(1.0d0 + x_prime_2)
            t              = 1.0d0 / sqrt_1_plus_x2
            !-----------------------------------------------------------
            
            v_power(0) = 1.0d0
            t_power(0) = 1.0d0
            
            do i = 1, 4
                v_power(i) = v * v_power(i-1)
            end do
            
            do i = 1, 4*3
                t_power(i) = t * t_power(i-1)
            end do

            log_Iv = LOG_INV_SQRT_2PI - 0.5d0*log(v) + v_eta - 0.25d0*log(1.0d0 + x_prime_2) + &
                log(abs(1.0d0 + u1(t_power)/v_power(1) + u2(t_power)/v_power(2) + &
                            u3(t_power)/v_power(3) + u4(t_power)/v_power(4)))
        
        ! Method 4: uk order 6
        else if ((x > 84.4153d0 .and. v > 0.46d0) .or. (v > 56.9971d0)) then
            
            !----- stable  v*(eta - x')  --------------------------------
            u = v / x
            if (u < 0.5d0) then
                asinh_u = asinh(u)
                t1      = u / (sqrt(1.0d0 + u*u) + 1.0d0)
                v_eta   = v * (t1 - asinh_u)
                x_prime   = 1.0d0 / u
                x_prime_2 = 1.0d0 / (u*u)
            else
                x_prime   = 1.0d0 / u
                x_prime_2 = x_prime * x_prime
                s         = sqrt(1.0d0 + x_prime_2)
                v_eta     = v * (s - x_prime + log(x_prime/(1.0d0 + s)))
            end if
            sqrt_1_plus_x2 = sqrt(1.0d0 + x_prime_2)
            t              = 1.0d0 / sqrt_1_plus_x2
            !-----------------------------------------------------------
            
            v_power(0) = 1.0d0
            t_power(0) = 1.0d0
            
            do i = 1, 6
                v_power(i) = v * v_power(i-1)
            end do
            
            do i = 1, 6*3
                t_power(i) = t * t_power(i-1)
            end do

            log_Iv = LOG_INV_SQRT_2PI - 0.5d0*log(v) + v_eta - 0.25d0*log(1.0d0 + x_prime_2) + &
                log(abs(1.0d0 + u1(t_power)/v_power(1) + u2(t_power)/v_power(2) + &
                            u3(t_power)/v_power(3) + u4(t_power)/v_power(4) + &
                            u5(t_power)/v_power(5) + u6(t_power)/v_power(6)))
        
        ! Method 5: uk order 9
        else if ((x > 35.9074d0 .and. v > 0.6d0) .or. (v > 20.1534d0)) then
            
            !----- stable  v*(eta - x')  --------------------------------
            u = v / x
            if (u < 0.5d0) then
                asinh_u = asinh(u)
                t1      = u / (sqrt(1.0d0 + u*u) + 1.0d0)
                v_eta   = v * (t1 - asinh_u)
                x_prime   = 1.0d0 / u
                x_prime_2 = 1.0d0 / (u*u)
            else
                x_prime   = 1.0d0 / u
                x_prime_2 = x_prime * x_prime
                s         = sqrt(1.0d0 + x_prime_2)
                v_eta     = v * (s - x_prime + log(x_prime/(1.0d0 + s)))
            end if
            sqrt_1_plus_x2 = sqrt(1.0d0 + x_prime_2)
            t              = 1.0d0 / sqrt_1_plus_x2
            !-----------------------------------------------------------
            
            v_power(0) = 1.0d0
            t_power(0) = 1.0d0
            
            do i = 1, 9
                v_power(i) = v * v_power(i-1)
            end do
            
            do i = 1, 9*3
                t_power(i) = t * t_power(i-1)
            end do
            
            log_Iv = LOG_INV_SQRT_2PI - 0.5d0*log(v) + v_eta - 0.25d0*log(1.0d0 + x_prime_2) + &
                log(abs(1.0d0 + u1(t_power)/v_power(1) + u2(t_power)/v_power(2) + &
                            u3(t_power)/v_power(3) + u4(t_power)/v_power(4) + &
                            u5(t_power)/v_power(5) + u6(t_power)/v_power(6) + &
                            u7(t_power)/v_power(7) + u8(t_power)/v_power(8) + &
                            u9(t_power)/v_power(9)))
        
        ! Method 6: uk order 13
        else if ((x > 19.6931d0 .and. v > 0.7d0) .or. (v > 12.6964d0)) then
            
            !----- stable  v*(eta - x')  --------------------------------
            u = v / x
            if (u < 0.5d0) then
                asinh_u = asinh(u)
                t1      = u / (sqrt(1.0d0 + u*u) + 1.0d0)
                v_eta   = v * (t1 - asinh_u)
                x_prime   = 1.0d0 / u
                x_prime_2 = 1.0d0 / (u*u)
            else
                x_prime   = 1.0d0 / u
                x_prime_2 = x_prime * x_prime
                s         = sqrt(1.0d0 + x_prime_2)
                v_eta     = v * (s - x_prime + log(x_prime/(1.0d0 + s)))
            end if
            sqrt_1_plus_x2 = sqrt(1.0d0 + x_prime_2)
            t              = 1.0d0 / sqrt_1_plus_x2
            !-----------------------------------------------------------
            
            v_power(0) = 1.0d0
            t_power(0) = 1.0d0
            
            do i = 1, 13
                v_power(i) = v * v_power(i-1)
            end do
            
            do i = 1, 13*3
                t_power(i) = t * t_power(i-1)
            end do
            
            log_Iv = LOG_INV_SQRT_2PI - 0.5d0*log(v) + v_eta - 0.25d0*log(1.0d0 + x_prime_2) + &
                log(abs(1.0d0 + u1(t_power)/v_power(1) + u2(t_power)/v_power(2) + &
                            u3(t_power)/v_power(3) + u4(t_power)/v_power(4) + &
                            u5(t_power)/v_power(5) + u6(t_power)/v_power(6) + &
                            u7(t_power)/v_power(7) + u8(t_power)/v_power(8) + &
                            u9(t_power)/v_power(9) + u10(t_power)/v_power(10) + &
                            u11(t_power)/v_power(11) + u12(t_power)/v_power(12) + &
                            u13(t_power)/v_power(13)))
        

        ! Method 7: Sum definition (default)
        else
            
            log_x2_over_4 = 2.0d0 * log_x - 1.38629436111989062d0
            x2_over_4 = x * x / 4.0d0
            inv_x2_over_4 = 1.0d0 / x2_over_4
            v_inv_x2_over_4 = v * inv_x2_over_4
            
            terms(0) = -log_gamma(v + 1.0d0)
            
            do k = 1, N_TERMS-1
                terms(k) = terms(k-1) - log(dble(k) * (dble(k)*inv_x2_over_4 + v_inv_x2_over_4))
            end do
            
            peak_k = floor((-v + dsqrt(v*v + x*x)) / 2.0d0)
            
            if (peak_k < 0) peak_k = 0
            
            sum_terms = 0.0d0
            do k = 0, N_TERMS-1
                sum_terms = sum_terms + exp(terms(k) - terms(peak_k))
            end do
            
            log_Iv = v * (log_x - LOG_2) + terms(peak_k) + log(sum_terms) - x 
            
        end if
 
      end function iv_log
    
      pure elemental function iv_scaled(v_in, x) result(Iv)

        double precision, intent(in) :: x
        integer, intent(in)          :: v_in
        double precision             :: Iv
        
          Iv = dexp(iv_log(dble(v_in),x))
          
      end function iv_scaled

      pure elemental function iv_log_scaled(v_in, x) result(Iv)

        double precision, intent(in) :: x
        integer, intent(in)          :: v_in
        double precision             :: Iv
        
          Iv = iv_log(dble(v_in),x)
          
      end function iv_log_scaled

end module bessel_functions
