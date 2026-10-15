subroutine integrate_ERI_sum(pattern_id,p,q,p_x,q_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq &
     &                                                                 ,ya,yb,yc,yd,yp,yq &
     &                                                                 ,za,zb,zc,zd,zp,zq,result)
      
      use quadpack , only : dqags
      use iso_c_binding
      use torus_init
      use bessel_functions
      use files
      use table_lookup_module
      use table_1d_lookup 
      !use gauss_legendre_quadrature 
      use precomputed_bessel

      use, intrinsic :: ieee_arithmetic

      implicit none

      ! Input parameters
      double precision, intent(in)       :: p_x , p
      double precision, intent(in)       :: q_x , q
      double precision, intent(in)       :: phi
      double precision, intent(in)       :: xpA , xpB 
      double precision, intent(in)       :: xqC , xqD
      double precision, intent(in)       :: xa, xb, xc, xd, xp, xq 
      double precision, intent(in)       :: ya, yb, yc, yd, yp, yq
      double precision, intent(in)       :: za, zb, zc, zd, zp, zq
      integer         , intent(in)       :: pattern_id
      

      ! Output parameters

      double precision, intent(out)      :: result
    
      ! Local variables

      double precision,parameter         :: epsabs = 1.0d-10 , epsrel = 1.0d-8
      integer,parameter                  :: inf = 1 
      double precision,parameter         :: bound = 0.0d0
      integer, parameter                 :: limit = 50
      integer, parameter                 :: lenw = limit*4
      integer                            :: ier, iwork(limit), last, neval
      double precision                   :: abserr, work(lenw)
      
      integer                            :: Nmax , Nmax_A , Nmax_B
      double precision,parameter         :: pi2  = pi * pi

      double precision                   :: inv_ax 
      double precision                   :: inv_ax2

      double precision                   :: cxpa   , sxpa
      double precision                   :: cxpb   , sxpb
      double precision                   :: cxqc   , sxqc
      double precision                   :: cxqd   , sxqd
      double precision                   :: c2xpab , s2xpab
      double precision                   :: c2xqcd , s2xqcd
      double precision                   :: ypq    , zpq 
      double precision                   :: ypa    , ypb 
      double precision                   :: yqc    , yqd 
      double precision                   :: zpa    , zpb 
      double precision                   :: zqc    , zqd 

      double precision                   :: A , B , C

      


      ! --------------------------------------------------------------- !
      ! backward ! 

      double precision                     :: I_A(0:50000)
      double precision                     :: I_B(0:50000)

      COMPLEX(KIND=KIND(1.0D0)), PARAMETER :: I_dp = (0.0D0, 1.0D0)
      integer                              :: i 
      COMPLEX(KIND=KIND(1.0D0))            :: exp_term(0:50000)

      logical                              :: no_converged_x = .true.

      ! --------------------------------------------------------------- !

      ! --------------------------------------------------------------- !

      inv_ax  = 1.d0/ax 
      inv_ax2 = inv_ax * inv_ax 

      cxpa  = dcos(xpA)
      sxpa  = dsin(xpA)

      cxpb  = dcos(xpB)
      sxpb  = dsin(xpB)
      
      cxqc  = dcos(xqC)
      sxqc  = dsin(xqC)

      cxqd  = dcos(xqD)
      sxqd  = dsin(xqD)

      c2xpab = dcos(ax*(2.d0*xp-xA-xB))
      s2xpab = dsin(ax*(2.d0*xp-xA-xB))

      c2xqcd = dcos(ax*(2.d0*xq-xC-xD))
      s2xqcd = dsin(ax*(2.d0*xq-xC-xD))


      ! /////////////////////////////////////////////////////////////// !
      ! Real Gaussian ! 

      ypq    = yp - yq 

      ypa    = yp - ya 
      ypb    = yp - yb 

      yqc    = yq - yc
      yqd    = yq - yd

      zpq    = zp - zq

      zpa    = zp - za 
      zpb    = zp - zb 

      zqc    = zq - zc
      zqd    = zq - zd

      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !

      
      A   = 2.d0  * p_x   * inv_ax2
      B   = 2.d0  * q_x   * inv_ax2

      Nmax_A = get_Nmax_from_bessel_1D(A)
      Nmax_B = get_Nmax_from_bessel_1D(B)
 
      call bessel_I_scaled_backward(Nmax_A+20, A, I_A)
      call bessel_I_scaled_backward(Nmax_B+20, B, I_B)

      do i = 1 , Nmax_A
        exp_term(i)    = exp(I_dp * i * phi)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=--=-=-=-=-=-=-=-=-=-=- !
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=--=-=-=-=-=-=-=-=-=-=- !

      !call dqags(transformed_integrand, 0.0d0, 1.0d0, epsabs, epsrel, &
      !          result, abserr, neval, ier, limit, lenw, last, iwork, work)

      !if (ier > 2) then
      !  write(*,'(A,I8,A)') 'Error code = ', ier
      !  write(outfile,'(A,I8,A)') 'you have an integral did not converge '
      !  !stop
      !end if

      !call gauss_legendre_auto(transformed_integrand, 0.0d0, 1.0d0, result, 2)
      
      !call gauss_legendre_64(transformed_integrand, 0.0d0, 1.0d0, result)

      contains
 
      function f_decay(t,i_quad) result(ft)
        double precision, intent(in) :: t
        integer, intent(in)          :: i_quad
        double precision             :: ft
         ft =  S(t,i_quad)
      end function f_decay

      ! Transformed integrand for u in [0,1)

      function transformed_integrand(u,i_quad) result(y)
        use precomputed_bessel
        double precision, intent(in) :: u
        integer, intent(in) :: i_quad
        double precision :: y, t
        
        if (u >= 1.0d0) then
          y = 0.0d0
          return
        endif
        
        ! Transformation: t = u/(1-u)
        t = u / (1.0d0 - u)
        
        ! Jacobian: dt/du = 1/(1-u)²
        y = S(t,i_quad) / ((1.0d0 - u) * (1.0d0 - u))
        
        ! Check for numerical issues
        if (.not. ieee_is_finite(y)) then
          y = 0.0d0
        endif

      end function transformed_integrand


      double precision function S(t,i_quad) result(sum)

      use bessel_derivatives
      use precomputed_bessel
      use functions

      implicit none

      double precision                     :: t
      integer, intent(in)                  :: i_quad
      double precision                     :: D, D2 , piD , piD2 
      double precision                     :: const
      integer                              :: n     
      COMPLEX(KIND=KIND(1.0D0))            :: termAn , termBn , term
      COMPLEX(KIND=KIND(1.0D0)), PARAMETER :: I_dp = (0.0D0, 1.0D0)
      COMPLEX(KIND=KIND(1.0D0))            :: expo_term , current_term

      
      double precision                     :: t2 , t3 , t4 
      double precision,parameter           :: eps = 2.22d-16


      double precision                     :: eta_t   , eta_p  , eta_q

      ! /////////////////////////////////////////////////////////////// !

      double precision                     :: Ireal_0_ppt , Ireal_0_upq
      double precision                     :: u

      ! /////////////////////////////////////////////////////////////// !


      t2 = t  * t 
      t3 = t2 * t 
      t4 = t3 * t 

      C    = 2.d0  * t * t * inv_ax2

      D    = 1.d0/(dsqrt(p*q+(p+q)*t2))

      D2   = D * D

      piD  = pi * D
      piD2 = piD * piD

      expo_term    = exp(I_dp * phi)
      current_term = expo_term

      ! /////////////////////////////////////////////////////////////// !

      eta_t  =  t2 / (p+t2)    
      eta_p  =  p*t2  * D2    
      eta_q  =  q*t2  * D2    
      
      u      = eta_t * p

      Ireal_0_ppt = Ireal(0,p+t2)
      Ireal_0_upq = Ireal(0,u+q )

      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !

      term    = 0.d0
      termAn  = 0.d0
      termBn  = 0.d0

      Nmax         = get_Nmax(A, B, C, Phi) + 5

      select case(pattern_id)
      
      ! /////////////////////////////////////////////////////////////// !

      case (0000) ! | s     s    s    s     ( 1 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do
          
      case (0001) ! | s     s    s    px    ( 2 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0002) ! | s     s    s    py    ( 3 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0003) ! | s     s    s    pz    ( 4 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0010) ! | s     s    px   s     ( 5 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0011) ! | s     s    px   px    ( 6 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0012) ! | s     s    px   py    ( 7 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0013) ! | s     s    px   pz    ( 8 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0020) ! | s     s    py   s     ( 9 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0021) ! | s     s    py   px    ( 10 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0022) ! | s     s    py   py    ( 11 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*ypq**2 + eta_p*ypq*(yqc + yqd) + yqc*yqd))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0023) ! | s     s    py   pz    ( 12 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)*(eta_p*zpq + zqd)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0030) ! | s     s    pz   s     ( 13 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0031) ! | s     s    pz   px    ( 14 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0032) ! | s     s    pz   py    ( 15 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)*(eta_p*zpq + zqc)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0033) ! | s     s    pz   pz    ( 16 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*zpq**2 + eta_p*zpq*(zqc + zqd) + zqc*zqd))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0100) ! | s     px   s    s     ( 17 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0101) ! | s     px   s    px    ( 18 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax**2*(cxpb*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0102) ! | s     px   s    py    ( 19 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0103) ! | s     px   s    pz    ( 20 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0110) ! | s     px   px   s     ( 21 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax**2*(cxpb*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0111) ! | s     px   px   px    ( 22 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax*inv_ax2*(cxpb*cxqc*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxpb*cxqc*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpb*cxqd*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpb*sxqc*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 2, n, I_B, B) + cxqc*cxqd*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxqc*sxpb*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxqd*sxpb*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + sxpb*sxqc*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 2, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0112) ! | s     px   px   py    ( 23 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)) 
      sum    = (inv_ax**2*(cxpb*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0113) ! | s     px   px   pz    ( 24 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)) 
      sum    = (inv_ax**2*(cxpb*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0120) ! | s     px   py   s     ( 25 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0121) ! | s     px   py   px    ( 26 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)) 
      sum    = (inv_ax**2*(cxpb*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0122) ! | s     px   py   py    ( 27 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*ypq**2 + eta_p*ypq*(yqc + yqd) + yqc*yqd))) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0123) ! | s     px   py   pz    ( 28 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)*(eta_p*zpq + zqd)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0130) ! | s     px   pz   s     ( 29 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0131) ! | s     px   pz   px    ( 30 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)) 
      sum    = (inv_ax**2*(cxpb*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0132) ! | s     px   pz   py    ( 31 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)*(eta_p*zpq + zqc)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0133) ! | s     px   pz   pz    ( 32 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*zpq**2 + eta_p*zpq*(zqc + zqd) + zqc*zqd))) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0200) ! | s     py   s    s     ( 33 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypb)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0201) ! | s     py   s    px    ( 34 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0202) ! | s     py   s    py    ( 35 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqd*(-eta_q*ypq + ypb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0203) ! | s     py   s    pz    ( 36 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)*(-eta_q*ypq + ypb)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0210) ! | s     py   px   s     ( 37 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0211) ! | s     py   px   px    ( 38 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0212) ! | s     py   px   py    ( 39 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqd*(-eta_q*ypq + ypb)))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0213) ! | s     py   px   pz    ( 40 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0220) ! | s     py   py   s     ( 41 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqc*(-eta_q*ypq + ypb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0221) ! | s     py   py   px    ( 42 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqc*(-eta_q*ypq + ypb)))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0222) ! | s     py   py   py    ( 43 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(-eta_q*ypq + eta_t*(2*eta_p*ypq + yqc + yqd) + ypb)*Ireal(2, q + u) + piD*(eta_p**2*(-eta_q*ypq**3 + ypb*ypq**2) + eta_p*(eta_q*ypq**2*(-yqc - yqd) + ypb*ypq*(yqc + yqd)) + yqc*yqd*(-eta_q*ypq + ypb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0223) ! | s     py   py   pz    ( 44 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqd)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqc*(-eta_q*ypq + ypb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0230) ! | s     py   pz   s     ( 45 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)*(-eta_q*ypq + ypb)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0231) ! | s     py   pz   px    ( 46 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0232) ! | s     py   pz   py    ( 47 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqc)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqd*(-eta_q*ypq + ypb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0233) ! | s     py   pz   pz    ( 48 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*zpq**2 + eta_p*zpq*(zqc + zqd) + zqc*zqd))*(-eta_q*ypq + ypb)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0300) ! | s     pz   s    s     ( 49 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpb)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0301) ! | s     pz   s    px    ( 50 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0302) ! | s     pz   s    py    ( 51 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)*(-eta_q*zpq + zpb)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0303) ! | s     pz   s    pz    ( 52 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqd*(-eta_q*zpq + zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0310) ! | s     pz   px   s     ( 53 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0311) ! | s     pz   px   px    ( 54 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0312) ! | s     pz   px   py    ( 55 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0313) ! | s     pz   px   pz    ( 56 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqd*(-eta_q*zpq + zpb)))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0320) ! | s     pz   py   s     ( 57 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)*(-eta_q*zpq + zpb)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0321) ! | s     pz   py   px    ( 58 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0322) ! | s     pz   py   py    ( 59 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*ypq**2 + eta_p*ypq*(yqc + yqd) + yqc*yqd))*(-eta_q*zpq + zpb)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0323) ! | s     pz   py   pz    ( 60 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqc)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqd*(-eta_q*zpq + zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0330) ! | s     pz   pz   s     ( 61 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqc*(-eta_q*zpq + zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0331) ! | s     pz   pz   px    ( 62 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqc*(-eta_q*zpq + zpb)))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0332) ! | s     pz   pz   py    ( 63 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqd)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqc*(-eta_q*zpq + zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (0333) ! | s     pz   pz   pz    ( 64 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(-eta_q*zpq + eta_t*(2*eta_p*zpq + zqc + zqd) + zpb)*Ireal(2, q + u) + piD*(eta_p**2*(-eta_q*zpq**3 + zpb*zpq**2) + eta_p*(eta_q*zpq**2*(-zqc - zqd) + zpb*zpq*(zqc + zqd)) + zqc*zqd*(-eta_q*zpq + zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1000) ! | px    s    s    s     ( 65 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1001) ! | px    s    s    px    ( 66 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax**2*(cxpa*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1002) ! | px    s    s    py    ( 67 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1003) ! | px    s    s    pz    ( 68 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1010) ! | px    s    px   s     ( 69 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax**2*(cxpa*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1011) ! | px    s    px   px    ( 70 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax*inv_ax2*(cxpa*cxqc*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxpa*cxqc*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpa*cxqd*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpa*sxqc*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 2, n, I_B, B) + cxqc*cxqd*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxqc*sxpa*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxqd*sxpa*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + sxpa*sxqc*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 2, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1012) ! | px    s    px   py    ( 71 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)) 
      sum    = (inv_ax**2*(cxpa*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1013) ! | px    s    px   pz    ( 72 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)) 
      sum    = (inv_ax**2*(cxpa*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1020) ! | px    s    py   s     ( 73 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1021) ! | px    s    py   px    ( 74 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)) 
      sum    = (inv_ax**2*(cxpa*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1022) ! | px    s    py   py    ( 75 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*ypq**2 + eta_p*ypq*(yqc + yqd) + yqc*yqd))) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1023) ! | px    s    py   pz    ( 76 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)*(eta_p*zpq + zqd)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1030) ! | px    s    pz   s     ( 77 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1031) ! | px    s    pz   px    ( 78 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)) 
      sum    = (inv_ax**2*(cxpa*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1032) ! | px    s    pz   py    ( 79 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)*(eta_p*zpq + zqc)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1033) ! | px    s    pz   pz    ( 80 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*zpq**2 + eta_p*zpq*(zqc + zqd) + zqc*zqd))) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1100) ! | px    px   s    s     ( 81 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1101) ! | px    px   s    px    ( 82 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax*inv_ax2*(cxpa*cxpb*cxqd*der_I_A(2, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*cxpb*sxqd*der_I_A(2, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxpa*cxqd*sxpb*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxpb*sxqd*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxpb*cxqd*sxpa*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxpa*sxqd*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpa*sxpb*der_I_A(0, 2, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxpb*sxqd*der_I_A(0, 2, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1102) ! | px    px   s    py    ( 83 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)) 
      sum    = (inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1103) ! | px    px   s    pz    ( 84 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)) 
      sum    = (inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1110) ! | px    px   px   s     ( 85 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax*inv_ax2*(cxpa*cxpb*cxqc*der_I_A(2, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*cxpb*sxqc*der_I_A(2, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxpa*cxqc*sxpb*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxpb*sxqc*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxpb*cxqc*sxpa*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxpa*sxqc*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpa*sxpb*der_I_A(0, 2, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxpb*sxqc*der_I_A(0, 2, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1111) ! | px    px   px   px    ( 86 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2) 
      sum    = (inv_ax2**2*(cxpa*cxpb*cxqc*cxqd*der_I_A(2, 0, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxpa*cxpb*cxqc*sxqd*der_I_A(2, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpa*cxpb*cxqd*sxqc*der_I_A(2, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpa*cxpb*sxqc*sxqd*der_I_A(2, 0, n, I_A, A)*der_I_B(0, 2, n, I_B, B) + cxpa*cxqc*cxqd*sxpb*der_I_A(1, 1, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxpa*cxqc*sxpb*sxqd*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpa*cxqd*sxpb*sxqc*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpa*sxpb*sxqc*sxqd*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 2, n, I_B, B) + cxpb*cxqc*cxqd*sxpa*der_I_A(1, 1, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxpb*cxqc*sxpa*sxqd*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpb*cxqd*sxpa*sxqc*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpb*sxpa*sxqc*sxqd*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 2, n, I_B, B) + cxqc*cxqd*sxpa*sxpb*der_I_A(0, 2, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxqc*sxpa*sxpb*sxqd*der_I_A(0, 2, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxqd*sxpa*sxpb*sxqc*der_I_A(0, 2, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + sxpa*sxpb*sxqc*sxqd*der_I_A(0, 2, n, I_A, A)*der_I_B(0, 2, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1112) ! | px    px   px   py    ( 87 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)) 
      sum    = (inv_ax*inv_ax2*(cxpa*cxpb*cxqc*der_I_A(2, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*cxpb*sxqc*der_I_A(2, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxpa*cxqc*sxpb*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxpb*sxqc*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxpb*cxqc*sxpa*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxpa*sxqc*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpa*sxpb*der_I_A(0, 2, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxpb*sxqc*der_I_A(0, 2, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1113) ! | px    px   px   pz    ( 88 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)) 
      sum    = (inv_ax*inv_ax2*(cxpa*cxpb*cxqc*der_I_A(2, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*cxpb*sxqc*der_I_A(2, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxpa*cxqc*sxpb*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxpb*sxqc*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxpb*cxqc*sxpa*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxpa*sxqc*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpa*sxpb*der_I_A(0, 2, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxpb*sxqc*der_I_A(0, 2, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1120) ! | px    px   py   s     ( 89 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)) 
      sum    = (inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1121) ! | px    px   py   px    ( 90 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)) 
      sum    = (inv_ax*inv_ax2*(cxpa*cxpb*cxqd*der_I_A(2, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*cxpb*sxqd*der_I_A(2, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxpa*cxqd*sxpb*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxpb*sxqd*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxpb*cxqd*sxpa*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxpa*sxqd*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpa*sxpb*der_I_A(0, 2, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxpb*sxqd*der_I_A(0, 2, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1122) ! | px    px   py   py    ( 91 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*ypq**2 + eta_p*ypq*(yqc + yqd) + yqc*yqd))) 
      sum    = (inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1123) ! | px    px   py   pz    ( 92 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)*(eta_p*zpq + zqd)) 
      sum    = (inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1130) ! | px    px   pz   s     ( 93 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)) 
      sum    = (inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1131) ! | px    px   pz   px    ( 94 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)) 
      sum    = (inv_ax*inv_ax2*(cxpa*cxpb*cxqd*der_I_A(2, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*cxpb*sxqd*der_I_A(2, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxpa*cxqd*sxpb*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxpb*sxqd*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxpb*cxqd*sxpa*der_I_A(1, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxpa*sxqd*der_I_A(1, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpa*sxpb*der_I_A(0, 2, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxpb*sxqd*der_I_A(0, 2, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1132) ! | px    px   pz   py    ( 95 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)*(eta_p*zpq + zqc)) 
      sum    = (inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1133) ! | px    px   pz   pz    ( 96 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*zpq**2 + eta_p*zpq*(zqc + zqd) + zqc*zqd))) 
      sum    = (inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax2*(cxpa*cxpb*der_I_A(2, 0, n, I_A, A) + cxpa*sxpb*der_I_A(1, 1, n, I_A, A) + cxpb*sxpa*der_I_A(1, 1, n, I_A, A) + sxpa*sxpb*der_I_A(0, 2, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1200) ! | px    py   s    s     ( 97 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1201) ! | px    py   s    px    ( 98 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax**2*(cxpa*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1202) ! | px    py   s    py    ( 99 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqd*(-eta_q*ypq + ypb)))) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1203) ! | px    py   s    pz    ( 100 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1210) ! | px    py   px   s     ( 101 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax**2*(cxpa*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1211) ! | px    py   px   px    ( 102 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax*inv_ax2*(cxpa*cxqc*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxpa*cxqc*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpa*cxqd*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpa*sxqc*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 2, n, I_B, B) + cxqc*cxqd*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxqc*sxpa*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxqd*sxpa*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + sxpa*sxqc*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 2, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1212) ! | px    py   px   py    ( 103 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqd*(-eta_q*ypq + ypb)))) 
      sum    = (inv_ax**2*(cxpa*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1213) ! | px    py   px   pz    ( 104 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax**2*(cxpa*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1220) ! | px    py   py   s     ( 105 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqc*(-eta_q*ypq + ypb)))) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1221) ! | px    py   py   px    ( 106 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqc*(-eta_q*ypq + ypb)))) 
      sum    = (inv_ax**2*(cxpa*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1222) ! | px    py   py   py    ( 107 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(-eta_q*ypq + eta_t*(2*eta_p*ypq + yqc + yqd) + ypb)*Ireal(2, q + u) + piD*(eta_p**2*(-eta_q*ypq**3 + ypb*ypq**2) + eta_p*(eta_q*ypq**2*(-yqc - yqd) + ypb*ypq*(yqc + yqd)) + yqc*yqd*(-eta_q*ypq + ypb)))) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1223) ! | px    py   py   pz    ( 108 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqd)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqc*(-eta_q*ypq + ypb)))) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1230) ! | px    py   pz   s     ( 109 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1231) ! | px    py   pz   px    ( 110 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax**2*(cxpa*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1232) ! | px    py   pz   py    ( 111 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqc)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqd*(-eta_q*ypq + ypb)))) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1233) ! | px    py   pz   pz    ( 112 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*zpq**2 + eta_p*zpq*(zqc + zqd) + zqc*zqd))*(-eta_q*ypq + ypb)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1300) ! | px    pz   s    s     ( 113 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1301) ! | px    pz   s    px    ( 114 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax**2*(cxpa*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1302) ! | px    pz   s    py    ( 115 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1303) ! | px    pz   s    pz    ( 116 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqd*(-eta_q*zpq + zpb)))) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1310) ! | px    pz   px   s     ( 117 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax**2*(cxpa*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1311) ! | px    pz   px   px    ( 118 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax*inv_ax2*(cxpa*cxqc*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxpa*cxqc*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpa*cxqd*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpa*sxqc*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 2, n, I_B, B) + cxqc*cxqd*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxqc*sxpa*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxqd*sxpa*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + sxpa*sxqc*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 2, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1312) ! | px    pz   px   py    ( 119 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax**2*(cxpa*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1313) ! | px    pz   px   pz    ( 120 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqd*(-eta_q*zpq + zpb)))) 
      sum    = (inv_ax**2*(cxpa*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1320) ! | px    pz   py   s     ( 121 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1321) ! | px    pz   py   px    ( 122 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax**2*(cxpa*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1322) ! | px    pz   py   py    ( 123 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*ypq**2 + eta_p*ypq*(yqc + yqd) + yqc*yqd))*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1323) ! | px    pz   py   pz    ( 124 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqc)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqd*(-eta_q*zpq + zpb)))) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1330) ! | px    pz   pz   s     ( 125 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqc*(-eta_q*zpq + zpb)))) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1331) ! | px    pz   pz   px    ( 126 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqc*(-eta_q*zpq + zpb)))) 
      sum    = (inv_ax**2*(cxpa*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpa*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpa*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpa*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1332) ! | px    pz   pz   py    ( 127 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqd)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqc*(-eta_q*zpq + zpb)))) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (1333) ! | px    pz   pz   pz    ( 128 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(-eta_q*zpq + eta_t*(2*eta_p*zpq + zqc + zqd) + zpb)*Ireal(2, q + u) + piD*(eta_p**2*(-eta_q*zpq**3 + zpb*zpq**2) + eta_p*(eta_q*zpq**2*(-zqc - zqd) + zpb*zpq*(zqc + zqd)) + zqc*zqd*(-eta_q*zpq + zpb)))) 
      sum    = (inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpa*der_I_A(1, 0, n, I_A, A) + sxpa*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2000) ! | py    s    s    s     ( 129 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypa)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2001) ! | py    s    s    px    ( 130 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2002) ! | py    s    s    py    ( 131 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqd*(-eta_q*ypq + ypa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2003) ! | py    s    s    pz    ( 132 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)*(-eta_q*ypq + ypa)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2010) ! | py    s    px   s     ( 133 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2011) ! | py    s    px   px    ( 134 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2012) ! | py    s    px   py    ( 135 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqd*(-eta_q*ypq + ypa)))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2013) ! | py    s    px   pz    ( 136 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2020) ! | py    s    py   s     ( 137 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqc*(-eta_q*ypq + ypa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2021) ! | py    s    py   px    ( 138 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqc*(-eta_q*ypq + ypa)))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2022) ! | py    s    py   py    ( 139 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(-eta_q*ypq + eta_t*(2*eta_p*ypq + yqc + yqd) + ypa)*Ireal(2, q + u) + piD*(eta_p**2*(-eta_q*ypq**3 + ypa*ypq**2) + eta_p*(eta_q*ypq**2*(-yqc - yqd) + ypa*ypq*(yqc + yqd)) + yqc*yqd*(-eta_q*ypq + ypa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2023) ! | py    s    py   pz    ( 140 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqd)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqc*(-eta_q*ypq + ypa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2030) ! | py    s    pz   s     ( 141 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)*(-eta_q*ypq + ypa)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2031) ! | py    s    pz   px    ( 142 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2032) ! | py    s    pz   py    ( 143 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqc)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqd*(-eta_q*ypq + ypa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2033) ! | py    s    pz   pz    ( 144 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*zpq**2 + eta_p*zpq*(zqc + zqd) + zqc*zqd))*(-eta_q*ypq + ypa)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2100) ! | py    px   s    s     ( 145 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2101) ! | py    px   s    px    ( 146 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax**2*(cxpb*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2102) ! | py    px   s    py    ( 147 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqd*(-eta_q*ypq + ypa)))) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2103) ! | py    px   s    pz    ( 148 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2110) ! | py    px   px   s     ( 149 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax**2*(cxpb*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2111) ! | py    px   px   px    ( 150 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax*inv_ax2*(cxpb*cxqc*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxpb*cxqc*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpb*cxqd*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpb*sxqc*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 2, n, I_B, B) + cxqc*cxqd*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxqc*sxpb*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxqd*sxpb*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + sxpb*sxqc*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 2, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2112) ! | py    px   px   py    ( 151 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqd*(-eta_q*ypq + ypa)))) 
      sum    = (inv_ax**2*(cxpb*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2113) ! | py    px   px   pz    ( 152 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqd)*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax**2*(cxpb*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2120) ! | py    px   py   s     ( 153 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqc*(-eta_q*ypq + ypa)))) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2121) ! | py    px   py   px    ( 154 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqc*(-eta_q*ypq + ypa)))) 
      sum    = (inv_ax**2*(cxpb*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2122) ! | py    px   py   py    ( 155 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(-eta_q*ypq + eta_t*(2*eta_p*ypq + yqc + yqd) + ypa)*Ireal(2, q + u) + piD*(eta_p**2*(-eta_q*ypq**3 + ypa*ypq**2) + eta_p*(eta_q*ypq**2*(-yqc - yqd) + ypa*ypq*(yqc + yqd)) + yqc*yqd*(-eta_q*ypq + ypa)))) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2123) ! | py    px   py   pz    ( 156 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqd)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqc*(-eta_q*ypq + ypa)))) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2130) ! | py    px   pz   s     ( 157 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2131) ! | py    px   pz   px    ( 158 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*zpq + zqc)*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax**2*(cxpb*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2132) ! | py    px   pz   py    ( 159 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqc)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqd*(-eta_q*ypq + ypa)))) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2133) ! | py    px   pz   pz    ( 160 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*zpq**2 + eta_p*zpq*(zqc + zqd) + zqc*zqd))*(-eta_q*ypq + ypa)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2200) ! | py    py   s    s     ( 161 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2201) ! | py    py   s    px    ( 162 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2202) ! | py    py   s    py    ( 163 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(eta_t**2*(eta_p*ypq + yqd) + eta_t*(-2*eta_q*ypq + ypa + ypb))*Ireal(2, q + u) + Ireal_0_upq*(eta_p*ypq + yqd)*Ireal(2, p + t2) + piD*(eta_p*(eta_q**2*ypq**3 + eta_q*ypq**2*(-ypa - ypb) + ypa*ypb*ypq) + yqd*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2203) ! | py    py   s    pz    ( 164 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqd)*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2210) ! | py    py   px   s     ( 165 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2211) ! | py    py   px   px    ( 166 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb))) 
      sum    = (inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2212) ! | py    py   px   py    ( 167 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(eta_t**2*(eta_p*ypq + yqd) + eta_t*(-2*eta_q*ypq + ypa + ypb))*Ireal(2, q + u) + Ireal_0_upq*(eta_p*ypq + yqd)*Ireal(2, p + t2) + piD*(eta_p*(eta_q**2*ypq**3 + eta_q*ypq**2*(-ypa - ypb) + ypa*ypb*ypq) + yqd*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb)))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2213) ! | py    py   px   pz    ( 168 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqd)*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2220) ! | py    py   py   s     ( 169 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(eta_t**2*(eta_p*ypq + yqc) + eta_t*(-2*eta_q*ypq + ypa + ypb))*Ireal(2, q + u) + Ireal_0_upq*(eta_p*ypq + yqc)*Ireal(2, p + t2) + piD*(eta_p*(eta_q**2*ypq**3 + eta_q*ypq**2*(-ypa - ypb) + ypa*ypb*ypq) + yqc*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2221) ! | py    py   py   px    ( 170 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(eta_t**2*(eta_p*ypq + yqc) + eta_t*(-2*eta_q*ypq + ypa + ypb))*Ireal(2, q + u) + Ireal_0_upq*(eta_p*ypq + yqc)*Ireal(2, p + t2) + piD*(eta_p*(eta_q**2*ypq**3 + eta_q*ypq**2*(-ypa - ypb) + ypa*ypb*ypq) + yqc*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb)))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2222) ! | py    py   py   py    ( 171 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(eta_t**2*((eta_p**2*ypq**2 + eta_p*ypq*(yqc + yqd) + yqc*yqd)*Ireal(2, q + u) + Ireal(4, q + u)) + (eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + eta_t*(eta_p*(-4*eta_q*ypq**2 + ypq*(2*ypa + 2*ypb)) + eta_q*ypq*(-2*yqc - 2*yqd) + ypa*(yqc + yqd) + ypb*(yqc + yqd)) + ypa*ypb)*Ireal(2, q + u)) + piD*(eta_p**2*(eta_q**2*ypq**4 + eta_q*ypq**3*(-ypa - ypb) + ypa*ypb*ypq**2) + eta_p*(eta_q**2*ypq**3*(yqc + yqd) + eta_q*ypq**2*(ypa*(-yqc - yqd) + ypb*(-yqc - yqd)) + ypa*ypb*ypq*(yqc + yqd)) + yqc*yqd*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb)) + (Ireal_0_upq*(eta_p**2*ypq**2 + eta_p*ypq*(yqc + yqd) + yqc*yqd) + Ireal(2, q + u))*Ireal(2, p + t2))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2223) ! | py    py   py   pz    ( 172 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqd)*(Ireal_0_ppt*(eta_t**2*(eta_p*ypq + yqc) + eta_t*(-2*eta_q*ypq + ypa + ypb))*Ireal(2, q + u) + Ireal_0_upq*(eta_p*ypq + yqc)*Ireal(2, p + t2) + piD*(eta_p*(eta_q**2*ypq**3 + eta_q*ypq**2*(-ypa - ypb) + ypa*ypb*ypq) + yqc*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2230) ! | py    py   pz   s     ( 173 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqc)*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2231) ! | py    py   pz   px    ( 174 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqc)*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2232) ! | py    py   pz   py    ( 175 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*zpq + zqc)*(Ireal_0_ppt*(eta_t**2*(eta_p*ypq + yqd) + eta_t*(-2*eta_q*ypq + ypa + ypb))*Ireal(2, q + u) + Ireal_0_upq*(eta_p*ypq + yqd)*Ireal(2, p + t2) + piD*(eta_p*(eta_q**2*ypq**3 + eta_q*ypq**2*(-ypa - ypb) + ypa*ypb*ypq) + yqd*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2233) ! | py    py   pz   pz    ( 176 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * ((Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*zpq**2 + eta_p*zpq*(zqc + zqd) + zqc*zqd))*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*ypq**2 + eta_q*ypq*(-ypa - ypb) + ypa*ypb))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2300) ! | py    pz   s    s     ( 177 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypa)*(-eta_q*zpq + zpb)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2301) ! | py    pz   s    px    ( 178 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypa)*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2302) ! | py    pz   s    py    ( 179 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*zpq + zpb)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqd*(-eta_q*ypq + ypa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2303) ! | py    pz   s    pz    ( 180 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*ypq + ypa)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqd*(-eta_q*zpq + zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2310) ! | py    pz   px   s     ( 181 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypa)*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2311) ! | py    pz   px   px    ( 182 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypa)*(-eta_q*zpq + zpb)) 
      sum    = (inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2312) ! | py    pz   px   py    ( 183 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*zpq + zpb)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqd*(-eta_q*ypq + ypa)))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2313) ! | py    pz   px   pz    ( 184 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*ypq + ypa)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqd*(-eta_q*zpq + zpb)))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2320) ! | py    pz   py   s     ( 185 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*zpq + zpb)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqc*(-eta_q*ypq + ypa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2321) ! | py    pz   py   px    ( 186 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*zpq + zpb)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqc*(-eta_q*ypq + ypa)))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2322) ! | py    pz   py   py    ( 187 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*zpq + zpb)*(Ireal_0_ppt*(-eta_q*ypq + eta_t*(2*eta_p*ypq + yqc + yqd) + ypa)*Ireal(2, q + u) + piD*(eta_p**2*(-eta_q*ypq**3 + ypa*ypq**2) + eta_p*(eta_q*ypq**2*(-yqc - yqd) + ypa*ypq*(yqc + yqd)) + yqc*yqd*(-eta_q*ypq + ypa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2323) ! | py    pz   py   pz    ( 188 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * ((Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqc*(-eta_q*ypq + ypa)))*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqd*(-eta_q*zpq + zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2330) ! | py    pz   pz   s     ( 189 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*ypq + ypa)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqc*(-eta_q*zpq + zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2331) ! | py    pz   pz   px    ( 190 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*ypq + ypa)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqc*(-eta_q*zpq + zpb)))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2332) ! | py    pz   pz   py    ( 191 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * ((Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypa*ypq) + yqd*(-eta_q*ypq + ypa)))*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpb*zpq) + zqc*(-eta_q*zpq + zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (2333) ! | py    pz   pz   pz    ( 192 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*ypq + ypa)*(Ireal_0_ppt*(-eta_q*zpq + eta_t*(2*eta_p*zpq + zqc + zqd) + zpb)*Ireal(2, q + u) + piD*(eta_p**2*(-eta_q*zpq**3 + zpb*zpq**2) + eta_p*(eta_q*zpq**2*(-zqc - zqd) + zpb*zpq*(zqc + zqd)) + zqc*zqd*(-eta_q*zpq + zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3000) ! | pz    s    s    s     ( 193 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpa)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3001) ! | pz    s    s    px    ( 194 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3002) ! | pz    s    s    py    ( 195 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)*(-eta_q*zpq + zpa)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3003) ! | pz    s    s    pz    ( 196 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqd*(-eta_q*zpq + zpa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3010) ! | pz    s    px   s     ( 197 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3011) ! | pz    s    px   px    ( 198 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3012) ! | pz    s    px   py    ( 199 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3013) ! | pz    s    px   pz    ( 200 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqd*(-eta_q*zpq + zpa)))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3020) ! | pz    s    py   s     ( 201 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)*(-eta_q*zpq + zpa)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3021) ! | pz    s    py   px    ( 202 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3022) ! | pz    s    py   py    ( 203 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*ypq**2 + eta_p*ypq*(yqc + yqd) + yqc*yqd))*(-eta_q*zpq + zpa)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3023) ! | pz    s    py   pz    ( 204 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqc)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqd*(-eta_q*zpq + zpa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3030) ! | pz    s    pz   s     ( 205 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqc*(-eta_q*zpq + zpa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3031) ! | pz    s    pz   px    ( 206 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqc*(-eta_q*zpq + zpa)))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3032) ! | pz    s    pz   py    ( 207 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqd)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqc*(-eta_q*zpq + zpa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3033) ! | pz    s    pz   pz    ( 208 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(-eta_q*zpq + eta_t*(2*eta_p*zpq + zqc + zqd) + zpa)*Ireal(2, q + u) + piD*(eta_p**2*(-eta_q*zpq**3 + zpa*zpq**2) + eta_p*(eta_q*zpq**2*(-zqc - zqd) + zpa*zpq*(zqc + zqd)) + zqc*zqd*(-eta_q*zpq + zpa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3100) ! | pz    px   s    s     ( 209 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3101) ! | pz    px   s    px    ( 210 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax**2*(cxpb*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3102) ! | pz    px   s    py    ( 211 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3103) ! | pz    px   s    pz    ( 212 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqd*(-eta_q*zpq + zpa)))) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3110) ! | pz    px   px   s     ( 213 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax**2*(cxpb*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3111) ! | pz    px   px   px    ( 214 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax*inv_ax2*(cxpb*cxqc*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxpb*cxqc*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpb*cxqd*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxpb*sxqc*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 2, n, I_B, B) + cxqc*cxqd*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(2, 0, n, I_B, B) + cxqc*sxpb*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + cxqd*sxpb*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 1, n, I_B, B) + sxpb*sxqc*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 2, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3112) ! | pz    px   px   py    ( 215 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqd)*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax**2*(cxpb*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3113) ! | pz    px   px   pz    ( 216 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqd*(-eta_q*zpq + zpa)))) 
      sum    = (inv_ax**2*(cxpb*cxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqc*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqc*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqc*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3120) ! | pz    px   py   s     ( 217 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3121) ! | pz    px   py   px    ( 218 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(eta_p*ypq + yqc)*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax**2*(cxpb*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3122) ! | pz    px   py   py    ( 219 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*ypq**2 + eta_p*ypq*(yqc + yqd) + yqc*yqd))*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3123) ! | pz    px   py   pz    ( 220 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqc)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqd*(-eta_q*zpq + zpa)))) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3130) ! | pz    px   pz   s     ( 221 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqc*(-eta_q*zpq + zpa)))) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3131) ! | pz    px   pz   px    ( 222 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqc*(-eta_q*zpq + zpa)))) 
      sum    = (inv_ax**2*(cxpb*cxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + cxpb*sxqd*der_I_A(1, 0, n, I_A, A)*der_I_B(0, 1, n, I_B, B) + cxqd*sxpb*der_I_A(0, 1, n, I_A, A)*der_I_B(1, 0, n, I_B, B) + sxpb*sxqd*der_I_A(0, 1, n, I_A, A)*der_I_B(0, 1, n, I_B, B))) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3132) ! | pz    px   pz   py    ( 223 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqd)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqc*(-eta_q*zpq + zpa)))) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3133) ! | pz    px   pz   pz    ( 224 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(-eta_q*zpq + eta_t*(2*eta_p*zpq + zqc + zqd) + zpa)*Ireal(2, q + u) + piD*(eta_p**2*(-eta_q*zpq**3 + zpa*zpq**2) + eta_p*(eta_q*zpq**2*(-zqc - zqd) + zpa*zpq*(zqc + zqd)) + zqc*zqd*(-eta_q*zpq + zpa)))) 
      sum    = (inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = inv_ax*(cxpb*der_I_A(1, 0, n, I_A, A) + sxpb*der_I_A(0, 1, n, I_A, A))
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3200) ! | pz    py   s    s     ( 225 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypb)*(-eta_q*zpq + zpa)) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3201) ! | pz    py   s    px    ( 226 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypb)*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3202) ! | pz    py   s    py    ( 227 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*zpq + zpa)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqd*(-eta_q*ypq + ypb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3203) ! | pz    py   s    pz    ( 228 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*ypq + ypb)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqd*(-eta_q*zpq + zpa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3210) ! | pz    py   px   s     ( 229 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypb)*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3211) ! | pz    py   px   px    ( 230 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD2*(-eta_q*ypq + ypb)*(-eta_q*zpq + zpa)) 
      sum    = (inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3212) ! | pz    py   px   py    ( 231 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*zpq + zpa)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqd*(-eta_q*ypq + ypb)))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3213) ! | pz    py   px   pz    ( 232 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*ypq + ypb)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqd*(-eta_q*zpq + zpa)))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3220) ! | pz    py   py   s     ( 233 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*zpq + zpa)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqc*(-eta_q*ypq + ypb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3221) ! | pz    py   py   px    ( 234 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*zpq + zpa)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqc*(-eta_q*ypq + ypb)))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3222) ! | pz    py   py   py    ( 235 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*zpq + zpa)*(Ireal_0_ppt*(-eta_q*ypq + eta_t*(2*eta_p*ypq + yqc + yqd) + ypb)*Ireal(2, q + u) + piD*(eta_p**2*(-eta_q*ypq**3 + ypb*ypq**2) + eta_p*(eta_q*ypq**2*(-yqc - yqd) + ypb*ypq*(yqc + yqd)) + yqc*yqd*(-eta_q*ypq + ypb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3223) ! | pz    py   py   pz    ( 236 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * ((Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqc*(-eta_q*ypq + ypb)))*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqd*(-eta_q*zpq + zpa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3230) ! | pz    py   pz   s     ( 237 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*ypq + ypb)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqc*(-eta_q*zpq + zpa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3231) ! | pz    py   pz   px    ( 238 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*ypq + ypb)*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqc*(-eta_q*zpq + zpa)))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3232) ! | pz    py   pz   py    ( 239 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * ((Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*ypq**2 + ypb*ypq) + yqd*(-eta_q*ypq + ypb)))*(Ireal_0_ppt*eta_t*Ireal(2, q + u) + piD*(eta_p*(-eta_q*zpq**2 + zpa*zpq) + zqc*(-eta_q*zpq + zpa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3233) ! | pz    py   pz   pz    ( 240 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(-eta_q*ypq + ypb)*(Ireal_0_ppt*(-eta_q*zpq + eta_t*(2*eta_p*zpq + zqc + zqd) + zpa)*Ireal(2, q + u) + piD*(eta_p**2*(-eta_q*zpq**3 + zpa*zpq**2) + eta_p*(eta_q*zpq**2*(-zqc - zqd) + zpa*zpq*(zqc + zqd)) + zqc*zqd*(-eta_q*zpq + zpa)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3300) ! | pz    pz   s    s     ( 241 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3301) ! | pz    pz   s    px    ( 242 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3302) ! | pz    pz   s    py    ( 243 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqd)*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3303) ! | pz    pz   s    pz    ( 244 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(eta_t**2*(eta_p*zpq + zqd) + eta_t*(-2*eta_q*zpq + zpa + zpb))*Ireal(2, q + u) + Ireal_0_upq*(eta_p*zpq + zqd)*Ireal(2, p + t2) + piD*(eta_p*(eta_q**2*zpq**3 + eta_q*zpq**2*(-zpa - zpb) + zpa*zpb*zpq) + zqd*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3310) ! | pz    pz   px   s     ( 245 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3311) ! | pz    pz   px   px    ( 246 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb))) 
      sum    = (inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax2*(cxqc*cxqd*der_I_B(2, 0, n, I_B, B) + cxqc*sxqd*der_I_B(1, 1, n, I_B, B) + cxqd*sxqc*der_I_B(1, 1, n, I_B, B) + sxqc*sxqd*der_I_B(0, 2, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3312) ! | pz    pz   px   py    ( 247 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqd)*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3313) ! | pz    pz   px   pz    ( 248 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(eta_t**2*(eta_p*zpq + zqd) + eta_t*(-2*eta_q*zpq + zpa + zpb))*Ireal(2, q + u) + Ireal_0_upq*(eta_p*zpq + zqd)*Ireal(2, p + t2) + piD*(eta_p*(eta_q**2*zpq**3 + eta_q*zpq**2*(-zpa - zpb) + zpa*zpb*zpq) + zqd*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb)))) 
      sum    = (inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqc*der_I_B(1, 0, n, I_B, B) + sxqc*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3320) ! | pz    pz   py   s     ( 249 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqc)*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3321) ! | pz    pz   py   px    ( 250 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqc)*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3322) ! | pz    pz   py   py    ( 251 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * ((Ireal_0_ppt*Ireal(2, q + u) + piD*(eta_p**2*ypq**2 + eta_p*ypq*(yqc + yqd) + yqc*yqd))*(Ireal_0_ppt*eta_t**2*Ireal(2, q + u) + Ireal_0_upq*Ireal(2, p + t2) + piD*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3323) ! | pz    pz   py   pz    ( 252 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqc)*(Ireal_0_ppt*(eta_t**2*(eta_p*zpq + zqd) + eta_t*(-2*eta_q*zpq + zpa + zpb))*Ireal(2, q + u) + Ireal_0_upq*(eta_p*zpq + zqd)*Ireal(2, p + t2) + piD*(eta_p*(eta_q**2*zpq**3 + eta_q*zpq**2*(-zpa - zpb) + zpa*zpb*zpq) + zqd*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3330) ! | pz    pz   pz   s     ( 253 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(eta_t**2*(eta_p*zpq + zqc) + eta_t*(-2*eta_q*zpq + zpa + zpb))*Ireal(2, q + u) + Ireal_0_upq*(eta_p*zpq + zqc)*Ireal(2, p + t2) + piD*(eta_p*(eta_q**2*zpq**3 + eta_q*zpq**2*(-zpa - zpb) + zpa*zpb*zpq) + zqc*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3331) ! | pz    pz   pz   px    ( 254 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(eta_t**2*(eta_p*zpq + zqc) + eta_t*(-2*eta_q*zpq + zpa + zpb))*Ireal(2, q + u) + Ireal_0_upq*(eta_p*zpq + zqc)*Ireal(2, p + t2) + piD*(eta_p*(eta_q**2*zpq**3 + eta_q*zpq**2*(-zpa - zpb) + zpa*zpb*zpq) + zqc*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb)))) 
      sum    = (inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))*I_A(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = inv_ax*(cxqd*der_I_B(1, 0, n, I_B, B) + sxqd*der_I_B(0, 1, n, I_B, B))
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3332) ! | pz    pz   pz   py    ( 255 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(eta_p*ypq + yqd)*(Ireal_0_ppt*(eta_t**2*(eta_p*zpq + zqc) + eta_t*(-2*eta_q*zpq + zpa + zpb))*Ireal(2, q + u) + Ireal_0_upq*(eta_p*zpq + zqc)*Ireal(2, p + t2) + piD*(eta_p*(eta_q**2*zpq**3 + eta_q*zpq**2*(-zpa - zpb) + zpa*zpb*zpq) + zqc*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb)))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case (3333) ! | pz    pz   pz   pz    ( 256 )
      n = 0
      const  = dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) * (piD*(Ireal_0_ppt*(eta_t**2*((eta_p**2*zpq**2 + eta_p*zpq*(zqc + zqd) + zqc*zqd)*Ireal(2, q + u) + Ireal(4, q + u)) + (eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + eta_t*(eta_p*(-4*eta_q*zpq**2 + zpq*(2*zpa + 2*zpb)) + eta_q*zpq*(-2*zqc - 2*zqd) + zpa*(zqc + zqd) + zpb*(zqc + zqd)) + zpa*zpb)*Ireal(2, q + u)) + piD*(eta_p**2*(eta_q**2*zpq**4 + eta_q*zpq**3*(-zpa - zpb) + zpa*zpb*zpq**2) + eta_p*(eta_q**2*zpq**3*(zqc + zqd) + eta_q*zpq**2*(zpa*(-zqc - zqd) + zpb*(-zqc - zqd)) + zpa*zpb*zpq*(zqc + zqd)) + zqc*zqd*(eta_q**2*zpq**2 + eta_q*zpq*(-zpa - zpb) + zpa*zpb)) + (Ireal_0_upq*(eta_p**2*zpq**2 + eta_p*zpq*(zqc + zqd) + zqc*zqd) + Ireal(2, q + u))*Ireal(2, p + t2))) 
      sum    = (I_A(n)*I_B(n)) * I_C_table_x(n, i_quad)
      do n   = 1 , Nmax
      termAn = I_A(n)
      termBn = I_B(n)
      term   = termAn * termBn * I_C_table_x(n, i_quad) * current_term
      sum = sum + 2.d0 * real(term)
      current_term = current_term * expo_term
      if (abs(term) < eps * dabs(sum) ) then
        no_converged_x = .false. 
        exit
      end if
      end do

      case default
        sum = 0.0d0

      end select

      sum = sum * const 

 end function S


end subroutine integrate_ERI_sum



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