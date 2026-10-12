subroutine integrate_ERI_sum_new(pattern_id,p,q,p_x,q_x,phi_x,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq &
     &                           ,ya,yb,yc,yd,yp,yq &
     &                           ,za,zb,zc,zd,zp,zq,result)
      
      use iso_c_binding
      use torus_init
      use bessel_functions
      use files
      use table_lookup_module
      use table_1d_lookup 
      use gauss_legendre_quadrature 
      use precomputed_bessel
      use functions

      use, intrinsic                     :: ieee_arithmetic

      implicit none

      ! Input parameters
      double precision, intent(in)       :: p_x , p
      double precision, intent(in)       :: q_x , q
      double precision, intent(in)       :: phi_x
      double precision, intent(in)       :: xpA , xpB 
      double precision, intent(in)       :: xqC , xqD
      double precision, intent(in)       :: xa, xb, xc, xd, xp, xq 
      double precision, intent(in)       :: ya, yb, yc, yd, yp, yq
      double precision, intent(in)       :: za, zb, zc, zd, zp, zq
      integer         , intent(in)       :: pattern_id

      ! Output parameters

      double precision, intent(out)      :: result
    
      ! Local variables
      
      integer                            :: Nmax_A , Nmax_B , Nmax_C
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

      double precision                   :: AAx , BBx , CCx

      ! --------------------------------------------------------------- !
      ! backward ! 

      double precision                   :: I_A_x(0:50000)
      double precision                   :: I_B_x(0:50000)
      double precision                   ::   F_x(0:50000)
      
      integer                            :: n

      ! --------------------------------------------------------------- !

      integer                              :: Nmax_x(64)
      double precision                     :: const(64)
      double precision                     :: t2(64), D(64), D2(64), piD(64), piD2(64)
      double precision                     :: eta_t(64), eta_p(64), eta_q(64), u(64)
      double precision                     :: Ireal_0_ppt(64), Ireal_0_upq(64)

      double precision                     :: real_An, real_Bn
      double precision                     :: imag_An, imag_Bn

      ! --------------------------------------------------------------- !

      ! =/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/= !
      !            calculate the integrals for the ERI sum              !
      ! =/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/= !


      ! calculate some constants for the x contribution !

      inv_ax  = 1.d0/ax 
      inv_ax2 = inv_ax * inv_ax 

      AAx   = 2.d0  * p_x   * inv_ax2
      BBx   = 2.d0  * q_x   * inv_ax2

      ! --------------------------------------------------------------- !

      Nmax_A = get_Nmax_from_bessel_1D(AAx)
      call bessel_I_scaled_backward(Nmax_A+20, AAx, I_A_x)

      Nmax_B = get_Nmax_from_bessel_1D(BBx)
      call bessel_I_scaled_backward(Nmax_B+20, BBx, I_B_x)


      ! /////////////////////////////////////////////////////////////// !

      ! Clifford Gaussians ! 

      cxpa  = dcos(xpA) ; cxpb  = dcos(xpB)
      sxpa  = dsin(xpA) ; sxpb  = dsin(xpB)

      cxqc  = dcos(xqC) ; cxqd  = dcos(xqD)
      sxqc  = dsin(xqC) ; sxqd  = dsin(xqD)

      c2xpab = dcos(ax*(2.d0*xp-xA-xB))
      s2xpab = dsin(ax*(2.d0*xp-xA-xB))

      c2xqcd = dcos(ax*(2.d0*xq-xC-xD))
      s2xqcd = dsin(ax*(2.d0*xq-xC-xD))

      ! /////////////////////////////////////////////////////////////// !

      ! Real Gaussians ! 

      ypq    = yp - yq ; zpq    = zp - zq

      ypa    = yp - ya ; zpa    = zp - za  
      ypb    = yp - yb ; zpb    = zp - zb  
      yqc    = yq - yc ; zqc    = zq - zc 
      yqd    = yq - yd ; zqd    = zq - zd 

      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !

      do n = 1, 64
       CCx = inv_ax2 * 2.d0 * t_val(n) * t_val(n)
       Nmax_x(n) = get_Nmax(AAx, BBx, CCx, phi_x) + 10
       !Nmax_C    = get_Nmax_from_bessel_1D(CCx)
       !Nmax_x(n) = min(Nmax_A, Nmax_B, Nmax_C) + 10
      end do

      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !

      !                     const for t constants                       ! 

      do n = 1, 64 
        t2(n)          = t_val(n) * t_val(n)
        D(n)           = 1.d0/(dsqrt(p*q+(p+q)*t2(n)))
        D2(n)          = D(n) * D(n)
        piD(n)         = pi * D(n)
        piD2(n)        = piD(n)   * piD(n)
        eta_t(n)       = t2(n) / (p+t2(n))    
        eta_p(n)       = p*t2(n)  * D2(n)
        eta_q(n)       = q*t2(n)  * D2(n)
        u(n)           = eta_t(n) * p
        Ireal_0_ppt(n) = Ireal(0,p+t2(n))
        Ireal_0_upq(n) = Ireal(0,u(n)+q)
      end do 

      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !

      !         calculate the constant for y,z contributions            !

      select case (pattern_id)

      case (0000) ! | s     s    s    s     ( 1 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (0001) ! | s     s    s    px    ( 2 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (0002) ! | s     s    s    py    ( 3 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)) 
      end do

      case (0003) ! | s     s    s    pz    ( 4 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)) 
      end do

      case (0010) ! | s     s    px   s     ( 5 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (0011) ! | s     s    px   px    ( 6 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (0012) ! | s     s    px   py    ( 7 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)) 
      end do

      case (0013) ! | s     s    px   pz    ( 8 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)) 
      end do

      case (0020) ! | s     s    py   s     ( 9 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)) 
      end do

      case (0021) ! | s     s    py   px    ( 10 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)) 
      end do

      case (0022) ! | s     s    py   py    ( 11 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))**2*ypq**2 + pi*D(n)*(eta_p(n))*ypq*yqc + pi*D(n)*(eta_p(n))*ypq*yqd + pi*D(n)*yqc*yqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (0023) ! | s     s    py   pz    ( 12 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)*((eta_p(n))*zpq + zqd)) 
      end do

      case (0030) ! | s     s    pz   s     ( 13 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)) 
      end do

      case (0031) ! | s     s    pz   px    ( 14 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)) 
      end do

      case (0032) ! | s     s    pz   py    ( 15 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)*((eta_p(n))*zpq + zqc)) 
      end do

      case (0033) ! | s     s    pz   pz    ( 16 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))**2*zpq**2 + pi*D(n)*(eta_p(n))*zpq*zqc + pi*D(n)*(eta_p(n))*zpq*zqd + pi*D(n)*zqc*zqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (0100) ! | s     px   s    s     ( 17 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (0101) ! | s     px   s    px    ( 18 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (0102) ! | s     px   s    py    ( 19 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)) 
      end do
  
      case (0103) ! | s     px   s    pz    ( 20 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)) 
      end do

      case (0110) ! | s     px   px   s     ( 21 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (0111) ! | s     px   px   px    ( 22 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (0112) ! | s     px   px   py    ( 23 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)) 
      end do

      case (0113) ! | s     px   px   pz    ( 24 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)) 
      end do

      case (0120) ! | s     px   py   s     ( 25 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)) 
      end do

      case (0121) ! | s     px   py   px    ( 26 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)) 
      end do

      case (0122) ! | s     px   py   py    ( 27 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))**2*ypq**2 + pi*D(n)*(eta_p(n))*ypq*yqc + pi*D(n)*(eta_p(n))*ypq*yqd + pi*D(n)*yqc*yqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (0123) ! | s     px   py   pz    ( 28 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)*((eta_p(n))*zpq + zqd)) 
      end do

      case (0130) ! | s     px   pz   s     ( 29 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)) 
      end do

      case (0131) ! | s     px   pz   px    ( 30 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)) 
      end do

      case (0132) ! | s     px   pz   py    ( 31 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)*((eta_p(n))*zpq + zqc)) 
      end do

      case (0133) ! | s     px   pz   pz    ( 32 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))**2*zpq**2 + pi*D(n)*(eta_p(n))*zpq*zqc + pi*D(n)*(eta_p(n))*zpq*zqd + pi*D(n)*zqc*zqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (0200) ! | s     py   s    s     ( 33 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (0201) ! | s     py   s    px    ( 34 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (0202) ! | s     py   s    py    ( 35 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypb*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (0203) ! | s     py   s    pz    ( 36 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (0210) ! | s     py   px   s     ( 37 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (0211) ! | s     py   px   px    ( 38 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (0212) ! | s     py   px   py    ( 39 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypb*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (0213) ! | s     py   px   pz    ( 40 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (0220) ! | s     py   py   s     ( 41 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypb*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (0221) ! | s     py   py   px    ( 42 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypb*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (0222) ! | s     py   py   py    ( 43 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))**2*(eta_q(n))*ypq**3 + pi*D(n)*(eta_p(n))**2*ypb*ypq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2*yqc - pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2*yqd + pi*D(n)*(eta_p(n))*ypb*ypq*yqc + pi*D(n)*(eta_p(n))*ypb*ypq*yqd - pi*D(n)*(eta_q(n))*ypq*yqc*yqd + pi*D(n)*ypb*yqc*yqd + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*ypq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*yqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*yqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*ypb*Ireal(2, q + u(n)))) 
      end do

      case (0223) ! | s     py   py   pz    ( 44 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqd)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypb*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (0230) ! | s     py   pz   s     ( 45 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (0231) ! | s     py   pz   px    ( 46 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (0232) ! | s     py   pz   py    ( 47 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqc)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypb*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (0233) ! | s     py   pz   pz    ( 48 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypb)*(pi*D(n)*(eta_p(n))**2*zpq**2 + pi*D(n)*(eta_p(n))*zpq*zqc + pi*D(n)*(eta_p(n))*zpq*zqd + pi*D(n)*zqc*zqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (0300) ! | s     pz   s    s     ( 49 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (0301) ! | s     pz   s    px    ( 50 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (0302) ! | s     pz   s    py    ( 51 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (0303) ! | s     pz   s    pz    ( 52 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpb*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (0310) ! | s     pz   px   s     ( 53 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (0311) ! | s     pz   px   px    ( 54 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (0312) ! | s     pz   px   py    ( 55 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (0313) ! | s     pz   px   pz    ( 56 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpb*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (0320) ! | s     pz   py   s     ( 57 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (0321) ! | s     pz   py   px    ( 58 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (0322) ! | s     pz   py   py    ( 59 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpb)*(pi*D(n)*(eta_p(n))**2*ypq**2 + pi*D(n)*(eta_p(n))*ypq*yqc + pi*D(n)*(eta_p(n))*ypq*yqd + pi*D(n)*yqc*yqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (0323) ! | s     pz   py   pz    ( 60 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqc)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpb*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (0330) ! | s     pz   pz   s     ( 61 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpb*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (0331) ! | s     pz   pz   px    ( 62 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpb*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (0332) ! | s     pz   pz   py    ( 63 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqd)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpb*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (0333) ! | s     pz   pz   pz    ( 64 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))**2*(eta_q(n))*zpq**3 + pi*D(n)*(eta_p(n))**2*zpb*zpq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2*zqc - pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2*zqd + pi*D(n)*(eta_p(n))*zpb*zpq*zqc + pi*D(n)*(eta_p(n))*zpb*zpq*zqd - pi*D(n)*(eta_q(n))*zpq*zqc*zqd + pi*D(n)*zpb*zqc*zqd + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*zpq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*zpb*Ireal(2, q + u(n)))) 
      end do

      case (1000) ! | px    s    s    s     ( 65 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (1001) ! | px    s    s    px    ( 66 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (1002) ! | px    s    s    py    ( 67 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)) 
      end do

      case (1003) ! | px    s    s    pz    ( 68 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)) 
      end do

      case (1010) ! | px    s    px   s     ( 69 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (1011) ! | px    s    px   px    ( 70 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (1012) ! | px    s    px   py    ( 71 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)) 
      end do

      case (1013) ! | px    s    px   pz    ( 72 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)) 
      end do

      case (1020) ! | px    s    py   s     ( 73 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)) 
      end do

      case (1021) ! | px    s    py   px    ( 74 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)) 
      end do

      case (1022) ! | px    s    py   py    ( 75 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))**2*ypq**2 + pi*D(n)*(eta_p(n))*ypq*yqc + pi*D(n)*(eta_p(n))*ypq*yqd + pi*D(n)*yqc*yqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (1023) ! | px    s    py   pz    ( 76 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)*((eta_p(n))*zpq + zqd)) 
      end do

      case (1030) ! | px    s    pz   s     ( 77 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)) 
      end do

      case (1031) ! | px    s    pz   px    ( 78 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)) 
      end do

      case (1032) ! | px    s    pz   py    ( 79 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)*((eta_p(n))*zpq + zqc)) 
      end do

      case (1033) ! | px    s    pz   pz    ( 80 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))**2*zpq**2 + pi*D(n)*(eta_p(n))*zpq*zqc + pi*D(n)*(eta_p(n))*zpq*zqd + pi*D(n)*zqc*zqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (1100) ! | px    px   s    s     ( 81 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (1101) ! | px    px   s    px    ( 82 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (1102) ! | px    px   s    py    ( 83 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)) 
      end do

      case (1103) ! | px    px   s    pz    ( 84 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)) 
      end do

      case (1110) ! | px    px   px   s     ( 85 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (1111) ! | px    px   px   px    ( 86 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2) 
      end do

      case (1112) ! | px    px   px   py    ( 87 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)) 
      end do

      case (1113) ! | px    px   px   pz    ( 88 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)) 
      end do

      case (1120) ! | px    px   py   s     ( 89 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)) 
      end do

      case (1121) ! | px    px   py   px    ( 90 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)) 
      end do

      case (1122) ! | px    px   py   py    ( 91 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))**2*ypq**2 + pi*D(n)*(eta_p(n))*ypq*yqc + pi*D(n)*(eta_p(n))*ypq*yqd + pi*D(n)*yqc*yqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (1123) ! | px    px   py   pz    ( 92 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)*((eta_p(n))*zpq + zqd)) 
      end do

      case (1130) ! | px    px   pz   s     ( 93 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)) 
      end do

      case (1131) ! | px    px   pz   px    ( 94 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)) 
      end do

      case (1132) ! | px    px   pz   py    ( 95 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)*((eta_p(n))*zpq + zqc)) 
      end do

      case (1133) ! | px    px   pz   pz    ( 96 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))**2*zpq**2 + pi*D(n)*(eta_p(n))*zpq*zqc + pi*D(n)*(eta_p(n))*zpq*zqd + pi*D(n)*zqc*zqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (1200) ! | px    py   s    s     ( 97 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (1201) ! | px    py   s    px    ( 98 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (1202) ! | px    py   s    py    ( 99 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypb*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (1203) ! | px    py   s    pz    ( 100 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (1210) ! | px    py   px   s     ( 101 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (1211) ! | px    py   px   px    ( 102 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (1212) ! | px    py   px   py    ( 103 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypb*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (1213) ! | px    py   px   pz    ( 104 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (1220) ! | px    py   py   s     ( 105 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypb*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (1221) ! | px    py   py   px    ( 106 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypb*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (1222) ! | px    py   py   py    ( 107 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))**2*(eta_q(n))*ypq**3 + pi*D(n)*(eta_p(n))**2*ypb*ypq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2*yqc - pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2*yqd + pi*D(n)*(eta_p(n))*ypb*ypq*yqc + pi*D(n)*(eta_p(n))*ypb*ypq*yqd - pi*D(n)*(eta_q(n))*ypq*yqc*yqd + pi*D(n)*ypb*yqc*yqd + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*ypq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*yqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*yqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*ypb*Ireal(2, q + u(n)))) 
      end do

      case (1223) ! | px    py   py   pz    ( 108 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqd)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypb*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (1230) ! | px    py   pz   s     ( 109 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (1231) ! | px    py   pz   px    ( 110 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)*(-(eta_q(n))*ypq + ypb)) 
      end do

      case (1232) ! | px    py   pz   py    ( 111 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqc)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypb*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (1233) ! | px    py   pz   pz    ( 112 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypb)*(pi*D(n)*(eta_p(n))**2*zpq**2 + pi*D(n)*(eta_p(n))*zpq*zqc + pi*D(n)*(eta_p(n))*zpq*zqd + pi*D(n)*zqc*zqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (1300) ! | px    pz   s    s     ( 113 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (1301) ! | px    pz   s    px    ( 114 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (1302) ! | px    pz   s    py    ( 115 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (1303) ! | px    pz   s    pz    ( 116 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpb*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (1310) ! | px    pz   px   s     ( 117 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (1311) ! | px    pz   px   px    ( 118 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (1312) ! | px    pz   px   py    ( 119 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (1313) ! | px    pz   px   pz    ( 120 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpb*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (1320) ! | px    pz   py   s     ( 121 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (1321) ! | px    pz   py   px    ( 122 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (1322) ! | px    pz   py   py    ( 123 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpb)*(pi*D(n)*(eta_p(n))**2*ypq**2 + pi*D(n)*(eta_p(n))*ypq*yqc + pi*D(n)*(eta_p(n))*ypq*yqd + pi*D(n)*yqc*yqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (1323) ! | px    pz   py   pz    ( 124 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqc)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpb*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (1330) ! | px    pz   pz   s     ( 125 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpb*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (1331) ! | px    pz   pz   px    ( 126 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpb*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (1332) ! | px    pz   pz   py    ( 127 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqd)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpb*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (1333) ! | px    pz   pz   pz    ( 128 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))**2*(eta_q(n))*zpq**3 + pi*D(n)*(eta_p(n))**2*zpb*zpq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2*zqc - pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2*zqd + pi*D(n)*(eta_p(n))*zpb*zpq*zqc + pi*D(n)*(eta_p(n))*zpb*zpq*zqd - pi*D(n)*(eta_q(n))*zpq*zqc*zqd + pi*D(n)*zpb*zqc*zqd + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*zpq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*zpb*Ireal(2, q + u(n)))) 
      end do

      case (2000) ! | py    s    s    s     ( 129 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2001) ! | py    s    s    px    ( 130 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2002) ! | py    s    s    py    ( 131 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypa*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2003) ! | py    s    s    pz    ( 132 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2010) ! | py    s    px   s     ( 133 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2011) ! | py    s    px   px    ( 134 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2012) ! | py    s    px   py    ( 135 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypa*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2013) ! | py    s    px   pz    ( 136 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2020) ! | py    s    py   s     ( 137 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypa*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2021) ! | py    s    py   px    ( 138 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypa*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2022) ! | py    s    py   py    ( 139 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))**2*(eta_q(n))*ypq**3 + pi*D(n)*(eta_p(n))**2*ypa*ypq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2*yqc - pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2*yqd + pi*D(n)*(eta_p(n))*ypa*ypq*yqc + pi*D(n)*(eta_p(n))*ypa*ypq*yqd - pi*D(n)*(eta_q(n))*ypq*yqc*yqd + pi*D(n)*ypa*yqc*yqd + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*ypq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*yqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*yqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*ypa*Ireal(2, q + u(n)))) 
      end do

      case (2023) ! | py    s    py   pz    ( 140 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqd)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypa*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2030) ! | py    s    pz   s     ( 141 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2031) ! | py    s    pz   px    ( 142 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2032) ! | py    s    pz   py    ( 143 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqc)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypa*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2033) ! | py    s    pz   pz    ( 144 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypa)*(pi*D(n)*(eta_p(n))**2*zpq**2 + pi*D(n)*(eta_p(n))*zpq*zqc + pi*D(n)*(eta_p(n))*zpq*zqd + pi*D(n)*zqc*zqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (2100) ! | py    px   s    s     ( 145 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2101) ! | py    px   s    px    ( 146 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2102) ! | py    px   s    py    ( 147 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypa*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2103) ! | py    px   s    pz    ( 148 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2110) ! | py    px   px   s     ( 149 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2111) ! | py    px   px   px    ( 150 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2112) ! | py    px   px   py    ( 151 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypa*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2113) ! | py    px   px   pz    ( 152 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqd)*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2120) ! | py    px   py   s     ( 153 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypa*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2121) ! | py    px   py   px    ( 154 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypa*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2122) ! | py    px   py   py    ( 155 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))**2*(eta_q(n))*ypq**3 + pi*D(n)*(eta_p(n))**2*ypa*ypq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2*yqc - pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2*yqd + pi*D(n)*(eta_p(n))*ypa*ypq*yqc + pi*D(n)*(eta_p(n))*ypa*ypq*yqd - pi*D(n)*(eta_q(n))*ypq*yqc*yqd + pi*D(n)*ypa*yqc*yqd + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*ypq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*yqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*yqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*ypa*Ireal(2, q + u(n)))) 
      end do

      case (2123) ! | py    px   py   pz    ( 156 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqd)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypa*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2130) ! | py    px   pz   s     ( 157 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2131) ! | py    px   pz   px    ( 158 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*zpq + zqc)*(-(eta_q(n))*ypq + ypa)) 
      end do

      case (2132) ! | py    px   pz   py    ( 159 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqc)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypa*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2133) ! | py    px   pz   pz    ( 160 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypa)*(pi*D(n)*(eta_p(n))**2*zpq**2 + pi*D(n)*(eta_p(n))*zpq*zqc + pi*D(n)*(eta_p(n))*zpq*zqd + pi*D(n)*zqc*zqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (2200) ! | py    py   s    s     ( 161 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_q(n))**2*ypq**2 - pi*D(n)*(eta_q(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypb*ypq + pi*D(n)*ypa*ypb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (2201) ! | py    py   s    px    ( 162 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_q(n))**2*ypq**2 - pi*D(n)*(eta_q(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypb*ypq + pi*D(n)*ypa*ypb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (2202) ! | py    py   s    py    ( 163 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))*(eta_q(n))**2*ypq**3 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypa*ypq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypb*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypb*ypq + pi*D(n)*(eta_q(n))**2*ypq**2*yqd - pi*D(n)*(eta_q(n))*ypa*ypq*yqd - pi*D(n)*(eta_q(n))*ypb*ypq*yqd + pi*D(n)*ypa*ypb*yqd + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*ypq*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*yqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypa*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))*ypq*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*yqd*Ireal(2, p + t2(n)))) 
      end do

      case (2203) ! | py    py   s    pz    ( 164 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqd)*(pi*D(n)*(eta_q(n))**2*ypq**2 - pi*D(n)*(eta_q(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypb*ypq + pi*D(n)*ypa*ypb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (2210) ! | py    py   px   s     ( 165 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_q(n))**2*ypq**2 - pi*D(n)*(eta_q(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypb*ypq + pi*D(n)*ypa*ypb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (2211) ! | py    py   px   px    ( 166 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_q(n))**2*ypq**2 - pi*D(n)*(eta_q(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypb*ypq + pi*D(n)*ypa*ypb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (2212) ! | py    py   px   py    ( 167 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))*(eta_q(n))**2*ypq**3 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypa*ypq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypb*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypb*ypq + pi*D(n)*(eta_q(n))**2*ypq**2*yqd - pi*D(n)*(eta_q(n))*ypa*ypq*yqd - pi*D(n)*(eta_q(n))*ypb*ypq*yqd + pi*D(n)*ypa*ypb*yqd + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*ypq*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*yqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypa*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))*ypq*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*yqd*Ireal(2, p + t2(n)))) 
      end do

      case (2213) ! | py    py   px   pz    ( 168 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqd)*(pi*D(n)*(eta_q(n))**2*ypq**2 - pi*D(n)*(eta_q(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypb*ypq + pi*D(n)*ypa*ypb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (2220) ! | py    py   py   s     ( 169 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))*(eta_q(n))**2*ypq**3 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypa*ypq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypb*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypb*ypq + pi*D(n)*(eta_q(n))**2*ypq**2*yqc - pi*D(n)*(eta_q(n))*ypa*ypq*yqc - pi*D(n)*(eta_q(n))*ypb*ypq*yqc + pi*D(n)*ypa*ypb*yqc + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*ypq*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*yqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypa*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))*ypq*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*yqc*Ireal(2, p + t2(n)))) 
      end do

      case (2221) ! | py    py   py   px    ( 170 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))*(eta_q(n))**2*ypq**3 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypa*ypq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypb*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypb*ypq + pi*D(n)*(eta_q(n))**2*ypq**2*yqc - pi*D(n)*(eta_q(n))*ypa*ypq*yqc - pi*D(n)*(eta_q(n))*ypb*ypq*yqc + pi*D(n)*ypa*ypb*yqc + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*ypq*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*yqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypa*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))*ypq*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*yqc*Ireal(2, p + t2(n)))) 
      end do

      case (2222) ! | py    py   py   py    ( 171 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))**2*(eta_q(n))**2*ypq**4 - pi*D(n)*(eta_p(n))**2*(eta_q(n))*ypa*ypq**3 - pi*D(n)*(eta_p(n))**2*(eta_q(n))*ypb*ypq**3 + pi*D(n)*(eta_p(n))**2*ypa*ypb*ypq**2 + pi*D(n)*(eta_p(n))*(eta_q(n))**2*ypq**3*yqc + pi*D(n)*(eta_p(n))*(eta_q(n))**2*ypq**3*yqd - pi*D(n)*(eta_p(n))*(eta_q(n))*ypa*ypq**2*yqc - pi*D(n)*(eta_p(n))*(eta_q(n))*ypa*ypq**2*yqd - pi*D(n)*(eta_p(n))*(eta_q(n))*ypb*ypq**2*yqc - pi*D(n)*(eta_p(n))*(eta_q(n))*ypb*ypq**2*yqd + pi*D(n)*(eta_p(n))*ypa*ypb*ypq*yqc + pi*D(n)*(eta_p(n))*ypa*ypb*ypq*yqd + pi*D(n)*(eta_q(n))**2*ypq**2*yqc*yqd - pi*D(n)*(eta_q(n))*ypa*ypq*yqc*yqd - pi*D(n)*(eta_q(n))*ypb*ypq*yqc*yqd + pi*D(n)*ypa*ypb*yqc*yqd + (Ireal_0_ppt(n))*(eta_p(n))**2*(eta_t(n))**2*ypq**2*Ireal(2, q + u(n)) - 4*(Ireal_0_ppt(n))*(eta_p(n))*(eta_q(n))*(eta_t(n))*ypq**2*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*ypq*yqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*ypq*yqd*Ireal(2, q + u(n)) + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*ypa*ypq*Ireal(2, q + u(n)) + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*ypb*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_q(n))**2*ypq**2*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*ypq*yqc*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*ypq*yqd*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*ypa*ypq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*ypb*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*yqc*yqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(4, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypa*yqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypa*yqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypb*yqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypb*yqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*ypa*ypb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))**2*ypq**2*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*(eta_p(n))*ypq*yqc*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*(eta_p(n))*ypq*yqd*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*yqc*yqd*Ireal(2, p + t2(n)) + Ireal(2, p + t2(n))*Ireal(2, q + u(n)))) 
      end do

      case (2223) ! | py    py   py   pz    ( 172 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqd)*(pi*D(n)*(eta_p(n))*(eta_q(n))**2*ypq**3 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypa*ypq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypb*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypb*ypq + pi*D(n)*(eta_q(n))**2*ypq**2*yqc - pi*D(n)*(eta_q(n))*ypa*ypq*yqc - pi*D(n)*(eta_q(n))*ypb*ypq*yqc + pi*D(n)*ypa*ypb*yqc + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*ypq*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*yqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypa*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))*ypq*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*yqc*Ireal(2, p + t2(n)))) 
      end do

      case (2230) ! | py    py   pz   s     ( 173 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqc)*(pi*D(n)*(eta_q(n))**2*ypq**2 - pi*D(n)*(eta_q(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypb*ypq + pi*D(n)*ypa*ypb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (2231) ! | py    py   pz   px    ( 174 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqc)*(pi*D(n)*(eta_q(n))**2*ypq**2 - pi*D(n)*(eta_q(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypb*ypq + pi*D(n)*ypa*ypb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (2232) ! | py    py   pz   py    ( 175 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*zpq + zqc)*(pi*D(n)*(eta_p(n))*(eta_q(n))**2*ypq**3 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypa*ypq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypb*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypb*ypq + pi*D(n)*(eta_q(n))**2*ypq**2*yqd - pi*D(n)*(eta_q(n))*ypa*ypq*yqd - pi*D(n)*(eta_q(n))*ypb*ypq*yqd + pi*D(n)*ypa*ypb*yqd + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*ypq*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*yqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypa*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*ypb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))*ypq*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*yqd*Ireal(2, p + t2(n)))) 
      end do

      case (2233) ! | py    py   pz   pz    ( 176 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * ((pi*D(n)*(eta_p(n))**2*zpq**2 + pi*D(n)*(eta_p(n))*zpq*zqc + pi*D(n)*(eta_p(n))*zpq*zqd + pi*D(n)*zqc*zqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))*(pi*D(n)*(eta_q(n))**2*ypq**2 - pi*D(n)*(eta_q(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypb*ypq + pi*D(n)*ypa*ypb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (2300) ! | py    pz   s    s     ( 177 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypa)*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (2301) ! | py    pz   s    px    ( 178 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypa)*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (2302) ! | py    pz   s    py    ( 179 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpb)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypa*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2303) ! | py    pz   s    pz    ( 180 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypa)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpb*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2310) ! | py    pz   px   s     ( 181 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypa)*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (2311) ! | py    pz   px   px    ( 182 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypa)*(-(eta_q(n))*zpq + zpb)) 
      end do

      case (2312) ! | py    pz   px   py    ( 183 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpb)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypa*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2313) ! | py    pz   px   pz    ( 184 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypa)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpb*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2320) ! | py    pz   py   s     ( 185 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpb)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypa*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2321) ! | py    pz   py   px    ( 186 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpb)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypa*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2322) ! | py    pz   py   py    ( 187 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpb)*(-pi*D(n)*(eta_p(n))**2*(eta_q(n))*ypq**3 + pi*D(n)*(eta_p(n))**2*ypa*ypq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2*yqc - pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2*yqd + pi*D(n)*(eta_p(n))*ypa*ypq*yqc + pi*D(n)*(eta_p(n))*ypa*ypq*yqd - pi*D(n)*(eta_q(n))*ypq*yqc*yqd + pi*D(n)*ypa*yqc*yqd + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*ypq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*yqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*yqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*ypa*Ireal(2, q + u(n)))) 
      end do

      case (2323) ! | py    pz   py   pz    ( 188 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * ((-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypa*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpb*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2330) ! | py    pz   pz   s     ( 189 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypa)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpb*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2331) ! | py    pz   pz   px    ( 190 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypa)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpb*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2332) ! | py    pz   pz   py    ( 191 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * ((-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypa*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypa*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpb*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpb*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (2333) ! | py    pz   pz   pz    ( 192 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypa)*(-pi*D(n)*(eta_p(n))**2*(eta_q(n))*zpq**3 + pi*D(n)*(eta_p(n))**2*zpb*zpq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2*zqc - pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2*zqd + pi*D(n)*(eta_p(n))*zpb*zpq*zqc + pi*D(n)*(eta_p(n))*zpb*zpq*zqd - pi*D(n)*(eta_q(n))*zpq*zqc*zqd + pi*D(n)*zpb*zqc*zqd + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*zpq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*zpb*Ireal(2, q + u(n)))) 
      end do

      case (3000) ! | pz    s    s    s     ( 193 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3001) ! | pz    s    s    px    ( 194 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3002) ! | pz    s    s    py    ( 195 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3003) ! | pz    s    s    pz    ( 196 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpa*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3010) ! | pz    s    px   s     ( 197 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3011) ! | pz    s    px   px    ( 198 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3012) ! | pz    s    px   py    ( 199 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3013) ! | pz    s    px   pz    ( 200 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpa*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3020) ! | pz    s    py   s     ( 201 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3021) ! | pz    s    py   px    ( 202 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3022) ! | pz    s    py   py    ( 203 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpa)*(pi*D(n)*(eta_p(n))**2*ypq**2 + pi*D(n)*(eta_p(n))*ypq*yqc + pi*D(n)*(eta_p(n))*ypq*yqd + pi*D(n)*yqc*yqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (3023) ! | pz    s    py   pz    ( 204 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqc)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpa*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3030) ! | pz    s    pz   s     ( 205 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpa*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3031) ! | pz    s    pz   px    ( 206 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpa*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3032) ! | pz    s    pz   py    ( 207 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqd)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpa*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3033) ! | pz    s    pz   pz    ( 208 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))**2*(eta_q(n))*zpq**3 + pi*D(n)*(eta_p(n))**2*zpa*zpq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2*zqc - pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2*zqd + pi*D(n)*(eta_p(n))*zpa*zpq*zqc + pi*D(n)*(eta_p(n))*zpa*zpq*zqd - pi*D(n)*(eta_q(n))*zpq*zqc*zqd + pi*D(n)*zpa*zqc*zqd + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*zpq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*zpa*Ireal(2, q + u(n)))) 
      end do

      case (3100) ! | pz    px   s    s     ( 209 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3101) ! | pz    px   s    px    ( 210 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3102) ! | pz    px   s    py    ( 211 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3103) ! | pz    px   s    pz    ( 212 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpa*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3110) ! | pz    px   px   s     ( 213 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3111) ! | pz    px   px   px    ( 214 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3112) ! | pz    px   px   py    ( 215 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqd)*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3113) ! | pz    px   px   pz    ( 216 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpa*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3120) ! | pz    px   py   s     ( 217 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3121) ! | pz    px   py   px    ( 218 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*((eta_p(n))*ypq + yqc)*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3122) ! | pz    px   py   py    ( 219 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpa)*(pi*D(n)*(eta_p(n))**2*ypq**2 + pi*D(n)*(eta_p(n))*ypq*yqc + pi*D(n)*(eta_p(n))*ypq*yqd + pi*D(n)*yqc*yqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))) 
      end do

      case (3123) ! | pz    px   py   pz    ( 220 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqc)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpa*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3130) ! | pz    px   pz   s     ( 221 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpa*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3131) ! | pz    px   pz   px    ( 222 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpa*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3132) ! | pz    px   pz   py    ( 223 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqd)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpa*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3133) ! | pz    px   pz   pz    ( 224 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-pi*D(n)*(eta_p(n))**2*(eta_q(n))*zpq**3 + pi*D(n)*(eta_p(n))**2*zpa*zpq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2*zqc - pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2*zqd + pi*D(n)*(eta_p(n))*zpa*zpq*zqc + pi*D(n)*(eta_p(n))*zpa*zpq*zqd - pi*D(n)*(eta_q(n))*zpq*zqc*zqd + pi*D(n)*zpa*zqc*zqd + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*zpq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*zpa*Ireal(2, q + u(n)))) 
      end do

      case (3200) ! | pz    py   s    s     ( 225 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypb)*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3201) ! | pz    py   s    px    ( 226 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypb)*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3202) ! | pz    py   s    py    ( 227 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpa)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypb*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3203) ! | pz    py   s    pz    ( 228 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypb)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpa*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3210) ! | pz    py   px   s     ( 229 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypb)*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3211) ! | pz    py   px   px    ( 230 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi2*D(n)**2*(-(eta_q(n))*ypq + ypb)*(-(eta_q(n))*zpq + zpa)) 
      end do

      case (3212) ! | pz    py   px   py    ( 231 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpa)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypb*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3213) ! | pz    py   px   pz    ( 232 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypb)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpa*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3220) ! | pz    py   py   s     ( 233 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpa)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypb*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3221) ! | pz    py   py   px    ( 234 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpa)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypb*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3222) ! | pz    py   py   py    ( 235 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*zpq + zpa)*(-pi*D(n)*(eta_p(n))**2*(eta_q(n))*ypq**3 + pi*D(n)*(eta_p(n))**2*ypb*ypq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2*yqc - pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2*yqd + pi*D(n)*(eta_p(n))*ypb*ypq*yqc + pi*D(n)*(eta_p(n))*ypb*ypq*yqd - pi*D(n)*(eta_q(n))*ypq*yqc*yqd + pi*D(n)*ypb*yqc*yqd + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*ypq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*ypq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*yqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*yqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*ypb*Ireal(2, q + u(n)))) 
      end do

      case (3223) ! | pz    py   py   pz    ( 236 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * ((-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqc + pi*D(n)*ypb*yqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqd + pi*D(n)*zpa*zqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3230) ! | pz    py   pz   s     ( 237 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypb)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpa*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3231) ! | pz    py   pz   px    ( 238 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypb)*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpa*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3232) ! | pz    py   pz   py    ( 239 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * ((-pi*D(n)*(eta_p(n))*(eta_q(n))*ypq**2 + pi*D(n)*(eta_p(n))*ypb*ypq - pi*D(n)*(eta_q(n))*ypq*yqd + pi*D(n)*ypb*yqd + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))*(-pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpq*zqc + pi*D(n)*zpa*zqc + (Ireal_0_ppt(n))*(eta_t(n))*Ireal(2, q + u(n)))) 
      end do

      case (3233) ! | pz    py   pz   pz    ( 240 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(-(eta_q(n))*ypq + ypb)*(-pi*D(n)*(eta_p(n))**2*(eta_q(n))*zpq**3 + pi*D(n)*(eta_p(n))**2*zpa*zpq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2*zqc - pi*D(n)*(eta_p(n))*(eta_q(n))*zpq**2*zqd + pi*D(n)*(eta_p(n))*zpa*zpq*zqc + pi*D(n)*(eta_p(n))*zpa*zpq*zqd - pi*D(n)*(eta_q(n))*zpq*zqc*zqd + pi*D(n)*zpa*zqc*zqd + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*zpq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*zpa*Ireal(2, q + u(n)))) 
      end do

      case (3300) ! | pz    pz   s    s     ( 241 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_q(n))**2*zpq**2 - pi*D(n)*(eta_q(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpb*zpq + pi*D(n)*zpa*zpb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (3301) ! | pz    pz   s    px    ( 242 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_q(n))**2*zpq**2 - pi*D(n)*(eta_q(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpb*zpq + pi*D(n)*zpa*zpb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (3302) ! | pz    pz   s    py    ( 243 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqd)*(pi*D(n)*(eta_q(n))**2*zpq**2 - pi*D(n)*(eta_q(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpb*zpq + pi*D(n)*zpa*zpb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (3303) ! | pz    pz   s    pz    ( 244 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))*(eta_q(n))**2*zpq**3 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpa*zpq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpb*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpb*zpq + pi*D(n)*(eta_q(n))**2*zpq**2*zqd - pi*D(n)*(eta_q(n))*zpa*zpq*zqd - pi*D(n)*(eta_q(n))*zpb*zpq*zqd + pi*D(n)*zpa*zpb*zqd + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*zpq*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*zqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpa*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))*zpq*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*zqd*Ireal(2, p + t2(n)))) 
      end do

      case (3310) ! | pz    pz   px   s     ( 245 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_q(n))**2*zpq**2 - pi*D(n)*(eta_q(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpb*zpq + pi*D(n)*zpa*zpb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (3311) ! | pz    pz   px   px    ( 246 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_q(n))**2*zpq**2 - pi*D(n)*(eta_q(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpb*zpq + pi*D(n)*zpa*zpb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (3312) ! | pz    pz   px   py    ( 247 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqd)*(pi*D(n)*(eta_q(n))**2*zpq**2 - pi*D(n)*(eta_q(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpb*zpq + pi*D(n)*zpa*zpb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (3313) ! | pz    pz   px   pz    ( 248 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))*(eta_q(n))**2*zpq**3 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpa*zpq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpb*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpb*zpq + pi*D(n)*(eta_q(n))**2*zpq**2*zqd - pi*D(n)*(eta_q(n))*zpa*zpq*zqd - pi*D(n)*(eta_q(n))*zpb*zpq*zqd + pi*D(n)*zpa*zpb*zqd + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*zpq*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*zqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpa*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))*zpq*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*zqd*Ireal(2, p + t2(n)))) 
      end do

      case (3320) ! | pz    pz   py   s     ( 249 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqc)*(pi*D(n)*(eta_q(n))**2*zpq**2 - pi*D(n)*(eta_q(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpb*zpq + pi*D(n)*zpa*zpb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (3321) ! | pz    pz   py   px    ( 250 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqc)*(pi*D(n)*(eta_q(n))**2*zpq**2 - pi*D(n)*(eta_q(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpb*zpq + pi*D(n)*zpa*zpb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (3322) ! | pz    pz   py   py    ( 251 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * ((pi*D(n)*(eta_p(n))**2*ypq**2 + pi*D(n)*(eta_p(n))*ypq*yqc + pi*D(n)*(eta_p(n))*ypq*yqd + pi*D(n)*yqc*yqd + (Ireal_0_ppt(n))*Ireal(2, q + u(n)))*(pi*D(n)*(eta_q(n))**2*zpq**2 - pi*D(n)*(eta_q(n))*zpa*zpq - pi*D(n)*(eta_q(n))*zpb*zpq + pi*D(n)*zpa*zpb + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*Ireal(2, p + t2(n)))) 
      end do

      case (3323) ! | pz    pz   py   pz    ( 252 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqc)*(pi*D(n)*(eta_p(n))*(eta_q(n))**2*zpq**3 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpa*zpq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpb*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpb*zpq + pi*D(n)*(eta_q(n))**2*zpq**2*zqd - pi*D(n)*(eta_q(n))*zpa*zpq*zqd - pi*D(n)*(eta_q(n))*zpb*zpq*zqd + pi*D(n)*zpa*zpb*zqd + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*zpq*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*zqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpa*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))*zpq*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*zqd*Ireal(2, p + t2(n)))) 
      end do

      case (3330) ! | pz    pz   pz   s     ( 253 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))*(eta_q(n))**2*zpq**3 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpa*zpq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpb*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpb*zpq + pi*D(n)*(eta_q(n))**2*zpq**2*zqc - pi*D(n)*(eta_q(n))*zpa*zpq*zqc - pi*D(n)*(eta_q(n))*zpb*zpq*zqc + pi*D(n)*zpa*zpb*zqc + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*zpq*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*zqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpa*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))*zpq*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*zqc*Ireal(2, p + t2(n)))) 
      end do

      case (3331) ! | pz    pz   pz   px    ( 254 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))*(eta_q(n))**2*zpq**3 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpa*zpq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpb*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpb*zpq + pi*D(n)*(eta_q(n))**2*zpq**2*zqc - pi*D(n)*(eta_q(n))*zpa*zpq*zqc - pi*D(n)*(eta_q(n))*zpb*zpq*zqc + pi*D(n)*zpa*zpb*zqc + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*zpq*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*zqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpa*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))*zpq*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*zqc*Ireal(2, p + t2(n)))) 
      end do

      case (3332) ! | pz    pz   pz   py    ( 255 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*((eta_p(n))*ypq + yqd)*(pi*D(n)*(eta_p(n))*(eta_q(n))**2*zpq**3 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpa*zpq**2 - pi*D(n)*(eta_p(n))*(eta_q(n))*zpb*zpq**2 + pi*D(n)*(eta_p(n))*zpa*zpb*zpq + pi*D(n)*(eta_q(n))**2*zpq**2*zqc - pi*D(n)*(eta_q(n))*zpa*zpq*zqc - pi*D(n)*(eta_q(n))*zpb*zpq*zqc + pi*D(n)*zpa*zpb*zqc + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*zpq*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*zqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpa*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))*zpq*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*zqc*Ireal(2, p + t2(n)))) 
      end do

      case (3333) ! | pz    pz   pz   pz    ( 256 )
      
      do n = 1, 64
        const(n)  = dexp(-p * q * t2(n) * D2(n) * ( ypq*ypq + zpq*zpq ) ) * (pi*D(n)*(pi*D(n)*(eta_p(n))**2*(eta_q(n))**2*zpq**4 - pi*D(n)*(eta_p(n))**2*(eta_q(n))*zpa*zpq**3 - pi*D(n)*(eta_p(n))**2*(eta_q(n))*zpb*zpq**3 + pi*D(n)*(eta_p(n))**2*zpa*zpb*zpq**2 + pi*D(n)*(eta_p(n))*(eta_q(n))**2*zpq**3*zqc + pi*D(n)*(eta_p(n))*(eta_q(n))**2*zpq**3*zqd - pi*D(n)*(eta_p(n))*(eta_q(n))*zpa*zpq**2*zqc - pi*D(n)*(eta_p(n))*(eta_q(n))*zpa*zpq**2*zqd - pi*D(n)*(eta_p(n))*(eta_q(n))*zpb*zpq**2*zqc - pi*D(n)*(eta_p(n))*(eta_q(n))*zpb*zpq**2*zqd + pi*D(n)*(eta_p(n))*zpa*zpb*zpq*zqc + pi*D(n)*(eta_p(n))*zpa*zpb*zpq*zqd + pi*D(n)*(eta_q(n))**2*zpq**2*zqc*zqd - pi*D(n)*(eta_q(n))*zpa*zpq*zqc*zqd - pi*D(n)*(eta_q(n))*zpb*zpq*zqc*zqd + pi*D(n)*zpa*zpb*zqc*zqd + (Ireal_0_ppt(n))*(eta_p(n))**2*(eta_t(n))**2*zpq**2*Ireal(2, q + u(n)) - 4*(Ireal_0_ppt(n))*(eta_p(n))*(eta_q(n))*(eta_t(n))*zpq**2*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*zpq*zqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))**2*zpq*zqd*Ireal(2, q + u(n)) + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*zpa*zpq*Ireal(2, q + u(n)) + 2*(Ireal_0_ppt(n))*(eta_p(n))*(eta_t(n))*zpb*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_q(n))**2*zpq**2*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*zpq*zqc*Ireal(2, q + u(n)) - 2*(Ireal_0_ppt(n))*(eta_q(n))*(eta_t(n))*zpq*zqd*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*zpa*zpq*Ireal(2, q + u(n)) - (Ireal_0_ppt(n))*(eta_q(n))*zpb*zpq*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*zqc*zqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))**2*Ireal(4, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpa*zqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpa*zqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpb*zqc*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*(eta_t(n))*zpb*zqd*Ireal(2, q + u(n)) + (Ireal_0_ppt(n))*zpa*zpb*Ireal(2, q + u(n)) + (Ireal_0_upq(n))*(eta_p(n))**2*zpq**2*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*(eta_p(n))*zpq*zqc*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*(eta_p(n))*zpq*zqd*Ireal(2, p + t2(n)) + (Ireal_0_upq(n))*zqc*zqd*Ireal(2, p + t2(n)) + Ireal(2, p + t2(n))*Ireal(2, q + u(n)))) 
      end do

      case default
        const(n) = 0.d0

      end select


      ! /////////////////////////////////////////////////////////////// !
      ! /////////////////////////////////////////////////////////////// !
      !         calulate the contributions of periodic x axis           !
      ! /////////////////////////////////////////////////////////////// !
      ! /////////////////////////////////////////////////////////////// !

      select case(pattern_id)
      
      case (0000) ! | s     s    s    s     ( 1 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0001) ! | s     s    s    px    ( 2 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0002) ! | s     s    s    py    ( 3 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0003) ! | s     s    s    pz    ( 4 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0010) ! | s     s    px   s     ( 5 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0011) ! | s     s    px   px    ( 6 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn  = 0.d0
        else
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0012) ! | s     s    px   py    ( 7 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0013) ! | s     s    px   pz    ( 8 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0020) ! | s     s    py   s     ( 9 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0021) ! | s     s    py   px    ( 10 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0022) ! | s     s    py   py    ( 11 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0023) ! | s     s    py   pz    ( 12 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0030) ! | s     s    pz   s     ( 13 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0031) ! | s     s    pz   px    ( 14 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0032) ! | s     s    pz   py    ( 15 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0033) ! | s     s    pz   pz    ( 16 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0100) ! | s     px   s    s     ( 17 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0101) ! | s     px   s    px    ( 18 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0102) ! | s     px   s    py    ( 19 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0103) ! | s     px   s    pz    ( 20 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0110) ! | s     px   px   s     ( 21 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0111) ! | s     px   px   px    ( 22 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn  = 0.d0
        else
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0112) ! | s     px   px   py    ( 23 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0113) ! | s     px   px   pz    ( 24 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0120) ! | s     px   py   s     ( 25 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0121) ! | s     px   py   px    ( 26 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0122) ! | s     px   py   py    ( 27 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0123) ! | s     px   py   pz    ( 28 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0130) ! | s     px   pz   s     ( 29 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0131) ! | s     px   pz   px    ( 30 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0132) ! | s     px   pz   py    ( 31 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0133) ! | s     px   pz   pz    ( 32 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (0200) ! | s     py   s    s     ( 33 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0201) ! | s     py   s    px    ( 34 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0202) ! | s     py   s    py    ( 35 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0203) ! | s     py   s    pz    ( 36 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0210) ! | s     py   px   s     ( 37 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0211) ! | s     py   px   px    ( 38 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn  = 0.d0
        else
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0212) ! | s     py   px   py    ( 39 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0213) ! | s     py   px   pz    ( 40 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0220) ! | s     py   py   s     ( 41 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0221) ! | s     py   py   px    ( 42 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0222) ! | s     py   py   py    ( 43 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0223) ! | s     py   py   pz    ( 44 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0230) ! | s     py   pz   s     ( 45 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0231) ! | s     py   pz   px    ( 46 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0232) ! | s     py   pz   py    ( 47 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0233) ! | s     py   pz   pz    ( 48 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0300) ! | s     pz   s    s     ( 49 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0301) ! | s     pz   s    px    ( 50 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0302) ! | s     pz   s    py    ( 51 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0303) ! | s     pz   s    pz    ( 52 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0310) ! | s     pz   px   s     ( 53 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0311) ! | s     pz   px   px    ( 54 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 

        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0312) ! | s     pz   px   py    ( 55 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0313) ! | s     pz   px   pz    ( 56 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0320) ! | s     pz   py   s     ( 57 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0321) ! | s     pz   py   px    ( 58 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0322) ! | s     pz   py   py    ( 59 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0323) ! | s     pz   py   pz    ( 60 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0330) ! | s     pz   pz   s     ( 61 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0331) ! | s     pz   pz   px    ( 62 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (0332) ! | s     pz   pz   py    ( 63 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (0333) ! | s     pz   pz   pz    ( 64 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (1000) ! | px    s    s    s     ( 65 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1001) ! | px    s    s    px    ( 66 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1002) ! | px    s    s    py    ( 67 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1003) ! | px    s    s    pz    ( 68 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1010) ! | px    s    px   s     ( 69 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1011) ! | px    s    px   px    ( 70 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        if (dabs(BBx) < 1.d-300) then
          imag_Bn = 0.d0
        else
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 

        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1012) ! | px    s    px   py    ( 71 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1013) ! | px    s    px   pz    ( 72 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1020) ! | px    s    py   s     ( 73 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1021) ! | px    s    py   px    ( 74 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1022) ! | px    s    py   py    ( 75 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1023) ! | px    s    py   pz    ( 76 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1030) ! | px    s    pz   s     ( 77 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1031) ! | px    s    pz   px    ( 78 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1032) ! | px    s    pz   py    ( 79 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1033) ! | px    s    pz   pz    ( 80 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1100) ! | px    px   s    s     ( 81 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1101) ! | px    px   s    px    ( 82 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1102) ! | px    px   s    py    ( 83 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1103) ! | px    px   s    pz    ( 84 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1110) ! | px    px   px   s     ( 85 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1111) ! | px    px   px   px    ( 86 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        if (dabs(BBx) < 1.d-300) then
          imag_Bn = 0.d0 
        else 
        imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1112) ! | px    px   px   py    ( 87 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1113) ! | px    px   px   pz    ( 88 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1120) ! | px    px   py   s     ( 89 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1121) ! | px    px   py   px    ( 90 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1122) ! | px    px   py   py    ( 91 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1123) ! | px    px   py   pz    ( 92 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1130) ! | px    px   pz   s     ( 93 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1131) ! | px    px   pz   px    ( 94 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1132) ! | px    px   pz   py    ( 95 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1133) ! | px    px   pz   pz    ( 96 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -0.5*inv_ax2*(-1.0*cxpa*cxpb*I_A_x(n)+0.5*cxpa*cxpb*I_A_x(n+2)+0.5*cxpa*cxpb*I_A_x(Abs(n-2))-1.0*sxpa*sxpb*I_A_x(n)-0.5*sxpa*sxpb*I_A_x(n+2)-0.5*sxpa*sxpb*I_A_x(Abs(n-2)))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then
          imag_An = 0.d0 
        else 
          imag_An  = 0.5*inv_ax2*(cxpa*sxpb+cxpb*sxpa)*(n*I_A_x(n+1)+n*I_A_x(Abs(n-1))+I_A_x(n+1)-I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1200) ! | px    py   s    s     ( 97 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1201) ! | px    py   s    px    ( 98 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1202) ! | px    py   s    py    ( 99 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1203) ! | px    py   s    pz    ( 100 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1210) ! | px    py   px   s     ( 101 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1211) ! | px    py   px   px    ( 102 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        if (dabs(BBx) < 1.d-300) then
          imag_Bn = 0.d0 
        else 
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1212) ! | px    py   px   py    ( 103 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1213) ! | px    py   px   pz    ( 104 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1220) ! | px    py   py   s     ( 105 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1221) ! | px    py   py   px    ( 106 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1222) ! | px    py   py   py    ( 107 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1223) ! | px    py   py   pz    ( 108 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1230) ! | px    py   pz   s     ( 109 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1231) ! | px    py   pz   px    ( 110 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1232) ! | px    py   pz   py    ( 111 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1233) ! | px    py   pz   pz    ( 112 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1300) ! | px    pz   s    s     ( 113 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1301) ! | px    pz   s    px    ( 114 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1302) ! | px    pz   s    py    ( 115 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1303) ! | px    pz   s    pz    ( 116 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1310) ! | px    pz   px   s     ( 117 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1311) ! | px    pz   px   px    ( 118 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        if (dabs(BBx) < 1.d-300) then
          imag_Bn = 0.d0 
        else 
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1312) ! | px    pz   px   py    ( 119 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1313) ! | px    pz   px   pz    ( 120 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1320) ! | px    pz   py   s     ( 121 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1321) ! | px    pz   py   px    ( 122 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1322) ! | px    pz   py   py    ( 123 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1323) ! | px    pz   py   pz    ( 124 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1330) ! | px    pz   pz   s     ( 125 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1331) ! | px    pz   pz   px    ( 126 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1332) ! | px    pz   pz   py    ( 127 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (1333) ! | px    pz   pz   pz    ( 128 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2000) ! | py    s    s    s     ( 129 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2001) ! | py    s    s    px    ( 130 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2002) ! | py    s    s    py    ( 131 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2003) ! | py    s    s    pz    ( 132 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2010) ! | py    s    px   s     ( 133 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2011) ! | py    s    px   px    ( 134 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        if (dabs(BBx) < 1.d-300) then
          imag_Bn = 0.d0 
        else 
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2012) ! | py    s    px   py    ( 135 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2013) ! | py    s    px   pz    ( 136 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2020) ! | py    s    py   s     ( 137 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2021) ! | py    s    py   px    ( 138 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2022) ! | py    s    py   py    ( 139 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2023) ! | py    s    py   pz    ( 140 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2030) ! | py    s    pz   s     ( 141 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2031) ! | py    s    pz   px    ( 142 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2032) ! | py    s    pz   py    ( 143 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2033) ! | py    s    pz   pz    ( 144 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2100) ! | py    px   s    s     ( 145 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2101) ! | py    px   s    px    ( 146 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2102) ! | py    px   s    py    ( 147 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2103) ! | py    px   s    pz    ( 148 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2110) ! | py    px   px   s     ( 149 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2111) ! | py    px   px   px    ( 150 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        if (dabs(BBx) < 1.d-300) then
          imag_Bn = 0.d0 
        else 
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2112) ! | py    px   px   py    ( 151 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2113) ! | py    px   px   pz    ( 152 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2120) ! | py    px   py   s     ( 153 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2121) ! | py    px   py   px    ( 154 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2122) ! | py    px   py   py    ( 155 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2123) ! | py    px   py   pz    ( 156 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2130) ! | py    px   pz   s     ( 157 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2131) ! | py    px   pz   px    ( 158 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2132) ! | py    px   pz   py    ( 159 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2133) ! | py    px   pz   pz    ( 160 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (2200) ! | py    py   s    s     ( 161 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2201) ! | py    py   s    px    ( 162 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2202) ! | py    py   s    py    ( 163 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2203) ! | py    py   s    pz    ( 164 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2210) ! | py    py   px   s     ( 165 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2211) ! | py    py   px   px    ( 166 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        if (dabs(BBx) < 1.d-300) then
          imag_Bn = 0.d0 
        else 
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2212) ! | py    py   px   py    ( 167 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2213) ! | py    py   px   pz    ( 168 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2220) ! | py    py   py   s     ( 169 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2221) ! | py    py   py   px    ( 170 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2222) ! | py    py   py   py    ( 171 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2223) ! | py    py   py   pz    ( 172 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2230) ! | py    py   pz   s     ( 173 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2231) ! | py    py   pz   px    ( 174 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2232) ! | py    py   pz   py    ( 175 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2233) ! | py    py   pz   pz    ( 176 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2300) ! | py    pz   s    s     ( 177 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2301) ! | py    pz   s    px    ( 178 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2302) ! | py    pz   s    py    ( 179 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2303) ! | py    pz   s    pz    ( 180 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2310) ! | py    pz   px   s     ( 181 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2311) ! | py    pz   px   px    ( 182 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        if (dabs(BBx) < 1.d-300) then
          imag_Bn = 0.d0 
        else 
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2312) ! | py    pz   px   py    ( 183 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2313) ! | py    pz   px   pz    ( 184 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2320) ! | py    pz   py   s     ( 185 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2321) ! | py    pz   py   px    ( 186 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2322) ! | py    pz   py   py    ( 187 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2323) ! | py    pz   py   pz    ( 188 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2330) ! | py    pz   pz   s     ( 189 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2331) ! | py    pz   pz   px    ( 190 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (2332) ! | py    pz   pz   py    ( 191 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (2333) ! | py    pz   pz   pz    ( 192 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3000) ! | pz    s    s    s     ( 193 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3001) ! | pz    s    s    px    ( 194 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3002) ! | pz    s    s    py    ( 195 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3003) ! | pz    s    s    pz    ( 196 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3010) ! | pz    s    px   s     ( 197 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3011) ! | pz    s    px   px    ( 198 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        if (dabs(BBx) < 1.d-300) then
          imag_Bn = 0.d0 
        else 
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3012) ! | pz    s    px   py    ( 199 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3013) ! | pz    s    px   pz    ( 200 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3020) ! | pz    s    py   s     ( 201 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3021) ! | pz    s    py   px    ( 202 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3022) ! | pz    s    py   py    ( 203 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3023) ! | pz    s    py   pz    ( 204 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3030) ! | pz    s    pz   s     ( 205 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3031) ! | pz    s    pz   px    ( 206 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3032) ! | pz    s    pz   py    ( 207 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3033) ! | pz    s    pz   pz    ( 208 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3100) ! | pz    px   s    s     ( 209 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3101) ! | pz    px   s    px    ( 210 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3102) ! | pz    px   s    py    ( 211 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3103) ! | pz    px   s    pz    ( 212 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3110) ! | pz    px   px   s     ( 213 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3111) ! | pz    px   px   px    ( 214 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        if (dabs(BBx) < 1.d-300) then
          imag_Bn = 0.d0 
        else 
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3112) ! | pz    px   px   py    ( 215 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3113) ! | pz    px   px   pz    ( 216 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3120) ! | pz    px   py   s     ( 217 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3121) ! | pz    px   py   px    ( 218 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3122) ! | pz    px   py   py    ( 219 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3123) ! | pz    px   py   pz    ( 220 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3130) ! | pz    px   pz   s     ( 221 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3131) ! | pz    px   pz   px    ( 222 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3132) ! | pz    px   pz   py    ( 223 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3133) ! | pz    px   pz   pz    ( 224 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      case (3200) ! | pz    py   s    s     ( 225 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3201) ! | pz    py   s    px    ( 226 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3202) ! | pz    py   s    py    ( 227 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3203) ! | pz    py   s    pz    ( 228 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3210) ! | pz    py   px   s     ( 229 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3211) ! | pz    py   px   px    ( 230 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else 
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3212) ! | pz    py   px   py    ( 231 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3213) ! | pz    py   px   pz    ( 232 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3220) ! | pz    py   py   s     ( 233 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3221) ! | pz    py   py   px    ( 234 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3222) ! | pz    py   py   py    ( 235 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3223) ! | pz    py   py   pz    ( 236 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3230) ! | pz    py   pz   s     ( 237 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3231) ! | pz    py   pz   px    ( 238 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3232) ! | pz    py   pz   py    ( 239 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3233) ! | pz    py   pz   pz    ( 240 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3300) ! | pz    pz   s    s     ( 241 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3301) ! | pz    pz   s    px    ( 242 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3302) ! | pz    pz   s    py    ( 243 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3303) ! | pz    pz   s    pz    ( 244 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3310) ! | pz    pz   px   s     ( 245 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3311) ! | pz    pz   px   px    ( 246 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -0.5*inv_ax2*(-1.0*cxqc*cxqd*I_B_x(n)+0.5*cxqc*cxqd*I_B_x(n+2)+0.5*cxqc*cxqd*I_B_x(Abs(n-2))-1.0*sxqc*sxqd*I_B_x(n)-0.5*sxqc*sxqd*I_B_x(n+2)-0.5*sxqc*sxqd*I_B_x(Abs(n-2)))
        if (dabs(BBx) < 1.d-300) then
          imag_Bn = 0.d0 
        else 
          imag_Bn  = -0.5*inv_ax2*(cxqc*sxqd+cxqd*sxqc)*(n*I_B_x(n+1)+n*I_B_x(Abs(n-1))+I_B_x(n+1)-I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3312) ! | pz    pz   px   py    ( 247 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3313) ! | pz    pz   px   pz    ( 248 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3320) ! | pz    pz   py   s     ( 249 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3321) ! | pz    pz   py   px    ( 250 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3322) ! | pz    pz   py   py    ( 251 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3323) ! | pz    pz   py   pz    ( 252 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3330) ! | pz    pz   pz   s     ( 253 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3331) ! | pz    pz   pz   px    ( 254 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      case (3332) ! | pz    pz   pz   py    ( 255 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case (3333) ! | pz    pz   pz   pz    ( 256 )
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      case default
        F_x(n) = 0.d0

      
      end select

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=--=-=-=-=-=-=-=-=-=-=- !
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=--=-=-=-=-=-=-=-=-=-=- ! 


      result = 0.0d0

      do n = 1, 64
        result = result + w_val(n) * const(n) * S(n) 
      end do

      contains
 
      double precision function S(i_quad) result(sum)

       use bessel_derivatives
       use precomputed_bessel

      implicit none

      ! --------------------------------------------------------------- !

      integer, intent(in)                  :: i_quad

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      sum = 2.d0 * dot_product(F_x(1:Nmax_x(i_quad)), I_C_table_x(1:Nmax_x(i_quad), i_quad))
      sum = sum + F_x(0) * I_C_table_x(0, i_quad)

 end function S


end subroutine integrate_ERI_sum_new