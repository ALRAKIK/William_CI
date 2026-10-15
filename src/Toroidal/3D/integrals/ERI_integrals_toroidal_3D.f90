subroutine ERI_integral_4_function_toroidal_3D(one,two,three,four,value)
      
      use torus_init    
      use classification_ERI
      use constants_module

      implicit none 

      !-----------------------------------------------------------------!

      type(ERI_function),intent(in)   :: one , two , three , four 

      double precision,intent(out)    :: value 

      integer                         ::  i , j   , k  , l
      character(LEN=2)                :: o1 , o2 , o3 , o4 
      
      double precision                :: xa  , ya , za 
      double precision                :: xb  , yb , zb
      double precision                :: xc  , yc , zc
      double precision                :: xd  , yd , zd
      double precision                :: xAB , xCD 
      double precision                :: yAB , yCD
      double precision                :: zAB , zCD
      double precision                :: xp  , xq 
      double precision                :: yp  , yq
      double precision                :: zp  , zq
      double precision                :: c1  , c2 , c3 , c4 
      double precision                :: alpha , beta , gamma , delta 
      double precision                :: const  
      double precision                :: value_s 
      double precision                :: mu_x, mu_y , mu_z
      double precision                :: nu_x, nu_y , nu_z
      double precision                :: mu  , nu 
      double precision                :: xpA , xpB , xqC , xqD , phix
      double precision                :: ypA , ypB , yqC , yqD , phiy
      double precision                :: zpA , zpB , zqC , zqD , phiz
      double precision                :: const_x , const_y , const_z 
      double precision                :: inv_ax2, inv_ay2, inv_az2
      integer                         :: pattern_id, encode_orbital_pattern


      inv_ax2 = 1.0d0/(ax*ax)
      inv_ay2 = 1.0d0/(ay*ay) 
      inv_az2 = 1.0d0/(az*az)

      !-----------------------------------------------------------------!

      value = 0.d0 
            
      xa =   one%x ; ya =   one%y ; za =   one%z 
      xb =   two%x ; yb =   two%y ; zb =   two%z
      xc = three%x ; yc = three%y ; zc = three%z
      xd =  four%x ; yd =  four%y ; zd =  four%z


      xAB = xa - xb ; xCD = xc - xd
      yAB = ya - yb ; yCD = yc - yd
      zAB = za - zb ; zCD = zc - zd 

      do i = 1 , size  (one%exponent)
        alpha = one%exponent(i)
        c1    = one%coefficient(i)
        o1    = one%orbital
        do j = 1 , size  (two%exponent)
          beta  = two%exponent(j)
          c2    = two%coefficient(j)
          o2    = two%orbital

          call bary_exponent_x(alpha,beta,XAB,mu_x)
          call bary_exponent_y(alpha,beta,YAB,mu_y)
          call bary_exponent_z(alpha,beta,ZAB,mu_z)

          call bary_center_toroidal_x(alpha,beta,xa,xb,xp)
          call bary_center_toroidal_y(alpha,beta,ya,yb,yp)
          call bary_center_toroidal_z(alpha,beta,za,zb,zp)

          mu = alpha+beta 
          
          do k = 1 , size(three%exponent)
            gamma = three%exponent(k)
            c3    = three%coefficient(k)
            o3    = three%orbital
            do l = 1 , size (four%exponent)
              delta = four%exponent(l)
              c4    = four%coefficient(l)
              o4    = four%orbital

              call bary_exponent_x(gamma,delta,XCD,nu_x)
              call bary_exponent_y(gamma,delta,YCD,nu_y)
              call bary_exponent_z(gamma,delta,ZCD,nu_z)

              call bary_center_toroidal_x(gamma,delta,xc,xd,xq)
              call bary_center_toroidal_y(gamma,delta,yc,yd,yq)
              call bary_center_toroidal_z(gamma,delta,zc,zd,zq)
               
              nu    = gamma + delta

              const   = (c1*c2*c3*c4) * 2.d0 /dsqrt(pi) * (Lx*Lx) * (Ly*Ly) * (Lz*Lz)  

              const_x = dexp(2.d0 * ( mu_x + nu_x - mu - nu )* inv_ax2 )
              const_y = dexp(2.d0 * ( mu_y + nu_y - mu - nu )* inv_ay2 )
              const_z = dexp(2.d0 * ( mu_z + nu_z - mu - nu )* inv_az2 )

              if ( dabs(const_x * const_y * const_z * const) < 2.22d-16 ) cycle

              xpA     = ax*(xp - xa) ; ypA     = ay*(yp - ya) ; zpA     = az*(zp - za)
              xpB     = ax*(xp - xb) ; ypB     = ay*(yp - yb) ; zpB     = az*(zp - zb) 
              xqC     = ax*(xq - xc) ; yqC     = ay*(yq - yc) ; zqC     = az*(zq - zc)
              xqD     = ax*(xq - xd) ; yqD     = ay*(yq - yd) ; zqD     = az*(zq - zd)
              phix    = ax*(xp - xq) ; phiy    = ay*(yp - yq) ; phiz    = az*(zp - zq)

              pattern_id = encode_orbital_pattern(o1, o2, o3, o4)
 
              call integrate_ERI_3D(pattern_id      ,mu_x,nu_x,phix,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,&
                                                     mu_y,nu_y,phiy,ypA,ypB,yqC,yqD,ya,yb,yc,yd,yp,yq,&
                                                     mu_z,nu_z,phiz,zpA,zpB,zqC,zqD,za,zb,zc,zd,zp,zq,value_s)

              value  = value    + const * const_x * const_y * const_z * value_s

            end do
          end do 
        end do 
      end do 

end subroutine ERI_integral_4_function_toroidal_3D

subroutine integrate_ERI_3D(pattern_id    ,p_x,q_x,phi_x,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,&
                                           p_y,q_y,phi_y,ypA,ypB,yqC,yqD,ya,yb,yc,yd,yp,yq,&
                                           p_z,q_z,phi_z,zpA,zpB,zqC,zqD,za,zb,zc,zd,zp,zq,result)
      
      use quadpack , only : dqagi , dqag , dqags
      use iso_c_binding
      use torus_init
      !use gauss_legendre_quadrature
      use table_1d_lookup
      use constants_module

      use, intrinsic :: ieee_arithmetic

      implicit none

      ! Input parameters

      double precision, intent(in)       :: p_x , q_x 
      double precision, intent(in)       :: p_y , q_y
      double precision, intent(in)       :: p_z , q_z 
      double precision, intent(in)       :: phi_x , phi_y , phi_z
      double precision, intent(in)       :: xpA , xpB , xqC , xqD
      double precision, intent(in)       :: ypA , ypB , yqC , yqD
      double precision, intent(in)       :: zpA , zpB , zqC , zqD
      double precision, intent(in)       :: xa , xb , xc ,xd ,xp ,xq 
      double precision, intent(in)       :: ya , yb , yc ,yd ,yp ,yq
      double precision, intent(in)       :: za , zb , zc ,zd ,zp ,zq
      integer         , intent(in)       :: pattern_id
      

      ! Output parameters

      double precision, intent(out)      :: result
    
      ! Local variables

      double precision,parameter         :: epsabs = 1.d-10 , epsrel = 1.d-8
      integer,parameter                  :: inf = 1 
      double precision,parameter         :: bound = 0.0d0
      integer, parameter                 :: limit = 50
      integer, parameter                 :: lenw = limit*4
      integer                            :: ier, iwork(limit), last, neval
      double precision                   :: abserr, work(lenw)
      double precision,parameter         :: pi2  = pi * pi


      double precision                   :: inv_ax  , inv_ay  , inv_az 
      double precision                   :: inv_ax2 , inv_ay2 , inv_az2 

      double precision                   :: cxpa    , sxpa
      double precision                   :: cypa    , sypa
      double precision                   :: czpa    , szpa

      double precision                   :: cxpb    , sxpb
      double precision                   :: cypb    , sypb
      double precision                   :: czpb    , szpb

      double precision                   :: cxqc    , sxqc
      double precision                   :: cyqc    , syqc
      double precision                   :: czqc    , szqc

      double precision                   :: cxqd    , sxqd
      double precision                   :: cyqd    , syqd
      double precision                   :: czqd    , szqd

      double precision                   :: c2xpab  , s2xpab
      double precision                   :: c2ypab  , s2ypab
      double precision                   :: c2zpab  , s2zpab

      double precision                   :: c2xqcd  , s2xqcd
      double precision                   :: c2yqcd  , s2yqcd
      double precision                   :: c2zqcd  , s2zqcd

      double precision                   :: AAx,  BBx, CCx
      double precision                   :: AAy,  BBy, CCy
      double precision                   :: AAz,  BBz, CCz

      integer                            :: Nmax_AAx , Nmax_BBx
      integer                            :: Nmax_AAy , Nmax_BBy
      integer                            :: Nmax_AAz , Nmax_BBz

       ! --------------------------------------------------------------- !
      double precision                     :: I_A_x(0:50000)
      double precision                     :: I_B_x(0:50000)
      double precision                     :: I_A_y(0:50000)
      double precision                     :: I_B_y(0:50000)
      double precision                     :: I_A_z(0:50000)
      double precision                     :: I_B_z(0:50000)

      ! COMPLEX(KIND=KIND(1.0D0)), PARAMETER :: I_dp = (0.0D0, 1.0D0)
      !COMPLEX(KIND=KIND(1.0D0))            :: expo_x(0:50000)
      !COMPLEX(KIND=KIND(1.0D0))            :: expo_y(0:50000)
      !COMPLEX(KIND=KIND(1.0D0))            :: expo_z(0:50000)
      ! integer                              :: i 

      !COMPLEX(KIND=KIND(1.0D0))            :: F_x(0:50000)
      !COMPLEX(KIND=KIND(1.0D0))            :: F_y(0:50000)
      !COMPLEX(KIND=KIND(1.0D0))            :: F_z(0:50000)

      ! double precision                     :: F_x(0:50000)
      ! double precision                     :: F_y(0:50000)
      ! double precision                     :: F_z(0:50000)

      ! double precision                     :: cos_phi_x(0:50000)
      ! double precision                     :: cos_phi_y(0:50000)
      ! double precision                     :: cos_phi_z(0:50000)
      ! --------------------------------------------------------------- !

      inv_ax  = 1.d0/ax 
      inv_ay  = 1.d0/ay 
      inv_az  = 1.d0/az 

      inv_ax2 = inv_ax * inv_ax
      inv_ay2 = inv_ay * inv_ay
      inv_az2 = inv_az * inv_az

      AAx   = inv_ax2 * 2.d0 * p_x
      BBx   = inv_ax2 * 2.d0 * q_x
      AAy   = inv_ay2 * 2.d0 * p_y
      BBy   = inv_ay2 * 2.d0 * q_y
      AAz   = inv_az2 * 2.d0 * p_z 
      BBz   = inv_az2 * 2.d0 * q_z

      ! --------------------------------------------------------------- !

      Nmax_AAx = get_Nmax_from_bessel_1D(AAx)
      Nmax_BBx = get_Nmax_from_bessel_1D(BBx)
 
      call bessel_I_scaled_backward(Nmax_AAx+20, AAx, I_A_x)
      call bessel_I_scaled_backward(Nmax_BBx+20, BBx, I_B_x)

      Nmax_AAy = get_Nmax_from_bessel_1D(AAy)
      Nmax_BBy = get_Nmax_from_bessel_1D(BBy)
 
      call bessel_I_scaled_backward(Nmax_AAy+20, AAy, I_A_y)
      call bessel_I_scaled_backward(Nmax_BBy+20, BBy, I_B_y)

      Nmax_AAz = get_Nmax_from_bessel_1D(AAz)
      Nmax_BBz = get_Nmax_from_bessel_1D(BBz)
 
      call bessel_I_scaled_backward(Nmax_AAz+20, AAz, I_A_z)
      call bessel_I_scaled_backward(Nmax_BBz+20, BBz, I_B_z)

      ! --------------------------------------------------------------- !

      cxpa  = dcos(xpA) ; cypa  = dcos(ypA) ; czpa  = dcos(zpA) 
      sxpa  = dsin(xpA) ; sypa  = dsin(ypA) ; szpa  = dsin(zpA) 

      cxpb  = dcos(xpB) ; cypb  = dcos(ypB) ; czpb  = dcos(zpB) 
      sxpb  = dsin(xpB) ; sypb  = dsin(ypB) ; szpb  = dsin(zpB) 
      
      cxqc  = dcos(xqC) ; cyqc  = dcos(yqC) ; czqc  = dcos(zqC) 
      sxqc  = dsin(xqC) ; syqc  = dsin(yqC) ; szqc  = dsin(zqC) 

      cxqd  = dcos(xqD) ; cyqd  = dcos(yqD) ; czqd  = dcos(zqD)
      sxqd  = dsin(xqD) ; syqd  = dsin(yqD) ; szqd  = dsin(zqD)

      c2xpab = dcos(ax*(2.d0*xp-xA-xB)) ; c2ypab = dcos(ay*(2.d0*yp-yA-yB)) ; c2zpab = dcos(az*(2.d0*zp-zA-zB))
      s2xpab = dsin(ax*(2.d0*xp-xA-xB)) ; s2ypab = dsin(ay*(2.d0*yp-yA-yB)) ; s2zpab = dsin(az*(2.d0*zp-zA-zB))

      c2xqcd = dcos(ax*(2.d0*xq-xC-xD)) ; c2yqcd = dcos(ay*(2.d0*yq-yC-yD)) ; c2zqcd = dcos(az*(2.d0*zq-zC-zD)) 
      s2xqcd = dsin(ax*(2.d0*xq-xC-xD)) ; s2yqcd = dsin(ay*(2.d0*yq-yC-yD)) ; s2zqcd = dsin(az*(2.d0*zq-zC-zD)) 
      
      !call gauss_legendre_64(transformed_integrand, 0.0d0, 1.0d0, result)

      ! if (ier /= 0) then
      !   write(*,'(A,I8,A)') 'Error code = ', ier
      ! end if

      contains

      function f_decay(t,i_quad) result(ft)
        double precision, intent(in) :: t
        integer, intent(in)          :: i_quad
        double precision             :: ft

        ft    =  S(t,i_quad)

      end function f_decay

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

      use bessel_functions
      use table_lookup_module
      use precomputed_bessel

      implicit none
      double precision, intent(in)         :: t
      integer, intent(in)                  :: i_quad
      COMPLEX(KIND=KIND(1.0D0)), PARAMETER :: I_dp = (0.0D0, 1.0D0)
      integer                              :: n  
      COMPLEX(KIND=KIND(1.0D0))            :: termAn , termBn
      !double precision                     :: termAn , termBn
      double precision                     :: termC 
      double precision                     :: sum1, sum2, sum3
      COMPLEX(KIND=KIND(1.0D0))            :: term

      ! --------------------------------------------------------------- !
      
      integer                              :: Nmax_x, Nmax_y, Nmax_z
      logical                              :: no_converged_x = .true. , no_converged_y = .true. , no_converged_z = .true.
      COMPLEX(KIND=KIND(1.0D0))            :: expo_term_x , expo_term_y , expo_term_z
      COMPLEX(KIND=KIND(1.0D0))            :: current_term_x , current_term_y , current_term_z
      double precision,parameter           :: eps = 2.22d-16


      CCx   = inv_ax2 * 2.d0 * t * t 
      CCy   = inv_ay2 * 2.d0 * t * t 
      CCz   = inv_az2 * 2.d0 * t * t 

      Nmax_x = get_Nmax(AAx, BBx, CCx, Phi_x) + 10
      Nmax_y = get_Nmax(AAy, BBy, CCy, Phi_y) + 10
      Nmax_z = get_Nmax(AAz, BBz, CCz, Phi_z) + 10


      expo_term_x    = exp(I_dp*phi_x)
      expo_term_y    = exp(I_dp*phi_y)
      expo_term_z    = exp(I_dp*phi_z)

      current_term_x = expo_term_x
      current_term_y = expo_term_y
      current_term_z = expo_term_z

      termAn    = 0.d0 
      termBn    = 0.d0 
      termC     = 0.d0

      sum1      = 0.d0
      sum2      = 0.d0
      sum3      = 0.d0
      sum       = 0.d0

      select case(pattern_id)

      case (0000) ! | s   s   s   s    ( 1 ) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0001) ! | s   s   s   px   ( 2 ) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0002) ! | s   s   s   py   ( 3 ) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0003) ! | s   s   s   pz   ( 4 ) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0010) ! | s   s   px  s    ( 5 ) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0011) ! | s   s   px  px   ( 6 ) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0012) ! | s   s   px  py   ( 7 ) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0013) ! | s   s   px  pz   ( 8 ) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0020) ! | s   s   py  s    ( 9 ) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0021) ! | s   s   py  px   ( 10) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0022) ! | s   s   py  py   ( 11) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0023) ! | s   s   py  pz   ( 12) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0030) ! | s   s   pz  s    ( 13) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0031) ! | s   s   pz  px   ( 14) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0032) ! | s   s   pz  py   ( 15) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0033) ! | s   s   pz  pz   ( 16) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0100) ! | s   px  s   s    ( 17) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0101) ! | s   px  s   px   ( 18) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0102) ! | s   px  s   py   ( 19) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0103) ! | s   px  s   pz   ( 20) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0110) ! | s   px  px  s    ( 21) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0111) ! | s   px  px  px   ( 22) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0112) ! | s   px  px  py   ( 23) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0113) ! | s   px  px  pz   ( 24) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0120) ! | s   px  py  s    ( 25) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0121) ! | s   px  py  px   ( 26) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0122) ! | s   px  py  py   ( 27) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0123) ! | s   px  py  pz   ( 28) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0130) ! | s   px  pz  s    ( 29) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0131) ! | s   px  pz  px   ( 30) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0132) ! | s   px  pz  py   ( 31) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0133) ! | s   px  pz  pz   ( 32) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0200) ! | s   py  s   s    ( 33) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0201) ! | s   py  s   px   ( 34) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0202) ! | s   py  s   py   ( 35) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0203) ! | s   py  s   pz   ( 36) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0210) ! | s   py  px  s    ( 37) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0211) ! | s   py  px  px   ( 38) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0212) ! | s   py  px  py   ( 39) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0213) ! | s   py  px  pz   ( 40) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0220) ! | s   py  py  s    ( 41) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0221) ! | s   py  py  px   ( 42) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0222) ! | s   py  py  py   ( 43) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0223) ! | s   py  py  pz   ( 44) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0230) ! | s   py  pz  s    ( 45) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0231) ! | s   py  pz  px   ( 46) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0232) ! | s   py  pz  py   ( 47) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0233) ! | s   py  pz  pz   ( 48) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0300) ! | s   pz  s   s    ( 49) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0301) ! | s   pz  s   px   ( 50) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0302) ! | s   pz  s   py   ( 51) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0303) ! | s   pz  s   pz   ( 52) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0310) ! | s   pz  px  s    ( 53) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0311) ! | s   pz  px  px   ( 54) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0312) ! | s   pz  px  py   ( 55) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0313) ! | s   pz  px  pz   ( 56) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0320) ! | s   pz  py  s    ( 57) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0321) ! | s   pz  py  px   ( 58) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0322) ! | s   pz  py  py   ( 59) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0323) ! | s   pz  py  pz   ( 60) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0330) ! | s   pz  pz  s    ( 61) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0331) ! | s   pz  pz  px   ( 62) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0332) ! | s   pz  pz  py   ( 63) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (0333) ! | s   pz  pz  pz   ( 64) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1000) ! | px  s   s   s    ( 65) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1001) ! | px  s   s   px   ( 66) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1002) ! | px  s   s   py   ( 67) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1003) ! | px  s   s   pz   ( 68) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1010) ! | px  s   px  s    ( 69) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1011) ! | px  s   px  px   ( 70) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1012) ! | px  s   px  py   ( 71) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1013) ! | px  s   px  pz   ( 72) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1020) ! | px  s   py  s    ( 73) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1021) ! | px  s   py  px   ( 74) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1022) ! | px  s   py  py   ( 75) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1023) ! | px  s   py  pz   ( 76) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1030) ! | px  s   pz  s    ( 77) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1031) ! | px  s   pz  px   ( 78) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1032) ! | px  s   pz  py   ( 79) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1033) ! | px  s   pz  pz   ( 80) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1100) ! | px  px  s   s    ( 81) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
        end if 
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1101) ! | px  px  s   px   ( 82) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1102) ! | px  px  s   py   ( 83) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1103) ! | px  px  s   pz   ( 84) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1110) ! | px  px  px  s    ( 85) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1111) ! | px  px  px  px   ( 86) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1112) ! | px  px  px  py   ( 87) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300 ) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1113) ! | px  px  px  pz   ( 88) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1120) ! | px  px  py  s    ( 89) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1121) ! | px  px  py  px   ( 90) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1122) ! | px  px  py  py   ( 91) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1123) ! | px  px  py  pz   ( 92) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1130) ! | px  px  pz  s    ( 93) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1131) ! | px  px  pz  px   ( 94) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300 ) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1132) ! | px  px  pz  py   ( 95) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1133) ! | px  px  pz  pz   ( 96) 

      n         = 0
      sum1      = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        if (dabs(AAx) < 1.d-300) then 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2))))
        else 
        termAn  = inv_ax2 * (cxpa * cxpb * I_A_x(n)-c2xpab*(0.25d0*(I_A_x(abs(n-2))+2.d0*I_A_x(n)+I_A_x(n+2)))+I_dp/AAx*n * s2xpab * (0.5d0*(I_A_x(abs(n-1))+I_A_x(n+1))))
         end if 
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1200) ! | px  py  s   s    ( 97) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1201) ! | px  py  s   px   ( 98) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1202) ! | px  py  s   py   ( 99) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1203) ! | px  py  s   pz   ( 100) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1210) ! | px  py  px  s    ( 101) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1211) ! | px  py  px  px   ( 102) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1212) ! | px  py  px  py   ( 103) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1)) - I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1)) - I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1213) ! | px  py  px  pz   ( 104) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1220) ! | px  py  py  s    ( 105) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1221) ! | px  py  py  px   ( 106) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1)) - I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1)) - I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1222) ! | px  py  py  py   ( 107) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1223) ! | px  py  py  pz   ( 108) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1230) ! | px  py  pz  s    ( 109) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1231) ! | px  py  pz  px   ( 110) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1232) ! | px  py  pz  py   ( 111) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1233) ! | px  py  pz  pz   ( 112) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1300) ! | px  pz  s   s    ( 113) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1301) ! | px  pz  s   px   ( 114) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1302) ! | px  pz  s   py   ( 115) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1303) ! | px  pz  s   pz   ( 116) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1310) ! | px  pz  px  s    ( 117) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1311) ! | px  pz  px  px   ( 118) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1312) ! | px  pz  px  py   ( 119) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1313) ! | px  pz  px  pz   ( 120) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1320) ! | px  pz  py  s    ( 121) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1321) ! | px  pz  py  px   ( 122) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1322) ! | px  pz  py  py   ( 123) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1323) ! | px  pz  py  pz   ( 124) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1330) ! | px  pz  pz  s    ( 125) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1331) ! | px  pz  pz  px   ( 126) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1332) ! | px  pz  pz  py   ( 127) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (1333) ! | px  pz  pz  pz   ( 128) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpa * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpa * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2000) ! | py  s   s   s    ( 129) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2001) ! | py  s   s   px   ( 130) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2002) ! | py  s   s   py   ( 131) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2003) ! | py  s   s   pz   ( 132) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2010) ! | py  s   px  s    ( 133) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2011) ! | py  s   px  px   ( 134) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2012) ! | py  s   px  py   ( 135) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2013) ! | py  s   px  pz   ( 136) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2020) ! | py  s   py  s    ( 137) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)   
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2021) ! | py  s   py  px   ( 138) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)   
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2022) ! | py  s   py  py   ( 139) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2023) ! | py  s   py  pz   ( 140) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)   
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2030) ! | py  s   pz  s    ( 141) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2031) ! | py  s   pz  px   ( 142) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2032) ! | py  s   pz  py   ( 143) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2033) ! | py  s   pz  pz   ( 144) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2100) ! | py  px  s   s    ( 145) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2101) ! | py  px  s   px   ( 146) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2102) ! | py  px  s   py   ( 147) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2103) ! | py  px  s   pz   ( 148) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2110) ! | py  px  px  s    ( 149) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2111) ! | py  px  px  px   ( 150) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2112) ! | py  px  px  py   ( 151) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2113) ! | py  px  px  pz   ( 152) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2120) ! | py  px  py  s    ( 153) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)   
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2121) ! | py  px  py  px   ( 154) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)   
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2122) ! | py  px  py  py   ( 155) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2123) ! | py  px  py  pz   ( 156) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)   
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2130) ! | py  px  pz  s    ( 157) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2131) ! | py  px  pz  px   ( 158) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2132) ! | py  px  pz  py   ( 159) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2133) ! | py  px  pz  pz   ( 160) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2200) ! | py  py  s   s    ( 161) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if 
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2201) ! | py  py  s   px   ( 162) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if        
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2202) ! | py  py  s   py   ( 163) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
      if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2203) ! | py  py  s   pz   ( 164) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2210) ! | py  py  px  s    ( 165) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2211) ! | py  py  px  px   ( 166) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2212) ! | py  py  px  py   ( 167) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2213) ! | py  py  px  pz   ( 168) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2220) ! | py  py  py  s    ( 169) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2221) ! | py  py  py  px   ( 170) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2222) ! | py  py  py  py   ( 171) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2223) ! | py  py  py  pz   ( 172) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2230) ! | py  py  pz  s    ( 173) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2231) ! | py  py  pz  px   ( 174) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2232) ! | py  py  pz  py   ( 175) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2233) ! | py  py  pz  pz   ( 176) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        if (dabs(AAy) < 1.d-300) then 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2))))
        else 
        termAn  = inv_ay2 * (cypa*cypb*I_A_y(n)- c2ypab * (0.25d0*(I_A_y(abs(n-2))+2.d0*I_A_y(n)+I_A_y(n+2)))+I_dp/AAy*n * s2ypab * (0.5d0*(I_A_y(abs(n-1))+I_A_y(n+1))))
         end if
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = I_A_z(n) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = I_A_z(n)
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2300) ! | py  pz  s   s    ( 177) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2301) ! | py  pz  s   px   ( 178) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2302) ! | py  pz  s   py   ( 179) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2303) ! | py  pz  s   pz   ( 180) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2310) ! | py  pz  px  s    ( 181) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2311) ! | py  pz  px  px   ( 182) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2312) ! | py  pz  px  py   ( 183) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2313) ! | py  pz  px  pz   ( 184) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2320) ! | py  pz  py  s    ( 185) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)   
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2321) ! | py  pz  py  px   ( 186) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)   
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2322) ! | py  pz  py  py   ( 187) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2323) ! | py  pz  py  pz   ( 188) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)   
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2330) ! | py  pz  pz  s    ( 189) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2331) ! | py  pz  pz  px   ( 190) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2332) ! | py  pz  pz  py   ( 191) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpb * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (2333) ! | py  pz  pz  pz   ( 192) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypa * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypa * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpb * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpb *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpb * (I_A_z(abs(n-1))+I_A_z(n+1)))
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3000) ! | pz  s   s   s    ( 193) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3001) ! | pz  s   s   px   ( 194) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3002) ! | pz  s   s   py   ( 195) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3003) ! | pz  s   s   pz   ( 196) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3010) ! | pz  s   px  s    ( 197) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3011) ! | pz  s   px  px   ( 198) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3012) ! | pz  s   px  py   ( 199) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3013) ! | pz  s   px  pz   ( 200) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3020) ! | pz  s   py  s    ( 201) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3021) ! | pz  s   py  px   ( 202) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3022) ! | pz  s   py  py   ( 203) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3023) ! | pz  s   py  pz   ( 204) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3030) ! | pz  s   pz  s    ( 205) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa  * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3031) ! | pz  s   pz  px   ( 206) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3032) ! | pz  s   pz  py   ( 207) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3033) ! | pz  s   pz  pz   ( 208) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3100) ! | pz  px  s   s    ( 209) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3101) ! | pz  px  s   px   ( 210) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3102) ! | pz  px  s   py   ( 211) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3103) ! | pz  px  s   pz   ( 212) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3110) ! | pz  px  px  s    ( 213) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3111) ! | pz  px  px  px   ( 214) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3112) ! | pz  px  px  py   ( 215) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3113) ! | pz  px  px  pz   ( 216) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3120) ! | pz  px  py  s    ( 217) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3121) ! | pz  px  py  px   ( 218) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3122) ! | pz  px  py  py   ( 219) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3123) ! | pz  px  py  pz   ( 220) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3130) ! | pz  px  pz  s    ( 221) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3131) ! | pz  px  pz  px   ( 222) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1)) - I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3132) ! | pz  px  pz  py   ( 223) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3133) ! | pz  px  pz  pz   ( 224) 

      n         = 0
      sum1      = inv_ax * 0.5d0 * (sxpb * (I_A_x(abs(n-1))+I_A_x(n+1))) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (I_A_x(abs(n-1)) - I_A_x(n+1)) + sxpb * (I_A_x(abs(n-1))+I_A_x(n+1)))
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3200) ! | pz  py  s   s    ( 225) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3201) ! | pz  py  s   px   ( 226) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3202) ! | pz  py  s   py   ( 227) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3203) ! | pz  py  s   pz   ( 228) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3210) ! | pz  py  px  s    ( 229) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3211) ! | pz  py  px  px   ( 230) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3212) ! | pz  py  px  py   ( 231) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3213) ! | pz  py  px  pz   ( 232) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3220) ! | pz  py  py  s    ( 233) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3221) ! | pz  py  py  px   ( 234) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3222) ! | pz  py  py  py   ( 235) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3223) ! | pz  py  py  pz   ( 236) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3230) ! | pz  py  pz  s    ( 237) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3231) ! | pz  py  pz  px   ( 238) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3232) ! | pz  py  pz  py   ( 239) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * ( I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * ( I_dp * czpa * (I_A_z(abs(n-1))-I_A_z(n+1)) + szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3233) ! | pz  py  pz  pz   ( 240) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = inv_ay * 0.5d0 * (sypb * (I_A_y(abs(n-1))+I_A_y(n+1))) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (I_A_y(abs(n-1))-I_A_y(n+1)) + sypb * (I_A_y(abs(n-1))+I_A_y(n+1)))
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az * 0.5d0 * (szpa * (I_A_z(abs(n-1))+I_A_z(n+1))) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        termAn  = inv_az * 0.5d0 * (I_dp * czpa *(I_A_z(abs(n-1))-I_A_z(n+1))+ szpa * (I_A_z(abs(n-1))+I_A_z(n+1)))
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3300) ! | pz  pz  s   s    ( 241) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if 
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3301) ! | pz  pz  s   px   ( 242) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3302) ! | pz  pz  s   py   ( 243) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3303) ! | pz  pz  s   pz   ( 244) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3310) ! | pz  pz  px  s    ( 245) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3311) ! | pz  pz  px  px   ( 246) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        if (dabs(BBx) < 1.d-300) then 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2))))
        else 
        termBn  = inv_ax2 * (cxqc * cxqd * I_B_x(n)-c2xqcd*(0.25d0*(I_B_x(abs(n-2))+2.d0*I_B_x(n)+I_B_x(n+2)))-I_dp/BBx*n * s2xqcd * (0.5d0*(I_B_x(abs(n-1))+I_B_x(n+1))))
         end if 
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3312) ! | pz  pz  px  py   ( 247) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3313) ! | pz  pz  px  pz   ( 248) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqc * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqc * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3320) ! | pz  pz  py  s    ( 249) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3321) ! | pz  pz  py  px   ( 250) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3322) ! | pz  pz  py  py   ( 251) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        if (dabs(BBy) < 1.d-300) then 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2))))
        else 
        termBn  = inv_ay2 * (cyqc*cyqd*I_B_y(n)- c2yqcd * (0.25d0*(I_B_y(abs(n-2))+2.d0*I_B_y(n)+I_B_y(n+2)))-I_dp/BBy*n * s2yqcd * (0.5d0*(I_B_y(abs(n-1))+I_B_y(n+1))))
         end if 
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * I_B_z(n) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = I_B_z(n)
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3323) ! | pz  pz  py  pz   ( 252) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqc * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqc * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * inv_az * 0.5d0 * (szqd * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqd * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3330) ! | pz  pz  pz  s    ( 253) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3331) ! | pz  pz  pz  px   ( 254) 

      n         = 0
      sum1      = I_A_x(n) * inv_ax * 0.5d0 * (sxqd * (I_B_x(abs(n-1))+I_B_x(n+1))) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (I_B_x(abs(n-1))-I_B_x(n+1)) + sxqd * (I_B_x(abs(n-1))+I_B_x(n+1)))
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3332) ! | pz  pz  pz  py   ( 255) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * inv_ay * 0.5d0 * (syqd * (I_B_y(abs(n-1))+I_B_y(n+1))) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (I_B_y(abs(n-1))-I_B_y(n+1)) + syqd * (I_B_y(abs(n-1))+I_B_y(n+1)))
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * inv_az * 0.5d0 * (szqc * (I_B_z(abs(n-1))+I_B_z(n+1))) * I_C_table_z(n, i_quad)  
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (I_B_z(abs(n-1))-I_B_z(n+1)) + szqc * (I_B_z(abs(n-1))+I_B_z(n+1)))
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do

      case (3333) ! | pz  pz  pz  pz   ( 256) 

      n         = 0
      sum1      = I_A_x(n) * I_B_x(n) * I_C_table_x(n, i_quad)  
      do n      = 1 , Nmax_x
        termAn  = I_A_x(n)
        termBn  = I_B_x(n)
        termc   = I_C_table_x(n, i_quad)  
        term    = current_term_x * termC * termAn * termBn
        sum1    = sum1 + 2.d0 * real(term)
        current_term_x = current_term_x * expo_term_x
        if (abs(term) < eps * dabs(sum1) ) then
          no_converged_x = .false.
          exit
        end if 
      end do

      n         = 0
      sum2      = I_A_y(n) * I_B_y(n) * I_C_table_y(n, i_quad)  
      do n      = 1 , Nmax_y
        termAn  = I_A_y(n)
        termBn  = I_B_y(n)
        termc   = I_C_table_y(n, i_quad)  
        term    = current_term_y * termC * termAn * termBn
        sum2    = sum2 + 2.d0 * real(term)
        current_term_y = current_term_y * expo_term_y
        if (abs(term) < eps * dabs(sum2) ) then
          no_converged_y = .false.
          exit
        end if 
      end do

      n         = 0
      sum3      = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))) * inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))) * I_C_table_z(n, i_quad) 
      do n      = 1 , Nmax_z
        if (dabs(AAz) < 1.d-300) then 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2))))
        else 
        termAn  = inv_az2 * (czpa*czpb*I_A_z(n)-c2zpab*(0.25d0*(I_A_z(abs(n-2))+2.d0*I_A_z(n)+I_A_z(n+2)))+I_dp/AAz*n * s2zpab * (0.5d0*(I_A_z(abs(n-1))+I_A_z(n+1))))
         end if
        if (dabs(BBz) < 1.d-300) then 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2))))
        else 
        termBn  = inv_az2 * (czqc*czqd*I_B_z(n)-c2zqcd*(0.25d0*(I_B_z(abs(n-2))+2.d0*I_B_z(n)+I_B_z(n+2)))-I_dp/BBz*n*s2zqcd*(0.5d0*(I_B_z(abs(n-1))+I_B_z(n+1))))
         end if 
        termc   = I_C_table_z(n, i_quad)  
        term    = current_term_z * termC * termAn * termBn
        sum3    = sum3 + 2.d0 * real(term)
        current_term_z = current_term_z * expo_term_z
        if (abs(term) < eps * dabs(sum3) ) then
          no_converged_z = .false.
          exit
        end if 
      end do


      case default
        sum = 0.0d0
      end select


      sum = sum1 * sum2 * sum3

      if ( (no_converged_x) .or. (no_converged_y) .or.  (no_converged_z) ) then
        print *, "Warning: The sum did not converge in x " ,no_converged_x, " or/and y " , no_converged_y, " or/and z " , no_converged_z, ", sum is " , sum, " for the integral " , pattern_id
        print *, "                                       " ,sum1, " or/and y " , sum2, " or/and z " , sum3
      end if

      end function S

end subroutine integrate_ERI_3D