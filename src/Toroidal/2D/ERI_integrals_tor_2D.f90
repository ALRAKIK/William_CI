subroutine ERI_integral_4_function_toroidal_2D(one,two,three,four,value)
      
      use torus_init    
      use classification_ERI
      use HeavisideModule

      implicit none 

      !-----------------------------------------------------------------!

      type(ERI_function),intent(in)   :: one , two , three , four 

      double precision,intent(out)    :: value 

      integer                         ::  i , j   , k  , l
      character(LEN=2)                :: o1 , o2 , o3 , o4 
      double precision,parameter      :: pi     = 3.14159265358979323846D00
      double precision,parameter      :: R2PI52 = 2.0d0*pi**(5.0d0/2.0d0)
      
      double precision                :: xa  , ya 
      double precision                :: xb  , yb  
      double precision                :: xc  , yc , zc
      double precision                :: xd  , yd , zd
      double precision                :: xAB , xCD 
      double precision                :: yAB , yCD
      double precision                :: xp  , xq 
      double precision                :: yp  , yq
      double precision                :: c1  , c2 , c3 , c4 
      double precision                :: alpha , beta , gamma , delta 
      double precision                :: const  
      double precision                :: value_s 
      double precision                :: eta = 1e-9
      double precision                :: mu_x, mu_y
      double precision                :: nu_x, nu_y 
      double precision                :: mu  , nu 



      !-----------------------------------------------------------------!



      value = 0.d0 
            
      xa =   one%x ; ya =   one%y 
      xb =   two%x ; yb =   two%y  
      xc = three%x ; yc = three%y ; zc = three%z
      xd =  four%x ; yd =  four%y ; zd =  four%z


      xAB = xa - xb ; xCD = xc - xd
      yAB = ya - yb ; yCD = yc - yd

      do i = 1 , size  (one%exponent)
        alpha = one%exponent(i)
        c1    = one%coefficient(i)
        o1    = one%orbital
        do j = 1 , size  (two%exponent)
          beta  = two%exponent(j)
          c2    = two%coefficient(j)
          o2    = two%orbital


          mu_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(XAB)))+eta
          mu_y  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(YAB)))+eta


          xp    = datan((alpha*dsin(ax*xa)+beta*dsin(ax*xb))/(alpha*dcos(ax*xa)+beta*dcos(ax*xb)))/ax + 0.5*Lx * Heaviside(-alpha*dcos(ax*xa)-beta*dcos(ax*xb))  
          yp    = datan((alpha*dsin(ay*ya)+beta*dsin(ay*yb))/(alpha*dcos(ay*ya)+beta*dcos(ay*yb)))/ay + 0.5*Ly * Heaviside(-alpha*dcos(ay*ya)-beta*dcos(ay*yb))


          mu = alpha+beta 

          
          do k = 1 , size(three%exponent)
            gamma = three%exponent(k)
            c3    = three%coefficient(k)
            o3    = three%orbital
            do l = 1 , size (four%exponent)
              delta = four%exponent(l)
              c4    = four%coefficient(l)
              o4    = four%orbital


              nu_x  = dsqrt(gamma**2+delta**2+2.d0*gamma*delta*dcos(ax*(XCD)))+eta
              nu_y  = dsqrt(gamma**2+delta**2+2.d0*gamma*delta*dcos(ay*(YCD)))+eta


              xq    = datan((gamma*dsin(ax*xC)+delta*dsin(ax*xD))/(gamma*dcos(ax*xC)+delta*dcos(ax*xD)))/ax + 0.5*Lx * Heaviside(-gamma*dcos(ax*xC)-delta*dcos(ax*xD))  
              yq    = datan((gamma*dsin(ay*yC)+delta*dsin(ay*yD))/(gamma*dcos(ay*yC)+delta*dcos(ay*yD)))/ay + 0.5*Ly * Heaviside(-gamma*dcos(ay*yC)-delta*dcos(ay*yD))
                             

              nu    = gamma + delta

              const   = (c1*c2*c3*c4) * (2.d0/ dsqrt(pi)) * (Lx*Lx*Ly*Ly)
              
              call integrate_ERI_2D(mu,nu,mu_x,nu_x,mu_y,nu_y,xp-xq,yp-yq,value_s)

              value = value + const * value_s
              
            end do
          end do 
        end do 
      end do 

end subroutine ERI_integral_4_function_toroidal_2D


subroutine integrate_ERI_2D(sigma,nu,sigma_x,nu_x,sigma_y,nu_y,xpq,ypq,result)
      
      use quadpack , only : dqagi
      use iso_c_binding
      use torus_init
      use bessel_functions
      use, intrinsic :: ieee_arithmetic

      implicit none

      ! Input parameters
      double precision, intent(in)  :: sigma_x , sigma_y, sigma 
      double precision, intent(in)  :: nu_x , nu_y ,  nu
      double precision, intent(in)  :: xpq , ypq 

      
      ! Output parameters

      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter         :: epsabs = 1.0e-8 , epsrel = 1.0e-6
      integer,parameter                  :: inf = 1 
      double precision,parameter         :: bound = 0.0d0
      integer, parameter                 :: limit = 100
      integer, parameter                 :: lenw = limit*4
      integer                            :: ier, iwork(limit), last, neval
      double precision                   :: abserr, work(lenw)
      integer,parameter                  :: Nmax = 150
      
      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      if (ier == 1) then
        write(*,'(A,I8,A)') 'Error code = ', ier
        stop 'Integration failed'
      end if

      contains

      function f_decay(t) result(ft)
        double precision, intent(in) :: t
        double precision             :: ft

        ft  = S(t)

      end function f_decay

      double precision function S(t) result(sum)

      use bessel_functions
      implicit none
      double precision, intent(in) :: t
      double precision             :: A, B, C, term
      double precision             :: sum1, sum2
      double precision             :: tol = 1.0d-12
      integer                      :: n
      
      A = 2.d0*sigma_x/(ax*ax)
      B = 2.d0*nu_x/(ax*ax)
      C = 2.d0*t*t/(ax*ax)

      sum1 = iv_scaled(0,A)*iv_scaled(0,B)*iv_scaled(0,C) * exp(A+B-2.d0*(sigma+nu)/(ax*ax))

      do n = 1, Nmax
         term = iv_scaled(n, A)*iv_scaled(n, B)*iv_scaled(n, C) * exp(A+B-2.d0*(sigma+nu)/(ax*ax))
         if (term < tol) exit
         sum1 = sum1 + term * 2.d0 * cos(dble(n)*ax*xpq) 
      end do

      A = 2.d0*sigma_y/(ay*ay)
      B = 2.d0*nu_y/(ay*ay)
      C = 2.d0*t*t/(ay*ay)

      sum2 = iv_scaled(0,A)*iv_scaled(0,B)*iv_scaled(0,C) * exp(A+B-2.d0*(sigma+nu)/(ay*ay))

      do n = 1, Nmax
         term = iv_scaled(n, A)*iv_scaled(n, B)*iv_scaled(n, C) * exp(A+B-2.d0*(sigma+nu)/(ay*ay))
         if (term < tol) exit
         sum2 = sum2 + term * 2.d0 * cos(dble(n)*ay*ypq) 
      end do

      sum = sum1 * sum2

end function S

end subroutine integrate_ERI_2D