subroutine overlap_integral_ss_toroidal_2D(r1,r2,AO1,AO2,S_ss_normal)

      use torus_init
      use classification_ERI
      use bessel_functions

      implicit none 

      !-----------------------------------------------------------------!

      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: S_ss_normal

      integer                          :: i , j 
      double precision,parameter       :: pi = 3.14159265358979323846D00
      double precision                 :: alpha , beta
      double precision                 :: c1    , c2 
      double precision                 :: x1 , x2 , y1 , y2 
      double precision                 :: X , Y
      double precision                 :: const 
      double precision                 :: overlap_x , overlap_y
      double precision                 :: gamma_x , gamma_y     
      double precision                 :: I_0_gamma_x , I_0_gamma_y 


      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)

      X        = (x1 - x2)
      Y        = (y1 - y2)

      !-----------------------------------------------------------------!
 
      S_ss_normal = 0.d0
      
      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        c1    = AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta = AO2%exponent(j)
          c2   = AO2%coefficient(j)

              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))
              gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(Y)))

              I_0_gamma_x = iv_scaled(0,2.d0*gamma_x/(ax*ax))
              I_0_gamma_y = iv_scaled(0,2.d0*gamma_y/(ay*ay))

              overlap_x   =  Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
              overlap_y   =  Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
              
              S_ss_normal =  S_ss_normal + c1 * c2 * overlap_x * overlap_y

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine overlap_integral_ss_toroidal_2D
