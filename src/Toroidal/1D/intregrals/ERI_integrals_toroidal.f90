subroutine ERI_integral_4_function_toroidal(one,two,three,four,value)
      
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
      double precision                :: xpA , xpB , xqC , xqD , phi 
      double precision                :: test 
      integer                         :: pattern_id, encode_orbital_pattern

      ! --------------------------------------------------------------- !
      double precision                :: albe     , gade
      double precision                :: inv_albe , inv_gade
      double precision                :: reduced_mu , reduced_nu



      !-----------------------------------------------------------------!

      value   = 0.d0 
      value_s = 0.d0 
            
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
          mu_y  = alpha + beta
          mu_z  = alpha + beta

          call bary_center_toroidal_x(alpha,beta,xa,xb,xp)

          !   Real Gaussian   !
      
          albe                = alpha + beta 
          inv_albe            = 1.d0 / albe 
          reduced_mu          = alpha * beta * inv_albe

          yp    = ( alpha * ya + beta * yb ) * inv_albe
          zp    = ( alpha * za + beta * zb ) * inv_albe
                    
          do k = 1 , size(three%exponent)
            gamma = three%exponent(k)
            c3    = three%coefficient(k)
            o3    = three%orbital
            do l = 1 , size (four%exponent)
              delta = four%exponent(l)
              c4    = four%coefficient(l)
              o4    = four%orbital

              call bary_exponent_x(gamma,delta,XCD,nu_x)
              nu_y  =  gamma+delta
              nu_z  =  gamma+delta

              call bary_center_toroidal_x(gamma,delta,xc,xd,xq)

              !   Real Gaussian   !
      
              gade                = gamma + delta 
              inv_gade            = 1.d0 / gade 
              reduced_nu          = gamma * delta * inv_gade

              yq     = (gamma * yc + delta * yd) * inv_gade
              zq     = (gamma * zc + delta * zd) * inv_gade

              const  = (c1*c2*c3*c4) * 2.d0 /dsqrt(pi) * (Lx * Lx) * dexp(- reduced_mu * (Yab*Yab + Zab*Zab) ) * dexp(- reduced_nu * (Ycd*Ycd + Zcd*Zcd) )

              test = dexp(-0.5d0*(alpha+beta+gamma+delta-mu_x-nu_x)*(Lx*Lx)/(pi*pi)) * const

              xpA     = ax*(xp - xa)
              xpB     = ax*(xp - xb)
              xqC     = ax*(xq - xc)
              xqD     = ax*(xq - xd)
              phi     = ax*(xp - xq)

              pattern_id = encode_orbital_pattern(o1, o2, o3, o4)

              if ( dabs(test) <= 2.22d-16 ) cycle

                call integrate_ERI_sum(pattern_id,albe,gade,mu_x,nu_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq&
                &                                                                        ,ya,yb,yc,yd,yp,yq&
                &                                                                        ,za,zb,zc,zd,zp,zq,value_s)
                value  = value    + test * value_s

              
            end do
          end do 
        end do 
      end do 

end subroutine ERI_integral_4_function_toroidal


subroutine ERI_integral_4_function_toroidal_new(one,two,three,four,value)
      
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
      double precision                :: xpA , xpB , xqC , xqD , phi 
      double precision                :: test 
      integer                         :: pattern_id, encode_orbital_pattern

      ! --------------------------------------------------------------- !
      double precision                :: albe     , gade
      double precision                :: inv_albe , inv_gade
      double precision                :: reduced_mu , reduced_nu
      !-----------------------------------------------------------------!

      value   = 0.d0 
      value_s = 0.d0 
            
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
          mu_y  = alpha + beta
          mu_z  = alpha + beta

          call bary_center_toroidal_x(alpha,beta,xa,xb,xp)

          !   Real Gaussian   !
      
          albe                = alpha + beta 
          inv_albe            = 1.d0 / albe 
          reduced_mu          = alpha * beta * inv_albe

          yp    = ( alpha * ya + beta * yb ) * inv_albe
          zp    = ( alpha * za + beta * zb ) * inv_albe
                    
          do k = 1 , size(three%exponent)
            gamma = three%exponent(k)
            c3    = three%coefficient(k)
            o3    = three%orbital
            do l = 1 , size (four%exponent)
              delta = four%exponent(l)
              c4    = four%coefficient(l)
              o4    = four%orbital

              call bary_exponent_x(gamma,delta,XCD,nu_x)
              nu_y  =  gamma+delta
              nu_z  =  gamma+delta

              call bary_center_toroidal_x(gamma,delta,xc,xd,xq)

              !   Real Gaussian   !
      
              gade                = gamma + delta 
              inv_gade            = 1.d0 / gade 
              reduced_nu          = gamma * delta * inv_gade

              yq     = (gamma * yc + delta * yd) * inv_gade
              zq     = (gamma * zc + delta * zd) * inv_gade

              const  = (c1*c2*c3*c4) * 2.d0 /dsqrt(pi) * (Lx * Lx) * dexp(- reduced_mu * (Yab*Yab) ) * dexp(- reduced_nu * (Ycd*Ycd) ) * dexp(- reduced_mu * (Zab*Zab) ) * dexp(- reduced_nu * (Zcd*Zcd) )

              test = dexp(-0.5d0*(alpha+beta+gamma+delta-mu_x-nu_x)*(Lx*Lx)/(pi*pi)) * const

              xpA     = ax*(xp - xa)
              xpB     = ax*(xp - xb)
              xqC     = ax*(xq - xc)
              xqD     = ax*(xq - xd)
              phi     = ax*(xp - xq)

              pattern_id = encode_orbital_pattern(o1, o2, o3, o4)

              if ( dabs(test) <= 2.22d-16 ) cycle

                call integrate_ERI_sum_new(pattern_id,albe,gade,mu_x,nu_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq&
                &                                                                        ,ya,yb,yc,yd,yp,yq&
                &                                                                        ,za,zb,zc,zd,zp,zq,value_s)
                value  = value    + test * value_s

              
            end do
          end do 
        end do 
      end do 

end subroutine ERI_integral_4_function_toroidal_new