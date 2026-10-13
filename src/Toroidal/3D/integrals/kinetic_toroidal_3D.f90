subroutine kinetic_matrix_toroidal_3D(number_of_atoms,number_of_functions,atoms,AO,kinetic)

      use files 
      use torus_init
      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      integer                      :: i , j , k , l 
      integer                      :: number_of_atoms
      integer                      :: number_of_functions
      integer                      :: fpuc

      type(atom)                   :: atoms(number_of_atoms)

      type(ERI_function)           :: AO (number_of_functions)
      type(ERI_function)           :: AO1 , AO2


      double precision             :: r1(3) , r2(3)

      ! output ! 

      double precision             :: kinetic_tmp(number_of_functions,number_of_functions)
      double precision,intent(out) :: kinetic(number_of_functions,number_of_functions)

      !-----------------------------------------------------------------!
      
      ! functions_per_unitcell ! 

      fpuc = 0 

      do i = 1 , number_of_atom_in_unitcell 
        fpuc = fpuc + atoms(i)%num_s_function +  3 * atoms(i)%num_p_function
      end do 

      !-----------------------------------------------------------------!


      kinetic(:,:) = 0.d0 

      do i = 1 , fpuc
      !do i = 1 , number_of_functions
        do j = 1 , number_of_functions

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z

          if (AO1%orbital =="s" .and. AO2%orbital == "s") then
        
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call kinetic_integral_ss_toroidal_3D(r1,r2,AO1,AO2,kinetic_tmp(i,j))
              end do 
            end do 

          end if        

          if (AO1%orbital =="s" .and. AO2%orbital(:1) == "p") then
            
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call kinetic_integral_sp_toroidal_3D(r1,r2,AO1,AO2,kinetic_tmp(i,j))
              end do 
            end do

          end if

          if (AO1%orbital(:1) == "p" .and. AO2%orbital == "s") then
                
            do k = 1, size(AO1%exponent)
              do l = 1, size(AO2%exponent)
                call kinetic_integral_sp_toroidal_3D(r2, r1, AO2, AO1, kinetic_tmp(i,j))
              end do 
            end do

          end if

          if (AO1%orbital(:1) == "p" .and. AO2%orbital(:1) == "p") then
                
            do k = 1, size(AO1%exponent)
              do l = 1, size(AO2%exponent)
                call kinetic_integral_pp_toroidal_3D(r1, r2, AO1, AO2, kinetic_tmp(i,j))
              end do 
            end do

          end if

        
        end do 
      end do 

!      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
!      !                    symmetry of the integrals                    !
!      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!


      call symmetry_of_integrals(number_of_functions,fpuc,kinetic_tmp,kinetic)

      ! do i = 1 , number_of_functions - 1 
      !   do j = 1 , number_of_functions
      !     if (abs(kinetic(i,j)) < 1e-15) kinetic(i,j) = 0.d0 
      !     kinetic(j,i) = kinetic(i,j)
      !   end do 
      ! end do


end subroutine kinetic_matrix_toroidal_3D


subroutine kinetic_integral_toroidal_3D(pattern_id,r1,r2,AO1,AO2,K_t)

      use torus_init
      use classification_ERI
      use bessel_functions
      use functions
      
      implicit none 

      !-----------------------------------------------------------------!

      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: K_t

      integer                          :: i , j 
      double precision,parameter       :: pi = 3.14159265358979323846D00
      double precision                 :: alpha  , beta   , beta2 
      double precision                 :: c1     , c2     , const
      double precision                 :: x1     , x2     , X  
      double precision                 :: y1     , y2     , Y
      double precision                 :: z1     , z2     , Z 
      double precision                 :: kx     , ky     , kz 
      double precision                 :: sx_int , sy_int , sz_int  
      double precision                 :: K 
       

      ! Clifford ! 

      double precision                 :: ax2 , inv_ax , inv_ax2 
      double precision                 :: ay2 , inv_ay , inv_ay2 
      double precision                 :: az2 , inv_az , inv_az2 
      double precision                 :: px  , py  , pz   
      double precision                 :: A_x , B_x
      double precision                 :: A_y , B_y
      double precision                 :: A_z , B_z
      double precision                 :: xp  , yp  , zp
      double precision                 :: spa_x  , spb_x
      double precision                 :: cpa_x  , cpb_x
      double precision                 :: spa_y  , spb_y
      double precision                 :: cpa_y  , cpb_y
      double precision                 :: spa_z  , spb_z
      double precision                 :: cpa_z  , cpb_z
      
      integer                          :: pattern_id

      double precision                 :: D_00_x , D_00_y , D_00_z
      double precision                 :: D_01_x , D_01_y , D_01_z
      double precision                 :: D_10_x , D_10_y , D_10_z
      double precision                 :: D_11_x , D_11_y , D_11_z
      double precision                 :: S_00_x , S_00_y , S_00_z
      double precision                 :: S_01_x , S_01_y , S_01_z
      double precision                 :: S_10_x , S_10_y , S_10_z
      double precision                 :: S_11_x , S_11_y , S_11_z



      !-----------------------------------------------------------------!
      !-----------------------------------------------------------------!

      x1  = r1(1) ; x2  = r2(1) 
      y1  = r1(2) ; y2  = r2(2)
      z1  = r1(3) ; z2  = r2(3)

      X   = (x1 - x2)
      Y   = (y1 - y2)
      Z   = (z1 - z2)

      ax2 = ax * ax
      ay2 = ay * ay
      az2 = az * az

      inv_ax  = 1.d0 / ax 
      inv_ay  = 1.d0 / ay 
      inv_az  = 1.d0 / az 

      inv_ax2 = inv_ax * inv_ax 
      inv_ay2 = inv_ay * inv_ay 
      inv_az2 = inv_az * inv_az 

      !-----------------------------------------------------------------!
 
      K   = 0.d0
      K_t = 0.d0 

      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        c1    = AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta = AO2%exponent(j)
          c2   = AO2%coefficient(j)

          const   =  c1 * c2

          ! ----------------- !
          ! Clifford Gaussian !
          ! ----------------- ! 

          call bary_exponent_x         (alpha,beta,X,px)
          call bary_center_toroidal_x  (alpha,beta,x1,x2,xp)

          call bary_exponent_y         (alpha,beta,Y,py)
          call bary_center_toroidal_y  (alpha,beta,y1,y2,yp)

          call bary_exponent_z         (alpha,beta,Z,pz)
          call bary_center_toroidal_z  (alpha,beta,z1,z2,zp)

          beta2  = beta * beta

          A_x      = 2.d0*px/ax2
          B_x      = 2.d0*beta/ax2
          A_y      = 2.d0*py/ay2
          B_y      = 2.d0*beta/ay2
          A_z      = 2.d0*pz/az2
          B_z      = 2.d0*beta/az2

          kx     = dexp(-2.d0 * ( alpha + beta - px ) / ax2)
          ky     = dexp(-2.d0 * ( alpha + beta - py ) / ay2)
          kz     = dexp(-2.d0 * ( alpha + beta - pz ) / az2)

          sx_int = Lx * iv_scaled(0, A_x)
          sy_int = Ly * iv_scaled(0, A_y)
          sz_int = Lz * iv_scaled(0, A_z)

          spa_x    = dsin(ax*(xp-x1)) ; cpa_x    = dcos(ax*(xp-x1))
          spb_x    = dsin(ax*(xp-x2)) ; cpb_x    = dcos(ax*(xp-x2))
          spa_y    = dsin(ay*(yp-y1)) ; cpa_y    = dcos(ay*(yp-y1))
          spb_y    = dsin(ay*(yp-y2)) ; cpb_y    = dcos(ay*(yp-y2))
          spa_z    = dsin(az*(zp-z1)) ; cpa_z    = dcos(az*(zp-z1))
          spb_z    = dsin(az*(zp-z2)) ; cpb_z    = dcos(az*(zp-z2)) 

      select case(pattern_id)

      case (00) ! | s     s     ( 1 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=0, n=0)
      
        S_00_x = sx_int
        S_00_y = sy_int
        S_00_z = sz_int
        D_00_x = Icliff(2,0,A_x)*B_x**2*ax2*cpb_x**2 + Icliff(1,1,A_x)*2*B_x**2*ax2*cpb_x*spb_x + Icliff(0,2,A_x)*B_x**2*ax2*spb_x**2 - Icliff(0,1,A_x)*B_x*ax2*cpb_x + Icliff(1,0,A_x)*B_x*ax2*spb_x
        D_00_y = Icliff(2,0,A_y)*B_y**2*ay2*cpb_y**2 + Icliff(1,1,A_y)*2*B_y**2*ay2*cpb_y*spb_y + Icliff(0,2,A_y)*B_y**2*ay2*spb_y**2 - Icliff(0,1,A_y)*B_y*ay2*cpb_y + Icliff(1,0,A_y)*B_y*ay2*spb_y
        D_00_z = Icliff(2,0,A_z)*B_z**2*az2*cpb_z**2 + Icliff(1,1,A_z)*2*B_z**2*az2*cpb_z*spb_z + Icliff(0,2,A_z)*B_z**2*az2*spb_z**2 - Icliff(0,1,A_z)*B_z*az2*cpb_z + Icliff(1,0,A_z)*B_z*az2*spb_z
      
        K = -0.5d0 * ( D_00_x * S_00_y * S_00_z + S_00_x * D_00_y * S_00_z + S_00_x * S_00_y * D_00_z )
      
      case (01) ! | s     px    ( 2 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=1, m=0, n=0)
      
        S_01_x = inv_ax * (Icliff(1,0,A_x)*cpb_x + Icliff(0,1,A_x)*spb_x)
        S_00_y = sy_int
        S_00_z = sz_int
        D_01_x = inv_ax * (Icliff(3,0,A_x)*B_x**2*ax2*cpb_x**3 + Icliff(2,1,A_x)*3*B_x**2*ax2*cpb_x**2*spb_x + Icliff(1,2,A_x)*3*B_x**2*ax2*cpb_x*spb_x**2 + Icliff(0,3,A_x)*B_x**2*ax2*spb_x**3 - Icliff(1,1,A_x)*3*B_x*ax2*cpb_x**2 - Icliff(0,2,A_x)*3*B_x*ax2*cpb_x*spb_x + Icliff(2,0,A_x)*3*B_x*ax2*cpb_x*spb_x + Icliff(1,1,A_x)*3*B_x*ax2*spb_x**2 - Icliff(1,0,A_x)*ax2*cpb_x - Icliff(0,1,A_x)*ax2*spb_x)
        D_00_y = Icliff(2,0,A_y)*B_y**2*ay2*cpb_y**2 + Icliff(1,1,A_y)*2*B_y**2*ay2*cpb_y*spb_y + Icliff(0,2,A_y)*B_y**2*ay2*spb_y**2 - Icliff(0,1,A_y)*B_y*ay2*cpb_y + Icliff(1,0,A_y)*B_y*ay2*spb_y
        D_00_z = Icliff(2,0,A_z)*B_z**2*az2*cpb_z**2 + Icliff(1,1,A_z)*2*B_z**2*az2*cpb_z*spb_z + Icliff(0,2,A_z)*B_z**2*az2*spb_z**2 - Icliff(0,1,A_z)*B_z*az2*cpb_z + Icliff(1,0,A_z)*B_z*az2*spb_z
    
        K = -0.5d0 * ( D_01_x * S_00_y * S_00_z + S_01_x * D_00_y * S_00_z + S_01_x * S_00_y * D_00_z )
      
      case (02) ! | s     py    ( 3 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        S_00_x = sx_int
        S_01_y = inv_ay * (Icliff(1,0,A_y)*cpb_y + Icliff(0,1,A_y)*spb_y)
        S_00_z = sz_int
        D_00_x = Icliff(2,0,A_x)*B_x**2*ax2*cpb_x**2 + Icliff(1,1,A_x)*2*B_x**2*ax2*cpb_x*spb_x + Icliff(0,2,A_x)*B_x**2*ax2*spb_x**2 - Icliff(0,1,A_x)*B_x*ax2*cpb_x + Icliff(1,0,A_x)*B_x*ax2*spb_x
        D_01_y = inv_ay * (Icliff(3,0,A_y)*B_y**2*ay2*cpb_y**3 + Icliff(2,1,A_y)*3*B_y**2*ay2*cpb_y**2*spb_y + Icliff(1,2,A_y)*3*B_y**2*ay2*cpb_y*spb_y**2 + Icliff(0,3,A_y)*B_y**2*ay2*spb_y**3 - Icliff(1,1,A_y)*3*B_y*ay2*cpb_y**2 - Icliff(0,2,A_y)*3*B_y*ay2*cpb_y*spb_y + Icliff(2,0,A_y)*3*B_y*ay2*cpb_y*spb_y + Icliff(1,1,A_y)*3*B_y*ay2*spb_y**2 - Icliff(1,0,A_y)*ay2*cpb_y - Icliff(0,1,A_y)*ay2*spb_y)
        D_00_z = Icliff(2,0,A_z)*B_z**2*az2*cpb_z**2 + Icliff(1,1,A_z)*2*B_z**2*az2*cpb_z*spb_z + Icliff(0,2,A_z)*B_z**2*az2*spb_z**2 - Icliff(0,1,A_z)*B_z*az2*cpb_z + Icliff(1,0,A_z)*B_z*az2*spb_z
      
        K = -0.5d0 * ( D_00_x * S_01_y * S_00_z + S_00_x * D_01_y * S_00_z + S_00_x * S_01_y * D_00_z )
      
      case (03) ! | s     pz    ( 4 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        S_00_x = sx_int
        S_00_y = sy_int
        S_01_z = inv_az * (Icliff(1,0,A_z)*cpb_z + Icliff(0,1,A_z)*spb_z)
        D_00_x = Icliff(2,0,A_x)*B_x**2*ax2*cpb_x**2 + Icliff(1,1,A_x)*2*B_x**2*ax2*cpb_x*spb_x + Icliff(0,2,A_x)*B_x**2*ax2*spb_x**2 - Icliff(0,1,A_x)*B_x*ax2*cpb_x + Icliff(1,0,A_x)*B_x*ax2*spb_x
        D_00_y = Icliff(2,0,A_y)*B_y**2*ay2*cpb_y**2 + Icliff(1,1,A_y)*2*B_y**2*ay2*cpb_y*spb_y + Icliff(0,2,A_y)*B_y**2*ay2*spb_y**2 - Icliff(0,1,A_y)*B_y*ay2*cpb_y + Icliff(1,0,A_y)*B_y*ay2*spb_y
        D_01_z = inv_az * (Icliff(3,0,A_z)*B_z**2*az2*cpb_z**3 + Icliff(2,1,A_z)*3*B_z**2*az2*cpb_z**2*spb_z + Icliff(1,2,A_z)*3*B_z**2*az2*cpb_z*spb_z**2 + Icliff(0,3,A_z)*B_z**2*az2*spb_z**3 - Icliff(1,1,A_z)*3*B_z*az2*cpb_z**2 - Icliff(0,2,A_z)*3*B_z*az2*cpb_z*spb_z + Icliff(2,0,A_z)*3*B_z*az2*cpb_z*spb_z + Icliff(1,1,A_z)*3*B_z*az2*spb_z**2 - Icliff(1,0,A_z)*az2*cpb_z - Icliff(0,1,A_z)*az2*spb_z)
      
        K = -0.5d0 * ( D_00_x * S_00_y * S_01_z + S_00_x * D_00_y * S_01_z + S_00_x * S_00_y * D_01_z )
      
      case (10) ! | px    s     ( 5 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=0)
      
        S_10_x = inv_ax * (Icliff(1,0,A_x)*cpa_x + Icliff(0,1,A_x)*spa_x)
        S_00_y = sy_int
        S_00_z = sz_int
        D_10_x = inv_ax * (Icliff(3,0,A_x)*B_x**2*ax2*cpa_x*cpb_x**2 + Icliff(2,1,A_x)*2*B_x**2*ax2*cpa_x*cpb_x*spb_x + Icliff(1,2,A_x)*B_x**2*ax2*cpa_x*spb_x**2 + Icliff(2,1,A_x)*B_x**2*ax2*cpb_x**2*spa_x + Icliff(1,2,A_x)*2*B_x**2*ax2*cpb_x*spa_x*spb_x + Icliff(0,3,A_x)*B_x**2*ax2*spa_x*spb_x**2 - Icliff(1,1,A_x)*B_x*ax2*cpa_x*cpb_x + Icliff(2,0,A_x)*B_x*ax2*cpa_x*spb_x - Icliff(0,2,A_x)*B_x*ax2*cpb_x*spa_x + Icliff(1,1,A_x)*B_x*ax2*spa_x*spb_x)
        D_00_y = Icliff(2,0,A_y)*B_y**2*ay2*cpb_y**2 + Icliff(1,1,A_y)*2*B_y**2*ay2*cpb_y*spb_y + Icliff(0,2,A_y)*B_y**2*ay2*spb_y**2 - Icliff(0,1,A_y)*B_y*ay2*cpb_y + Icliff(1,0,A_y)*B_y*ay2*spb_y
        D_00_z = Icliff(2,0,A_z)*B_z**2*az2*cpb_z**2 + Icliff(1,1,A_z)*2*B_z**2*az2*cpb_z*spb_z + Icliff(0,2,A_z)*B_z**2*az2*spb_z**2 - Icliff(0,1,A_z)*B_z*az2*cpb_z + Icliff(1,0,A_z)*B_z*az2*spb_z
      
        K = -0.5d0 * ( D_10_x * S_00_y * S_00_z + S_10_x * D_00_y * S_00_z + S_10_x * S_00_y * D_00_z )
      
      case (11) ! | px    px    ( 6 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=1, m=0, n=0)
      
        S_11_x = inv_ax2 * (Icliff(2,0,A_x)*cpa_x*cpb_x + Icliff(1,1,A_x)*cpa_x*spb_x + Icliff(1,1,A_x)*cpb_x*spa_x + Icliff(0,2,A_x)*spa_x*spb_x)
        S_00_y = sy_int
        S_00_z = sz_int
        D_11_x = inv_ax2 * (Icliff(4,0,A_x)*B_x**2*ax2*cpa_x*cpb_x**3 + Icliff(3,1,A_x)*3*B_x**2*ax2*cpa_x*cpb_x**2*spb_x + Icliff(2,2,A_x)*3*B_x**2*ax2*cpa_x*cpb_x*spb_x**2 + Icliff(1,3,A_x)*B_x**2*ax2*cpa_x*spb_x**3 + Icliff(3,1,A_x)*B_x**2*ax2*cpb_x**3*spa_x + Icliff(2,2,A_x)*3*B_x**2*ax2*cpb_x**2*spa_x*spb_x + Icliff(1,3,A_x)*3*B_x**2*ax2*cpb_x*spa_x*spb_x**2 + Icliff(0,4,A_x)*B_x**2*ax2*spa_x*spb_x**3 - Icliff(2,1,A_x)*3*B_x*ax2*cpa_x*cpb_x**2 - Icliff(1,2,A_x)*3*B_x*ax2*cpa_x*cpb_x*spb_x + Icliff(3,0,A_x)*3*B_x*ax2*cpa_x*cpb_x*spb_x + Icliff(2,1,A_x)*3*B_x*ax2*cpa_x*spb_x**2 - Icliff(1,2,A_x)*3*B_x*ax2*cpb_x**2*spa_x - Icliff(0,3,A_x)*3*B_x*ax2*cpb_x*spa_x*spb_x + Icliff(2,1,A_x)*3*B_x*ax2*cpb_x*spa_x*spb_x + Icliff(1,2,A_x)*3*B_x*ax2*spa_x*spb_x**2 - Icliff(2,0,A_x)*ax2*cpa_x*cpb_x - Icliff(1,1,A_x)*ax2*cpa_x*spb_x - Icliff(1,1,A_x)*ax2*cpb_x*spa_x - Icliff(0,2,A_x)*ax2*spa_x*spb_x)
        D_00_y = Icliff(2,0,A_y)*B_y**2*ay2*cpb_y**2 + Icliff(1,1,A_y)*2*B_y**2*ay2*cpb_y*spb_y + Icliff(0,2,A_y)*B_y**2*ay2*spb_y**2 - Icliff(0,1,A_y)*B_y*ay2*cpb_y + Icliff(1,0,A_y)*B_y*ay2*spb_y
        D_00_z = Icliff(2,0,A_z)*B_z**2*az2*cpb_z**2 + Icliff(1,1,A_z)*2*B_z**2*az2*cpb_z*spb_z + Icliff(0,2,A_z)*B_z**2*az2*spb_z**2 - Icliff(0,1,A_z)*B_z*az2*cpb_z + Icliff(1,0,A_z)*B_z*az2*spb_z
      
        K = -0.5d0 * ( D_11_x * S_00_y * S_00_z + S_11_x * D_00_y * S_00_z + S_11_x * S_00_y * D_00_z )
      
      case (12) ! | px    py    ( 7 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        S_10_x = inv_ax * (Icliff(1,0,A_x)*cpa_x + Icliff(0,1,A_x)*spa_x)
        S_01_y = inv_ay * (Icliff(1,0,A_y)*cpb_y + Icliff(0,1,A_y)*spb_y)
        S_00_z = sz_int
        D_10_x = inv_ax * (Icliff(3,0,A_x)*B_x**2*ax2*cpa_x*cpb_x**2 + Icliff(2,1,A_x)*2*B_x**2*ax2*cpa_x*cpb_x*spb_x + Icliff(1,2,A_x)*B_x**2*ax2*cpa_x*spb_x**2 + Icliff(2,1,A_x)*B_x**2*ax2*cpb_x**2*spa_x + Icliff(1,2,A_x)*2*B_x**2*ax2*cpb_x*spa_x*spb_x + Icliff(0,3,A_x)*B_x**2*ax2*spa_x*spb_x**2 - Icliff(1,1,A_x)*B_x*ax2*cpa_x*cpb_x + Icliff(2,0,A_x)*B_x*ax2*cpa_x*spb_x - Icliff(0,2,A_x)*B_x*ax2*cpb_x*spa_x + Icliff(1,1,A_x)*B_x*ax2*spa_x*spb_x)
        D_01_y = inv_ay * (Icliff(3,0,A_y)*B_y**2*ay2*cpb_y**3 + Icliff(2,1,A_y)*3*B_y**2*ay2*cpb_y**2*spb_y + Icliff(1,2,A_y)*3*B_y**2*ay2*cpb_y*spb_y**2 + Icliff(0,3,A_y)*B_y**2*ay2*spb_y**3 - Icliff(1,1,A_y)*3*B_y*ay2*cpb_y**2 - Icliff(0,2,A_y)*3*B_y*ay2*cpb_y*spb_y + Icliff(2,0,A_y)*3*B_y*ay2*cpb_y*spb_y + Icliff(1,1,A_y)*3*B_y*ay2*spb_y**2 - Icliff(1,0,A_y)*ay2*cpb_y - Icliff(0,1,A_y)*ay2*spb_y)
        D_00_z = Icliff(2,0,A_z)*B_z**2*az2*cpb_z**2 + Icliff(1,1,A_z)*2*B_z**2*az2*cpb_z*spb_z + Icliff(0,2,A_z)*B_z**2*az2*spb_z**2 - Icliff(0,1,A_z)*B_z*az2*cpb_z + Icliff(1,0,A_z)*B_z*az2*spb_z
      
        K = -0.5d0 * ( D_10_x * S_01_y * S_00_z + S_10_x * D_01_y * S_00_z + S_10_x * S_01_y * D_00_z )
      
      case (13) ! | px    pz    ( 8 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        S_10_x = inv_ax * (Icliff(1,0,A_x)*cpa_x + Icliff(0,1,A_x)*spa_x)
        S_00_y = sy_int
        S_01_z = inv_az * (Icliff(1,0,A_z)*cpa_z + Icliff(0,1,A_z)*spa_z)
        D_10_x = inv_ax * (Icliff(3,0,A_x)*B_x**2*ax2*cpa_x*cpb_x**2 + Icliff(2,1,A_x)*2*B_x**2*ax2*cpa_x*cpb_x*spb_x + Icliff(1,2,A_x)*B_x**2*ax2*cpa_x*spb_x**2 + Icliff(2,1,A_x)*B_x**2*ax2*cpb_x**2*spa_x + Icliff(1,2,A_x)*2*B_x**2*ax2*cpb_x*spa_x*spb_x + Icliff(0,3,A_x)*B_x**2*ax2*spa_x*spb_x**2 - Icliff(1,1,A_x)*B_x*ax2*cpa_x*cpb_x + Icliff(2,0,A_x)*B_x*ax2*cpa_x*spb_x - Icliff(0,2,A_x)*B_x*ax2*cpb_x*spa_x + Icliff(1,1,A_x)*B_x*ax2*spa_x*spb_x)
        D_00_y = Icliff(2,0,A_y)*B_y**2*ay2*cpb_y**2 + Icliff(1,1,A_y)*2*B_y**2*ay2*cpb_y*spb_y + Icliff(0,2,A_y)*B_y**2*ay2*spb_y**2 - Icliff(0,1,A_y)*B_y*ay2*cpb_y + Icliff(1,0,A_y)*B_y*ay2*spb_y
        D_01_z = inv_az * (Icliff(3,0,A_z)*B_z**2*az2*cpb_z**3 + Icliff(2,1,A_z)*3*B_z**2*az2*cpb_z**2*spb_z + Icliff(1,2,A_z)*3*B_z**2*az2*cpb_z*spb_z**2 + Icliff(0,3,A_z)*B_z**2*az2*spb_z**3 - Icliff(1,1,A_z)*3*B_z*az2*cpb_z**2 - Icliff(0,2,A_z)*3*B_z*az2*cpb_z*spb_z + Icliff(2,0,A_z)*3*B_z*az2*cpb_z*spb_z + Icliff(1,1,A_z)*3*B_z*az2*spb_z**2 - Icliff(1,0,A_z)*az2*cpb_z - Icliff(0,1,A_z)*az2*spb_z)
      
        K = -0.5d0 * ( D_10_x * S_00_y * S_01_z + S_10_x * D_00_y * S_01_z + S_10_x * S_00_y * D_01_z )
      
      case (20) ! | py    s     ( 9 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=0)
      
        S_00_x = sx_int
        S_10_y = inv_ay * (Icliff(1,0,A_y)*cpa_y + Icliff(0,1,A_y)*spa_y)
        S_00_z = sz_int
        D_00_x = Icliff(2,0,A_x)*B_x**2*ax2*cpb_x**2 + Icliff(1,1,A_x)*2*B_x**2*ax2*cpb_x*spb_x + Icliff(0,2,A_x)*B_x**2*ax2*spb_x**2 - Icliff(0,1,A_x)*B_x*ax2*cpb_x + Icliff(1,0,A_x)*B_x*ax2*spb_x
        D_10_y = inv_ay * (Icliff(3,0,A_y)*B_y**2*ay2*cpa_y*cpb_y**2 + Icliff(2,1,A_y)*2*B_y**2*ay2*cpa_y*cpb_y*spb_y + Icliff(1,2,A_y)*B_y**2*ay2*cpa_y*spb_y**2 + Icliff(2,1,A_y)*B_y**2*ay2*cpb_y**2*spa_y + Icliff(1,2,A_y)*2*B_y**2*ay2*cpb_y*spa_y*spb_y + Icliff(0,3,A_y)*B_y**2*ay2*spa_y*spb_x**2 - Icliff(1,1,A_y)*B_y*ay2*cpa_y*cpb_y + Icliff(2,0,A_y)*B_y*ay2*cpa_y*spb_y - Icliff(0,2,A_y)*B_y*ay2*cpb_y*spa_y + Icliff(1,1,A_y)*B_y*ay2*spa_y*spb_y)
        D_00_z = Icliff(2,0,A_z)*B_z**2*az2*cpb_z**2 + Icliff(1,1,A_z)*2*B_z**2*az2*cpb_z*spb_z + Icliff(0,2,A_z)*B_z**2*az2*spb_z**2 - Icliff(0,1,A_z)*B_z*az2*cpb_z + Icliff(1,0,A_z)*B_z*az2*spb_z
      
        K = -0.5d0 * ( D_00_x * S_10_y * S_00_z + S_00_x * D_10_y * S_00_z + S_00_x * S_10_y * D_00_z )
      
      case (21) ! | py    px    ( 10 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=1, m=0, n=0)
      
        S_01_x = inv_ax * (Icliff(1,0,A_x)*cpb_x + Icliff(0,1,A_x)*spb_x)
        S_10_y = inv_ay * (Icliff(1,0,A_y)*cpa_y + Icliff(0,1,A_y)*spa_y)
        S_00_z = sz_int
        D_01_x = inv_ax * (Icliff(3,0,A_x)*B_x**2*ax2*cpb_x**3 + Icliff(2,1,A_x)*3*B_x**2*ax2*cpb_x**2*spb_x + Icliff(1,2,A_x)*3*B_x**2*ax2*cpb_x*spb_x**2 + Icliff(0,3,A_x)*B_x**2*ax2*spb_x**3 - Icliff(1,1,A_x)*3*B_x*ax2*cpb_x**2 - Icliff(0,2,A_x)*3*B_x*ax2*cpb_x*spb_x + Icliff(2,0,A_x)*3*B_x*ax2*cpb_x*spb_x + Icliff(1,1,A_x)*3*B_x*ax2*spb_x**2 - Icliff(1,0,A_x)*ax2*cpb_x - Icliff(0,1,A_x)*ax2*spb_x)
        D_10_y = inv_ay * (Icliff(3,0,A_y)*B_y**2*ay2*cpa_y*cpb_y**2 + Icliff(2,1,A_y)*2*B_y**2*ay2*cpa_y*cpb_y*spb_y + Icliff(1,2,A_y)*B_y**2*ay2*cpa_y*spb_y**2 + Icliff(2,1,A_y)*B_y**2*ay2*cpb_y**2*spa_y + Icliff(1,2,A_y)*2*B_y**2*ay2*cpb_y*spa_y*spb_y + Icliff(0,3,A_y)*B_y**2*ay2*spa_y*spb_x**2 - Icliff(1,1,A_y)*B_y*ay2*cpa_y*cpb_y + Icliff(2,0,A_y)*B_y*ay2*cpa_y*spb_y - Icliff(0,2,A_y)*B_y*ay2*cpb_y*spa_y + Icliff(1,1,A_y)*B_y*ay2*spa_y*spb_y)
        D_00_z = Icliff(2,0,A_z)*B_z**2*az2*cpb_z**2 + Icliff(1,1,A_z)*2*B_z**2*az2*cpb_z*spb_z + Icliff(0,2,A_z)*B_z**2*az2*spb_z**2 - Icliff(0,1,A_z)*B_z*az2*cpb_z + Icliff(1,0,A_z)*B_z*az2*spb_z
      
        K = -0.5d0 * ( D_01_x * S_10_y * S_00_z + S_01_x * D_10_y * S_00_z + S_01_x * S_10_y * D_00_z )
      
      case (22) ! | py    py    ( 11 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=1, n=0)
      
        S_00_x = sx_int
        S_11_y = inv_ay2 * (Icliff(2,0,A_y)*cpa_y*cpb_y + Icliff(1,1,A_y)*cpa_y*spb_y + Icliff(1,1,A_y)*cpb_y*spa_y + Icliff(0,2,A_y)*spa_y*spb_y)
        S_00_z = sz_int
        D_00_x = Icliff(2,0,A_x)*B_x**2*ax2*cpb_x**2 + Icliff(1,1,A_x)*2*B_x**2*ax2*cpb_x*spb_x + Icliff(0,2,A_x)*B_x**2*ax2*spb_x**2 - Icliff(0,1,A_x)*B_x*ax2*cpb_x + Icliff(1,0,A_x)*B_x*ax2*spb_x
        D_11_y = inv_ay2 * (Icliff(4,0,A_y)*B_y**2*ay2*cpa_y*cpb_y**3 + Icliff(3,1,A_y)*3*B_y**2*ay2*cpa_y*cpb_y**2*spb_y + Icliff(2,2,A_y)*3*B_y**2*ay2*cpa_y*cpb_y*spb_y**2 + Icliff(1,3,A_y)*B_y**2*ay2*cpa_y*spb_y**3 + Icliff(3,1,A_y)*B_y**2*ay2*cpb_y**3*spa_y + Icliff(2,2,A_y)*3*B_y**2*ay2*cpb_y**2*spa_y*spb_y + Icliff(1,3,A_y)*3*B_y**2*ay2*cpb_y*spa_y*spb_y**2 + Icliff(0,4,A_y)*B_y**2*ay2*spa_y*spb_y**3 - Icliff(2,1,A_y)*3*B_y*ay2*cpa_y*cpb_y**2 - Icliff(1,2,A_y)*3*B_y*ay2*cpa_y*cpb_y*spb_y + Icliff(3,0,A_y)*3*B_y*ay2*cpa_y*cpb_y*spb_y + Icliff(2,1,A_y)*3*B_y*ay2*cpa_y*spb_y**2 - Icliff(1,2,A_y)*3*B_y*ay2*cpb_y**2*spa_y - Icliff(0,3,A_y)*3*B_y*ay2*cpb_y*spa_y*spb_y + Icliff(2,1,A_y)*3*B_y*ay2*cpb_y*spa_y*spb_y + Icliff(1,2,A_y)*3*B_y*ay2*spa_y*spb_y**2 - Icliff(2,0,A_y)*ay2*cpa_y*cpb_y - Icliff(1,1,A_y)*ay2*cpa_y*spb_y - Icliff(1,1,A_y)*ay2*cpb_y*spa_y - Icliff(0,2,A_y)*ay2*spa_y*spb_y)
        D_00_z = Icliff(2,0,A_z)*B_z**2*az2*cpb_z**2 + Icliff(1,1,A_z)*2*B_z**2*az2*cpb_z*spb_z + Icliff(0,2,A_z)*B_z**2*az2*spb_z**2 - Icliff(0,1,A_z)*B_z*az2*cpb_z + Icliff(1,0,A_z)*B_z*az2*spb_z
      
        K = -0.5d0 * ( D_00_x * S_11_y * S_00_z + S_00_x * D_11_y * S_00_z + S_00_x * S_11_y * D_00_z )
      
      case (23) ! | py    pz    ( 12 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=1)
      
        S_00_x = sx_int
        S_10_y = inv_ay * (Icliff(1,0,A_y)*cpa_y + Icliff(0,1,A_y)*spa_y)
        S_01_z = inv_az * (Icliff(1,0,A_z)*cpa_z + Icliff(0,1,A_z)*spa_z)
        D_00_x = Icliff(2,0,A_x)*B_x**2*ax2*cpb_x**2 + Icliff(1,1,A_x)*2*B_x**2*ax2*cpb_x*spb_x + Icliff(0,2,A_x)*B_x**2*ax2*spb_x**2 - Icliff(0,1,A_x)*B_x*ax2*cpb_x + Icliff(1,0,A_x)*B_x*ax2*spb_x
        D_10_y = inv_ay * (Icliff(3,0,A_y)*B_y**2*ay2*cpa_y*cpb_y**2 + Icliff(2,1,A_y)*2*B_y**2*ay2*cpa_y*cpb_y*spb_y + Icliff(1,2,A_y)*B_y**2*ay2*cpa_y*spb_y**2 + Icliff(2,1,A_y)*B_y**2*ay2*cpb_y**2*spa_y + Icliff(1,2,A_y)*2*B_y**2*ay2*cpb_y*spa_y*spb_y + Icliff(0,3,A_y)*B_y**2*ay2*spa_y*spb_x**2 - Icliff(1,1,A_y)*B_y*ay2*cpa_y*cpb_y + Icliff(2,0,A_y)*B_y*ay2*cpa_y*spb_y - Icliff(0,2,A_y)*B_y*ay2*cpb_y*spa_y + Icliff(1,1,A_y)*B_y*ay2*spa_y*spb_y)
        D_01_z = inv_az * (Icliff(3,0,A_z)*B_z**2*az2*cpb_z**3 + Icliff(2,1,A_z)*3*B_z**2*az2*cpb_z**2*spb_z + Icliff(1,2,A_z)*3*B_z**2*az2*cpb_z*spb_z**2 + Icliff(0,3,A_z)*B_z**2*az2*spb_z**3 - Icliff(1,1,A_z)*3*B_z*az2*cpb_z**2 - Icliff(0,2,A_z)*3*B_z*az2*cpb_z*spb_z + Icliff(2,0,A_z)*3*B_z*az2*cpb_z*spb_z + Icliff(1,1,A_z)*3*B_z*az2*spb_z**2 - Icliff(1,0,A_z)*az2*cpb_z - Icliff(0,1,A_z)*az2*spb_z)
      
        K = -0.5d0 * ( D_00_x * S_10_y * S_01_z + S_00_x * D_10_y * S_01_z + S_00_x * S_10_y * D_01_z )
      
      case (30) ! | pz    s     ( 13 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=0)
      
        S_00_x = sx_int
        S_00_y = sy_int
        S_10_z = inv_az * (Icliff(1,0,A_z)*cpa_z + Icliff(0,1,A_z)*spa_z)
        D_00_x = Icliff(2,0,A_x)*B_x**2*ax2*cpb_x**2 + Icliff(1,1,A_x)*2*B_x**2*ax2*cpb_x*spb_x + Icliff(0,2,A_x)*B_x**2*ax2*spb_x**2 - Icliff(0,1,A_x)*B_x*ax2*cpb_x + Icliff(1,0,A_x)*B_x*ax2*spb_x
        D_00_y = Icliff(2,0,A_y)*B_y**2*ay2*cpb_y**2 + Icliff(1,1,A_y)*2*B_y**2*ay2*cpb_y*spb_y + Icliff(0,2,A_y)*B_y**2*ay2*spb_y**2 - Icliff(0,1,A_y)*B_y*ay2*cpb_y + Icliff(1,0,A_y)*B_y*ay2*spb_y
        D_10_z = inv_az * (Icliff(3,0,A_z)*B_z**2*az2*cpa_z*cpb_z**2 + Icliff(2,1,A_z)*2*B_z**2*az2*cpa_z*cpb_z*spb_z + Icliff(1,2,A_z)*B_z**2*az2*cpa_z*spb_z**2 + Icliff(2,1,A_z)*B_z**2*az2*cpb_z**2*spa_z + Icliff(1,2,A_z)*2*B_z**2*az2*cpb_z*spa_z*spb_z + Icliff(0,3,A_z)*B_z**2*az2*spa_z*spb_x**2 - Icliff(1,1,A_z)*B_z*az2*cpa_z*cpb_z + Icliff(2,0,A_z)*B_z*az2*cpa_z*spb_z - Icliff(0,2,A_z)*B_z*az2*cpb_z*spa_z + Icliff(1,1,A_z)*B_z*az2*spa_z*spb_z)
      
        K = -0.5d0 * ( D_00_x * S_00_y * S_10_z + S_00_x * D_00_y * S_10_z + S_00_x * S_00_y * D_10_z )
      
      case (31) ! | pz    px    ( 14 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=1, m=0, n=0)
      
        S_01_x = inv_ax * (Icliff(1,0,A_x)*cpb_x + Icliff(0,1,A_x)*spb_x)
        S_00_y = sy_int
        S_10_z = inv_az * (Icliff(1,0,A_z)*cpa_z + Icliff(0,1,A_z)*spa_z)
        D_01_x = inv_ax * (Icliff(3,0,A_x)*B_x**2*ax2*cpb_x**3 + Icliff(2,1,A_x)*3*B_x**2*ax2*cpb_x**2*spb_x + Icliff(1,2,A_x)*3*B_x**2*ax2*cpb_x*spb_x**2 + Icliff(0,3,A_x)*B_x**2*ax2*spb_x**3 - Icliff(1,1,A_x)*3*B_x*ax2*cpb_x**2 - Icliff(0,2,A_x)*3*B_x*ax2*cpb_x*spb_x + Icliff(2,0,A_x)*3*B_x*ax2*cpb_x*spb_x + Icliff(1,1,A_x)*3*B_x*ax2*spb_x**2 - Icliff(1,0,A_x)*ax2*cpb_x - Icliff(0,1,A_x)*ax2*spb_x)
        D_00_y = Icliff(2,0,A_y)*B_y**2*ay2*cpb_y**2 + Icliff(1,1,A_y)*2*B_y**2*ay2*cpb_y*spb_y + Icliff(0,2,A_y)*B_y**2*ay2*spb_y**2 - Icliff(0,1,A_y)*B_y*ay2*cpb_y + Icliff(1,0,A_y)*B_y*ay2*spb_y
        D_10_z = inv_az * (Icliff(3,0,A_z)*B_z**2*az2*cpa_z*cpb_z**2 + Icliff(2,1,A_z)*2*B_z**2*az2*cpa_z*cpb_z*spb_z + Icliff(1,2,A_z)*B_z**2*az2*cpa_z*spb_z**2 + Icliff(2,1,A_z)*B_z**2*az2*cpb_z**2*spa_z + Icliff(1,2,A_z)*2*B_z**2*az2*cpb_z*spa_z*spb_z + Icliff(0,3,A_z)*B_z**2*az2*spa_z*spb_x**2 - Icliff(1,1,A_z)*B_z*az2*cpa_z*cpb_z + Icliff(2,0,A_z)*B_z*az2*cpa_z*spb_z - Icliff(0,2,A_z)*B_z*az2*cpb_z*spa_z + Icliff(1,1,A_z)*B_z*az2*spa_z*spb_z)
      
        K = -0.5d0 * ( D_01_x * S_00_y * S_10_z + S_01_x * D_00_y * S_10_z + S_01_x * S_00_y * D_10_z )
      
      case (32) ! | pz    py    ( 15 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=1, n=0)
      
        S_00_x = sx_int
        S_01_y = inv_ay * (Icliff(1,0,A_y)*cpb_y + Icliff(0,1,A_y)*spb_y)
        S_10_z = inv_az * (Icliff(1,0,A_z)*cpa_z + Icliff(0,1,A_z)*spa_z)
        D_00_x = Icliff(2,0,A_x)*B_x**2*ax2*cpb_x**2 + Icliff(1,1,A_x)*2*B_x**2*ax2*cpb_x*spb_x + Icliff(0,2,A_x)*B_x**2*ax2*spb_x**2 - Icliff(0,1,A_x)*B_x*ax2*cpb_x + Icliff(1,0,A_x)*B_x*ax2*spb_x
        D_01_y = inv_ay * (Icliff(3,0,A_y)*B_y**2*ay2*cpb_y**3 + Icliff(2,1,A_y)*3*B_y**2*ay2*cpb_y**2*spb_y + Icliff(1,2,A_y)*3*B_y**2*ay2*cpb_y*spb_y**2 + Icliff(0,3,A_y)*B_y**2*ay2*spb_y**3 - Icliff(1,1,A_y)*3*B_y*ay2*cpb_y**2 - Icliff(0,2,A_y)*3*B_y*ay2*cpb_y*spb_y + Icliff(2,0,A_y)*3*B_y*ay2*cpb_y*spb_y + Icliff(1,1,A_y)*3*B_y*ay2*spb_y**2 - Icliff(1,0,A_y)*ay2*cpb_y - Icliff(0,1,A_y)*ay2*spb_y)
        D_10_z = inv_az * (Icliff(3,0,A_z)*B_z**2*az2*cpa_z*cpb_z**2 + Icliff(2,1,A_z)*2*B_z**2*az2*cpa_z*cpb_z*spb_z + Icliff(1,2,A_z)*B_z**2*az2*cpa_z*spb_z**2 + Icliff(2,1,A_z)*B_z**2*az2*cpb_z**2*spa_z + Icliff(1,2,A_z)*2*B_z**2*az2*cpb_z*spa_z*spb_z + Icliff(0,3,A_z)*B_z**2*az2*spa_z*spb_x**2 - Icliff(1,1,A_z)*B_z*az2*cpa_z*cpb_z + Icliff(2,0,A_z)*B_z*az2*cpa_z*spb_z - Icliff(0,2,A_z)*B_z*az2*cpb_z*spa_z + Icliff(1,1,A_z)*B_z*az2*spa_z*spb_z)
      
        K = -0.5d0 * ( D_00_x * S_01_y * S_10_z + S_00_x * D_01_y * S_10_z + S_00_x * S_01_y * D_10_z )
      
      case (33) ! | pz    pz    ( 16 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=1)
      
        S_00_x = sx_int
        S_00_y = sy_int
        S_11_z = inv_az2 * (Icliff(2,0,A_z)*cpa_z*cpb_z + Icliff(1,1,A_z)*cpa_z*spb_z + Icliff(1,1,A_z)*cpb_z*spa_z + Icliff(0,2,A_z)*spa_z*spb_z)
        D_00_x = Icliff(2,0,A_x)*B_x**2*ax2*cpb_x**2 + Icliff(1,1,A_x)*2*B_x**2*ax2*cpb_x*spb_x + Icliff(0,2,A_x)*B_x**2*ax2*spb_x**2 - Icliff(0,1,A_x)*B_x*ax2*cpb_x + Icliff(1,0,A_x)*B_x*ax2*spb_x
        D_00_y = Icliff(2,0,A_y)*B_y**2*ay2*cpb_y**2 + Icliff(1,1,A_y)*2*B_y**2*ay2*cpb_y*spb_y + Icliff(0,2,A_y)*B_y**2*ay2*spb_y**2 - Icliff(0,1,A_y)*B_y*ay2*cpb_y + Icliff(1,0,A_y)*B_y*ay2*spb_y
        D_11_z = inv_az2 * (Icliff(4,0,A_z)*B_z**2*az2*cpa_z*cpb_z**3 + Icliff(3,1,A_z)*3*B_z**2*az2*cpa_z*cpb_z**2*spb_z + Icliff(2,2,A_z)*3*B_z**2*az2*cpa_z*cpb_z*spb_z**2 + Icliff(1,3,A_z)*B_z**2*az2*cpa_z*spb_z**3 + Icliff(3,1,A_z)*B_z**2*az2*cpb_z**3*spa_z + Icliff(2,2,A_z)*3*B_z**2*az2*cpb_z**2*spa_z*spb_z + Icliff(1,3,A_z)*3*B_z**2*az2*cpb_z*spa_z*spb_z**2 + Icliff(0,4,A_z)*B_z**2*az2*spa_z*spb_z**3 - Icliff(2,1,A_z)*3*B_z*az2*cpa_z*cpb_z**2 - Icliff(1,2,A_z)*3*B_z*az2*cpa_z*cpb_z*spb_z + Icliff(3,0,A_z)*3*B_z*az2*cpa_z*cpb_z*spb_z + Icliff(2,1,A_z)*3*B_z*az2*cpa_z*spb_z**2 - Icliff(1,2,A_z)*3*B_z*az2*cpb_z**2*spa_z - Icliff(0,3,A_z)*3*B_z*az2*cpb_z*spa_z*spb_z + Icliff(2,1,A_z)*3*B_z*az2*cpb_z*spa_z*spb_z + Icliff(1,2,A_z)*3*B_z*az2*spa_z*spb_z**2 - Icliff(2,0,A_z)*az2*cpa_z*cpb_z - Icliff(1,1,A_z)*az2*cpa_z*spb_z - Icliff(1,1,A_z)*az2*cpb_z*spa_z - Icliff(0,2,A_z)*az2*spa_z*spb_z)
      
        K = -0.5d0 * ( D_00_x * S_00_y * S_11_z + S_00_x * D_00_y * S_11_z + S_00_x * S_00_y * D_11_z )

      case default
        K = 0.0d0

      end select

        K_t  =  K_t + const * Kx * Ky * Kz * K 

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine kinetic_integral_toroidal_3D