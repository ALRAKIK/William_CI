subroutine kinetic_matrix_toroidal_1D(nBas,AO,kinetic)

      use files 
      use torus_init
      use atom_basis
      use classification_ERI
      use keywords
      use unitcell_module

      implicit none 

      !-----------------------------------------------------------------!

      integer                      :: i , j
      integer                      :: nBas

      type(ERI_function)           :: AO (nBas)
      type(ERI_function)           :: AO1 , AO2


      double precision             :: r1(3) , r2(3)

      integer                      :: pattern_id, encode_orbital_pattern_AO

      ! output ! 

      double precision,intent(out) :: kinetic(nBas,nBas)

      !-----------------------------------------------------------------!

      kinetic(:,:) = 0.d0

      
      do i = 1 , nfuc
        do j = 1 , nBas

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z


          pattern_id = encode_orbital_pattern_AO(AO1%orbital, AO2%orbital)

          call Kinetic_integral_toroidal_1D(pattern_id,r1,r2,AO1,AO2,kinetic(i,j))
          
        end do 
      end do 

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      do i = 1 , nfuc
        do j = 1 , nBas
          if (abs(kinetic(i,j)) < 1e-14) kinetic(i,j) = 0.d0 
        end do 
      end do 

      call fill_with_symmetry(nBas, kinetic)

      if (c_K) then 
        write(outfile,'(a)') ""
        write(outfile,'(a)') "------------------"
        write(outfile,'(a)') "The kinetic matrix"
        write(outfile,'(a)') "------------------"
        write(outfile,'(a)') ""
        do i = 1 , nBas - 1 
          do j = i , nBas
            write(outfile,'(i3,i3, 6x,1000(f16.12,2x))') i , j, kinetic(i,j)
          end do
        end do 
      end if 

end subroutine kinetic_matrix_toroidal_1D


subroutine kinetic_integral_toroidal_1D(pattern_id,r1,r2,AO1,AO2,K_t)

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
      double precision                 :: px
      double precision                 :: A , B , Br
      double precision                 :: xp
      double precision                 :: spa  , spb
      double precision                 :: cpa  , cpb

      !   Real   ! 

      double precision                 :: albe, inv_albe , mu 
      double precision                 :: yp   , zp
      double precision                 :: ypa  , ypb
      double precision                 :: zpa  , zpb

      
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

      inv_ax  = 1.d0 / ax 
      inv_ax2 = inv_ax * inv_ax 

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

          beta2  = beta * beta

          A      = 2.d0*px/ax2
          B      = 2.d0*beta/ax2

          kx     = dexp(-2.d0 * ( alpha + beta - px ) / ax2)

          sx_int = Lx * iv_scaled(0, A)

          spa    = dsin(ax*(xp-x1)) ; cpa    = dcos(ax*(xp-x1))
          spb    = dsin(ax*(xp-x2)) ; cpb    = dcos(ax*(xp-x2)) 

          ! ----------------- !
          !   Real Gaussian   !
          ! ----------------- !

          albe     = alpha + beta
          inv_albe = 1.d0 / albe
          mu       = alpha * beta * inv_albe
          yp       = (alpha * y1 + beta * y2) * inv_albe
          zp       = (alpha * z1 + beta * z2) * inv_albe

          ky       = dexp(- mu * ( Y * Y ))
          kz       = dexp(- mu * ( Z * Z ))

          sy_int   = dsqrt(pi * inv_albe)
          sz_int   = dsqrt(pi * inv_albe)

          Br       = beta

          ypa  = yp - y1
          ypb  = yp - y2

          zpa  = zp - z1 
          zpb  = zp - z2

      select case(pattern_id)

      case (00) ! | s     s     ( 1 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=0, n=0)
      
        S_00_x = sx_int
        S_00_y = sy_int
        S_00_z = sz_int
        D_00_x = Icliff(2,0,A)*B**2*ax2*cpb**2 + Icliff(1,1,A)*2*B**2*ax2*cpb*spb + Icliff(0,2,A)*B**2*ax2*spb**2 - Icliff(0,1,A)*B*ax2*cpb + Icliff(1,0,A)*B*ax2*spb
        D_00_y = sy_int*4*Br**2*ypb**2 + Ireal(1,albe)*8*Br**2*ypb + Ireal(2,albe)*4*Br**2 - sy_int*2*Br
        D_00_z = sz_int*4*Br**2*zpb**2 + Ireal(1,albe)*8*Br**2*zpb + Ireal(2,albe)*4*Br**2 - sz_int*2*Br
      
        K = -0.5d0 * ( D_00_x * S_00_y * S_00_z + S_00_x * D_00_y * S_00_z + S_00_x * S_00_y * D_00_z )
      
      case (01) ! | s     px    ( 2 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=1, m=0, n=0)
      
        S_01_x = inv_ax * (Icliff(1,0,A)*cpb + Icliff(0,1,A)*spb)
        S_00_y = sy_int
        S_00_z = sz_int
        D_01_x = inv_ax * (Icliff(3,0,A)*B**2*ax2*cpb**3 + Icliff(2,1,A)*3*B**2*ax2*cpb**2*spb + Icliff(1,2,A)*3*B**2*ax2*cpb*spb**2 + Icliff(0,3,A)*B**2*ax2*spb**3 - Icliff(1,1,A)*3*B*ax2*cpb**2 - Icliff(0,2,A)*3*B*ax2*cpb*spb + Icliff(2,0,A)*3*B*ax2*cpb*spb + Icliff(1,1,A)*3*B*ax2*spb**2 - Icliff(1,0,A)*ax2*cpb - Icliff(0,1,A)*ax2*spb)
        D_00_y = sy_int*4*Br**2*ypb**2 + Ireal(1,albe)*8*Br**2*ypb + Ireal(2,albe)*4*Br**2 - sy_int*2*Br
        D_00_z = sz_int*4*Br**2*zpb**2 + Ireal(1,albe)*8*Br**2*zpb + Ireal(2,albe)*4*Br**2 - sz_int*2*Br
      
        K = -0.5d0 * ( D_01_x * S_00_y * S_00_z + S_01_x * D_00_y * S_00_z + S_01_x * S_00_y * D_00_z )
      
      case (02) ! | s     py    ( 3 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        S_00_x = sx_int
        S_01_y = sy_int*ypb + Ireal(1,albe)
        S_00_z = sz_int
        D_00_x = Icliff(2,0,A)*B**2*ax2*cpb**2 + Icliff(1,1,A)*2*B**2*ax2*cpb*spb + Icliff(0,2,A)*B**2*ax2*spb**2 - Icliff(0,1,A)*B*ax2*cpb + Icliff(1,0,A)*B*ax2*spb
        D_01_y = sy_int*4*Br**2*ypb**3 + Ireal(1,albe)*12*Br**2*ypb**2 + Ireal(2,albe)*12*Br**2*ypb + Ireal(3,albe)*4*Br**2 - sy_int*6*Br*ypb - Ireal(1,albe)*6*Br
        D_00_z = sz_int*4*Br**2*zpb**2 + Ireal(1,albe)*8*Br**2*zpb + Ireal(2,albe)*4*Br**2 - sz_int*2*Br
      
        K = -0.5d0 * ( D_00_x * S_01_y * S_00_z + S_00_x * D_01_y * S_00_z + S_00_x * S_01_y * D_00_z )
      
      case (03) ! | s     pz    ( 4 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        S_00_x = sx_int
        S_00_y = sy_int
        S_01_z = sz_int*zpb + Ireal(1,albe)
        D_00_x = Icliff(2,0,A)*B**2*ax2*cpb**2 + Icliff(1,1,A)*2*B**2*ax2*cpb*spb + Icliff(0,2,A)*B**2*ax2*spb**2 - Icliff(0,1,A)*B*ax2*cpb + Icliff(1,0,A)*B*ax2*spb
        D_00_y = sy_int*4*Br**2*ypb**2 + Ireal(1,albe)*8*Br**2*ypb + Ireal(2,albe)*4*Br**2 - sy_int*2*Br
        D_01_z = sz_int*4*Br**2*zpb**3 + Ireal(1,albe)*12*Br**2*zpb**2 + Ireal(2,albe)*12*Br**2*zpb + Ireal(3,albe)*4*Br**2 - sz_int*6*Br*zpb - Ireal(1,albe)*6*Br
      
        K = -0.5d0 * ( D_00_x * S_00_y * S_01_z + S_00_x * D_00_y * S_01_z + S_00_x * S_00_y * D_01_z )
      
      case (10) ! | px    s     ( 5 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=0)
      
        S_10_x = inv_ax * (Icliff(1,0,A)*cpa + Icliff(0,1,A)*spa)
        S_00_y = sy_int
        S_00_z = sz_int
        D_10_x = inv_ax * (Icliff(3,0,A)*B**2*ax2*cpa*cpb**2 + Icliff(2,1,A)*2*B**2*ax2*cpa*cpb*spb + Icliff(1,2,A)*B**2*ax2*cpa*spb**2 + Icliff(2,1,A)*B**2*ax2*cpb**2*spa + Icliff(1,2,A)*2*B**2*ax2*cpb*spa*spb + Icliff(0,3,A)*B**2*ax2*spa*spb**2 - Icliff(1,1,A)*B*ax2*cpa*cpb + Icliff(2,0,A)*B*ax2*cpa*spb - Icliff(0,2,A)*B*ax2*cpb*spa + Icliff(1,1,A)*B*ax2*spa*spb)
        D_00_y = sy_int*4*Br**2*ypb**2 + Ireal(1,albe)*8*Br**2*ypb + Ireal(2,albe)*4*Br**2 - sy_int*2*Br
        D_00_z = sz_int*4*Br**2*zpb**2 + Ireal(1,albe)*8*Br**2*zpb + Ireal(2,albe)*4*Br**2 - sz_int*2*Br
      
        K = -0.5d0 * ( D_10_x * S_00_y * S_00_z + S_10_x * D_00_y * S_00_z + S_10_x * S_00_y * D_00_z )
      
      case (11) ! | px    px    ( 6 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=1, m=0, n=0)
      
        S_11_x = inv_ax2 * (Icliff(2,0,A)*cpa*cpb + Icliff(1,1,A)*cpa*spb + Icliff(1,1,A)*cpb*spa + Icliff(0,2,A)*spa*spb)
        S_00_y = sy_int
        S_00_z = sz_int
        D_11_x = inv_ax2 * (Icliff(4,0,A)*B**2*ax2*cpa*cpb**3 + Icliff(3,1,A)*3*B**2*ax2*cpa*cpb**2*spb + Icliff(2,2,A)*3*B**2*ax2*cpa*cpb*spb**2 + Icliff(1,3,A)*B**2*ax2*cpa*spb**3 + Icliff(3,1,A)*B**2*ax2*cpb**3*spa + Icliff(2,2,A)*3*B**2*ax2*cpb**2*spa*spb + Icliff(1,3,A)*3*B**2*ax2*cpb*spa*spb**2 + Icliff(0,4,A)*B**2*ax2*spa*spb**3 - Icliff(2,1,A)*3*B*ax2*cpa*cpb**2 - Icliff(1,2,A)*3*B*ax2*cpa*cpb*spb + Icliff(3,0,A)*3*B*ax2*cpa*cpb*spb + Icliff(2,1,A)*3*B*ax2*cpa*spb**2 - Icliff(1,2,A)*3*B*ax2*cpb**2*spa - Icliff(0,3,A)*3*B*ax2*cpb*spa*spb + Icliff(2,1,A)*3*B*ax2*cpb*spa*spb + Icliff(1,2,A)*3*B*ax2*spa*spb**2 - Icliff(2,0,A)*ax2*cpa*cpb - Icliff(1,1,A)*ax2*cpa*spb - Icliff(1,1,A)*ax2*cpb*spa - Icliff(0,2,A)*ax2*spa*spb)
        D_00_y = sy_int*4*Br**2*ypb**2 + Ireal(1,albe)*8*Br**2*ypb + Ireal(2,albe)*4*Br**2 - sy_int*2*Br
        D_00_z = sz_int*4*Br**2*zpb**2 + Ireal(1,albe)*8*Br**2*zpb + Ireal(2,albe)*4*Br**2 - sz_int*2*Br
      
        K = -0.5d0 * ( D_11_x * S_00_y * S_00_z + S_11_x * D_00_y * S_00_z + S_11_x * S_00_y * D_00_z )
      
      case (12) ! | px    py    ( 7 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        S_10_x = inv_ax * (Icliff(1,0,A)*cpa + Icliff(0,1,A)*spa)
        S_01_y = sy_int*ypb + Ireal(1,albe)
        S_00_z = sz_int
        D_10_x = inv_ax * (Icliff(3,0,A)*B**2*ax2*cpa*cpb**2 + Icliff(2,1,A)*2*B**2*ax2*cpa*cpb*spb + Icliff(1,2,A)*B**2*ax2*cpa*spb**2 + Icliff(2,1,A)*B**2*ax2*cpb**2*spa + Icliff(1,2,A)*2*B**2*ax2*cpb*spa*spb + Icliff(0,3,A)*B**2*ax2*spa*spb**2 - Icliff(1,1,A)*B*ax2*cpa*cpb + Icliff(2,0,A)*B*ax2*cpa*spb - Icliff(0,2,A)*B*ax2*cpb*spa + Icliff(1,1,A)*B*ax2*spa*spb)
        D_01_y = sy_int*4*Br**2*ypb**3 + Ireal(1,albe)*12*Br**2*ypb**2 + Ireal(2,albe)*12*Br**2*ypb + Ireal(3,albe)*4*Br**2 - sy_int*6*Br*ypb - Ireal(1,albe)*6*Br
        D_00_z = sz_int*4*Br**2*zpb**2 + Ireal(1,albe)*8*Br**2*zpb + Ireal(2,albe)*4*Br**2 - sz_int*2*Br
      
        K = -0.5d0 * ( D_10_x * S_01_y * S_00_z + S_10_x * D_01_y * S_00_z + S_10_x * S_01_y * D_00_z )
      
      case (13) ! | px    pz    ( 8 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        S_10_x = inv_ax * (Icliff(1,0,A)*cpa + Icliff(0,1,A)*spa)
        S_00_y = sy_int
        S_01_z = sz_int*zpb + Ireal(1,albe)
        D_10_x = inv_ax * (Icliff(3,0,A)*B**2*ax2*cpa*cpb**2 + Icliff(2,1,A)*2*B**2*ax2*cpa*cpb*spb + Icliff(1,2,A)*B**2*ax2*cpa*spb**2 + Icliff(2,1,A)*B**2*ax2*cpb**2*spa + Icliff(1,2,A)*2*B**2*ax2*cpb*spa*spb + Icliff(0,3,A)*B**2*ax2*spa*spb**2 - Icliff(1,1,A)*B*ax2*cpa*cpb + Icliff(2,0,A)*B*ax2*cpa*spb - Icliff(0,2,A)*B*ax2*cpb*spa + Icliff(1,1,A)*B*ax2*spa*spb)
        D_00_y = sy_int*4*Br**2*ypb**2 + Ireal(1,albe)*8*Br**2*ypb + Ireal(2,albe)*4*Br**2 - sy_int*2*Br
        D_01_z = sz_int*4*Br**2*zpb**3 + Ireal(1,albe)*12*Br**2*zpb**2 + Ireal(2,albe)*12*Br**2*zpb + Ireal(3,albe)*4*Br**2 - sz_int*6*Br*zpb - Ireal(1,albe)*6*Br
      
        K = -0.5d0 * ( D_10_x * S_00_y * S_01_z + S_10_x * D_00_y * S_01_z + S_10_x * S_00_y * D_01_z )
      
      case (20) ! | py    s     ( 9 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=0)
      
        S_00_x = sx_int
        S_10_y = sy_int*ypa + Ireal(1,albe)
        S_00_z = sz_int
        D_00_x = Icliff(2,0,A)*B**2*ax2*cpb**2 + Icliff(1,1,A)*2*B**2*ax2*cpb*spb + Icliff(0,2,A)*B**2*ax2*spb**2 - Icliff(0,1,A)*B*ax2*cpb + Icliff(1,0,A)*B*ax2*spb
        D_10_y = sy_int*4*Br**2*ypa*ypb**2 + Ireal(1,albe)*8*Br**2*ypa*ypb + Ireal(2,albe)*4*Br**2*ypa + Ireal(1,albe)*4*Br**2*ypb**2 + Ireal(2,albe)*8*Br**2*ypb + Ireal(3,albe)*4*Br**2 - sy_int*2*Br*ypa - Ireal(1,albe)*2*Br
        D_00_z = sz_int*4*Br**2*zpb**2 + Ireal(1,albe)*8*Br**2*zpb + Ireal(2,albe)*4*Br**2 - sz_int*2*Br
      
        K = -0.5d0 * ( D_00_x * S_10_y * S_00_z + S_00_x * D_10_y * S_00_z + S_00_x * S_10_y * D_00_z )
      
      case (21) ! | py    px    ( 10 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=1, m=0, n=0)
      
        S_01_x = inv_ax * (Icliff(1,0,A)*cpb + Icliff(0,1,A)*spb)
        S_10_y = sy_int*ypa + Ireal(1,albe)
        S_00_z = sz_int
        D_01_x = inv_ax * (Icliff(3,0,A)*B**2*ax2*cpb**3 + Icliff(2,1,A)*3*B**2*ax2*cpb**2*spb + Icliff(1,2,A)*3*B**2*ax2*cpb*spb**2 + Icliff(0,3,A)*B**2*ax2*spb**3 - Icliff(1,1,A)*3*B*ax2*cpb**2 - Icliff(0,2,A)*3*B*ax2*cpb*spb + Icliff(2,0,A)*3*B*ax2*cpb*spb + Icliff(1,1,A)*3*B*ax2*spb**2 - Icliff(1,0,A)*ax2*cpb - Icliff(0,1,A)*ax2*spb)
        D_10_y = sy_int*4*Br**2*ypa*ypb**2 + Ireal(1,albe)*8*Br**2*ypa*ypb + Ireal(2,albe)*4*Br**2*ypa + Ireal(1,albe)*4*Br**2*ypb**2 + Ireal(2,albe)*8*Br**2*ypb + Ireal(3,albe)*4*Br**2 - sy_int*2*Br*ypa - Ireal(1,albe)*2*Br
        D_00_z = sz_int*4*Br**2*zpb**2 + Ireal(1,albe)*8*Br**2*zpb + Ireal(2,albe)*4*Br**2 - sz_int*2*Br
      
        K = -0.5d0 * ( D_01_x * S_10_y * S_00_z + S_01_x * D_10_y * S_00_z + S_01_x * S_10_y * D_00_z )
      
      case (22) ! | py    py    ( 11 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=1, n=0)
      
        S_00_x = sx_int
        S_11_y = sy_int*ypa*ypb + Ireal(1,albe)*ypa + Ireal(1,albe)*ypb + Ireal(2,albe)
        S_00_z = sz_int
        D_00_x = Icliff(2,0,A)*B**2*ax2*cpb**2 + Icliff(1,1,A)*2*B**2*ax2*cpb*spb + Icliff(0,2,A)*B**2*ax2*spb**2 - Icliff(0,1,A)*B*ax2*cpb + Icliff(1,0,A)*B*ax2*spb
        D_11_y = sy_int*4*Br**2*ypa*ypb**3 + Ireal(1,albe)*12*Br**2*ypa*ypb**2 + Ireal(2,albe)*12*Br**2*ypa*ypb + Ireal(3,albe)*4*Br**2*ypa + Ireal(1,albe)*4*Br**2*ypb**3 + Ireal(2,albe)*12*Br**2*ypb**2 + Ireal(3,albe)*12*Br**2*ypb + Ireal(4,albe)*4*Br**2 - sy_int*6*Br*ypa*ypb - Ireal(1,albe)*6*Br*ypa - Ireal(1,albe)*6*Br*ypb - Ireal(2,albe)*6*Br
        D_00_z = sz_int*4*Br**2*zpb**2 + Ireal(1,albe)*8*Br**2*zpb + Ireal(2,albe)*4*Br**2 - sz_int*2*Br
      
        K = -0.5d0 * ( D_00_x * S_11_y * S_00_z + S_00_x * D_11_y * S_00_z + S_00_x * S_11_y * D_00_z )
      
      case (23) ! | py    pz    ( 12 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=1)
      
        S_00_x = sx_int
        S_10_y = sy_int*ypa + Ireal(1,albe)
        S_01_z = sz_int*zpb + Ireal(1,albe)
        D_00_x = Icliff(2,0,A)*B**2*ax2*cpb**2 + Icliff(1,1,A)*2*B**2*ax2*cpb*spb + Icliff(0,2,A)*B**2*ax2*spb**2 - Icliff(0,1,A)*B*ax2*cpb + Icliff(1,0,A)*B*ax2*spb
        D_10_y = sy_int*4*Br**2*ypa*ypb**2 + Ireal(1,albe)*8*Br**2*ypa*ypb + Ireal(2,albe)*4*Br**2*ypa + Ireal(1,albe)*4*Br**2*ypb**2 + Ireal(2,albe)*8*Br**2*ypb + Ireal(3,albe)*4*Br**2 - sy_int*2*Br*ypa - Ireal(1,albe)*2*Br
        D_01_z = sz_int*4*Br**2*zpb**3 + Ireal(1,albe)*12*Br**2*zpb**2 + Ireal(2,albe)*12*Br**2*zpb + Ireal(3,albe)*4*Br**2 - sz_int*6*Br*zpb - Ireal(1,albe)*6*Br
      
        K = -0.5d0 * ( D_00_x * S_10_y * S_01_z + S_00_x * D_10_y * S_01_z + S_00_x * S_10_y * D_01_z )
      
      case (30) ! | pz    s     ( 13 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=0)
      
        S_00_x = sx_int
        S_00_y = sy_int
        S_10_z = sz_int*zpa + Ireal(1,albe)
        D_00_x = Icliff(2,0,A)*B**2*ax2*cpb**2 + Icliff(1,1,A)*2*B**2*ax2*cpb*spb + Icliff(0,2,A)*B**2*ax2*spb**2 - Icliff(0,1,A)*B*ax2*cpb + Icliff(1,0,A)*B*ax2*spb
        D_00_y = sy_int*4*Br**2*ypb**2 + Ireal(1,albe)*8*Br**2*ypb + Ireal(2,albe)*4*Br**2 - sy_int*2*Br
        D_10_z = sz_int*4*Br**2*zpa*zpb**2 + Ireal(1,albe)*8*Br**2*zpa*zpb + Ireal(2,albe)*4*Br**2*zpa + Ireal(1,albe)*4*Br**2*zpb**2 + Ireal(2,albe)*8*Br**2*zpb + Ireal(3,albe)*4*Br**2 - sz_int*2*Br*zpa - Ireal(1,albe)*2*Br
      
        K = -0.5d0 * ( D_00_x * S_00_y * S_10_z + S_00_x * D_00_y * S_10_z + S_00_x * S_00_y * D_10_z )
      
      case (31) ! | pz    px    ( 14 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=1, m=0, n=0)
      
        S_01_x = inv_ax * (Icliff(1,0,A)*cpb + Icliff(0,1,A)*spb)
        S_00_y = sy_int
        S_10_z = sz_int*zpa + Ireal(1,albe)
        D_01_x = inv_ax * (Icliff(3,0,A)*B**2*ax2*cpb**3 + Icliff(2,1,A)*3*B**2*ax2*cpb**2*spb + Icliff(1,2,A)*3*B**2*ax2*cpb*spb**2 + Icliff(0,3,A)*B**2*ax2*spb**3 - Icliff(1,1,A)*3*B*ax2*cpb**2 - Icliff(0,2,A)*3*B*ax2*cpb*spb + Icliff(2,0,A)*3*B*ax2*cpb*spb + Icliff(1,1,A)*3*B*ax2*spb**2 - Icliff(1,0,A)*ax2*cpb - Icliff(0,1,A)*ax2*spb)
        D_00_y = sy_int*4*Br**2*ypb**2 + Ireal(1,albe)*8*Br**2*ypb + Ireal(2,albe)*4*Br**2 - sy_int*2*Br
        D_10_z = sz_int*4*Br**2*zpa*zpb**2 + Ireal(1,albe)*8*Br**2*zpa*zpb + Ireal(2,albe)*4*Br**2*zpa + Ireal(1,albe)*4*Br**2*zpb**2 + Ireal(2,albe)*8*Br**2*zpb + Ireal(3,albe)*4*Br**2 - sz_int*2*Br*zpa - Ireal(1,albe)*2*Br
      
        K = -0.5d0 * ( D_01_x * S_00_y * S_10_z + S_01_x * D_00_y * S_10_z + S_01_x * S_00_y * D_10_z )
      
      case (32) ! | pz    py    ( 15 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=1, n=0)
      
        S_00_x = sx_int
        S_01_y = sy_int*ypb + Ireal(1,albe)
        S_10_z = sz_int*zpa + Ireal(1,albe)
        D_00_x = Icliff(2,0,A)*B**2*ax2*cpb**2 + Icliff(1,1,A)*2*B**2*ax2*cpb*spb + Icliff(0,2,A)*B**2*ax2*spb**2 - Icliff(0,1,A)*B*ax2*cpb + Icliff(1,0,A)*B*ax2*spb
        D_01_y = sy_int*4*Br**2*ypb**3 + Ireal(1,albe)*12*Br**2*ypb**2 + Ireal(2,albe)*12*Br**2*ypb + Ireal(3,albe)*4*Br**2 - sy_int*6*Br*ypb - Ireal(1,albe)*6*Br
        D_10_z = sz_int*4*Br**2*zpa*zpb**2 + Ireal(1,albe)*8*Br**2*zpa*zpb + Ireal(2,albe)*4*Br**2*zpa + Ireal(1,albe)*4*Br**2*zpb**2 + Ireal(2,albe)*8*Br**2*zpb + Ireal(3,albe)*4*Br**2 - sz_int*2*Br*zpa - Ireal(1,albe)*2*Br
      
        K = -0.5d0 * ( D_00_x * S_01_y * S_10_z + S_00_x * D_01_y * S_10_z + S_00_x * S_01_y * D_10_z )
      
      case (33) ! | pz    pz    ( 16 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=1)
      
        S_00_x = sx_int
        S_00_y = sy_int
        S_11_z = sz_int*zpa*zpb + Ireal(1,albe)*zpa + Ireal(1,albe)*zpb + Ireal(2,albe)
        D_00_x = Icliff(2,0,A)*B**2*ax2*cpb**2 + Icliff(1,1,A)*2*B**2*ax2*cpb*spb + Icliff(0,2,A)*B**2*ax2*spb**2 - Icliff(0,1,A)*B*ax2*cpb + Icliff(1,0,A)*B*ax2*spb
        D_00_y = sy_int*4*Br**2*ypb**2 + Ireal(1,albe)*8*Br**2*ypb + Ireal(2,albe)*4*Br**2 - sy_int*2*Br
        D_11_z = sz_int*4*Br**2*zpa*zpb**3 + Ireal(1,albe)*12*Br**2*zpa*zpb**2 + Ireal(2,albe)*12*Br**2*zpa*zpb + Ireal(3,albe)*4*Br**2*zpa + Ireal(1,albe)*4*Br**2*zpb**3 + Ireal(2,albe)*12*Br**2*zpb**2 + Ireal(3,albe)*12*Br**2*zpb + Ireal(4,albe)*4*Br**2 - sz_int*6*Br*zpa*zpb - Ireal(1,albe)*6*Br*zpa - Ireal(1,albe)*6*Br*zpb - Ireal(2,albe)*6*Br
      
        K = -0.5d0 * ( D_00_x * S_00_y * S_11_z + S_00_x * D_00_y * S_11_z + S_00_x * S_00_y * D_11_z )

      case default
        K = 0.0d0

      end select

        K_t  =  K_t + const * Kx * Ky * Kz * K 

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine kinetic_integral_toroidal_1D