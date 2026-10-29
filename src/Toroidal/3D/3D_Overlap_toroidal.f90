subroutine overlap_matrix_toroidal_3D(nBas,AO,overlap)

      use files
      use atom_basis
      use torus_init
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

      double precision,intent(out) :: overlap(nBas,nBas)

      !-----------------------------------------------------------------!

      overlap(:,:) = 0.d0
      
      ! --------------------------------------------------------------- !

      do i = 1 , nfuc
        do j = 1 , nBas

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z


          pattern_id = encode_orbital_pattern_AO(AO1%orbital, AO2%orbital)

          call overlap_integral_toroidal_3D(pattern_id,r1,r2,AO1,AO2,overlap(i,j))
          
        end do 
      end do 

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
       
      do i = 1 , nfuc
        do j = 1 , nBas
          if (abs(overlap(i,j)) < 1e-14) overlap(i,j) = 0.d0 
        end do 
      end do 

      call fill_with_symmetry(nBas, overlap)
 
      if (c_OV) then 
        write(outfile,'(a)') ""
        write(outfile,'(a)') "------------------"
        write(outfile,'(a)') "The overlap matrix"
        write(outfile,'(a)') "------------------"
        write(outfile,'(a)') ""
        do i = 1 , nBas - 1 
          do j = i , nBas
            write(outfile,'(i3,i3, 6x,1000(f16.12,2x))') i , j, overlap(i,j)
          end do
        end do 
      end if 

end subroutine overlap_matrix_toroidal_3D


subroutine overlap_integral_toroidal_3D(pattern_id,r1,r2,AO1,AO2,S)


      use torus_init
      use classification_ERI
      use bessel_functions
      use functions

      implicit none 

      !-----------------------------------------------------------------!

      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: S

      integer                          :: i , j 
      double precision,parameter       :: pi = 3.14159265358979323846D00
      double precision                 :: alpha  , beta
      double precision                 :: c1     , c2     , const
      double precision                 :: x1     , x2     , X  
      double precision                 :: y1     , y2     , Y
      double precision                 :: z1     , z2     , Z 
      double precision                 :: kx     , ky     , kz
      double precision                 :: sx     , sy     , sz  
      double precision                 :: sx_int , sy_int , sz_int

      ! Clifford ! 

      double precision                 :: ax2 , inv_ax , inv_ax2
      double precision                 :: ay2 , inv_ay , inv_ay2 
      double precision                 :: az2 , inv_az , inv_az2

      double precision                 :: px , xp
      double precision                 :: py , yp
      double precision                 :: pz , zp

      double precision                 :: A_x , A_y , A_z
      double precision                 :: spa_x  , spb_x
      double precision                 :: cpa_x  , cpb_x
      double precision                 :: spa_y  , spb_y
      double precision                 :: cpa_y  , cpb_y
      double precision                 :: spa_z  , spb_z
      double precision                 :: cpa_z  , cpb_z

      integer                          :: pattern_id

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

      s  = 0.d0 
      sx = 0.d0
      sy = 0.d0 
      sz = 0.d0 

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

          A_x    = 2.d0*px/ax2
          A_y    = 2.d0*py/ay2
          A_z    = 2.d0*pz/az2

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
      
        sx   = Sx_int
        sy   = Sy_int
        sz   = Sz_int
      
      case (01) ! | s     px    ( 2 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=1, m=0, n=0)
      
        sx   = inv_ax * (Icliff(1,0,A_x)*cpb_x + Icliff(0,1,A_x)*spb_x)
        sy   = Sy_int
        sz   = Sz_int
      
      case (02) ! | s     py    ( 3 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        sx   = Sx_int
        sy   = inv_ay * (Icliff(1,0,A_y)*cpb_y + Icliff(0,1,A_y)*spb_y)
        sz   = Sz_int
      
      case (03) ! | s     pz    ( 4 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        sx   = Sx_int
        sy   = Sy_int
        sz   = inv_az * (Icliff(1,0,A_z)*cpb_z + Icliff(0,1,A_z)*spb_z)
      
      case (10) ! | px    s     ( 5 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=0)
      
        sx   = inv_ax * (Icliff(1,0,A_x)*cpa_x + Icliff(0,1,A_x)*spa_x)
        sy   = Sy_int
        sz   = Sz_int
      
      case (11) ! | px    px    ( 6 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=1, m=0, n=0)
      
        sx   = inv_ax2 * (Icliff(2,0,A_x)*cpa_x*cpb_x + Icliff(1,1,A_x)*cpa_x*spb_x + Icliff(1,1,A_x)*cpb_x*spa_x + Icliff(0,2,A_x)*spa_x*spb_x)
        sy   = Sy_int
        sz   = Sz_int
      
      case (12) ! | px    py    ( 7 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        sx   = inv_ax * (Icliff(1,0,A_x)*cpa_x + Icliff(0,1,A_x)*spa_x)
        sy   = inv_ay * (Icliff(1,0,A_y)*cpb_y + Icliff(0,1,A_y)*spb_y)
        sz   = Sz_int
      
      case (13) ! | px    pz    ( 8 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        sx   = inv_ax * (Icliff(1,0,A_x)*cpa_x + Icliff(0,1,A_x)*spa_x)
        sy   = Sy_int
        sz   = inv_az * (Icliff(1,0,A_z)*cpb_z + Icliff(0,1,A_z)*spb_z)
      
      case (20) ! | py    s     ( 9 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=0)
      
        sx   = Sx_int
        sy   = inv_ay * (Icliff(1,0,A_y)*cpa_y + Icliff(0,1,A_y)*spa_y)
        sz   = Sz_int
      
      case (21) ! | py    px    ( 10 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=1, m=0, n=0)
      
        sx   = inv_ax * (Icliff(1,0,A_x)*cpb_x + Icliff(0,1,A_x)*spb_x)
        sy   = inv_ay * (Icliff(1,0,A_y)*cpa_y + Icliff(0,1,A_y)*spa_y)
        sz   = Sz_int
      
      case (22) ! | py    py    ( 11 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=1, n=0)
      
        sx   = Sx_int
        sy   = inv_ay2 * (Icliff(2,0,A_y)*cpa_y*cpb_y + Icliff(1,1,A_y)*cpa_y*spb_y + Icliff(1,1,A_y)*cpb_y*spa_y + Icliff(0,2,A_y)*spa_y*spb_y)
        sz   = Sz_int
      
      case (23) ! | py    pz    ( 12 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=1)
      
        sx   = Sx_int
        sy   = inv_ay * (Icliff(1,0,A_y)*cpa_y + Icliff(0,1,A_y)*spa_y)
        sz   = inv_az * (Icliff(1,0,A_z)*cpb_z + Icliff(0,1,A_z)*spb_z)
      
      case (30) ! | pz    s     ( 13 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=0)
      
        sx   = Sx_int
        sy   = Sy_int
        sz   = inv_az * (Icliff(1,0,A_z)*cpa_z + Icliff(0,1,A_z)*spa_z)
      
      case (31) ! | pz    px    ( 14 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=1, m=0, n=0)
      
        sx   = inv_ax * (Icliff(1,0,A_x)*cpb_x + Icliff(0,1,A_x)*spb_x)
        sy   = Sy_int
        sz   = inv_az * (Icliff(1,0,A_z)*cpa_z + Icliff(0,1,A_z)*spa_z)
      
      case (32) ! | pz    py    ( 15 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=1, n=0)
      
        sx   = Sx_int
        sy   = inv_ay * (Icliff(1,0,A_y)*cpb_y + Icliff(0,1,A_y)*spb_y)
        sz   = inv_az * (Icliff(1,0,A_z)*cpa_z + Icliff(0,1,A_z)*spa_z)
      
      case (33) ! | pz    pz    ( 16 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=1)
      
        sx   = Sx_int
        sy   = Sy_int
        sz   = inv_az2 * (Icliff(2,0,A_z)*cpa_z*cpb_z + Icliff(1,1,A_z)*cpa_z*spb_z + Icliff(1,1,A_z)*cpb_z*spa_z + Icliff(0,2,A_z)*spa_z*spb_z)

      case default

        s = 0.0d0

      end select

        s  =  s + const * kx * ky * kz * sx * sy * sz

        end do 
      end do

      
end subroutine overlap_integral_toroidal_3D