subroutine overlap_matrix_toroidal(nBas,AO,overlap)

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

          call overlap_integral_toroidal_1D(pattern_id,r1,r2,AO1,AO2,overlap(i,j))
          
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

end subroutine overlap_matrix_toroidal

      !-----------------------------------------------------------------!

subroutine overlap_integral_toroidal_1D(pattern_id,r1,r2,AO1,AO2,S)

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
      double precision                 :: px , xp
      double precision                 :: A
      double precision                 :: spa  , spb
      double precision                 :: cpa  , cpb

      !   Real   ! 

      double precision                 :: albe, inv_albe , mu 
      double precision                 :: yp   , zp
      double precision                 :: ypa  , ypb 
      double precision                 :: zpa  , zpb 
      
      integer                          :: pattern_id

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
 
      s       = 0.d0

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

          A      = 2.d0*px/ax2

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

          ypa = yp - y1
          ypb = yp - y2

          zpa = zp - z1 
          zpb = zp - z2
          

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
      
        sx   = inv_ax * (Icliff(1,0,A)*cpb + Icliff(0,1,A)*spb)
        sy   = Sy_int
        sz   = Sz_int
      
      case (02) ! | s     py    ( 3 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        sx   = Sx_int
        sy   = Sy_int*ypb + Ireal(1,albe)
        sz   = Sz_int
      
      case (03) ! | s     pz    ( 4 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        sx   = Sx_int
        sy   = Sy_int
        sz   = Sz_int*zpb + Ireal(1,albe)
      
      case (10) ! | px    s     ( 5 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=0)
      
        sx   = inv_ax * (Icliff(1,0,A)*cpa + Icliff(0,1,A)*spa)
        sy   = Sy_int
        sz   = Sz_int
      
      case (11) ! | px    px    ( 6 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=1, m=0, n=0)
      
        sx   = inv_ax2 * (Icliff(2,0,A)*cpa*cpb + Icliff(1,1,A)*cpa*spb + Icliff(1,1,A)*cpb*spa + Icliff(0,2,A)*spa*spb)
        sy   = Sy_int
        sz   = Sz_int
      
      case (12) ! | px    py    ( 7 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        sx   = inv_ax * (Icliff(1,0,A)*cpa + Icliff(0,1,A)*spa)
        sy   = Sy_int*ypb + Ireal(1,albe)
        sz   = Sz_int
      
      case (13) ! | px    pz    ( 8 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        sx   = inv_ax * (Icliff(1,0,A)*cpa + Icliff(0,1,A)*spa)
        sy   = Sy_int
        sz   = Sz_int*zpb + Ireal(1,albe)
      
      case (20) ! | py    s     ( 9 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=0)
      
        sx   = Sx_int
        sy   = Sy_int*ypa + Ireal(1,albe)
        sz   = Sz_int
      
      case (21) ! | py    px    ( 10 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=1, m=0, n=0)
      
        sx   = inv_ax * (Icliff(1,0,A)*cpb + Icliff(0,1,A)*spb)
        sy   = Sy_int*ypa + Ireal(1,albe)
        sz   = Sz_int
      
      case (22) ! | py    py    ( 11 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=1, n=0)
      
        sx   = Sx_int
        sy   = Sy_int*ypa*ypb + Ireal(1,albe)*ypa + Ireal(1,albe)*ypb + Ireal(2,albe)
        sz   = Sz_int
      
      case (23) ! | py    pz    ( 12 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=1)
      
        sx   = Sx_int
        sy   = Sy_int*ypa + Ireal(1,albe)
        sz   = Sz_int*zpb + Ireal(1,albe)
      
      case (30) ! | pz    s     ( 13 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=0)
      
        sx   = Sx_int
        sy   = Sy_int
        sz   = Sz_int*zpa + Ireal(1,albe)
      
      case (31) ! | pz    px    ( 14 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=1, m=0, n=0)
      
        sx   = inv_ax * (Icliff(1,0,A)*cpb + Icliff(0,1,A)*spb)
        sy   = Sy_int
        sz   = Sz_int*zpa + Ireal(1,albe)
      
      case (32) ! | pz    py    ( 15 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=1, n=0)
      
        sx   = Sx_int
        sy   = Sy_int*ypb + Ireal(1,albe)
        sz   = Sz_int*zpa + Ireal(1,albe)
      
      case (33) ! | pz    pz    ( 16 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=1)
      
        sx   = Sx_int
        sy   = Sy_int
        sz   = Sz_int*zpa*zpb + Ireal(1,albe)*zpa + Ireal(1,albe)*zpb + Ireal(2,albe)

      case default

        s = 0.0d0

      end select

        s  =  s + const * kx * ky * kz * sx * sy * sz

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine overlap_integral_toroidal_1D