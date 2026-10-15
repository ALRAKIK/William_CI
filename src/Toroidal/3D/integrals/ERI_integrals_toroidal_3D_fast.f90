subroutine ERI_integral_4_function_toroidal_3D_fast(one,two,three,four,value)
      
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
 
              call integrate_ERI_3D_fast(pattern_id ,mu_x,nu_x,phix,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,&
                                                     mu_y,nu_y,phiy,ypA,ypB,yqC,yqD,ya,yb,yc,yd,yp,yq,&
                                                     mu_z,nu_z,phiz,zpA,zpB,zqC,zqD,za,zb,zc,zd,zp,zq,value_s)

              value  = value    + const * const_x * const_y * const_z * value_s

            end do
          end do 
        end do 
      end do 

end subroutine ERI_integral_4_function_toroidal_3D_fast

subroutine integrate_ERI_3D_fast(pattern_id, p_x,q_x,phi_x,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,&
                                             p_y,q_y,phi_y,ypA,ypB,yqC,yqD,ya,yb,yc,yd,yp,yq,&
                                             p_z,q_z,phi_z,zpA,zpB,zqC,zqD,za,zb,zc,zd,zp,zq,result)

      use iso_c_binding
      use torus_init
      use precomputed_bessel
      use table_1d_lookup
      use constants_module
      use table_lookup_module

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

      double precision                     :: F_x(0:50000)
      double precision                     :: F_y(0:50000)
      double precision                     :: F_z(0:50000)

      double precision                     :: real_An, real_Bn
      double precision                     :: imag_An, imag_Bn
      integer                              :: n


      integer                              :: Nmax_x(nquad_point)
      integer                              :: Nmax_y(nquad_point)
      integer                              :: Nmax_z(nquad_point)

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

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      do n = 1, nquad_point
        CCx = inv_ax2 * 2.d0 * t_val(n) * t_val(n) 
        CCy = inv_ay2 * 2.d0 * t_val(n) * t_val(n) 
        CCz = inv_az2 * 2.d0 * t_val(n) * t_val(n) 
        Nmax_x(n) = get_Nmax(AAx, BBx, CCx, phi_x) + 10
        Nmax_y(n) = get_Nmax(AAy, BBy, CCy, phi_y) + 10
        Nmax_z(n) = get_Nmax(AAz, BBz, CCz, phi_z) + 10
      end do

      select case (pattern_id)

      case (0000) ! | s   s   s   s    ( 1 ) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0001) ! | s   s   s   px   ( 2 ) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0002) ! | s   s   s   py   ( 3 ) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0003) ! | s   s   s   pz   ( 4 ) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0010) ! | s   s   px  s    ( 5 ) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0011) ! | s   s   px  px   ( 6 ) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0012) ! | s   s   px  py   ( 7 ) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0013) ! | s   s   px  pz   ( 8 ) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0020) ! | s   s   py  s    ( 9 ) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0021) ! | s   s   py  px   ( 10) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0022) ! | s   s   py  py   ( 11) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0023) ! | s   s   py  pz   ( 12) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0030) ! | s   s   pz  s    ( 13) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0031) ! | s   s   pz  px   ( 14) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0032) ! | s   s   pz  py   ( 15) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0033) ! | s   s   pz  pz   ( 16) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0100) ! | s   px  s   s    ( 17) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0101) ! | s   px  s   px   ( 18) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0102) ! | s   px  s   py   ( 19) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0103) ! | s   px  s   pz   ( 20) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0110) ! | s   px  px  s    ( 21) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0111) ! | s   px  px  px   ( 22) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0112) ! | s   px  px  py   ( 23) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0113) ! | s   px  px  pz   ( 24) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0120) ! | s   px  py  s    ( 25) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0121) ! | s   px  py  px   ( 26) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0122) ! | s   px  py  py   ( 27) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0123) ! | s   px  py  pz   ( 28) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0130) ! | s   px  pz  s    ( 29) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0131) ! | s   px  pz  px   ( 30) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0132) ! | s   px  pz  py   ( 31) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0133) ! | s   px  pz  pz   ( 32) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0200) ! | s   py  s   s    ( 33) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0201) ! | s   py  s   px   ( 34) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0202) ! | s   py  s   py   ( 35) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0203) ! | s   py  s   pz   ( 36) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0210) ! | s   py  px  s    ( 37) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0211) ! | s   py  px  px   ( 38) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0212) ! | s   py  px  py   ( 39) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0213) ! | s   py  px  pz   ( 40) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0220) ! | s   py  py  s    ( 41) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0221) ! | s   py  py  px   ( 42) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0222) ! | s   py  py  py   ( 43) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0223) ! | s   py  py  pz   ( 44) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0230) ! | s   py  pz  s    ( 45) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0231) ! | s   py  pz  px   ( 46) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0232) ! | s   py  pz  py   ( 47) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0233) ! | s   py  pz  pz   ( 48) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0300) ! | s   pz  s   s    ( 49) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0301) ! | s   pz  s   px   ( 50) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0302) ! | s   pz  s   py   ( 51) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0303) ! | s   pz  s   pz   ( 52) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0310) ! | s   pz  px  s    ( 53) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0311) ! | s   pz  px  px   ( 54) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0312) ! | s   pz  px  py   ( 55) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0313) ! | s   pz  px  pz   ( 56) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0320) ! | s   pz  py  s    ( 57) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0321) ! | s   pz  py  px   ( 58) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0322) ! | s   pz  py  py   ( 59) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0323) ! | s   pz  py  pz   ( 60) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0330) ! | s   pz  pz  s    ( 61) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0331) ! | s   pz  pz  px   ( 62) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0332) ! | s   pz  pz  py   ( 63) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (0333) ! | s   pz  pz  pz   ( 64) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1000) ! | px  s   s   s    ( 65) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1001) ! | px  s   s   px   ( 66) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1002) ! | px  s   s   py   ( 67) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1003) ! | px  s   s   pz   ( 68) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1010) ! | px  s   px  s    ( 69) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1011) ! | px  s   px  px   ( 70) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1012) ! | px  s   px  py   ( 71) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1013) ! | px  s   px  pz   ( 72) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1020) ! | px  s   py  s    ( 73) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1021) ! | px  s   py  px   ( 74) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1022) ! | px  s   py  py   ( 75) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1023) ! | px  s   py  pz   ( 76) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1030) ! | px  s   pz  s    ( 77) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1031) ! | px  s   pz  px   ( 78) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1032) ! | px  s   pz  py   ( 79) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1033) ! | px  s   pz  pz   ( 80) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1100) ! | px  px  s   s    ( 81) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (1101) ! | px  px  s   px   ( 82) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1102) ! | px  px  s   py   ( 83) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1103) ! | px  px  s   pz   ( 84) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1110) ! | px  px  px  s    ( 85) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1111) ! | px  px  px  px   ( 86) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1112) ! | px  px  px  py   ( 87) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1113) ! | px  px  px  pz   ( 88) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1120) ! | px  px  py  s    ( 89) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1121) ! | px  px  py  px   ( 90) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1122) ! | px  px  py  py   ( 91) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1123) ! | px  px  py  pz   ( 92) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1130) ! | px  px  pz  s    ( 93) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1131) ! | px  px  pz  px   ( 94) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1132) ! | px  px  pz  py   ( 95) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1133) ! | px  px  pz  pz   ( 96) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = -1.0*inv_ax2*(0.5*c2xpab*I_A_x(n)+0.25*c2xpab*I_A_x(n+2)+0.25*c2xpab*I_A_x(Abs(n-2))-1.0*cxpa*cxpb*I_A_x(n))
        real_Bn  = I_B_x(n)
        if (dabs(AAx) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ax2*n*s2xpab*(I_A_x(n+1)+I_A_x(Abs(n-1)))/AAx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1200) ! | px  py  s   s    ( 97) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1201) ! | px  py  s   px   ( 98) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1202) ! | px  py  s   py   ( 99) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1203) ! | px  py  s   pz   ( 100) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1210) ! | px  py  px  s    ( 101) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1211) ! | px  py  px  px   ( 102) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1212) ! | px  py  px  py   ( 103) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1213) ! | px  py  px  pz   ( 104) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1220) ! | px  py  py  s    ( 105) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1221) ! | px  py  py  px   ( 106) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1222) ! | px  py  py  py   ( 107) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1223) ! | px  py  py  pz   ( 108) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1230) ! | px  py  pz  s    ( 109) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1231) ! | px  py  pz  px   ( 110) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1232) ! | px  py  pz  py   ( 111) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1233) ! | px  py  pz  pz   ( 112) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1300) ! | px  pz  s   s    ( 113) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1301) ! | px  pz  s   px   ( 114) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1302) ! | px  pz  s   py   ( 115) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1303) ! | px  pz  s   pz   ( 116) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1310) ! | px  pz  px  s    ( 117) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1311) ! | px  pz  px  px   ( 118) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1312) ! | px  pz  px  py   ( 119) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1313) ! | px  pz  px  pz   ( 120) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1320) ! | px  pz  py  s    ( 121) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
          
      case (1321) ! | px  pz  py  px   ( 122) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1322) ! | px  pz  py  py   ( 123) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1323) ! | px  pz  py  pz   ( 124) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1330) ! | px  pz  pz  s    ( 125) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1331) ! | px  pz  pz  px   ( 126) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1332) ! | px  pz  pz  py   ( 127) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (1333) ! | px  pz  pz  pz   ( 128) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpa*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpa*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2000) ! | py  s   s   s    ( 129) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2001) ! | py  s   s   px   ( 130) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2002) ! | py  s   s   py   ( 131) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2003) ! | py  s   s   pz   ( 132) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2010) ! | py  s   px  s    ( 133) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2011) ! | py  s   px  px   ( 134) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2012) ! | py  s   px  py   ( 135) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2013) ! | py  s   px  pz   ( 136) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2020) ! | py  s   py  s    ( 137) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2021) ! | py  s   py  px   ( 138) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2022) ! | py  s   py  py   ( 139) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2023) ! | py  s   py  pz   ( 140) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2030) ! | py  s   pz  s    ( 141) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2031) ! | py  s   pz  px   ( 142) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2032) ! | py  s   pz  py   ( 143) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2033) ! | py  s   pz  pz   ( 144) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2100) ! | py  px  s   s    ( 145) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2101) ! | py  px  s   px   ( 146) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2102) ! | py  px  s   py   ( 147) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2103) ! | py  px  s   pz   ( 148) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2110) ! | py  px  px  s    ( 149) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2111) ! | py  px  px  px   ( 150) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2112) ! | py  px  px  py   ( 151) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2113) ! | py  px  px  pz   ( 152) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2120) ! | py  px  py  s    ( 153) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2121) ! | py  px  py  px   ( 154) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2122) ! | py  px  py  py   ( 155) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2123) ! | py  px  py  pz   ( 156) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2130) ! | py  px  pz  s    ( 157) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2131) ! | py  px  pz  px   ( 158) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2132) ! | py  px  pz  py   ( 159) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2133) ! | py  px  pz  pz   ( 160) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !
      
      case (2200) ! | py  py  s   s    ( 161) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do
      
      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = I_B_y(n)
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do
      
      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2201) ! | py  py  s   px   ( 162) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = I_B_y(n)
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2202) ! | py  py  s   py   ( 163) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2203) ! | py  py  s   pz   ( 164) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = I_B_y(n)
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2210) ! | py  py  px  s    ( 165) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = I_B_y(n)
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2211) ! | py  py  px  px   ( 166) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = I_B_y(n)
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2212) ! | py  py  px  py   ( 167) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2213) ! | py  py  px  pz   ( 168) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = I_B_y(n)
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2220) ! | py  py  py  s    ( 169) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2221) ! | py  py  py  px   ( 170) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2222) ! | py  py  py  py   ( 171) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = I_B_z(n)
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2223) ! | py  py  py  pz   ( 172) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2230) ! | py  py  pz  s    ( 173) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = I_B_y(n)
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2231) ! | py  py  pz  px   ( 174) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = I_B_y(n)
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2232) ! | py  py  pz  py   ( 175) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2233) ! | py  py  pz  pz   ( 176) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = -1.0*inv_ay2*(0.5*c2ypab*I_A_y(n)+0.25*c2ypab*I_A_y(n+2)+0.25*c2ypab*I_A_y(Abs(n-2))-1.0*cypa*cypb*I_A_y(n))
        real_Bn  = I_B_y(n)
        if (dabs(AAy) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_ay2*n*s2ypab*(I_A_y(n+1)+I_A_y(Abs(n-1)))/AAy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = I_A_z(n)
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( real_An * imag_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2300) ! | py  pz  s   s    ( 177) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2301) ! | py  pz  s   px   ( 178) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2302) ! | py  pz  s   py   ( 179) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2303) ! | py  pz  s   pz   ( 180) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2310) ! | py  pz  px  s    ( 181) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2311) ! | py  pz  px  px   ( 182) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2312) ! | py  pz  px  py   ( 183) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2313) ! | py  pz  px  pz   ( 184) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2320) ! | py  pz  py  s    ( 185) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2321) ! | py  pz  py  px   ( 186) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2322) ! | py  pz  py  py   ( 187) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2323) ! | py  pz  py  pz   ( 188) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2330) ! | py  pz  pz  s    ( 189) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2331) ! | py  pz  pz  px   ( 190) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2332) ! | py  pz  pz  py   ( 191) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (2333) ! | py  pz  pz  pz   ( 192) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypa*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypa*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpb*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        imag_An  = 0.5*czpb*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3000) ! | pz  s   s   s    ( 193) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3001) ! | pz  s   s   px   ( 194) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3002) ! | pz  s   s   py   ( 195) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3003) ! | pz  s   s   pz   ( 196) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3010) ! | pz  s   px  s    ( 197) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3011) ! | pz  s   px  px   ( 198) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3012) ! | pz  s   px  py   ( 199) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3013) ! | pz  s   px  pz   ( 200) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3020) ! | pz  s   py  s    ( 201) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3021) ! | pz  s   py  px   ( 202) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3022) ! | pz  s   py  py   ( 203) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if  
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3023) ! | pz  s   py  pz   ( 204) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3030) ! | pz  s   pz  s    ( 205) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3031) ! | pz  s   pz  px   ( 206) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3032) ! | pz  s   pz  py   ( 207) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3033) ! | pz  s   pz  pz   ( 208) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3100) ! | pz  px  s   s    ( 209) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3101) ! | pz  px  s   px   ( 210) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3102) ! | pz  px  s   py   ( 211) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3103) ! | pz  px  s   pz   ( 212) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3110) ! | pz  px  px  s    ( 213) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3111) ! | pz  px  px  px   ( 214) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3112) ! | pz  px  px  py   ( 215) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3113) ! | pz  px  px  pz   ( 216) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3120) ! | pz  px  py  s    ( 217) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3121) ! | pz  px  py  px   ( 218) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3122) ! | pz  px  py  py   ( 219) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3123) ! | pz  px  py  pz   ( 220) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3130) ! | pz  px  pz  s    ( 221) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3131) ! | pz  px  pz  px   ( 222) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_x) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3132) ! | pz  px  pz  py   ( 223) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3133) ! | pz  px  pz  pz   ( 224) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = 0.5*inv_ax*sxpb*(I_A_x(n+1)+I_A_x(Abs(n-1)))
        real_Bn  = I_B_x(n)
        imag_An  = 0.5*cxpb*inv_ax*(-I_A_x(n+1)+I_A_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( imag_An * real_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3200) ! | pz  py  s   s    ( 225) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3201) ! | pz  py  s   px   ( 226) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3202) ! | pz  py  s   py   ( 227) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3203) ! | pz  py  s   pz   ( 228) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3210) ! | pz  py  px  s    ( 229) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3211) ! | pz  py  px  px   ( 230) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3212) ! | pz  py  px  py   ( 231) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3213) ! | pz  py  px  pz   ( 232) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3220) ! | pz  py  py  s    ( 233) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3221) ! | pz  py  py  px   ( 234) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3222) ! | pz  py  py  py   ( 235) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = I_B_z(n)
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3223) ! | pz  py  py  pz   ( 236) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3230) ! | pz  py  pz  s    ( 237) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3231) ! | pz  py  pz  px   ( 238) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3232) ! | pz  py  pz  py   ( 239) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_y) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3233) ! | pz  py  pz  pz   ( 240) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = 0.5*inv_ay*sypb*(I_A_y(n+1)+I_A_y(Abs(n-1)))
        real_Bn  = I_B_y(n)
        imag_An  = 0.5*cypb*inv_ay*(-I_A_y(n+1)+I_A_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( imag_An * real_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = 0.5*inv_az*szpa*(I_A_z(n+1)+I_A_z(Abs(n-1)))
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        imag_An  = 0.5*czpa*inv_az*(-I_A_z(n+1)+I_A_z(Abs(n-1)))
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3300) ! | pz  pz  s   s    ( 241) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = I_B_z(n)
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3301) ! | pz  pz  s   px   ( 242) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = I_B_z(n)
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3302) ! | pz  pz  s   py   ( 243) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = I_B_z(n)
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3303) ! | pz  pz  s   pz   ( 244) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3310) ! | pz  pz  px  s    ( 245) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = I_B_z(n)
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3311) ! | pz  pz  px  px   ( 246) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = -1.0*inv_ax2*(0.5*c2xqcd*I_B_x(n)+0.25*c2xqcd*I_B_x(n+2)+0.25*c2xqcd*I_B_x(Abs(n-2))-1.0*cxqc*cxqd*I_B_x(n))
        if (dabs(BBx) < 1.d-300) then 
          imag_Bn = 0.d0 
        else   
          imag_Bn  = -0.5*inv_ax2*n*s2xqcd*(I_B_x(n+1)+I_B_x(Abs(n-1)))/BBx
        end if 
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = I_B_z(n)
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3312) ! | pz  pz  px  py   ( 247) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = I_B_z(n)
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3313) ! | pz  pz  px  pz   ( 248) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqc*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqc*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3320) ! | pz  pz  py  s    ( 249) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = I_B_z(n)
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3321) ! | pz  pz  py  px   ( 250) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = I_B_z(n)
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3322) ! | pz  pz  py  py   ( 251) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = -1.0*inv_ay2*(0.5*c2yqcd*I_B_y(n)+0.25*c2yqcd*I_B_y(n+2)+0.25*c2yqcd*I_B_y(Abs(n-2))-1.0*cyqc*cyqd*I_B_y(n))
        if (dabs(BBy) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_ay2*n*s2yqcd*(I_B_y(n+1)+I_B_y(Abs(n-1)))/BBy
        end if 
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = I_B_z(n)
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        F_z(n)   = ( real_An * real_Bn ) * dcos(n * phi_z) - ( imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3323) ! | pz  pz  py  pz   ( 252) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqc*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqc*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = 0.5*inv_az*szqd*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        imag_Bn  = -0.5*czqd*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3330) ! | pz  pz  pz  s    ( 253) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3331) ! | pz  pz  pz  px   ( 254) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = 0.5*inv_ax*sxqd*(I_B_x(n+1)+I_B_x(Abs(n-1)))
        imag_Bn  = -0.5*cxqd*inv_ax*(-I_B_x(n+1)+I_B_x(Abs(n-1)))
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x) - ( real_An * imag_Bn ) * dsin(n * phi_x) 
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3332) ! | pz  pz  pz  py   ( 255) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = 0.5*inv_ay*syqd*(I_B_y(n+1)+I_B_y(Abs(n-1)))
        imag_Bn  = -0.5*cyqd*inv_ay*(-I_B_y(n+1)+I_B_y(Abs(n-1)))
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y) - ( real_An * imag_Bn ) * dsin(n * phi_y) 
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = 0.5*inv_az*szqc*(I_B_z(n+1)+I_B_z(Abs(n-1)))
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        imag_Bn  = -0.5*czqc*inv_az*(-I_B_z(n+1)+I_B_z(Abs(n-1)))
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      case (3333) ! | pz  pz  pz  pz   ( 256) 
      
      do n       = 0 , maxval(Nmax_x)
        real_An  = I_A_x(n)
        real_Bn  = I_B_x(n)
        F_x(n)   = ( real_An * real_Bn ) * dcos(n * phi_x)
      end do

      do n       = 0 , maxval(Nmax_y)
        real_An  = I_A_y(n)
        real_Bn  = I_B_y(n)
        F_y(n)   = ( real_An * real_Bn ) * dcos(n * phi_y)
      end do

      do n       = 0 , maxval(Nmax_z)
        real_An  = -1.0*inv_az2*(0.5*c2zpab*I_A_z(n)+0.25*c2zpab*I_A_z(n+2)+0.25*c2zpab*I_A_z(Abs(n-2))-1.0*czpa*czpb*I_A_z(n))
        real_Bn  = -1.0*inv_az2*(0.5*c2zqcd*I_B_z(n)+0.25*c2zqcd*I_B_z(n+2)+0.25*c2zqcd*I_B_z(Abs(n-2))-1.0*czqc*czqd*I_B_z(n))
        if (dabs(AAz) < 1.d-300) then 
          imag_An = 0.d0 
        else
          imag_An  = 0.5*inv_az2*n*s2zpab*(I_A_z(n+1)+I_A_z(Abs(n-1)))/AAz
        end if 
        if (dabs(BBz) < 1.d-300) then 
          imag_Bn = 0.d0 
        else
          imag_Bn  = -0.5*inv_az2*n*s2zqcd*(I_B_z(n+1)+I_B_z(Abs(n-1)))/BBz
        end if 
        F_z(n)   = ( real_An * real_Bn - imag_An * imag_Bn) * dcos(n * phi_z) - ( real_An * imag_Bn + imag_An * real_Bn ) * dsin(n * phi_z) 
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !


      case default
        F_x(n) = 0.d0
        F_y(n) = 0.d0
        F_z(n) = 0.d0  
          
      end select


      result = 0.0d0

      do n = 1, nquad_point
        result = result + w_val(n) * S(n)
      end do


      contains 

      double precision function S(i_quad) result(sum)

      use bessel_functions
      use table_lookup_module
      use precomputed_bessel
      
      ! --------------------------------------------------------------- !

      implicit none

      integer, intent(in)                  :: i_quad
      double precision                     :: sum1, sum2, sum3

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      sum1 = 2.d0 * dot_product(F_x(1:Nmax_x(i_quad)), I_C_table_x(1:Nmax_x(i_quad), i_quad))
      sum1 = sum1 + F_x(0) * I_C_table_x(0, i_quad)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      sum2 = 2.d0 * dot_product(F_y(1:Nmax_y(i_quad)), I_C_table_y(1:Nmax_y(i_quad), i_quad))
      sum2 = sum2 + F_y(0) * I_C_table_y(0, i_quad)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      sum3 = 2.d0 * dot_product(F_z(1:Nmax_z(i_quad)), I_C_table_z(1:Nmax_z(i_quad), i_quad))
      sum3 = sum3 + F_z(0) * I_C_table_z(0, i_quad)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !

      sum = sum1 * sum2 * sum3

      end function S


end subroutine integrate_ERI_3D_fast