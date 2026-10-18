subroutine nuclear_attraction_matrix_toroidal_1D(number_of_atoms,number_of_functions,geometry,atoms,AO,NA)

      use files
      use atom_basis
      use torus_init
      use classification_ERI
      use omp_lib
      use keywords

      implicit none 

      !-----------------------------------------------------------------!

      integer                       :: i , j
      integer                       :: index_atom1 , index_sym
      integer                       :: index_unitcell
      integer                       :: number_of_atoms
      integer                       :: number_of_functions

      type(atom)                    :: atoms(number_of_atoms)

      type(ERI_function)            :: AO (number_of_functions)
      type(ERI_function)            :: AO1 , AO2

      
      double precision              :: r1(3) , r2(3)
      double precision  ,intent(in) :: geometry(number_of_atoms,3)

      integer                       :: pattern_id, encode_orbital_pattern_AO

      double precision              :: temp_integral

      ! output !

      double precision,intent(out)  :: NA(number_of_functions,number_of_functions)

      !-----------------------------------------------------------------!
      ! Parallelization variables
      integer                        :: total_ij_pairs, ij_index
      integer, allocatable           :: i_index(:), j_index(:)
      integer                        :: num_threads, optimal_chunk_size
      !-----------------------------------------------------------------!

      !-----------------------------------------------------------------!

      NA(:,:) = 0.d0

      index_atom1 = atoms(1)%num_s_function + 3*atoms(1)%num_p_function
      
      index_unitcell = 0

      do i = 1 , number_of_atom_in_unitcell
        index_unitcell = index_unitcell  + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do

      index_sym   = 0

      do i = 1 , number_of_atoms/2 + 1
        index_sym = index_sym + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do

      ! --------------------------------------------------------------- !

      !-----------------------------------------------------------------!
      ! Setup OpenMP
      !-----------------------------------------------------------------!
      call omp_set_dynamic(.false.)
      call omp_set_num_threads(omp_get_max_threads())

      !$omp parallel
      if (omp_get_thread_num() == 0) then
       print *, "Nuclear Attraction: Running with", omp_get_num_threads(), "threads"
      endif
      !$omp end parallel
      
      !$omp parallel
       !$omp single
         num_threads = omp_get_num_threads()
       !$omp end single
      !$omp end parallel
      
      !Calculate optimal chunk size based on available threads
      if (num_threads <= 16) then
       optimal_chunk_size = 16
      else if (num_threads <= 64) then
       optimal_chunk_size = 8
      else 
       optimal_chunk_size = 1
      end if

      !-----------------------------------------------------------------!
      ! Precompute all i-j pairs
      !-----------------------------------------------------------------!
      total_ij_pairs = index_unitcell * number_of_functions
      allocate(i_index(total_ij_pairs), j_index(total_ij_pairs))

      ij_index = 0
      do i = 1, index_unitcell
       do j = 1, number_of_functions
         ij_index = ij_index + 1
         i_index(ij_index) = i
         j_index(ij_index) = j
       end do
      end do

      !-----------------------------------------------------------------!
      !Parallel computation
      !-----------------------------------------------------------------!

      !$omp parallel do private(ij_index,i,j,AO1,AO2,r1,r2, pattern_id, temp_integral) &
      !$omp shared(NA, AO, i_index, j_index, number_of_atoms, geometry, atoms, total_ij_pairs) &
      !$omp schedule(dynamic,optimal_chunk_size)


      do ij_index = 1, total_ij_pairs
        i = i_index(ij_index)
        j = j_index(ij_index)

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z

          pattern_id = encode_orbital_pattern_AO(AO1%orbital, AO2%orbital)

          call nuclear_attraction_integral_toroidal_1D(pattern_id,number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,temp_integral)
          NA(i,j) = temp_integral 
          
        end do 
      
      !$omp end parallel do
      
      deallocate(i_index, j_index)

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!


      do i = 1 , index_unitcell
       do j = 1 , number_of_functions
         if (abs(NA(i,j)) < 1e-15) NA(i,j) = 0.d0 
       end do 
      end do 

      do i = index_unitcell + 1   , number_of_functions
       do j = index_unitcell + 1 , number_of_functions
         NA(i,j) = NA(i-index_unitcell,j-index_unitcell)
       end do 
      end do 
 
      do i = 1 , number_of_functions - 1 
       do j = i , number_of_functions
         NA(j,i) = NA(i,j)
       end do 
      end do

      if (c_NA) then 
        write(outfile,'(a)') ""
        write(outfile,'(a)') "------------------"
        write(outfile,'(a)') "The nuclear attraction matrix"
        write(outfile,'(a)') "------------------"
        write(outfile,'(a)') ""
        do i = 1 , number_of_functions - 1 
          do j = i , number_of_functions
            write(outfile,'(i3,i3, 6x,1000(f16.12,2x))') i , j, NA(i,j)
          end do
        end do 
      end if 

end subroutine nuclear_attraction_matrix_toroidal_1D

      !-----------------------------------------------------------------!

subroutine nuclear_attraction_integral_toroidal_1D(pattern_id,n_atoms,geometry,atoms,r1,r2,AO1,AO2,NA_prime)

      use quadpack, only : dqagi
      use torus_init
      use classification_ERI
      use bessel_functions
      use atom_basis

      implicit none 

      !-----------------------------------------------------------------!

      
      type(ERI_function),intent(in)    :: AO1 , AO2
      integer           ,intent(in)    :: n_atoms
      type(atom)        ,intent(in)    :: atoms(n_atoms)
      double precision  ,intent(in)    :: r1(3) , r2(3)
      double precision  ,intent(in)    :: geometry(n_atoms,3)
      double precision  ,intent(out)   :: NA_prime

      integer                          :: i , j , k
      integer                          :: charge_atom
      double precision,parameter       :: pi = 3.14159265358979323846D00
      double precision                 :: alpha  , beta
      double precision                 :: c1     , c2     , const
      double precision                 :: x1     , x2     , X  
      double precision                 :: y1     , y2     , Y
      double precision                 :: z1     , z2     , Z  
      double precision                 :: kx     , ky     , kz 
      double precision                 :: Ix     , Iy     , Iz  
      double precision                 :: Ix_int , Iy_int , Iz_int  
      double precision                 :: NA
       
      !-----------------------------------------------------------------!

      ! Clifford ! 

      double precision                 :: ax2 , inv_ax , inv_ax2 
      double precision                 :: px
      double precision                 :: A
      double precision                 :: xp
      double precision                 :: sda   , sdb
      double precision                 :: cda   , cdb

      !-----------------------------------------------------------------!

      !   Real   ! 

      double precision                 :: albe, inv_albe , mu 
      double precision                 :: yp   , zp
      double precision                 :: ypa  , ypb 
      double precision                 :: zpa  , zpb

      !-----------------------------------------------------------------!

      !  Atoms ! 

      double precision                 :: xc , yc , zc
      double precision                 :: xpc ,  ypc,  zpc
      double precision                 :: ypc2, zpc2
      double precision                 :: kc
      double precision                 :: dx , xD

      !-----------------------------------------------------------------!
      
      ! Local variables for numerical integration
      
      double precision,parameter    :: epsabs = 1.d-10 , epsrel = 1.d-8
      integer,parameter             :: inf = 1 
      double precision,parameter    :: bound = 0.0d0
      integer, parameter            :: limit = 100
      integer, parameter            :: lenw = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)
      double precision              :: result

      !-----------------------------------------------------------------!

      integer                       :: pattern_id

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
 
      NA       = 0.d0
      NA_prime = 0.d0 

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

          !kx     = dexp(-2.d0 * ( alpha + beta ) / ax2)
          kx     = 1.d0 

          ! care here the did not add the px in the exponent as 
          ! we dont have to calculat it, but we calculate dx 
          ! inside the integral 

          ! Note all the standard integrals should go 
          ! inside the integration for now 

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

          ypa = yp - y1
          ypb = yp - y2

          zpa = zp - z1 
          zpb = zp - z2

          do k = 1 , n_atoms

            xc  = geometry(k,1) 
            yc  = geometry(k,2)
            zc  = geometry(k,3)
            kc  = 2.d0 / sqrt(pi)
            xpc = xp - xc 
            ypc = yp - yc 
            zpc = zp - zc 

            ypc2 = ypc * ypc 
            zpc2 = zpc * zpc 
              
            charge_atom = (-1)*atoms(k)%charge

            call dqagi(f_decay, bound, inf, epsabs, epsrel, result, abserr, neval, ier, Limit, Lenw, Last, Iwork, Work)

            NA_prime = NA_prime + charge_atom * result * kx * ky * kz * kc * const 

            if (ier /= 0) then
              write(*,'(A,//)')   'Error code from ss NA :'
              write(*,'(A,I8,A)') 'Error code = ', ier
            end if
            
          end do 

        end do 
      end do

      contains

      function f_decay(t) result(I_t)

      use functions

      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !
      double precision, intent(in) :: t
      double precision             :: I_t
      double precision             :: t2
      double precision             :: albept , inv_albept
      double precision             :: eta_t
      double precision             :: k_na_x , k_na_y, k_na_z
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

      t2         = t  * t 
      albept     = albe + t2
      inv_albept = 1.d0 / albept
      eta_t      = t2 * inv_albept
      
      ! Clifford ! 

      call bary_exponent_x(px,t2,xpc,dx)
      call bary_center_toroidal_x(px,t2,xp,xc,xd)

      dx = 2.d0 * dx  * inv_ax2

      A = dx 

      Ix_int = Lx * iv_scaled(0, A)

      k_na_x = dexp( - 2.d0 * (t2 + albe)  * inv_ax2 + dx )

      sda    = dsin(ax*(xd-x1))
      cda    = dcos(ax*(xd-x1))

      sdb    = dsin(ax*(xd-x2))
      cdb    = dcos(ax*(xd-x2))

      ! Real ! 

      k_na_y = dexp(- albe * eta_t * ( yPC * yPC ) )
      k_na_z = dexp(- albe * eta_t * ( zPC * zPC ) )

      Iy_int = dsqrt(pi * inv_albept)
      Iz_int = dsqrt(pi * inv_albept)

      select case(pattern_id)

      case (00) ! | s     s     ( 1 )
              
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=0, n=0)
              
        Ix   = Ix_int
        Iy   = Iy_int
        Iz   = Iz_int
              
      case (01) ! | s     px    ( 2 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=1, m=0, n=0)
      
        Ix   = inv_ax * (Icliff(1,0,A)*cdb + Icliff(0,1,A)*sdb)
        Iy   = Iy_int
        Iz   = Iz_int
      
      case (02) ! | s     py    ( 3 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        Ix   = Ix_int
        Iy   = -Iy_int*eta_t*ypc + Iy_int*ypb + Ireal(1,albept)
        Iz   = Iz_int
      
      case (03) ! | s     pz    ( 4 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        Ix   = Ix_int
        Iy   = Iy_int
        Iz   = -Iz_int*eta_t*zpc + Iz_int*zpb + Ireal(1,albept)
      
      case (10) ! | px    s     ( 5 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=0)
      
        Ix   = inv_ax * (Icliff(1,0,A)*cda + Icliff(0,1,A)*sda)
        Iy   = Iy_int
        Iz   = Iz_int
      
      case (11) ! | px    px    ( 6 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=1, m=0, n=0)
      
        Ix   = inv_ax2 * (Icliff(2,0,A)*cda*cdb + Icliff(1,1,A)*cda*sdb + Icliff(1,1,A)*cdb*sda + Icliff(0,2,A)*sda*sdb)
        Iy   = Iy_int
        Iz   = Iz_int
      
      case (12) ! | px    py    ( 7 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        Ix   = inv_ax * (Icliff(1,0,A)*cda + Icliff(0,1,A)*sda)
        Iy   = -Iy_int*eta_t*ypc + Iy_int*ypb + Ireal(1,albept)
        Iz   = Iz_int
      
      case (13) ! | px    pz    ( 8 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        Ix   = inv_ax * (Icliff(1,0,A)*cda + Icliff(0,1,A)*sda)
        Iy   = Iy_int
        Iz   = -Iz_int*eta_t*zpc + Iz_int*zpb + Ireal(1,albept)
      
      case (20) ! | py    s     ( 9 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=0)
      
        Ix   = Ix_int
        Iy   = -Iy_int*eta_t*ypc + Iy_int*ypa + Ireal(1,albept)
        Iz   = Iz_int
      
      case (21) ! | py    px    ( 10 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=1, m=0, n=0)
      
        Ix   = inv_ax * (Icliff(1,0,A)*cdb + Icliff(0,1,A)*sdb)
        Iy   = -Iy_int*eta_t*ypc + Iy_int*ypa + Ireal(1,albept)
        Iz   = Iz_int
      
      case (22) ! | py    py    ( 11 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=1, n=0)
      
        Ix   = Ix_int
        Iy   = Iy_int*eta_t**2*ypc**2 - Iy_int*eta_t*ypa*ypc - Iy_int*eta_t*ypb*ypc - Ireal(1,albept)*2*eta_t*ypc + Iy_int*ypa*ypb + Ireal(1,albept)*ypa + Ireal(1,albept)*ypb + Ireal(2,albept)
        Iz   = Iz_int
      
      case (23) ! | py    pz    ( 12 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=1)
      
        Ix   = Ix_int
        Iy   = -Iy_int*eta_t*ypc + Iy_int*ypa + Ireal(1,albept)
        Iz   = -Iz_int*eta_t*zpc + Iz_int*zpb + Ireal(1,albept)
      
      case (30) ! | pz    s     ( 13 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=0)
      
        Ix   = Ix_int
        Iy   = Iy_int
        Iz   = -Iz_int*eta_t*zpc + Iz_int*zpa + Ireal(1,albept)
      
      case (31) ! | pz    px    ( 14 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=1, m=0, n=0)
      
        Ix   = inv_ax * (Icliff(1,0,A)*cdb + Icliff(0,1,A)*sdb)
        Iy   = Iy_int
        Iz   = -Iz_int*eta_t*zpc + Iz_int*zpa + Ireal(1,albept)
      
      case (32) ! | pz    py    ( 15 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=1, n=0)
      
        Ix   = Ix_int
        Iy   = -Iy_int*eta_t*ypc + Iy_int*ypb + Ireal(1,albept)
        Iz   = -Iz_int*eta_t*zpc + Iz_int*zpa + Ireal(1,albept)
      
      case (33) ! | pz    pz    ( 16 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=1)
      
        Ix   = Ix_int
        Iy   = Iy_int
        Iz   = Iz_int*eta_t**2*zpc**2 - Iz_int*eta_t*zpa*zpc - Iz_int*eta_t*zpb*zpc - Ireal(1,albept)*2*eta_t*zpc + Iz_int*zpa*zpb + Ireal(1,albept)*zpa + Ireal(1,albept)*zpb + Ireal(2,albept)

      case default
        I_t = 0.0d0

      end select

        I_t  = Ix * Iy * Iz * k_na_x * k_na_y * k_na_z

      !-----------------------------------------------------------------!

      end function f_decay

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_toroidal_1D