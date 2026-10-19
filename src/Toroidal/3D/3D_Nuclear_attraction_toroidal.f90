subroutine nuclear_attraction_matrix_toroidal_3D(number_of_atoms,nBas,geometry,atoms,AO,NA)

      use files
      use atom_basis
      use torus_init
      use classification_ERI
      use omp_lib
      use keywords
      use unitcell_module

      implicit none 

      !-----------------------------------------------------------------!

      integer                       :: i , j
      integer                       :: number_of_atoms
      integer                       :: nBas

      type(atom)                    :: atoms(number_of_atoms)

      type(ERI_function)            :: AO (nBas)
      type(ERI_function)            :: AO1 , AO2

      
      double precision              :: r1(3) , r2(3)
      double precision  ,intent(in) :: geometry(number_of_atoms,3)

      integer                       :: pattern_id, encode_orbital_pattern_AO

      double precision              :: temp_integral

      ! output !

      double precision,intent(out)  :: NA(nBas,nBas)

      !-----------------------------------------------------------------!
      ! Parallelization variables
      integer                        :: total_ij_pairs, ij_index
      integer, allocatable           :: i_index(:), j_index(:)
      integer                        :: num_threads, optimal_chunk_size
      !-----------------------------------------------------------------!

      !-----------------------------------------------------------------!

      NA(:,:) = 0.d0

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
      total_ij_pairs = nfuc * nBas

      allocate(i_index(total_ij_pairs), j_index(total_ij_pairs))

      ij_index = 0
      do i = 1, nfuc
       do j = 1, nBas
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

          call nuclear_attraction_integral_toroidal_3D(pattern_id,number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,temp_integral)

          NA(i,j) = temp_integral 
          
        end do 
      
      !$omp end parallel do
      
      deallocate(i_index, j_index)

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      do i = 1 , nfuc
       do j = 1 , nBas
         if (abs(NA(i,j)) < 1e-14) NA(i,j) = 0.d0 
       end do 
      end do

      call fill_with_symmetry(nBas, NA)

      if (c_NA) then 
        write(outfile,'(a)') ""
        write(outfile,'(a)') "------------------"
        write(outfile,'(a)') "The nuclear attraction matrix"
        write(outfile,'(a)') "------------------"
        write(outfile,'(a)') ""
        do i = 1 , nBas - 1 
          do j = i , nBas
            write(outfile,'(i3,i3, 6x,1000(f16.12,2x))') i , j, NA(i,j)
          end do
        end do 
      end if 

end subroutine nuclear_attraction_matrix_toroidal_3D


subroutine nuclear_attraction_integral_toroidal_3D(pattern_id,n_atoms,geometry,atoms,r1,r2,AO1,AO2,NA_prime)

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
      double precision                 :: ay2 , inv_ay , inv_ay2 
      double precision                 :: az2 , inv_az , inv_az2 
      double precision                 :: px  , py     , pz 
      double precision                 :: A_x , A_y    , A_z
      double precision                 :: xp  , yp     , zp 
      double precision                 :: sda_x   , sdb_x
      double precision                 :: cda_x   , cdb_x
      double precision                 :: sda_y   , sdb_y
      double precision                 :: cda_y   , cdb_y
      double precision                 :: sda_z   , sdb_z
      double precision                 :: cda_z   , cdb_z

      !-----------------------------------------------------------------!

      double precision                 :: albe

      !-----------------------------------------------------------------!

      !  Atoms ! 

      double precision                 :: xc , yc , zc
      double precision                 :: xpc ,  ypc,  zpc
      double precision                 :: ypc2, zpc2
      double precision                 :: kc
      double precision                 :: dx , xD
      double precision                 :: dy , yD
      double precision                 :: dz , zD

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
      ay2 = ay * ay
      az2 = az * az

      inv_ax  = 1.d0 / ax 
      inv_ay  = 1.d0 / ay 
      inv_az  = 1.d0 / az 

      inv_ax2 = inv_ax * inv_ax 
      inv_ay2 = inv_ay * inv_ay 
      inv_az2 = inv_az * inv_az 

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

          call bary_exponent_y         (alpha,beta,Y,py)
          call bary_center_toroidal_y  (alpha,beta,y1,y2,yp)

          call bary_exponent_z         (alpha,beta,Z,pz)
          call bary_center_toroidal_z  (alpha,beta,z1,z2,zp)

          kx     = 1.d0 
          ky     = 1.d0 
          kz     = 1.d0 

          albe   = alpha + beta

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
      double precision             :: k_na_x , k_na_y, k_na_z
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

      t2         = t  * t 
      
      ! Clifford ! 

      call bary_exponent_x(px,t2,xpc,dx)
      call bary_center_toroidal_x(px,t2,xp,xc,xd)

      call bary_exponent_y(py,t2,ypc,dy)
      call bary_center_toroidal_y(py,t2,yp,yc,yd)

      call bary_exponent_z(pz,t2,zpc,dz)
      call bary_center_toroidal_z(pz,t2,zp,zc,zd)

      dx = 2.d0 * dx  * inv_ax2
      dy = 2.d0 * dy  * inv_ay2
      dz = 2.d0 * dz  * inv_az2

      A_x = dx
      A_y = dy
      A_z = dz

      Ix_int = Lx * iv_scaled(0, A_x)
      Iy_int = Ly * iv_scaled(0, A_y)
      Iz_int = Lz * iv_scaled(0, A_z)

      k_na_x = dexp( - 2.d0 * (t2 + albe)  * inv_ax2 + dx )
      k_na_y = dexp( - 2.d0 * (t2 + albe)  * inv_ay2 + dy )
      k_na_z = dexp( - 2.d0 * (t2 + albe)  * inv_az2 + dz )

      sda_x    = dsin(ax*(xd-x1))
      cda_x    = dcos(ax*(xd-x1))

      sdb_x    = dsin(ax*(xd-x2))
      cdb_x    = dcos(ax*(xd-x2))

      sda_y    = dsin(ay*(yd-y1))
      cda_y    = dcos(ay*(yd-y1))

      sdb_y    = dsin(ay*(yd-y2))
      cdb_y    = dcos(ay*(yd-y2))

      sda_z    = dsin(az*(zd-z1))
      cda_z    = dcos(az*(zd-z1))

      sdb_z    = dsin(az*(zd-z2))
      cdb_z    = dcos(az*(zd-z2))

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
      
        Ix   = inv_ax * (Icliff(1,0,A_x)*cdb_x + Icliff(0,1,A_x)*sdb_x)
        Iy   = Iy_int
        Iz   = Iz_int
      
      case (02) ! | s     py    ( 3 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        Ix   = Ix_int
        Iy   = inv_ay * (Icliff(1,0,A_y)*cdb_y + Icliff(0,1,A_y)*sdb_y)
        Iz   = Iz_int
      
      case (03) ! | s     pz    ( 4 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        Ix   = Ix_int
        Iy   = Iy_int
        Iz   = inv_az * (Icliff(1,0,A_z)*cdb_z + Icliff(0,1,A_z)*sdb_z)
      
      case (10) ! | px    s     ( 5 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=0)
      
        Ix   = inv_ax * (Icliff(1,0,A_x)*cda_x + Icliff(0,1,A_x)*sda_x)
        Iy   = Iy_int
        Iz   = Iz_int
      
      case (11) ! | px    px    ( 6 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=1, m=0, n=0)
      
        Ix   = inv_ax2 * (Icliff(2,0,A_x)*cda_x*cdb_x + Icliff(1,1,A_x)*cda_x*sdb_x + Icliff(1,1,A_x)*cdb_x*sda_x + Icliff(0,2,A_x)*sda_x*sdb_x)
        Iy   = Iy_int
        Iz   = Iz_int
      
      case (12) ! | px    py    ( 7 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        Ix   = inv_ax * (Icliff(1,0,A_x)*cda_x + Icliff(0,1,A_x)*sda_x)
        Iy   = inv_ay * (Icliff(1,0,A_y)*cdb_y + Icliff(0,1,A_y)*sdb_y)
        Iz   = Iz_int
      
      case (13) ! | px    pz    ( 8 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        Ix   = inv_ax * (Icliff(1,0,A_x)*cda_x + Icliff(0,1,A_x)*sda_x)
        Iy   = Iy_int
        Iz   = inv_az * (Icliff(1,0,A_z)*cdb_z + Icliff(0,1,A_z)*sdb_z)
      
      case (20) ! | py    s     ( 9 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=0)
      
        Ix   = Ix_int
        Iy   = inv_ay * (Icliff(1,0,A_y)*cda_y + Icliff(0,1,A_y)*sda_y)
        Iz   = Iz_int
      
      case (21) ! | py    px    ( 10 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=1, m=0, n=0)
      
        Ix   = inv_ax * (Icliff(1,0,A_x)*cdb_x + Icliff(0,1,A_x)*sdb_x)
        Iy   = inv_ay * (Icliff(1,0,A_y)*cda_y + Icliff(0,1,A_y)*sda_y)
        Iz   = Iz_int
      
      case (22) ! | py    py    ( 11 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=1, n=0)
      
        Ix   = Ix_int
        Iy   = inv_ay2 * (Icliff(2,0,A_y)*cda_y*cdb_y + Icliff(1,1,A_y)*cda_y*sdb_y + Icliff(1,1,A_y)*cdb_y*sda_y + Icliff(0,2,A_y)*sda_y*sdb_y)
        Iz   = Iz_int
      
      case (23) ! | py    pz    ( 12 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=1)
      
        Ix   = Ix_int
        Iy   = inv_ay * (Icliff(1,0,A_y)*cda_y + Icliff(0,1,A_y)*sda_y)
        Iz   = inv_az * (Icliff(1,0,A_z)*cdb_z + Icliff(0,1,A_z)*sdb_z)
      
      case (30) ! | pz    s     ( 13 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=0)
      
        Ix   = Ix_int
        Iy   = Iy_int
        Iz   = inv_az * (Icliff(1,0,A_z)*cda_z + Icliff(0,1,A_z)*sda_z)
      
      case (31) ! | pz    px    ( 14 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=1, m=0, n=0)
      
        Ix   = inv_ax * (Icliff(1,0,A_x)*cdb_x + Icliff(0,1,A_x)*sdb_x)
        Iy   = Iy_int
        Iz   = inv_az * (Icliff(1,0,A_z)*cda_z + Icliff(0,1,A_z)*sda_z)
      
      case (32) ! | pz    py    ( 15 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=1, n=0)
      
        Ix   = Ix_int
        Iy   = inv_ay * (Icliff(1,0,A_y)*cdb_y + Icliff(0,1,A_y)*sdb_y)
        Iz   = inv_az * (Icliff(1,0,A_z)*cda_z + Icliff(0,1,A_z)*sda_z)
      
      case (33) ! | pz    pz    ( 16 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=1)
      
        Ix   = Ix_int
        Iy   = Iy_int
        Iz   = inv_az2 * (Icliff(2,0,A_z)*cda_z*cdb_z + Icliff(1,1,A_z)*cda_z*sdb_z + Icliff(1,1,A_z)*cdb_z*sda_z + Icliff(0,2,A_z)*sda_z*sdb_z)

      case default

        I_t = 0.0d0

      end select

        I_t  = Ix * Iy * Iz * k_na_x * k_na_y * k_na_z

      !-----------------------------------------------------------------!

      end function f_decay

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_toroidal_3D