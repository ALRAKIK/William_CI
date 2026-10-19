subroutine ERI_integral_toroidal_2D(number_of_atoms,geometry,atoms)

      use files
      use omp_lib
      use files
      use torus_init
      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      integer                        :: i , j , k , l , num_int , num_total_int
      integer                        :: number_of_atoms
      type(atom)                     :: atoms(number_of_atoms)
      type(ERI_function),allocatable :: ERI  (:)

      double precision               :: geometry(number_of_atoms,3)
      double precision,allocatable   ::       two_electron(:,:,:,:)
      double precision,allocatable   ::            two_eri(:,:,:,:)
      double precision               :: value
      double precision               :: start_time, end_time
      integer                        :: number_of_functions
      integer                        :: number_of_functions_per_unitcell 
      integer                        :: index_sym

      integer                        :: days, hours, minutes, seconds , t 


      !-----------------------------------------------------------------!

      call omp_set_dynamic(.false.)
      call omp_set_num_threads(omp_get_max_threads())

      number_of_functions_per_unitcell = 0 
      do i = 1 , number_of_atom_in_unitcell
        number_of_functions_per_unitcell = number_of_functions_per_unitcell + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      number_of_functions = 0 
      do i = 1 , number_of_atoms
        number_of_functions = number_of_functions + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do

      allocate(ERI(number_of_functions))
      allocate(two_electron(number_of_functions,number_of_functions,number_of_functions,number_of_functions))
      allocate(two_eri(number_of_functions,number_of_functions,number_of_functions,number_of_functions))

      call classification(number_of_atoms,number_of_functions,geometry,atoms,ERI)

      index_sym   = 0 
    
      do i = 1 , number_of_atoms/2 + 1 
        index_sym = index_sym + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do 

      ! 2-fold symmetry implementation (k,l permutation only)

      !$omp parallel
      if (omp_get_thread_num() == 0) then
        print *, "Running with", omp_get_num_threads(), "threads"
      endif
      !$omp end parallel

      start_time = omp_get_wtime()


      num_int = 0
      num_total_int = number_of_functions_per_unitcell * (number_of_functions) * (number_of_functions) * (number_of_functions+1) / 2 

      !$omp parallel do collapse(3) private(j, k, l, value) shared(two_electron, ERI,num_int) schedule(dynamic,16)

            do i = 1, number_of_functions_per_unitcell
              do j = 1, number_of_functions
                  do k = 1, number_of_functions
                      do l = k, number_of_functions
                          call ERI_integral_4_function_toroidal_2D(ERI(i),ERI(j),ERI(k),ERI(l), value)
                          
                          !$omp atomic
                          num_int = num_int + 1
                          two_electron(i,j,k,l) = value
                          two_electron(i,j,l,k) = value

                      end do
                  end do
              end do
            end do

      !$omp end parallel do

          

      end_time = omp_get_wtime()

      write(outfile,"(a)") "Translation symmetry applied to integrals"
      write(outfile,"(a)") "" 

      t = int(end_time - start_time)
      days= (t/86400)
      hours=mod(t,86400)/3600
      minutes=mod(mod(t,86400),3600)/60
      seconds=mod(mod(mod(t,86400),3600),60)

      write(outfile,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') '2 el-integrals calculation time = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec     "
      write(outfile,"(a)") "" 

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      call shift_integrals(two_electron,two_eri,number_of_functions,number_of_functions_per_unitcell)

      !open(1,file="./tmp/ERI.dat")
      open(1,file=trim(tmp_file_name)//"/ERI.dat")
        do i = 1, number_of_functions
          do j = 1 , number_of_functions
            do k = 1 , number_of_functions
              do l = 1 , number_of_functions
                if (abs(two_eri(i,j,k,l)) > 1e-8 ) write(1,*) i , j , k , l , two_eri(i,j,k,l)  
              end do 
            end do 
          end do 
        end do 
      close(1)

      deallocate(ERI)
      deallocate(two_electron)
      deallocate(two_eri)

end subroutine ERI_integral_toroidal_2D
!-----------------------------------------------------------------!