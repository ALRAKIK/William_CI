subroutine nuclear_attraction_matrix_toroidal_SAO_1D(number_of_atoms,nBas,geometry,atoms,AO,NA)

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

          call nuclear_attraction_integral_toroidal_1D(pattern_id,number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,temp_integral)

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

end subroutine nuclear_attraction_matrix_toroidal_SAO_1D