subroutine overlap_matrix_toroidal_1D_sym(nBas,AO,overlap,phase)

      use files
      use torus_init
      use classification_ERI
      use keywords
      use unitcell_module

      implicit none

      integer           ,intent(in)  :: nBas
      type(ERI_function),intent(in)  :: AO(nBas)
      integer           ,intent(in)  :: phase(nfuc)


      double precision  ,intent(out) :: overlap(nBas,nBas)


      type(ERI_function)            :: AO1, AO2
      double precision              :: r1(3), r2(3)

      integer                       :: i, j, jj
      integer                       :: pattern_id, encode_orbital_pattern_AO
      integer                       :: cx, cy, cz, icx, icy, icz, orb
      integer                       :: cell_index, inv_cell_index

      

      !-----------------------------------------------------------------!

      overlap(:,:) = 0.d0

      !=================================================================!
      !  Pass 1: compute half of the cells (cell_index <= inv)          !
      !=================================================================!

      do cz = 0, nz - 1
         icz = modulo(-cz, nz)
         do cy = 0, ny - 1
            icy = modulo(-cy, ny)
            do cx = 0, nx - 1
               icx = modulo(-cx, nx)
            
               cell_index     =  cx + nx*cy     + nx*ny*cz
               inv_cell_index = icx + nx*icy    + nx*ny*icz
            
               if (cell_index > inv_cell_index) cycle
            
               do i = 1, nfuc
                  AO1 = AO(i)
                  do orb = 1, nfuc
                     j    = cell_index*nfuc + orb
                     AO2  = AO(j)
                  
                     r1(1)=AO1%x ; r2(1)=AO2%x
                     r1(2)=AO1%y ; r2(2)=AO2%y
                     r1(3)=AO1%z ; r2(3)=AO2%z
                  
                     pattern_id = encode_orbital_pattern_AO(AO1%orbital, AO2%orbital)
                     call overlap_integral_toroidal_1D(pattern_id, r1, r2, AO1, AO2, overlap(i,j))
                  end do
               end do
            end do
         end do
      end do

      !=================================================================!
      !  Pass 2: fill the other half via inversion                      !
      !                                                                 !
      !  O(i,j) = phase(i) * phase(orb) * O(i, inv(j))                  !
      !=================================================================!

      do cz = 0, nz - 1
         icz = modulo(-cz, nz)
         do cy = 0, ny - 1
            icy = modulo(-cy, ny)
            do cx = 0, nx - 1
               icx = modulo(-cx, nx)
            
               cell_index     =  cx + nx*cy  + nx*ny*cz
               inv_cell_index = icx + nx*icy + nx*ny*icz
            
               if (cell_index <= inv_cell_index) cycle
            
               do i = 1, nfuc
                  do orb = 1, nfuc
                     j  = cell_index*nfuc     + orb
                     jj = inv_cell_index*nfuc + inv_perm(orb)
                     overlap(i,j) = phase(i) * phase(orb) * overlap(inv_perm(i),jj)
                  end do
               end do
            end do
         end do
      end do

      !-----------------------------------------------------------------!
      !  Clean tiny values                                              !
      !-----------------------------------------------------------------!

      do i = 1, nfuc
         do j = 1, nBas
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

end subroutine overlap_matrix_toroidal_1D_sym

subroutine kinetic_matrix_toroidal_1D_sym(nBas,AO,kinetic,phase)

      use files
      use torus_init
      use classification_ERI
      use keywords
      use unitcell_module

      implicit none


      integer           ,intent(in)  :: nBas
      type(ERI_function),intent(in)  :: AO(nBas)
      integer           ,intent(in)  :: phase(nfuc)

      double precision  ,intent(out) :: kinetic(nBas,nBas)
      
      type(ERI_function)             :: AO1, AO2
      double precision               :: r1(3), r2(3)
      integer                        :: pattern_id, encode_orbital_pattern_AO
      integer                        :: i, j, jj

      integer                        :: cx, cy, cz, icx, icy, icz, orb
      integer                        :: cell_index, inv_cell_index

      

      !-----------------------------------------------------------------!

      kinetic(:,:) = 0.d0
   
      !=================================================================!
      !  Pass 1: compute half of the cells (cell_index <= inv)          !
      !=================================================================!

      do cz = 0, nz - 1
         icz = modulo(-cz, nz)
         do cy = 0, ny - 1
            icy = modulo(-cy, ny)
            do cx = 0, nx - 1
               icx = modulo(-cx, nx)
            
               cell_index     =  cx + nx*cy     + nx*ny*cz
               inv_cell_index = icx + nx*icy    + nx*ny*icz
            
               if (cell_index > inv_cell_index) cycle
            
               do i = 1, nfuc
                  AO1 = AO(i)
                  do orb = 1, nfuc
                     j    = cell_index*nfuc + orb

                     AO2  = AO(j)
                  
                     r1(1)=AO1%x ; r2(1)=AO2%x
                     r1(2)=AO1%y ; r2(2)=AO2%y
                     r1(3)=AO1%z ; r2(3)=AO2%z
                  
                     pattern_id = encode_orbital_pattern_AO(AO1%orbital, AO2%orbital)

                     call Kinetic_integral_toroidal_1D(pattern_id,r1,r2,AO1,AO2,kinetic(i,j))

                  end do
               end do
            end do
         end do
      end do

      !=================================================================!
      !  Pass 2: fill the other half via inversion                      !
      !                                                                 !
      !  O(i,j) = phase(i) * phase(orb) * O(i, inv(j))                  !
      !=================================================================!

      do cz = 0, nz - 1
         icz = modulo(-cz, nz)
         do cy = 0, ny - 1
            icy = modulo(-cy, ny)
            do cx = 0, nx - 1
               icx = modulo(-cx, nx)
            
               cell_index     =  cx + nx*cy  + nx*ny*cz
               inv_cell_index = icx + nx*icy + nx*ny*icz
            
               if (cell_index <= inv_cell_index) cycle
            
               do i = 1, nfuc
                  do orb = 1, nfuc
                     j  = cell_index*nfuc     + orb
                     jj = inv_cell_index*nfuc + inv_perm(orb)
                     kinetic(i,j) = phase(i) * phase(orb) * kinetic(inv_perm(i),jj)
                  end do
               end do
            end do
         end do
      end do

      !-----------------------------------------------------------------!
      !  Clean tiny values                                              !
      !-----------------------------------------------------------------!
      do i = 1, nfuc
         do j = 1, nBas
            if (abs(kinetic(i,j)) < 1e-14) kinetic(i,j) = 0.d0
         end do
      end do

      !-----------------------------------------------------------------!
      ! Clean tiny values + translational symmetry                      !
      !-----------------------------------------------------------------!
      
      do i = 1 , nfuc
        do j = 1 , nBas
          if (abs(kinetic(i,j)) < 1e-14) kinetic(i,j) = 0.d0 
        end do 
      end do 

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

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
   
end subroutine kinetic_matrix_toroidal_1D_sym


subroutine nuclear_attraction_matrix_toroidal_1D_sym(number_of_atoms,nBas,geometry,atoms,AO,NA,phase)

      use files
      use atom_basis
      use torus_init
      use classification_ERI
      use omp_lib
      use keywords
      use unitcell_module

      implicit none

      integer           ,intent(in)  :: nBas
      type(ERI_function),intent(in)  :: AO(nBas)
      integer           ,intent(in)  :: phase(nfuc)
      double precision,intent(in)    :: geometry(number_of_atoms,3)

      double precision  ,intent(out) :: NA(nBas,nBas)



      integer                       :: i, j, jj
      integer                       :: number_of_atoms

      type(atom)                    :: atoms(number_of_atoms)
      type(ERI_function)            :: AO1, AO2

      double precision              :: r1(3), r2(3)
      

      integer                       :: pattern_id, encode_orbital_pattern_AO
      double precision              :: temp_integral

      !--- parallelisation
      integer                       :: total_ij_pairs, ij_index
      integer, allocatable          :: i_index(:), j_index(:)
      integer                       :: num_threads, optimal_chunk_size

      !--- inversion symmetry
      
      integer, allocatable          :: cell_list(:)
      integer                       :: n_computed_cells, idx, cidx
      integer                       :: cx, cy, cz, icx, icy, icz, orb
      integer                       :: cell_index, inv_cell_index

      !-----------------------------------------------------------------!

      NA(:,:) = 0.d0

      !-----------------------------------------------------------------!
      ! Build list of "computed" cells: cell_index <= inv_cell_index    !
      !-----------------------------------------------------------------!

      allocate(cell_list(N_cell))
      n_computed_cells = 0
      do cz = 0, nz - 1
         icz = modulo(-cz, nz)
         do cy = 0, ny - 1
            icy = modulo(-cy, ny)
            do cx = 0, nx - 1
               icx = modulo(-cx, nx)
               cell_index     =  cx + nx*cy     + nx*ny*cz
               inv_cell_index = icx + nx*icy    + nx*ny*icz
               if (cell_index <= inv_cell_index) then
                  n_computed_cells = n_computed_cells + 1
                  cell_list(n_computed_cells) = cell_index
               end if
            end do
         end do
      end do

      !-----------------------------------------------------------------!
      ! OpenMP setup (unchanged)                                        !
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

      if (num_threads <= 16) then
         optimal_chunk_size = 16
      else if (num_threads <= 64) then
         optimal_chunk_size = 8
      else
         optimal_chunk_size = 1
      end if

      !-----------------------------------------------------------------!
      ! Precompute (i,j) pairs, j in a computed cell only               !
      ! total = nfuc * nfuc * n_computed_cells (~ half of old nfuc*nBas)!
      !-----------------------------------------------------------------!

      total_ij_pairs = nfuc * nfuc * n_computed_cells
      allocate(i_index(total_ij_pairs), j_index(total_ij_pairs))

      ij_index = 0
      do idx = 1, n_computed_cells
         cidx = cell_list(idx)
         do orb = 1, nfuc
            j = cidx*nfuc + orb
            do i = 1, nfuc
               ij_index = ij_index + 1
               i_index(ij_index) = i
               j_index(ij_index) = j
            end do
         end do
      end do

      !-----------------------------------------------------------------!
      ! Parallel computation of the "computed" half                     !
      !-----------------------------------------------------------------!

      !$omp parallel do private(ij_index,i,j,AO1,AO2,r1,r2,pattern_id,temp_integral) &
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

         call nuclear_attraction_integral_toroidal_1D(pattern_id,number_of_atoms,geometry,atoms, &
                                                      r1,r2,AO1,AO2,temp_integral)
         NA(i,j) = temp_integral
      end do
      !$omp end parallel do

      deallocate(i_index, j_index, cell_list)

      !-----------------------------------------------------------------!
      ! Fill the other half via inversion                               !
      !                                                                 !
      ! NA(i,j) = phase(i) * phase(orb) * NA(i, inv(j))                 !
      !-----------------------------------------------------------------!

      do cz = 0, nz - 1
         icz = modulo(-cz, nz)
         do cy = 0, ny - 1
            icy = modulo(-cy, ny)
            do cx = 0, nx - 1
               icx = modulo(-cx, nx)
               cell_index     =  cx + nx*cy  + nx*ny*cz
               inv_cell_index = icx + nx*icy + nx*ny*icz
               if (cell_index <= inv_cell_index) cycle
               do i = 1, nfuc
                  do orb = 1, nfuc
                     j  = cell_index*nfuc     + orb
                     jj = inv_cell_index*nfuc + inv_perm(orb)
                     NA(i,j) = phase(i) * phase(orb) * NA(inv_perm(i),jj)
                  end do
               end do
            end do
         end do
      end do

      !-----------------------------------------------------------------!
      ! Clean tiny values                                               !
      !-----------------------------------------------------------------!

      do i = 1, nfuc
         do j = 1, nBas
            if (abs(NA(i,j)) < 1e-14) NA(i,j) = 0.d0
         end do
      end do

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

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

end subroutine nuclear_attraction_matrix_toroidal_1D_sym

subroutine ERI_integral_toroidal_1D_sym(number_of_atoms,geometry,nBas,atoms,two_electron_integrals)

      use files
      use keywords
      use omp_lib
      use files
      use torus_init
      use atom_basis
      use classification_ERI
      use unitcell_module
      use precomputed_bessel

      implicit none 

      !-----------------------------------------------------------------!

      integer                        :: i, j, k, l
      integer                        :: number_of_atoms
      integer                        :: nBas
      type(atom)                     :: atoms(number_of_atoms)
      type(ERI_function),allocatable :: ERI(:)

      double precision               :: geometry(number_of_atoms,3)
      double precision,allocatable   :: two_electron(:,:,:,:)
      double precision               :: num_total_int
      double precision               :: num_total_int_without

      ! --------------------------------------------------------------- !
      
      double precision,intent(out)   :: two_electron_integrals(ERI_size,nBas,nBas,nBas)

      ! ----------------------    Time     ---------------------------- !
      integer                        :: days, hours, minutes, seconds
      integer                        :: t , time 
      double precision               :: start,end
      double precision               :: start_time, end_time
      !-----------------------------------------------------------------!
      integer                        :: total_ij_pairs, ij_index
      integer                        :: total_ij_pairs_without
      integer, allocatable           :: i_index(:), j_index(:)
      integer, allocatable           :: i_index_without(:), j_index_without(:)
      integer                        :: num_threads, optimal_chunk_size
      !-----------------------------------------------------------------!

      integer                        :: current_pct 
      integer                        :: last_pct = -1
      double precision               :: integrals_done = 0.d0
      double precision               :: val

      !-----------------------------------------------------------------!
      !                      Schwarz screening                          !
      double precision, allocatable  :: Q_schwarz(:,:)
      double precision               :: val_iijj
      double precision, parameter    :: schwarz_thresh = 1.d-14
      integer(8)                     :: n_screened, n_phase1_computed

      !-- inversion symmetry
      integer                        :: phase(nfuc)
      logical, allocatable           :: cell_in_C(:)
      integer                        :: ii, jj, kk, ll, orb_j, orb_k, orb_l
      double precision               :: ph
      logical                        :: found
      integer                        :: a1, b1, c1, d1, ip
      integer, external              :: inv_index

      !-----------------------------------------------------------------!
        
      !-----------------------------------------------------------------!

      !                   Precompute Bessel functions                   !

      call initialize_bessel_table()
      call gauss_legendre_points()
      
      !-----------------------------------------------------------------!

      call omp_set_dynamic(.false.)
      call omp_set_num_threads(omp_get_max_threads())

      allocate(ERI(nBas))
      
      

      if (.not. c_SAO) allocate(two_electron(nfuc,nBas,nBas,nBas))
      call classification_tor(number_of_atoms,nBas,geometry,atoms,ERI)

      !-- set up inversion
      call one_electron_inv_sym(nfuc, ERI, phase)
      call setup_cell_inversion(nx, ny, nz, cell_in_C)

      !$omp parallel
      if (omp_get_thread_num() == 0) then
        write(*,'(A)')    ""
        write(*,'(A,I0)') "Number of threads available: ", omp_get_num_threads()
      endif
      !$omp end parallel

      !$omp parallel
        !$omp single
          num_threads = omp_get_num_threads()
        !$omp end single
      !$omp end parallel

      ! Calculate optimal chunk size based on available threads

      if (num_threads <= 16) then
        optimal_chunk_size = 16           ! Good for laptops (8-16 cores)
      else if (num_threads <= 64) then
        optimal_chunk_size = 8            ! Good for medium clusters
      else 
        optimal_chunk_size = 1            ! Good for large clusters (128+ cores)
      end if

      start_time = omp_get_wtime()

      !-----------------------------------------------------------------!
      ! Pair list for Phase 1: only cell_j in the computed set C         !
      !-----------------------------------------------------------------!
      total_ij_pairs = 0
      total_ij_pairs_without = 0
      do i = 1, nfuc
        do j = i, nBas
          total_ij_pairs_without = total_ij_pairs_without + 1
          if (cell_in_C((j-1)/nfuc)) total_ij_pairs = total_ij_pairs + 1
        end do
      end do

      allocate(i_index(total_ij_pairs), j_index(total_ij_pairs))
      allocate(i_index_without(total_ij_pairs_without), j_index_without(total_ij_pairs_without))

      total_ij_pairs         = 0
      total_ij_pairs_without = 0
      do i = 1, nfuc
        do j = i, nBas
            total_ij_pairs_without = total_ij_pairs_without + 1 
            i_index_without(total_ij_pairs_without) = i
            j_index_without(total_ij_pairs_without) = j
          if (cell_in_C((j-1)/nfuc)) then
            total_ij_pairs = total_ij_pairs + 1
            i_index(total_ij_pairs) = i
            j_index(total_ij_pairs) = j
          end if
        end do
      end do

      !-- count Phase-1 integrals (for progress bar)

      num_total_int         = 0.d0
      do ij_index = 1, total_ij_pairs
        i = i_index(ij_index)
        j = j_index(ij_index)
        do k = 1, nBas
          do l = k, nBas
            if (i <= k .or. (i == k .and. j <= l)) then
              num_total_int = num_total_int + 1.d0
            end if
          end do
        end do
      end do

      num_total_int_without = 0.d0
      do ij_index = 1, total_ij_pairs_without
        i = i_index_without(ij_index)
        j = j_index_without(ij_index)
        do k = 1, nBas
          do l = k, nBas
            if (i <= k .or. (i == k .and. j <= l)) then
              num_total_int_without = num_total_int_without + 1.d0
            end if
          end do
        end do
      end do

      !-----------------------------------------------------------------!
      !              Schwarz screening precomputation                   !
      !   Q_schwarz(i,j) = sqrt(|(ij|ij)|)  for all unique pairs        !
      !   Done ONCE here, reused for every screening check below        !
      !-----------------------------------------------------------------!
      allocate(Q_schwarz(nBas, nBas))
      Q_schwarz = 0.d0
      write(*,'(A)')          'Precomputing Schwarz screening matrix...'
      write(outfile,'(6x,A)') 'Precomputing Schwarz screening matrix...'
      flush(outfile)

      do i = 1, nBas
        do j = i, nBas
          call ERI_integral_4_function_toroidal(ERI(i), ERI(j), ERI(i), ERI(j), val_iijj)
          Q_schwarz(i,j) = dsqrt(dabs(val_iijj))
          Q_schwarz(j,i) = Q_schwarz(i,j)
        end do
      end do

      write(*,'(A)')          'Schwarz matrix done'
      write(outfile,'(6x,A)') 'Schwarz matrix done'
      flush(outfile)

      !-- Phase-1 diagnostics
      
      write(*,'(A)') 'Sample Q_schwarz values:'
      write(*,'(A,ES12.4)') 'Max Q value = ', maxval(Q_schwarz)
      write(*,'(A,ES12.4)') 'Min Q value = ', minval(Q_schwarz)
      write(*,'(A,ES12.4)') 'Threshold   = ', schwarz_thresh

      ! --- Diagnostic: count how many integrals survive ---
      n_screened = 0
      do ij_index = 1, total_ij_pairs
        i = i_index(ij_index)
        j = j_index(ij_index)
        do k = 1, nBas
          do l = k, nBas
            if (i <= k .or. (i == k .and. j <= l)) then
              if (Q_schwarz(i,j) * Q_schwarz(k,l) < schwarz_thresh) then
                n_screened = n_screened + 1
              end if
            end if
          end do
        end do
      end do

      write(*,'(A)') ''

      n_phase1_computed = int(num_total_int,8) - n_screened


      write(*,'(6x,A,F8.2,A)') ''
      write(*,'(6x,A,F8.2,A)') '------------------------------------------------------------'
      write(*,'(6x,A,I0)')   'Total number of unique integrals without symmetry : ', int(num_total_int_without,16)
      write(*,'(6x,A,I0)')   'Total number of unique integrals with    symmetry : ', int(num_total_int,16)
      write(*,'(6x,A,I0)')   'Screened out by Schwarz + inversion   ', n_screened
      write(*,'(6x,A,I0)')   'Total unique integrals will computed: ', n_phase1_computed
      write(*,'(6x,A,F8.2,A)') 'Reduction with inversion            ', &
        100.d0 * dble(n_screened) / num_total_int + 100.d0 * dble(num_total_int /num_total_int_without) , ' %'
      write(*,'(6x,A,F8.2,A)') '------------------------------------------------------------'
      write(*,'(6x,A,F8.2,A)') ''
      flush(6)

      write(outfile,'(6x,A,F8.2,A)') ''
      write(outfile,'(6x,A,F8.2,A)') '------------------------------------------------------------'
      write(outfile,'(6x,A,I0)')   'Total number of unique integrals without symmetry : ', int(num_total_int_without,16)
      write(outfile,'(6x,A,I0)')   'Total number of unique integrals with    symmetry : ', int(num_total_int,16)
      write(outfile,'(6x,A,I0)')   'Screened out by Schwarz + inversion   ', n_screened
      write(outfile,'(6x,A,I0)')   'Total unique integrals will computed: ', n_phase1_computed
      write(outfile,'(6x,A,F8.2,A)') 'Reduction with inversion            ', &
        100.d0 * dble(n_screened) / num_total_int + 100.d0 * dble(num_total_int /num_total_int_without) , ' %'
      write(outfile,'(6x,A,F8.2,A)') '------------------------------------------------------------'
      write(outfile,'(6x,A,F8.2,A)') ''
      flush(outfile)
        
      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !

      if (c_SAO) then
        two_electron_integrals = 0.d0
      else
        two_electron = 0.d0
      end if

      !-----------------------------------------------------------------!
      ! Phase 1: compute all (i,j,k,l) with cell_j in computed set C     !
      !-----------------------------------------------------------------!

      if (c_SAO) then

      !$omp parallel do private(ij_index,i,j,k,l,val) &
      !$omp shared(ERI, i_index, j_index, Q_schwarz, two_electron_integrals) &
      !$omp schedule(dynamic,optimal_chunk_size)

      do ij_index = 1, total_ij_pairs
        i = i_index(ij_index)
        j = j_index(ij_index)
        do k = 1, nBas
          do l = k, nBas
            if (i <= k .or. (i == k .and. j <= l)) then

             ! ---- Schwarz screen: ONE multiply + ONE compare ----
             if (Q_schwarz(i,j) * Q_schwarz(k,l) < schwarz_thresh) cycle

              call ERI_integral_4_function_toroidal(ERI(i),ERI(j),ERI(k),ERI(l), val)         ! chemist notation

              if (c_physicist) then 
                ! Store in physicist notation:
                if (i <= nfuc) two_electron_integrals(i,k,j,l) = val
                if (j <= nfuc) two_electron_integrals(j,k,i,l) = val
                if (i <= nfuc) two_electron_integrals(i,l,j,k) = val
                if (j <= nfuc) two_electron_integrals(j,l,i,k) = val
                if (k <= nfuc) two_electron_integrals(k,i,l,j) = val
                if (l <= nfuc) two_electron_integrals(l,i,k,j) = val
                if (k <= nfuc) two_electron_integrals(k,j,l,i) = val
                if (l <= nfuc) two_electron_integrals(l,j,k,i) = val
              else
                if (i <= nfuc) two_electron_integrals(i,j,k,l) = val
                if (j <= nfuc) two_electron_integrals(j,i,k,l) = val
                if (i <= nfuc) two_electron_integrals(i,j,l,k) = val
                if (j <= nfuc) two_electron_integrals(j,i,l,k) = val
                if (k <= nfuc) two_electron_integrals(k,l,i,j) = val
                if (l <= nfuc) two_electron_integrals(l,k,i,j) = val
                if (k <= nfuc) two_electron_integrals(k,l,j,i) = val
                if (l <= nfuc) two_electron_integrals(l,k,j,i) = val
              end if

              !$omp critical
              integrals_done = integrals_done + 1.d0
              current_pct = int((integrals_done * 100.0d0) / dble(num_total_int - n_screened))
              if (current_pct > last_pct .and. mod(current_pct, 5) == 0) then
                  write(*,'(I3,"% done")') current_pct
                  flush(6)
                  last_pct = current_pct
              endif
              !$omp end critical
            end if
          end do
        end do
      end do
      !$omp end parallel do

      else 

      !$omp parallel do private(ij_index,i,j,k,l,val) &
      !$omp shared(two_electron, ERI, i_index, j_index, Q_schwarz) &
      !$omp schedule(dynamic,optimal_chunk_size)

      do ij_index = 1, total_ij_pairs
          i = i_index(ij_index)
          j = j_index(ij_index)
        do k = 1, nBas 
          do l = k, nBas
            if (i <= k .or. (i == k .and. j <= l)) then

             ! ---- Schwarz screen: ONE multiply + ONE compare ----
             if (Q_schwarz(i,j) * Q_schwarz(k,l) < schwarz_thresh) cycle

              call ERI_integral_4_function_toroidal(ERI(i),ERI(j),ERI(k),ERI(l), val)         ! chemist notation

              if (c_physicist) then 
                ! Store in physicist notation:
                if (i <= nfuc) two_electron(i,k,j,l) = val
                if (j <= nfuc) two_electron(j,k,i,l) = val
                if (i <= nfuc) two_electron(i,l,j,k) = val
                if (j <= nfuc) two_electron(j,l,i,k) = val
                if (k <= nfuc) two_electron(k,i,l,j) = val
                if (l <= nfuc) two_electron(l,i,k,j) = val
                if (k <= nfuc) two_electron(k,j,l,i) = val
                if (l <= nfuc) two_electron(l,j,k,i) = val
              else
                if (i <= nfuc) two_electron(i,j,k,l) = val
                if (j <= nfuc) two_electron(j,i,k,l) = val
                if (i <= nfuc) two_electron(i,j,l,k) = val
                if (j <= nfuc) two_electron(j,i,l,k) = val
                if (k <= nfuc) two_electron(k,l,i,j) = val
                if (l <= nfuc) two_electron(l,k,i,j) = val
                if (k <= nfuc) two_electron(k,l,j,i) = val
                if (l <= nfuc) two_electron(l,k,j,i) = val
              end if

              !$omp critical
              integrals_done = integrals_done + 1.d0
              current_pct = int((integrals_done * 100.0d0) / dble(num_total_int - n_screened))
              if (current_pct > last_pct .and. mod(current_pct, 5) == 0) then
                  write(*,'(I3,"% done")') current_pct
                  flush(6)
                  last_pct = current_pct
              endif
              !$omp end critical
            end if
          end do
        end do
      end do
      !$omp end parallel do

      end if

      deallocate(i_index, j_index)

      write(*,'(A)') ''
      write(*,'(A)') 'Phase 1 done. Filling remaining half via inversion symmetry...'
      write(*,'(A)') ''
      flush(6)

      !=================================================================!
      ! Phase 2: fill all (i,j,k,l) with cell_j NOT in C, via inversion !
      !   (ij|kl) = s_i s_j s_k s_l * (inv_i,inv_j|inv_k,inv_l)          !
      !=================================================================!

      if (c_SAO) then

      !$omp parallel do &
      !$omp private(i,j,k,l,ii,jj,kk,ll,orb_j,orb_k,orb_l,ph,found,a1,b1,c1,d1,ip,val) &
      !$omp shared(ERI, two_electron_integrals, phase, cell_in_C, Q_schwarz) &
      !$omp schedule(dynamic,optimal_chunk_size)

      do i = 1, nfuc
        do j = i, nBas
          if (cell_in_C((j-1)/nfuc)) cycle
          orb_j = mod(j-1, nfuc) + 1
          jj = inv_index(j, nfuc, nx, ny, nz)

          do k = 1, nBas
            orb_k = mod(k-1, nfuc) + 1
            kk = inv_index(k, nfuc, nx, ny, nz)
            do l = k, nBas
              if (.not. (i <= k .or. (i == k .and. j <= l))) cycle

              orb_l = mod(l-1, nfuc) + 1
              ll = inv_index(l, nfuc, nx, ny, nz)
              ii = inv_perm(i)

              ph = dble(phase(i)) * dble(phase(orb_j)) * &
                   dble(phase(orb_k)) * dble(phase(orb_l))

              found = .false.
              val   = 0.d0

              do ip = 1, 8
                select case(ip)
                case(1); a1=ii; b1=jj; c1=kk; d1=ll
                case(2); a1=ii; b1=jj; c1=ll; d1=kk
                case(3); a1=jj; b1=ii; c1=kk; d1=ll
                case(4); a1=jj; b1=ii; c1=ll; d1=kk
                case(5); a1=kk; b1=ll; c1=ii; d1=jj
                case(6); a1=kk; b1=ll; c1=jj; d1=ii
                case(7); a1=ll; b1=kk; c1=ii; d1=jj
                case(8); a1=ll; b1=kk; c1=jj; d1=ii
                end select
                if (a1 <= nfuc .and. b1 >= a1 .and. d1 >= c1 .and. &
                    a1 <= c1 .and. cell_in_C((b1-1)/nfuc)) then
                  if (c_physicist) then
                    val = ph * two_electron_integrals(a1,c1,b1,d1)
                  else
                    val = ph * two_electron_integrals(a1,b1,c1,d1)
                  end if
                  found = .true.
                  exit
                end if
              end do

              if (.not. found) then
                if (Q_schwarz(i,j) * Q_schwarz(k,l) < schwarz_thresh) cycle
                call ERI_integral_4_function_toroidal(ERI(i),ERI(j),ERI(k),ERI(l), val)
              end if

              if (c_physicist) then
                if (i <= nfuc) two_electron_integrals(i,k,j,l) = val
                if (j <= nfuc) two_electron_integrals(j,k,i,l) = val
                if (i <= nfuc) two_electron_integrals(i,l,j,k) = val
                if (j <= nfuc) two_electron_integrals(j,l,i,k) = val
                if (k <= nfuc) two_electron_integrals(k,i,l,j) = val
                if (l <= nfuc) two_electron_integrals(l,i,k,j) = val
                if (k <= nfuc) two_electron_integrals(k,j,l,i) = val
                if (l <= nfuc) two_electron_integrals(l,j,k,i) = val
              else
                if (i <= nfuc) two_electron_integrals(i,j,k,l) = val
                if (j <= nfuc) two_electron_integrals(j,i,k,l) = val
                if (i <= nfuc) two_electron_integrals(i,j,l,k) = val
                if (j <= nfuc) two_electron_integrals(j,i,l,k) = val
                if (k <= nfuc) two_electron_integrals(k,l,i,j) = val
                if (l <= nfuc) two_electron_integrals(l,k,i,j) = val
                if (k <= nfuc) two_electron_integrals(k,l,j,i) = val
                if (l <= nfuc) two_electron_integrals(l,k,j,i) = val
              end if
            end do
          end do
        end do
      end do
      !$omp end parallel do

      else

      !$omp parallel do &
      !$omp private(i,j,k,l,ii,jj,kk,ll,orb_j,orb_k,orb_l,ph,found,a1,b1,c1,d1,ip,val) &
      !$omp shared(ERI, two_electron, phase, cell_in_C, Q_schwarz) &
      !$omp schedule(dynamic,optimal_chunk_size)

      do i = 1, nfuc
        do j = i, nBas
          if (cell_in_C((j-1)/nfuc)) cycle
          orb_j = mod(j-1, nfuc) + 1
          jj = inv_index(j, nfuc, nx, ny, nz)

          do k = 1, nBas
            orb_k = mod(k-1, nfuc) + 1
            kk = inv_index(k, nfuc, nx, ny, nz)
            do l = k, nBas
              if (.not. (i <= k .or. (i == k .and. j <= l))) cycle

              orb_l = mod(l-1, nfuc) + 1
              ll = inv_index(l, nfuc, nx, ny, nz)
              ii = i

              ph = dble(phase(i)) * dble(phase(orb_j)) * &
                   dble(phase(orb_k)) * dble(phase(orb_l))

              found = .false.
              val   = 0.d0

              do ip = 1, 8
                select case(ip)
                case(1); a1=ii; b1=jj; c1=kk; d1=ll
                case(2); a1=ii; b1=jj; c1=ll; d1=kk
                case(3); a1=jj; b1=ii; c1=kk; d1=ll
                case(4); a1=jj; b1=ii; c1=ll; d1=kk
                case(5); a1=kk; b1=ll; c1=ii; d1=jj
                case(6); a1=kk; b1=ll; c1=jj; d1=ii
                case(7); a1=ll; b1=kk; c1=ii; d1=jj
                case(8); a1=ll; b1=kk; c1=jj; d1=ii
                end select
                if (a1 <= nfuc .and. b1 >= a1 .and. d1 >= c1 .and. &
                    a1 <= c1 .and. cell_in_C((b1-1)/nfuc)) then
                  if (c_physicist) then
                    val = ph * two_electron(a1,c1,b1,d1)
                  else
                    val = ph * two_electron(a1,b1,c1,d1)
                  end if
                  found = .true.
                  exit
                end if
              end do

              if (.not. found) then
                if (Q_schwarz(i,j) * Q_schwarz(k,l) < schwarz_thresh) cycle
                call ERI_integral_4_function_toroidal(ERI(i),ERI(j),ERI(k),ERI(l), val)
              end if

              if (c_physicist) then
                if (i <= nfuc) two_electron(i,k,j,l) = val
                if (j <= nfuc) two_electron(j,k,i,l) = val
                if (i <= nfuc) two_electron(i,l,j,k) = val
                if (j <= nfuc) two_electron(j,l,i,k) = val
                if (k <= nfuc) two_electron(k,i,l,j) = val
                if (l <= nfuc) two_electron(l,i,k,j) = val
                if (k <= nfuc) two_electron(k,j,l,i) = val
                if (l <= nfuc) two_electron(l,j,k,i) = val
              else
                if (i <= nfuc) two_electron(i,j,k,l) = val
                if (j <= nfuc) two_electron(j,i,k,l) = val
                if (i <= nfuc) two_electron(i,j,l,k) = val
                if (j <= nfuc) two_electron(j,i,l,k) = val
                if (k <= nfuc) two_electron(k,l,i,j) = val
                if (l <= nfuc) two_electron(l,k,i,j) = val
                if (k <= nfuc) two_electron(k,l,j,i) = val
                if (l <= nfuc) two_electron(l,k,j,i) = val
              end if
            end do
          end do
        end do
      end do
      !$omp end parallel do

      end if

      write(*,'(a)') ""
      write(*,'(a)') "*************************************************"
      write(*,'(a)') "* Two-electron integrals calculation completed  *"
      write(*,'(a)') "*************************************************"
      flush(6)
      write(*,'(a)') ""

      deallocate(Q_schwarz)
      deallocate(cell_in_C)

      end_time = omp_get_wtime()

      t = int(end_time - start_time)
      days= (t/86400)
      hours=mod(t,86400)/3600
      minutes=mod(mod(t,86400),3600)/60
      seconds=mod(mod(mod(t,86400),3600),60)

      write(outfile,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') '2 el-integrals calculation time = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec"
      write(outfile,"(a)") "" 
      FLUSH(outfile)

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      
      if (c_SAO) go to 200

      call cpu_time(start)
        call symmetry_of_integrals_ERI(nBas,nfuc,two_electron,two_electron_integrals)
      call cpu_time(end)

      time = int(end - start)
      days = (time/86400)
      hours=mod(time,86400)/3600
      minutes=mod(mod(time,86400),3600)/60
      seconds=mod(mod(mod(time,86400),3600),60)

      write(outfile,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') 'CPU time for Translational and 8 fold symmetry = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec"
      write(outfile,"(a)") ""
      write(outfile,"(a)") "8-fold symmetry with translation applied to integrals"
      write(outfile,"(a)") ""

      deallocate(two_electron)

200   continue
      deallocate(ERI)

contains

      !-----------------------------------------------------------------!
      ! Build "computed cell" table: cell_in_C(c) = .true. iff          !
      !   c <= inv_c   (so half + self-inverse cells are computed)      !
      !-----------------------------------------------------------------!
      subroutine setup_cell_inversion(nx_loc, ny_loc, nz_loc, cell_in_C_loc)
        integer, intent(in)              :: nx_loc, ny_loc, nz_loc
        logical, intent(out), allocatable :: cell_in_C_loc(:)
        integer :: ncell_loc, c, cx, cy, cz, icx, icy, icz, inv_c
        ncell_loc = nx_loc*ny_loc*nz_loc
        allocate(cell_in_C_loc(0:ncell_loc-1))
        do c = 0, ncell_loc-1
          cx = mod(c, nx_loc)
          cy = mod(c/nx_loc, ny_loc)
          cz = c/(nx_loc*ny_loc)
          icx = modulo(-cx, nx_loc)
          icy = modulo(-cy, ny_loc)
          icz = modulo(-cz, nz_loc)
          inv_c = icx + nx_loc*icy + nx_loc*ny_loc*icz
          cell_in_C_loc(c) = (c <= inv_c)
        end do
      end subroutine setup_cell_inversion

end subroutine ERI_integral_toroidal_1D_sym