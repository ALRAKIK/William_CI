subroutine ERI_integral(number_of_atoms,number_of_functions,geometry,atoms,two_electron)

      use files 
      use omp_lib
      use atom_basis
      use classification_ERI
      use keywords

      implicit none 

      !-----------------------------------------------------------------!

      type(atom),intent(in)          :: atoms(number_of_atoms)
      integer,intent(in)             :: number_of_atoms


      type(ERI_function),allocatable :: ERI  (:)
      integer                        :: i , j , k , l
      integer                        :: number_of_functions
      double precision               :: geometry(number_of_atoms,3)
      double precision               :: two_electron(number_of_functions,number_of_functions,number_of_functions,number_of_functions)
      double precision               :: value 

      double precision               :: num_total_int


      ! ----------------------    Time     ---------------------------- !
      double precision               :: start_time, end_time
      integer                        :: days, hours, minutes, seconds , t 
      !-----------------------------------------------------------------!
      integer                        :: total_ij_pairs, ij_index
      integer, allocatable           :: i_index(:), j_index(:)
      integer                        :: num_threads, optimal_chunk_size
      !-----------------------------------------------------------------!
      integer                        :: current_pct 
      integer                        :: last_pct = -1
      double precision               :: integrals_done = 0.d0
      !-----------------------------------------------------------------!
      !                      Schwarz screening                          !
      
      double precision, allocatable :: Q_schwarz(:,:)
      double precision              :: val_iijj
      double precision, parameter   :: schwarz_thresh = 1.d-14
      integer                       :: n_screened, n_computed
      
      !-----------------------------------------------------------------!

      call omp_set_dynamic(.false.)
      call omp_set_num_threads(omp_get_max_threads())
      
      allocate(ERI(number_of_functions))
      
      call classification(number_of_atoms,number_of_functions,geometry,atoms,ERI)

      !$omp parallel
      if (omp_get_thread_num() == 0) then
        write(*,'(A)')    ""
        write(*,'(A,I0)') "Number of threads available: ", omp_get_num_threads()
      endif
      !$omp end parallel

      !$omp    parallel
          !$omp     single
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

      ! Precompute all i-j pairs for manual collapse

      total_ij_pairs = 0
      do i = 1, number_of_functions
        do j = i, number_of_functions
          total_ij_pairs = total_ij_pairs + 1
        end do
      end do

      allocate(i_index(total_ij_pairs), j_index(total_ij_pairs))

      total_ij_pairs = 0
      do i = 1, number_of_functions
        do j = i, number_of_functions
          total_ij_pairs = total_ij_pairs + 1
          i_index(total_ij_pairs) = i
          j_index(total_ij_pairs) = j
        end do
      end do

      num_total_int = 0.d0
      do ij_index = 1, total_ij_pairs
        i = i_index(ij_index)
        j = j_index(ij_index)
        do k = 1, number_of_functions
          do l = k, number_of_functions
            if (i <= k .or. (i == k .and. j <= l)) then
              num_total_int = num_total_int + 1.d0
            end if 
          end do
        end do
      end do

      !-----------------------------------------------------------------!
      !              Schwarz screening precomputation                   !
      !   Q_schwarz(i,j) = sqrt(|(ij|ij)|)  for all unique pairs        !
      !   Done ONCE here, reused for every screening check below        !
      !-----------------------------------------------------------------!

      allocate(Q_schwarz(number_of_functions, number_of_functions))
      Q_schwarz = 0.d0
      write(*,'(A)')          'Precomputing Schwarz screening matrix...'
      write(outfile,'(6x,A)') 'Precomputing Schwarz screening matrix...'


      do i = 1, number_of_functions
        do j = i, number_of_functions
          call ERI_integral_4_function(ERI(i), ERI(j), ERI(i), ERI(j), val_iijj)
          Q_schwarz(i,j) = dsqrt(dabs(val_iijj))
          Q_schwarz(j,i) = Q_schwarz(i,j)
        end do
      end do

      write(*,'(A)')       'Schwarz matrix done'
      write(outfile,'(6x,A)') 'Schwarz matrix done'
      flush(outfile)

      write(*,'(A)') 'Sample Q_schwarz values:'
      write(*,'(A,ES12.4)') 'Max Q value = ', maxval(Q_schwarz)
      write(*,'(A,ES12.4)') 'Min Q value = ', minval(Q_schwarz)
      write(*,'(A,ES12.4)') 'Threshold   = ', schwarz_thresh

      ! --- Diagnostic: count how many integrals survive ---

      n_screened = 0
      n_computed = 0
      do ij_index = 1, total_ij_pairs
        i = i_index(ij_index)
        j = j_index(ij_index)
        do k = 1, number_of_functions
          do l = k, number_of_functions
            if (i <= k .or. (i == k .and. j <= l)) then
              if (Q_schwarz(i,j) * Q_schwarz(k,l) < schwarz_thresh) then
                n_screened = n_screened + 1
              else
                n_computed = n_computed + 1
              end if
            end if
          end do
        end do
      end do
      
      write(*,'(A,I0)')   'Total number of unique integrals :    ', int(num_total_int,16)
      write(*,'(A,I0)')   'Screened out by Schwarz:              ', n_screened
      write(*,'(A,I0)')   'Total unique integrals will computed: ', n_computed
      write(*,'(A,F8.2,A)') 'Reduction:                 ', &
        100.d0 * dble(n_screened) / num_total_int, ' %'

      write(*,'(A)') ''

      write(outfile,'(6x,A,I0)')   'Total number of unique integrals :    ', int(num_total_int,16)
      write(outfile,'(6x,A,I0)')   'Screened out by Schwarz:              ', n_screened
      write(outfile,'(6x,A,I0)')   'Total unique integrals will computed: ', n_computed
      write(outfile,'(6x,A,F8.2,A)') 'Reduction:                 ', &
        100.d0 * dble(n_screened) / num_total_int, ' %'

      write(outfile,'(A)') ''
      flush(outfile)
        
      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !

      two_electron = 0.d0 

      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !

      !$omp parallel do private(ij_index,i,j,k,l,value) &
      !$omp shared(two_electron, ERI, i_index, j_index, Q_schwarz) &
      !$omp schedule(dynamic,optimal_chunk_size)

      do ij_index = 1, total_ij_pairs
          i = i_index(ij_index)
          j = j_index(ij_index)
        do k = 1, number_of_functions
          do l = k, number_of_functions

            if (i <= k .or. (i == k .and. j <= l)) then

                ! ---- Schwarz screen: ONE multiply + ONE compare ----

                if (Q_schwarz(i,j) * Q_schwarz(k,l) < schwarz_thresh) cycle

                call ERI_integral_4_function(ERI(i),ERI(j),ERI(k),ERI(l),value)

                if (c_physicist) then 

                  two_electron(i,k,j,l) = value
                  two_electron(j,k,i,l) = value
                  two_electron(i,l,j,k) = value
                  two_electron(j,l,i,k) = value
                  two_electron(k,i,l,j) = value
                  two_electron(l,i,k,j) = value
                  two_electron(k,j,l,i) = value
                  two_electron(l,j,k,i) = value

                else 

                  two_electron(i,j,k,l) = value
                  two_electron(j,i,k,l) = value
                  two_electron(i,j,l,k) = value
                  two_electron(j,i,l,k) = value
                  two_electron(k,l,i,j) = value
                  two_electron(l,k,i,j) = value
                  two_electron(k,l,j,i) = value
                  two_electron(l,k,j,i) = value

                end if 

              !$omp critical
              integrals_done = integrals_done + 1.d0
              current_pct = int((integrals_done * 100.0d0) / n_computed)
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

      write(*,'(a)') ""
      write(*,'(a)') "*************************************************"
      write(*,'(a)') "* Two-electron integrals calculation completed  *"
      write(*,'(a)') "*************************************************"
      flush(6)
      write(*,'(a)') ""

      deallocate(i_index, j_index)
      deallocate(Q_schwarz)

      end_time = omp_get_wtime()
    
      t = int(end_time - start_time)
      days= (t/86400)
      hours=mod(t,86400)/3600
      minutes=mod(mod(t,86400),3600)/60
      seconds=mod(mod(mod(t,86400),3600),60)

      write(outfile,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') '2 el-integrals calculation time = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec"
      write(outfile,"(a)") "" 
      FLUSH(outfile)


end subroutine