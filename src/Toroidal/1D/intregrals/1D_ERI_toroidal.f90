subroutine ERI_integral_toroidal_1D(number_of_atoms,geometry,nBas,atoms,two_electron_integrals)

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

      integer                        :: i , j , k , l
      integer                        :: number_of_atoms
      integer                        :: nBas
      type(atom)                     :: atoms(number_of_atoms)
      type(ERI_function),allocatable :: ERI  (:)

      double precision               :: geometry(number_of_atoms,3)
      double precision,allocatable   :: two_electron(:,:,:,:)
      double precision               :: num_total_int

      ! --------------------------------------------------------------- !
      
      double precision,intent(out)   :: two_electron_integrals(nBas,nBas,nBas,nBas)

      ! ----------------------    Time     ---------------------------- !
      integer                        :: days, hours, minutes, seconds , t , time 
      double precision               :: start,end
      double precision               :: start_time, end_time
      !-----------------------------------------------------------------!
      integer                        :: total_ij_pairs, ij_index
      integer, allocatable           :: i_index(:), j_index(:)
      integer                        :: num_threads, optimal_chunk_size
      !-----------------------------------------------------------------!

      integer                        :: current_pct 
      integer                        :: last_pct = -1
      double precision               :: integrals_done = 0.d0
      double precision               :: val

      !-----------------------------------------------------------------!
      !                      Schwarz screening                          !
      
      double precision, allocatable :: Q_schwarz(:,:)
      double precision              :: val_iijj
      double precision, parameter   :: schwarz_thresh = 1.d-14
      integer(8)                    :: n_screened, n_computed

      !-----------------------------------------------------------------!
        
      !-----------------------------------------------------------------!

      !                   Precompute Bessel functions                   !

      call initialize_bessel_table()
      call gauss_legendre_points()
      
      !-----------------------------------------------------------------!

      call omp_set_dynamic(.false.)
      call omp_set_num_threads(omp_get_max_threads())


      allocate(ERI(nBas))
      allocate(two_electron(nfuc,nBas,nBas,nBas))

      call classification_tor(number_of_atoms,nBas,geometry,atoms,ERI)

      !$omp parallel
      if (omp_get_thread_num() == 0) then
        print *, "Running with", omp_get_num_threads(), "threads"
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
      do i = 1, nfuc
        do j = i, nBas
          total_ij_pairs = total_ij_pairs + 1
        end do
      end do

      allocate(i_index(total_ij_pairs), j_index(total_ij_pairs))

      total_ij_pairs = 0
      do i = 1, nfuc
        do j = i, nBas
          total_ij_pairs = total_ij_pairs + 1
          i_index(total_ij_pairs) = i
          j_index(total_ij_pairs) = j
        end do
      end do

      num_total_int = 0.d0
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

      write(outfile,'(a,I0,a)') 'Need to  compute ', int(num_total_int,16), ' unique integrals (Two electron Integrals)'
      write(*      ,'(a,I0,a)') 'Need to  compute ', int(num_total_int,16), ' unique integrals (Two electron Integrals)'
      write(outfile,*) ''
      flush(outfile)


      !-----------------------------------------------------------------!
      !              Schwarz screening precomputation                   !
      !   Q_schwarz(i,j) = sqrt(|(ij|ij)|)  for all unique pairs        !
      !   Done ONCE here, reused for every screening check below        !
      !-----------------------------------------------------------------!

      allocate(Q_schwarz(nBas, nBas))
      Q_schwarz = 0.d0
      write(*,'(A)') 'Precomputing Schwarz screening matrix...'
      write(outfile,'(A)') 'Precomputing Schwarz screening matrix...'


      do i = 1, nBas
        do j = i, nBas
          call ERI_integral_4_function_toroidal(ERI(i), ERI(j), ERI(i), ERI(j), val_iijj)
          Q_schwarz(i,j) = dsqrt(dabs(val_iijj))
          Q_schwarz(j,i) = Q_schwarz(i,j)
        end do
      end do

      write(*,'(A)') 'Schwarz matrix done.'
      write(outfile,'(A)') 'Schwarz matrix done.'
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
        do k = 1, nBas
          do l = k, nBas
            if (i <= k .or. (i == k .and. j <= l)) then
              n_computed = n_computed + 1
              if (Q_schwarz(i,j) * Q_schwarz(k,l) < schwarz_thresh) then
                n_screened = n_screened + 1
              end if
            end if
          end do
        end do
      end do

      write(*,'(A,I0)')   'Total unique integrals:    ', n_computed
      write(*,'(A,I0)')   'Screened out by Schwarz:   ', n_screened
      write(*,'(A,F8.2,A)') 'Reduction:                 ', &
        100.d0 * dble(n_screened) / dble(n_computed), ' %'
      write(outfile,'(A,I0)')   'Total unique integrals:    ', n_computed
      write(outfile,'(A,I0)')   'Screened out by Schwarz:   ', n_screened
      write(outfile,'(A,F8.2,A)') 'Reduction:                 ', &
        100.d0 * dble(n_screened) / dble(n_computed), ' %'
      write(outfile,'(A)') ''
      flush(outfile)
      !-----------------------------------------------------------------!
      
      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !

      two_electron = 0.d0 

      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !

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
             if (Q_schwarz(i,j) * Q_schwarz(k,l) < schwarz_thresh) then
               two_electron(i,j,k,l) = 0.d0   ! screened integrals are zero
               cycle
             end if

              call ERI_integral_4_function_toroidal(ERI(i),ERI(j),ERI(k),ERI(l), val)

              if (c_physicist) then 
                ! Store in physicist notation:
                two_electron(i,k,j,l) = val    ! <ik|jl> = (ij|kl)_chem
                two_electron(i,l,j,k) = val    ! <il|jk> = (ij|lk)_chem = (ij|kl)_chem
              else 
                two_electron(i,j,k,l) = val
                two_electron(i,j,l,k) = val
              end if

               !$omp critical
               integrals_done = integrals_done + 1.d0
               current_pct = int((integrals_done * 100.0d0) / num_total_int)
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

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      if (c_save) then 
        call save_two_electron_integrals(nfuc,nBas,two_electron)
      end if

      call cpu_time(start)
        call symmetry_of_integrals_ERI_mod(nBas,nfuc,two_electron,two_electron_integrals)
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

      deallocate(ERI)
      deallocate(two_electron)

end subroutine ERI_integral_toroidal_1D