module table_1d_lookup
      implicit none
      
      ! 1D table data
      double precision, allocatable      :: A_table(:)
      integer, allocatable               :: Nmax_table(:)
      logical                            :: table_loaded = .false.
      integer                            :: total_records
      
contains

      subroutine load_table_1d(filename)
            
            implicit none
            character(len=*), intent(in) :: filename
            integer                      :: file_unit, ios, i
            real                         :: start_time, end_time
            
            if (table_loaded) then
                  print*, '1D table already loaded!'
                  return
            end if
            
            call cpu_time(start_time)
            file_unit = 40  ! Different unit number from your 4D table
            
            open(unit=file_unit, file=filename, status='old', &
                 form='unformatted', access='stream', iostat=ios)
            
            if (ios /= 0) then
                  write(*,*) 'ERROR: Cannot open 1D binary file ', trim(filename)
                  write(*,*) 'Please generate the table first using number_of_terms program'
                  stop
            end if
            
            write(*,'(A)') ''
            write(*,'(A)') '************************************************************************'
            write(*,'(A,4x,A)') 'Loading 1D table from: ', trim(filename)
            write(*,'(A)') ''
            
            ! Get file size and calculate number of records
            inquire(unit=file_unit, size=total_records)
            ! Each record: 8 bytes (A) + 4 bytes (Nmax) = 12 bytes       
            total_records = total_records / 12
            
            if (total_records == 0) then
                  write(*,*) 'ERROR: File is empty!'
                  close(file_unit)
                  stop
            end if
            
            write(*,'(A,4x,I0,4x,A)') '  Found', total_records, 'records in table'
            
            ! Allocate arrays
            allocate(A_table(total_records))
            allocate(Nmax_table(total_records))
            
            ! Read all data into arrays
            rewind(file_unit)
            do i = 1, total_records
                  read(file_unit, iostat=ios) A_table(i), Nmax_table(i)
                  if (ios /= 0) then
                        print*, 'ERROR: Failed to read record', i
                        close(file_unit)
                        stop
                  end if
            end do
            
            close(file_unit)
            table_loaded = .true.
            
            call cpu_time(end_time)
            
            
            write(*,'(A,4x,ES10.0,4x,A,ES10.0)') '  A range:', A_table(1), '  to  ', A_table(total_records)
            write(*,'(A,F12.8,4x,A)') '  Load time   :', end_time - start_time, 'seconds'
            write(*,'(A,F12.8,4x,A)') '  Memory usage:', (total_records * 12.d0) / (1024.d0*1024.d0), 'MB'
            write(*,'(A)') ''
            write(*,'(A)') '1D table loaded successfully!'
            write(*,'(A)') '************************************************************************'
            write(*,'(A)') ''
            
      end subroutine load_table_1d
      
      
      integer function get_Nmax_from_bessel_1D(A_target)
            
            implicit none 
            double precision, intent(in) :: A_target
            integer                      :: low, high, mid
            
            if (.not. table_loaded) then
                  print*, 'ERROR: 1D table not loaded! Call load_table_1d() first.'
                  get_Nmax_from_bessel_1D = -1
                  return
            end if
            
            ! Check boundaries
            if (A_target <= A_table(1)) then
                  get_Nmax_from_bessel_1D = Nmax_table(1)
                  if (A_target < A_table(1)) then
                        print*, 'Warning: A=', A_target, 'below table range. Using smallest A=', A_table(1)
                  end if
                  return
            end if
            
            if (A_target >= A_table(total_records)) then
                  get_Nmax_from_bessel_1D = Nmax_table(total_records)
                  if (A_target > A_table(total_records)) then
                        print*, 'Warning: A=', A_target, 'above table range. Using largest A=', A_table(total_records)
                  end if
                  return
            end if
            
            ! Binary search
            low = 1
            high = total_records
            
            do while (high - low > 1)
                  mid = (low + high) / 2
                  if (A_table(mid) <= A_target) then
                        low = mid
                  else
                        high = mid
                  end if
            end do
            
            ! Return the closest value
            if (abs(A_table(low) - A_target) <= abs(A_table(high) - A_target)) then
                  get_Nmax_from_bessel_1D = Nmax_table(low)
            else
                  get_Nmax_from_bessel_1D = Nmax_table(high)
            end if
            
      end function get_Nmax_from_bessel_1D
      
end module table_1d_lookup