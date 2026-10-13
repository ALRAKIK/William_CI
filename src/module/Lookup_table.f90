module table_lookup_module

      use constants_module
      
      implicit none
      
      ! Grid parameters
      integer                     :: n_points
      double precision            :: A_min, A_max
      double precision            :: log_min, log_max, log_step
      integer, parameter          :: n_phi = 20
      
      ! 4D lookup table (allocatable for flexibility)
      integer, allocatable        :: Nmax_table(:,:,:,:)
      logical                     :: table_loaded = .false.
      
contains

      subroutine load_table(filename)

      character(len=*), intent(in) :: filename
      integer                      :: file_unit, ios
      integer                      :: i, j, k, l, Nmax_val
      real                         :: start_time, end_time
            
      if (table_loaded) then
            print*, 'Table already loaded!'
            return
      end if
            
      call cpu_time(start_time)
            
      file_unit = 30
            
      ! =============================================================== !
      ! OPEN BINARY FILE
      ! =============================================================== !


      open(unit=file_unit, file=filename, status='old', &
           form='unformatted', access='stream', iostat=ios)
      
      if (ios /= 0) then
            print*, 'ERROR: Cannot open binary file ', trim(filename)
            print*, 'Make sure you generated table_data.bin'
            stop
      end if
      
      write(*,"(A)") ""
      write(*,"(A)") "************************************************************************" 
      write(*,"(A,A,A)")  'Loading binary table from: ', trim(filename) , " For the S function maximum terms"
      write(*,"(A)") ""
            
      ! ============================================================
      ! READ HEADER
      ! ============================================================
      read(file_unit, iostat=ios) n_points
      read(file_unit, iostat=ios) A_min
      read(file_unit, iostat=ios) A_max
      read(file_unit, iostat=ios) log_min
      read(file_unit, iostat=ios) log_max
      read(file_unit, iostat=ios) log_step
            
      if (ios /= 0) then
            print*, 'ERROR: Failed to read header'
            close(file_unit)
            stop
      end if
            
      write(*,"(A)") "The table is loaded with the following parameters:"
      write(*,"(A)")         '  '
      write(*,"(A,4x,I0)")     '  n_points =', n_points
      write(*,"(A,ES10.0)")    '  A_Max    =', A_max
      write(*,"(A,ES10.0)")    '  log_step =', log_step
      write(*,"(A)")         '  '
            
      ! =============================================================== !
      ! ALLOCATE TABLE
      ! =============================================================== !

      allocate(Nmax_table(0:n_points, 0:n_points, 0:n_points, 0:n_phi-1))
            
      ! =============================================================== !
      ! READ ALL DATA IN ONE SHOT (SUPER FAST!)
      ! =============================================================== !

      do i = 0, n_points
        do j = 0, n_points
          do k = 0, n_points
            do l = 0, n_phi-1
              read(file_unit, iostat=ios) Nmax_val
              
              if (ios /= 0) then
                    print*, 'ERROR: Unexpected end of file at', i, j, k, l
                    close(file_unit)
                    stop
              end if
              
              Nmax_table(i, j, k, l) = max(Nmax_val, 1)
            end do
          end do
        end do
      end do
            
      close(file_unit)

      table_loaded = .true.
            
      call cpu_time(end_time)
            
      
      write(*,"(A,4x,I0)")      '  Entries:', (n_points+1)**3 * n_phi
      write(*,"(A,F12.8,4x,A)") '  Load time   :', end_time - start_time, 'seconds'
      write(*,"(A,F12.8,4x,A)") '  Memory usage:', &
            (n_points+1)**3 * n_phi * 4 / (1024.0*1024.0), 'MB'
      write(*,"(A)") ''
      write(*,"(A)")         '  Binary table loaded successfully!'
      write(*,"(A)") ''
      write(*,"(A)") '************************************************************************'
            
      end subroutine load_table
      
      
      integer function get_Nmax(A, B, C, phi)

      double precision, intent(in) :: A, B, C, phi
      integer                      :: idx_A, idx_B, idx_C, idx_phi
      double precision             :: A_safe, B_safe, C_safe, phi_safe
            
      if (.not. table_loaded) then
            print*, 'ERROR: Table not loaded! Call load_table() first.'
            get_Nmax = 100
            return
      end if
            
      ! Clamp to valid range
      A_safe = max(A_min, min(A_max, A))
      B_safe = max(A_min, min(A_max, B))
      C_safe = max(A_min, min(A_max, C))
      phi_safe = mod(phi, 2.0d0 * pi)

      if (phi_safe < 0.0d0) phi_safe = phi_safe + 2.0d0 * pi
      
      ! Calculate indices
      idx_A = nint((log10(A_safe) - log_min) / log_step)
      idx_B = nint((log10(B_safe) - log_min) / log_step)
      idx_C = nint((log10(C_safe) - log_min) / log_step)
      idx_phi = nint(phi_safe / (pi / 10.0d0))
      
      ! Clamp to array bounds
      idx_A = max(0, min(n_points, idx_A))
      idx_B = max(0, min(n_points, idx_B))
      idx_C = max(0, min(n_points, idx_C))
      idx_phi = max(0, min(n_phi-1, idx_phi))
      
      ! Lookup
      get_Nmax = Nmax_table(idx_A, idx_B, idx_C, idx_phi)
      
      if (get_Nmax < 1) get_Nmax = 10
            
      end function get_Nmax
      
      subroutine cleanup_table()
            ! Call this at program end if needed
            if (allocated(Nmax_table)) deallocate(Nmax_table)
            table_loaded = .false.
      end subroutine cleanup_table
      
end module table_lookup_module