subroutine read_basis_class_tor(atom_type)
        
      use files
      use atom_basis
      
      implicit none

      integer             :: i , j 
      type(atom)          :: atom_type
      character(len=100)  :: lines
      character(len=10)   :: atom_type_charge
      logical             :: found_s, found_p , found_d  , found_f , found_g
      logical             :: basis_problem
      
      write(atom_type_charge,'(A,I0)') "A ", atom_type%charge 

      atom_type%num_s_function = 0 ; atom_type%num_exponent_s = 0
      atom_type%num_p_function = 0 ; atom_type%num_exponent_p = 0
      atom_type%num_d_function = 0 ; atom_type%num_exponent_d = 0
      atom_type%num_f_function = 0 ; atom_type%num_exponent_f = 0
      atom_type%num_g_function = 0 ; atom_type%num_exponent_g = 0
      
      found_s                  = .false.
      found_p                  = .false.
      found_d                  = .false.
      found_f                  = .false.
      found_g                  = .false.
    
      open(1,file=trim(tmp_file_name)//"/Basis_normalized")

        do 

          read(1,'(A)',end=3) lines 

          !------------------------------------------------
          ! S-TYPE
          !------------------------------------------------

          if (lines == atom_type_charge) then 
            do 
              read(1,'(A)',end=3) lines 
              if (lines == "$ S-TYPE FUNCTIONS") then
                found_s = .true.
                read(1,*) atom_type%num_exponent_s , atom_type%num_s_function
                allocate(atom_type%exponent_s(atom_type%num_exponent_s))
                allocate(atom_type%coefficient_s(atom_type%num_exponent_s,atom_type%num_s_function))
               do i = 1 , atom_type%num_exponent_s
                 read(1,*) atom_type%exponent_s(i) , (atom_type%coefficient_s(i,j),j=1,atom_type%num_s_function)
               end do 
              end if 
              !------------------------------------------------
              ! PX-TYPE
              !------------------------------------------------
              if (lines == "$ PX-TYPE FUNCTIONS") then
                found_p = .true.
                read(1,*) atom_type%num_exponent_p , atom_type%num_p_function
                allocate(atom_type%exponent_p(atom_type%num_exponent_p))
                allocate(atom_type%coefficient_px(atom_type%num_exponent_p,atom_type%num_p_function))
               do i = 1 , atom_type%num_exponent_p
                 read(1,*) atom_type%exponent_p(i) , (atom_type%coefficient_px(i,j),j=1,atom_type%num_p_function)
               end do
              end if
              !------------------------------------------------
              ! PY-TYPE
              !------------------------------------------------
              if (lines == "$ PY-TYPE FUNCTIONS") then
                found_p = .true.
                read(1,*) atom_type%num_exponent_p , atom_type%num_p_function
                allocate(atom_type%coefficient_py(atom_type%num_exponent_p,atom_type%num_p_function))
               do i = 1 , atom_type%num_exponent_p
                 read(1,*) atom_type%exponent_p(i) , (atom_type%coefficient_py(i,j),j=1,atom_type%num_p_function)
               end do
              end if
              !------------------------------------------------
              ! PZ-TYPE
              !------------------------------------------------
              if (lines == "$ PZ-TYPE FUNCTIONS") then
                found_p = .true.
                read(1,*) atom_type%num_exponent_p , atom_type%num_p_function
                allocate(atom_type%coefficient_pz(atom_type%num_exponent_p,atom_type%num_p_function))
               do i = 1 , atom_type%num_exponent_p
                 read(1,*) atom_type%exponent_p(i) , (atom_type%coefficient_pz(i,j),j=1,atom_type%num_p_function)
               end do
              end if

              if (lines == "$ D-TYPE FUNCTIONS") then
                found_d = .true.
                read(1,*) atom_type%num_exponent_d , atom_type%num_d_function
                allocate(atom_type%exponent_d(atom_type%num_exponent_d))
                allocate(atom_type%coefficient_d(atom_type%num_exponent_d,atom_type%num_d_function))
               do i = 1 , atom_type%num_exponent_d
                 read(1,*) atom_type%exponent_d(i) , (atom_type%coefficient_d(i,j),j=1,atom_type%num_d_function)
               end do
              end if
              if (lines == "$ F-TYPE FUNCTIONS") then
                found_f = .true.
                read(1,*) atom_type%num_exponent_f , atom_type%num_f_function
                allocate(atom_type%exponent_f(atom_type%num_exponent_f))
                allocate(atom_type%coefficient_f(atom_type%num_exponent_f,atom_type%num_f_function))
               do i = 1 , atom_type%num_exponent_f
                 read(1,*) atom_type%exponent_f(i) , (atom_type%coefficient_f(i,j),j=1,atom_type%num_f_function)
               end do 
              end if
              if (lines == "$ G-TYPE FUNCTIONS") then
                found_g = .true.
                read(1,*) atom_type%num_exponent_g , atom_type%num_g_function
                allocate(atom_type%exponent_g(atom_type%num_exponent_g))
                allocate(atom_type%coefficient_g(atom_type%num_exponent_g,atom_type%num_g_function))
               do i = 1 , atom_type%num_exponent_g
                 read(1,*) atom_type%exponent_g(i) , (atom_type%coefficient_g(i,j),j=1,atom_type%num_g_function)
               end do 
              end if
              if (adjustl(lines(1:1))=="A" .or. adjustl(lines(1:1))=="a") exit
            end do 

            if (.not. found_p) then
              atom_type%num_exponent_p = 0
              atom_type%num_p_function = 0
            end if

            if (.not. found_d) then
              atom_type%num_exponent_d = 0
              atom_type%num_d_function = 0
            end if

            if (.not. found_f) then
              atom_type%num_exponent_f = 0
              atom_type%num_f_function = 0
            end if

            if (.not. found_g) then
              atom_type%num_exponent_g = 0
              atom_type%num_g_function = 0
            end if

            exit

          end if 

        end do 

        basis_problem = .not. (found_s .and. found_p .and. found_d .and. found_f .and. found_g) 

        if (basis_problem) then 
          write(outfile,'(6x,A)') "Problem in the basis set:"
          write(outfile,'(6x,A)') " * -  check the basis  "
          write(outfile,'(6x,A)') " * -  check the syntax "

          write(*,'(6x,A)') "Problem in the basis set:"
          write(*,'(6x,A)') " * -  check the basis  "
          write(*,'(6x,A)') " * -  check the syntax "
        end if 

3      close(1)

end subroutine read_basis_class_tor