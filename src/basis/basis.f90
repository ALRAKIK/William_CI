subroutine basis(number_of_atoms,charge,atoms)

      use files 
      use atom_basis

      implicit none 

      integer,intent(in)      :: number_of_atoms
      integer,intent(in)      :: charge(number_of_atoms)
      type(atom),intent(out)  :: atoms(number_of_atoms)

      integer                 :: i 

      do i = 1, number_of_atoms
        atoms(i)%charge = charge(i)
        call read_atom_shells_gen(atoms(i))
      end do

end subroutine basis 

subroutine read_atom_shells_gen(at)

      use atom_basis

      implicit none

      type(atom), intent(inout) :: at

      character(len=1000) :: line, key
      character(len=10)   :: tag1, tag2
      character(len=1)    :: shl
      integer             :: ios, ne, nc, i, j, u
      logical             :: inside, found
      double precision, allocatable :: ex(:), co(:,:)

      write(tag1,'(A,I0)') "A ", at%charge
      write(tag2,'(A,I0)') "a ", at%charge

      at%num_s_function = 0 ; at%num_exponent_s = 0
      at%num_p_function = 0 ; at%num_exponent_p = 0
      at%num_d_function = 0 ; at%num_exponent_d = 0
      at%num_f_function = 0 ; at%num_exponent_f = 0
      at%num_g_function = 0 ; at%num_exponent_g = 0

      open(newunit=u, file="Basis", status='old', action='read', iostat=ios)
      if (ios /= 0) then
        write(*,*) 'basis_tor_gen: cannot open the file Basis'
        stop
      end if

      inside = .false. ; found = .false.
      do
        read(u,'(A)',iostat=ios) line
        if (ios /= 0) exit
        if (.not. inside) then
          if (trim(adjustl(line)) == trim(tag1) .or. trim(adjustl(line)) == trim(tag2)) then
            inside = .true. ; found = .true.
          end if
          cycle
        end if
        key = adjustl(line)
        if (key(1:1) == 'A' .or. key(1:1) == 'a') exit           ! next atom
        if (key(1:1) /= '$') cycle
        if (index(key, '-TYPE FUNCTIONS') == 0) cycle             ! comment line
        shl = key(3:3)
        read(u,*) ne , nc
        allocate(ex(ne), co(ne,nc))
        do i = 1 , ne
          read(u,*) ex(i) , (co(i,j), j=1,nc)
        end do
        select case (shl)
        case ('S','s')
          at%num_exponent_s = ne ; at%num_s_function = nc
          at%exponent_s = ex     ; at%coefficient_s = co
        case ('P','p')
          at%num_exponent_p = ne ; at%num_p_function = nc
          at%exponent_p = ex     ; at%coefficient_p = co
        case ('D','d')
          at%num_exponent_d = ne ; at%num_d_function = nc
          at%exponent_d = ex     ; at%coefficient_d = co
        case ('F','f')
          at%num_exponent_f = ne ; at%num_f_function = nc
          at%exponent_f = ex     ; at%coefficient_f = co
        case ('G','g')
          at%num_exponent_g = ne ; at%num_g_function = nc
          at%exponent_g = ex     ; at%coefficient_g = co
        case default
          write(*,*) 'basis_tor_gen: shell type not supported: ', trim(key)
          stop
        end select
        deallocate(ex, co)
      end do
      close(u)

      ! absent shells -> zero-size arrays (size() is then defined, e.g. in lib_int)
      if (.not. allocated(at%exponent_s)) allocate(at%exponent_s(0), at%coefficient_s(0,0))
      if (.not. allocated(at%exponent_p)) allocate(at%exponent_p(0), at%coefficient_p(0,0))
      if (.not. allocated(at%exponent_d)) allocate(at%exponent_d(0), at%coefficient_d(0,0))
      if (.not. allocated(at%exponent_f)) allocate(at%exponent_f(0), at%coefficient_f(0,0))
      if (.not. allocated(at%exponent_g)) allocate(at%exponent_g(0), at%coefficient_g(0,0))

      if (.not. found) then
        write(*,'(A,I0,A)') ' basis_tor_gen: no basis for atom "A ', at%charge, '" in the file Basis'
        stop
      end if

end subroutine read_atom_shells_gen