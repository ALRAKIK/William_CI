subroutine extract_basis(number_of_atoms,charge)

      use files
      implicit none 


      integer              :: number_of_atoms
      integer              :: charge(number_of_atoms)


      integer              :: i 
      character(len=1000)  :: lines

      character(len=10)    :: atom_type1(number_of_atoms)
      character(len=10)    :: atom_type2(number_of_atoms)
      
      ! --------------------------------------------------------------- !

      do i = 1 , number_of_atoms
        write(atom_type1(i),'(A,I0)') "A ", charge(i)
      end do 

      do i = 1 , number_of_atoms
        write(atom_type2(i),'(A,I0)') "a ", charge(i)
      end do 
      
      open(1,file="Basis")
      open(2,file=trim(tmp_file_name)//"/Basis_scratch")

      do
        read(1,'(A)',end=3) lines
        do i = 1 , number_of_atoms 
          if (lines == atom_type1(i) .or. lines == atom_type2(i) ) then 
            write(2,'(A)') trim(lines)
            do
              read (1,'(A)',end=3) lines
              if (adjustl(lines(1:1))=="A" .or. adjustl(lines(1:1))=="a") exit
              write(2,'(A)') trim(lines)
            end do 
          end if 
        end do 
      end do 

3      close(1)

close(2)

      ! --------------------------------------------------------------- !

end subroutine