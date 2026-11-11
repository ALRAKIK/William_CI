subroutine print_orbital_table(ERI, nBas)
      
      use files
      use classification_ERI
      use keywords

      implicit none

      integer, intent(in)             :: nBas
      type(ERI_function) , intent(in) :: ERI(nBas)
      integer                         :: i, j

      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !

      call header_under("Atomic Orbitals", -1)
      ! Print header
      write(outfile,'(A)') "|-----------------------------------------------------------------------------------------|"
      write(outfile,'(A)') "|  Idx         X         Y        Z      Type           Exponent        Coefficient       |"
      write(outfile,'(A)') "|-----------------------------------------------------------------------------------------|"
      do i = 1, nBas
        write(outfile,'(a1,I5,3x,3f10.6,3x,a5,3x,f16.10,2x,f16.10,4x,a3)')  "|", i, ERI(i)%x, ERI(i)%y, ERI(i)%z, ERI(i)%orbital, ERI(i)%exponent(1), ERI(i)%coefficient(1), "  |"
        do j = 2, size(ERI(i)%exponent)
          write(outfile,'(a1,49x,f16.10,2x,f16.10,4x,a3)') "|", ERI(i)%exponent(j), ERI(i)%coefficient(j), "  |"
        end do
        write(outfile,'(A)') "|-----------------------------------------------------------------------------------------|"
      end do
      write(outfile,'(A)') ""
      FLUSH(outfile)

      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !

      if (c_Orbitals) then 

      write(orbfile,'(A)') "|-----------------------------------------------------------------------------------------|"
      write(orbfile,'(A)') "|                                     The atomic orbitals                                 |"
      write(orbfile,'(A)') "|-----------------------------------------------------------------------------------------|"

      write(orbfile,'(A)') "|-----------------------------------------------------------------------------------------|"
      write(orbfile,'(A)') "|  Idx         X         Y        Z      Type           Exponent        Coefficient       |"
      write(orbfile,'(A)') "|-----------------------------------------------------------------------------------------|"
      do i = 1, nBas
        write(orbfile,'(a1,I5,3x,3f10.6,3x,a5,3x,f16.10,2x,f16.10,4x,a3)')  "|", i, ERI(i)%x, ERI(i)%y, ERI(i)%z, ERI(i)%orbital, ERI(i)%exponent(1), ERI(i)%coefficient(1), "  |"
        do j = 2, size(ERI(i)%exponent)
          write(orbfile,'(a1,49x,f16.10,2x,f16.10,4x,a3)') "|", ERI(i)%exponent(j), ERI(i)%coefficient(j), "  |"
        end do
        write(orbfile,'(A)') "|-----------------------------------------------------------------------------------------|"
      end do
      write(orbfile,'(A)') ""
      FLUSH(orbfile)

      end if


end subroutine print_orbital_table