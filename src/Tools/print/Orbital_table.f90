
!---------------------------------------------------------------------!
! Orbital table per SHELL: s, p, d, f, g (the Cartesian components    !
! px, py, ... are listed only as an index range)                      !
!---------------------------------------------------------------------!

subroutine print_shell_table(AO, nBas)

      use files
      use classification_ERI
      use keywords

      implicit none

      integer           ,intent(in) :: nBas
      type(ERI_function),intent(in) :: AO(nBas)

      integer :: i, j, ncomp, nsh, u, pass
      character(len=1) :: shl

      call header_under("Atomic Orbitals", -1)
      do pass = 1 , 2
        if (pass == 1) then
          u = outfile
        else
          if (.not. c_Orbitals) exit
          u = orbfile
        end if
        write(u,'(a)') "|"//repeat("-",96)//"|"
        write(u,'(a1,2x,a6,2x,a8,9x,a1,11x,a1,11x,a1,6x,a4,7x,a8,7x,a11,a2)') "|", "Shell", "AO range ", "X", "Y", "Z", "Type", &
                                                          "Exponent", "Coefficient", " |"
        write(u,'(a)') "|"//repeat("-",96)//"|"
        i = 1 ; nsh = 0
        do while (i <= nBas)
          shl = AO(i)%orbital(1:1)
          select case (shl)
          case ('s') ; ncomp = 1
          case ('p') ; ncomp = 3
          case ('d') ; ncomp = 6
          case ('f') ; ncomp = 10
          case ('g') ; ncomp = 15
          case default ; ncomp = 1
          end select
          nsh = nsh + 1
          write(u,'(a1,i5,2x,i5,a1,i5,1x,3f12.6,4x,a1,4x,f15.8,f15.8,a3)') "|", nsh, i, " ", i+ncomp-1, &
               AO(i)%x, AO(i)%y, AO(i)%z, shl, AO(i)%exponent(1), AO(i)%coefficient(1), "  |"
          do j = 2 , size(AO(i)%exponent)
            write(u,'(a1,64x,f15.8,f15.8,a3)') "|", AO(i)%exponent(j), AO(i)%coefficient(j), "  |"
          end do
          write(u,'(a)') "|"//repeat("-",96)//"|"
          i = i + ncomp
        end do
        write(u,'(A)') ""
        flush(u)
      end do

end subroutine print_shell_table
