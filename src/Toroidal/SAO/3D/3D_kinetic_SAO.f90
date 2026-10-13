subroutine kinetic_matrix_toroidal_SAO_3D(nBas,AO,kinetic)

      use files 
      use torus_init
      use atom_basis
      use classification_ERI
      use keywords
      use unitcell_module

      implicit none 

      !-----------------------------------------------------------------!

      integer                      :: i , j
      integer                      :: nBas

      type(ERI_function)           :: AO (nBas)
      type(ERI_function)           :: AO1 , AO2


      double precision             :: r1(3) , r2(3)

      integer                      :: pattern_id, encode_orbital_pattern_AO

      ! output ! 

      double precision,intent(out) :: kinetic(nBas,nBas)

      !-----------------------------------------------------------------!

      kinetic(:,:) = 0.d0
      
      do i = 1 , nfuc
        do j = 1 , nBas

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z


          pattern_id = encode_orbital_pattern_AO(AO1%orbital, AO2%orbital)

          call Kinetic_integral_toroidal_3D(pattern_id,r1,r2,AO1,AO2,kinetic(i,j))
          
        end do 
      end do 

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      do i = 1 , nfuc
        do j = 1 , nBas
          if (abs(kinetic(i,j)) < 1e-14) kinetic(i,j) = 0.d0 
        end do 
      end do 

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



end subroutine kinetic_matrix_toroidal_SAO_3D