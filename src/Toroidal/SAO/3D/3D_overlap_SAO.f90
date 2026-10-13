subroutine overlap_matrix_toroidal_SAO_3D(nBas,AO,overlap)

      use files
      use atom_basis
      use torus_init
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

      double precision,intent(out) :: overlap(nBas,nBas)

      !-----------------------------------------------------------------!

      overlap(:,:) = 0.d0
      
      ! --------------------------------------------------------------- !

      do i = 1 , nfuc
        do j = 1 , nBas

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z


          pattern_id = encode_orbital_pattern_AO(AO1%orbital, AO2%orbital)

          call overlap_integral_toroidal_3D(pattern_id,r1,r2,AO1,AO2,overlap(i,j))
          
        end do 
      end do 

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
       
      do i = 1 , nfuc
        do j = 1 , nBas
          if (abs(overlap(i,j)) < 1e-14) overlap(i,j) = 0.d0 
        end do 
      end do 

      call fill_with_symmetry(nBas, overlap)
 
      if (c_OV) then 
        write(outfile,'(a)') ""
        write(outfile,'(a)') "------------------"
        write(outfile,'(a)') "The overlap matrix"
        write(outfile,'(a)') "------------------"
        write(outfile,'(a)') ""
        do i = 1 , nBas - 1 
          do j = i , nBas
            write(outfile,'(i3,i3, 6x,1000(f16.12,2x))') i , j, overlap(i,j)
          end do
        end do 
      end if 

end subroutine overlap_matrix_toroidal_SAO_3D