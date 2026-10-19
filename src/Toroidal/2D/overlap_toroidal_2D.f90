subroutine overlap_matrix_toroidal_2D(number_of_atoms,number_of_functions,atoms,AO)

      use files 
      use atom_basis
      use torus_init
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      integer                      ::           i , j , k , l 
      integer                      :: index_atom1 , index_sym
      integer                      ::          index_unitcell
      integer                      ::         number_of_atoms
      integer                      ::     number_of_functions

      type(atom)                   ::   atoms(number_of_atoms)

      type(ERI_function)           :: AO (number_of_functions)
      type(ERI_function)           ::                AO1 , AO2

      double precision,allocatable ::             overlap(:,:)
      double precision             ::            r1(3) , r2(3)


      double precision,parameter   :: pi = 3.14159265358979323846D00

      !-----------------------------------------------------------------!

      allocate(overlap(number_of_functions,number_of_functions))

      overlap(:,:) = 0.d0

      index_atom1 = atoms(1)%num_s_function + 3*atoms(1)%num_p_function
      
      index_unitcell = 0

      do i = 1 , number_of_atom_in_unitcell
        index_unitcell = index_unitcell  + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do

      index_sym   = 0

      do i = 1 , number_of_atoms/2 + 1
        index_sym = index_sym + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do

      ! --------------------------------------------------------------- !

      do i = 1 , index_unitcell
        do j = i , number_of_functions

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z
          
          if (AO1%orbital =="s" .and. AO2%orbital == "s") then
            
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call overlap_integral_ss_toroidal_2D(r1,r2,AO1,AO2,overlap(i,j))
              end do 
            end do 

          end if

        end do 
      end do 

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      do i = index_unitcell + 1   , number_of_functions
        do j = index_unitcell + 1 , number_of_functions
          overlap(i,j) = overlap(i-index_unitcell,j-index_unitcell)
        end do 
      end do 
 
      do i = 1 , number_of_functions - 1 
        do j = i , number_of_functions
          overlap(j,i) = overlap(i,j)
        end do 
      end do 

      !open(1,file="./tmp/OV.dat")
      open(1,file=trim(tmp_file_name)//"/OV.dat ")
        do i = 1 , size(overlap,1)
          do j = i , size(overlap,1)
            !write(1,'(I5,I5,f24.16)') i , j , overlap(i,j)
            write(1,*) i , j , overlap(i,j)
          end do 
        end do 
      close(1)

      !open(1,file="./tmp/OV_matrix.dat")
      open(1,file=trim(tmp_file_name)//"/OV_matrix.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(overlap,1))
      do i = 1 , size(overlap,1)
        write(1,'(i3,6x,1000(f16.12,2x))') i ,  (overlap(i,j),j=1,size(overlap,1))
      end do 
      close(1)

      deallocate(overlap)


end subroutine overlap_matrix_toroidal_2D

      !-----------------------------------------------------------------!