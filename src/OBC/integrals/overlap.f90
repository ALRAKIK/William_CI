subroutine overlap_matrix(number_of_atoms,number_of_functions,geometry,atoms,overlap)

      use files
      use atom_basis
      use keywords

      implicit none 

      !-----------------------------------------------------------------!

      integer                      :: i , j , k , l 
      integer                      :: number_of_functions
      integer                      :: func_index(1:number_of_atoms)
      integer                      :: idx_s  , idx_p 
      integer                      :: idx_p1 , idx_p2 
      integer                      :: number_of_atoms
      type(atom)                   :: atoms(number_of_atoms)
      type(atom)                   :: atom1 , atom2 

      double precision             :: geometry(number_of_atoms,3)

      double precision             :: overlap(number_of_functions,number_of_functions)
      double precision             :: r1(3) , r2(3)

      double precision,parameter   :: pi = dacos(-1.d0)


      double precision             :: SS 
      double precision             :: SP(3) , PS(3)
      double precision             :: PP(3,3)

      !-----------------------------------------------------------------!

      overlap(:,:) = 0.d0 

      func_index(1) = 1
      do i = 2, number_of_atoms
        func_index(i) = func_index(i-1) + atoms(i-1)%num_s_function + 3 * atoms(i-1)%num_p_function
      end do

      do i = 1 , number_of_atoms
        do j = i , number_of_atoms
          
          r1(1) = geometry(i,1) ; r2(1) = geometry(j,1)
          r1(2) = geometry(i,2) ; r2(2) = geometry(j,2)
          r1(3) = geometry(i,3) ; r2(3) = geometry(j,3)

          atom1 = atoms(i)
          atom2 = atoms(j)

          do k = 1 , atom1%num_s_function
            do l = 1 , atom2%num_s_function
              call overlap_integral_ss(r1,r2,atom1,atom2,k,l,SS)
              overlap(func_index(i) + k - 1, func_index(j) + l - 1) = SS    ! SS 
              overlap(func_index(j) + l - 1, func_index(i) + k - 1) = SS    ! SS 
            end do 
          end do 

          do k = 1 , atom1%num_s_function
            do l = 1 , atom2%num_p_function
              call overlap_integral_sp(r1,r2,atom1,atom2,k,l,SP)
              idx_s = func_index(i) + k - 1 
              idx_p = func_index(j) + atom2%num_s_function - 1 + (l-1) * 3 
              overlap(  idx_s, idx_p+1) = SP(1)                             !  SPx
              overlap(idx_p+1,   idx_s) = SP(1)
              overlap(  idx_s, idx_p+2) = SP(2)                             !  SPy
              overlap(idx_p+2,   idx_s) = SP(2)  
              overlap(  idx_s, idx_p+3) = SP(3)                             !  SPz
              overlap(idx_p+3,   idx_s) = SP(3)  
            end do 
          end do 

          do k = 1 , atom1%num_p_function
            do l = 1 , atom2%num_s_function
              call overlap_integral_ps(r1,r2,atom1,atom2,k,l,PS)
              idx_p = func_index(i) + atom1%num_s_function - 1 + (k-1) * 3  
              idx_s = func_index(j) + l - 1                                 
              overlap(idx_p+1,   idx_s) = PS(1)                               ! PxS
              overlap(  idx_s, idx_p+1) = PS(1)                               
              overlap(idx_p+2,   idx_s) = PS(2)                               ! PyS
              overlap(  idx_s, idx_p+2) = PS(2)                               
              overlap(idx_p+3,   idx_s) = PS(3)                               ! PzS
              overlap(  idx_s, idx_p+3) = PS(3)                               
            end do 
          end do

          do k = 1 , atom1%num_p_function
            do l = 1 , atom2%num_p_function
              call overlap_integral_pp(r1,r2,atom1,atom2,k,l,PP)
              idx_p1 = func_index(i) + atom1%num_s_function - 1 + (k-1) * 3  ! p-function index in atom1
              idx_p2 = func_index(j) + atom2%num_s_function - 1 + (l-1) * 3  ! p-function index in atom2
              
              overlap(idx_p1+1, idx_p2+1) = PP(1,1)                            ! PxPx
              overlap(idx_p2+1, idx_p1+1) = PP(1,1)                            ! PxPx
              overlap(idx_p1+1, idx_p2+2) = PP(1,2)                            ! PxPy
              overlap(idx_p2+1, idx_p1+2) = PP(1,2)                            ! PxPy
              overlap(idx_p1+1, idx_p2+3) = PP(1,3)                            ! PxPz
              overlap(idx_p2+1, idx_p1+3) = PP(1,3)                            ! PxPz

              overlap(idx_p1+2, idx_p2+1) = PP(2,1)                            ! PyPx
              overlap(idx_p2+2, idx_p1+1) = PP(2,1)                            ! PyPx
              overlap(idx_p1+2, idx_p2+2) = PP(2,2)                            ! PyPy
              overlap(idx_p2+2, idx_p1+2) = PP(2,2)                            ! PyPy
              overlap(idx_p1+2, idx_p2+3) = PP(2,3)                            ! PyPz
              overlap(idx_p2+2, idx_p1+3) = PP(2,3)                            ! PyPz

              overlap(idx_p1+3, idx_p2+1) = PP(3,1)                            ! PzPx
              overlap(idx_p2+3, idx_p1+1) = PP(3,1)                            ! PzPx
              overlap(idx_p1+3, idx_p2+2) = PP(3,2)                            ! PzPy
              overlap(idx_p2+3, idx_p1+2) = PP(3,2)                            ! PzPy
              overlap(idx_p1+3, idx_p2+3) = PP(3,3)                            ! PzPz
              overlap(idx_p2+3, idx_p1+3) = PP(3,3)                            ! PzPz

            end do 
          end do 

        end do 
      end do 

      if (c_details) then 
      open(1,file=trim(tmp_file_name)//"/OV.dat")
        do i = 1 , size(overlap,1)
          do j = i , size(overlap,1)
            if (abs(overlap(i,j)) > 1e-15 ) write(1,*) i , j , overlap(i,j)
          end do 
        end do 
      close(1)
      end if 

end subroutine
