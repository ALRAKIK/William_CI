subroutine get_g3space_spin_adapted(nBas, nOcc, nSpace, g03space)
   
      ! Input variables
      integer, intent(in)           :: nBas, nOcc,  nSpace

      ! Output 
      integer, intent(out)          :: g03space(:,:)
         
      ! Local variables
      integer                       :: i, j, l, idx
      double precision              :: fi, fj, fl
      integer, allocatable          :: tempSpace(:,:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '         |  Auxiliary Subspace Creation  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      ! Temporary workspace (maximum possible size)
      allocate(tempSpace(3, nSpace))
         
      ! Selection algorithm
      idx = 1
      do i = 1, nBas
         fi = 1.0d0

         if (i > nOcc) fi = 0.0d0
            
         do j = 1, nBas

            fj = 1.0d0
            if (j > nOcc) fj = 0.0d0
                  
            if (i >= j) then
               do l = 1, nBas
                  fl = 1.0d0
                        
                  if (l > nOcc) fl = 0.0d0

                  ! Makes sure it's either 2e1h or 2h1e
                  if ((fi - fl) * (fj - fl) /= 0.0d0) then
                     tempSpace(1, idx) = i
                     tempSpace(2, idx) = j
                     tempSpace(3, idx) = l
                     idx = idx + 1
                  end if
               end do
            end if
         end do
      end do 

      ! Copy to output list
      g03space = tempSpace(:, 1:nSpace)

      deallocate(tempSpace)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '   |  Auxiliary Subspace Spin Adapted Finished  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

end subroutine get_g3space_spin_adapted


