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

subroutine get_g4space(nBas, nOcc, nSpace, g04space)
   
      ! Input variables
      integer, intent(in)           :: nBas, nOcc,  nSpace

      ! Output 
      integer, intent(out)          :: g04space(:,:)
         
      ! Local variables
      integer                       :: i, j, l, n, idx, nBasZ, nOccZ
      double precision              :: fi, fj, fl, fn
      integer, allocatable          :: tempSpace(:,:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '         |  Auxiliary Subspace Creation  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      ! Temporary workspace (maximum possible size)
      allocate(tempSpace(4, nSpace))

      nOccZ = 2 * nOcc
      nBasZ = 2 * nBas
         
      ! Selection algorithm
      idx = 1
      do i = 1, nBasZ
         fi = 1.0d0
         ! nO times 2 because we use the spin orbital basis set
         if (i > nOccZ) fi = 0.0d0
            
         do j = 1, nBasZ

            fj = 1.0d0
            if (j > nOccZ) fj = 0.0d0
                  
            if (i > j) then
               do l = 1, nBasZ
                  fl = 1.0d0
                        
                  if (l > nOccZ) fl = 0.0d0

                  do n = 1, nBasZ

                     fn = 1.0d0
                     if (n > nOccZ) fn = 0.0d0

                     if (l > n) then

                        ! Makes sure it's 2e2h
                        if ((fi - fl) * (fi - fl)*(fj - fn) /= 0.0d0) then
                           tempSpace(1, idx) = i
                           tempSpace(2, idx) = j
                           tempSpace(3, idx) = l
                           tempSpace(4, idx) = n
                           idx = idx + 1
                        end if
                     end if
                  end do
               end do
            end if
         end do
      end do 

      ! Copy to output list
      g04space = tempSpace(:, 1:nSpace)

      deallocate(tempSpace)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '    |  Auxiliary Subspace Creation Finished  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)
   end subroutine get_g4space


subroutine get_g2space_4_2(nBas, nOcc, nSpace, g02space)
   
      ! Input variables
      integer, intent(in)           :: nBas, nOcc,  nSpace

      ! Output 
      integer, intent(out)          :: g02space(:,:)
         
      ! Local variables
      integer                       :: j, l, idx, nBasZ, nOccZ
      double precision              :: fj, fl
      integer, allocatable          :: tempSpace(:,:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '         |  Auxiliary Subspace Creation  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      ! Temporary workspace (maximum possible size)
      allocate(tempSpace(2, nSpace))

      nOccZ = 2 * nOcc
      nBasZ = 2 * nBas
      write(*,*)'nBas = ', nBasZ   
      ! Selection algorithm
      idx = 1
      do j = 1, nBasZ
         fj = 1.0d0
         ! nO times 2 because we use the spin orbital basis set
         if (j > nOccZ) fj = 0.0d0
            
         do l = 1, nBasZ

            fl = 1.0d0
            if (l > nOccZ) fl = 0.0d0

            if (j > l) then

                        ! Makes sure it's 2e
               if (fj+fl-1 /= 0.0d0) then
                  tempSpace(1, idx) = j
                  tempSpace(2, idx) = l
                  idx = idx + 1
               end if               
            end if
         end do
      end do


      ! Copy to output list
      g02space = tempSpace(:, 1:nSpace)

      deallocate(tempSpace)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '    |  Auxiliary Subspace Creation Finished  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

end subroutine get_g2space_4_2


subroutine get_g4space_4_2(nBas, nOcc, nSpace, g04space)
   
      ! Input variables
      integer, intent(in)           :: nBas, nOcc, nSpace

      ! Output 
      integer, intent(out)          :: g04space(:,:)
         
      ! Local variables
      integer                       :: i, j, l, n, idx, nBasZ, nOccZ
      double precision              :: fi, fj, fl, fn
      integer, allocatable          :: tempSpace(:,:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '         |  Auxiliary Subspace Creation  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      ! Temporary workspace (maximum possible size)
      allocate(tempSpace(4, nSpace))
      write(*,*)"nspace4 =", nspace

      nOccZ = 2 * nOcc
      nBasZ = 2 * nBas
         
      ! Selection algorithm
      idx = 1
      do i = 1, nBasZ
         fi = 1.0d0
         ! nO times 2 because we use the spin orbital basis set
         if (i > nOccZ) fi = 0.0d0
            
         do j = 1, nBasZ

            fj = 1.0d0
            if (j > nOccZ) fj = 0.0d0
                  
            if (i > j) then
               do l = 1, nBasZ
                  fl = 1.0d0
                        
                  if (l > nOccZ) fl = 0.0d0

                  do n = 1, nBasZ

                     fn = 1.0d0
                     if (n > nOccZ) fn = 0.0d0

                     if (j > l) then

                        ! Makes sure it's 3e1h
                        if ((fi - fn) * (fl - fn)*(fj - fn) /= 0.0d0) then
                           tempSpace(1, idx) = i
                           tempSpace(2, idx) = j
                           tempSpace(3, idx) = l
                           tempSpace(4, idx) = n
                           idx = idx + 1
                        end if
                     end if
                  end do
               end do
            end if
         end do
      end do 

      ! Copy to output list

      g04space = tempSpace(:, 1:nSpace)
      deallocate(tempSpace)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '    |  Auxiliary Subspace Creation Finished  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)
   end subroutine get_g4space_4_2