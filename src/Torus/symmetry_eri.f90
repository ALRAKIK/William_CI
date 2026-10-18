subroutine shift_integrals(known_integrals, full_integrals, nbasis,nbasis_unitcell)
      
      IMPLICIT NONE
  
      ! Input/Output declarations

      integer          , intent(in)  :: nbasis                                        ! Number of basis functions
      integer          , intent(in)  :: nbasis_unitcell                               ! Number of basis functions_per_unitcell
      double precision , intent(in)  :: known_integrals(nbasis,nbasis,nbasis,nbasis)  ! Integrals involving index 1
      double precision , intent(out) :: full_integrals(nbasis,nbasis,nbasis,nbasis)   ! Complete set of integrals
      

      ! Local variables

      integer                        :: index_i , index_j , index_k , index_l , shift 
      integer                        :: i, j, k, l , atom_i 


      full_integrals = 0.0D0
  
      ! First, copy all known integrals to the full array

      do i = 1, nbasis
        do j = 1, nbasis
          do k = 1, nbasis
            do l = 1, nbasis
              if (dabs(known_integrals(i,j,k,l)) >  1.d-30) full_integrals(i,j,k,l) = known_integrals(i,j,k,l)
            end do 
          end do 
        end do 
      end do 

      ! Process all orbitals beyond atom 1

      do i = nbasis_unitcell + 1 , nbasis

          ! Calculate which atom this orbital belongs to

          atom_i = (i-1)/nbasis_unitcell + 1

          ! Calculate the proper shift based on atom number

          shift = (atom_i - 1) * nbasis_unitcell

          do j = 1, nbasis
              do k = 1, nbasis
                  do l = 1, nbasis
                      ! Apply the correct shift with proper wrapping
                      index_i = i - shift
                      if (index_i <= 0) index_i = index_i + nbasis

                      index_j = j - shift
                      if (index_j <= 0) index_j = index_j + nbasis

                      index_k = k - shift
                      if (index_k <= 0) index_k = index_k + nbasis

                      index_l = l - shift
                      if (index_l <= 0) index_l = index_l + nbasis

                      full_integrals(i,j,k,l) = known_integrals(index_i, index_j, index_k, index_l)
                  end do
              end do
          end do
      end do
  
end subroutine shift_integrals

