subroutine check_symmetric_matrix(nBas,S,T,V,HC)

      implicit none 

      integer, intent(in) :: nBas
      integer             :: i , j

      double precision, intent(in) :: S (nBas,nBas)
      double precision, intent(in) :: T (nBas,nBas)
      double precision, intent(in) :: V (nBas,nBas)
      double precision, intent(in) :: HC(nBas,nBas) 


      do i = 1 , nBas-1
        do j = i , nBas 
          if (abs(S(i,j)-S(j,i)) > 1.0d-10) then 
            print * , "S is not symmetric"
            stop 
          end if 
          if (abs(T(i,j)-T(j,i)) > 1.0d-10) then 
            print * , "T is not symmetric"
            stop 
          end if 
          if (abs(V(i,j)-V(j,i)) > 1.0d-10) then 
            print * , "V is not symmetric"
            stop 
          end if 
          if (abs(HC(i,j)-HC(j,i)) > 1.0d-10) then 
            print * , "HC is not symmetric"
            stop 
          end if
        end do 
      end do 

end subroutine check_symmetric_matrix

subroutine fill_with_symmetry(nBas, matrix)

      use constants_module
      use unitcell_module
      implicit none

      integer, intent(in)     :: nBas
      real(dp), intent(inout) :: matrix(nBas, nBas)

      integer :: i, j, a, b, Ri, Rj
      integer :: Rxi, Ryi, Rzi, Rxj, Ryj, Rzj
      integer :: Rdx, Rdy, Rdz, src_cell, j_src

      ! Loop over all rows beyond the first unit cell
      do i = nfuc + 1, nBas

        ! Decode row index i into (cell index Ri, orbital a)
        Ri = (i-1)/nfuc + 1
        a  = mod(i-1, nfuc) + 1

        ! Decode cell Ri into (Rxi, Ryi, Rzi) — x fastest
        Rzi = (Ri-1)/(nx*ny) + 1
        Ryi = mod((Ri-1)/nx, ny) + 1
        Rxi = mod(Ri-1, nx) + 1

        do j = 1, nBas

          ! Decode column index j into (cell index Rj, orbital b)
          Rj = (j-1)/nfuc + 1
          b  = mod(j-1, nfuc) + 1

          ! Decode cell Rj into (Rxj, Ryj, Rzj)
          Rzj = (Rj-1)/(nx*ny) + 1
          Ryj = mod((Rj-1)/nx, ny) + 1
          Rxj = mod(Rj-1, nx) + 1

          ! Relative cell shift (modular, 0-based)
          Rdx = mod(Rxj - Rxi + nx, nx)
          Rdy = mod(Ryj - Ryi + ny, ny)
          Rdz = mod(Rzj - Rzi + nz, nz)

          ! Corresponding cell in the reference row (cell 0)
          src_cell = Rdz*nx*ny + Rdy*nx + Rdx + 1
          j_src = b + (src_cell-1)*nfuc

          ! Copy the value from the first nfuc rows
          matrix(i,j) = matrix(a, j_src)

        end do
      end do

end subroutine fill_with_symmetry