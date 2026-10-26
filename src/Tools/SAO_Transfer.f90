subroutine transform_AO_to_SAO(nBas, AO_matrix, SAO_matrix)

      use constants_module
      use unitcell_module
      
      implicit none

      ! input
      integer, intent(in)       :: nBas
      real(dp), intent(in)      :: AO_matrix(nBas, nBas)

      ! output
      complex(dpc), intent(out) :: SAO_matrix(nBas, nBas)

      ! local
      integer                   :: kx, ky, kz 
      integer                   :: Rx, Ry, Rz
      integer                   :: a, b, j
      integer                   :: mu_sao, nu_sao
      integer                   :: k_index, cell_index
      real(dp)                  :: theta
      complex(dpc)              :: val

      integer     , allocatable :: coord(:,:)
      complex(dpc), allocatable :: phase(:,:)
      integer                   :: ix, iy, iz, idx
      integer                   :: k1 , r2

      SAO_matrix = (0.0_dp, 0.0_dp)


      allocate(coord(3,N_cell))
      allocate(phase(N_cell,N_cell))

      ! Zero-based coordinates, x fastest.
      ! The same indexing is used for cells and k points.

      do iz = 0, nz-1
        do iy = 0, ny-1
          do ix = 0, nx-1
            idx = iz*nx*ny + iy*nx + ix + 1
            coord(:,idx) = [ix,iy,iz]
          end do
        end do
      end do

      ! phase(k,R) = exp(+i k.R)

      do k1 = 1, N_cell
        do r2 = 1, N_cell
          theta = 2.0_dp*pi * ( &
            real(coord(1,k1),dp)*real(coord(1,r2),dp)/real(nx,dp) &
          + real(coord(2,k1),dp)*real(coord(2,r2),dp)/real(ny,dp) &
          + real(coord(3,k1),dp)*real(coord(3,r2),dp)/real(nz,dp))
          theta = theta - 2.0_dp*pi * anint(theta/(2.0_dp*pi))
          phase(k1,r2) = exp(cmplx(0.0_dp,theta,kind=dpc))
        end do
      end do

      ! loop over all k-points (reciprocal cells)

      do kz = 1, nz
        do ky = 1, ny
          do kx = 1, nx

              ! linear index of this k-point (1-based, x fastest)
              k_index = (kz-1)*nx*ny + (ky-1)*nx + kx

              do a = 1, nfuc
                mu_sao = (k_index-1)*nfuc + a

                do b = 1, nfuc
                  nu_sao = (k_index-1)*nfuc + b

                  val    = (0.0_dp, 0.0_dp)

                  ! sum over all real-space cells
                  do Rz = 1, nz
                    do Ry = 1, ny
                      do Rx = 1, nx

                        ! linear cell index (1-based, x fastest)
                        cell_index = (Rz-1)*nx*ny + (Ry-1)*nx + Rx

                        ! column index in AO_matrix for basis function b in cell (Rx,Ry,Rz)
                        j = b + (cell_index-1)*nfuc

                        val = val + phase(k_index,cell_index) * AO_matrix(a, j)

                      end do
                    end do
                  end do

                  SAO_matrix(mu_sao, nu_sao) = val

                end do
              end do

          end do
        end do
      end do

      !call clean_and_symmetrize(nBas, SAO_matrix)

end subroutine transform_AO_to_SAO


subroutine transform_SAO_to_AO(nBas, M_SAO, M_AO)

      use constants_module
      use unitcell_module
      implicit none

      ! input
      integer, intent(in)            :: nBas
      complex(dpc), intent(in)       :: M_SAO(nBas, nBas)

      ! output
      double precision, intent(out)  :: M_AO(nBas, nBas)

      ! local
      integer                        :: kx, ky, kz
      integer                        :: Rx, Ry, Rz, Rpx, Rpy, Rpz
      integer                        :: a, b, mu, nu, mu_sao, nu_sao
      integer                        :: k_index, cell_index, cellp_index
      real(dp)                       :: theta
      complex(dpc)                   :: val
      integer                        :: Nc
      complex(dpc)                   :: M_AO_temp(nBas, nBas)

      integer     , allocatable      :: coord(:,:)
      complex(dpc), allocatable      :: phase(:,:)
      integer                        :: ix, iy, iz, idx
      integer                        :: k1 , r2
      integer                        :: dRx, dRy, dRz, dR_index



      Nc        = nx * ny * nz
      M_AO_temp = (0.0_dp, 0.0_dp)
      M_AO      = 0.0_dp

      allocate(coord(3,N_cell))
      allocate(phase(N_cell,N_cell))


      ! Zero-based coordinates, x fastest.
      ! The same indexing is used for cells and k points.

      do iz = 0, nz-1
        do iy = 0, ny-1
          do ix = 0, nx-1
            idx = iz*nx*ny + iy*nx + ix + 1
            coord(:,idx) = [ix,iy,iz]
          end do
        end do
      end do

      ! phase(k,R) = exp(+i k.R)

      do k1 = 1, N_cell
        do r2 = 1, N_cell
          theta = 2.0_dp*pi * ( &
            real(coord(1,k1),dp)*real(coord(1,r2),dp)/real(nx,dp) &
          + real(coord(2,k1),dp)*real(coord(2,r2),dp)/real(ny,dp) &
          + real(coord(3,k1),dp)*real(coord(3,r2),dp)/real(nz,dp))
          theta = theta - 2.0_dp*pi * anint(theta/(2.0_dp*pi))
          phase(k1,r2) = exp(cmplx(0.0_dp,theta,kind=dpc))
        end do
      end do


      ! loop over all k-points (reciprocal cells)
      do kz = 1, nz
        do ky = 1, ny
          do kx = 1, nx

            ! linear index of this k-point (1-based, x fastest)
            k_index = (kz-1)*nx*ny + (ky-1)*nx + kx

            do a = 1, nfuc
              mu_sao = (k_index-1)*nfuc + a

              do b = 1, nfuc
                nu_sao = (k_index-1)*nfuc + b   ! block-diagonal: same k_index
                val = M_SAO(mu_sao, nu_sao)

                if (abs(val) < 1.0e-15_dp) cycle   ! skip zero blocks

                ! sum over real-space cells R (for orbital a) and R' (for orbital b)
                do Rz = 1, nz
                  do Ry = 1, ny
                    do Rx = 1, nx

                      cell_index = (Rz-1)*nx*ny + (Ry-1)*nx + Rx
                      mu = a + (cell_index-1)*nfuc

                      do Rpz = 1, nz
                        do Rpy = 1, ny
                          do Rpx = 1, nx

                            cellp_index = (Rpz-1)*nx*ny + (Rpy-1)*nx + Rpx
                            nu = b + (cellp_index-1)*nfuc

                            dRx = mod(Rx - Rpx + nx, nx)
                            dRy = mod(Ry - Rpy + ny, ny)
                            dRz = mod(Rz - Rpz + nz, nz)
                            dR_index = dRz*nx*ny + dRy*nx + dRx + 1

                            M_AO_temp(mu,nu) = M_AO_temp(mu,nu) + phase(k_index, dR_index) * val

                          end do
                        end do
                      end do

                    end do
                  end do
                end do

              end do
            end do

          end do
        end do
      end do

      M_AO_temp = M_AO_temp / dble(Nc)   ! inverse Fourier transform normalisation

      ! call clean_and_symmetrize(nBas, M_AO_temp)

      M_AO(:,:) = dble(M_AO_temp(:,:))

end subroutine transform_SAO_to_AO

subroutine clean_and_symmetrize(N, A)

      use constants_module

      implicit none

      integer, intent(in)         :: N
      complex(dpc), intent(inout) :: A(N,N)
      integer                     :: i, j
      real(dp), parameter         :: tol = 1.0d-15

      do i = 1, N
        do j = 1, i-1
          ! Zero out tiny elements (both real and imag)
          if (abs(A(i,j)) < tol) A(i,j) = (0.0_dp, 0.0_dp)
          ! Enforce Hermitian: upper = conjg(lower)
          A(j,i) = conjg(A(i,j))
        end do
        ! Diagonal must be real
        A(i,i) = dble(A(i,i))
      end do
      
end subroutine clean_and_symmetrize


subroutine sort_array(arr, n)

      use constants_module
      implicit none
      integer, intent(in)     :: n
      real(dp), intent(inout) :: arr(n)
      integer                 :: i, j
      real(dp)                :: tmp

      do i = 1, n-1
        do j = i+1, n
          if (arr(i) > arr(j)) then
            tmp = arr(i)
            arr(i) = arr(j)
            arr(j) = tmp
          end if
        end do
      end do

end subroutine sort_array