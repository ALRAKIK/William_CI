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
      integer                   :: kx, ky, kz, Rx, Ry, Rz
      integer                   :: a, b, mu_sao, nu_sao, j
      integer                   :: k_index, cell_index
      real(dp)                  :: theta
      complex(dpc)              :: val

      SAO_matrix = (0.0_dp, 0.0_dp)

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

                        ! phase factor exp(i k·R)
                        theta = 2.0_dp * pi * ( dble(kx-1)*dble(Rx-1)/dble(nx) &
                                              + dble(ky-1)*dble(Ry-1)/dble(ny) &
                                              + dble(kz-1)*dble(Rz-1)/dble(nz) )

                        val = val + exp(dcmplx(0.0_dp, theta)) * AO_matrix(a, j)

                      end do
                    end do
                  end do

                  SAO_matrix(mu_sao, nu_sao) = val

                end do
              end do

          end do
        end do
      end do

      call clean_and_symmetrize(nBas, SAO_matrix)

end subroutine transform_AO_to_SAO


subroutine transform_SAO_to_AO(nBas, M_SAO, M_AO)

      use constants_module
      use unitcell_module
      implicit none

      ! input
      integer, intent(in)       :: nBas
      complex(dpc), intent(in)  :: M_SAO(nBas, nBas)

      ! output
      complex(dpc), intent(out) :: M_AO(nBas, nBas)

      ! local
      integer                   :: kx, ky, kz
      integer                   :: Rx, Ry, Rz, Rpx, Rpy, Rpz
      integer                   :: a, b, mu, nu, mu_sao, nu_sao
      integer                   :: k_index, cell_index, cellp_index
      real(dp)                  :: theta, thetap
      complex(dpc)              :: phase, phasep, val
      integer                   :: Nc

      Nc = nx * ny * nz
      M_AO = (0.0_dp, 0.0_dp)

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

                if (abs(val) < 1.0e-14_dp) cycle   ! skip zero blocks

                ! sum over real-space cells R (for orbital a) and R' (for orbital b)
                do Rz = 1, nz
                  do Ry = 1, ny
                    do Rx = 1, nx

                      cell_index = (Rz-1)*nx*ny + (Ry-1)*nx + Rx
                      mu = a + (cell_index-1)*nfuc

                      theta = 2.0_dp * pi * ( dble(kx-1)*dble(Rx-1)/dble(nx) &
                                            + dble(ky-1)*dble(Ry-1)/dble(ny) &
                                            + dble(kz-1)*dble(Rz-1)/dble(nz) )
                      phase = exp(dcmplx(0.0_dp, theta))

                        do Rpz = 1, nz
                          do Rpy = 1, ny
                            do Rpx = 1, nx

                              cellp_index = (Rpz-1)*nx*ny + (Rpy-1)*nx + Rpx
                              nu = b + (cellp_index-1)*nfuc

                              thetap = 2.0_dp * pi * ( dble(kx-1)*dble(Rpx-1)/dble(nx) &
                                                     + dble(ky-1)*dble(Rpy-1)/dble(ny) &
                                                     + dble(kz-1)*dble(Rpz-1)/dble(nz) )
                              phasep = exp(dcmplx(0.0_dp, thetap))

                              M_AO(mu, nu) = M_AO(mu, nu) + phase * val * conjg(phasep)

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

      M_AO = M_AO / dble(Nc)   ! inverse Fourier transform normalisation

      call clean_and_symmetrize(nBas, M_AO)

end subroutine transform_SAO_to_AO

subroutine clean_and_symmetrize(N, A)

      use constants_module

      implicit none

      integer, intent(in)         :: N
      complex(dpc), intent(inout) :: A(N,N)
      integer                     :: i, j
      real(dp), parameter         :: tol = 1.0d-14

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