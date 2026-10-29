subroutine transform_ERI_AO_to_SAO(nBas, ERI_AO, ERI_SAO)

      use constants_module, only: dp, dpc, pi
      use unitcell_module , only: nfuc, nx, ny, nz, N_cell

      implicit none

      integer       , intent(in)  :: nBas
      real(dp)      , intent(in)  :: ERI_AO(nfuc,nBas,nBas,nBas)
      complex(dpc)  , intent(out) :: ERI_SAO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      integer                     :: c, d, j3, j4, lo, hi
      integer                     :: k1, k2, k3, k4
      integer                     :: r2, r3, r4
      integer                     :: ix, iy, iz
      integer                     :: k4x, k4y, k4z
      integer                     :: idx

      integer     , allocatable   :: coord(:,:)
      complex(dpc), allocatable   :: phase(:,:)
      complex(dpc), allocatable   :: X1(:,:,:,:,:)     ! (a, j2, j3, d, k4)
      complex(dpc), allocatable   :: X2(:,:,:,:,:,:)   ! (a, j2, c, d, k3, k4)

      complex(dpc)                :: ph
      real(dp)                    :: theta, invN

      if (min(nx,ny,nz) < 1)       error stop 'transform_ERI_AO_to_SAO: invalid cell grid'
      if (N_cell /= nx*ny*nz)      error stop 'transform_ERI_AO_to_SAO: inconsistent N_cell'
      if (nBas /= nfuc*N_cell)     error stop 'transform_ERI_AO_to_SAO: nBas /= nfuc*N_cell'

      allocate(coord(3,N_cell), phase(N_cell,N_cell))

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
          phase(k1,r2) = exp(cmplx(0.0_dp,theta,kind=dpc))
        end do
      end do
      invN = 1.0_dp/real(N_cell,dp)

      ! ---- Step 1: R4 -> k4   X1(a,j2,j3,d,k4) = sum_R4 ERI_AO(a,j2,j3,d R4) e^{+i k4 R4}

      allocate(X1(nfuc,nBas,nBas,nfuc,N_cell))
      !$omp parallel do collapse(2) schedule(static) default(shared) private(k4,d,r4,j4,ph)
      do k4 = 1, N_cell
        do d = 1, nfuc
          X1(:,:,:,d,k4) = (0.0_dp,0.0_dp)
          do r4 = 1, N_cell
            j4 = d + (r4-1)*nfuc
            ph = phase(k4,r4)
            X1(:,:,:,d,k4) = X1(:,:,:,d,k4) + ERI_AO(:,:,:,j4)*ph
          end do
        end do
      end do
      !$omp end parallel do

      ! ---- Step 2: R3 -> k3   X2(a,j2,c,d,k3,k4) = sum_R3 X1(a,j2,c R3,d,k4) e^{-i k3 R3}
      allocate(X2(nfuc,nBas,nfuc,nfuc,N_cell,N_cell))
      !$omp parallel do collapse(3) schedule(static) default(shared) private(k4,k3,d,c,r3,j3,ph)
      do k4 = 1, N_cell
        do k3 = 1, N_cell
          do d = 1, nfuc
            do c = 1, nfuc
              X2(:,:,c,d,k3,k4) = (0.0_dp,0.0_dp)
              do r3 = 1, N_cell
                j3 = c + (r3-1)*nfuc
                ph = conjg(phase(k3,r3))
                X2(:,:,c,d,k3,k4) = X2(:,:,c,d,k3,k4) + X1(:,:,j3,d,k4)*ph
              end do
            end do
          end do
        end do
      end do
      !$omp end parallel do
      deallocate(X1)

      ! ---- Step 3: R2 -> k2, pick k4 = k1 - k2 + k3, apply 1/N
      !$omp parallel do collapse(3) schedule(static) default(shared) &
      !$omp   private(k1,k2,k3,k4,k4x,k4y,k4z,d,c,r2,lo,hi,ph)
      do k1 = 1, N_cell
        do k2 = 1, N_cell
          do k3 = 1, N_cell
            k4x = modulo(coord(1,k1)-coord(1,k2)+coord(1,k3),nx)
            k4y = modulo(coord(2,k1)-coord(2,k2)+coord(2,k3),ny)
            k4z = modulo(coord(3,k1)-coord(3,k2)+coord(3,k3),nz)
            k4  = k4z*nx*ny + k4y*nx + k4x + 1
            do d = 1, nfuc
              do c = 1, nfuc
                ERI_SAO(:,:,c,d,k1,k2,k3) = (0.0_dp,0.0_dp)
                do r2 = 1, N_cell
                  lo = (r2-1)*nfuc + 1
                  hi = r2*nfuc
                  ph = phase(k2,r2)*invN
                  ERI_SAO(:,:,c,d,k1,k2,k3) = ERI_SAO(:,:,c,d,k1,k2,k3) &
                                            + X2(:,lo:hi,c,d,k3,k4)*ph
                end do
              end do
            end do
          end do
        end do
      end do
      !$omp end parallel do

      deallocate(X2, coord, phase)

end subroutine transform_ERI_AO_to_SAO


subroutine transform_ERI_SAO_to_kMO(nBas, coeff_SAO, ERI_SAO, ERI_kMO)

      use constants_module, only: dp, dpc
      use unitcell_module , only: nfuc, nx, ny, nz, N_cell

      implicit none

      integer     , intent(in)  :: nBas
      complex(dpc), intent(in)  :: coeff_SAO(nBas,nBas)
      complex(dpc), intent(in)  :: ERI_SAO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: ERI_kMO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      integer                   :: k1, k2, k3, k4
      integer                   :: ix, iy, iz
      integer                   :: k4x, k4y, k4z
      integer                   :: idx 
      integer                   :: off1, off2, off3, off4
      integer                   :: a, b, c, d
      integer                   :: p, q, r, s
      integer                   :: row, col

      integer     , allocatable :: coord(:,:)
      complex(dpc), allocatable :: scr1(:,:,:,:), scr2(:,:,:,:), scr3(:,:,:,:)


      real(dp)                  :: max_offblock, coeff_scale
      real(dp), parameter       :: block_tol = 1.0e-10_dp

      if (min(nx,ny,nz) < 1)   error stop 'transform_ERI_SAO_to_kMO: invalid grid'
      if (N_cell /= nx*ny*nz)  error stop 'transform_ERI_SAO_to_kMO: inconsistent N_cell'
      if (nBas /= nfuc*N_cell) error stop 'transform_ERI_SAO_to_kMO: inconsistent nBas'

      max_offblock = 0.0_dp
      coeff_scale  = maxval(abs(coeff_SAO))
      do col = 1, nBas
        do row = 1, nBas
          if ((row-1)/nfuc /= (col-1)/nfuc) max_offblock = max(max_offblock,abs(coeff_SAO(row,col)))
        end do
      end do
      if (max_offblock > block_tol*max(1.0_dp,coeff_scale)) then
        write(*,'(A,ES16.8)') 'Maximum off-k-block coefficient: ', max_offblock
        error stop 'transform_ERI_SAO_to_kMO: coefficients not k-block ordered'
      end if

      allocate(coord(3,N_cell))
      do iz = 0, nz-1
        do iy = 0, ny-1
          do ix = 0, nx-1
            idx = iz*nx*ny + iy*nx + ix + 1
            coord(:,idx) = [ix,iy,iz]
          end do
        end do
      end do

      !$omp parallel default(shared) &
      !$omp   private(k1,k2,k3,k4,k4x,k4y,k4z,off1,off2,off3,off4, &
      !$omp           a,b,c,d,p,q,r,s,scr1,scr2,scr3)
      allocate(scr1(nfuc,nfuc,nfuc,nfuc), scr2(nfuc,nfuc,nfuc,nfuc), scr3(nfuc,nfuc,nfuc,nfuc))

      !$omp do collapse(3) schedule(static)
      do k1 = 1, N_cell
        do k2 = 1, N_cell
          do k3 = 1, N_cell
            off1 = (k1-1)*nfuc
            off2 = (k2-1)*nfuc
            off3 = (k3-1)*nfuc
            k4x = modulo(coord(1,k1)-coord(1,k2)+coord(1,k3),nx)
            k4y = modulo(coord(2,k1)-coord(2,k2)+coord(2,k3),ny)
            k4z = modulo(coord(3,k1)-coord(3,k2)+coord(3,k3),nz)
            k4  = k4z*nx*ny + k4y*nx + k4x + 1
            off4 = (k4-1)*nfuc

            ! 4th index: d -> s
            scr1 = (0.0_dp,0.0_dp)
            do s = 1, nfuc
              do d = 1, nfuc
                do c = 1, nfuc
                  do b = 1, nfuc
                    do a = 1, nfuc
                      scr1(a,b,c,s) = scr1(a,b,c,s) + ERI_SAO(a,b,c,d,k1,k2,k3)*coeff_SAO(off4+d,off4+s)
                    end do
                  end do
                end do
              end do
            end do
            ! 3rd index: c -> r (conjugated)
            scr2 = (0.0_dp,0.0_dp)
            do s = 1, nfuc
              do r = 1, nfuc
                do c = 1, nfuc
                  do b = 1, nfuc
                    do a = 1, nfuc
                      scr2(a,b,r,s) = scr2(a,b,r,s) + scr1(a,b,c,s)*conjg(coeff_SAO(off3+c,off3+r))
                    end do
                  end do
                end do
              end do
            end do
            ! 2nd index: b -> q
            scr3 = (0.0_dp,0.0_dp)
            do s = 1, nfuc
              do r = 1, nfuc
                do q = 1, nfuc
                  do b = 1, nfuc
                    do a = 1, nfuc
                      scr3(a,q,r,s) = scr3(a,q,r,s) + scr2(a,b,r,s)*coeff_SAO(off2+b,off2+q)
                    end do
                  end do
                end do
              end do
            end do
            ! 1st index: a -> p (conjugated); each thread owns this slab
            ERI_kMO(:,:,:,:,k1,k2,k3) = (0.0_dp,0.0_dp)
            do s = 1, nfuc
              do r = 1, nfuc
                do q = 1, nfuc
                  do p = 1, nfuc
                    do a = 1, nfuc
                      ERI_kMO(p,q,r,s,k1,k2,k3) = ERI_kMO(p,q,r,s,k1,k2,k3) &
                        + conjg(coeff_SAO(off1+a,off1+p))*scr3(a,q,r,s)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      !$omp end do

      deallocate(scr1, scr2, scr3)
      !$omp end parallel
      deallocate(coord)

end subroutine transform_ERI_SAO_to_kMO
