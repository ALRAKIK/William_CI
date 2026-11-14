!=====================================================================!
!                                                                     !
!     SAO path : real-space AO integrals  ->  k-space MO integrals    !
!                                                                     !
!   ERI_AO(a, j2, j3, j4)       a in the unit cell, j = band + cell   !
!   ERI_out(p,q,r,s, k1,k2,k3)  = (p k1, q k2 | r k3, s k4)           !
!                                 k4 = k1 - k2 + k3                   !
!                                                                     !
!   Done k4 by k4 (in small batches), never storing ERI_SAO:          !
!     1) R4 -> k4 :  X1(a,j2,j3,d)    = sum_R4 ERI_AO(a,j2,j3,d R4)   !
!                                        e^{+i k4 R4}                 !
!     2) R3 -> k3 :  X2(a,j2,c,d,k3)  = sum_R3 X1(a,j2,c R3,d)        !
!                                        e^{-i k3 R3}                 !
!     3) for every (k1,k2) with this k4 (k3 = k4 - k1 + k2):           !
!        R2 -> k2 :  S(a,b,c,d)       = 1/N sum_R2 X2(a,b R2,c,d,k3)  !
!                                        e^{+i k2 R2}                 !
!        4 quarter transforms with the k-block coefficients            !
!        -> ERI_out(:,:,:,:,k1,k2,k3)                                 !
!   Every element is summed in the same order as the former           !
!   transform_ERI_AO_to_SAO + transform_ERI_SAO_to_kMO (same numbers).!
!                                                                     !
!   Output index ranges (custom bounds):                              !
!     full  (CCSD)  : p,q,r,s = 1..nfuc                               !
!     (ov|ov) (MP2) : p,r = 1..o_hi (occ) ,  q,s = v_lo..nfuc (vir)    !
!                                                                     !
!   Memory: ERI_AO + ERI_out + X1 + nb * X2  (X1, X2 = 1/N_cell of    !
!   ERI_AO, complex) + f^4 work slabs:                                !
!     task mode  : one (k1,k2,k3) block per thread, 4 slabs PER THREAD !
!     shared mode: all threads on the same block, 4 slabs in total    !
!   task mode is used when its slabs are <= 5% of ERI_out, otherwise  !
!   shared mode (large nfuc, many threads). Same numbers in both.     !
!=====================================================================!

subroutine transform_ERI_AO_to_kMO(nBas, ERI_AO, coeff_SAO, p_lo, p_hi, q_lo, q_hi, ERI_out)

      use constants_module, only: dp, dpc, pi
      use unitcell_module , only: nfuc, nx, ny, nz, N_cell
      use omp_lib

      implicit none

      ! - input - !

      integer       , intent(in)  :: nBas
      real(dp)      , intent(in)  :: ERI_AO(nfuc,nBas,nBas,nBas)
      complex(dpc)  , intent(in)  :: coeff_SAO(nBas,nBas)
      integer       , intent(in)  :: p_lo, p_hi                 ! range of p and r
      integer       , intent(in)  :: q_lo, q_hi                 ! range of q and s

      ! - output - !

      complex(dpc)  , intent(out) :: ERI_out(p_lo:p_hi,q_lo:q_hi,p_lo:p_hi,q_lo:q_hi,N_cell,N_cell,N_cell)

      ! - local - !

      integer                     :: a, b, c, d, p, q, r, s
      integer                     :: j3, j4, lo, hi
      integer                     :: k1, k2, k3, k4, kb, ib, nb, nthreads, ntask, it
      integer                     :: r2, r3, r4
      integer                     :: ix, iy, iz, idx
      integer                     :: k3x, k3y, k3z
      integer                     :: off1, off2, off3, off4
      integer                     :: row, col
      logical                     :: task_mode
      double precision            :: slab_bytes, out_bytes
      integer     , allocatable   :: coord(:,:), batch(:)
      complex(dpc), allocatable   :: phase(:,:)
      complex(dpc), allocatable   :: X1(:,:,:,:)          ! (a, j2, j3, d)        one k4
      complex(dpc), allocatable   :: X2(:,:,:,:,:,:)      ! (a, j2, c, d, k3, ib) batch of k4
      complex(dpc), allocatable   :: slab(:,:,:,:), scr1(:,:,:,:), scr2(:,:,:,:), scr3(:,:,:,:)
      complex(dpc)                :: ph
      real(dp)                    :: theta, invN, max_offblock, coeff_scale
      real(dp)    , parameter     :: block_tol = 1.0e-10_dp

      !-----------------------------------------------------------------!

      if (min(nx,ny,nz) < 1)       error stop 'transform_ERI_AO_to_kMO: invalid cell grid'
      if (N_cell /= nx*ny*nz)      error stop 'transform_ERI_AO_to_kMO: inconsistent N_cell'
      if (nBas /= nfuc*N_cell)     error stop 'transform_ERI_AO_to_kMO: nBas /= nfuc*N_cell'

      ! ---- the coefficients must be k-block diagonal ---- !
      max_offblock = 0.0_dp
      coeff_scale  = maxval(abs(coeff_SAO))
      do col = 1, nBas
        do row = 1, nBas
          if ((row-1)/nfuc /= (col-1)/nfuc) max_offblock = max(max_offblock,abs(coeff_SAO(row,col)))
        end do
      end do
      if (max_offblock > block_tol*max(1.0_dp,coeff_scale)) then
        write(*,'(A,ES16.8)') 'Maximum off-k-block coefficient: ', max_offblock
        error stop 'transform_ERI_AO_to_kMO: coefficients not k-block ordered'
      end if

      ! ---- cell coordinates and phases e^{+i k.R} ---- !
      allocate(coord(3,N_cell), phase(N_cell,N_cell))
      do iz = 0, nz-1
        do iy = 0, ny-1
          do ix = 0, nx-1
            idx = iz*nx*ny + iy*nx + ix + 1
            coord(:,idx) = [ix,iy,iz]
          end do
        end do
      end do
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

      ! ---- batch of k4: enough (k4,k1,k2) tasks to keep all threads busy ---- !
      nthreads = omp_get_max_threads()
      nb = min(N_cell, max(1, (4*nthreads + N_cell*N_cell - 1)/(N_cell*N_cell)))
      allocate(batch(nb))
      allocate(X1(nfuc,nBas,nBas,nfuc))
      allocate(X2(nfuc,nBas,nfuc,nfuc,N_cell,nb))

      ! per-thread slabs only if they stay small compared with ERI_out
      slab_bytes = 16.d0*dble(nthreads)*( dble(nfuc)**4 + dble(nfuc)**3*dble(q_hi-q_lo+1)            &
                 + dble(nfuc)**2*dble(p_hi-p_lo+1)*dble(q_hi-q_lo+1)                             &
                 + dble(nfuc)*dble(p_hi-p_lo+1)*dble(q_hi-q_lo+1)**2 )
      out_bytes  = 16.d0*(dble(p_hi-p_lo+1)*dble(q_hi-q_lo+1))**2*dble(N_cell)**3
      task_mode  = (nthreads == 1) .or. (slab_bytes <= 0.05d0*out_bytes)
      if (.not. task_mode) then
        allocate(slab(nfuc,nfuc,nfuc,nfuc))
        allocate(scr1(nfuc,nfuc,nfuc,q_lo:q_hi))
        allocate(scr2(nfuc,nfuc,p_lo:p_hi,q_lo:q_hi))
        allocate(scr3(nfuc,q_lo:q_hi,p_lo:p_hi,q_lo:q_hi))
      end if

      do k4 = 1, N_cell, nb

        ! ------------- steps 1 and 2 for every k4 of the batch ------------- !
        do ib = 1, min(nb, N_cell-k4+1)
          kb = k4 + ib - 1
          batch(ib) = kb

          ! step 1 : R4 -> k4
          !$omp parallel do collapse(2) schedule(static) default(shared) private(d,j3,r4,j4,ph)
          do d = 1, nfuc
            do j3 = 1, nBas
              X1(:,:,j3,d) = (0.0_dp,0.0_dp)
              do r4 = 1, N_cell
                j4 = d + (r4-1)*nfuc
                ph = phase(kb,r4)
                X1(:,:,j3,d) = X1(:,:,j3,d) + ERI_AO(:,:,j3,j4)*ph
              end do
            end do
          end do
          !$omp end parallel do

          ! step 2 : R3 -> k3
          !$omp parallel do collapse(3) schedule(static) default(shared) private(k3,d,c,r3,j3,ph)
          do k3 = 1, N_cell
            do d = 1, nfuc
              do c = 1, nfuc
                X2(:,:,c,d,k3,ib) = (0.0_dp,0.0_dp)
                do r3 = 1, N_cell
                  j3 = c + (r3-1)*nfuc
                  ph = conjg(phase(k3,r3))
                  X2(:,:,c,d,k3,ib) = X2(:,:,c,d,k3,ib) + X1(:,:,j3,d)*ph
                end do
              end do
            end do
          end do
          !$omp end parallel do
        end do

        ! ------------- step 3 + quarter transforms : (ib, k1, k2) ------------- !
        ntask = min(nb, N_cell-k4+1) * N_cell * N_cell

        if (.not. task_mode) then
          do it = 0, ntask-1
            call shared_block(it)
          end do
          cycle
        end if

        !$omp parallel default(shared) &
        !$omp   private(it,ib,kb,k1,k2,k3,k3x,k3y,k3z,off1,off2,off3,off4, &
        !$omp           a,b,c,d,p,q,r,s,r2,lo,hi,ph,slab,scr1,scr2,scr3)
        allocate(slab(nfuc,nfuc,nfuc,nfuc))
        allocate(scr1(nfuc,nfuc,nfuc,q_lo:q_hi))
        allocate(scr2(nfuc,nfuc,p_lo:p_hi,q_lo:q_hi))
        allocate(scr3(nfuc,q_lo:q_hi,p_lo:p_hi,q_lo:q_hi))

        !$omp do schedule(dynamic)
        do it = 0, ntask-1
          ib = it / (N_cell*N_cell) + 1
          k1 = mod(it / N_cell, N_cell) + 1
          k2 = mod(it, N_cell) + 1
          kb = batch(ib)

          ! k3 = k4 - k1 + k2  (k4 = k1 - k2 + k3)
          k3x = modulo(coord(1,kb)-coord(1,k1)+coord(1,k2),nx)
          k3y = modulo(coord(2,kb)-coord(2,k1)+coord(2,k2),ny)
          k3z = modulo(coord(3,kb)-coord(3,k1)+coord(3,k2),nz)
          k3  = k3z*nx*ny + k3y*nx + k3x + 1

          off1 = (k1-1)*nfuc
          off2 = (k2-1)*nfuc
          off3 = (k3-1)*nfuc
          off4 = (kb-1)*nfuc

          ! step 3 : R2 -> k2 , 1/N
          do d = 1, nfuc
            do c = 1, nfuc
              slab(:,:,c,d) = (0.0_dp,0.0_dp)
              do r2 = 1, N_cell
                lo = (r2-1)*nfuc + 1
                hi = r2*nfuc
                ph = phase(k2,r2)*invN
                slab(:,:,c,d) = slab(:,:,c,d) + X2(:,lo:hi,c,d,k3,ib)*ph
              end do
            end do
          end do

          ! 4th index : d -> s
          scr1 = (0.0_dp,0.0_dp)
          do s = q_lo, q_hi
            do d = 1, nfuc
              do c = 1, nfuc
                do b = 1, nfuc
                  do a = 1, nfuc
                    scr1(a,b,c,s) = scr1(a,b,c,s) + slab(a,b,c,d)*coeff_SAO(off4+d,off4+s)
                  end do
                end do
              end do
            end do
          end do
          ! 3rd index : c -> r (conjugated)
          scr2 = (0.0_dp,0.0_dp)
          do s = q_lo, q_hi
            do r = p_lo, p_hi
              do c = 1, nfuc
                do b = 1, nfuc
                  do a = 1, nfuc
                    scr2(a,b,r,s) = scr2(a,b,r,s) + scr1(a,b,c,s)*conjg(coeff_SAO(off3+c,off3+r))
                  end do
                end do
              end do
            end do
          end do
          ! 2nd index : b -> q
          scr3 = (0.0_dp,0.0_dp)
          do s = q_lo, q_hi
            do r = p_lo, p_hi
              do q = q_lo, q_hi
                do b = 1, nfuc
                  do a = 1, nfuc
                    scr3(a,q,r,s) = scr3(a,q,r,s) + scr2(a,b,r,s)*coeff_SAO(off2+b,off2+q)
                  end do
                end do
              end do
            end do
          end do
          ! 1st index : a -> p (conjugated) ; each task owns this slab
          ERI_out(:,:,:,:,k1,k2,k3) = (0.0_dp,0.0_dp)
          do s = q_lo, q_hi
            do r = p_lo, p_hi
              do q = q_lo, q_hi
                do p = p_lo, p_hi
                  do a = 1, nfuc
                    ERI_out(p,q,r,s,k1,k2,k3) = ERI_out(p,q,r,s,k1,k2,k3) &
                      + conjg(coeff_SAO(off1+a,off1+p))*scr3(a,q,r,s)
                  end do
                end do
              end do
            end do
          end do
        end do
        !$omp end do

        deallocate(slab, scr1, scr2, scr3)
        !$omp end parallel

      end do

      if (.not. task_mode) deallocate(slab, scr1, scr2, scr3)
      deallocate(X1, X2, batch, coord, phase)

contains

      !-----------------------------------------------------------------!
      ! shared mode: one (k1,k2,k3) block, all threads split its loops  !
      ! (same operations and summation order as the task mode)          !
      !-----------------------------------------------------------------!
      subroutine shared_block(itask)

        integer, intent(in) :: itask
        integer :: ib_, kb_, k1_, k2_, k3_, o1, o2, o3, o4, a_, b_, c_, d_, p_, q_, r_, s_, r2_, lo_, hi_
        complex(dpc) :: ph_

        ib_ = itask / (N_cell*N_cell) + 1
        k1_ = mod(itask / N_cell, N_cell) + 1
        k2_ = mod(itask, N_cell) + 1
        kb_ = batch(ib_)
        k3_ = modulo(coord(3,kb_)-coord(3,k1_)+coord(3,k2_),nz)*nx*ny &
            + modulo(coord(2,kb_)-coord(2,k1_)+coord(2,k2_),ny)*nx    &
            + modulo(coord(1,kb_)-coord(1,k1_)+coord(1,k2_),nx) + 1
        o1 = (k1_-1)*nfuc ; o2 = (k2_-1)*nfuc ; o3 = (k3_-1)*nfuc ; o4 = (kb_-1)*nfuc

        !$omp parallel default(shared) private(a_,b_,c_,d_,p_,q_,r_,s_,r2_,lo_,hi_,ph_)

        ! step 3 : R2 -> k2 , 1/N
        !$omp do collapse(2) schedule(static)
        do d_ = 1, nfuc
          do c_ = 1, nfuc
            slab(:,:,c_,d_) = (0.0_dp,0.0_dp)
            do r2_ = 1, N_cell
              lo_ = (r2_-1)*nfuc + 1
              hi_ = r2_*nfuc
              ph_ = phase(k2_,r2_)*invN
              slab(:,:,c_,d_) = slab(:,:,c_,d_) + X2(:,lo_:hi_,c_,d_,k3_,ib_)*ph_
            end do
          end do
        end do
        !$omp end do

        ! 4th index : d -> s
        !$omp do collapse(2) schedule(static)
        do s_ = q_lo, q_hi
          do c_ = 1, nfuc
            scr1(:,:,c_,s_) = (0.0_dp,0.0_dp)
            do d_ = 1, nfuc
              do b_ = 1, nfuc
                do a_ = 1, nfuc
                  scr1(a_,b_,c_,s_) = scr1(a_,b_,c_,s_) + slab(a_,b_,c_,d_)*coeff_SAO(o4+d_,o4+s_)
                end do
              end do
            end do
          end do
        end do
        !$omp end do

        ! 3rd index : c -> r (conjugated)
        !$omp do collapse(2) schedule(static)
        do s_ = q_lo, q_hi
          do r_ = p_lo, p_hi
            scr2(:,:,r_,s_) = (0.0_dp,0.0_dp)
            do c_ = 1, nfuc
              do b_ = 1, nfuc
                do a_ = 1, nfuc
                  scr2(a_,b_,r_,s_) = scr2(a_,b_,r_,s_) + scr1(a_,b_,c_,s_)*conjg(coeff_SAO(o3+c_,o3+r_))
                end do
              end do
            end do
          end do
        end do
        !$omp end do

        ! 2nd index : b -> q
        !$omp do collapse(2) schedule(static)
        do s_ = q_lo, q_hi
          do r_ = p_lo, p_hi
            scr3(:,:,r_,s_) = (0.0_dp,0.0_dp)
            do q_ = q_lo, q_hi
              do b_ = 1, nfuc
                do a_ = 1, nfuc
                  scr3(a_,q_,r_,s_) = scr3(a_,q_,r_,s_) + scr2(a_,b_,r_,s_)*coeff_SAO(o2+b_,o2+q_)
                end do
              end do
            end do
          end do
        end do
        !$omp end do

        ! 1st index : a -> p (conjugated)
        !$omp do collapse(2) schedule(static)
        do s_ = q_lo, q_hi
          do r_ = p_lo, p_hi
            ERI_out(:,:,r_,s_,k1_,k2_,k3_) = (0.0_dp,0.0_dp)
            do q_ = q_lo, q_hi
              do p_ = p_lo, p_hi
                do a_ = 1, nfuc
                  ERI_out(p_,q_,r_,s_,k1_,k2_,k3_) = ERI_out(p_,q_,r_,s_,k1_,k2_,k3_) &
                    + conjg(coeff_SAO(o1+a_,o1+p_))*scr3(a_,q_,r_,s_)
                end do
              end do
            end do
          end do
        end do
        !$omp end do

        !$omp end parallel

      end subroutine shared_block

end subroutine transform_ERI_AO_to_kMO
