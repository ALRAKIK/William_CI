!=====================================================================!
!                                                                     !
!   [Nov 2026] keyword direct: k-space MO integrals for MP2 WITHOUT   !
!   the stored AO integrals (integral-direct).                        !
!                                                                     !
!   Same result as transform_ERI_AO_to_kMO (AO_to_MO_SAO.f90):        !
!     ERI_out(p,q,r,s,k1,k2,k3) = (p k1 q k2 | r k3 s k4),            !
!     k4 = k1 - k2 + k3 , p,r in p_lo..p_hi , q,s in q_lo..q_hi       !
!                                                                     !
!   The transform is linear in the 4th AO index d, so it is done one  !
!   unit-cell shell SD at a time (d in SD, every cell R4):            !
!     slab(a,j2,j3,dd,R4) = (a0 j2|j3, d R4)   (eri_direct:direct_slab)!
!     for every k4:                                                   !
!       X1(a,j2,j3,dd)    = sum_R4 slab e^{+i k4 R4}                  !
!       X2(a,j2,c,dd,k3)  = sum_R3 X1(a,j2,c R3,dd) e^{-i k3 R3}      !
!       for every (k1,k2), k3 = k4 - k1 + k2 :                        !
!         S(a,b,c,dd)     = 1/N sum_R2 X2(a,b R2,c,dd,k3) e^{+i k2 R2} !
!         ERI_out += C1^* C2 C3^* C4(dd,s) S  (4 quarter transforms,  !
!                                              zgemm)                 !
!   Memory: ERI_out + slab (= stored SAO integrals * n_SD / nfuc,     !
!   real) + X1 + X2 (complex, n_SD/N_cell of that).                   !
!   Integral work: every (SA in cell 0, SJ | SK, SD R4) shell quartet !
!   once per run, without the 8-fold symmetry.                        !
!=====================================================================!

subroutine transform_ERI_direct_kMO(nBas, coeff_SAO, p_lo, p_hi, q_lo, q_hi, ERI_out, iout)

      use constants_module, only: dp, dpc, pi
      use unitcell_module , only: nfuc, nx, ny, nz, N_cell
      use eri_direct      , only: direct_slab, d_nshuc, d_sh_first, d_sh_n
      use omp_lib

      implicit none

      integer       , intent(in)  :: nBas, iout
      complex(dpc)  , intent(in)  :: coeff_SAO(nBas,nBas)
      integer       , intent(in)  :: p_lo, p_hi, q_lo, q_hi
      complex(dpc)  , intent(out) :: ERI_out(p_lo:p_hi,q_lo:q_hi,p_lo:p_hi,q_lo:q_hi,N_cell,N_cell,N_cell)

      integer :: SD, nd, d0, dd, k1, k2, k3, k4, r2, r3, r4, j3, c, b, ix, iy, iz, idx, it, ntask
      integer :: off1, off2, off3, off4, no, nv, lo, hi, r, s
      integer :: coord(3,N_cell)
      real(dp) :: theta, invN, t_int, t_tr, w0
      complex(dpc) :: phase(N_cell,N_cell), ph
      real(dp)    , allocatable :: slab(:,:,:,:,:)
      complex(dpc), allocatable :: X1(:,:,:,:), X2(:,:,:,:,:)
      complex(dpc), allocatable :: Sb(:,:,:,:), scr1(:,:,:,:), scr2(:,:,:,:), scr3(:,:,:,:), Cd(:,:), C3(:,:), tmp(:,:)
      complex(dpc), parameter :: one=(1.0_dp,0.0_dp), zero=(0.0_dp,0.0_dp)

      if (N_cell /= nx*ny*nz)  error stop 'transform_ERI_direct_kMO: inconsistent N_cell'
      if (nBas /= nfuc*N_cell) error stop 'transform_ERI_direct_kMO: nBas /= nfuc*N_cell'
      no = p_hi - p_lo + 1 ; nv = q_hi - q_lo + 1

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

      ERI_out = zero
      t_int = 0.0_dp ; t_tr = 0.0_dp
      ntask = N_cell*N_cell

      do SD = 1, d_nshuc
        nd = d_sh_n(SD) ; d0 = d_sh_first(SD)

        ! ---- the integrals of this shell (all cells R4) ---- !
        w0 = omp_get_wtime()
        allocate(slab(nfuc,nBas,nBas,nd,N_cell))
        call direct_slab(SD, nfuc, nBas, N_cell, slab)
        t_int = t_int + omp_get_wtime() - w0

        w0 = omp_get_wtime()
        allocate(X1(nfuc,nBas,nBas,nd), X2(nfuc,nBas,nfuc,nd,N_cell))
        do k4 = 1, N_cell

          ! step 1 : R4 -> k4
          !$omp parallel do collapse(2) schedule(static) private(dd,j3,r4)
          do dd = 1, nd
            do j3 = 1, nBas
              X1(:,:,j3,dd) = zero
              do r4 = 1, N_cell
                X1(:,:,j3,dd) = X1(:,:,j3,dd) + slab(:,:,j3,dd,r4)*phase(k4,r4)
              end do
            end do
          end do
          !$omp end parallel do

          ! step 2 : R3 -> k3
          !$omp parallel do collapse(3) schedule(static) private(k3,dd,c,r3,j3)
          do k3 = 1, N_cell
            do dd = 1, nd
              do c = 1, nfuc
                X2(:,:,c,dd,k3) = zero
                do r3 = 1, N_cell
                  j3 = c + (r3-1)*nfuc
                  X2(:,:,c,dd,k3) = X2(:,:,c,dd,k3) + X1(:,:,j3,dd)*conjg(phase(k3,r3))
                end do
              end do
            end do
          end do
          !$omp end parallel do

          ! step 3 + quarter transforms, one (k1,k2) block per task
          !$omp parallel default(shared) private(it,k1,k2,k3,off1,off2,off3,off4,Sb,scr1,scr2,scr3,Cd,C3,tmp,dd,c,b,r2,lo,hi,r,s)
          allocate(Sb(nfuc,nfuc,nfuc,nd), scr1(nfuc,nfuc,nfuc,q_lo:q_hi), scr2(nfuc,nfuc,p_lo:p_hi,q_lo:q_hi), &
                   scr3(nfuc,q_lo:q_hi,p_lo:p_hi,q_lo:q_hi), Cd(nd,q_lo:q_hi), C3(nfuc,p_lo:p_hi), tmp(nfuc,q_lo:q_hi))
          !$omp do schedule(dynamic)
          do it = 0, ntask-1
            k1 = it/N_cell + 1
            k2 = mod(it,N_cell) + 1
            k3 = modulo(coord(3,k4)-coord(3,k1)+coord(3,k2),nz)*nx*ny &
               + modulo(coord(2,k4)-coord(2,k1)+coord(2,k2),ny)*nx    &
               + modulo(coord(1,k4)-coord(1,k1)+coord(1,k2),nx) + 1
            off1 = (k1-1)*nfuc ; off2 = (k2-1)*nfuc ; off3 = (k3-1)*nfuc ; off4 = (k4-1)*nfuc

            ! step 3 : R2 -> k2 , 1/N
            do dd = 1, nd
              do c = 1, nfuc
                Sb(:,:,c,dd) = zero
                do r2 = 1, N_cell
                  lo = (r2-1)*nfuc + 1 ; hi = r2*nfuc
                  Sb(:,:,c,dd) = Sb(:,:,c,dd) + X2(:,lo:hi,c,dd,k3)*(phase(k2,r2)*invN)
                end do
              end do
            end do
            ! 4th index : d -> s   scr1(a,b,c,s) = Sb(a,b,c,dd) C4(d,s)
            Cd = coeff_SAO(off4+d0:off4+d0+nd-1, off4+q_lo:off4+q_hi)
            call zgemm('N','N', nfuc**3, nv, nd, one, Sb, nfuc**3, Cd, nd, zero, scr1, nfuc**3)
            ! 3rd index : c -> r (conjugated)   scr2(a,b,r,s) = scr1(a,b,c,s) C3^*(c,r)
            C3 = conjg(coeff_SAO(off3+1:off3+nfuc, off3+p_lo:off3+p_hi))
            do s = q_lo, q_hi
              call zgemm('N','N', nfuc*nfuc, no, nfuc, one, scr1(1,1,1,s), nfuc*nfuc, C3, nfuc, zero, &
                         scr2(1,1,p_lo,s), nfuc*nfuc)
            end do
            ! 2nd index : b -> q   scr3(a,q,r,s) = scr2(a,b,r,s) C2(b,q)
            do s = q_lo, q_hi
              do r = p_lo, p_hi
                call zgemm('N','N', nfuc, nv, nfuc, one, scr2(1,1,r,s), nfuc, coeff_SAO(off2+1,off2+q_lo), nBas, &
                           zero, tmp, nfuc)
                scr3(:,:,r,s) = tmp
              end do
            end do
            ! 1st index : a -> p (conjugated), accumulated over the shells SD
            call zgemm('C','N', no, nv*no*nv, nfuc, one, coeff_SAO(off1+1,off1+p_lo), nBas, scr3, nfuc, one, &
                       ERI_out(p_lo,q_lo,p_lo,q_lo,k1,k2,k3), no)
          end do
          !$omp end do
          deallocate(Sb, scr1, scr2, scr3, Cd, C3, tmp)
          !$omp end parallel

        end do
        deallocate(slab, X1, X2)
        t_tr = t_tr + omp_get_wtime() - w0
      end do

      write(iout,'(a,f10.2,a,f10.2,a)') '  [direct] MP2 integrals: kernel ', t_int, ' s , transformation ', t_tr, ' s'
      write(iout,'(a,f9.3,a,f9.3,a)') '  [direct] MP2 memory: (ov|ov) ', &
           16.d0*dble(no*nv)**2*dble(N_cell)**3/1024.d0**3, ' GB , largest shell slab ', &
           8.d0*dble(nfuc)*dble(nBas)**2*dble(maxval(d_sh_n(1:d_nshuc)))*dble(N_cell)/1024.d0**3, ' GB'
      flush(iout)

end subroutine transform_ERI_direct_kMO
