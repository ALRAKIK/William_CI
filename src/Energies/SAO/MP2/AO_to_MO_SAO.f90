subroutine transform_ERI_AO_to_SAO(nBas, ERI_AO, ERI_SAO)

      use constants_module, only: dp, dpc, pi
      use unitcell_module, only: nfuc, nx, ny, nz, N_cell
      implicit none

      ! Input: first AO index belongs to the reference cell.
      integer, intent(in) :: nBas
      real(dp), intent(in) :: ERI_AO(nfuc,nBas,nBas,nBas)

      ! Output:
      ! ERI_SAO(a,b,c,d,k1,k2,k3) = (a k1,b k2 | c k3,d k4)
      ! k4 is determined by momentum conservation.
      complex(dpc), intent(out) :: ERI_SAO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      integer :: a, b, c, d
      integer :: k1, k2, k3, k4
      integer :: r2, r3, r4
      integer :: j2, j3, j4
      integer :: ix, iy, iz, idx
      integer :: k4x, k4y, k4z

      integer, allocatable      :: coord(:,:)
      complex(dpc), allocatable :: phase(:,:)

      real(dp) :: theta
      complex(dpc) :: val, phase23

      ! Check consistency of the cell grid and AO ordering.
      if (min(nx,ny,nz) < 1) then
        error stop 'transform_ERI_AO_to_SAO: invalid cell grid'
      end if

      if (N_cell /= nx*ny*nz) then
        error stop 'transform_ERI_AO_to_SAO: inconsistent N_cell'
      end if

      if (nBas /= nfuc*N_cell) then
        error stop 'transform_ERI_AO_to_SAO: nBas /= nfuc*N_cell'
      end if

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

          phase(k1,r2) = exp(cmplx(0.0_dp,theta,kind=dpc))
        end do
      end do

      do k1 = 1, N_cell
        do k2 = 1, N_cell
          do k3 = 1, N_cell

            ! Component-wise momentum conservation.
            k4x = modulo(coord(1,k1)-coord(1,k2)+coord(1,k3),nx)
            k4y = modulo(coord(2,k1)-coord(2,k2)+coord(2,k3),ny)
            k4z = modulo(coord(3,k1)-coord(3,k2)+coord(3,k3),nz)

            k4 = k4z*nx*ny + k4y*nx + k4x + 1

            do d = 1, nfuc
              do c = 1, nfuc
                do b = 1, nfuc
                  do a = 1, nfuc

                    val = cmplx(0.0_dp,0.0_dp,kind=dpc)

                    do r4 = 1, N_cell
                      j4 = d + (r4-1)*nfuc

                      do r3 = 1, N_cell
                        j3 = c + (r3-1)*nfuc

                        do r2 = 1, N_cell
                          j2 = b + (r2-1)*nfuc

                          phase23 = phase(k2,r2)*conjg(phase(k3,r3))

                          val = val + ERI_AO(a,j2,j3,j4) &
                                    * phase23 * phase(k4,r4)
                        end do
                      end do
                    end do

                    ERI_SAO(a,b,c,d,k1,k2,k3) = &
                      val / real(N_cell,dp)

                  end do
                end do
              end do
            end do

          end do
        end do
      end do

      deallocate(coord,phase)

end subroutine transform_ERI_AO_to_SAO

subroutine transform_ERI_SAO_to_kMO(nBas, coeff_SAO, ERI_SAO, ERI_kMO)

      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, nx, ny, nz, N_cell
      implicit none

      integer, intent(in) :: nBas

      ! Full MO coefficients expressed in the SAO basis.
      ! Rows and columns must both be ordered by k block.
      complex(dpc), intent(in) :: coeff_SAO(nBas,nBas)

      complex(dpc), intent(in) :: &
        ERI_SAO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      ! All orbital combinations, not just occupied-virtual.
      complex(dpc), intent(out) :: &
        ERI_kMO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      integer :: k1, k2, k3, k4
      integer :: ix, iy, iz, idx
      integer :: k4x, k4y, k4z
      integer :: off1, off2, off3, off4
      integer :: a, b, c, d, p, q, r, s
      integer :: row, col

      integer, allocatable :: coord(:,:)
      complex(dpc), allocatable :: scr1(:,:,:,:)
      complex(dpc), allocatable :: scr2(:,:,:,:)
      complex(dpc), allocatable :: scr3(:,:,:,:)

      real(dp) :: max_offblock, coeff_scale
      real(dp), parameter :: block_tol = 1.0e-10_dp

      if (min(nx,ny,nz) < 1) then
        error stop 'transform_ERI_SAO_to_kMO: invalid grid'
      end if

      if (N_cell /= nx*ny*nz) then
        error stop 'transform_ERI_SAO_to_kMO: inconsistent N_cell'
      end if

      if (nBas /= nfuc*N_cell) then
        error stop 'transform_ERI_SAO_to_kMO: inconsistent nBas'
      end if

      ! Reject significant mixing between different k blocks.
      max_offblock = 0.0_dp
      coeff_scale = maxval(abs(coeff_SAO))

      do col = 1, nBas
        do row = 1, nBas
          if ((row-1)/nfuc /= (col-1)/nfuc) then
            max_offblock = max(max_offblock,abs(coeff_SAO(row,col)))
          end if
        end do
      end do

      if (max_offblock > block_tol*max(1.0_dp,coeff_scale)) then
        write(*,'(A,ES16.8)') 'Maximum off-k-block coefficient: ', &
          max_offblock
        error stop 'transform_ERI_SAO_to_kMO: coefficients not k-block ordered'
      end if

      allocate(coord(3,N_cell))
      allocate(scr1(nfuc,nfuc,nfuc,nfuc))
      allocate(scr2(nfuc,nfuc,nfuc,nfuc))
      allocate(scr3(nfuc,nfuc,nfuc,nfuc))

      ! Zero-based coordinates, x fastest.
      do iz = 0, nz-1
        do iy = 0, ny-1
          do ix = 0, nx-1
            idx = iz*nx*ny + iy*nx + ix + 1
            coord(:,idx) = [ix,iy,iz]
          end do
        end do
      end do

      do k1 = 1, N_cell
        off1 = (k1-1)*nfuc

        do k2 = 1, N_cell
          off2 = (k2-1)*nfuc

          do k3 = 1, N_cell
            off3 = (k3-1)*nfuc

            k4x = modulo(coord(1,k1)-coord(1,k2)+coord(1,k3),nx)
            k4y = modulo(coord(2,k1)-coord(2,k2)+coord(2,k3),ny)
            k4z = modulo(coord(3,k1)-coord(3,k2)+coord(3,k3),nz)

            k4 = k4z*nx*ny + k4y*nx + k4x + 1
            off4 = (k4-1)*nfuc

            ! Fourth index: d -> s, no conjugation.
            scr1 = cmplx(0.0_dp,0.0_dp,kind=dpc)

            do s = 1, nfuc
              do d = 1, nfuc
                do c = 1, nfuc
                  do b = 1, nfuc
                    do a = 1, nfuc
                      scr1(a,b,c,s) = scr1(a,b,c,s) &
                        + ERI_SAO(a,b,c,d,k1,k2,k3) &
                        * coeff_SAO(off4+d,off4+s)
                    end do
                  end do
                end do
              end do
            end do

            ! Third index: c -> r, conjugated.
            scr2 = cmplx(0.0_dp,0.0_dp,kind=dpc)

            do s = 1, nfuc
              do r = 1, nfuc
                do c = 1, nfuc
                  do b = 1, nfuc
                    do a = 1, nfuc
                      scr2(a,b,r,s) = scr2(a,b,r,s) &
                        + scr1(a,b,c,s) &
                        * conjg(coeff_SAO(off3+c,off3+r))
                    end do
                  end do
                end do
              end do
            end do

            ! Second index: b -> q, no conjugation.
            scr3 = cmplx(0.0_dp,0.0_dp,kind=dpc)

            do s = 1, nfuc
              do r = 1, nfuc
                do q = 1, nfuc
                  do b = 1, nfuc
                    do a = 1, nfuc
                      scr3(a,q,r,s) = scr3(a,q,r,s) &
                        + scr2(a,b,r,s) &
                        * coeff_SAO(off2+b,off2+q)
                    end do
                  end do
                end do
              end do
            end do

            ! First index: a -> p, conjugated.
            ERI_kMO(:,:,:,:,k1,k2,k3) = &
              cmplx(0.0_dp,0.0_dp,kind=dpc)

            do s = 1, nfuc
              do r = 1, nfuc
                do q = 1, nfuc
                  do p = 1, nfuc
                    do a = 1, nfuc
                      ERI_kMO(p,q,r,s,k1,k2,k3) = &
                        ERI_kMO(p,q,r,s,k1,k2,k3) &
                        + conjg(coeff_SAO(off1+a,off1+p)) &
                        * scr3(a,q,r,s)
                    end do
                  end do
                end do
              end do
            end do

          end do
        end do
      end do

      deallocate(coord,scr1,scr2,scr3)

end subroutine transform_ERI_SAO_to_kMO

subroutine check_ERI_SAO_roundtrip(nBas, ERI_AO, ERI_SAO)

      use constants_module, only: dp, dpc, pi
      use unitcell_module, only: nfuc, nx, ny, nz, N_cell
      implicit none

      integer, intent(in) :: nBas
      real(dp), intent(in) :: ERI_AO(nfuc,nBas,nBas,nBas)
      complex(dpc), intent(in) :: &
        ERI_SAO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      integer :: a, b, c, d
      integer :: r2, r3, r4, j2, j3, j4
      integer :: k1, k2, k3, k4
      integer :: ix, iy, iz, idx
      integer :: k4x, k4y, k4z
      integer :: worst(4)

      integer, allocatable :: coord(:,:)
      complex(dpc), allocatable :: phase(:,:)

      real(dp) :: theta, err
      real(dp) :: max_error, max_imag, ao_scale
      real(dp) :: original_worst
      complex(dpc) :: reconstructed, reconstructed_worst

      if (min(nx,ny,nz) < 1) then
        error stop 'check_ERI_SAO_roundtrip: invalid grid'
      end if
      if (N_cell /= nx*ny*nz .or. nBas /= nfuc*N_cell) then
        error stop 'check_ERI_SAO_roundtrip: inconsistent dimensions'
      end if

      allocate(coord(3,N_cell), phase(N_cell,N_cell))

      do iz = 0, nz-1
        do iy = 0, ny-1
          do ix = 0, nx-1
            idx = iz*nx*ny + iy*nx + ix + 1
            coord(:,idx) = [ix,iy,iz]
          end do
        end do
      end do

      ! Positive phases; the inverse uses their conjugates.
      do k1 = 1, N_cell
        do r2 = 1, N_cell
          theta = 2.0_dp*pi * ( &
            real(coord(1,k1),dp)*real(coord(1,r2),dp)/real(nx,dp) &
          + real(coord(2,k1),dp)*real(coord(2,r2),dp)/real(ny,dp) &
          + real(coord(3,k1),dp)*real(coord(3,r2),dp)/real(nz,dp))

          phase(k1,r2) = exp(cmplx(0.0_dp,theta,kind=dpc))
        end do
      end do

      max_error = 0.0_dp
      max_imag  = 0.0_dp
      ao_scale  = maxval(abs(ERI_AO))
      worst = [1,1,1,1]
      original_worst = ERI_AO(1,1,1,1)
      reconstructed_worst = cmplx(0.0_dp,0.0_dp,kind=dpc)

      do r4 = 1, N_cell
        do r3 = 1, N_cell
          do r2 = 1, N_cell
            do d = 1, nfuc
              j4 = d + (r4-1)*nfuc
              do c = 1, nfuc
                j3 = c + (r3-1)*nfuc
                do b = 1, nfuc
                  j2 = b + (r2-1)*nfuc
                  do a = 1, nfuc

                    reconstructed = cmplx(0.0_dp,0.0_dp,kind=dpc)

                    do k3 = 1, N_cell
                      do k2 = 1, N_cell
                        do k1 = 1, N_cell
                          k4x = modulo(coord(1,k1)-coord(1,k2) &
                                     + coord(1,k3),nx)
                          k4y = modulo(coord(2,k1)-coord(2,k2) &
                                     + coord(2,k3),ny)
                          k4z = modulo(coord(3,k1)-coord(3,k2) &
                                     + coord(3,k3),nz)

                          k4 = k4z*nx*ny + k4y*nx + k4x + 1

                          reconstructed = reconstructed &
                            + ERI_SAO(a,b,c,d,k1,k2,k3) &
                            * conjg(phase(k2,r2)) &
                            * phase(k3,r3) &
                            * conjg(phase(k4,r4))
                        end do
                      end do
                    end do

                    reconstructed = reconstructed / &
                      (real(N_cell,dp)*real(N_cell,dp))

                    err = abs(reconstructed - ERI_AO(a,j2,j3,j4))
                    max_imag = max(max_imag,abs(aimag(reconstructed)))

                    if (err > max_error) then
                      max_error = err
                      worst = [a,j2,j3,j4]
                      original_worst = ERI_AO(a,j2,j3,j4)
                      reconstructed_worst = reconstructed
                    end if

                  end do
                end do
              end do
            end do
          end do
        end do
      end do

      write(*,'(/,A)') 'ERI AO -> SAO -> AO round-trip check'
      write(*,'(A,ES16.8)') 'Maximum AO magnitude       : ', ao_scale
      write(*,'(A,ES16.8)') 'Maximum reconstruction error: ', max_error
      write(*,'(A,ES16.8)') 'Maximum reconstructed Im    : ', max_imag

      if (ao_scale > 0.0_dp) then
        write(*,'(A,ES16.8)') 'Error / maximum AO magnitude: ', &
          max_error/ao_scale
      end if

      write(*,'(A,4I6)')    'Worst AO indices            : ', worst
      write(*,'(A,ES24.16)') 'Original AO integral        : ', original_worst
      write(*,'(A,2ES24.16)') 'Reconstructed AO (Re, Im)   : ', &
        real(reconstructed_worst,dp), aimag(reconstructed_worst)

      deallocate(coord,phase)

end subroutine check_ERI_SAO_roundtrip


subroutine check_ERI_SAO_symmetries(ERI_SAO)

      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, nx, ny, nz, N_cell
      implicit none

      complex(dpc), intent(in) :: &
        ERI_SAO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      integer :: a, b, c, d
      integer :: k1, k2, k3, k4
      integer :: ix, iy, iz, idx
      integer :: k4x, k4y, k4z
      integer, allocatable :: coord(:,:)

      real(dp) :: err_perm, max_err_perm
      real(dp) :: err_conj, max_err_conj
      real(dp) :: max_imag_diag, min_real_diag
      complex(dpc) :: val1, val_perm, val_conj

      allocate(coord(3,N_cell))

      do iz = 0, nz-1
        do iy = 0, ny-1
          do ix = 0, nx-1
            idx = iz*nx*ny + iy*nx + ix + 1
            coord(:,idx) = [ix,iy,iz]
          end do
        end do
      end do

      max_err_perm  = 0.0_dp
      max_err_conj  = 0.0_dp
      max_imag_diag = 0.0_dp
      min_real_diag = 1.0d30

      do k1 = 1, N_cell
        do k2 = 1, N_cell
          do k3 = 1, N_cell

            ! Component-wise momentum conservation: k4 = k1 - k2 + k3
            k4x = modulo(coord(1,k1)-coord(1,k2)+coord(1,k3), nx)
            k4y = modulo(coord(2,k1)-coord(2,k2)+coord(2,k3), ny)
            k4z = modulo(coord(3,k1)-coord(3,k2)+coord(3,k3), nz)
            k4  = k4z*nx*ny + k4y*nx + k4x + 1

            do d = 1, nfuc
              do c = 1, nfuc
                do b = 1, nfuc
                  do a = 1, nfuc

                    val1 = ERI_SAO(a,b,c,d, k1,k2,k3)

                    ! 1. Permutation of electron 1 and electron 2:
                    !    (a k1, b k2 | c k3, d k4) = (c k3, d k4 | a k1, b k2)
                    !    For (c k3, d k4 | a k1, b k2), the first three momenta are (k3, k4, k1)
                    val_perm = ERI_SAO(c,d,a,b, k3,k4,k1)
                    err_perm = abs(val1 - val_perm)
                    max_err_perm = max(max_err_perm, err_perm)

                    ! 2. Complex conjugation:
                    !    (a k1, b k2 | c k3, d k4)* = (b k2, a k1 | d k4, c k3)
                    !    For (b k2, a k1 | d k4, c k3), the first three momenta are (k2, k1, k4)
                    val_conj = ERI_SAO(b,a,d,c, k2,k1,k4)
                    err_conj = abs(conjg(val1) - val_conj)
                    max_err_conj = max(max_err_conj, err_conj)

                  end do
                end do
              end do
            end do

          end do
        end do
      end do

      ! 3. Check diagonal self-repulsion: (a k1, a k1 | a k1, a k1)
      do k1 = 1, N_cell
        do a = 1, nfuc
          val1 = ERI_SAO(a,a,a,a, k1,k1,k1)
          max_imag_diag = max(max_imag_diag, abs(aimag(val1)))
          min_real_diag = min(min_real_diag, real(val1, dp))
        end do
      end do

      write(*,'(/,A)') '=================================================='
      write(*,'(A)')   '         ERI_SAO Physical Symmetry Checks         '
      write(*,'(A)')   '=================================================='
      write(*,'(A,ES16.8)') 'Max Electron Permutation Error : ', max_err_perm
      write(*,'(A,ES16.8)') 'Max Complex Conjugation Error  : ', max_err_conj
      write(*,'(A,ES16.8)') 'Max Imag part on Diagonal      : ', max_imag_diag
      write(*,'(A,ES16.8)') 'Min Real part on Diagonal      : ', min_real_diag
      write(*,'(A,/)') '=================================================='

      deallocate(coord)

end subroutine check_ERI_SAO_symmetries





subroutine check_ERI_kMO_symmetry(N_cell, blk, ERI_kMO, max_occ, max_vir)

      use constants_module, only: dp, dpc
      use unitcell_module, only: nx, ny, nz
      use k_blocks
      implicit none

      integer, intent(in) :: N_cell, max_occ, max_vir
      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc), intent(in) :: ERI_kMO(max_occ, max_vir, max_occ, max_vir, N_cell, N_cell, N_cell)

      integer :: ki, ka, kj, kb
      integer :: i, a, j, b
      integer :: nocc_i, nocc_j, nvir_a, nvir_b
      integer :: kbx, kby, kbz, idx, ix, iy, iz
      integer, allocatable :: coord(:,:)
      real(dp) :: err, max_err
      complex(dpc) :: g1, g2

      allocate(coord(3, N_cell))
      do iz = 0, nz-1
        do iy = 0, ny-1
          do ix = 0, nx-1
            idx = iz*nx*ny + iy*nx + ix + 1
            coord(:,idx) = [ix, iy, iz]
          end do
        end do
      end do

      max_err = 0.0_dp

      do ki = 1, N_cell
        nocc_i = count(blk(ki)%occupation == 2)
        do ka = 1, N_cell
          nvir_a = count(blk(ka)%occupation == 0)
          do kj = 1, N_cell
            nocc_j = count(blk(kj)%occupation == 2)

            kbx = modulo(coord(1,ki) - coord(1,ka) + coord(1,kj), nx)
            kby = modulo(coord(2,ki) - coord(2,ka) + coord(2,kj), ny)
            kbz = modulo(coord(3,ki) - coord(3,ka) + coord(3,kj), nz)
            kb  = kbz*nx*ny + kby*nx + kbx + 1
            nvir_b = count(blk(kb)%occupation == 0)

            do i = 1, nocc_i
              do a = 1, nvir_a
                do j = 1, nocc_j
                  do b = 1, nvir_b
                    ! g1 = (i ki, a ka | j kj, b kb)
                    g1 = ERI_kMO(i, a, j, b, ki, ka, kj)
                    ! g2 = (j kj, b kb | i ki, a ka)
                    g2 = ERI_kMO(j, b, i, a, kj, kb, ki)

                    err = abs(g1 - g2)
                    max_err = max(max_err, err)
                  end do
                end do
              end do
            end do

          end do
        end do
      end do

      write(*,'(/,A)') '=================================================='
      write(*,'(A)')   '         ERI_kMO Physical Symmetry Check          '
      write(*,'(A)')   '=================================================='
      write(*,'(A,ES16.8)') 'Max (ia|jb) - (jb|ia) error : ', max_err
      write(*,'(A,/)') '=================================================='

      deallocate(coord)

end subroutine check_ERI_kMO_symmetry


subroutine check_permutation_defect_propagation(nBas, coeff_SAO, &
                                                ERI_SAO, ERI_kMO)

      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, nx, ny, nz, N_cell
      implicit none

      integer, intent(in) :: nBas

      complex(dpc), intent(in) :: coeff_SAO(nBas,nBas)

      complex(dpc), intent(in) :: &
        ERI_SAO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      complex(dpc), intent(in) :: &
        ERI_kMO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      integer :: a, b, c, d
      integer :: p, q, r, s
      integer :: k1, k2, k3, k4
      integer :: ix, iy, iz, idx
      integer :: k4x, k4y, k4z

      integer, allocatable :: coord(:,:)

      complex(dpc), allocatable :: defect_SAO(:,:,:,:,:,:,:)
      complex(dpc), allocatable :: defect_predicted(:,:,:,:,:,:,:)

      complex(dpc) :: actual_defect
      real(dp) :: max_predicted
      real(dp) :: max_actual
      real(dp) :: max_difference
      real(dp) :: difference
      real(dp) :: relative_difference

      integer :: worst_p, worst_q, worst_r, worst_s
      integer :: worst_k1, worst_k2, worst_k3
      complex(dpc) :: predicted_at_worst
      complex(dpc) :: actual_at_worst

      if (nBas /= nfuc*N_cell) then
        error stop 'check_permutation_defect_propagation: inconsistent nBas'
      end if

      if (N_cell /= nx*ny*nz) then
        error stop 'check_permutation_defect_propagation: inconsistent grid'
      end if

      allocate(coord(3,N_cell))
      allocate(defect_SAO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell))
      allocate(defect_predicted(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell))

      ! Zero-based reciprocal-grid coordinates; x varies fastest.
      do iz = 0, nz-1
        do iy = 0, ny-1
          do ix = 0, nx-1
            idx = iz*nx*ny + iy*nx + ix + 1
            coord(:,idx) = [ix,iy,iz]
          end do
        end do
      end do

      ! --------------------------------------------------------------- !
      ! Build D_SAO = ERI_SAO - electron-pair-permuted ERI_SAO.        !
      ! --------------------------------------------------------------- !
      defect_SAO = cmplx(0.0_dp,0.0_dp,kind=dpc)

      do k1 = 1, N_cell
        do k2 = 1, N_cell
          do k3 = 1, N_cell

            k4x = modulo(coord(1,k1) - coord(1,k2) + coord(1,k3), nx)
            k4y = modulo(coord(2,k1) - coord(2,k2) + coord(2,k3), ny)
            k4z = modulo(coord(3,k1) - coord(3,k2) + coord(3,k3), nz)
            k4  = k4z*nx*ny + k4y*nx + k4x + 1

            do d = 1, nfuc
              do c = 1, nfuc
                do b = 1, nfuc
                  do a = 1, nfuc

                    defect_SAO(a,b,c,d,k1,k2,k3) = &
                      ERI_SAO(a,b,c,d,k1,k2,k3) &
                    - ERI_SAO(c,d,a,b,k3,k4,k1)

                  end do
                end do
              end do
            end do

          end do
        end do
      end do

      ! --------------------------------------------------------------- !
      ! Apply exactly the same SAO -> kMO transformation to D_SAO.    !
      ! This is the predicted MO-space permutation defect.             !
      ! --------------------------------------------------------------- !
      call transform_ERI_SAO_to_kMO(nBas, coeff_SAO, defect_SAO, &
                                    defect_predicted)

      ! --------------------------------------------------------------- !
      ! Construct actual kMO defect and compare to predicted defect.   !
      ! --------------------------------------------------------------- !
      max_predicted  = 0.0_dp
      max_actual     = 0.0_dp
      max_difference = 0.0_dp

      worst_p  = 1
      worst_q  = 1
      worst_r  = 1
      worst_s  = 1
      worst_k1 = 1
      worst_k2 = 1
      worst_k3 = 1

      predicted_at_worst = cmplx(0.0_dp,0.0_dp,kind=dpc)
      actual_at_worst    = cmplx(0.0_dp,0.0_dp,kind=dpc)

      do k1 = 1, N_cell
        do k2 = 1, N_cell
          do k3 = 1, N_cell

            k4x = modulo(coord(1,k1) - coord(1,k2) + coord(1,k3), nx)
            k4y = modulo(coord(2,k1) - coord(2,k2) + coord(2,k3), ny)
            k4z = modulo(coord(3,k1) - coord(3,k2) + coord(3,k3), nz)
            k4  = k4z*nx*ny + k4y*nx + k4x + 1

            do s = 1, nfuc
              do r = 1, nfuc
                do q = 1, nfuc
                  do p = 1, nfuc

                    ! (p k1, q k2 | r k3, s k4)
                    ! - (r k3, s k4 | p k1, q k2)
                    actual_defect = ERI_kMO(p,q,r,s,k1,k2,k3) &
                                  - ERI_kMO(r,s,p,q,k3,k4,k1)

                    difference = abs(actual_defect - &
                                     defect_predicted(p,q,r,s,k1,k2,k3))

                    max_predicted = max(max_predicted, &
                      abs(defect_predicted(p,q,r,s,k1,k2,k3)))

                    max_actual = max(max_actual,abs(actual_defect))

                    if (difference > max_difference) then
                      max_difference = difference

                      worst_p  = p
                      worst_q  = q
                      worst_r  = r
                      worst_s  = s
                      worst_k1 = k1
                      worst_k2 = k2
                      worst_k3 = k3

                      predicted_at_worst = &
                        defect_predicted(p,q,r,s,k1,k2,k3)

                      actual_at_worst = actual_defect
                    end if

                  end do
                end do
              end do
            end do

          end do
        end do
      end do

      relative_difference = max_difference / &
        max(1.0_dp,max(max_predicted,max_actual))

      write(*,'(/,A)') '=================================================='
      write(*,'(A)')   '   SAO -> kMO Permutation-Defect Propagation    '
      write(*,'(A)')   '=================================================='
      write(*,'(A,ES16.8)') 'Max predicted kMO defect : ', max_predicted
      write(*,'(A,ES16.8)') 'Max actual kMO defect    : ', max_actual
      write(*,'(A,ES16.8)') 'Max prediction difference: ', max_difference
      write(*,'(A,ES16.8)') 'Relative difference      : ', relative_difference
      write(*,'(A,4I5)')    'Worst MO indices p,q,r,s : ', &
        worst_p,worst_q,worst_r,worst_s
      write(*,'(A,3I5)')    'Worst k indices k1,k2,k3 : ', &
        worst_k1,worst_k2,worst_k3
      write(*,'(A,2ES20.10)') 'Predicted defect (Re,Im): ', &
        real(predicted_at_worst,dp), aimag(predicted_at_worst)
      write(*,'(A,2ES20.10)') 'Actual defect    (Re,Im): ', &
        real(actual_at_worst,dp), aimag(actual_at_worst)
      write(*,'(A)')   '=================================================='

      deallocate(coord,defect_SAO,defect_predicted)

end subroutine check_permutation_defect_propagation