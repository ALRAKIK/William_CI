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


subroutine symmetry_of_integrals_ERI(nf, fpuc, eri_tmp, eri)

      use unitcell_module
      use omp_lib

      implicit none 

      integer                      :: i, j, k, l
      integer                      :: nf, fpuc
      integer                      :: func_i
      integer                      :: idx, ii, m, n, p
      double precision             :: eri_tmp(fpuc, nf, nf, nf)
      double precision,intent(out) ::     eri(nf, nf, nf, nf)

      ! Precomputation arrays
      integer                      :: icell(nf), jcell(nf), kcell(nf), func_idx(nf)
      integer                      :: equiv(nf, nf)
      integer                      :: eqv(nf)
      integer                      :: el , ek
      
      ! Timing

      !-----------------------------------------------------------------!

      eri(:,:,:,:) = 0.d0
      
      ! Precompute cell and function indices
      do idx = 1, nf
        func_idx(idx) = mod(idx - 1, fpuc) + 1
        ii = (idx - 1) / fpuc
        kcell(idx) = mod(ii, nz)
        ii = ii / nz
        jcell(idx) = mod(ii, ny)
        icell(idx) = ii / ny
      end do
    
      
      eri(1:nfuc,:,:,:) = eri_tmp(:,:,:,:)

      ! Precompute equivalence mapping table
      do i = 1, nf
        do j = 1, nf
          m = modulo(icell(j) - icell(i), nx)
          n = modulo(jcell(j) - jcell(i), ny)
          p = modulo(kcell(j) - kcell(i), nz)
          equiv(i,j) = m * ny * nz * fpuc + n * nz * fpuc + p * fpuc + func_idx(j)
        end do
      end do

      !$OMP PARALLEL DO PRIVATE(i,func_i,eqv,l,el,k,ek,j) SCHEDULE(static)
      do i = fpuc + 1, nf
        func_i = func_idx(i)
        ! Copy equiv for this i into local array
        do j = 1, nf
          eqv(j) = equiv(i, j)
        end do
        do l = 1, nf
          el = eqv(l)
          do k = 1, nf
            ek = eqv(k)
            ! Innermost loop over j (contiguous stride 1)
            do j = 1, nf
              eri(i, j, k, l) = eri(func_i, eqv(j), ek, el)
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO
      !-----------------------------------------------------------------!

end subroutine symmetry_of_integrals_ERI