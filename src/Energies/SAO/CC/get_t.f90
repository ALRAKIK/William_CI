! This is The verified MP2 starting point, NOT a CCSD solver.
! Spatial restricted amplitudes: t1(i,a,ki), t2(i,j,a,b,ki,kj,ka).
! Bare chemist-order ERIs. Do not antisymmetrize the driving term.
! All energies below are total supercell correlation energies.

subroutine get_t_term(blk,g,t2)

      use constants_module, only: dp, dpc
      use unitcell_module , only: nfuc, N_cell
      use k_blocks        , only: K_block

      implicit none 

      ! input ! 

      type(K_block), intent(in)  :: blk(N_cell)
      complex(dpc) , intent(in)  :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      ! local ! 

      integer                    :: ki,kj,ka,kb
      integer                    :: i,j,a,b
      integer                    :: k,nocc
      integer                    :: fourth_k
      real(dp)                   :: gap
      complex(dpc)               :: gd

      ! output ! 

      complex(dpc), intent(out)  :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      nocc = 0
      do k=1,N_cell
        nocc = nocc + count(blk(k)%occupation == 2)
      end do

      t2 = cmplx(0.0_dp,0.0_dp,kind=dpc)

      do ka=1,N_cell
        do kj=1,N_cell
          do ki=1,N_cell
            kb = fourth_k(ki,kj,ka)
            do b=1,nfuc
              if (blk(kb)%occupation(b) /= 0) cycle
              do a=1,nfuc
                if (blk(ka)%occupation(a) /= 0) cycle
                do j=1,nfuc
                  if (blk(kj)%occupation(j) /= 2) cycle
                  do i=1,nfuc
                    if (blk(ki)%occupation(i) /= 2) cycle
                    gap = blk(ka)%energy(a)+blk(kb)%energy(b) &
                        - blk(ki)%energy(i)-blk(kj)%energy(j)
                    if (gap <= 1.0e-12_dp) error stop 'CC baseline: nonpositive excitation gap'
                    gd = g(i,a,j,b,ki,ka,kj)
                    t2(i,j,a,b,ki,kj,ka) = -conjg(gd)/gap
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      
end subroutine get_t_term