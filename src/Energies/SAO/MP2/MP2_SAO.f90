subroutine kMP2_calc(N_cell, blk, ERI_kMO, ENuc, EHF)

      use files
      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, nx, ny, nz
      use k_blocks
      implicit none

      ! Inputs
      integer         , intent(in)   :: N_cell
      type(K_block)   , intent(in)   :: blk(N_cell)
      complex(dpc)    , intent(in)   :: ERI_kMO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      real(dp)        , intent(in)   :: ENuc, EHF

      ! Local variables
      integer                        :: ki, ka, kj, kb
      integer                        :: i, a, j, b
      integer                        :: kbx, kby, kbz, idx, ix, iy, iz
      integer, allocatable           :: coord(:,:)
      
      real(dp)                       :: Dijab
      real(dp)                       :: E2d, E2x, EcMP2
      complex(dpc)                   :: g_direct, g_exch

      ! Build coordinate lookup for momentum conservation

      allocate(coord(3, N_cell))
      do iz = 0, nz-1
        do iy = 0, ny-1
          do ix = 0, nx-1
            idx = iz*nx*ny + iy*nx + ix + 1
            coord(:,idx) = [ix, iy, iz]
          end do
        end do
      end do

      E2d = 0.0_dp
      E2x = 0.0_dp

      ! Loop over independent momentum vectors
      do ki = 1, N_cell
        do kj = 1, N_cell
          do ka = 1, N_cell

            ! Momentum conservation: kb = ki - ka + kj  (mod lattice)
            kbx = modulo(coord(1,ki) - coord(1,ka) + coord(1,kj), nx)
            kby = modulo(coord(2,ki) - coord(2,ka) + coord(2,kj), ny)
            kbz = modulo(coord(3,ki) - coord(3,ka) + coord(3,kj), nz)
            kb  = kbz*nx*ny + kby*nx + kbx + 1

            ! Loop over orbitals within each block
            do i = 1, blk(ki)%dim
              if (blk(ki)%occupation(i) /= 2) cycle  ! Occupied only

              do j = 1, blk(kj)%dim
                if (blk(kj)%occupation(j) /= 2) cycle  ! Occupied only

                do a = 1, blk(ka)%dim
                  if (blk(ka)%occupation(a) /= 0) cycle  ! Virtual only

                  do b = 1, blk(kb)%dim
                    if (blk(kb)%occupation(b) /= 0) cycle  ! Virtual only

                    ! Denominator: e_a(ka) + e_b(kb) - e_i(ki) - e_j(kj)
                    Dijab = blk(ka)%energy(a) + blk(kb)%energy(b) &
                          - blk(ki)%energy(i) - blk(kj)%energy(j)

                    ! Direct: (i ki, a ka | j kj, b kb)
                    g_direct = ERI_kMO(i, a, j, b, ki, ka, kj)

                    ! Exchange: (i ki, b kb | j kj, a ka)
                    g_exch   = ERI_kMO(i, b, j, a, ki, kb, kj)

                    ! Direct sum: |g_direct|^2 / Dijab
                    E2d = E2d - real(g_direct * conjg(g_direct), dp) / Dijab

                    ! Exchange sum: Re(g_direct * conjg(g_exch)) / Dijab
                    E2x = E2x - real(g_direct * conjg(g_exch), dp) / Dijab

                  end do
                end do

              end do
            end do

          end do
        end do
      end do

      EcMP2  = 2.0_dp * E2d - E2x
    
      write(outfile,*)
      write(outfile,'(A32)')           '--------------------------'
      write(outfile,'(A32)')           ' K-MP2 calculation        '
      write(outfile,'(A32)')           '--------------------------'
      write(outfile,'(A32,1X,F24.16)') ' MP2 correlation energy = ',EcMP2
      write(outfile,'(A32,1X,F24.16)') ' Direct part            = ',2d0*E2d
      write(outfile,'(A32,1X,F24.16)') ' Exchange part          = ',-E2x
      write(outfile,'(A32)')           '--------------------------'
      write(outfile,'(A32,1X,F24.16)') ' Hartree Fock    energy = ',EHF + ENuc
      write(outfile,'(A32,1X,F24.16)') ' MP2 electronic  energy = ',EcMP2
      write(outfile,'(A32,1X,F24.16)') ' MP2 total       energy = ',ENuc + EHF + EcMP2
      write(outfile,'(A32)')           '--------------------------'
      write(outfile,*)

      deallocate(coord)

end subroutine kMP2_calc
