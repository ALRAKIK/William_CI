subroutine k_ccsd_T(blk,ERI_kMO,EHF,Enuc,EMP2)

      use files
      use constants_module, only: dp, dpc
      use unitcell_module , only: nfuc, N_cell
      use k_blocks        , only: K_block

      implicit none 

      ! input !

      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: ERI_kMO(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      real(dp)     , intent(in) :: EHF , Enuc , EMP2

      ! local ! 

      integer ,parameter        :: max_iter      = 100
      real(dp),parameter        :: damping       = 0.7_dp
      real(dp),parameter        :: energy_tol    = 1.0e-10_dp
      real(dp),parameter        :: residual_tol  = 1.0e-9_dp
      integer ,parameter        :: ndiis         = 6
      
      complex(dpc),allocatable  :: t1(:,:,:)
      complex(dpc),allocatable  :: t2(:,:,:,:,:,:,:)

      logical                   :: cc_converged

      ! output ! 

      real(dp)                  :: E_corr , E_corrt 


      write(*,'(a)') ""
      write(*,'(a)') "╔───────────────────────────────────────────────╗"
      write(*,'(a)') "│               CCSD(T) Started                 │"
      write(*,'(a)') "╚───────────────────────────────────────────────╝"
      flush(6)
      write(*,'(a)') ""


      allocate(t1(nfuc,nfuc,N_cell))
      allocate(t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell))

      t1 = cmplx(0.0_dp,0.0_dp,kind=dpc)

      call get_t_term(blk,ERI_kMO,t2)

      call solve_kccsd(blk, ERI_kMO, t1, t2, cc_converged, max_iter, damping, energy_tol, residual_tol, ndiis,E_corr)

      if (cc_converged) then
         call build_ccsd_t(blk, ERI_kMO, t1, t2,E_corrt)
      end if

      write(outfile,*)
      write(outfile,'(4x,A32)')        '------------------------------'
      write(outfile,'(4x,A32)')           ' K-CCSD(T) calculation        '
      write(outfile,'(4x,A32)')           '------------------------------'
      write(outfile,'(4x,A32,1X,F24.16)') ' Hartree Fock       energy = ',EHF + ENuc
      write(outfile,'(4x,A32,1X,F24.16)') ' MP2 correlation    energy = ',EMP2
      write(outfile,'(4x,A32)')           '------------------------------'
      write(outfile,'(4x,A32,1X,F24.16)') ' CCSD electronic    energy = ',E_corr
      write(outfile,'(4x,A32,1X,F24.16)') ' HF + CCSD          energy = ',ENuc + EHF + E_corr
      write(outfile,'(4x,A32,1X,F24.16)') ' CCSD(T) electronic energy = ',E_corrt
      write(outfile,'(4x,A32,1X,F24.16)') ' HF + CCSD(T)       energy = ',ENuc + EHF + E_corr + E_corrt
      write(outfile,'(4x,A32)')           '------------------------------'
      write(outfile,*)


      write(*,'(a)') ""
      write(*,'(a)') "╔───────────────────────────────────────────────╗"
      write(*,'(a)') "│               CCSD(T) completed               │"
      write(*,'(a)') "╚───────────────────────────────────────────────╝"
      flush(6)
      write(*,'(a)') ""



end subroutine