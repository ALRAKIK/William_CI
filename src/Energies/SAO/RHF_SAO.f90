subroutine RHF_SAO(nBas,nO,S,T,V,Hc,ERI,X,ENuc,EHF,e,coeff_SAO,blk,natom,Zat)

      ! Perform a restricted Hartree-Fock calculation

      use files 
      use keywords
      use eri_direct, only : direct_JK          ! [Nov 2026] keyword direct
      use df_toroidal, only : df_JK             ! [Nov 2026] keyword density_fitting
      use constants_module
      use unitcell_module
      use k_blocks
      use debug_tools
      use HF_parameter
      use cc_dims , only : use_ksym      ! [Nov 2026]

      implicit none

      ! Input variables

      integer,intent(in)            ::                    nBas
  
      integer,intent(in)            ::                      nO
      double precision,intent(in)   ::             S(nBas,nBas)
      double precision,intent(in)   ::             T(nBas,nBas)
      double precision,intent(in)   ::             V(nBas,nBas)
      double precision,intent(in)   ::            Hc(nBas,nBas) 
      double precision,intent(in)   ::             X(nBas,nBas) 
      double precision,intent(in)   :: ERI(nfuc,nBas,nBas,nBas)

      double precision,intent(in)   ::                     ENuc
  
      ! Local variables
  
      
      integer                       ::                     nSCF
      double precision              ::                     Conv
      double precision              ::                      Gap
      double precision              ::                    ET,EV
      double precision              ::                    EJ,EK
      double precision              ::                  EHF_old
      double precision,allocatable  ::                   P(:,:)
      double precision,allocatable  ::                   J(:,:)
      double precision,allocatable  ::                   K(:,:)
      double precision,allocatable  ::                   F(:,:)
      double precision,allocatable  ::               error(:,:)
      double precision,allocatable  ::             error_o(:,:)   ! [Nov 2026] X^T error X
      double precision,allocatable  ::          wA(:,:), wB(:,:)  ! [Nov 2026] work for the dgemm products
      logical                       ::               converged     ! [Nov 2026]
      logical                       ::               ksym_ok       ! [Nov 2026] fast_ccsdt gauge
      integer                       ::               fourth_k      ! [Nov 2026]
      integer                       ::        kb , b1 , b2 , nf    ! [Nov 2026] k-block bounds
      double precision,external     ::             trace_matrix

      integer                       ::                   n_diis
      integer                       ::                   i  , o
      integer                       ::                  mu , nu
      double precision              ::                      sum 
      double precision              ::              sum_o(nbas)
      double precision              ::                    rcond
      double precision,allocatable  ::            err_diis(:,:)
      double precision,allocatable  ::              F_diis(:,:)
  
      ! Output variables
  
      double precision,intent(out)  ::                      EHF 
      double precision,intent(out)  ::                  e(nBas)
      complex(dpc)    ,intent(out)  ::     coeff_SAO(nBas,nBas)


      ! ----------------------    SAO     ----------------------------- !

      complex(dpc),allocatable      ::              S_SAO(:,:)
      complex(dpc),allocatable      ::              X_SAO(:,:)

      complex(dpc),allocatable      ::              F_SAO(:,:)
      complex(dpc),allocatable      ::             Fp_SAO(:,:)
      complex(dpc),allocatable      ::              P_SAO(:,:)


      integer,allocatable           ::                  idx(:) 
      integer                       ::            tmp, occ_col
      double precision,allocatable  ::             e_sorted(:)

      type(K_block)                 ::             blk(N_cell)

      integer         , intent(in)  :: natom
      integer         , intent(in)  :: Zat(natom)
      
      !-----------------------------------------------------------------!

      if (c_details) then 
        open(HFfile,file=trim(tmp_file_name)//"/RHF.out")
      end if 

      write(outfile,*)
      write(outfile,*)'******************************************************************************************'
      write(outfile,*)'|                          Restricted Hartree-Fock calculation                           |'
      write(outfile,*)'******************************************************************************************'
      write(outfile,*)
  
      ! Memory allocation
  
      allocate(P(nBas,nBas) ,error(nBas,nBas), error_o(nBas,nBas))
      allocate(wA(nBas,nBas), wB(nBas,nBas))
      allocate(J(nBas,nBas) ,    K(nBas,nBas))
      allocate(F(nBas,nBas))

      allocate(err_diis(nBas*nBas,max_diis))
      allocate(  F_diis(nBas*nBas,max_diis))

      allocate(S_SAO(nBas,nBas) ,      X_SAO(nBas,nBas))
      allocate(P_SAO(nBas,nBas))
      allocate(F_SAO(nBas,nBas) ,     Fp_SAO(nBas,nBas))

      allocate(idx(nBas) ,   e_sorted(nBas))
       
      if (c_Huckel) then 
        call header_under("Huckel guess",-1)
        call guess_Huckel_RHF_SAO(nBas,c_details,nO,Hc,X,ENuc,S,T,V,P)
      else 
        call header_under("AO guess",-1)
        call guess_AO_RHF_SAO(nBas,c_details,nO,HC,X,ENuc,S,T,V,P)
      end if 
         
      ! --------------------------------------------------------------- !
      !              check that P_{mu nu} S_{mu nu} = N 
      ! --------------------------------------------------------------- !

      sum = 0.d0 
      do mu = 1 , nbas 
        do nu = 1 , nbas 
          sum = sum + P(mu,nu)*S(mu,nu)
        end do 
      end do 

      call header_under("P_{mu nu} S_{mu nu} = N",-1)

      write(outfile,"(a,f16.8)") "P_{mu nu} S_{mu nu} =  " , sum 

      ! --------------------------------------------------------------- !

      ! --------------------------------------------------------------- !
      !              apply the density on the functions 
      ! [Nov 2026] printed only with the keyword Details (nBas lines)
      ! --------------------------------------------------------------- !

      if (c_details) then
        sum_o(:) = 0.d0 

        do o = 1 , nbas
            sum_o(o) = dot_product(S(:,o), matmul(P, S(:,o)))
        end do 

        call header_under("P_{mu nu} S_{mu s} S_{nu s}",-1)

        write(outfile,"(a)") "P_{mu nu} S_{mu s} S_{nu s} :  "
        do i = 1 , nbas 
          write(outfile,"(40x,i3,1x,f16.8)") i , sum_o(i)
        end do 
      end if

      ! --------------------------------------------------------------- !

      call transform_AO_to_SAO(nBas,S,S_SAO)

      ! [Nov 2026] S_SAO is block diagonal in k: Lowdin S^-1/2 computed block
      ! by block (same X, cost N_cell*nfuc^3 instead of nBas^3)
      ! (was: call get_X_from_overlap_complex(nBas,S_SAO,X_SAO))
      call get_X_from_overlap_complex_blocks(nBas,nfuc,N_cell,S_SAO,X_SAO)
      
      ! --------------------------------------------------------------- !

      ! Initialization
  
      nSCF = 0
      Conv = 1d0

      if (c_DIIS) then 
        n_diis        = 0
        F_diis(:,:)   = 0d0
        err_diis(:,:) = 0d0
        rcond         = 0d0
      end if 

      EHF_old         = 0d0 
      EHF             = 0d0 
      converged       = .false.
  
      !------------------------------------------------------------------------
      ! Main SCF loop
      !------------------------------------------------------------------------
  
      write(outfile,*)
      write(outfile,*) repeat('-', 110)
      write(outfile,*) "|",repeat(' ', 47),"RHF calculation",repeat(' ', 46),"|"
      write(outfile,*) repeat('-', 110)

      write(outfile,'(1x,a1,a2,1x,a1,a,a1,a,a1,a,a1,a,a1,a,a1,a,a1)') &
      "|","#","|","      HF energy  ","|", "    Conv   ","|","   HL Gap  ","|", "    T contribution  ",&
      "|","   V contribution   ","|","   Two contribution ","|"
      write(outfile,*) repeat('-', 110)
      flush(outfile)

      write(*,*)
      write(*,*) repeat('-', 46)
      write(*,*) "|",repeat(' ', 15),"RHF calculation",repeat(' ', 15),"|"
      write(*,*) repeat('-', 46)

      write(*,'(1x,a1,a2,1x,a1,a,a1,a,a1,a,a1,a,a1,a,a1,a,a1)') &
      "|","#","|","      HF energy  ","|", "    Conv   ","|","   HL Gap  ","|"
      write(*,*) repeat('-', 46)
      flush(6)
      
      do while(nSCF < maxSCF)
  
      !   Increment 
  
      nSCF = nSCF + 1
  
      !   [Nov 2026] Coulomb J and exchange K (rows mu <= nfuc) in one
      !   parallel pass (was: hartree_potential_SAO + exchange_potential_SAO)
      if (c_df) then
        call df_JK(nBas,P,J,K)                      ! [Nov 2026] keyword density_fitting: fitted J and K
      else if (c_direct) then
        call direct_JK(nBas,P,J,K)                  ! [Nov 2026] keyword direct: integrals recomputed
      else
      call hartree_exchange_potential(nfuc,nBas,P,ERI,J,K)
      end if
      ! ****************** !

      !   Build Fock operator
    
      F(:,:) = Hc(:,:) + J(:,:) + K(:,:)

      if (c_debug) call check_scf_state(nBas, F, P, S, X, outfile)

      call transform_AO_to_SAO(nBas,F,F_SAO)
      call transform_SAO_to_AO(nBas,F_SAO,F)

      !---------------------------------------------------------------- !
      !   Compute HF energy
      !   [Nov 2026] computed BEFORE the convergence test, so the final
      !   energy belongs to the final density (before: one iteration late)
      !---------------------------------------------------------------- !
  
      !   [Nov 2026] Tr(P A) = sum_ij P(i,j) A(j,i): O(n^2), no matrix product
      !   (was: trace_matrix(nBas,matmul(P,A)), O(n^3) each)

      !   Compute the kinetic energy
      ET = dot_product(reshape(P,[nBas*nBas]), reshape(transpose(T),[nBas*nBas]))

      !   Compute the potential energy
      EV = dot_product(reshape(P,[nBas*nBas]), reshape(transpose(V),[nBas*nBas]))

      !   Compute the Hartree energy
      EJ = 0.5d0*dot_product(reshape(P,[nBas*nBas]), reshape(transpose(J),[nBas*nBas])) * N_cell

      !   Compute the exchange energy
      EK = 0.5d0*dot_product(reshape(P,[nBas*nBas]), reshape(transpose(K),[nBas*nBas])) * N_cell

      ! ****************** !

      !   Total HF energy
  
      EHF = ET + EV + EJ + EK

      ! Compute the error vector and extract the convergence criterion  !
      ! [Nov 2026] error in the orthogonal basis, X^T (FPS - SPF) X : used
      ! for the convergence test AND for DIIS (standard, better conditioned)

      ! [Nov 2026] same products with BLAS dgemm (was: matmul, ~10x slower)
      call dgemm('N','N',nBas,nBas,nBas,1.d0,P,nBas,S,nBas,0.d0,wA,nBas)        ! P S
      call dgemm('N','N',nBas,nBas,nBas,1.d0,F,nBas,wA,nBas,0.d0,error,nBas)    ! F P S
      call dgemm('N','N',nBas,nBas,nBas,1.d0,S,nBas,P,nBas,0.d0,wA,nBas)        ! S P
      call dgemm('N','N',nBas,nBas,nBas,-1.d0,wA,nBas,F,nBas,1.d0,error,nBas)   ! - S P F
      call dgemm('N','N',nBas,nBas,nBas,1.d0,error,nBas,X,nBas,0.d0,wB,nBas)    ! error X
      call dgemm('T','N',nBas,nBas,nBas,1.d0,X,nBas,wB,nBas,0.d0,error_o,nBas)  ! X^T error X
      if(nSCF > 1) Conv = maxval(abs(error_o))

      ! [Nov 2026] stop ONLY when the orbitals are converged. Before, the
      ! loop also stopped when |EHF - EHF_old| < thresh (energy converges
      ! much faster than the orbitals -> MP2/CC errors 1e-8 .. 1e-7)
      if (Conv < thresh) then
        converged = .true.
        exit
      end if

      EHF_old = EHF

      ! ****************** !

      if (c_details) call details_print_HF(nBas,nSCF,J,K,F,1)

      ! DIIS extrapolation ! 

      if (c_DIIS) then 
        if(max_diis > 1) then
          n_diis = min(n_diis+1,max_diis)
          call DIIS_extrapolation(rcond,nBas*nBas,nBas*nBas,n_diis,err_diis,F_diis,error_o,F)   ! [Nov 2026] error_o
        end if
      end if

      ! [Nov 2026] BUG FIX: F_SAO must be rebuilt from the DIIS-extrapolated F.
      ! Before, F_SAO was the one built before DIIS, so DIIS had no effect at
      ! all in the SAO path (identical runs with and without the keyword DIIS).
      call transform_AO_to_SAO(nBas,F,F_SAO)

      !   Transform for the Fock matrix F in the orthogonal basis 
      !   [Nov 2026] block by block (F_SAO and X_SAO are block diagonal in k):
      !   N_cell products of nfuc x nfuc instead of one nBas x nBas product
      !   (was: Fp_SAO = matmul(conjg(transpose(X_SAO)), matmul(F_SAO, X_SAO)))

      Fp_SAO = (0.0_dp, 0.0_dp)
      do kb = 1, N_cell
        b1 = (kb-1)*nfuc + 1 ; b2 = kb*nfuc
        Fp_SAO(b1:b2,b1:b2) = matmul(conjg(transpose(X_SAO(b1:b2,b1:b2))), &
                                     matmul(F_SAO(b1:b2,b1:b2), X_SAO(b1:b2,b1:b2)))
      end do
            
      !   Diagonalize F' to get MO coefficients (eigenvectors in the orthogonal basis) c' and MO energies (eigenvalues) e
      
      ! -------------------- Full diagonalization --------------------- !

      !call diagonalize_matrix_complex(nBas,Fp_SAO,e)

      ! ---------------- Block by Block Diagonalization --------------- !

      call diagonalize_block_diagonal_complex_para(nBas, nfuc, N_cell, Fp_SAO, e)

      do i = 1, nBas
        idx(i) = i
      end do

      do mu = 1, nBas-1
        do nu = mu+1, nBas
          if (e(idx(mu)) > e(idx(nu))) then
            tmp = idx(mu)
            idx(mu) = idx(nu)
            idx(nu) = tmp
          end if
        end do
      end do      

      !   Back-transform the MO coefficients c in the original non-orthogonal basis

      !   [Nov 2026] block by block (was: coeff_SAO = matmul(X_SAO, Fp_SAO))
      coeff_SAO = (0.0_dp, 0.0_dp)
      do kb = 1, N_cell
        b1 = (kb-1)*nfuc + 1 ; b2 = kb*nfuc
        coeff_SAO(b1:b2,b1:b2) = matmul(X_SAO(b1:b2,b1:b2), Fp_SAO(b1:b2,b1:b2))
      end do

      !   Compute the density matrix P  !
      !   [Nov 2026] an orbital of block kb has coefficients only in rows of
      !   block kb, so its contribution 2 c c^* is limited to that block
      !   (was: a loop over all nBas x nBas elements for every orbital)

      P_SAO = (0.0_dp, 0.0_dp)
      nf    = nfuc

      do i = 1, nO
        occ_col = idx(i)
        kb = (occ_col-1)/nf + 1
        b1 = (kb-1)*nf + 1 ; b2 = kb*nf
        do nu = b1, b2
          do mu = b1, b2
            P_SAO(mu, nu) = P_SAO(mu, nu) + 2.0_dp * coeff_SAO(mu, occ_col) * conjg(coeff_SAO(nu, occ_col))
          end do
        end do
      end do

      call transform_SAO_to_AO(nBas,P_SAO,P)

      !   Compute HOMO-LUMO gap

      if (nBas > nO) then
        Gap = e(idx(nO+1)) - e(idx(nO))
      else
        Gap = 0d0
      endif

      if (c_details) then 

      call header_HF("The Real orbital coeff c ", -1)
      write(HFfile,'(17x,*(i3,15x))') (i,i=1,nBas)
      do i = 1 , nBas
        write(HFfile,'(i3,6x,*(f16.10,2x))') i ,  (real(coeff_SAO(o,i)),o=1,nBas)
      end do 
      write(HFfile,'(a)') ""

      call header_HF("The imag orbital coeff c ", -1)
      write(HFfile,'(17x,*(i3,15x))') (i,i=1,nBas)
      do i = 1 , nBas
        write(HFfile,'(i3,6x,*(f16.10,2x))') i ,  (aimag(coeff_SAO(o,i)),o=1,nBas)
      end do 
      write(HFfile,'(a)') ""
      
      call header_HF("Density matrix P = 2 c c^t", -1)
      write(HFfile,'(17x,*(i3,15x))') (i,i=1,nfuc)
      do i = 1 , size(P,1)
        write(HFfile,'(i3,6x,*(f16.10,2x))') i ,  (P(o,i),o=1,nfuc)
      end do 
      write(HFfile,'(a)') ""

      call header_HF("HF Orbital Energies", -1)
      do i = 1 , size(e)
        write(HFfile,'(i3,6x,*(f16.10,2x))') i ,  e(i)
      end do

      write(HFfile,'(a)') ""
      write(HFfile,'(a,f16.10)')   " The Kinetic   Energy    = ", ET
      write(HFfile,'(a,f16.10)')   " The Potential Energy    = ", EV
      write(HFfile,'(a,f16.10)')   " The Hartree   Energy    = ", EJ
      write(HFfile,'(a,f16.10,a)') " The Exchange  Energy    = ", EK  , "    +"
      write(HFfile,'(a,f16.10,a)') " The Nuclear   Energy    = ", ENuc 
      write(HFfile,'(a,f16.10)')   "-----------------------------------"
      write(HFfile,'(a,f16.10)')   " The HF Energy        = ", EHF+ENuc
      write(HFfile,'(a)') ""
      write(HFfile,'(a)') "" 
      
      write(HFfile,'(2a)')  repeat('*_',36) , "*"
      write(HFfile,'(a)')   repeat('_',73)

      end if 

      !   Dump results
   
      write(outfile,"(1x,a1,i2,1x,a1,f16.8,1x,a1,f10.6,1x,a1,f10.6,1x,a1,f16.8,4x,a1,f16.8,4x,a1,f16.8,4x,a1)")      & 
      "|",nSCF,"|",EHF+ENuc,"|",Conv,"|",Gap,"|",ET,"|",EV,"|",EJ+EK,"|"
      write(outfile,*) repeat('-', 110)
      flush(outfile)

      write(*,"(1x,a1,i2,1x,a1,f16.8,1x,a1,f10.8,1x,a1,f10.6,1x,a1,f16.8,4x,a1,f16.8,4x,a1,f16.8,4x,a1)")      & 
      "|",nSCF,"|",EHF+ENuc,"|",Conv,"|",Gap,"|"
      write(*,*) repeat('-', 46)
      flush(6)
      
      enddo

      if (c_details) then 
        close(HFfile)
      end if 

      call print_band_structure(nBas, nfuc, N_cell, e, 'band_structure.dat')

      ! [Nov 2026] fast_ccsdt: time-reversal gauge C(-k) = C(k)^* (ksym_gauge.f90)
      !            (same density and energy; lets the fast CCSD(T) skip half of the k blocks)
      ksym_ok = .false.
      if (c_fast_ccsdt) call ksym_gauge(nBas, nfuc, N_cell, S_SAO, F_SAO, coeff_SAO, e, ksym_ok)

      call build_k_blocks(nBas,nO, N_cell, coeff_SAO, e, blk)

      ! [Nov 2026] the symmetry also needs the same occupations at k and -k
      use_ksym = .false.
      if (c_fast_ccsdt .and. ksym_ok) then
        use_ksym = .true.
        do kb = 1, N_cell
          if (any(blk(kb)%occupation /= blk(fourth_k(1,1,kb))%occupation)) use_ksym = .false.
        end do
        if (.not. use_ksym) write(outfile,'(a)') &
          '  [fast_ccsdt] WARNING: occupations differ at k and -k, k <-> -k symmetry not used'
      end if
      
      !------------------------------------------------------------------------
      ! End of SCF loop
      !------------------------------------------------------------------------
  
      if(.not. converged) then          ! [Nov 2026] was nSCF == maxSCF
  
        write(outfile,*)
        write(outfile,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(outfile,*)'                 Convergence failed                 '
        write(outfile,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(outfile,*)

        call print_RHF_SAO(nBas,nO,e,coeff_SAO,ENuc,ET,EV,EJ,EK,EHF,blk)

        deallocate(P,error,error_o,J,K,F,err_diis,F_diis,S_SAO,X_SAO,P_SAO,F_SAO,Fp_SAO,idx,e_sorted)
  
        stop
  
      endif
  
        ! Compute final HF energy

        call print_RHF_SAO(nBas,nO,e,coeff_SAO,ENuc,ET,EV,EJ,EK,EHF,blk)

        call print_mulliken_charges(nBas, natom, S, P, Zat)   ! [Nov 2026] real number of functions per atom

        ! Memory deallocation

        deallocate(P,error,error_o,J,K,F,err_diis,F_diis,S_SAO,X_SAO,P_SAO,F_SAO,Fp_SAO,idx,e_sorted)

        write(*,'(a)') ""
        write(*,'(a)') "╔───────────────────────────────────────────────╗"
        write(*,'(a)') "│             Hartree-Fock completed            │"
        write(*,'(a)') "╚───────────────────────────────────────────────╝"
        flush(6)
        write(*,'(a)') ""
        
end subroutine RHF_SAO