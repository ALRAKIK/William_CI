subroutine RHF_SAO(nBas,nO,S,T,V,Hc,ERI,X,ENuc,EHF,e,coeff_SAO)

      ! Perform a restricted Hartree-Fock calculation

      use files 
      use keywords
      use constants_module
      use unitcell_module
      use k_blocks

      implicit none

      ! Input variables

      integer,intent(in)            :: nBas
  
      integer,intent(in)            ::                      nO
      double precision,intent(in)   ::             S(nBas,nBas)
      double precision,intent(in)   ::             T(nBas,nBas)
      double precision,intent(in)   ::             V(nBas,nBas)
      double precision,intent(in)   ::            Hc(nBas,nBas) 
      double precision,intent(in)   ::             X(nBas,nBas) 
      double precision,intent(in)   :: ERI(nfuc,nBas,nBas,nBas)

      double precision,intent(in)   :: ENuc
  
      ! Local variables
  
      integer,parameter             :: maxSCF = 100
      double precision,parameter    :: thresh = 1d-8
      integer                       :: nSCF
      double precision              :: Conv
      double precision              :: Gap
      double precision              :: ET,EV,EJ , EHF_old
      double precision              :: EK
      double precision,allocatable  :: cp(:,:)
      double precision,allocatable  :: P(:,:)
      double precision,allocatable  :: J(:,:)
      double precision,allocatable  :: K(:,:)
      double precision,allocatable  :: F(:,:),Fp(:,:)
      double precision,allocatable  :: error(:,:)
      double precision,external     :: trace_matrix

      integer                       :: max_diis = 5
      integer                       :: n_diis
      integer                       :: i  , o
      integer                       :: mu , nu 
      double precision              :: sum 
      double precision              :: sum_o(nbas)
      double precision              :: rcond
      double precision,allocatable  :: err_diis(:,:)
      double precision,allocatable  :: F_diis(:,:)
  
      ! Output variables
  
      double precision,intent(out)  :: EHF 
      double precision,intent(out)  :: e(nBas)
      complex(dpc)    ,intent(out)  :: coeff_SAO(nBas,nBas)


      ! ----------------------    SAO     ----------------------------- !

      complex(dpc)                  ::     S_SAO(nBas,nBas)
      complex(dpc)                  ::     X_SAO(nBas,nBas)

      complex(dpc)                  ::     F_SAO(nBas,nBas)
      complex(dpc)                  ::    Fp_SAO(nBas,nBas)
      complex(dpc)                  ::     P_SAO(nBas,nBas)
      complex(dpc)                  ::  P_SAO_AO(nBas,nBas)

      integer                       :: idx(nBas) , tmp , occ_col
      double precision              :: e_sorted(nBas)

      type(K_block), allocatable    :: blk(:)

      
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
  
      allocate(cp(nBas,nBas), P(nBas,nBas) ,error(nBas,nBas))
      allocate(J(nBas,nBas) , K(nBas,nBas))
      allocate(F(nBas,nBas) ,Fp(nBas,nBas))

      allocate(err_diis(nBas*nBas,max_diis))
      allocate(  F_diis(nBas*nBas,max_diis))
       
      if (c_Huckel) then 
        call header_under("Huckel guess",-1)
        call guess_Huckel_RHF_SAO(nBas,c_details,nO,Hc,X,ENuc,S,T,V,P,ERI)
      else 
        call header_under("AO guess",-1)
        call guess_AO_RHF_SAO(nBas,c_details,nO,HC,X,ENuc,S,T,V,P,ERI)
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
      ! --------------------------------------------------------------- !

      sum_o(:) = 0.d0 

      do o = 1 , nbas
        do mu = 1 , nbas 
          do nu = 1 , nbas 
            sum_o(o) = sum_o(o) + P(mu,nu)*S(mu,o)*S(nu,o)
          end do 
        end do 
      end do 

      call header_under("P_{mu nu} S_{mu s} S_{nu s}",-1)

      write(outfile,"(a)") "P_{mu nu} S_{mu s} S_{nu s} :  "
      do i = 1 , nbas 
        write(outfile,"(40x,i3,1x,f16.8)") i , sum_o(i)
      end do 

      ! --------------------------------------------------------------- !

      call transform_AO_to_SAO(nBas,S,S_SAO)

      call get_X_from_overlap_complex(nBas,S_SAO,X_SAO)

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
      
      do while(nSCF < maxSCF)
  
      !   Increment 
  
      nSCF = nSCF + 1
  
      !   Compute the Hartree potential J
      call hartree_potential_SAO(nBas,P,ERI,J)
      ! ****************** !
  
      !   Compute the exchange potential K
      call exchange_potential_SAO(nBas,P,ERI,K)
      ! ****************** !

      !   Build Fock operator
    
      F(:,:) = Hc(:,:) + J(:,:) + K(:,:)

      call transform_AO_to_SAO(nBas,F,F_SAO)

      !   Compute the error vector and extract the convergence criterion

      error = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)
      if(nSCF > 1) Conv = maxval(abs(error))
      if (Conv < thresh) exit 

      ! ****************** !

      if (c_details) then 

      write(HFfile,'(a)') ""
      write(HFfile,'(3a)') "!" ,repeat('-',44), "!"
      write(HFfile,'(a,I3)') "                 Iter  = ", nSCF
      write(HFfile,'(3a)') "!" ,repeat('-',44), "!"
      write(HFfile,'(a)') ""

    
      call header_HF("Hartree potential J = P(la,si)*ERI(mu,nu,la,si)",-1)
      write(HFfile,'(17x,1000(i3,15x))') (i,i=1,nfuc)
      do i = 1 , size(j,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (j(o,i),o=1,nfuc)
      end do 
      write(HFfile,'(a)') ""

      call header_HF("Exchange potential K = -0.5 P(la,si)*ERI(mu,la,si,nu)",-1)
      write(HFfile,'(17x,1000(i3,15x))') (i,i=1,nfuc)
      do i = 1 , size(k,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (k(o,i),o=1,nfuc)
      end do 
      write(HFfile,'(a)') ""

      call header_HF("Fock matrix F = Hc + J + K", -1)
      write(HFfile,'(17x,1000(i3,15x))') (i,i=1,nfuc)
      do i = 1 , size(F,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (F(o,i),o=1,nfuc)
      end do 
      write(HFfile,'(a)') ""

      end if 

      !------------------------------------------------------------------------
      !   Compute HF energy
      !------------------------------------------------------------------------
  
      !   Compute the kinetic energy
      ET = trace_matrix(nBas,matmul(P,T))
      ! ****************** !

      !   Compute the potential energy
      EV = trace_matrix(nBas,matmul(P,V))
      ! ****************** !

      !   Compute the Hartree energy
      EJ = 0.5d0*trace_matrix(nBas,matmul(P,J)) * N_cell

      ! ****************** !

      !   Compute the exchange energy
      EK = 0.5d0*trace_matrix(nBas,matmul(P,K)) * N_cell

      ! ****************** !

      !   Total HF energy
  
      EHF = ET + EV + EJ + EK

      if ((abs(EHF - EHF_old) < thresh ) .and. nSCF > 2 ) exit

      if (nSCF > 2) then 
        EHF_old = EHF 
      end if 

      !   Compute HOMO-LUMO gap

      if (nBas > nO) then
        e_sorted = e
        call sort_array(e_sorted, nBas)
        Gap = e_sorted(nO+1) - e_sorted(nO)
      else
        Gap = 0.0_dp
      endif

      ! ****************** !

      ! DIIS extrapolation ! 

      if (c_DIIS) then 
        if(max_diis > 1) then
          n_diis = min(n_diis+1,max_diis)
          call DIIS_extrapolation(rcond,nBas*nBas,nBas*nBas,n_diis,err_diis,F_diis,error,F)
        end if
      end if

      !   Transform for the Fock matrix F in the orthogonal basis 
      
      Fp_SAO = matmul(conjg(transpose(X_SAO)), matmul(F_SAO, X_SAO)) 

            
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

      coeff_SAO = matmul(X_SAO, Fp_SAO)

      !   Compute the density matrix P  !

      P_SAO = (0.0_dp, 0.0_dp)

      do i = 1, nO
        occ_col = idx(i)
        do mu = 1, nBas
          do nu = 1, nBas
            P_SAO(mu, nu) = P_SAO(mu, nu) + 2.0_dp * coeff_SAO(mu, occ_col) * conjg(coeff_SAO(nu, occ_col))
          end do
        end do
      end do

      call transform_SAO_to_AO(nBas,P_SAO,P_SAO_AO)

      P(:,:) = dble(P_SAO_AO(:,:))

      if (c_details) then 

      call header_HF("The Real orbital coeff c ", -1)
      write(HFfile,'(17x,1000(i3,15x))') (i,i=1,nBas)
      do i = 1 , nBas
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (real(coeff_SAO(o,i)),o=1,nBas)
      end do 
      write(HFfile,'(a)') ""

      call header_HF("The imag orbital coeff c ", -1)
      write(HFfile,'(17x,1000(i3,15x))') (i,i=1,nBas)
      do i = 1 , nBas
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (aimag(coeff_SAO(o,i)),o=1,nBas)
      end do 
      write(HFfile,'(a)') ""
      
      call header_HF("Density matrix P = 2 c c^t", -1)
      write(HFfile,'(17x,1000(i3,15x))') (i,i=1,nfuc)
      do i = 1 , size(P,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (P(o,i),o=1,nfuc)
      end do 
      write(HFfile,'(a)') ""

      call header_HF("HF Orbital Energies", -1)
      do i = 1 , size(e)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  e(i)
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
      
      enddo

      if (c_details) then 
        close(HFfile)
      end if 

      call print_band_structure(nBas, nfuc, N_cell, e, 'band_structure.dat')

      call build_k_blocks(nBas,nO, N_cell, coeff_SAO, e, blk)

      !------------------------------------------------------------------------
      ! End of SCF loop
      !------------------------------------------------------------------------
  
      if(nSCF == maxSCF) then
  
        write(outfile,*)
        write(outfile,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(outfile,*)'                 Convergence failed                 '
        write(outfile,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(outfile,*)

        call print_RHF_SAO(nBas,nO,e,coeff_SAO,ENuc,ET,EV,EJ,EK,EHF,blk)
  
        stop
  
      endif
  
        ! Compute final HF energy

        call print_RHF_SAO(nBas,nO,e,coeff_SAO,ENuc,ET,EV,EJ,EK,EHF,blk)

        
end subroutine RHF_SAO