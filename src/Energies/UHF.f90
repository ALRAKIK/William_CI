subroutine UHF(nBas,c_details,n_alpha,n_beta,S,T,V,Hc,ERI,X,ENuc,EHF,e,c)

      ! Perform a restricted Hartree-Fock calculation

      use files 
      implicit none

      ! Input variables

      integer         ,intent(in)   :: nBas
  
      integer         ,intent(in)   :: n_alpha , n_beta
      double precision,intent(in)   :: S(nBas,nBas)
      double precision,intent(in)   :: T(nBas,nBas)
      double precision,intent(in)   :: V(nBas,nBas)
      double precision,intent(in)   :: Hc(nBas,nBas) 
      double precision,intent(in)   :: X(nBas,nBas) 
      double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

      double precision,intent(in)   :: ENuc

      logical         ,intent(in)   :: c_details
  
      ! Local variables
  
      integer,parameter             :: maxSCF = 100
      double precision,parameter    :: thresh = 1d-5
      integer                       :: nSCF
      double precision              :: Conv
      double precision              :: Gap
      double precision              :: Gap_alpha , Gap_beta
      double precision              :: ET,EV , EHF_old
      double precision              :: EK,EJ
      double precision              :: EK_alpha
      double precision              :: EK_beta

      double precision,allocatable  ::    cp_alpha(:,:) ,    cp_beta(:,:)
      double precision,allocatable  ::     P_alpha(:,:) ,     P_beta(:,:)
      double precision,allocatable  ::           J(:,:)
      double precision,allocatable  ::     K_alpha(:,:) ,     K_beta(:,:)
      double precision,allocatable  ::     F_alpha(:,:) ,     F_beta(:,:) 
      double precision,allocatable  ::    Fp_alpha(:,:) ,    Fp_beta(:,:)
      double precision,allocatable  :: error_alpha(:,:) , error_beta(:,:)


      double precision,allocatable  ::     P_total(:,:)

      double precision,external     :: trace_matrix

      integer                       :: max_diis = 0
      integer                       :: n_diis_alpha
      integer                       ::  n_diis_beta 
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
      double precision,intent(out)  :: c(nBas,nBas)

      !-----------------------------------------------------------------!

      !--------------------------  UHF      ----------------------------!

      double precision              :: e_alpha(nBas)
      double precision              :: e_beta(nBas)
      double precision              :: c_alpha(nBas,nBas)
      double precision              :: c_beta(nBas,nBas)
      character(len=10)             :: type_calc
      double precision              :: ET1, EV1, ET2, EV2

      !-----------------------------------------------------------------!

      
      if (c_details) then 
        open(HFfile,file=trim(tmp_file_name)//"/RHF.out")
      end if 

      write(outfile,*)
      write(outfile,*)'******************************************************************************************'
      write(outfile,*)'|                         Un-Restricted Hartree-Fock calculation                          |'
      write(outfile,*)'******************************************************************************************'
      write(outfile,*)
  
      ! Memory allocation
  
      allocate(cp_alpha(nBas,nBas), P_alpha(nBas,nBas) ,error_alpha(nBas,nBas))
      allocate(K_alpha(nBas,nBas))
      allocate(F_alpha(nBas,nBas) ,Fp_alpha(nBas,nBas))


      allocate(cp_beta(nBas,nBas), P_beta(nBas,nBas) ,error_beta(nBas,nBas))
      allocate(K_beta(nBas,nBas))
      allocate(F_beta(nBas,nBas) ,Fp_beta(nBas,nBas))


      allocate(P_total(nBas,nBas))
      allocate(J(nBas,nBas))


      allocate(err_diis(nBas*nBas,max_diis))
      allocate(F_diis(nBas*nBas,max_diis))
       
      type_calc = "alpha"
      call    guess_UHF(nBas,c_details,n_alpha,HC,X,ENuc,T,V,P_alpha,type_calc,ET1,EV1,ET2,EV2)
      type_calc = "beta"
      call    guess_UHF(nBas,c_details,n_beta ,HC,X,ENuc,T,V,P_beta,type_calc,ET1,EV1,ET2,EV2)
         
      write(outfile,'(a)') ""
      write(outfile,'(a,f16.10)')   "      The Kinetic   Energy            = ", ET1+ET2
      write(outfile,'(a,f16.10)')   "      The Potential Energy            = ", EV1+EV2
      write(outfile,'(a,f16.10,a)') "      The Nuclear   Energy            = ", ENuc , "  +"
      write(outfile,'(a,f16.10)')   "      -----------------------------------------"
      write(outfile,'(a,f16.10)')   "      The guess Energy                = ", ET1+ET2+EV1+EV2+ENuc
      write(outfile,'(a)') ""

      ! --------------------------------------------------------------- !
      !              check that P_{mu nu} S_{mu nu} = N 
      ! --------------------------------------------------------------- !

      sum = 0.d0 
      do mu = 1 , nbas 
        do nu = 1 , nbas 
          sum = sum + P_alpha(mu,nu)*S(mu,nu) + P_beta(mu,nu)*S(mu,nu)
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
            sum_o(o) = sum_o(o) + P_alpha(mu,nu)*S(mu,o)*S(nu,o) + P_beta(mu,nu)*S(mu,o)*S(nu,o)
          end do 
        end do 
      end do 

      call header_under("P_{mu nu} S_{mu s} S_{nu s}",-1)

      write(outfile,"(a)") "P_{mu nu} S_{mu s} S_{nu s} :  "
      do i = 1 , nbas 
        write(outfile,"(40x,i3,1x,f16.8)") i , sum_o(i)
      end do 

      ! --------------------------------------------------------------- !

      ! Initialization
  
      nSCF = 0
      Conv = 1d0

      n_diis_alpha  = 0
      n_diis_beta   = 0
      F_diis(:,:)   = 0d0
      err_diis(:,:) = 0d0
      rcond         = 0d0

      EHF_old       = 0d0 
      EHF           = 0d0 
  
      Gap_alpha     = 0.d0
      Gap_beta      = 0.d0

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
      
      do while(nSCF < maxSCF)
  
      !   Increment 
  
      nSCF = nSCF + 1
  
      !   Compute the Hartree potential J

      call hartree_potential(nBas,P_alpha+P_beta,ERI,J)

      ! ****************** !
  
      !   Compute the exchange potential K
      call exchange_potential_UHF(nBas,P_alpha,ERI,K_alpha)
      call exchange_potential_UHF(nBas,P_beta ,ERI,K_beta)

      ! ****************** !

      !   Build Fock operator
    
      F_alpha(:,:) = Hc(:,:) + J(:,:) - K_alpha(:,:)
      F_beta (:,:) = Hc(:,:) + J(:,:) -  K_beta(:,:)
       
      !   Compute the error vector and extract the convergence criterion

      error_alpha = matmul(F_alpha,matmul(P_alpha,S)) - matmul(matmul(S,P_alpha),F_alpha)
      error_beta  = matmul(F_beta , matmul(P_beta,S)) - matmul(matmul(S,P_beta) , F_beta)

      if(nSCF > 1) Conv = max(maxval(abs(error_alpha)), maxval(abs(error_beta)))
      
      if (Conv < thresh) exit 

      ! ****************** !
  
      !---------------------------------------------------------------- !
      !   Compute HF energy
      !---------------------------------------------------------------- !

      P_total(:,:) = P_alpha(:,:) + P_beta(:,:)
  
      !   Compute the kinetic energy
      ET = trace_matrix(nBas,matmul(P_total,T))
      ! ****************** !

      !   Compute the potential energy
      EV = trace_matrix(nBas,matmul(P_total,V))
      ! ****************** !

      !   Compute the Hartree energy

      EJ = 0.5d0*trace_matrix(nBas,matmul(P_total,J))

      ! ****************** !

      !   Compute the exchange energy
      EK_alpha = -0.5d0 * trace_matrix(nBas,matmul(P_alpha,K_alpha))
      EK_beta  = -0.5d0 * trace_matrix(nBas,matmul(P_beta ,K_beta))

      EK = EK_alpha + EK_beta

      ! ****************** !

      !   Total HF energy
  
      EHF = ET + EV + EJ + EK

      if ((abs(EHF - EHF_old) < thresh ) .and. nSCF > maxSCF/4 ) exit

      if (nSCF > 2) then 
        EHF_old = EHF 
      end if 


      ! ****************** !
      !   Compute HOMO-LUMO gap
      ! ****************** !

      if (n_alpha > 0 ) then 
        if(nBas > n_alpha) then
          Gap_alpha = e_alpha(n_alpha+1) - e_alpha(n_alpha)
        else
          Gap_alpha = 0d0
        endif
      end if 

      if (n_beta > 0 ) then 
        if(nBas > n_beta) then
          Gap_beta = e_beta(n_beta+1) - e_beta(n_beta)
        else
          Gap_beta = 0d0
        endif
      end if 

      if (n_alpha > 0 .and. n_beta > 0) then 
        Gap = min(Gap_alpha, Gap_beta)
      else if (n_alpha <= 0) then 
        Gap = Gap_beta
      elseif (n_beta <= 0) then
        Gap = Gap_alpha
      end if 

      ! ****************** !

      ! DIIS extrapolation ! 

      if(max_diis > 1) then
        n_diis_alpha = min(n_diis_alpha+1,max_diis)
        n_diis_beta  = min(n_diis_beta+1,max_diis)
        call DIIS_extrapolation(rcond,nBas*nBas,nBas*nBas,n_diis_alpha,err_diis,F_diis,error_alpha,F_alpha)
        call DIIS_extrapolation(rcond,nBas*nBas,nBas*nBas, n_diis_beta,err_diis,F_diis,error_beta , F_beta)
      end if

      !   Transform for the Fock matrix F in the orthogonal basis 
      
      Fp_alpha = matmul(transpose(X),matmul(F_alpha,X))
      Fp_beta  = matmul(transpose(X),matmul(F_beta,X))
            
      !   Diagonalize F' to get MO coefficients (eigenvectors in the orthogonal basis) c' and MO energies (eigenvalues) e

      cp_alpha(:,:) = Fp_alpha(:,:)
      cp_beta (:,:) =  Fp_beta(:,:)

      call diagonalize_matrix(nBas,cp_alpha,e_alpha)
      call diagonalize_matrix(nBas,cp_beta ,e_beta)

      !   Back-transform the MO coefficients c in the original non-orthogonal basis

      c_alpha = matmul(X,cp_alpha)
      c_beta  = matmul(X,cp_beta)

      !   Compute the density matrix P 

      P_alpha(:,:) = matmul(c_alpha(:,1:n_alpha),transpose(c_alpha(:,1:n_alpha)))
      P_beta(:,:)  = matmul(  c_beta(:,1:n_beta),transpose(c_beta (:,1:n_beta )))

      !   Dump results
   
      write(outfile,"(1x,a1,i2,1x,a1,f16.8,1x,a1,f10.6,1x,a1,f10.6,1x,a1,f16.8,4x,a1,f16.8,4x,a1,f16.8,4x,a1)")      & 
      "|",nSCF,"|",EHF+ENuc,"|",Conv,"|",Gap,"|",ET,"|",EV,"|",Ej+EK,"|"
      write(outfile,*) repeat('-', 110)
      
      enddo

      !------------------------------------------------------------------------
      ! End of SCF loop
      !------------------------------------------------------------------------
  
      if(nSCF == maxSCF) then
  
        write(outfile,*)
        write(outfile,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(outfile,*)'                 Convergence failed                 '
        write(outfile,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(outfile,*)

        call print_UHF(nBas,n_alpha,n_beta,e_alpha,e_beta,C_alpha,C_beta,ENuc,ET,EV,EJ,EK,EHF)
  
        stop
  
      endif
  
      ! Compute final HF energy

        call print_UHF(nBas,n_alpha,n_beta,e_alpha,e_beta,C_alpha,C_beta,ENuc,ET,EV,EJ,EK,EHF)
  
end subroutine UHF