subroutine guess_AO_RHF_SAO(nBas,c_details,nO,HC,X,ENuc,S,T,V,P)

      use files
      use unitcell_module
      use constants_module


      implicit none 

      ! --------------------------------------------------------------- !
      integer         , intent(in)  ::                nBas , nO 
      double precision, intent(in)  ::             HC(nBas,nBas)
      double precision, intent(in)  ::              X(nBas,nBas)
      double precision, intent(in)  ::              T(nBas,nBas)
      double precision, intent(in)  ::              V(nBas,nBas)
      double precision, intent(in)  ::              S(nBas,nBas)
      double precision, intent(in)  ::                     ENuc
      logical         , intent(in)  ::                c_details
      double precision,external     ::             trace_matrix
      ! --------------------------------------------------------------- !
      double precision,allocatable  :: cp(:,:), c(:,:) , e(:) , ct(:,:)
      integer                       :: i      , o
      double precision              :: ET     , EV
      ! --------------------------------------------------------------- !
      double precision, intent(out) :: P(nbas,nbas)
      ! --------------------------------------------------------------- !

      allocate(cp(nbas,nbas),c(nbas,nbas),ct(nbas,nbas), e(nbas))

      cp(:,:) = 0.d0 

      cp(:,:) = matmul(transpose(X(:,:)), matmul(HC(:,:), X(:,:)))

      call diagonalize_matrix(nbas, cp, e)

      c(:,:)  = matmul(X(:,:), cp(:,:))

      ct(:,:) = transpose(c)

      P(:,:) = 2d0 * matmul(c(:,1:nO), transpose(c(:,1:nO)))

      ! --------------------------------------------------------------- !
      !                          Print                                  !
      ! --------------------------------------------------------------- !

      ET = trace_matrix(nBas,matmul(P,T))
      EV = trace_matrix(nBas,matmul(P,V))

      if (c_details) then 

      write(HFfile,'(3a)') "!" ,repeat('-',44), "!"
      write(HFfile,'(a)') "                 Iter  =   0"
      write(HFfile,'(3a)') "!" ,repeat('-',44), "!"
      write(HFfile,'(a)') ""
      
      call header_HF("Guess F Matrix, F = HC ", -1)

      write(HFfile,'(17x,1000(i3,15x))') (i,i=1,nfuc)
      do i = 1 , size(HC,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (HC(o,i),o=1,nfuc)
      end do 
      write(HFfile,'(a)') ""

      call header_HF("Guess X Matrix, X = S^(-1/2)", -1)
      write(HFfile,'(17x,1000(i3,15x))') (i,i=1,nfuc)
      do i = 1 , size(X,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (X(o,i),o=1,nfuc)
      end do 
      write(HFfile,'(a)') ""

      call header_HF("Guess Density Matrix, P = 2 C C^t", -1)
      write(HFfile,'(17x,1000(i3,15x))') (i,i=1,nfuc)
      do i = 1 , size(P,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (P(o,i),o=1,nfuc)
      end do 
      write(HFfile,'(a)') ""
      
      call header_HF("Guess Orbital Energies", -1)

      do i = 1 , size(e)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  e(i)
      end do 
      
      write(HFfile,'(a)') ""

      write(HFfile,'(a)') ""
      write(HFfile,'(a,f24.16)')   "      The Kinetic   Energy    = ", ET
      write(HFfile,'(a,f24.16)')   "      The Potential Energy    = ", EV
      write(HFfile,'(a,f24.16,a)') "      The Nuclear   Energy    = ", ENuc
      write(HFfile,'(a,f24.16)')   "-----------------------------------"
      write(HFfile,'(a,f24.16)') " The Energy of the guess  = ", ET+EV+ENuc
      write(HFfile,'(a)') ""

      write(HFfile,'(2a)')  repeat('*_',36) , "*"
      write(HFfile,'(a)')   repeat('_',73)

      end if 

      write(outfile,'(a)') ""
      write(outfile,'(a,f24.16)')   "      The Kinetic   Energy    = ", ET
      write(outfile,'(a,f24.16)')   "      The Potential Energy    = ", EV
      write(outfile,'(a,f24.16,a)') "      The Nuclear   Energy    = ", ENuc
      write(outfile,'(a,f24.16)')   "      -----------------------------------"
      write(outfile,'(a,f24.16)')   "      The guess Energy        = ", ET+EV+ENuc
      write(outfile,'(a)') ""


      ! --------------------------------------------------------------- !
      deallocate(cp,c,ct,e)
      ! --------------------------------------------------------------- !


end subroutine guess_AO_RHF_SAO