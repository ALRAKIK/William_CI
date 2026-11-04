subroutine AO_to_MO_ERI_chemist(nBas,c,ERI_AO,ERI_MO)

      ! Expression of bi-electronic integrals in the MO basis set

      use files 
      use keywords

      implicit none

      ! Input variables
      
      integer,intent(in)            :: nBas
      double precision,intent(in)   :: c(nBas,nBas)
      double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
      
      ! Local variables
      
      integer                       :: mu,nu,la,si
      integer                       :: p,q,r,s
      double precision,allocatable  :: scr(:,:,:,:)

      ! Output variables

      double precision,intent(out)  :: ERI_MO(nBas,nBas,nBas,nBas)

      ! Memory allocation

      allocate(scr(nBas,nBas,nBas,nBas))

  
      !--------------------------------------
      ! AO to MO transformation starts here !
      !--------------------------------------

      ! Transform 4th index

      scr(:,:,:,:) = 0d0

      do mu=1,nBas
         do nu=1,nBas
           do la=1,nBas
             do si=1,nBas
            
               do s=1,nBas
              
                 scr(mu,nu,la,s) = scr(mu,nu,la,s) + c(si,s)*ERI_AO(mu,nu,la,si)
              
               enddo    
             
             enddo
          enddo
        enddo
      enddo

      ! Transform 3rd index

      ERI_MO(:,:,:,:) = 0d0

      do mu=1,nBas
        do nu=1,nBas
          do la=1,nBas
          
            do r=1,nBas
              do s=1,nBas
              
                ERI_MO(mu,nu,r,s) = ERI_MO(mu,nu,r,s) + c(la,r)*scr(mu,nu,la,s)
              
              enddo
            enddo    
          
          enddo
        enddo
      enddo

      ! Transform 2nd index

      scr(:,:,:,:) = 0d0

      do mu=1,nBas
        do nu=1,nBas
        
          do q=1,nBas
            do r=1,nBas
              do s=1,nBas
              
                scr(mu,q,r,s) = scr(mu,q,r,s) + c(nu,q)*ERI_MO(mu,nu,r,s)
              
              enddo
            enddo
          enddo    
        
        enddo
      enddo

      ! Transform 1st (and last) index

      ERI_MO(:,:,:,:) = 0d0

      do mu=1,nBas
      
        do p=1,nBas
          do q=1,nBas
            do r=1,nBas
              do s=1,nBas
              
                ERI_MO(p,q,r,s) = ERI_MO(p,q,r,s) + c(mu,p)*scr(mu,q,r,s)
              
              enddo
            enddo
          enddo
        enddo
      
      enddo

end subroutine AO_to_MO_ERI_chemist


subroutine AO_to_MO_ERI_physicist(nBas,c,ERI_AO,ERI_MO)

      ! AO -> MO transformation for physicist-notation ERIs
      ! Input  : ERI_AO(mu,nu,la,si) = (mu nu | la si) in chemist notation
      ! Output : ERI_MO(p,q,r,s)     = <p q | r s> = (p r | q s)

      use files
      use keywords

      implicit none

      ! Input variables
      integer,intent(in)            :: nBas
      double precision,intent(in)   :: c(nBas,nBas)
      double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

      ! Output variables
      double precision,intent(out)  :: ERI_MO(nBas,nBas,nBas,nBas)

      ! Local variables
      integer                       :: mu,nu,la,si
      integer                       :: p,q,r,s
      double precision,allocatable  :: scr(:,:,:,:)

      allocate(scr(nBas,nBas,nBas,nBas))

      !-------------------------------------------------
      ! Transform 4th AO index: si -> s
      ! scr(mu,nu,la,s) = sum_si C(si,s) ERI_AO(mu,nu,la,si)
      !-------------------------------------------------
      scr(:,:,:,:) = 0d0

      do mu=1,nBas
        do nu=1,nBas
          do la=1,nBas
            do s=1,nBas
              do si=1,nBas
                scr(mu,nu,la,s) = scr(mu,nu,la,s) &
                                + c(si,s)*ERI_AO(mu,nu,la,si)
              enddo
            enddo
          enddo
        enddo
      enddo

      !-------------------------------------------------
      ! Transform 3rd AO index: la -> q
      ! ERI_MO(mu,nu,q,s) = sum_la C(la,q) scr(mu,nu,la,s)
      !-------------------------------------------------
      ERI_MO(:,:,:,:) = 0d0

      do mu=1,nBas
        do nu=1,nBas
          do q=1,nBas
            do s=1,nBas
              do la=1,nBas
                ERI_MO(mu,nu,q,s) = ERI_MO(mu,nu,q,s) &
                                  + c(la,q)*scr(mu,nu,la,s)
              enddo
            enddo
          enddo
        enddo
      enddo

      !-------------------------------------------------
      ! Transform 2nd AO index: nu -> r
      ! scr(mu,r,q,s) = sum_nu C(nu,r) ERI_MO(mu,nu,q,s)
      !-------------------------------------------------
      scr(:,:,:,:) = 0d0

      do mu=1,nBas
        do r=1,nBas
          do q=1,nBas
            do s=1,nBas
              do nu=1,nBas
                scr(mu,r,q,s) = scr(mu,r,q,s) &
                              + c(nu,r)*ERI_MO(mu,nu,q,s)
              enddo
            enddo
          enddo
        enddo
      enddo

      !-------------------------------------------------
      ! Transform 1st AO index: mu -> p
      ! ERI_MO(p,q,r,s) = sum_mu C(mu,p) scr(mu,r,q,s)
      !-------------------------------------------------
      ERI_MO(:,:,:,:) = 0d0

      do p=1,nBas
        do q=1,nBas
          do r=1,nBas
            do s=1,nBas
              do mu=1,nBas
                ERI_MO(p,q,r,s) = ERI_MO(p,q,r,s) &
                                + c(mu,p)*scr(mu,r,q,s)
              enddo
            enddo
          enddo
        enddo
      enddo

      deallocate(scr)

end subroutine AO_to_MO_ERI_physicist