subroutine Hartree_Fock_SAO(nBas,nO,S,T,V,Hc,ERI,X,Enuc,EHF,e,c,coeff_SAO &
                       ,n_alpha, n_beta,blk,natom,Zat)

      use keywords
      use constants_module
      use unitcell_module
      use k_blocks

      implicit none

      ! input ! 

      integer,intent(in)            ::                         nBas
  
      integer,intent(in)            ::                          nO
      double precision,intent(in)   ::                 S(nBas,nBas)
      double precision,intent(in)   ::                 T(nBas,nBas)
      double precision,intent(in)   ::                 V(nBas,nBas)
      double precision,intent(in)   ::                Hc(nBas,nBas)
      double precision,intent(in)   ::                 X(nBas,nBas)
      double precision,intent(in)   :: ERI(ERI_size,nBas,nBas,nBas)

      integer        ,intent(in)    :: n_alpha , n_beta 

      double precision,intent(in)   :: ENuc

      integer         , intent(in)  :: natom
      integer         , intent(in)  :: Zat(natom)


      ! local ! 

      ! output ! 

      double precision,intent(out)  :: EHF 
      double precision,intent(out)  :: e(nBas)
      double precision,intent(out)  :: c(nBas,nBas)
      complex(dpc),    intent(out)  :: coeff_SAO(nBas,nBas)
      type(K_block),   intent(out)  :: blk(N_cell)
       

      if (c_UHF) then 
        call UHF(nBas,c_details,n_alpha,n_beta,S,T,V,Hc,ERI,X,Enuc,EHF,e,c)
        go to 100
      end if

      call RHF_SAO(nBas,nO,S,T,V,Hc,ERI,X,Enuc,EHF,e,coeff_SAO,blk,natom,Zat)

100 continue

end subroutine 


subroutine Hartree_Fock(nBas,nO,S,T,V,Hc,ERI,X,Enuc,EHF,e,c,n_alpha, n_beta,natom,Zat)

      use keywords
      use constants_module
      use unitcell_module

      implicit none

      ! input ! 

      integer,intent(in)            ::                         nBas
  
      integer,intent(in)            ::                          nO
      double precision,intent(in)   ::                 S(nBas,nBas)
      double precision,intent(in)   ::                 T(nBas,nBas)
      double precision,intent(in)   ::                 V(nBas,nBas)
      double precision,intent(in)   ::                Hc(nBas,nBas)
      double precision,intent(in)   ::                 X(nBas,nBas)
      double precision,intent(in)   :: ERI(ERI_size,nBas,nBas,nBas)

      integer        ,intent(in)    :: n_alpha , n_beta 

      double precision,intent(in)   :: ENuc

      integer         , intent(in)  :: natom
      integer         , intent(in)  :: Zat(natom)


      ! local ! 

      ! output ! 

      double precision,intent(out)  :: EHF 
      double precision,intent(out)  :: e(nBas)
      double precision,intent(out)  :: c(nBas,nBas)
       
      if (c_UHF) then 
        call UHF(nBas,c_details,n_alpha,n_beta,S,T,V,Hc,ERI,X,Enuc,EHF,e,c)
        go to 100
      end if

      call RHF(nBas,nO,S,T,V,Hc,ERI,X,Enuc,EHF,e,c,natom,Zat)

100 continue

end subroutine 