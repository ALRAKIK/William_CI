subroutine exchange_potential_UHF(nBas,P,ERI,K)
      use keywords      ! [Nov 2026] c_physicist

      ! Compute exchange matrix in the AO basis

      implicit none

      ! Input variables

      integer,intent(in)            :: nBas
      double precision,intent(in)   :: P(nBas,nBas)
      double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

      ! Local variables

      integer                       :: mu,nu,la,si

      ! Output variables

      double precision,intent(out)  :: K(nBas,nBas)
      
      K(:,:) = 0d0

      ! [Nov 2026] physicist-ordered integrals (keyword physicist or MCDE):
      ! before, this routine always assumed chemist order while the
      ! Coulomb part followed the keyword -> wrong UHF Fock matrix
      if (c_physicist) then
        do nu=1,nBas
          do si=1,nBas
            do la=1,nBas
              do mu=1,nBas
                K(mu,nu) = K(mu,nu) + P(la,si)*ERI(mu,si,la,nu)
              end do
            end do
          end do
        end do
        return
      end if

      do nu=1,nBas
        do si=1,nBas
          do la=1,nBas
            do mu=1,nBas
              K(mu,nu) = K(mu,nu) + P(la,si)*ERI(mu,la,si,nu)
            end do
          end do
        end do
      end do

end subroutine exchange_potential_UHF