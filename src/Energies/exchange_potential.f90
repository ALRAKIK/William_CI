subroutine exchange_potential(nBas,P,ERI,K)

      ! Compute exchange matrix in the AO basis
      use keywords

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
      
      if (c_physicist .or. c_MCDE) then 

      do nu=1,nBas
        do si=1,nBas
          do la=1,nBas
            do mu=1,nBas
              K(mu,nu) = K(mu,nu) - 0.5d0 * P(la,si)*ERI(mu,si,la,nu)  ! physicist notation
            end do
          end do
        end do
      end do

      else 

      do nu=1,nBas
        do si=1,nBas
          do la=1,nBas
            do mu=1,nBas
              K(mu,nu) = K(mu,nu) - 0.5d0 * P(la,si)*ERI(mu,la,si,nu)   !   chemist notation
            end do
          end do
        end do
      end do

      end if 

      







end subroutine exchange_potential