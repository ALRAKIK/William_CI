subroutine hartree_potential(nBas,P,ERI,J)

      ! Compute Coulomb matrix

      use keywords
      implicit none

      ! Input variables

      integer,intent(in)            :: nBas
      double precision,intent(in)   :: P(nBas,nBas)
      double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)


      ! Local variables

      integer                       :: mu,nu,la,si


      ! Output variables

      double precision,intent(out)  :: J(nBas,nBas)

      J(:,:) = 0d0


      if (c_physicist) then 

      do si=1,nBas
        do la=1,nBas
          do nu=1,nBas
            do mu=1,nBas
              J(mu,nu) = J(mu,nu) + P(la,si)*ERI(mu,la,nu,si)  ! physicist notation  
            end do
          end do
        end do
      end do

      else 

      do si=1,nBas
        do la=1,nBas
          do nu=1,nBas
            do mu=1,nBas
              J(mu,nu) = J(mu,nu) + P(la,si)*ERI(mu,nu,la,si)   !   chemist notation  
            end do
          end do
        end do
      end do


      end if 

      



end subroutine hartree_potential