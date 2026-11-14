!=====================================================================!
!                                                                     !
!   [Nov 2026] Coulomb J and exchange K in ONE pass over the ERIs,    !
!   parallel with OpenMP. Replaces, in RHF and RHF_SAO, the separate  !
!   serial routines hartree_potential + exchange_potential and        !
!   hartree_potential_SAO + exchange_potential_SAO (same formulas).   !
!                                                                     !
!   n1 = nBas : normal path, ERI(nBas,nBas,nBas,nBas)                 !
!   n1 = nfuc : SAO path,    ERI(nfuc,nBas,nBas,nBas), only the rows  !
!               mu <= nfuc of J and K are computed (the others are 0, !
!               as before: the Fock matrix is then symmetrized        !
!               through the SAO transformation)                       !
!                                                                     !
!   chemist   ERI(a,b,c,d) = (ab|cd)                                  !
!      J(mu,nu) =        sum P(la,si) (mu nu|la si)                   !
!      K(mu,nu) = -1/2 * sum P(la,si) (mu la|si nu)                   !
!      one element E = ERI(a,b,c,d) gives                              !
!      J(a,b) += E P(c,d)      K(a,d) += -1/2 E P(b,c)                 !
!                                                                     !
!   physicist ERI(a,b,c,d) = (ac|bd)                                  !
!      J(a,c) += E P(b,d)      K(a,d) += -1/2 E P(c,b)                 !
!                                                                     !
!   Parallel over the last index d: column K(:,d) belongs to one      !
!   thread; J is accumulated per thread and summed at the end.        !
!   Why faster: the ERIs (the big array) are read ONCE instead of     !
!   twice, the inner loop runs over the contiguous first index, and   !
!   all cores work (the old routines were serial).                    !
!=====================================================================!

subroutine hartree_exchange_potential(n1, nBas, P, ERI, J, K)

      use keywords , only : c_physicist

      implicit none

      ! - input - !

      integer          , intent(in)  :: n1 , nBas
      double precision , intent(in)  :: P(nBas,nBas)
      double precision , intent(in)  :: ERI(n1,nBas,nBas,nBas)

      ! - output - !

      double precision , intent(out) :: J(nBas,nBas)
      double precision , intent(out) :: K(nBas,nBas)

      ! - local - !

      integer                        :: b , c , d
      double precision               :: pcd , pbc
      double precision , allocatable :: Jloc(:,:) , kcol(:)

      !-----------------------------------------------------------------!

      J(:,:) = 0.d0
      K(:,:) = 0.d0

      !$omp parallel default(shared) private(b,c,d,pcd,pbc,Jloc,kcol)
      allocate(Jloc(n1,nBas), kcol(n1))
      Jloc = 0.d0

      if (.not. c_physicist) then

        !$omp do schedule(dynamic)
        do d = 1 , nBas
          kcol = 0.d0
          do c = 1 , nBas
            pcd = P(c,d)
            do b = 1 , nBas
              pbc = P(b,c)
              Jloc(:,b) = Jloc(:,b) + ERI(:,b,c,d)*pcd
              kcol(:)   = kcol(:)   + ERI(:,b,c,d)*pbc
            end do
          end do
          K(1:n1,d) = -0.5d0*kcol(:)
        end do
        !$omp end do

      else

        !$omp do schedule(dynamic)
        do d = 1 , nBas
          kcol = 0.d0
          do c = 1 , nBas
            do b = 1 , nBas
              Jloc(:,c) = Jloc(:,c) + ERI(:,b,c,d)*P(b,d)
              kcol(:)   = kcol(:)   + ERI(:,b,c,d)*P(c,b)
            end do
          end do
          K(1:n1,d) = -0.5d0*kcol(:)
        end do
        !$omp end do

      end if

      !$omp critical
      J(1:n1,:) = J(1:n1,:) + Jloc(:,:)
      !$omp end critical

      deallocate(Jloc, kcol)
      !$omp end parallel

end subroutine hartree_exchange_potential
