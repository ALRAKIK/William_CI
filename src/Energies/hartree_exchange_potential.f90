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
      use iso_c_binding, only : c_loc, c_intptr_t

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

      ! [Nov 2026] skip the columns ERI(:,b,c,d) that are exactly zero
      ! (screened integrals; ~96% of them for a 50-cell chain). The list of
      ! non-zero columns is built once per ERI array and reused in every SCF
      ! iteration. Same sums in the same order as the dense loop, so the
      ! result is bit-for-bit identical. Kept dense when > 50% are non-zero.
      integer,    save, allocatable :: dptr(:), colb(:), colc(:)
      integer(8), save              :: key_loc = 0
      integer,    save              :: key_n1 = -1, key_nb = -1
      double precision, save        :: key_fp = 0.d0
      logical,    save              :: use_sparse = .false.
      integer(8)                    :: nnz, ii, loc_now
      integer,    allocatable       :: cnt(:)
      double precision              :: fp
      integer                       :: ip

      !-----------------------------------------------------------------!

      J(:,:) = 0.d0
      K(:,:) = 0.d0

      ! ---- (re)build the column list when the ERI array is a new one ---- !
      loc_now = int(loc(ERI), 8)
      fp = 0.d0                                     ! fingerprint: a sample of values
      do ii = 1 , 997
        d = int(mod(ii*7919_8, int(nBas,8))) + 1 ; c = int(mod(ii*104729_8, int(nBas,8))) + 1
        b = int(mod(ii*1299709_8, int(nBas,8))) + 1
        fp = fp + ERI(1+mod(int(ii),n1),b,c,d)*dble(ii)
      end do
      if (loc_now /= key_loc .or. n1 /= key_n1 .or. nBas /= key_nb .or. fp /= key_fp) then
        if (allocated(dptr)) deallocate(dptr, colb, colc)
        allocate(cnt(nBas))
        !$omp parallel do default(shared) private(b,c,d) schedule(dynamic)
        do d = 1 , nBas
          cnt(d) = 0
          do c = 1 , nBas
            do b = 1 , nBas
              if (any(ERI(:,b,c,d) /= 0.d0)) cnt(d) = cnt(d) + 1
            end do
          end do
        end do
        !$omp end parallel do
        nnz = sum(int(cnt,8))
        use_sparse = (nnz <= int(nBas,8)*int(nBas,8)*int(nBas,8)/2)
        if (use_sparse) then
          allocate(dptr(nBas+1), colb(max(nnz,1_8)), colc(max(nnz,1_8)))
          dptr(1) = 1
          do d = 1 , nBas
            dptr(d+1) = dptr(d) + cnt(d)
          end do
          !$omp parallel do default(shared) private(b,c,d,ip) schedule(dynamic)
          do d = 1 , nBas
            ip = dptr(d)
            do c = 1 , nBas
              do b = 1 , nBas
                if (any(ERI(:,b,c,d) /= 0.d0)) then
                  colb(ip) = b ; colc(ip) = c ; ip = ip + 1
                end if
              end do
            end do
          end do
          !$omp end parallel do
        end if
        deallocate(cnt)
        key_loc = loc_now ; key_n1 = n1 ; key_nb = nBas ; key_fp = fp
      end if

      if (use_sparse) then
        !$omp parallel default(shared) private(b,c,d,ip,Jloc,kcol)
        allocate(Jloc(n1,nBas), kcol(n1))
        Jloc = 0.d0
        if (.not. c_physicist) then
          !$omp do schedule(dynamic)
          do d = 1 , nBas
            kcol = 0.d0
            do ip = dptr(d) , dptr(d+1)-1
              b = colb(ip) ; c = colc(ip)
              Jloc(:,b) = Jloc(:,b) + ERI(:,b,c,d)*P(c,d)
              kcol(:)   = kcol(:)   + ERI(:,b,c,d)*P(b,c)
            end do
            K(1:n1,d) = -0.5d0*kcol(:)
          end do
          !$omp end do
        else
          !$omp do schedule(dynamic)
          do d = 1 , nBas
            kcol = 0.d0
            do ip = dptr(d) , dptr(d+1)-1
              b = colb(ip) ; c = colc(ip)
              Jloc(:,c) = Jloc(:,c) + ERI(:,b,c,d)*P(b,d)
              kcol(:)   = kcol(:)   + ERI(:,b,c,d)*P(c,b)
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
        return
      end if

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
