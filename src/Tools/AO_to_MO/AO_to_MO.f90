!=====================================================================!
!                                                                     !
!        AO -> MO transformation of the two-electron integrals        !
!                         IN PLACE (one N^4 array)                    !
!                                                                     !
!   chemist   : ERI(p,q,r,s) = sum C(mu,p) C(nu,q) C(la,r) C(si,s)    !
!                                   ERI(mu,nu,la,si)                  !
!   physicist : the same transformation of the stored array, then    !
!               indices 2 and 3 swapped:  ERI_MO(p,q,r,s) = X(p,r,q,s)!
!               (exactly what the former AO_to_MO_ERI_physicist did)  !
!                                                                     !
!   step A : indices 1,2 slice by slice   ERI(:,:,la,si) = C^T M C    !
!   step B : indices 3,4 on blocks of p   T(p,la,si) -> T(p,r,s)      !
!   step C : (physicist only) swap of indices 2 and 3, in place      !
!                                                                     !
!   Extra memory: per thread 2 N^2 (step A) and 2 pb N^2 (step B),    !
!   pb = block of p (<= 32). The O(N^5) work is done with matmul.     !
!=====================================================================!

subroutine AO_to_MO_ERI_inplace(nBas, c, ERI, physicist)

      implicit none

      ! - input - !

      integer          , intent(in)    :: nBas
      double precision , intent(in)    :: c(nBas,nBas)
      logical          , intent(in)    :: physicist

      ! - input / output - !

      double precision , intent(inout) :: ERI(nBas,nBas,nBas,nBas)

      ! - local - !

      integer          , parameter     :: pb_max = 32
      integer                          :: n, la, si, q, r, s, ib, nblk, p1, p2, m
      double precision , allocatable   :: ct(:,:), W(:,:), T(:,:,:), U(:,:,:), tmp(:)

      !-----------------------------------------------------------------!

      n  = nBas
      allocate(ct(n,n))
      ct = transpose(c)
      nblk = (n + pb_max - 1) / pb_max

      !$omp parallel default(shared) private(la,si,q,r,s,ib,p1,p2,m,W,T,U,tmp)

      ! ---- step A : (mu nu | la si) -> (p q | la si) ---- !
      allocate(W(n,n))
      !$omp do collapse(2) schedule(static)
      do si = 1 , n
        do la = 1 , n
          W = matmul(ct, ERI(:,:,la,si))
          ERI(:,:,la,si) = matmul(W, c)
        end do
      end do
      !$omp end do
      deallocate(W)

      ! ---- step B : (p q | la si) -> (p q | r s), blocks of p for fixed q ---- !
      !$omp do collapse(2) schedule(dynamic)
      do q = 1 , n
        do ib = 1 , nblk
          p1 = (ib-1)*pb_max + 1
          p2 = min(n, ib*pb_max)
          m  = p2 - p1 + 1
          allocate(T(m,n,n), U(m,n,n))
          T(:,:,:) = ERI(p1:p2,q,:,:)
          do si = 1 , n
            U(:,:,si) = matmul(T(:,:,si), c)            ! la -> r
          end do
          call mat_times(m*n, n, n, U, c, T)            ! si -> s
          ERI(p1:p2,q,:,:) = T(:,:,:)
          deallocate(T, U)
        end do
      end do
      !$omp end do

      ! ---- step C : physicist order, swap indices 2 and 3 ---- !
      if (physicist) then
        allocate(tmp(n))
        !$omp do schedule(dynamic)
        do s = 1 , n
          do r = 2 , n
            do q = 1 , r - 1
              tmp(:)        = ERI(:,q,r,s)
              ERI(:,q,r,s)  = ERI(:,r,q,s)
              ERI(:,r,q,s)  = tmp(:)
            end do
          end do
        end do
        !$omp end do
        deallocate(tmp)
      end if

      !$omp end parallel

      deallocate(ct)

contains

      ! C(m,n) = A(m,k) B(k,n) on explicit-shape views (no copies)
      subroutine mat_times(m, k, nn, A, B, Cm)
        integer          , intent(in)  :: m, k, nn
        double precision , intent(in)  :: A(m,k), B(k,nn)
        double precision , intent(out) :: Cm(m,nn)
        Cm = matmul(A, B)
      end subroutine mat_times

end subroutine AO_to_MO_ERI_inplace
