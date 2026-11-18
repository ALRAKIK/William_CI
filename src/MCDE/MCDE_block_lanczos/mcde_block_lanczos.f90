! ===================================================================== !
!  [Nov 2026] Block (multi-vector) Lanczos spectral function for MCDE   !
!                                                                       !
!  New, self-contained routine. It does not replace or call the old     !
!  Lanczos_Haydock / diagonalize_matrix_d / spectral_function.          !
!                                                                       !
!  It computes the same quantity as                                     !
!     diagonalize_matrix_d + spectral_function :                        !
!     A(w) = sum_{i=1..nStart} sum_j |U(i,j)|^2 eta/pi /((E_j-w)^2+eta^2)
!          = -1/pi Im sum_{i=1..nStart} [ (w + i eta - H)^-1 ]_{ii}     !
!  without diagonalizing H.                                             !
!                                                                       !
!  Method: one Haydock (Lanczos) chain per start vector e_i,            !
!  i = 1..nStart, all chains advanced together so that every step is   !
!  ONE matrix-block product  W = H * V  (dgemm, N x N times N x nStart),!
!  which uses all threads of the BLAS library. Each chain gives a       !
!  tridiagonal (a_k, b_k); G_ii(z) is the continued fraction            !
!     G_ii(z) = 1/(z - a_1 - b_2^2/(z - a_2 - b_3^2/(...)))            !
!  Memory: 3 blocks of N x nStart doubles on top of H.                  !
!                                                                       !
!  H must be real symmetric (as for dsyevd). H is not modified.         !
! ===================================================================== !

subroutine mcde_blanczos_spectrum(N, nStart, H, nIter, e_min, e_max, delta, eta, A_w)

      use omp_lib
      implicit none

      ! Input variables
      integer         , intent(in)    :: N, nStart, nIter
      double precision, intent(in)    :: H(N,N)
      double precision, intent(in)    :: e_min, e_max, delta, eta

      ! Output: A_w(:,1) = omega, A_w(:,2) = spectral function (same grid as spectral_function)
      double precision, intent(inout) :: A_w(int((e_max-e_min)/delta),2)

      ! Local variables
      double precision, allocatable   :: Vold(:,:), V(:,:), W(:,:)
      double precision, allocatable   :: a(:,:), b(:,:)
      integer         , allocatable   :: nlen(:)
      logical         , allocatable   :: alive(:)
      double precision, parameter     :: PI  = 4.0d0*atan(1.0d0)
      double precision, parameter     :: tol = 1.0d-10
      integer                         :: k, c, iw, nw, nKry, nAlive
      double precision                :: s, omega, Gsum, t0, t1
      complex(kind=8)                 :: z, g

      nKry = min(nIter, N)
      nw   = int((e_max-e_min)/delta)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '        |  Block Lanczos spectral function  |'
      write(*,*) '---------------------------------------------------'
      write(*,'(a,i12)') '  Matrix dimension      = ', N
      write(*,'(a,i12)') '  Start vectors (chains)= ', nStart
      write(*,'(a,i12)') '  Lanczos steps         = ', nKry
      write(*,'(a,f12.2,a)') '  Work memory (3 blocks)= ', 3.0d0*8.0d0*dble(N)*dble(nStart)/1.0d9, ' GB'
      write(*,*)

      allocate(Vold(N,nStart), V(N,nStart), W(N,nStart))
      allocate(a(nKry,nStart), b(nKry+1,nStart), nlen(nStart), alive(nStart))

      ! start vectors e_1 ... e_nStart
      Vold = 0.0d0
      V    = 0.0d0
      do c = 1, nStart
         V(c,c) = 1.0d0
      end do
      a     = 0.0d0
      b     = 0.0d0
      nlen  = nKry
      alive = .true.

      t0 = omp_get_wtime()

      do k = 1, nKry

         ! W = H * V   (the expensive, threaded step)
         call dgemm('N','N', N, nStart, N, 1.0d0, H, N, V, N, 0.0d0, W, N)

         !$omp parallel do default(shared) private(c, s) schedule(dynamic)
         do c = 1, nStart
            if (.not. alive(c)) then
               W(:,c) = 0.0d0
               cycle
            end if
            a(k,c) = dot_product(V(:,c), W(:,c))
            W(:,c) = W(:,c) - a(k,c)*V(:,c) - b(k,c)*Vold(:,c)
            s = norm2(W(:,c))
            b(k+1,c) = s
            if (s < tol .or. k == nKry) then
               nlen(c)  = k
               alive(c) = .false.
               b(k+1,c) = 0.0d0
               W(:,c)   = 0.0d0
            else
               W(:,c) = W(:,c) / s
            end if
            Vold(:,c) = V(:,c)
            V(:,c)    = W(:,c)
         end do
         !$omp end parallel do

         nAlive = count(alive)
         if (mod(k,10) == 0 .or. nAlive == 0 .or. k == 1) then
            t1 = omp_get_wtime()
            write(*,'(a,i6,a,i6,a,f10.1,a)') '  step ', k, '   chains alive ', nAlive, &
                 '   elapsed ', t1-t0, ' s'
            flush(6)
         end if
         if (nAlive == 0) exit

      end do

      deallocate(Vold, V, W)

      ! spectral function from the continued fractions
      !$omp parallel do default(shared) private(iw, omega, z, c, k, g, Gsum) schedule(static)
      do iw = 1, nw
         omega = (iw-1)*delta + e_min
         z     = cmplx(omega, eta, kind=8)
         Gsum  = 0.0d0
         do c = 1, nStart
            g = (0.0d0, 0.0d0)
            do k = nlen(c), 1, -1
               g = 1.0d0 / (z - a(k,c) - b(k+1,c)**2 * g)
            end do
            Gsum = Gsum - aimag(g)/PI
         end do
         A_w(iw,1) = omega
         A_w(iw,2) = Gsum
      end do
      !$omp end parallel do

      write(*,'(a,f10.1,a)') '  Block Lanczos finished in ', omp_get_wtime()-t0, ' s'
      write(*,*)

      deallocate(a, b, nlen, alive)

end subroutine mcde_blanczos_spectrum
