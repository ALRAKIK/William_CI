subroutine DIIS_extrapolation(rcond, n_err, n_e, n_diis, error, e, error_in, e_inout)

      implicit none
      integer, intent(in)             :: n_err, n_e
      double precision, intent(in)    :: error_in(n_err)
      double precision, intent(inout) :: error(n_err, n_diis)
      double precision, intent(inout) :: e(n_e, n_diis)
      double precision, intent(out)   :: rcond
      integer, intent(inout)          :: n_diis
      double precision, intent(inout) :: e_inout(n_e)

      double precision, allocatable   :: A(:,:), b(:), w(:)
      integer                         :: info = 0

      allocate(A(n_diis+1, n_diis+1), b(n_diis+1), w(n_diis+1))

      ! Update history
      call prepend(n_err, n_diis, error, error_in)
      call prepend(n_e,   n_diis, e,     e_inout)

      ! Build A matrix
      
      call dgemm('T', 'N', n_diis, n_diis, n_err, 1.0d0, error, n_err, error, n_err, 0.0d0, A, n_diis+1)

      A(1:n_diis, n_diis+1) = -1.0d0
      A(n_diis+1, 1:n_diis) = -1.0d0
      A(n_diis+1, n_diis+1) =  0.0d0

      b(1:n_diis) = 0.0d0
      b(n_diis+1) = -1.0d0

      call linear_solve(n_diis+1, A, b, w, rcond)

      ! If solver fails (info != 0) or rcond too small, reset DIIS and skip extrapolation

      if (info /= 0 .or. rcond < 1.0d-14) then
        ! Reset DIIS history
        n_diis = 0
        error  = 0.0d0
        e      = 0.0d0
        deallocate(A, b, w)
        return
      end if

      ! Extrapolate Fock matrix

      call dgemm('N', 'N', n_e, 1, n_diis, 1.0d0, e, n_e, w, n_diis, 0.0d0, e_inout, n_e)

      deallocate(A, b, w)

end subroutine DIIS_extrapolation

subroutine linear_solve(N,A,b,x,rcond)

      ! Solve the linear system A.x = b where A is a NxN matrix
      ! and x and x are vectors of size N
    
      implicit none
    
      integer,intent(in)             :: N
      double precision,intent(inout) :: A(N,N)
      double precision,intent(out)   :: b(N),rcond
      double precision,intent(out)   :: x(N)
    
      integer                        :: info,lwork
      double precision               :: ferr,berr
      integer,allocatable            :: ipiv(:),iwork(:)
      double precision,allocatable   :: AF(:,:),work(:)
    
      ! Find optimal size for temporary arrays
    
      allocate(work(1))
      allocate(AF(N,N),ipiv(N),iwork(N))
    
      lwork = -1
      call dsysvx('N','U',N,1,A,N,AF,N,ipiv,b,N,x,N,rcond,ferr,berr,work,lwork,iwork,info)
      lwork = int(work(1))
    
      deallocate(work)
    
      allocate(work(lwork))
    
      call dsysvx('N','U',N,1,A,N,AF,N,ipiv,b,N,x,N,rcond,ferr,berr,work,lwork,iwork,info)
    
end subroutine

subroutine prepend(N,M,A,b)

      ! Prepend the vector b of size N into the matrix A of size NxM
  
      implicit none
  
      ! Input variables
  
      integer,intent(in)            :: N,M
      double precision,intent(in)   :: b(N)
  
      ! Local viaruabkes
  
      integer                       :: i,j
  
      ! Output variables
  
      double precision,intent(inout)  :: A(N,M)
  
  
      ! print*,'b in append'
      ! call matout(N,1,b)
  
      do i=1,N
        do j=M-1,1,-1
          A(i,j+1) = A(i,j)
        end do
        A(i,1) = b(i)
      end do
  
end subroutine


!==============================================================================
!  Complex DIIS extrapolation for the k-CCSD/k-LCCD amplitude iterations.
!
!  Mirrors the structure of the (real) HF DIIS you already use, but:
!    * works on complex(dpc) vectors,
!    * extrapolates the AMPLITUDE vector t2 from the residual/error history.
!
!  Conventions (matching solve_klccd):
!    * error vector  = the residual  R2 = -val_res / Dijab   (complex)
!    * vector stored = the amplitude t2                      (complex)
!
!  Usage inside solve_klccd, replacing the line
!        t2 = t2 + damping * R2
!  with
!        call diis_extrapolate_complex( namp, t2_flat, R2_flat, nhist, ndiis, &
!                                       t2_hist, res_hist, ok )
!
!  where t2_flat and R2_flat are flat complex(dpc) arrays of length namp.
!  On return, t2_flat holds the DIIS-extrapolated amplitudes and ok=.true.
!  If the solve failed or history was reset, ok=.false. and you should fall
!  back to the plain damped update t2 = t2 + damping * R2 for this step.
!
!  NOTE on the B-matrix sign convention: this uses the same convention as the
!  HF routine (B entries = <e_i|e_j>, off-diagonal = -1, rhs = [0,...,-1]),
!  which yields weights w that sum to 1.
!==============================================================================
subroutine diis_extrapolate_complex( namp, tvec, evec, nhist_inout, ndiis, &
                                     t_hist, e_hist, ok )
      use constants_module, only: dp, dpc
      implicit none
      ! Input
      integer,       intent(in)    :: namp              ! length of amplitude/residual
      integer,       intent(in)    :: ndiis             ! max history length
      complex(dpc),  intent(in)    :: evec(namp)        ! current residual (this iteration's error)
      complex(dpc),  intent(inout) :: tvec(namp)        ! in: updated amplitude; out: extrapolated
      complex(dpc),  intent(inout) :: t_hist(namp, ndiis)  ! history of amplitudes
      complex(dpc),  intent(inout) :: e_hist(namp, ndiis)  ! history of residuals
      integer,       intent(inout) :: nhist_inout       ! current number stored (0..ndiis)
      logical,       intent(out)   :: ok
      ! Local
      integer                      :: i, j, n
      integer                      :: nhist
      complex(dpc), allocatable    :: Bmat(:,:), rhs(:)
      integer,      allocatable    :: ipiv(:)
      integer                      :: info
      complex(dpc)                 :: zsum
      real(dp)                     :: bnorm_thresh

      ok = .false.
      bnorm_thresh = 1.0e-14_dp

      !---- Advance the rolling history: shift old columns right, store current at col 1
      nhist = min(nhist_inout + 1, ndiis)

      do j = ndiis, 2, -1
        e_hist(:, j) = e_hist(:, j-1)
        t_hist(:, j) = t_hist(:, j-1)
      end do
      e_hist(:, 1) = evec(:)
      t_hist(:, 1) = tvec(:)

      nhist_inout = nhist

      !---- After the FIRST iteration we only have one vector: no mix possible ----
      if (nhist < 2) then
        ok = .true.          ! nothing to extrapolate; keep tvec as-is
        return
      end if

      !---- Build augmented DIIS linear system of size (nhist+1) ----
      allocate(Bmat(nhist+1, nhist+1), rhs(nhist+1), ipiv(nhist+1))

      Bmat = cmplx(0.0_dp, 0.0_dp, kind=dpc)
      rhs  = cmplx(0.0_dp, 0.0_dp, kind=dpc)

      do i = 1, nhist
        do j = 1, nhist
          zsum = cmplx(0.0_dp, 0.0_dp, kind=dpc)
          do n = 1, namp
            zsum = zsum + conjg(e_hist(n, i)) * e_hist(n, j)
          end do
          Bmat(i, j) = zsum
        end do
        Bmat(i,       nhist+1) = cmplx(-1.0_dp, 0.0_dp, kind=dpc)
        Bmat(nhist+1, i)       = cmplx(-1.0_dp, 0.0_dp, kind=dpc)
      end do
      Bmat(nhist+1, nhist+1) = cmplx(0.0_dp, 0.0_dp, kind=dpc)

      rhs(nhist+1) = cmplx(-1.0_dp, 0.0_dp, kind=dpc)

      !---- Solve B*w = rhs (complex, general) ----
      call zgesv(nhist+1, 1, Bmat, nhist+1, ipiv, rhs, nhist+1, info)

      if (info /= 0) then
        ! Singular or failed: reset history, skip extrapolation
        nhist_inout = 0
        e_hist = cmplx(0.0_dp, 0.0_dp, kind=dpc)
        t_hist = cmplx(0.0_dp, 0.0_dp, kind=dpc)
        deallocate(Bmat, rhs, ipiv)
        return
      end if

      !---- Extrapolate: t_out = sum_j w_j * t_hist(:, j)  (weights = rhs(1:nhist))
      tvec(:) = cmplx(0.0_dp, 0.0_dp, kind=dpc)
      do j = 1, nhist
        do n = 1, namp
          tvec(n) = tvec(n) + rhs(j) * t_hist(n, j)
        end do
      end do

      ok = .true.
      deallocate(Bmat, rhs, ipiv)
end subroutine diis_extrapolate_complex