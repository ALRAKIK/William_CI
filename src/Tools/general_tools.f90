! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
subroutine diagonalize_matrix(N,A,e)

      ! =============================================================== !

      ! Diagonalize a square matrix

      use files
      implicit none

      ! =============================================================== !

      ! Input variables

      integer,intent(in)            :: N
      double precision,intent(inout):: A(N,N)
      double precision,intent(out)  :: e(N)

      ! =============================================================== !

      ! Local variables

      integer                       :: lwork,info
      integer                       :: i
      double precision,allocatable  :: work(:)

      ! =============================================================== !

      ! Memory allocation

      allocate(work(1))
      lwork = -1 
      call dsyev('V','U',N,A,N,e,work,lwork,info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))

      ! diagonalization 

      call dsyev('V','U',N,A,N,e,work,lwork,info)

      if(info /= 0) then
        write(outfile,'(a,I0)') 'Problem in diagonalize_matrix (dsyev)!!  --> ' , info 
        stop
      endif

      do i = 1 , N
        if (abs(e(i)) < 1d-15) e(i) = 0.d0
      end do

      deallocate(work)

end subroutine diagonalize_matrix

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine diagonalize_matrix_d(N,A,e)

      ! =============================================================== !

      ! Diagonalize a symmetric square matrix using DSYEVD
      ! (Divide & Conquer — scales much better with threads than DSYEV)

      use files
      implicit none

      ! =============================================================== !

      ! Input variables

      integer,intent(in)            :: N
      double precision,intent(inout):: A(N,N)
      double precision,intent(out)  :: e(N)

      ! =============================================================== !

      ! Local variables

      integer                       :: lwork, liwork, info
      integer                       :: i
      double precision,allocatable  :: work(:)
      integer,allocatable           :: iwork(:)

      ! =============================================================== !

      ! Workspace query

      allocate(work(1))
      allocate(iwork(1))
      lwork  = -1
      liwork = -1
      call dsyevd('V','U',N,A,N,e,work,lwork,iwork,liwork,info)

      if(info /= 0) then
        write(outfile,'(a,I0)') 'Problem in workspace query (dsyevd)!!  --> ',info
        stop
      endif

      lwork  = int(work(1))
      liwork = iwork(1)
      deallocate(work)
      deallocate(iwork)
      allocate(work(lwork))
      allocate(iwork(liwork))

      ! =============================================================== !

      ! Diagonalization

      call dsyevd('V','U',N,A,N,e,work,lwork,iwork,liwork,info)

      if(info /= 0) then
        write(outfile,'(a,I0)') 'Problem in diagonalize_matrix (dsyevd)!!  --> ',info
        stop
      endif

      ! =============================================================== !

      ! Zero out numerical noise (parallelized loop)

      !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
      do i = 1 , N
        if (abs(e(i)) < 1d-15) e(i) = 0.d0
      end do
      !$OMP END PARALLEL DO

      deallocate(work)
      deallocate(iwork)

end subroutine diagonalize_matrix_d 

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !



subroutine diagonalize_matrix_complex(N,A,e)

      ! =============================================================== !

      ! Diagonalize a complex square matrix

      use constants_module
      use files
      implicit none

      ! =============================================================== !

      ! Input variables

      integer,intent(in)            :: N
      complex(dpc),intent(inout)    :: A(N,N)
      double precision,intent(out)  :: e(N)

      ! =============================================================== !

      ! Local variables

      integer                       :: lwork,info
      integer                       :: i
      
      complex(dpc), allocatable     :: work(:)
      double precision, allocatable :: rwork(:)

      ! =============================================================== !

      ! Memory allocation

      allocate(rwork(max(1, 3*N-2)))
      allocate(work(1))

      lwork = -1
      call zheev('V', 'U', N, A, N, e, work, lwork, rwork, info)
      lwork = int(real(work(1))) + 1 
      deallocate(work)
      allocate(work(lwork))

      ! Diagonalization 

      call zheev('V', 'U', N, A, N, e, work, lwork, rwork, info)

      if(info /= 0) then
        write(outfile,'(a,I0)') 'Problem in diagonalize_matrix_complex (zheev)!!  --> ' , info 
        stop
      endif

      do i = 1 , N
        if (abs(e(i)) < 1d-15) e(i) = 0.d0
      end do

      deallocate(rwork, work)

end subroutine diagonalize_matrix_complex

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine check_the_overlap(N,over)

      use files 
      use keywords, only : s_eig_min      ! [Nov 2026] used by density fitting
      implicit none 

      ! input !

      integer,intent(in)            :: N
      double precision, intent(in)  :: over(N,N)

      ! local !

      integer                       :: i
      double precision              :: Evec(N,N)
      double precision              :: Eval(N)
      double precision,parameter    :: thresh = 1d-6       ! = lindep threshold of the X builders
      double precision,parameter    :: neg_tol = 1d-10     ! [Nov 2026] relative tolerance for "negative"
      double precision,parameter    :: machine_eps = 2.22E-16
      integer                       :: n_small
      double precision              :: cond
      double precision              :: double_digits, digits_lost, reliable_digits


      ! output !
      
      Evec(:,:) = over(:,:)

      call diagonalize_matrix(N,Evec,Eval)

      ! /////////////////////////////////////////////////////////////// !

      call header_under("The Overlap Eigenvalues",-1)

      write(outfile, "(5f16.8)") (Eval(i), i=1,n)

      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !

      s_eig_min       = MINVAL(Eval)
      double_digits   = log10(1.0 / machine_eps)
      ! [Nov 2026] a zero or round-off negative eigenvalue gave cond = ****
      ! and NaN digits: cond is now taken as infinite (singular) in that case
      if (MINVAL(Eval) > 0.d0) then
        cond            = MAXVAL(Eval)/MINVAL(Eval)
        digits_lost     = log10(cond)
      else
        cond            = huge(1.d0)
        digits_lost     = double_digits
      end if
      reliable_digits = double_digits - digits_lost
      n_small         = count(Eval <= thresh)


      write(outfile,*) ""
      write(outfile,*) ""
      write(outfile,'(5X,a,es16.6)') " The smallest eigenvalue   : " , MINVAL(Eval)
      write(outfile,'(5X,a,f16.8)')  " The Largest  eigenvalue   : " , MAXVAL(Eval)
      if (MINVAL(Eval) > 0.d0) then
        write(outfile,'(5X,a,es16.6)') " The conditional number    : " , cond
      else
        write(outfile,'(5X,a,a16)')    " The conditional number    : " , "infinite"
      end if
      write(outfile,'(5X,a,f16.4)') " Estimated digits lost     : " , digits_lost
      write(outfile,'(5X,a,f16.4)') " Estimated reliable digits : " , reliable_digits
      write(outfile,'(5X,a,i16)')   " Eigenvalues <= 1e-6       : " , n_small
      write(outfile,*) ""
      if (n_small > 0) then
        write(outfile,'(a)') " *** NOTE: near linear dependence. The eigenvectors of S below 1e-6 are removed   ***"
        write(outfile,'(a)') " ***       by canonical orthogonalization when X is built (see the HF output).    ***"
        write(outfile,*) ""
      end if
      if (cond > 1.0E12) then
        write(outfile,'(a)') " *** WARNING: Matrix is ill-conditioned! Results may be inaccurate. ***"
      end if
      if (cond > 1.0E15) then
        write(outfile,'(a)') " *** CRITICAL: Matrix is numerically singular! Results are likely garbage. ***"
      end if
      FLUSH(outfile)

      ! [Nov 2026] stop only for a TRULY negative eigenvalue. A tiny negative
      ! one (|e| ~ 1e-14, round-off of dsyev on a numerically singular S,
      ! e.g. aug basis on a dense chain) is a zero eigenvalue: it is removed
      ! by canonical orthogonalization in get_X_from_overlap*.
      if (Minval(Eval) < -neg_tol*max(1.d0,MAXVAL(Eval))) then 
        call header("Error",-1)
          write(outfile,'(a80)')'*******************************************************************************************'
          write(outfile,'(a80)')           "* The smallest eigenvalue of the overlap is negative, exiting the program      *"
          write(outfile,'(a80)')'*******************************************************************************************'
          FLUSH(outfile)
          stop            
      end if
      ! /////////////////////////////////////////////////////////////// !

      ! --------------------------------------------------------------- !
  
end subroutine check_the_overlap

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine diagonalize_block_diagonal_complex_para(nBas, nfuc, Nc, A, e)

      use constants_module
      implicit none
      integer, intent(in)         :: nBas, nfuc, Nc
      complex(dpc), intent(inout) :: A(nBas, nBas)
      real(dp), intent(out)       :: e(nBas)
      complex(dpc), allocatable   :: A_block(:,:)
      real(dp), allocatable       :: e_block(:)
      integer                     :: k, i, j


      !$omp parallel private(k, i, j, A_block, e_block)
      allocate(A_block(nfuc, nfuc), e_block(nfuc))
      !$omp do schedule(static)
        do k = 1, Nc
          i = (k-1)*nfuc + 1
          j = k*nfuc

          A_block(:,:) = A(i:i+nfuc-1, i:i+nfuc-1)

          call clean_and_symmetrize(nfuc, A_block)

          call diagonalize_matrix_complex(nfuc, A_block, e_block)

          e(i:i+nfuc-1) = e_block
          A(i:i+nfuc-1, i:i+nfuc-1) = A_block
        end do
      !$omp end do

      deallocate(A_block, e_block)
      !$omp end parallel
      

end subroutine diagonalize_block_diagonal_complex_para

subroutine print_band_structure(nBas, nfuc, Nc, e, filename)
  
      use constants_module
      implicit none

      integer, intent(in) :: nBas, nfuc, Nc
      real(dp), intent(in) :: e(nBas)
      character(len=*), intent(in) :: filename
      integer :: k, i, unit
      real(dp) :: k_value

      open(newunit=unit, file=filename, status='replace')
      
      do k = 0, Nc-1
        k_value = dble(k) / dble(Nc)  
        write(unit, '(F10.6, 100F16.10)') k_value, (e(k*nfuc + i), i = 1, nfuc)
      end do

  close(unit)

  write(*,*) 'Band structure written to: ', trim(filename)

end subroutine print_band_structure
