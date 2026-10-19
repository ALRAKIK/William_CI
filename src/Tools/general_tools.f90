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
      implicit none 

      ! input !

      integer,intent(in)            :: N
      double precision, intent(in)  :: over(N,N)

      ! local !

      integer                       :: i
      double precision              :: Evec(N,N)
      double precision              :: Eval(N)
      double precision,parameter    :: thresh = 1d-6
      double precision,parameter    :: machine_eps = 2.22E-16
      double precision              :: cond
      double precision              :: double_digits, digits_lost, reliable_digits


      ! output !
      
      Evec(:,:) = over(:,:)

      call diagonalize_matrix(N,Evec,Eval)

      ! /////////////////////////////////////////////////////////////// !

      call header_under("The Overlap Eigenvalues",-1)

      write(outfile, "(5f16.8)") (Eval(i), i=1,n)

      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !

      cond            = MAXVAL(Eval)/MINVAL(Eval)
      double_digits   = log10(1.0 / machine_eps)
      digits_lost     = log10(cond)
      reliable_digits = double_digits - digits_lost


      write(outfile,*) ""
      write(outfile,*) ""
      write(outfile,'(5X,a,f16.8)') " The smallest eigenvalue   : " , MINVAL(Eval)
      write(outfile,'(5X,a,f16.8)') " The Largest  eigenvalue   : " , MAXVAL(Eval)
      write(outfile,'(5X,a,f16.8)') " The conditional number    : " , cond
      write(outfile,'(5X,a,f16.4)') " Estimated digits lost     : " , digits_lost
      write(outfile,'(5X,a,f16.4)') " Estimated reliable digits : " , reliable_digits
      write(outfile,*) ""
      if (cond > 1.0E12) then
        write(outfile,'(a)') " *** WARNING: Matrix is ill-conditioned! Results may be inaccurate. ***"
      end if
      if (cond > 1.0E15) then
        write(outfile,'(a)') " *** CRITICAL: Matrix is numerically singular! Results are likely garbage. ***"
      end if
      FLUSH(outfile)

      if (Minval(Eval) < 0.d0) then 
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

subroutine check_the_overlap_SAO(N,over)

      use files
      use unitcell_module
      use constants_module

      implicit none 

      ! input !

      integer,intent(in)            :: N
      double precision, intent(in)  :: over(N,N)

      ! local !

      integer                       :: i
      double precision              :: Evec(N,N)
      double precision              :: Eval(N)
      double precision,parameter    :: thresh = 1d-6
      double precision,parameter    :: machine_eps = 2.22E-16
      double precision              :: cond
      double precision              :: double_digits, digits_lost, reliable_digits

      complex(dpc)                  :: over_k(N,N)


      ! output !


      call transform_AO_to_SAO(N,nfuc,N_cell,over,over_k)
      
      call diagonalize_matrix_complex(N,over_k,Eval)

      ! /////////////////////////////////////////////////////////////// !

      call header_under("The Overlap Eigenvalues",-1)

      write(outfile, "(5f16.8)") (Eval(i), i=1,n)

      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !

      cond            = MAXVAL(Eval)/MINVAL(Eval)
      double_digits   = log10(1.0 / machine_eps)
      digits_lost     = log10(cond)
      reliable_digits = double_digits - digits_lost


      write(outfile,*) ""
      write(outfile,*) ""
      write(outfile,'(5X,a,f16.8)') " The smallest eigenvalue   : " , MINVAL(Eval)
      write(outfile,'(5X,a,f16.8)') " The Largest  eigenvalue   : " , MAXVAL(Eval)
      write(outfile,'(5X,a,f16.8)') " The conditional number    : " , cond
      write(outfile,'(5X,a,f16.4)') " Estimated digits lost     : " , digits_lost
      write(outfile,'(5X,a,f16.4)') " Estimated reliable digits : " , reliable_digits
      write(outfile,*) ""
      if (cond > 1.0E12) then
        write(outfile,'(a)') " *** WARNING: Matrix is ill-conditioned! Results may be inaccurate. ***"
      end if
      if (cond > 1.0E15) then
        write(outfile,'(a)') " *** CRITICAL: Matrix is numerically singular! Results are likely garbage. ***"
      end if
      FLUSH(outfile)

      if (Minval(Eval) < 0.d0) then 
        call header("Error",-1)
          write(outfile,'(a80)')'*******************************************************************************************'
          write(outfile,'(a80)')           "* The smallest eigenvalue of the overlap is negative, exiting the program      *"
          write(outfile,'(a80)')'*******************************************************************************************'
          FLUSH(outfile)
          stop            
      end if
      ! /////////////////////////////////////////////////////////////// !

      ! --------------------------------------------------------------- !
  
end subroutine check_the_overlap_SAO

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine diagonalize_block_diagonal_complex(nBas, nfuc, Nc, A, e)

      use constants_module
      implicit none
      integer, intent(in)         :: nBas, nfuc, Nc
      complex(dpc), intent(inout) :: A(nBas, nBas)                       ! block-diagonal; on output eigenvectors per block
      real(dp), intent(out)       :: e(nBas)                             ! eigenvalues in block order
      complex(dpc), allocatable   :: A_block(:,:)
      real(dp), allocatable       :: e_block(:)

      integer                     :: k, i, j

      do k = 1, Nc
        ! Indices for the k-th block
        i = (k-1)*nfuc + 1
        j = k*nfuc

        allocate(A_block(nfuc, nfuc), e_block(nfuc))

        A_block(:,:) = A(i:i+nfuc-1,i:i+nfuc-1)

        call clean_and_symmetrize(nfuc, A_block)

        call diagonalize_matrix_complex(nfuc, A_block, e_block)

        e(i:i+nfuc-1) = e_block
        A(i:i+nfuc-1, i:i+nfuc-1) = A_block

        deallocate(A_block, e_block)

      end do
      
end subroutine diagonalize_block_diagonal_complex

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

subroutine chemist_to_physics_notation(nBas, ERI)

      implicit none

      integer, intent(in)             :: nBas
      double precision, intent(inout) :: ERI(nBas, nBas, nBas, nBas)

      ! Local variables
      integer                         :: i, j, k, l
      double precision, allocatable   :: tmp(:,:,:,:)

      ! Allocate temporary array
      allocate(tmp(nBas, nBas, nBas, nBas))

      ! Copy original integrals
      tmp = ERI

      ! Convert:  <i j | k l> = (i k | j l)
      ! i.e. swap the second and third indices

      do i = 1, nBas
        do j = 1, nBas
          do k = 1, nBas
            do l = 1, nBas
              ERI(i, j, k, l) = tmp(i, k, j, l)
            end do
          end do
        end do
      end do

      deallocate(tmp)

end subroutine chemist_to_physics_notation