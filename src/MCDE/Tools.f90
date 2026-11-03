integer function d(a,b)

      implicit none

      integer                                   :: a,b
   
      if (a == b) then
         d = 1.0d0
      else
         d = 0.0d0
      end if
    
      return

end function


subroutine inv_S(N,eigve, S)

      implicit none

      ! Input variables
      integer         , intent(in)            :: N
      double precision, intent(in)           :: eigve(N,N)
      double precision, intent(inout)        :: S(N,N)

      ! Local variable
      double precision, allocatable          :: Piv(:,:), LDA(:,:)
      integer                                :: i, j, length, err


      length = size(eigve(1,:))
      allocate(Piv(length,length))
      allocate(LDA(length,length))

      do i=1, length
         do j=1, length
            S(i,j)=dot_product(eigve(i,:),eigve(j,:))
         enddo
      enddo

      call dgetrf(length, length, S, length, Piv, err)

      call dgetri(length, S, length, Piv, LDA, length, err)

end subroutine inv_S


! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine diagonalize_matrix_nonh(N,A,e,c,VR,VL)

      ! =============================================================== !

      ! Diagonalize a square matrix

      use files
      use sorts
      implicit none

      ! =============================================================== !

      ! Input variables

      integer,intent(in)            :: N
      double precision,intent(inout):: A(N,N), VR(N,N), VL(N,N)
      double precision,intent(out)  :: e(N)
      double precision,intent(out)  :: c(N)

      ! =============================================================== !

      ! Local variables

      integer                       :: lwork,info
      integer                       :: i
      double precision,allocatable  :: work(:)

      ! =============================================================== !

      ! Memory allocation

      allocate(work(1))
      lwork = -1 
      call dgeev('V','V',N,A,N,e,c,VL,N,VR,N,work,lwork,info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))

      ! diagonalization 

      call dgeev('V','V',N,A,N,e,c,VL,N,VR,N,work,lwork,info)
      ! call bubble_sort(e)

      if(info /= 0) then
        write(outfile,'(a,I0)') 'Problem in diagonalize_matrix (dgeev)!!  --> ' , info 
        stop
      endif

      do i = 1 , N
        if (abs(e(i)) < 1d-15) e(i) = 0.d0
      end do

      deallocate(work)

end subroutine diagonalize_matrix_nonh

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !