! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine read_integrals(nBas,S,T,V,Hc,G)

      ! Read one- and two-electron integrals from files

      use files
      implicit none

      ! Input variables

      integer,intent(in)            :: nBas

      ! Local variables

      integer                       :: mu,nu,la,si
      double precision              :: Ov,Kin,Nuc,ERI

      double precision,intent(out)  :: S(nBas,nBas)
      double precision,intent(out)  :: T(nBas,nBas)
      double precision,intent(out)  :: V(nBas,nBas)
      double precision,intent(out)  :: Hc(nBas,nBas)
      double precision,intent(out)  :: G(nBas,nBas,nBas,nBas)

      open( 8,file=trim(tmp_file_name)// "/OV.dat" )
      open( 9,file=trim(tmp_file_name)// "/KI.dat" )
      open(10,file=trim(tmp_file_name)// "/NA.dat" )
      open(11,file=trim(tmp_file_name)// "/ERI.dat")


      S(:,:) = 0d0
      do 
        read(8,*,end=8) mu,nu,Ov
        S(mu,nu) = Ov
        S(nu,mu) = Ov
      enddo
      8 close(8)

      ! Read kinetic integrals

      T(:,:) = 0d0
      do 
        read(9,*,end=9) mu,nu,Kin
        T(mu,nu) = Kin
        T(nu,mu) = Kin
      enddo
      9 close(9)

      ! Read nuclear integrals

      V(:,:) = 0d0
      do 
        read(10,*,end=10) mu,nu,Nuc
        V(mu,nu) = Nuc
        V(nu,mu) = Nuc
      enddo
      10 close(unit=10)

      ! Define core Hamiltonian

      Hc(:,:) = T(:,:) + V(:,:)

      ! Read ERI integrals

      G(:,:,:,:) = 0.d0

      do 
      read(11,*,end=11) mu,nu,la,si,ERI
        G(mu, nu, la, si) = ERI
      enddo
      11 close(unit=11)


end subroutine read_integrals

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine read_integrals_from_file(nBas,S,T,V,Hc,G)

      ! Read one- and two-electron integrals from files

      implicit none

      ! Input variables

      integer,intent(in)            :: nBas

      ! Local variables

      integer                       :: mu,nu,la,si
      double precision              :: Ov,Kin,Nuc,ERI

      double precision,intent(out)  :: S(nBas,nBas)
      double precision,intent(out)  :: T(nBas,nBas)
      double precision,intent(out)  :: V(nBas,nBas)
      double precision,intent(out)  :: Hc(nBas,nBas)
      double precision,intent(out)  :: G(nBas,nBas,nBas,nBas)

      open( 8,file='./Read/OV.dat' )
      open( 9,file='./Read/KI.dat' )
      open(10,file='./Read/NA.dat' )
      open(11,file='./Read/ERI.dat')

      S(:,:) = 0d0
      do 
        read(8,*,end=8) mu,nu,Ov
        S(mu,nu) = Ov
        S(nu,mu) = Ov
      enddo
      8 close(8)

      ! Read kinetic integrals

      T(:,:) = 0d0
      do 
        read(9,*,end=9) mu,nu,Kin
        T(mu,nu) = Kin
        T(nu,mu) = Kin
      enddo
      9 close(9)

      ! Read nuclear integrals

      V(:,:) = 0d0
      do 
        read(10,*,end=10) mu,nu,Nuc
        V(mu,nu) = Nuc
        V(nu,mu) = Nuc
      enddo
      10 close(unit=10)

      ! Define core Hamiltonian

      Hc(:,:) = T(:,:) + V(:,:)

      ! Read ERI integrals

      G(:,:,:,:) = 0.d0

      do 
      read(11,*,end=11) mu,nu,la,si,ERI
        G(mu, nu, la, si) = ERI
      enddo
      11 close(unit=11)


end subroutine read_integrals_from_file

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine get_X_from_overlap(N,over,X)

      use files 
      implicit none 

      ! input !

      integer,intent(in)            :: N
      double precision, intent(in)  :: over(N,N)

      ! local !

      integer                       :: i , o
      double precision              :: Evec(N,N)
      double precision              :: Eval(N)
      double precision,parameter    :: thresh = 1d-6

      double precision              :: Xt(N,N)

      ! output !

      double precision, intent(out) :: X(N,N)
      

      Evec(:,:) = over(:,:)

      call diagonalize_matrix(N,Evec,Eval)

      ! //////////////////////////////////////////////////////////////// !
      !call header_under("The Overlap Eigenvalues",-1)
      !write(outfile, "(5f16.8)") (Eval(i), i=1,n)
      !write(outfile,*) ""
      !write(outfile,*) ""
      !write(outfile,'(a,g16.8)') " The smallest eigenvalue : " , MINVAL(Eval)
      !write(outfile,*) ""
      !if (Minval(Eval) < 0.d0) then 
      !  call header("Error",-1)
      !    write(outfile,'(a80)')'*******************************************************************************************'
      !    write(outfile,'(a80)')           "* The smallest eigenvalue of the overlap is negative, exiting the program      *"
      !    write(outfile,'(a80)')'*******************************************************************************************'
      !    stop            
      !end if
      ! //////////////////////////////////////////////////////////////// !

      do i = 1 , N
          if (Eval(i) > thresh) then
            Eval(i) = dsqrt(1.d0 / Eval(i))
          else
            write(outfile,"(a,I3,a,4x,E16.10)") 'Eigenvalue ',i,' is too small for Lowdin orthogonalization',Eval(i)
          end if
      end do
    
      call ADAt(N,Evec,Eval,X)

      open(1,file=trim(tmp_file_name)//"/X.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(X,1))
      do i = 1 , size(X,1)
        write(1,'(i3,6x,1000(f16.10,2x))') i ,  (X(i,o),o=1,size(X,1))
      end do 
      close(1)

      Xt = matmul(transpose(X),(matmul(over,X)))

      ! --------------------------------------------------------------- !

      !call header("check Unitary matrix X^(t) S X  = 1 ",-1)
      !call matout(N,N,Xt)

      ! --------------------------------------------------------------- !
  
end subroutine get_X_from_overlap

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine get_X_from_overlap_2D(N,over,X)

      use files 
      implicit none 

      ! input !

      integer,intent(in)            :: N
      double precision, intent(in)  :: over(N,N)

      ! local !

      integer                       :: i , o
      double precision              :: Evec(N,N)
      double precision              :: Eval(N)
      double precision,parameter    :: thresh = 1d-6

      double precision              :: Xt(N,N)
      double precision              :: a, b , c , d 
      double precision              :: trace_val , det_val , discriminant
      double precision              :: lambda1 , lambda2 
      double precision              :: eigenvector1(2) , eigenvector2(2)

      ! output !

      double precision, intent(out) :: X(N,N)

      Evec(:,:) = over(:,:)


      a = Evec(1,1)
      b = Evec(1,2)
      c = Evec(2,1)
      d = Evec(2,2)

      trace_val = a + d
      det_val   = a*d - b*c

      discriminant = trace_val**2 - 4.0d0 * det_val

      lambda1 = (trace_val + SQRT(discriminant)) / 2.0d0
      lambda2 = (trace_val - SQRT(discriminant)) / 2.0d0

      Eval(1) = lambda1 
      Eval(2) = lambda2 

      Evec(1,1) =  b  
      Evec(2,1) = -(a - lambda1)

      eigenvector1 = (/ b, -(a - lambda1) /)
      eigenvector1 = eigenvector1 / NORM2(eigenvector1)

      eigenvector2 = (/ b, -(a - lambda2) /)
      eigenvector2 = eigenvector2 / NORM2(eigenvector2)

      Evec(1,1) = eigenvector1(1)
      Evec(2,1) = eigenvector1(2)
      Evec(1,2) = eigenvector2(1)
      Evec(2,2) = eigenvector2(2)

      !call diagonalize_matrix(N,Evec,Eval)

      ! //////////////////////////////////////////////////////////////// !
      call header_under("The Overlap Eigenvalues",-1)

      write(outfile, "(5f16.8)") (Eval(i), i=1,n)

      write(outfile,*) ""
      write(outfile,*) ""
      write(outfile,'(a,g16.8)') " The smallest eigenvalue : " , MINVAL(Eval)
      write(outfile,*) ""

      if (Minval(Eval) < 0.d0) then 
        call header("Error",-1)
          write(outfile,'(a80)')'*******************************************************************************************'
          write(outfile,'(a80)')           "* The smallest eigenvalue of the overlap is negative, exiting the program      *"
          write(outfile,'(a80)')'*******************************************************************************************'
          stop            
      end if
      ! //////////////////////////////////////////////////////////////// !

      do i = 1 , N
          if (Eval(i) > thresh) then
            Eval(i) = dsqrt(1.d0 / Eval(i))
          else
            write(outfile,"(a,I3,a,4x,E16.10)") 'Eigenvalue ',i,' is too small for Lowdin orthogonalization',Eval(i)
            Eval(i) = dsqrt(1.d0 / Eval(i))
          end if
      end do
    
      call ADAt(N,Evec,Eval,X)

      open(1,file=trim(tmp_file_name)//"/X.dat ")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(X,1))
      do i = 1 , size(X,1)
        write(1,'(i3,6x,1000(f16.10,2x))') i ,  (X(i,o),o=1,size(X,1))
      end do 
      close(1)

      open(2,file=trim(tmp_file_name)//"/eigen_values.dat ")
      do i = 1 , N
      write(2,*) Eval(i)
      end do 
      close(2)

      Xt = matmul(transpose(X),(matmul(over,X)))

      ! --------------------------------------------------------------- !

!      call header("check Unitary matrix X^(t) S X  = 1 ",-1)
!      call matout(N,N,Xt)

      ! --------------------------------------------------------------- !
  
end subroutine get_X_from_overlap_2D

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine AtDA(N,A,D,B)

      ! Perform B = At.D.A where A is a NxN matrix and D is a diagonal matrix given
      ! as a vector of length N

      implicit none

      ! Input variables

      integer,intent(in)            :: N
      double precision,intent(in)   :: A(N,N),D(N)

      ! Local viaruabkes

      integer                       :: i,j,k

      ! Output variables
  
      double precision,intent(out)  :: B(N,N)
  
      B = 0d0
  
      do i=1,N
        do j=1,N
          do k=1,N
            B(i,k) = B(i,k) + A(j,i)*D(j)*A(j,k)
          enddo
        enddo
      enddo
  
end subroutine AtDA
  
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine ADAt(N,A,D,B)
  
      ! Perform B = A.D.At where A is a NxN matrix and D is a diagonal matrix given 
      ! as a vector of length N

      implicit none

      ! Input variables
  
      integer,intent(in)            :: N
      double precision,intent(in)   :: A(N,N),D(N)
    
      ! Local viaruabkes
  
      integer                       :: i,j,k
  
      ! Output variables
  
      double precision,intent(out)  :: B(N,N)
  
      B = 0d0

      do i=1,N
        do j=1,N
          do k=1,N
            B(i,k) = B(i,k) + A(i,j)*D(j)*A(k,j)
          enddo
        enddo
      enddo
  
end subroutine ADAt

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine matout(m,n,A)

      ! Print the MxN array A
  
      use files
      use sorts
      implicit none

      integer,parameter             :: ncol = 5
      double precision,parameter    :: small = 1d-10
      integer,intent(in)            :: m,n
      double precision,intent(in)   :: A(m,n)
      double precision              :: B(ncol)
      integer                       :: ilower,iupper,num,i,j


      do ilower=1,n,ncol
        iupper = min(ilower + ncol - 1,n)
        num = iupper - ilower + 1
        write(outfile,'(3X,10(9X,I6))') (j,j=ilower,iupper)
        do i=1,m
          do j=ilower,iupper
            B(j-ilower+1) = A(i,j)
          enddo
          ! call bubble_sort(B)
          do j=1,num
            if(abs(B(j)) < small) B(j) = 0d0
          enddo
          write(outfile,'(I7,10F15.8)') i,(B(j),j=1,num)
        enddo
      enddo
  
end subroutine matout

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

function trace_matrix(n,A) result(Tr)

      ! Calculate the trace of a square matrix A
  
      implicit none
  
      ! Input variables
  
      integer,intent(in)            :: n
      double precision,intent(in)   :: A(n,n)
  
      ! Local variables
  
      integer                       :: i
  
      ! Output variables

      double precision              :: Tr
  
      Tr = 0d0
      do i=1,n
        Tr = Tr + A(i,i)
      enddo
  
end function trace_matrix

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !


function trace_complex(N, A) result(Tr)

      ! Calculate the trace of a complex square matrix A

      use constants_module
      implicit none
      
      ! Input variables

      integer, intent(in)      :: n
      complex(dpc), intent(in) :: A(n,n)

      ! Local variables

      integer                  :: i

      ! Output variables

      complex(dpc)             :: Tr
      

      Tr = (0.0_dp, 0.0_dp)

      do i = 1, n
        Tr = Tr + A(i,i)
      end do

end function trace_complex

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine get_X_from_overlap_complex(N,over,X)

      use files 
      use constants_module
      implicit none 

      ! input !

      integer,intent(in)            :: N
      complex(dpc), intent(in)      :: over(N,N)

      ! local !

      integer                       :: i , o
      complex(dpc)                  :: Evec(N,N)
      double precision              :: Eval(N)
      double precision,parameter    :: thresh = 1d-6

      complex(dpc)                  :: Xt(N,N)

      ! output !

      complex(dpc), intent(out)     :: X(N,N)
      
      Evec(:,:) = over(:,:)

      call diagonalize_matrix_complex(N,Evec,Eval)

      do i = 1 , N
          if (Eval(i) > thresh) then
            Eval(i) = dsqrt(1.d0 / Eval(i))
          else
            write(outfile,"(a,I3,a,4x,E16.10)") 'Eigenvalue ',i,' is too small for Lowdin orthogonalization',Eval(i)
          end if
      end do
    
      call ADAt_complex(N,Evec,Eval,X)

      open(1,file=trim(tmp_file_name)//"/X.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(X,1))
      do i = 1 , size(X,1)
        write(1,'(i3,6x,1000(f16.10,2x))') i ,  (X(i,o),o=1,size(X,1))
      end do 
      close(1)

      Xt = matmul(conjg(transpose(X)),(matmul(over,X)))

      ! --------------------------------------------------------------- !

!      call header("check Unitary matrix X^(t) S X  = 1 ",-1)
!      call matout(N,N,Xt)

      ! --------------------------------------------------------------- !
  
end subroutine get_X_from_overlap_complex

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine ADAt_complex(N,A,D,B)
  
      ! Perform B = A.D.At where A is a NxN matrix and D is a diagonal matrix given 
      ! as a vector of length N

      use constants_module
      implicit none

      ! Input variables
  
      integer,intent(in)            :: N
      complex(dpc),intent(in)       :: A(N,N)
      double precision,intent(in)   :: D(N)
    
      ! Local viaruabkes
  
      integer                       :: i,j,k
  
      ! Output variables
  
      complex(dpc),intent(out)      :: B(N,N)
  
      B = 0d0

      do i=1,N
        do j=1,N
          do k=1,N
            B(i,k) = B(i,k) + A(i,j)*D(j)*conjg(A(k,j))
          enddo
        enddo
      enddo
  
end subroutine ADAt_complex

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

subroutine check_lowdin(N, S, X)
  use constants_module
  use files 
  implicit none
  integer, intent(in) :: N
  complex(dpc), intent(in) :: S(N,N), X(N,N)
  complex(dpc), allocatable :: Xt(:,:), diff(:,:)
  real(dp) :: max_diff
  integer :: i

  allocate(Xt(N,N), diff(N,N))

  ! Compute X† S X
  Xt = matmul(conjg(transpose(X)), matmul(S, X))

  diff = Xt
  do i = 1, N
    diff(i,i) = diff(i,i) - (1.0_dp, 0.0_dp)
  end do

  max_diff = maxval(abs(diff))
  write(outfile, '(a, ES15.7)') 'Max deviation from identity (X† S X):', max_diff

  if (max_diff < 1.0d-10) then
    write(outfile, '(a)') '✓ Lowdin orthogonalisation is CORRECT'
  else
    write(outfile, '(a)') '✗ Lowdin orthogonalisation is WRONG – check your conjugation'
  end if

  deallocate(Xt, diff)
  
end subroutine check_lowdin
