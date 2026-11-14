! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

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

      do i = 1 , N
          if (Eval(i) > thresh) then
            Eval(i) = dsqrt(1.d0 / Eval(i))
          else
            write(outfile,"(a,I3,a,4x,E16.10)") 'Eigenvalue ',i,' is too small for Lowdin orthogonalization',Eval(i)
            !Eval(i) = 0.d0
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

! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !
! --------------------------------------------------------------------- !

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
          do j=1,num
            if(abs(B(j)) < small) B(j) = 0d0
          enddo
          write(outfile,'(6x,I4,10F15.8)') i,(B(j),j=1,num)
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
