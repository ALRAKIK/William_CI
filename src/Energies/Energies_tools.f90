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

      integer                       :: i , o , ndrop
      double precision              :: Evec(N,N)
      double precision              :: Eval(N)
      double precision,parameter    :: thresh = 1d-6

      double precision              :: Xt(N,N)

      ! output !

      double precision, intent(out) :: X(N,N)
      

      Evec(:,:) = over(:,:)

      call diagonalize_matrix(N,Evec,Eval)

      ! [Nov 2026] canonical orthogonalization when S is near singular:
      ! eigenvectors with eigenvalue <= thresh are removed (zero column of X),
      ! the others give X = U s^-1/2. Before, the tiny eigenvalue was kept
      ! (only a warning), so X had huge, noisy columns. No eigenvalue below
      ! thresh -> Lowdin X = S^-1/2 exactly as before.
      call canonical_X_real(N,Evec,Eval,thresh,X,ndrop)
      if (ndrop > 0) then
        write(outfile,'(a,i6,a,i6,a,es12.4,a)') ' Canonical orthogonalization (AO): ', ndrop, ' of ', N, &
             ' vectors removed (S eigenvalue <= ', thresh, ')'
      end if

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

      integer                       :: i , o , ndrop
      complex(dpc)                  :: Evec(N,N)
      double precision              :: Eval(N)
      double precision,parameter    :: thresh = 1d-6

      complex(dpc)                  :: Xt(N,N)

      ! output !

      complex(dpc), intent(out)     :: X(N,N)
      
      Evec(:,:) = over(:,:)

      call diagonalize_matrix_complex(N,Evec,Eval)

      ! [Nov 2026] canonical orthogonalization (see get_X_from_overlap)
      call canonical_X_complex(N,Evec,Eval,thresh,X,ndrop)
      if (ndrop > 0) then
        write(outfile,'(a,i6,a,i6,a,es12.4,a)') ' Canonical orthogonalization: ', ndrop, ' of ', N, &
             ' vectors removed (S eigenvalue <= ', thresh, ')'
      end if

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

!=====================================================================!
! [Nov 2026] Lowdin X = S^-1/2 of a k-block-diagonal complex overlap,  !
! block by block (same result as get_X_from_overlap_complex on the     !
! full matrix, since S^-1/2 of a block-diagonal matrix is the          !
! block-diagonal matrix of the S_k^-1/2). Cost N_cell * nfuc^3 instead !
! of nBas^3. Same threshold and message for small eigenvalues. X.dat  !
! is written only with the keyword Details.                           !
!=====================================================================!
subroutine get_X_from_overlap_complex_blocks(N,nf,nk,over,X)

      use files
      use keywords , only : c_details
      use constants_module

      implicit none

      integer         , intent(in)  :: N , nf , nk
      complex(dpc)    , intent(in)  :: over(N,N)
      complex(dpc)    , intent(out) :: X(N,N)

      integer                       :: k , b1 , b2 , i , o
      double precision, parameter   :: thresh = 1d-6
      complex(dpc)    , allocatable :: Evec(:,:) , Xb(:,:)
      double precision, allocatable :: Eval(:)
      integer                       :: ndrop(nk)
      double precision              :: emin(nk)

      X = (0.0_dp, 0.0_dp)

      ! [Nov 2026] canonical orthogonalization per k block (see
      ! get_X_from_overlap): vectors with S_k eigenvalue <= thresh are
      ! removed (zero column of X); blocks without such vectors keep the
      ! Lowdin X_k = S_k^-1/2 exactly as before. The dropped orbitals are
      ! moved to the top of the spectrum in RHF_SAO (shift_dropped_complex).
      !$omp parallel default(shared) private(k,b1,b2,i,Evec,Xb,Eval)
      allocate(Evec(nf,nf), Xb(nf,nf), Eval(nf))
      !$omp do schedule(static)
      do k = 1 , nk
        b1 = (k-1)*nf + 1 ; b2 = k*nf
        Evec(:,:) = over(b1:b2,b1:b2)
        call diagonalize_matrix_complex(nf,Evec,Eval)
        emin(k) = Eval(1)
        call canonical_X_complex(nf,Evec,Eval,thresh,Xb,ndrop(k))
        X(b1:b2,b1:b2) = Xb
      end do
      !$omp end do
      deallocate(Evec, Xb, Eval)
      !$omp end parallel

      if (sum(ndrop) > 0) then
        write(outfile,'(a)') ''
        write(outfile,'(a,es10.2,a)') ' Canonical orthogonalization (S_k eigenvalues <= ', thresh, ' removed)'
        do k = 1 , nk
          if (ndrop(k) > 0) write(outfile,'(a,i5,a,i4,a,i4,a,es12.4)') '   k block ', k, ' : ', ndrop(k), &
                                    ' of ', nf, ' removed, smallest eigenvalue ', emin(k)
        end do
        write(outfile,'(a,i8,a,i8)') '   total removed : ', sum(ndrop), ' of ', N
        write(outfile,'(a)') '   (the removed orbitals have zero coefficients and are placed at the top of the'
        write(outfile,'(a)') '    spectrum with an orbital energy of 1.0E+04, they never get occupied)'
        write(outfile,'(a)') ''
      end if

      if (c_details) then
        open(1,file=trim(tmp_file_name)//"/X.dat")
        write(1,'(15x,1000(i3,15x))') (i,i=1,size(X,1))
        do i = 1 , size(X,1)
          write(1,'(i3,6x,1000(f16.10,2x))') i ,  (X(i,o),o=1,size(X,1))
        end do
        close(1)
      end if

end subroutine get_X_from_overlap_complex_blocks

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


!=====================================================================!
! [Nov 2026] Canonical orthogonalization helpers.                      !
!                                                                     !
!  U, e : eigenvectors / eigenvalues of S (e ascending, from dsyev or  !
!         zheev).                                                     !
!  If no e <= thresh : X = U e^-1/2 U^T (Lowdin, unchanged behaviour). !
!  Else              : X(:,j) = U(:,j)/sqrt(e_j) for e_j > thresh and  !
!                      X(:,j) = 0 for the removed vectors.            !
!  X^T S X = 1 on the kept space, 0 on the removed one. With a zero    !
!  column j, F' = X^T F X has a zero row and column j; shift_dropped_* !
!  puts lindep_shift on F'(j,j) so the removed orbital is the highest  !
!  one (C = X e_j = 0, it is never occupied and gives zero integrals). !
!=====================================================================!
subroutine canonical_X_real(N,U,e,thresh,X,ndrop)

      implicit none
      integer         , intent(in)  :: N
      double precision, intent(in)  :: U(N,N) , e(N) , thresh
      double precision, intent(out) :: X(N,N)
      integer         , intent(out) :: ndrop
      double precision              :: d(N)
      integer                       :: j

      ndrop = count(e <= thresh)
      if (ndrop == 0) then
        d = 1.d0 / dsqrt(e)
        call ADAt(N,U,d,X)
      else
        do j = 1 , N
          if (e(j) > thresh) then
            X(:,j) = U(:,j) / dsqrt(e(j))
          else
            X(:,j) = 0.d0
          end if
        end do
      end if

end subroutine canonical_X_real

subroutine canonical_X_complex(N,U,e,thresh,X,ndrop)

      use constants_module
      implicit none
      integer         , intent(in)  :: N
      complex(dpc)    , intent(in)  :: U(N,N)
      double precision, intent(in)  :: e(N) , thresh
      complex(dpc)    , intent(out) :: X(N,N)
      integer         , intent(out) :: ndrop
      double precision              :: d(N)
      integer                       :: j

      ndrop = count(e <= thresh)
      if (ndrop == 0) then
        d = 1.d0 / dsqrt(e)
        call ADAt_complex(N,U,d,X)
      else
        do j = 1 , N
          if (e(j) > thresh) then
            X(:,j) = U(:,j) / dsqrt(e(j))
          else
            X(:,j) = (0.0_dp, 0.0_dp)
          end if
        end do
      end if

end subroutine canonical_X_complex

subroutine shift_dropped_real(N,X,Fp)

      implicit none
      integer         , intent(in)    :: N
      double precision, intent(in)    :: X(N,N)
      double precision, intent(inout) :: Fp(N,N)
      double precision, parameter     :: lindep_shift = 1.d4
      integer                         :: j

      do j = 1 , N
        if (all(X(:,j) == 0.d0)) Fp(j,j) = lindep_shift
      end do

end subroutine shift_dropped_real

subroutine shift_dropped_complex(N,X,Fp)

      use constants_module
      implicit none
      integer         , intent(in)    :: N
      complex(dpc)    , intent(in)    :: X(N,N)
      complex(dpc)    , intent(inout) :: Fp(N,N)
      double precision, parameter     :: lindep_shift = 1.d4
      integer                         :: j

      do j = 1 , N
        if (all(X(:,j) == (0.0_dp, 0.0_dp))) Fp(j,j) = cmplx(lindep_shift, 0.0_dp, kind=dpc)
      end do

end subroutine shift_dropped_complex
