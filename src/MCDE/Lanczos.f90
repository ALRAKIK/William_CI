subroutine Lanczos_Haydock(heff_matrix, nSpace1, nSpace3, heffbas, it, a, b)

      ! Input variables
      integer         , intent(in)                        :: nSpace1, nSpace3, heffbas
      double precision, intent(in)                        :: heff_matrix(heffbas, heffbas)
      integer         , intent(in)                        :: it

      ! Output variables     
      double precision, allocatable, intent(out)          :: a(:)
      double precision, allocatable, intent(out)          :: b(:)

      ! Local variables
      double precision, allocatable          :: s0(:,:), w0(:,:), w(:,:), alpha(:,:), beta(:,:), v_prec(:,:), v_curr(:,:)
      double precision, allocatable          :: Krylov(:,:,:)
      double precision, allocatable          :: s(:,:)
      integer                                :: i,j,f

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '             |  Lanczos algorithm  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)
      
      allocate(s0(heffbas,1))
      s0=0
      allocate(w(heffbas,1))
      allocate(alpha(heffbas,1))
      allocate(beta(heffbas,1))
      do i = 1, nSpace1
         s0(i,1) = 1.0d0
      enddo
      allocate(Krylov(heffbas,heffbas,1))
      s0=s0/norm2(s0)
      
      Krylov(1,:,1)=s0(:,1)
      w0=Transpose(MATMUL(heff_matrix,s0))
      ! write(*,*)'s0 = ',s0
      ! write(*,*)
      ! write(*,*)'w0 = ',w0
      ! write(*,*)

      alpha(1,1) = dot_product(s0(:,1),w0(1,:))
      beta(1,1)=0.0d0
      w=TRANSPOSE(w0)-alpha(1,1)*s0
      beta(2,1)=NORM2(w)
      Krylov(2,:,:)=w/beta(2,1)

      do i = 2, (it-1)
         v_prec=Krylov(i-1,:,:)
         v_curr=Krylov(i,:,:)
         w0=Transpose(matmul(heff_matrix,v_curr))
         alpha(i,1)=dot_product((v_curr(:,1)),w0(1,:)) 
         w=TRANSPOSE(w0)-(alpha(i,1)*v_curr+beta(i,1)*v_prec)
         do j = 1, i
            w = w - matmul(matmul(Krylov(j,:,:),Transpose(w)),Krylov(j,:,:))
         enddo
         beta(i+1,1) = norm2(w)
         Krylov(i+1,:,:) = w/(beta(i+1,1))
         write(*,*)'b = ', beta(i+1,1)
         if (abs(beta(i+1,1)) < 1.0d-12) then
            EXIT 
         endif
      enddo
      f = i

      s=matmul(matmul(TRANSPOSE(Krylov(f,:,:)),heff_matrix),(Krylov(f,:,:)))
      alpha(f,1)=s(1,1)
      allocate(a(f))
      allocate(b(f))
      a=alpha(:f,1)
      b=beta(2:f+1,1)     
   
end subroutine Lanczos_Haydock

subroutine Lanczos_spectral_function(e_min,e_max,delta,a,b,A_w,eta)

      ! Input variables
      double precision, intent(in)           :: e_min, e_max, delta
      double precision, intent(in)           :: a(:), b(:)
      double precision, intent(in)           :: eta

      ! Output variables     
      double precision, allocatable, intent(out)          :: A_w(:,:)

      ! Local variables
      integer                                             :: steps
      integer                                             :: i
      double precision, allocatable                       :: omega(:)
      complex, allocatable                                :: z(:), result(:)
      double precision, parameter                         :: PI = 4*atan(1.0d0)

      steps = int((e_max-e_min)/delta)
      allocate(omega(steps))
      allocate(z(steps))
      allocate(result(steps))
      allocate(A_w(steps,2))
      do i = 1,steps
         omega(i)=(i-1)*delta+e_min
      enddo
      z = cmplx(omega,eta)
      write(*,*)'z = ', z
      result = z-a(size(a))
      do i = size(b), 1, -1
         result = z - a(i) - b(i)**2 / result
      end do
      result = -PI*AIMAG(1/result)
      A_w(:,1) = omega
      A_w(:,2) = result

   end subroutine Lanczos_spectral_function