subroutine spectral_function(A_w, eigve, eta, eigva, e_min, e_max, delta, nBas, heffbas)

      implicit none
      
      ! Input variables
      double precision     , intent(in)               :: eigve(heffbas, heffbas), eigva(heffbas), eta, e_min, e_max, delta
      double precision     , intent(inout)            :: A_w(int((e_max-e_min)/delta),2)
      integer              , intent(in)               :: nBas, heffbas

      ! Local variables
      double precision     , allocatable              :: omega(:), eigve_squared(:,:), result(:)
      complex              , allocatable              :: z(:)
      double precision     , parameter                :: PI = 4*atan(1.0d0)
      integer                                         :: i, j, step

      step = int((e_max-e_min)/delta)
      allocate(omega(step))
      allocate(z(step))
      allocate(result(step))
      allocate(eigve_squared(Size(eigve(1,:)),size(eigve(1,:))))

      do i = 1,step
         omega(i)=(i-1)*delta+e_min
      enddo
      eigve_squared = eigve**2d0
      A_w = 0.0d0
      
      do i = 1, nBas
         do j = 1, heffbas
            A_w(:,2) = A_w(:,2) + 1/PI*eta*eigve_squared(i,j)/((eigva(j)-omega)**2+eta**2)
         enddo
      enddo
      A_w(:,1) = omega

end subroutine spectral_function