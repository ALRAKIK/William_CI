integer function fourth_k(k1,k2,k3) result(k4)

      use unitcell_module
      implicit none 
      integer, intent(in) :: k1,k2,k3
      integer             :: x1,y1,z1,x2,y2,z2,x3,y3,z3

      x1 = modulo(k1-1,nx)
      y1 = modulo((k1-1)/nx,ny)
      z1 = (k1-1)/(nx*ny)

      x2 = modulo(k2-1,nx)
      y2 = modulo((k2-1)/nx,ny)
      z2 = (k2-1)/(nx*ny)

      x3 = modulo(k3-1,nx)
      y3 = modulo((k3-1)/nx,ny)
      z3 = (k3-1)/(nx*ny)

      k4 = modulo(x1+x2-x3,nx) + nx*modulo(y1+y2-y3,ny) + nx*ny*modulo(z1+z2-z3,nz) + 1

end function fourth_k


subroutine diis_coefficients(errors,c,ok)
          
      use constants_module, only: dp, dpc
      
      implicit none 

      ! input !
      complex(dpc), intent(in) :: errors(:,:)

      ! output ! 
      real(dp)    , intent(out) :: c(:)
      logical     , intent(out) :: ok
      
      ! local ! 
      real(dp) :: mat(size(c)+1,size(c)+1),rhs(size(c)+1),row(size(c)+1)
      real(dp) :: factor,scale,tmp
      integer  :: m,i,j,k,p



      m   = size(c)
      ok  = .false.
      c   = 0.0_dp
      mat = 0.0_dp
      rhs = 0.0_dp

      do i=1,m
        do j=1,m
          mat(i,j)=real(sum(conjg(errors(:,i))*errors(:,j)),dp)
        end do
      end do

      scale=maxval(abs(mat(:m,:m)))


      if (scale<=tiny(1.0_dp)) return

      mat(:m,:m)=mat(:m,:m)/scale

      do i=1,m
        mat(i,i)=mat(i,i)+1.0e-12_dp
      end do

      mat(:m,m+1)=-1.0_dp
      mat(m+1,:m)=-1.0_dp
      rhs(m+1)   =-1.0_dp

      do k=1,m+1
        p=k-1+maxloc(abs(mat(k:m+1,k)),dim=1)
        if (abs(mat(p,k))<1.0e-14_dp) return
        row=mat(k,:); mat(k,:)=mat(p,:); mat(p,:)=row
        tmp=rhs(k); rhs(k)=rhs(p); rhs(p)=tmp
        do i=k+1,m+1
          factor=mat(i,k)/mat(k,k)
          mat(i,k:m+1)=mat(i,k:m+1)-factor*mat(k,k:m+1)
          rhs(i)=rhs(i)-factor*rhs(k)
        end do
      end do

      do i=m+1,1,-1
        rhs(i)=(rhs(i)-sum(mat(i,i+1:m+1)*rhs(i+1:m+1)))/mat(i,i)
      end do

      c=rhs(:m)
      ok=all(abs(c)<100.0_dp).and.abs(sum(c)-1.0_dp)<1.0e-8_dp

end subroutine diis_coefficients