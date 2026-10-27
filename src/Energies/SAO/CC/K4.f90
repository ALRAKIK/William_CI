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


! CHEMIST integral order (p q | r s): ks = kp-kq+kr.
integer function cc_integral_k4(kp,kq,kr) result(ks)
    
      implicit none 
      integer, intent(in) :: kp,kq,kr
      integer             :: fourth_k

      ks = fourth_k(kp,kr,kq)

end function cc_integral_k4