subroutine bary_center_toroidal_x(e1,e2,r1,r2,rp)

      use torus_init
      use HeavisideModule

      implicit none
      double precision, intent(in)  :: e1, e2
      double precision, intent(in)  :: r1, r2
      double precision, intent(out) :: rp

      ! local !

      double precision,parameter    :: epsilon = 1.d-6
      double precision              :: val 

      ! ----------------------------------------------------------------!

      if (dabs(r1-r2) < 0.5d0*Lx + epsilon .and. dabs(r1-r2) > 0.5d0*Lx - epsilon) then
        if (dabs(e1-e2) < 1.d-10) then 
          rp = 0.5d0 * ( r1 + r2 )
        else 
          if (e1 > e2) then 
            rp = r1 
          else if (e2 > e1) then 
            rp = r2 
          else 
            rp = 0.5d0 * ( r1 + r2 )
          end if 
        end if
    
      else if (dabs(e1-e2) < 1.d-10) then 
        if (dabs(r1-r2) > 0.5*Lx + epsilon ) then
          rp = 0.5d0 * ( r1 + r2 ) - 0.5 * Lx
        else 
          rp = 0.5d0 * ( r1 + r2 )
        end if 
      else 
        val = e1*dcos(ax*r1)+e2*dcos(ax*r2)
        rp = datan((e1*dsin(ax*r1)+e2*dsin(ax*r2))/val)/ax + 0.5d0 * Lx * Heaviside(-val)
      end if
    
end subroutine bary_center_toroidal_x

      ! ----------------------------------------------------------------!
      ! ----------------------------------------------------------------!

subroutine bary_center_toroidal_y(e1,e2,r1,r2,rp)

      use torus_init
      use HeavisideModule

      implicit none
      double precision, intent(in)  :: e1, e2
      double precision, intent(in)  :: r1, r2
      double precision, intent(out) :: rp

      ! local !

      double precision,parameter    :: epsilon = 1.d-6
      double precision              :: val 

      ! ----------------------------------------------------------------!

      if (dabs(r1-r2) < 0.5d0*Ly + epsilon .and. dabs(r1-r2) > 0.5d0*Ly - epsilon) then
        if (dabs(e1-e2) < 1.d-10) then 
          rp = 0.5d0 * ( r1 + r2 )
        else 
          if (e1 > e2) then 
            rp = r1 
          else if (e2 > e1) then 
            rp = r2 
          else 
            rp = 0.5d0 * ( r1 + r2 )
          end if 
        end if
    
      else if (dabs(e1-e2) < 1.d-10) then 
        if (dabs(r1-r2) > 0.5*Ly + epsilon ) then
          rp = 0.5d0 * ( r1 + r2 ) - 0.5 * Ly
        else 
          rp = 0.5d0 * ( r1 + r2 )
        end if 
      else 
        val = e1*dcos(ay*r1)+e2*dcos(ay*r2)
        rp = datan((e1*dsin(ay*r1)+e2*dsin(ay*r2))/val)/ay + 0.5d0 * Ly * Heaviside(-val)
      end if
    
end subroutine bary_center_toroidal_y

      ! ----------------------------------------------------------------!
      ! ----------------------------------------------------------------!

subroutine bary_center_toroidal_z(e1,e2,r1,r2,rp)

      use torus_init
      use HeavisideModule

      implicit none
      double precision, intent(in)  :: e1, e2
      double precision, intent(in)  :: r1, r2
      double precision, intent(out) :: rp

      ! local !

      double precision,parameter    :: epsilon = 1.d-6
      double precision              :: val 

      ! ----------------------------------------------------------------!

      if (dabs(r1-r2) < 0.5d0*Lz + epsilon .and. dabs(r1-r2) > 0.5d0*Lz - epsilon) then
        if (dabs(e1-e2) < 1.d-10) then 
          rp = 0.5d0 * ( r1 + r2 )
        else 
          if (e1 > e2) then 
            rp = r1 
          else if (e2 > e1) then 
            rp = r2 
          else 
            rp = 0.5d0 * ( r1 + r2 )
          end if 
        end if
    
      else if (dabs(e1-e2) < 1.d-10) then 
        if (dabs(r1-r2) > 0.5*Lz + epsilon ) then
          rp = 0.5d0 * ( r1 + r2 ) - 0.5 * Lz
        else 
          rp = 0.5d0 * ( r1 + r2 )
        end if 
      else 
        val = e1*dcos(az*r1)+e2*dcos(az*r2)
        rp = datan((e1*dsin(az*r1)+e2*dsin(az*r2))/val)/az + 0.5d0 * Lz * Heaviside(-val)
      end if
    
end subroutine bary_center_toroidal_z

      ! ----------------------------------------------------------------!
      ! ----------------------------------------------------------------!


subroutine bary_exponent_x(e1,e2,x,p)

      use torus_init
      use HeavisideModule

      implicit none
      double precision, intent(in)  :: e1, e2
      double precision, intent(in)  :: x
      double precision, intent(out) :: p

      ! local !

      double precision,parameter    :: epsilon = 1.d-6
      double precision              :: e12 , e22

      ! ----------------------------------------------------------------!


      e12 = e1 * e1 
      e22 = e2 * e2 

      
      if (dabs(x) > 0.5d0*Lx-epsilon .and. dabs(x) < 0.5d0*Lx+epsilon) then 
        if (dabs(e1 - e2) < 1.d-8) then
          p = 0.d0 
        else 
          p = dabs(e1 - e2)
        end if 
      else 
        p = e12 + e22 + 2.d0 * e1 * e2 * dcos(ax*x)
        p = dsqrt(dabs(p))
      end if

      if (p > e1 + e2) p = e1+e2

      if (e1+e2-p < 1.d-14) p = e1 + e2 

end subroutine bary_exponent_x

      ! ----------------------------------------------------------------!
      ! ----------------------------------------------------------------!

subroutine bary_exponent_y(e1,e2,y,p)

      use torus_init
      use HeavisideModule

      implicit none
      double precision, intent(in)  :: e1, e2
      double precision, intent(in)  :: y
      double precision, intent(out) :: p

      ! local !

      double precision,parameter    :: epsilon = 1.d-6
      double precision              :: e12 , e22

      ! ----------------------------------------------------------------!


      e12 = e1 * e1 
      e22 = e2 * e2 

      
      if (dabs(y) > 0.5d0*Ly-epsilon .and. dabs(y) < 0.5d0*Ly+epsilon) then 
        if (dabs(e1 - e2) < 1.d-8) then
          p = 0.d0 
        else 
          p = dabs(e1 - e2)
        end if 
      else 
        p = e12 + e22 + 2.d0 * e1 * e2 * dcos(ay*y)
        p = dsqrt(dabs(p))
      end if

      if (p > e1 + e2) p = e1+e2

      if (e1+e2-p < 1.d-14) p = e1 + e2 

end subroutine bary_exponent_y

      ! ----------------------------------------------------------------!
      ! ----------------------------------------------------------------!

subroutine bary_exponent_z(e1,e2,z,p)

      use torus_init
      use HeavisideModule

      implicit none
      double precision, intent(in)  :: e1, e2
      double precision, intent(in)  :: z
      double precision, intent(out) :: p

      ! local !

      double precision,parameter    :: epsilon = 1.d-6
      double precision              :: e12 , e22

      ! ----------------------------------------------------------------!


      e12 = e1 * e1 
      e22 = e2 * e2 

      
      if (dabs(z) > 0.5d0*Lz-epsilon .and. dabs(z) < 0.5d0*Lz+epsilon) then 
        if (dabs(e1 - e2) < 1.d-8) then
          p = 0.d0 
        else 
          p = dabs(e1 - e2)
        end if 
      else 
        p = e12 + e22 + 2.d0 * e1 * e2 * dcos(az*z)
        p = dsqrt(dabs(p))
      end if

      if (p > e1 + e2) p = e1+e2

      if (e1+e2-p < 1.d-14) p = e1 + e2 
      

end subroutine bary_exponent_z

      !-----------------------------------------------------------------!
      !-----------------------------------------------------------------!