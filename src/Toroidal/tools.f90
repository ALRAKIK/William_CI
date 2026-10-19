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


subroutine symmetry_of_integrals_ERI(nf, fpuc, eri_tmp, eri)

      use unitcell_module
      use omp_lib

      implicit none 

      integer                      :: i, j, k, l
      integer                      :: nf, fpuc
      integer                      :: func_i
      integer                      :: idx, ii, m, n, p
      double precision             :: eri_tmp(fpuc, nf, nf, nf)
      double precision,intent(out) ::     eri(nf, nf, nf, nf)

      ! Precomputation arrays
      integer                      :: icell(nf), jcell(nf), kcell(nf), func_idx(nf)
      integer                      :: equiv(nf, nf)
      integer                      :: eqv(nf)
      integer                      :: el , ek
      
      ! Timing

      !-----------------------------------------------------------------!

      eri(:,:,:,:) = 0.d0
      
      ! Precompute cell and function indices
      do idx = 1, nf
        func_idx(idx) = mod(idx - 1, fpuc) + 1
        ii = (idx - 1) / fpuc
        kcell(idx) = mod(ii, nz)
        ii = ii / nz
        jcell(idx) = mod(ii, ny)
        icell(idx) = ii / ny
      end do
    
      
      eri(1:nfuc,:,:,:) = eri_tmp(:,:,:,:)

      ! Precompute equivalence mapping table
      do i = 1, nf
        do j = 1, nf
          m = modulo(icell(j) - icell(i), nx)
          n = modulo(jcell(j) - jcell(i), ny)
          p = modulo(kcell(j) - kcell(i), nz)
          equiv(i,j) = m * ny * nz * fpuc + n * nz * fpuc + p * fpuc + func_idx(j)
        end do
      end do

      !$OMP PARALLEL DO PRIVATE(i,func_i,eqv,l,el,k,ek,j) SCHEDULE(static)
      do i = fpuc + 1, nf
        func_i = func_idx(i)
        ! Copy equiv for this i into local array
        do j = 1, nf
          eqv(j) = equiv(i, j)
        end do
        do l = 1, nf
          el = eqv(l)
          do k = 1, nf
            ek = eqv(k)
            ! Innermost loop over j (contiguous stride 1)
            do j = 1, nf
              eri(i, j, k, l) = eri(func_i, eqv(j), ek, el)
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO
      !-----------------------------------------------------------------!

end subroutine symmetry_of_integrals_ERI

subroutine save_two_electron_integrals(fpuc, nf, int)
    
      implicit none
      integer, intent(in)          :: fpuc, nf
      double precision, intent(in) :: int(fpuc, nf, nf, nf)

      integer                      :: i, j, k, l

      ! Open file as unformatted, stream access → pure binary, no record overhead

      open(10, file="ERI_unique", form='unformatted', access='stream', status='replace')

      do i = 1, fpuc
          do j = 1, nf
              do k = 1, nf
                  do l = 1, nf
                      if (dabs(int(i,j,k,l)) > 1.d-12) then
                          ! Write indices and value directly as binary
                          write(10) i, j, k, l, int(i,j,k,l)
                      end if
                  end do
              end do
          end do
      end do

      close(10)

end subroutine save_two_electron_integrals