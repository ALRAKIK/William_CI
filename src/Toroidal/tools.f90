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

      ! ----------------------------------------------------------------!
      ! ----------------------------------------------------------------!










subroutine symmetry_of_integrals(nf,fpuc,mat_tmp,mat)

      use unitcell_module

      implicit none 

      integer                      :: i , j
      integer                      :: i_cell  , j_cell  , k_cell 
      integer                      :: i_shift , j_shift , k_shift 
      integer                      :: m , n , l
      integer                      :: nf
      integer                      :: fpuc
      integer                      :: func_i , func_j
      integer                      :: equivalent_j
      
      double precision             :: mat_tmp(nf,nf)
      double precision,intent(out) :: mat(nf,nf)

      !-----------------------------------------------------------------!

      mat(:,:) = 0.d0

      ! First, copy what you already computed (first fpuc rows)
      do i = 1, fpuc
        do j = 1, nf
          mat(i,j) = mat_tmp(i,j)
        end do
      end do

      ! Now fill the rest using symmetry
      
      do i = fpuc + 1, nf
        ! Find which unit cell and function this corresponds to
        call find_cell_and_function(i, fpuc, i_cell, j_cell, k_cell, func_i)
        
        do j = 1, nf
          ! Find which unit cell and function j corresponds to
          call find_cell_and_function(j, fpuc, i_shift, j_shift, k_shift, func_j)

          ! The overlap(i,j) should equal overlap(func_i, equivalent_j)
          ! where equivalent_j is the function in the reference unit cell that corresponds to j
          ! relative to i's unit cell
          
          ! Calculate the relative shift
          m = modulo(i_shift - i_cell, nx)
          n = modulo(j_shift - j_cell, ny)  
          l = modulo(k_shift - k_cell, nz)
          
          ! Find the equivalent function in the first fpuc rows
          equivalent_j = m * ny * nz * fpuc + n * nz * fpuc + l * fpuc + func_j
          
          ! Use the symmetry: S(i,j) = S(func_i, equivalent_j)

          if (equivalent_j <= nf) then
            mat(i,j) = mat_tmp(func_i, equivalent_j)
          endif
        end do
      end do

      !-----------------------------------------------------------------!

end subroutine 

subroutine symmetry_of_integrals_ERI(nf, fpuc, eri_tmp, eri)

      use unitcell_module
      use omp_lib

      implicit none 

      integer                      :: i, j, k, l
      integer                      :: i_cell, j_cell, k_cell
      integer                      :: i_shift, j_shift, k_shift
      integer                      :: m, n, p
      integer                      :: nf
      integer                      :: fpuc
      integer                      :: func_i, func_j, func_k, func_l
      integer                      :: equiv_j, equiv_k, equiv_l

      double precision             :: eri_tmp(fpuc, nf, nf, nf)
      double precision,intent(out) ::     eri(nf, nf, nf, nf)
      double precision             :: val 

      ! ----------------------    Time     ---------------------------- !
      integer                        :: days, hours, minutes, seconds , time 
      double precision               :: start,end
      !-----------------------------------------------------------------!


      !-----------------------------------------------------------------!

      eri(:,:,:,:) = 0.d0
    
      ! First, copy what you already computed (first index up to fpuc)

      call cpu_time(start)

      do i = 1, fpuc
        do j = 1, nf
          do k = 1, nf
            do l = 1, nf
              if (dabs(eri_tmp(i,j,k,l)) > 1.d-12) eri(i,j,k,l) = eri_tmp(i,j,k,l)
            end do
          end do
        end do
      end do

      call cpu_time(end)

      time = int(end - start)
      days = (time/86400)
      hours=mod(time,86400)/3600
      minutes=mod(mod(time,86400),3600)/60
      seconds=mod(mod(mod(time,86400),3600),60)

      write(*,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') 'CPU time for The first loop of The Translational and 8 fold symmetry = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec"

      call cpu_time(start)

      !$OMP PARALLEL DO PRIVATE(i,j,k,l,val) SCHEDULE(static)
      do i = 1, nf
        do j = i, nf
          do k = 1, nf
            do l = k, nf
              val = eri(i,j,k,l)
              if (val /= 0.d0) then
                eri(i,j,l,k) = val
                eri(j,i,k,l) = val
                eri(j,i,l,k) = val
                eri(k,l,i,j) = val
                eri(k,l,j,i) = val
                eri(l,k,i,j) = val
                eri(l,k,j,i) = val
              end if
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      call cpu_time(end)

      time = int(end - start)
      days = (time/86400)
      hours=mod(time,86400)/3600
      minutes=mod(mod(time,86400),3600)/60
      seconds=mod(mod(mod(time,86400),3600),60)

      write(*,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') 'CPU time for The second loop of The Translational and 8 fold symmetry = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec"

      call cpu_time(start)

      ! Now fill the rest using symmetry
    
      do i = fpuc + 1, nf
        ! Find which unit cell and function i corresponds to
        call find_cell_and_function(i, fpuc, i_cell, j_cell, k_cell, func_i)
        
        do j = 1, nf
          ! Find which unit cell and function j corresponds to  
          call find_cell_and_function(j, fpuc, i_shift, j_shift, k_shift, func_j)
          
          ! Calculate relative shift for j
          m = modulo(i_shift - i_cell, nx)
          n = modulo(j_shift - j_cell, ny)  
          p = modulo(k_shift - k_cell, nz)
          
          ! Find equivalent j in reference frame
          equiv_j = m * ny * nz * fpuc + n * nz * fpuc + p * fpuc + func_j
          
          do k = 1, nf
            ! Find which unit cell and function k corresponds to  
            call find_cell_and_function(k, fpuc, i_shift, j_shift, k_shift, func_k)
            
            ! Calculate relative shift for k
            m = modulo(i_shift - i_cell, nx)
            n = modulo(j_shift - j_cell, ny)  
            p = modulo(k_shift - k_cell, nz)
            
            ! Find equivalent k in reference frame
            equiv_k = m * ny * nz * fpuc + n * nz * fpuc + p * fpuc + func_k
            
            do l = 1, nf

              ! Find which unit cell and function l corresponds to  

              call find_cell_and_function(l, fpuc, i_shift, j_shift, k_shift, func_l)
              
              ! Calculate relative shift for l
              m = modulo(i_shift - i_cell, nx)
              n = modulo(j_shift - j_cell, ny)  
              p = modulo(k_shift - k_cell, nz)
              
              ! Find equivalent l in reference frame

              equiv_l = m * ny * nz * fpuc + n * nz * fpuc + p * fpuc + func_l
              
              ! Use the symmetry: (i,j,k,l) = (func_i, equiv_j, equiv_k, equiv_l)

              if (equiv_j <= nf .and. &
                  equiv_k <= nf .and. &
                  equiv_l <= nf) then
                  eri(i,j,k,l) = eri(func_i, equiv_j, equiv_k, equiv_l)
              endif

            end do
          end do
        end do
      end do

      call cpu_time(end)

      time = int(end - start)
      days = (time/86400)
      hours=mod(time,86400)/3600
      minutes=mod(mod(time,86400),3600)/60
      seconds=mod(mod(mod(time,86400),3600),60)

      write(*,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') 'CPU time for The Third loop of The Translational and 8 fold symmetry = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec"

      !-----------------------------------------------------------------!

end subroutine symmetry_of_integrals_ERI

subroutine find_cell_and_function(index, fpuc, i_cell, j_cell, k_cell, func_index)

      use unitcell_module
      implicit none

      integer, intent(in) :: index, fpuc
      integer, intent(out) :: i_cell, j_cell, k_cell, func_index
  
      integer :: temp
  
      ! Function index within unit cell (1 to fpuc)
      func_index = mod(index - 1, fpuc) + 1
  
      ! Remove function part and get cell index
      temp = (index - 1) / fpuc
  
      ! Decompose cell index into i,j,k coordinates
      k_cell = mod(temp, nz)
      temp = temp / nz
      j_cell = mod(temp, ny) 
      i_cell = temp / ny
  
end subroutine find_cell_and_function

      !-----------------------------------------------------------------!
      !-----------------------------------------------------------------!


subroutine symmetry_of_integrals_ERI_mod(nf, fpuc, eri_tmp, eri)

      use unitcell_module
      use omp_lib

      implicit none 

      integer                      :: i, j, k, l
      integer                      :: nf, fpuc
      integer                      :: func_i
      integer                      :: idx, ii, m, n, p
      double precision             :: eri_tmp(fpuc, nf, nf, nf)
      double precision,intent(out) ::     eri(nf, nf, nf, nf)
      double precision             :: val 

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
    
      ! Copy computed integrals (first index up to fpuc)
      do i = 1, fpuc
        do j = 1, nf
          do k = 1, nf
            do l = 1, nf
              if (dabs(eri_tmp(i,j,k,l)) > 1.d-12) eri(i,j,k,l) = eri_tmp(i,j,k,l)
            end do
          end do
        end do
      end do

      ! Precompute equivalence mapping table
      do i = 1, nf
        do j = 1, nf
          m = modulo(icell(j) - icell(i), nx)
          n = modulo(jcell(j) - jcell(i), ny)
          p = modulo(kcell(j) - kcell(i), nz)
          equiv(i,j) = m * ny * nz * fpuc + n * nz * fpuc + p * fpuc + func_idx(j)
        end do
      end do

      ! 8-fold symmetry expansion (keep as is)
      !$OMP PARALLEL DO PRIVATE(i,j,k,l,val) SCHEDULE(static)
      do i = 1, nf
        do j = i, nf
          do k = 1, nf
            do l = k, nf
              val = eri(i,j,k,l)
              if (val /= 0.d0) then
                eri(i,j,l,k) = val
                eri(j,i,k,l) = val
                eri(j,i,l,k) = val
                eri(k,l,i,j) = val
                eri(k,l,j,i) = val
                eri(l,k,i,j) = val
                eri(l,k,j,i) = val
              end if
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO


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

end subroutine symmetry_of_integrals_ERI_mod

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