!=====================================================================!
!                                                                     !
!                 1D toroidal : OVERLAP INTEGRALS                     !
!                                                                     !
!   x : Clifford Gaussian   g^i = [sin(a(x-A))/a]^i                   !
!                                 exp(-2alpha/a^2 (1-cos a(x-A)))     !
!   y,z : real Gaussians                                              !
!                                                                     !
!   S = sum_prim c_i c_j  kx Sx * ky Sy * kz Sz                        !
!                                                                     !
!   overlap_integral_toroidal_1D (l1,l2,r1,r2,AO1,AO2,S)  one pair    !
!   overlap_matrix_toroidal_1D   (nBas,AO,overlap)        full matrix !
!                                                                     !
!   any l (s ... g); l1, l2 = (lx,ly,lz) from orbital_to_lxyz         !
!=====================================================================!

subroutine overlap_integral_toroidal_1D(l1,l2,r1,r2,AO1,AO2,S)

      use torus_init
      use classification_ERI
      use clifford_tools

      implicit none

      integer           ,intent(in)    :: l1(3) , l2(3)     ! (lx,ly,lz) of AO1, AO2
      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: S

      integer                          :: i , j , nmax
      double precision                 :: alpha , beta , const
      double precision                 :: x1 , x2 , y1 , y2 , z1 , z2 , X , Y , Z
      double precision                 :: kx , ky , kz , sx , sy , sz
      double precision                 :: ax2 , inv_ax , px , xp , A
      double precision                 :: spa , spb , cpa , cpb
      double precision                 :: R(0:NMAX_CLIFF,0:NMAX_CLIFF)
      double precision                 :: albe , inv_albe , mu , yp , zp

      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1)
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X  = x1 - x2 ; Y = y1 - y2 ; Z = z1 - z2

      ax2    = ax * ax
      inv_ax = 1.d0 / ax

      ! highest Bessel order needed: s + m <= lx1 + lx2
      nmax = min(NMAX_CLIFF, l1(1) + l2(1))

      S = 0.d0

      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        do j = 1 , size(AO2%exponent)
          beta  = AO2%exponent(j)
          const = AO1%coefficient(i) * AO2%coefficient(j)

          ! ---------- Clifford Gaussian (x) ---------- !
          call bary_exponent_x        (alpha,beta,X,px)
          call bary_center_toroidal_x (alpha,beta,x1,x2,xp)

          A   = 2.d0*px/ax2
          kx  = dexp(-2.d0*(alpha + beta - px)/ax2)

          spa = dsin(ax*(xp-x1)) ; cpa = dcos(ax*(xp-x1))
          spb = dsin(ax*(xp-x2)) ; cpb = dcos(ax*(xp-x2))

          call clifford_bessel_ratios(nmax, A, R)
          sx = overlap_1D_clifford(l1(1), l2(1), spa, cpa, spb, cpb, inv_ax, R, Lx)

          ! ---------- Real Gaussians (y, z) ---------- !
          albe     = alpha + beta
          inv_albe = 1.d0 / albe
          mu       = alpha * beta * inv_albe
          yp       = (alpha*y1 + beta*y2) * inv_albe
          zp       = (alpha*z1 + beta*z2) * inv_albe

          ky = dexp(-mu*Y*Y)
          kz = dexp(-mu*Z*Z)

          sy = overlap_1D_real(l1(2), l2(2), albe, yp-y1, yp-y2)
          sz = overlap_1D_real(l1(3), l2(3), albe, zp-z1, zp-z2)

          S = S + const * kx * ky * kz * sx * sy * sz

        end do
      end do

end subroutine overlap_integral_toroidal_1D

!=====================================================================!
!  Overlap matrix: rows i = 1..nfuc (unit cell), all columns j,       !
!  the rest from the translational symmetry (fill_with_symmetry)      !
!=====================================================================!
subroutine overlap_matrix_toroidal_1D(nBas,AO,overlap)

      use files
      use torus_init
      use classification_ERI
      use keywords
      use unitcell_module
      use clifford_tools, only : init_clifford, orbital_to_lxyz

      implicit none

      integer           ,intent(in)  :: nBas
      type(ERI_function),intent(in)  :: AO(nBas)
      double precision  ,intent(out) :: overlap(nBas,nBas)

      integer                        :: i , j
      integer                        :: lxyz(3,nBas)
      double precision               :: r1(3) , r2(3)

      !-----------------------------------------------------------------!

      call init_clifford()

      do i = 1 , nBas
        call orbital_to_lxyz(AO(i)%orbital, lxyz(:,i))
      end do

      overlap(:,:) = 0.d0

      do i = 1 , nfuc
        r1 = (/ AO(i)%x , AO(i)%y , AO(i)%z /)
        do j = 1 , nBas
          r2 = (/ AO(j)%x , AO(j)%y , AO(j)%z /)
          call overlap_integral_toroidal_1D(lxyz(:,i),lxyz(:,j),r1,r2,AO(i),AO(j),overlap(i,j))
        end do
      end do

      do i = 1 , nfuc
        do j = 1 , nBas
          if (abs(overlap(i,j)) < 1e-14) overlap(i,j) = 0.d0
        end do
      end do

      call fill_with_symmetry(nBas, overlap)

      if (c_OV) then
        write(outfile,'(a)') ""
        write(outfile,'(a)') "------------------------------------"
        write(outfile,'(a)') "The overlap matrix                  "
        write(outfile,'(a)') "------------------------------------"
        write(outfile,'(a)') ""
        do i = 1 , nBas - 1
          do j = i , nBas
            write(outfile,'(i3,i3, 6x,1000(f16.12,2x))') i , j, overlap(i,j)
          end do
        end do
      end if

end subroutine overlap_matrix_toroidal_1D
