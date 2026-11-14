!=====================================================================!
!                                                                     !
!              2D toroidal : KINETIC ENERGY INTEGRALS                 !
!                                                                     !
!   x,y : Clifford Gaussians (periodic),  z : real Gaussian           !
!   T = -1/2 ( Dx Sy Sz + Sx Dy Sz + Sx Sy Dz )   (Laplacian on ket)  !
!   periodic d : kinetic_1D_clifford,  real d : kinetic_1D_real       !
!                                                                     !
!   kinetic_integral_toroidal_2D (l1,l2,r1,r2,AO1,AO2,K_t) one pair   !
!   kinetic_matrix_toroidal_2D   (nBas,AO,K)               matrix     !
!=====================================================================!

subroutine kinetic_integral_toroidal_2D(l1,l2,r1,r2,AO1,AO2,K_t)

      use torus_init
      use classification_ERI
      use clifford_tools

      implicit none

      integer           ,intent(in)  :: l1(3) , l2(3)
      double precision  ,intent(in)  :: r1(3) , r2(3)
      type(ERI_function),intent(in)  :: AO1 , AO2
      double precision  ,intent(out) :: K_t

      integer          :: i , j , d
      double precision :: alpha , beta , pe(2) , cen(2) , av(2) , Lcell(2) , Aarg , kd
      double precision :: sd(3) , dd(3) , spa , cpa , spb , cpb , esum , mu , zP , kz
      double precision :: R(0:NMAX_CLIFF,0:NMAX_CLIFF)

      av    = (/ ax , ay /)
      Lcell = (/ Lx , Ly /)
      K_t = 0.d0
      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        do j = 1 , size(AO2%exponent)
          beta = AO2%exponent(j)
          call clifford_pair(2, alpha, beta, r1, r2, pe, cen)
          do d = 1 , 2
            Aarg = 2.d0*pe(d)/(av(d)*av(d))
            kd   = dexp(-2.d0*(alpha + beta - pe(d))/(av(d)*av(d)))
            spa = dsin(av(d)*(cen(d)-r1(d))) ; cpa = dcos(av(d)*(cen(d)-r1(d)))
            spb = dsin(av(d)*(cen(d)-r2(d))) ; cpb = dcos(av(d)*(cen(d)-r2(d)))
            call clifford_bessel_ratios(min(NMAX_CLIFF, l1(d)+l2(d)+2), Aarg, R)
            sd(d) = kd * overlap_1D_clifford(l1(d), l2(d), spa, cpa, spb, cpb, 1.d0/av(d), R, Lcell(d))
            dd(d) = kd * kinetic_1D_clifford(l1(d), l2(d), beta, av(d), spa, cpa, spb, cpb, R, Lcell(d))
          end do
          esum = alpha + beta ; mu = alpha*beta/esum
          zP   = (alpha*r1(3) + beta*r2(3))/esum
          kz   = dexp(-mu*(r1(3)-r2(3))**2)
          sd(3) = kz * overlap_1D_real(l1(3), l2(3), esum, zP-r1(3), zP-r2(3))
          dd(3) = kz * kinetic_1D_real(l1(3), l2(3), beta, esum, zP-r1(3), zP-r2(3))
          K_t = K_t - 0.5d0 * AO1%coefficient(i)*AO2%coefficient(j) * &
                ( dd(1)*sd(2)*sd(3) + sd(1)*dd(2)*sd(3) + sd(1)*sd(2)*dd(3) )
        end do
      end do

end subroutine kinetic_integral_toroidal_2D

subroutine kinetic_matrix_toroidal_2D(nBas,AO,K)

      use files
      use classification_ERI
      use keywords
      use unitcell_module
      use clifford_tools, only : init_clifford, orbital_to_lxyz

      implicit none

      integer           ,intent(in)  :: nBas
      type(ERI_function),intent(in)  :: AO(nBas)
      double precision  ,intent(out) :: K(nBas,nBas)

      integer          :: i , j , lxyz(3,nBas)
      double precision :: r1(3) , r2(3)

      call init_clifford()
      do i = 1 , nBas
        call orbital_to_lxyz(AO(i)%orbital, lxyz(:,i))
      end do

      K = 0.d0
      !$omp parallel do default(shared) private(i,j,r1,r2) schedule(dynamic,4)
      do i = 1 , nfuc
        r1 = (/ AO(i)%x , AO(i)%y , AO(i)%z /)
        do j = 1 , nBas
          r2 = (/ AO(j)%x , AO(j)%y , AO(j)%z /)
          call kinetic_integral_toroidal_2D(lxyz(:,i),lxyz(:,j),r1,r2,AO(i),AO(j),K(i,j))
        end do
      end do
      !$omp end parallel do

      do i = 1 , nfuc
        do j = 1 , nBas
          if (abs(K(i,j)) < 1e-14) K(i,j) = 0.d0
        end do
      end do
      call fill_with_symmetry(nBas, K)

      if (c_K) then
        write(outfile,'(a)') ""
        write(outfile,'(a)') "------------------------------------------"
        write(outfile,'(a)') "The kinetic matrix (2D)   "
        write(outfile,'(a)') "------------------------------------------"
        do i = 1 , nBas - 1
          do j = i , nBas
            write(outfile,'(i3,i3, 6x,1000(f16.12,2x))') i , j, K(i,j)
          end do
        end do
      end if

end subroutine kinetic_matrix_toroidal_2D
