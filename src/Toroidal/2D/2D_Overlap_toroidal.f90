!=====================================================================!
!                                                                     !
!                 2D toroidal : OVERLAP INTEGRALS                     !
!                                                                     !
!   x,y : Clifford Gaussians (periodic),  z : real Gaussian           !
!   S = sum_prim c_i c_j  prod_d k_d S_d                              !
!                                                                     !
!   overlap_integral_toroidal_2D (l1,l2,r1,r2,AO1,AO2,S)  one pair    !
!   overlap_matrix_toroidal_2D   (nBas,AO,OV)             full matrix !
!=====================================================================!

subroutine overlap_integral_toroidal_2D(l1,l2,r1,r2,AO1,AO2,S)

      use torus_init
      use classification_ERI
      use clifford_tools

      implicit none

      integer           ,intent(in)  :: l1(3) , l2(3)
      double precision  ,intent(in)  :: r1(3) , r2(3)
      type(ERI_function),intent(in)  :: AO1 , AO2
      double precision  ,intent(out) :: S

      integer          :: i , j , d
      double precision :: alpha , beta , pe(2) , cen(2) , av(2) , Lcell(2) , Aarg , kd , sd(3)
      double precision :: esum , mu , zP
      double precision :: R(0:NMAX_CLIFF,0:NMAX_CLIFF)

      av    = (/ ax , ay /)
      Lcell = (/ Lx , Ly /)
      S = 0.d0
      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        do j = 1 , size(AO2%exponent)
          beta = AO2%exponent(j)
          call clifford_pair(2, alpha, beta, r1, r2, pe, cen)
          do d = 1 , 2
            Aarg = 2.d0*pe(d)/(av(d)*av(d))
            kd   = dexp(-2.d0*(alpha + beta - pe(d))/(av(d)*av(d)))
            call clifford_bessel_ratios(min(NMAX_CLIFF, l1(d)+l2(d)), Aarg, R)
            sd(d) = kd * overlap_1D_clifford(l1(d), l2(d),                                   &
                         dsin(av(d)*(cen(d)-r1(d))), dcos(av(d)*(cen(d)-r1(d))),           &
                         dsin(av(d)*(cen(d)-r2(d))), dcos(av(d)*(cen(d)-r2(d))),           &
                         1.d0/av(d), R, Lcell(d))
          end do
          esum = alpha + beta ; mu = alpha*beta/esum
          zP   = (alpha*r1(3) + beta*r2(3))/esum
          sd(3) = dexp(-mu*(r1(3)-r2(3))**2) * overlap_1D_real(l1(3), l2(3), esum, zP-r1(3), zP-r2(3))
          S = S + AO1%coefficient(i)*AO2%coefficient(j) * sd(1)*sd(2)*sd(3)
        end do
      end do

end subroutine overlap_integral_toroidal_2D

subroutine overlap_matrix_toroidal_2D(nBas,AO,OV)

      use files
      use classification_ERI
      use keywords
      use unitcell_module
      use clifford_tools, only : init_clifford, orbital_to_lxyz

      implicit none

      integer           ,intent(in)  :: nBas
      type(ERI_function),intent(in)  :: AO(nBas)
      double precision  ,intent(out) :: OV(nBas,nBas)

      integer          :: i , j , lxyz(3,nBas)
      double precision :: r1(3) , r2(3)

      call init_clifford()
      do i = 1 , nBas
        call orbital_to_lxyz(AO(i)%orbital, lxyz(:,i))
      end do

      OV = 0.d0
      !$omp parallel do default(shared) private(i,j,r1,r2) schedule(dynamic,4)
      do i = 1 , nfuc
        r1 = (/ AO(i)%x , AO(i)%y , AO(i)%z /)
        do j = 1 , nBas
          r2 = (/ AO(j)%x , AO(j)%y , AO(j)%z /)
          call overlap_integral_toroidal_2D(lxyz(:,i),lxyz(:,j),r1,r2,AO(i),AO(j),OV(i,j))
        end do
      end do
      !$omp end parallel do

      do i = 1 , nfuc
        do j = 1 , nBas
          if (abs(OV(i,j)) < 1e-14) OV(i,j) = 0.d0
        end do
      end do
      call fill_with_symmetry(nBas, OV)

      if (c_OV) then
        write(outfile,'(a)') ""
        write(outfile,'(a)') "------------------------------------------"
        write(outfile,'(a)') "The overlap matrix (2D)   "
        write(outfile,'(a)') "------------------------------------------"
        do i = 1 , nBas - 1
          do j = i , nBas
            write(outfile,'(i3,i3, 6x,1000(f16.12,2x))') i , j, OV(i,j)
          end do
        end do
      end if

end subroutine overlap_matrix_toroidal_2D
