!=====================================================================!
!                                                                     !
!          2D toroidal : NUCLEAR ATTRACTION INTEGRALS                !
!                                                                     !
!   x,y : Clifford Gaussians (periodic),  z : real Gaussian           !
!   1/|r-C|_E = 2/sqrt(pi) Int_0^inf dt prod_d exp(-t^2 |x_d-C_d|^2)  !
!   periodic d (once per primitive pair, t-independent):              !
!      G^d_n = sum_k d_k Ĩ_{|n+k|}(A_d)                               !
!   per nucleus  F^d_n = Re[ e^{i n a_d (P_d - C_d)} G^d_n ]          !
!   per t        I_d(t) = F^d_0 Ĩ_0(C_d) + 2 sum_n F^d_n Ĩ_n(C_d)      !
!                (one dot product with the column t_bessel(:,t,d))    !
!   real d : Gaussian product with exp(-t^2 (z-Cz)^2)                 !
!                                                                     !
!   nuclear_attraction_integral_toroidal_2D  one pair, all nuclei     !
!   nuclear_attraction_matrix_toroidal_2D    full matrix (OpenMP)     !
!=====================================================================!

subroutine nuclear_attraction_integral_toroidal_2D(l1,l2,n_atoms,geometry,atoms,r1,r2,AO1,AO2,V)

      use torus_init
      use classification_ERI
      use atom_basis
      use clifford_tools

      implicit none

      integer           ,intent(in)  :: l1(3) , l2(3) , n_atoms
      double precision  ,intent(in)  :: geometry(n_atoms,3)
      type(atom)        ,intent(in)  :: atoms(n_atoms)
      double precision  ,intent(in)  :: r1(3) , r2(3)
      type(ERI_function),intent(in)  :: AO1 , AO2
      double precision  ,intent(out) :: V

      double precision, parameter :: pi = 3.14159265358979323846d0
      integer          :: i , j , d , n , k , c , kk , M(2) , Np(2) , nw , nt , Mmax
      double precision :: alpha , beta , pe(2) , cen(2) , av(2) , Lcell(2) , Aarg , pref
      double precision :: phi , c1 , s1 , cn , sn , tmp , res , prod , sd
      double precision :: esum , mu , zP , zpa , zpb , zpc , t2 , ppt , eta
      complex(8)       :: acc
      complex(8)      , allocatable :: dco(:) , G(:,:)
      double precision, allocatable :: Iv(:) , F(:,:)

      av    = (/ ax , ay /)
      Lcell = (/ Lx , Ly /)
      M     = l1(1:2) + l2(1:2)
      Mmax  = maxval(M)
      nw    = ubound(t_bessel,1) + Mmax + 2
      allocate(dco(-Mmax:Mmax), Iv(0:nw), G(0:nw,2), F(0:nw,2))

      V = 0.d0
      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        do j = 1 , size(AO2%exponent)
          beta = AO2%exponent(j)
          call clifford_pair(2, alpha, beta, r1, r2, pe, cen)
          esum = alpha + beta ; mu = alpha*beta/esum
          zP   = (alpha*r1(3) + beta*r2(3))/esum
          zpa  = zP - r1(3) ; zpb = zP - r2(3)

          pref = AO1%coefficient(i)*AO2%coefficient(j) * 2.d0/dsqrt(pi) * dexp(-mu*(r1(3)-r2(3))**2)
          do d = 1 , 2
            Aarg  = 2.d0*pe(d)/(av(d)*av(d))
            pref  = pref * dexp(-2.d0*(alpha + beta - pe(d))/(av(d)*av(d))) * Lcell(d) / av(d)**M(d)
            Np(d) = min(bessel_nmax_estimate(Aarg), nw - Mmax - 2)
            call fourier_sin_product(l1(d), l2(d), av(d)*(cen(d)-r1(d)), av(d)*(cen(d)-r2(d)), dco(-M(d):M(d)))
            call bessel_scaled_vector_fast(Np(d)+M(d), Aarg, Iv(0:Np(d)+M(d)))
            do n = 0 , Np(d)
              acc = (0.d0,0.d0)
              do kk = -M(d) , M(d)
                acc = acc + dco(kk) * Iv(abs(n+kk))
              end do
              G(n,d) = acc
            end do
          end do

          do c = 1 , n_atoms
            do d = 1 , 2
              phi = av(d) * (cen(d) - geometry(c,d))
              c1 = dcos(phi) ; s1 = dsin(phi)
              cn = 1.d0      ; sn = 0.d0
              do n = 0 , Np(d)
                if (n > 0) then
                  if (mod(n,32) == 0) then
                    cn = dcos(n*phi) ; sn = dsin(n*phi)
                  else
                    tmp = cn*c1 - sn*s1
                    sn  = sn*c1 + cn*s1
                    cn  = tmp
                  end if
                end if
                F(n,d) = dble(G(n,d))*cn - aimag(G(n,d))*sn
              end do
            end do
            zpc = zP - geometry(c,3)

            res = 0.d0
            do k = 1 , nq_t
              prod = 1.d0
              do d = 1 , 2
                nt = min(Np(d), t_nmax(k,d))
                sd = F(0,d) * t_bessel(0,k,d)
                if (nt >= 1) sd = sd + 2.d0 * dot_product(F(1:nt,d), t_bessel(1:nt,k,d))
                prod = prod * sd
              end do
              t2  = t_node(k)*t_node(k)
              ppt = esum + t2
              eta = t2/ppt
              prod = prod * dexp(-esum*eta*zpc*zpc) * overlap_1D_real(l1(3), l2(3), ppt, zpa - eta*zpc, zpb - eta*zpc)
              res = res + t_weight(k) * prod
            end do

            V = V - dble(atoms(c)%charge) * pref * res
          end do

        end do
      end do

      deallocate(dco, Iv, G, F)

end subroutine nuclear_attraction_integral_toroidal_2D

subroutine nuclear_attraction_matrix_toroidal_2D(number_of_atoms,nBas,geometry,atoms,AO,NA)

      use files
      use atom_basis
      use torus_init
      use classification_ERI
      use keywords
      use unitcell_module
      use clifford_tools, only : init_clifford, orbital_to_lxyz, setup_t_tables, bessel_nmax_estimate

      implicit none

      integer           ,intent(in)  :: number_of_atoms , nBas
      double precision  ,intent(in)  :: geometry(number_of_atoms,3)
      type(atom)        ,intent(in)  :: atoms(number_of_atoms)
      type(ERI_function),intent(in)  :: AO(nBas)
      double precision  ,intent(out) :: NA(nBas,nBas)

      integer          :: i , j , lxyz(3,nBas) , lmax , ncap
      double precision :: r1(3) , r2(3) , emax , amin

      call init_clifford()
      lmax = 0 ; emax = 0.d0
      do i = 1 , nBas
        call orbital_to_lxyz(AO(i)%orbital, lxyz(:,i))
        lmax = max(lmax, maxval(lxyz(1:2,i)))
        emax = max(emax, maxval(AO(i)%exponent))
      end do
      amin = min(ax, ay)
      ncap = bessel_nmax_estimate(4.d0*emax/(amin*amin)) + 2*lmax + 1
      call setup_t_tables(ncap, 2)

      NA = 0.d0
      !$omp parallel do default(shared) private(i,j,r1,r2) schedule(dynamic,4) collapse(2)
      do i = 1 , nfuc
        do j = 1 , nBas
          r1 = (/ AO(i)%x , AO(i)%y , AO(i)%z /)
          r2 = (/ AO(j)%x , AO(j)%y , AO(j)%z /)
          call nuclear_attraction_integral_toroidal_2D(lxyz(:,i),lxyz(:,j),number_of_atoms,geometry,atoms, &
                                                           r1,r2,AO(i),AO(j),NA(i,j))
        end do
      end do
      !$omp end parallel do

      do i = 1 , nfuc
        do j = 1 , nBas
          if (abs(NA(i,j)) < 1e-14) NA(i,j) = 0.d0
        end do
      end do
      call fill_with_symmetry(nBas, NA)

      if (c_NA) then
        write(outfile,'(a)') ""
        write(outfile,'(a)') "------------------------------------------------------"
        write(outfile,'(a)') "The nuclear attraction matrix (2D)    "
        write(outfile,'(a)') "------------------------------------------------------"
        do i = 1 , nBas - 1
          do j = i , nBas
            write(outfile,'(i3,i3, 6x,1000(f16.12,2x))') i , j, NA(i,j)
          end do
        end do
      end if

end subroutine nuclear_attraction_matrix_toroidal_2D
