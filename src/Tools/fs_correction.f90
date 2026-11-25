!=====================================================================!
!                                                                     !
!   [Nov 2026] keyword fs_correction (Tori3D, SAO RHF)                !
!                                                                     !
!   Leading (1/K^2) finite-size correction of the electrostatic       !
!   energy on the Clifford torus (doc/madelung_clifford_c2.pdf):      !
!                                                                     !
!     dE_ES = 1/(2 Omega) sum_{G != 0} |rho(G)|^2 dv(G)   per cell    !
!     dv(G) = pi sum_i a_i^2 [8 G_i^4/G^8 - 8 G_i^2/G^6 + 1/G^4]      !
!     a_i = 2 pi / L_i (L_i = supercell length), G = reciprocal       !
!     vectors of the UNIT cell (G_i = 2 pi n_i / c_i)                 !
!                                                                     !
!   rho(G) = rho_nuc(G) - rho_el(G), per unit cell:                   !
!     rho_nuc(G) = sum_{j in cell 0} Z_j exp(-i G.R_j)                !
!     rho_el(G)  = sum_{mu in cell 0, nu} P(mu,nu) (mu nu)(G)         !
!     (mu nu)(G) = int over the torus  mu(r) nu(r) exp(-i G.r)        !
!   Clifford pair product along d (p = barycentric exponent, P its    !
!   centre, A = 2p/a^2, u = a(x - P)):                                !
!     mu nu = a^-M e^{-2(al+be-p)/a^2} sum_k d_k e^{iku}              !
!             * e^{A} sum_n I~_n(A) e^{inu}                           !
!     => FT at g = m a :  L e^{-2(al+be-p)/a^2} a^-M c_m e^{-i m a P}  !
!        c_m = sum_k d_k I~_{|m-k|}(A)                                !
!   (same pair quantities as the nuclear attraction; d_k from         !
!   fourier_sin_product). A unit-cell G has m_d = n_d * K_d.          !
!                                                                     !
!   The three parts are printed separately:                           !
!     nuc-nuc   : point nuclei (lattice sum, converges ~1/nmax:       !
!                 several nmax, extrapolated)                         !
!     nuc-el    : cross term,  el-el : electrons (converge fast)      !
!   Check: rho_el(0) = number of electrons per unit cell.             !
!=====================================================================!
subroutine fs_correction_3D(nBas, AO, n_atoms, geometry, charge, P, EHF)

      use files,              only : outfile
      use classification_ERI, only : ERI_function
      use unitcell_module,    only : nfuc, nx, ny, nz
      use torus_init,         only : Lx, Ly, Lz, ax, ay, az
      use clifford_tools,     only : init_clifford, orbital_to_lxyz, bessel_scaled_vector_fast, &
                                     fourier_sin_product, clifford_pair
      use omp_lib

      implicit none

      integer           , intent(in) :: nBas, n_atoms
      type(ERI_function), intent(in) :: AO(nBas)
      double precision  , intent(in) :: geometry(n_atoms,3)
      integer           , intent(in) :: charge(n_atoms)
      double precision  , intent(in) :: P(nBas,nBas)
      double precision  , intent(in) :: EHF      ! total HF energy of the supercell (electronic + nuclear)

      integer, parameter :: NG_EL = 10                 ! electrons: |n_i| <= NG_EL
      integer, parameter :: NGN(5) = (/ 30, 40, 50, 60, 70 /)   ! nuclei: several cutoffs, extrapolated
      type(ERI_function), allocatable :: AOn(:)
      integer :: nbes, Kc(3), lq(3,nBas), mu, nu, ip, jp, d, k, n, m, Mm, nmx, i1, i2, i3, nat_cell, j, nG
      double precision :: Lv(3), av(3), cell(3), Omega, pe(3), cen(3), r1(3), r2(3), pref, ex1, ex2, A
      double precision :: G(3), G2, dv, e_nn, e_ne, e_ee, nel0, wmax, w0
      double precision :: e_nn_cut(5), xf(3), coef(3), r_nn, r_ne, r_ee
      double precision :: Rg(0:NG_EL,0:NG_EL,0:NG_EL)
      complex(8), allocatable :: rho_el(:,:,:), rho_nuc(:,:,:), cm(:,:), dco(:)
      double precision, allocatable :: Iv(:)
      complex(8) :: zn, ze, ph

      w0 = omp_get_wtime()
      Kc = (/ max(nx,1), max(ny,1), max(nz,1) /)
      Lv = (/ Lx, Ly, Lz /) ; av = (/ ax, ay, az /)
      cell = Lv / dble(Kc)
      Omega = cell(1)*cell(2)*cell(3)
      nat_cell = n_atoms / (Kc(1)*Kc(2)*Kc(3))

      ! normalized AO list, as in Tori3D
      allocate(AOn(nBas))
      AOn = AO
      call normalize_AO_toroidal_3D(nBas, AOn, .false.)
      call init_clifford()
      do mu = 1, nBas
        call orbital_to_lxyz(AOn(mu)%orbital, lq(:,mu))
      end do

      ! ---------------- electrons: rho_el(n1,n2,n3) ---------------- !
      allocate(rho_el(-NG_EL:NG_EL,-NG_EL:NG_EL,-NG_EL:NG_EL))
      rho_el = (0.d0,0.d0)
      nmx = NG_EL*maxval(Kc) + 12
      !$omp parallel default(shared) private(nbes,mu,nu,ip,jp,d,k,n,m,Mm,pe,cen,r1,r2,pref,ex1,ex2,A,cm,dco,Iv,i1,i2,i3,ph,coef)
      allocate(cm(-NG_EL:NG_EL,3), dco(-8:8), Iv(0:4000))
      !$omp do schedule(dynamic,8) reduction(+:rho_el)
      do mu = 1, nfuc
        do nu = 1, nBas
          if (abs(P(mu,nu)) < 1.d-12) cycle
          r1 = (/ AOn(mu)%x, AOn(mu)%y, AOn(mu)%z /)
          r2 = (/ AOn(nu)%x, AOn(nu)%y, AOn(nu)%z /)
          do ip = 1, size(AOn(mu)%exponent)
            ex1 = AOn(mu)%exponent(ip)
            do jp = 1, size(AOn(nu)%exponent)
              ex2 = AOn(nu)%exponent(jp)
              call clifford_pair(3, ex1, ex2, r1, r2, pe, cen)
              pref = P(mu,nu)*AOn(mu)%coefficient(ip)*AOn(nu)%coefficient(jp)
              do d = 1, 3
                Mm = lq(d,mu) + lq(d,nu)
                coef(d) = Lv(d)*dexp(-2.d0*(ex1+ex2-pe(d))/av(d)**2)/av(d)**Mm
                A = 2.d0*pe(d)/av(d)**2
                call fourier_sin_product(lq(d,mu), lq(d,nu), av(d)*(cen(d)-r1(d)), av(d)*(cen(d)-r2(d)), dco(-Mm:Mm))
                nbes = min(4000, max(nmx+Mm, Mm + 40 + int(10.d0*dsqrt(A))))   ! [Nov 2026] backward recurrence start
                call bessel_scaled_vector_fast(nbes, A, Iv(0:nbes))
                do n = -NG_EL, NG_EL
                  m = n*Kc(d)
                  cm(n,d) = (0.d0,0.d0)
                  do k = -Mm, Mm
                    cm(n,d) = cm(n,d) + dco(k)*Iv(abs(m-k))
                  end do
                  cm(n,d) = coef(d)*cm(n,d)*exp(dcmplx(0.d0, -dble(m)*av(d)*cen(d)))
                end do
              end do
              if (pref*coef(1)*coef(2)*coef(3) == 0.d0) cycle
              do i3 = -NG_EL, NG_EL
                do i2 = -NG_EL, NG_EL
                  ph = pref*cm(i2,2)*cm(i3,3)
                  do i1 = -NG_EL, NG_EL
                    rho_el(i1,i2,i3) = rho_el(i1,i2,i3) + ph*cm(i1,1)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      !$omp end do
      deallocate(cm, dco, Iv)
      !$omp end parallel
      nel0 = dble(rho_el(0,0,0))

      ! ---------------- nuclei ---------------- !
      ! cross and el-el terms with |n| <= NG_EL ; nuc-nuc with large cutoffs (extrapolated)
      e_ne = 0.d0 ; e_ee = 0.d0
      do i3 = -NG_EL, NG_EL
        do i2 = -NG_EL, NG_EL
          do i1 = -NG_EL, NG_EL
            if (i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
            G = 2.d0*acos(-1.d0)*(/ dble(i1), dble(i2), dble(i3) /)/cell
            dv = dvG(G)
            zn = (0.d0,0.d0)
            do j = 1, nat_cell
              zn = zn + dble(charge(j))*exp(dcmplx(0.d0, -sum(G*geometry(j,:))))
            end do
            ze = rho_el(i1,i2,i3)
            e_ne = e_ne + dv*(-2.d0*dble(zn*conjg(ze)))
            e_ee = e_ee + dv*dble(ze*conjg(ze))
          end do
        end do
      end do
      do k = 1, 5
        nG = NGN(k)
        e_nn_cut(k) = 0.d0
        !$omp parallel do default(shared) private(i1,i2,i3,G,zn,j) reduction(+:e_nn_cut) collapse(2)
        do i3 = -nG, nG
          do i2 = -nG, nG
            do i1 = -nG, nG
              if (i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
              G = 2.d0*acos(-1.d0)*(/ dble(i1), dble(i2), dble(i3) /)/cell
              zn = (0.d0,0.d0)
              do j = 1, nat_cell
                zn = zn + dble(charge(j))*exp(dcmplx(0.d0, -sum(G*geometry(j,:))))
              end do
              e_nn_cut(k) = e_nn_cut(k) + dvG(G)*dble(zn*conjg(zn))
            end do
          end do
        end do
        !$omp end parallel do
      end do
      e_nn = extrap(e_nn_cut)
      e_nn = e_nn/(2.d0*Omega) ; e_ne = e_ne/(2.d0*Omega) ; e_ee = e_ee/(2.d0*Omega)

      ! [Nov 2026] higher orders of the kernel: R(G) = dv_exact(G) - dv_a2(G), |n_i| <= NG_EL
      ! (tail |n_i| > NG_EL of R ~ a^4/G^6 neglected)
      call es_remainder(Rg)
      r_nn = 0.d0 ; r_ne = 0.d0 ; r_ee = 0.d0
      do i3 = -NG_EL, NG_EL
        do i2 = -NG_EL, NG_EL
          do i1 = -NG_EL, NG_EL
            if (i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
            G = 2.d0*acos(-1.d0)*(/ dble(i1), dble(i2), dble(i3) /)/cell
            zn = (0.d0,0.d0)
            do j = 1, nat_cell
              zn = zn + dble(charge(j))*exp(dcmplx(0.d0, -sum(G*geometry(j,:))))
            end do
            ze = rho_el(i1,i2,i3)
            dv = Rg(abs(i1),abs(i2),abs(i3))
            r_nn = r_nn + dv*dble(zn*conjg(zn))
            r_ne = r_ne + dv*(-2.d0*dble(zn*conjg(ze)))
            r_ee = r_ee + dv*dble(ze*conjg(ze))
          end do
        end do
      end do
      r_nn = r_nn/(2.d0*Omega) ; r_ne = r_ne/(2.d0*Omega) ; r_ee = r_ee/(2.d0*Omega)

      write(outfile,'(a)') ''
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      write(outfile,'(a)') '  Finite-size correction (keyword fs_correction), electrostatics'
      write(outfile,'(a)') '  dE_ES = 1/(2 Omega) sum_G |rho(G)|^2 dv(G), leading order 1/K^2'
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      write(outfile,'(a,3i4)')        '  cells per direction K           : ', Kc
      write(outfile,'(a,f14.8)')      '  electrons per cell (rho_el(0))  : ', nel0
      write(outfile,'(a,i4)')         '  electrons: G cutoff |n_i| <=    : ', NG_EL
      write(outfile,'(a,5f14.8)')     '  nuc-nuc, cutoffs 30..70 (x 1/2Omega):', e_nn_cut/(2.d0*Omega)
      write(outfile,'(a,f18.10,a)')   '  dE nuc-nuc (extrapolated)       : ', e_nn, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE nuc-el                       : ', e_ne, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE el-el                        : ', e_ee, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_ES total                     : ', e_nn+e_ne+e_ee, ' Eh/cell'
      write(outfile,'(a)')            '  higher orders of the kernel (exact - a^2 term, |n_i| <= 10):'
      write(outfile,'(a,f18.10,a)')   '    nuc-nuc                       : ', r_nn, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '    nuc-el                        : ', r_ne, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '    el-el                         : ', r_ee, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_ES all orders                : ', e_nn+e_ne+e_ee+r_nn+r_ne+r_ee, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  E_HF per cell                   : ', EHF/dble(Kc(1)*Kc(2)*Kc(3)), ' Eh'
      write(outfile,'(a,f18.10,a)')   '  E_HF per cell - dE_ES (a^2)     : ', EHF/dble(Kc(1)*Kc(2)*Kc(3)) - (e_nn+e_ne+e_ee), ' Eh'
      write(outfile,'(a,f10.2,a)')    '  wall time                       : ', omp_get_wtime()-w0, ' s'
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      flush(outfile)
      deallocate(rho_el)

      ! [Nov 2026] exchange part (localized orbitals, Gaussian model)
      call fs_exchange_3D(nBas, AOn, lq, P, Lv, av, Kc, EHF, e_nn+e_ne+e_ee, e_nn+e_ne+e_ee+r_nn+r_ne+r_ee)
      deallocate(AOn)

contains

      ! [Nov 2026] R(G) = dv_exact(G) - dv_a2(G) for G = 2 pi n / cell, 0 <= n_i <= NG_EL
      !   dv_exact = (2/sqrt pi) int dt [ prod_d L_d I~_{m_d}(z_d) - P0 ],  P0 = pi^3/2 t^-3 exp(-G^2/4t^2)
      !   dv_a2    = (2/sqrt pi) int dt  P0 sum_d eps_d                     (analytic: dvG)
      !   eps_d    = a_d^2/(48 t^2) (G_d^4/(4 t^4) - 3 G_d^2/t^2 + 3),  m_d = n_d K_d,  z_d = 2 t^2/a_d^2
      ! t = exp(s), trapezoid in s.
      subroutine es_remainder(R)
        double precision, intent(out) :: R(0:NG_EL,0:NG_EL,0:NG_EL)
        double precision, allocatable :: Ivt(:), fac(:,:)
        double precision :: gmin, gmax, smin, smax, ds, t, wt, z, P0, eps, gv(3), g2, pi
        integer :: ns, is, dd, nb, j1, j2, j3
        pi = acos(-1.d0)
        gmin = 2.d0*pi/maxval(cell)
        gmax = 2.d0*pi*dble(NG_EL)*dsqrt(sum(1.d0/cell**2))
        smin = dlog(gmin/40.d0) ; smax = dlog(20.d0*gmax)
        ds = 0.01d0 ; ns = int((smax-smin)/ds) + 1
        R = 0.d0
        !$omp parallel default(shared) private(is,t,wt,dd,z,nb,Ivt,fac,j1,j2,j3,gv,g2,P0,eps)
        allocate(fac(0:NG_EL,3))
        !$omp do schedule(dynamic,4) reduction(+:R)
        do is = 0, ns
          t = dexp(smin + dble(is)*ds)
          wt = t*ds ; if (is == 0 .or. is == ns) wt = 0.5d0*wt
          do dd = 1, 3
            z = 2.d0*t*t/av(dd)**2
            nb = NG_EL*Kc(dd) + 40 + int(10.d0*dsqrt(z))
            allocate(Ivt(0:nb))
            call bessel_scaled_vector_fast(nb, z, Ivt)
            do j1 = 0, NG_EL
              fac(j1,dd) = Lv(dd)*Ivt(j1*Kc(dd))
            end do
            deallocate(Ivt)
          end do
          do j3 = 0, NG_EL
            do j2 = 0, NG_EL
              do j1 = 0, NG_EL
                if (j1 == 0 .and. j2 == 0 .and. j3 == 0) cycle
                gv = 2.d0*pi*(/ dble(j1), dble(j2), dble(j3) /)/cell
                g2 = sum(gv**2)
                P0 = pi**1.5d0/t**3*dexp(-g2/(4.d0*t*t))
                eps = sum(av**2/(48.d0*t*t)*(gv**4/(4.d0*t**4) - 3.d0*gv**2/(t*t) + 3.d0))
                R(j1,j2,j3) = R(j1,j2,j3) + wt*(fac(j1,1)*fac(j2,2)*fac(j3,3) - P0*(1.d0+eps))
              end do
            end do
          end do
        end do
        !$omp end do
        deallocate(fac)
        !$omp end parallel
        R = 2.d0/dsqrt(pi)*R
      end subroutine es_remainder

      double precision function dvG(Gv)
        double precision, intent(in) :: Gv(3)
        double precision :: g2
        integer :: i
        g2 = sum(Gv**2)
        dvG = 0.d0
        do i = 1, 3
          dvG = dvG + av(i)**2*(8.d0*Gv(i)**4/g2**4 - 8.d0*Gv(i)**2/g2**3 + 1.d0/g2**2)
        end do
        dvG = acos(-1.d0)*dvG
      end function dvG

      double precision function extrap(v)      ! fit c0 + c1/n + c2/n^2 over the 5 cutoffs
        double precision, intent(in) :: v(5)
        double precision :: Am(5,3), b(5), At(3,3), rhs(3), x(3), f
        integer :: r, c, i, jj
        do r = 1, 5
          Am(r,:) = (/ 1.d0, 1.d0/dble(NGN(r)), 1.d0/dble(NGN(r))**2 /)
          b(r) = v(r)
        end do
        At = matmul(transpose(Am), Am) ; rhs = matmul(transpose(Am), b)
        do i = 1, 3                                    ! Gauss elimination 3x3
          do jj = i+1, 3
            f = At(jj,i)/At(i,i)
            At(jj,:) = At(jj,:) - f*At(i,:) ; rhs(jj) = rhs(jj) - f*rhs(i)
          end do
        end do
        do i = 3, 1, -1
          x(i) = (rhs(i) - sum(At(i,i+1:3)*x(i+1:3)))/At(i,i)
        end do
        extrap = x(1)
      end function extrap

end subroutine fs_correction_3D


!=====================================================================!
!                                                                     !
!   [Nov 2026] keyword fs_correction: EXCHANGE part (Tori3D, RHF)     !
!                                                                     !
!   On the Clifford torus the Coulomb kernel is 1/d, with             !
!     d^2 = sum_i (2/a_i^2)(1 - cos a_i x_i)  ~  r^2 - sum a_i^2 x_i^4/12 !
!   so for |r| << L                                                    !
!     1/d - 1/r = sum_i a_i^2 x_i^4 / (24 r^3)          (= dv(r))      !
!   There is NO 1/L (Madelung) term in the exchange; the leading      !
!   error is 1/L^2, set by the size of the exchange hole.             !
!                                                                     !
!   Model: E_X = - sum_{ij} (ij|ji). With well localized orbitals the !
!   self terms dominate; the correction kept here is                  !
!     dE_X = - sum_i  int int rho_i(1) rho_i(2) dv(r_12)              !
!   with rho_i = |phi_i|^2 a Gaussian of variances sigma_{i,d}^2      !
!   (r_12 Gaussian with variances 2 sigma_{i,d}^2):                    !
!     dE_X = - sum_i sum_d (a_d^2/24) < x_d^4 / r^3 >                  !
!     < x_d^4/r^3 > = (4/sqrt(pi)) int_0^inf t^2 3 s_d^4               !
!          (1+2t^2 s_d^2)^(-5/2) prod_{k/=d} (1+2t^2 s_k^2)^(-1/2) dt  !
!          s_k^2 = 2 sigma_{i,k}^2                                    !
!   Isotropic check: - sum_i a^2 sigma_i / (10 sqrt(pi))               !
!   (verified against the exact torus Fourier sum of 1/d, ~1%).       !
!   Neglected: inter-orbital terms (ij|ji), i /= j, and the deviation !
!   of rho_i from a Gaussian.                                          !
!                                                                     !
!   Orbitals: occupied space from the final SCF density P             !
!     (Loewdin: eigenvectors of S^1/2 (P/2) S^1/2 with eigenvalue 1,  !
!      C = S^-1/2 U), then localized by Jacobi rotations maximizing  !
!     Omega = sum_d sum_i |z_ii^d|^2 ,  z^d = C^T Z^d C,              !
!     Z^d_{mu nu} = <mu| exp(-i a_d x_d) |nu>   (Resta / Silvestrelli) !
!   Spreads: sigma_{i,d}^2 = - ln |z_ii^d|^2 / a_d^2                   !
!   S and Z^d from the same Clifford pair Fourier transform as rho_el !
!   (m = 0 and m = 1).                                                !
!=====================================================================!
subroutine fs_exchange_3D(nBas, AOn, lq, P, Lv, av, Kc, EHF, dE_ES, dE_ES_all)

      use files,              only : outfile
      use classification_ERI, only : ERI_function
      use clifford_tools,     only : bessel_scaled_vector_fast, fourier_sin_product, clifford_pair
      use omp_lib

      implicit none

      integer           , intent(in) :: nBas, lq(3,nBas), Kc(3)
      type(ERI_function), intent(in) :: AOn(nBas)
      double precision  , intent(in) :: P(nBas,nBas), Lv(3), av(3), EHF, dE_ES, dE_ES_all

      integer, parameter :: MAXSWEEP = 500
      double precision, allocatable :: S(:,:), Zi(:,:,:), ev(:), Sh(:,:), Smh(:,:), &
                                       C(:,:), T1(:,:), T2(:,:), sig2(:,:), work(:), Rx(:)
      complex(8), allocatable :: z(:,:,:), zci(:), zcj(:)
      complex(8) :: u, vv, w
      integer :: d, k, nocc, lwork, info, i, j, isweep, ncls, ic
      double precision :: ntot, w0, &
                          B, Cc, th, cs, sn, gain, gsum, Om0, Om1, dEx, s2(3), dEi, sgm, ocmin, ocmax, dExa, ncell
      double precision :: cls_s(200), cls_e(200)
      integer :: cls_n(200)
      double precision, parameter :: pi = 3.14159265358979323846d0

      w0 = omp_get_wtime()

      ! ---------- overlap S over the whole supercell ---------- !
      allocate(S(nBas,nBas))
      call pair_matrix(0, S)

      ! ---------- occupied orbitals from P (Loewdin) ---------- !
      ! memory: at most 4 nBas^2 at a time
      ntot = sum(P*S)                                  ! number of electrons in the supercell
      nocc = nint(0.5d0*ntot)
      allocate(ev(nBas))
      lwork = -1 ; allocate(work(1))
      call dsyev('V', 'U', nBas, S, nBas, ev, work, lwork, info)       ! S <- eigenvectors V
      lwork = int(work(1)) ; deallocate(work) ; allocate(work(lwork))
      call dsyev('V', 'U', nBas, S, nBas, ev, work, lwork, info)
      allocate(T1(nBas,nBas), Sh(nBas,nBas), Smh(nBas,nBas))
      do i = 1, nBas
        T1(:,i) = S(:,i)*dsqrt(ev(i))
      end do
      call dgemm('N', 'T', nBas, nBas, nBas, 1.d0, T1, nBas, S, nBas, 0.d0, Sh, nBas)     ! S^1/2
      do i = 1, nBas
        T1(:,i) = S(:,i)/dsqrt(ev(i))
      end do
      call dgemm('N', 'T', nBas, nBas, nBas, 1.d0, T1, nBas, S, nBas, 0.d0, Smh, nBas)    ! S^-1/2
      call dgemm('N', 'N', nBas, nBas, nBas, 0.5d0, P, nBas, Sh, nBas, 0.d0, T1, nBas)     ! (P/2) S^1/2
      call dgemm('N', 'N', nBas, nBas, nBas, 1.d0, Sh, nBas, T1, nBas, 0.d0, S, nBas)      ! S <- M
      deallocate(Sh, T1)
      deallocate(work) ; lwork = -1 ; allocate(work(1))
      call dsyev('V', 'U', nBas, S, nBas, ev, work, lwork, info)
      lwork = int(work(1)) ; deallocate(work) ; allocate(work(lwork))
      call dsyev('V', 'U', nBas, S, nBas, ev, work, lwork, info)
      ocmin = ev(nBas-nocc+1)                          ! smallest occupied eigenvalue (should be 1)
      ocmax = 0.d0 ; if (nBas > nocc) ocmax = ev(nBas-nocc)   ! largest virtual one (should be 0)
      allocate(C(nBas,nocc))
      call dgemm('N', 'N', nBas, nocc, nBas, 1.d0, Smh, nBas, S(:,nBas-nocc+1:nBas), nBas, 0.d0, C, nBas)
      deallocate(work, ev, Smh)

      ! ---------- z^d = C^T Z^d C, one direction at a time ---------- !
      allocate(z(nocc,nocc,3), T1(nBas,nocc), T2(nocc,nocc), Zi(nBas,nBas,1))
      do d = 1, 3
        call pair_matrix(d, S, Zi(:,:,1))              ! S <- Re Z^d , Zi <- Im Z^d
        call dgemm('N', 'N', nBas, nocc, nBas, 1.d0, S, nBas, C, nBas, 0.d0, T1, nBas)
        call dgemm('T', 'N', nocc, nocc, nBas, 1.d0, C, nBas, T1, nBas, 0.d0, T2, nocc)
        z(:,:,d) = dcmplx(T2, 0.d0)
        call dgemm('N', 'N', nBas, nocc, nBas, 1.d0, Zi(:,:,1), nBas, C, nBas, 0.d0, T1, nBas)
        call dgemm('T', 'N', nocc, nocc, nBas, 1.d0, C, nBas, T1, nBas, 0.d0, T2, nocc)
        z(:,:,d) = z(:,:,d) + dcmplx(0.d0, T2)
      end do
      deallocate(S, Zi, T1, T2, C)

      ! ---------- Jacobi localization: maximize sum_d sum_i |z_ii^d|^2 ---------- !
      Om0 = omega_val()
      allocate(zci(nocc), zcj(nocc))
      do isweep = 1, MAXSWEEP
        gsum = 0.d0
        do i = 1, nocc-1
          do j = i+1, nocc
            ! Omega(theta) = const + B cos(4 theta) + Cc sin(4 theta)
            B = 0.d0 ; Cc = 0.d0
            do d = 1, 3
              vv = 0.5d0*(z(i,i,d) - z(j,j,d)) ; w = z(i,j,d)
              B  = B  + (abs(vv)**2 - abs(w)**2)
              Cc = Cc + 2.d0*dble(vv*conjg(w))
            end do
            gain = dsqrt(B*B + Cc*Cc) - B
            if (gain < 1.d-14) cycle
            gsum = gsum + gain
            th = 0.25d0*atan2(Cc, B)
            cs = cos(th) ; sn = sin(th)
            do d = 1, 3                                 ! z <- R^T z R   (phi_i' = c phi_i + s phi_j)
              zci = z(:,i,d) ; zcj = z(:,j,d)
              z(:,i,d) =  cs*zci + sn*zcj
              z(:,j,d) = -sn*zci + cs*zcj
              zci = z(i,:,d) ; zcj = z(j,:,d)
              z(i,:,d) =  cs*zci + sn*zcj
              z(j,:,d) = -sn*zci + cs*zcj
            end do
          end do
        end do
        if (gsum < 1.d-11*dble(nocc)) exit
      end do
      Om1 = omega_val()

      ! ---------- spreads and dE_X ---------- !
      allocate(sig2(3,nocc))
      dEx = 0.d0 ; ncls = 0
      do i = 1, nocc
        do d = 1, 3
          sig2(d,i) = -dlog(abs(z(i,i,d))**2)/av(d)**2
        end do
        s2 = 2.d0*sig2(:,i)
        dEi = 0.d0
        do d = 1, 3
          dEi = dEi + av(d)**2/24.d0*x4_over_r3(s2, d)
        end do
        dEx = dEx - dEi
        ! classes of orbitals (same mean spread and same dE_i to 1e-4)
        sgm = dsqrt(sum(sig2(:,i))/3.d0)
        ic = 0
        do k = 1, ncls
          if (abs(cls_s(k)-sgm) < 1.d-3*max(1.d0,sgm) .and. abs(cls_e(k)+dEi) < 1.d-4*abs(dEi)+1.d-10) then
            ic = k ; exit
          end if
        end do
        if (ic == 0 .and. ncls < 200) then
          ncls = ncls + 1 ; ic = ncls ; cls_s(ic) = sgm ; cls_e(ic) = -dEi ; cls_n(ic) = 0
        end if
        if (ic > 0) cls_n(ic) = cls_n(ic) + 1
      end do

      ! [Nov 2026] all orders of the kernel, same Gaussian holes
      allocate(Rx(nocc))
      call x_remainder(Rx)
      dExa = dEx - sum(Rx)
      ncell = dble(Kc(1)*Kc(2)*Kc(3))

      write(outfile,'(a)') ''
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      write(outfile,'(a)') '  Finite-size correction (keyword fs_correction), exchange'
      write(outfile,'(a)') '  dE_X = - sum_i sum_d (a_d^2/24) < x_d^4/r^3 >_i   (localized orbitals,'
      write(outfile,'(a)') '  Gaussian model, self terms only), leading order 1/K^2, no 1/K term'
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      write(outfile,'(a,f14.6,i8)')   '  electrons in supercell, n_occ   : ', ntot, nocc
      write(outfile,'(a,2es12.3)')    '  Loewdin check (occ min, virt max): ', ocmin, ocmax
      write(outfile,'(a,i6)')         '  Jacobi sweeps                   : ', min(isweep, MAXSWEEP)
      write(outfile,'(a,2f14.6)')     '  Omega = sum |z_ii|^2 (start/end): ', Om0, Om1
      write(outfile,'(a,f12.6,a,es14.6)') '  isotropic self-check, sigma = 1 : ', &
           x4_over_r3((/2.d0,2.d0,2.d0/),1)*3.d0/24.d0, '  exact 1/(10 sqrt(pi)) = ', 1.d0/(10.d0*dsqrt(pi))
      write(outfile,'(a)')            '  orbital classes:  count   sigma_rms (bohr)   dE per orbital (Eh)'
      do k = 1, ncls
        write(outfile,'(a,i8,f16.6,es22.10)') '                 ', cls_n(k), cls_s(k), cls_e(k)
      end do
      write(outfile,'(a,3f12.6)')     '  sigma_d mean over orbitals (x,y,z):', (dsqrt(sum(sig2(d,:))/nocc), d=1,3)
      write(outfile,'(a,f18.10,a)')   '  dE_X  (a^2 term)                : ', dEx/ncell, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_X  (all orders)              : ', dExa/ncell, ' Eh/cell'
      write(outfile,'(a)') ''
      write(outfile,'(a)') '  Summary (Eh per cell)                 a^2 term          all orders'
      write(outfile,'(a,2f18.10)')    '  dE_ES                           : ', dE_ES, dE_ES_all
      write(outfile,'(a,2f18.10)')    '  dE_X                            : ', dEx/ncell, dExa/ncell
      write(outfile,'(a,2f18.10)')    '  dE_ES + dE_X                    : ', dE_ES + dEx/ncell, dE_ES_all + dExa/ncell
      write(outfile,'(a,f18.10)')     '  E_HF per cell                   : ', EHF/ncell
      write(outfile,'(a,2f18.10)')    '  E_HF per cell - dE_ES - dE_X    : ', &
           (EHF - dEx)/ncell - dE_ES, (EHF - dExa)/ncell - dE_ES_all
      write(outfile,'(a,f10.2,a)')    '  wall time                       : ', omp_get_wtime()-w0, ' s'
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      flush(outfile)
      deallocate(z, zci, zcj, sig2, Rx)

contains

      ! [Nov 2026] supercell AO matrices from the Clifford pair Fourier transform
      !   dir = 0 : Ar = S_{mu nu}                       (m = 0 in all directions)
      !   dir = d : Ar + i Ai = <mu| exp(-i a_d x_d) |nu> (m = 1 along d, 0 elsewhere)
      ! pairs whose prefactor (Gaussian overlap factor) is below 1e-14 are skipped
      subroutine pair_matrix(dir, Ar, Ai)
        integer, intent(in) :: dir
        double precision, intent(out) :: Ar(nBas,nBas)
        double precision, intent(out), optional :: Ai(nBas,nBas)
        double precision, allocatable :: Iv(:)
        complex(8) :: dco(-8:8), cc(3), val
        integer :: mu, nu, ip, jp, dd, kk, Mm, nbes, mdir
        double precision :: r1(3), r2(3), pe(3), cen(3), ex1, ex2, A, coef(3), pref
        Ar = 0.d0
        if (present(Ai)) Ai = 0.d0
        !$omp parallel default(shared) private(mu,nu,ip,jp,dd,kk,Mm,nbes,mdir,r1,r2,pe,cen,ex1,ex2,A,coef,pref,dco,cc,val,Iv)
        allocate(Iv(0:4000))
        !$omp do schedule(dynamic,4)
        do mu = 1, nBas
          r1 = (/ AOn(mu)%x, AOn(mu)%y, AOn(mu)%z /)
          do nu = mu, nBas
            r2 = (/ AOn(nu)%x, AOn(nu)%y, AOn(nu)%z /)
            val = (0.d0,0.d0)
            do ip = 1, size(AOn(mu)%exponent)
              ex1 = AOn(mu)%exponent(ip)
              do jp = 1, size(AOn(nu)%exponent)
                ex2 = AOn(nu)%exponent(jp)
                call clifford_pair(3, ex1, ex2, r1, r2, pe, cen)
                pref = AOn(mu)%coefficient(ip)*AOn(nu)%coefficient(jp)
                do dd = 1, 3
                  coef(dd) = Lv(dd)*dexp(-2.d0*(ex1+ex2-pe(dd))/av(dd)**2)/av(dd)**(lq(dd,mu)+lq(dd,nu))
                end do
                if (abs(pref*coef(1)*coef(2)*coef(3)) < 1.d-14) cycle
                do dd = 1, 3
                  Mm = lq(dd,mu) + lq(dd,nu)
                  mdir = 0 ; if (dd == dir) mdir = 1
                  A = 2.d0*pe(dd)/av(dd)**2
                  call fourier_sin_product(lq(dd,mu), lq(dd,nu), av(dd)*(cen(dd)-r1(dd)), av(dd)*(cen(dd)-r2(dd)), dco(-Mm:Mm))
                  nbes = min(4000, Mm + 40 + int(10.d0*dsqrt(A)))   ! backward recurrence: start well above sqrt(A)
                  call bessel_scaled_vector_fast(nbes, A, Iv(0:nbes))
                  cc(dd) = (0.d0,0.d0)
                  do kk = -Mm, Mm
                    cc(dd) = cc(dd) + dco(kk)*Iv(abs(mdir-kk))
                  end do
                  cc(dd) = coef(dd)*cc(dd)
                  if (mdir == 1) cc(dd) = cc(dd)*exp(dcmplx(0.d0, -av(dd)*cen(dd)))
                end do
                val = val + pref*cc(1)*cc(2)*cc(3)
              end do
            end do
            Ar(mu,nu) = dble(val) ; Ar(nu,mu) = dble(val)
            if (present(Ai)) then
              Ai(mu,nu) = aimag(val) ; Ai(nu,mu) = aimag(val)
            end if
          end do
        end do
        !$omp end do
        deallocate(Iv)
        !$omp end parallel
      end subroutine pair_matrix

      ! [Nov 2026] all orders minus a^2, per orbital (Gaussian hole, variances s_d^2 = 2 sigma_d^2)
      !   <1/d> - <1/r> = (2/sqrt pi) int dt [ prod_d F_d - prod_d B_d ]
      !   F_d = sum_k I~_k(z_d) exp(-k^2 a_d^2 s_d^2/2), z_d = 2t^2/a_d^2  (exact torus, wrapped Gaussian)
      !   B_d = (1+2t^2 s_d^2)^-1/2
      !   a^2 part of the integrand: t^2 sum_d (a_d^2/12) 3 s_d^4 (1+2t^2s_d^2)^-5/2 prod_{k/=d} B_k
      subroutine x_remainder(R)
        double precision, intent(out) :: R(nocc)
        double precision, allocatable :: Ivt(:,:), s2a(:,:)
        double precision :: smin, smax, ds, t, wt, zz, F(3), Bq(3), a2t, smn, smx
        integer :: ns, is, dd, nb(3), io, kk, kmx, nbm
        allocate(s2a(3,nocc))
        s2a = 2.d0*sig2
        smn = dsqrt(minval(s2a)) ; smx = dsqrt(maxval(s2a))
        smin = dlog(1.d-3*minval(av)) ; smax = dlog(200.d0/smn)
        ds = 0.01d0 ; ns = int((smax-smin)/ds) + 1
        R = 0.d0
        !$omp parallel default(shared) private(is,t,wt,dd,zz,nb,nbm,Ivt,io,F,Bq,a2t,kk,kmx)
        !$omp do schedule(dynamic,4) reduction(+:R)
        do is = 0, ns
          t = dexp(smin + dble(is)*ds)
          wt = t*ds ; if (is == 0 .or. is == ns) wt = 0.5d0*wt
          do dd = 1, 3
            zz = 2.d0*t*t/av(dd)**2
            nb(dd) = int(12.d0/(av(dd)*smn)) + 40 + int(10.d0*dsqrt(zz))
          end do
          nbm = maxval(nb)
          allocate(Ivt(0:nbm,3))
          do dd = 1, 3
            call bessel_scaled_vector_fast(nb(dd), 2.d0*t*t/av(dd)**2, Ivt(0:nb(dd),dd))
          end do
          do io = 1, nocc
            do dd = 1, 3
              kmx = min(nb(dd), int(12.d0/(av(dd)*dsqrt(s2a(dd,io)))) + 30)
              F(dd) = Ivt(0,dd)
              do kk = 1, kmx
                F(dd) = F(dd) + 2.d0*Ivt(kk,dd)*dexp(-0.5d0*dble(kk*kk)*av(dd)**2*s2a(dd,io))
              end do
              Bq(dd) = 1.d0/dsqrt(1.d0 + 2.d0*t*t*s2a(dd,io))
            end do
            a2t = 0.d0
            do dd = 1, 3
              a2t = a2t + av(dd)**2/12.d0*3.d0*s2a(dd,io)**2*Bq(dd)**5*product(Bq)/Bq(dd)
            end do
            R(io) = R(io) + wt*(F(1)*F(2)*F(3) - Bq(1)*Bq(2)*Bq(3) - t*t*a2t)
          end do
          deallocate(Ivt)
        end do
        !$omp end do
        !$omp end parallel
        R = 2.d0/dsqrt(pi)*R
        deallocate(s2a)
      end subroutine x_remainder

      double precision function omega_val()
        integer :: ii, dd
        omega_val = 0.d0
        do dd = 1, 3
          do ii = 1, nocc
            omega_val = omega_val + abs(z(ii,ii,dd))**2
          end do
        end do
      end function omega_val

      ! < x_d^4 / r^3 > for r Gaussian with variances s2(1:3)
      ! (4/sqrt(pi)) int_0^inf t^2 3 s_d^4 (1+2t^2 s_d^2)^-5/2 prod_{k/=d}(1+2t^2 s_k^2)^-1/2 dt
      ! t = tau u/(1-u), composite Simpson on u in [0,1]
      double precision function x4_over_r3(s2v, dd)
        double precision, intent(in) :: s2v(3)
        integer, intent(in) :: dd
        integer, parameter :: NQ = 4000
        integer :: q, kk
        double precision :: tau, uq, t, f, wq, ssum
        tau = 1.d0/dsqrt(sum(s2v)/3.d0)
        ssum = 0.d0
        do q = 1, NQ-1                                  ! end points: integrand 0
          uq = dble(q)/dble(NQ)
          t = tau*uq/(1.d0-uq)
          f = t*t*3.d0*s2v(dd)**2*(1.d0+2.d0*t*t*s2v(dd))**(-2.5d0)
          do kk = 1, 3
            if (kk /= dd) f = f/dsqrt(1.d0+2.d0*t*t*s2v(kk))
          end do
          f = f*tau/(1.d0-uq)**2
          wq = 2.d0 ; if (mod(q,2) == 1) wq = 4.d0
          ssum = ssum + wq*f
        end do
        x4_over_r3 = 4.d0/dsqrt(pi)*ssum/(3.d0*dble(NQ))
      end function x4_over_r3

end subroutine fs_exchange_3D

