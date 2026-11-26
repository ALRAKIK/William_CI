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

      ! [Nov 2026] adaptive G cutoff for the electrons: a primitive pair (exponent p) contributes
      ! up to |G_d| = 10 sqrt(p) (its transform ~ exp(-G^2/4p) < 1e-11); NG = max over pairs,
      ! capped at NG_CAP. Before: fixed |n| <= 10, too small for large cells (e.g. atom in a box).
      integer, parameter :: NG_CAP = 200
      integer, parameter :: NR_MAX = 16                ! exact-kernel remainder: |n_i| <= min(NG, NR_MAX)
      type(ERI_function), allocatable :: AOn(:)
      integer :: nbes, Kc(3), lq(3,nBas), mu, nu, ip, jp, d, k, n, m, Mm, nmx, i1, i2, i3, nat_cell, j
      integer :: NG, NR, NGN(5), nl(3), npair_cap, nthr
      double precision :: pmax, neglect
      double precision :: Lv(3), av(3), cell(3), Omega, pe(3), cen(3), r1(3), r2(3), pref, ex1, ex2, A
      double precision :: G(3), dv, e_nn, e_ne, e_ee, nel0, w0, e_nn_in, e_nn_tail
      double precision :: e_nn_cut(5), coef(3), r_nn, r_ne, r_ee
      double precision, allocatable :: Rg(:,:,:)
      complex(8), allocatable :: rho_el(:,:,:), cm(:,:), dco(:)
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

      ! ---------------- G cutoff for the electrons ---------------- !
      NG = 4 ; npair_cap = 0 ; pmax = 0.d0
      do mu = 1, nfuc
        do ip = 1, size(AOn(mu)%exponent)
          do jp = 1, size(AOn(mu)%exponent)        ! tightest pairs are on the same function / atom
            ex1 = AOn(mu)%exponent(ip) + AOn(mu)%exponent(jp)
            pmax = max(pmax, ex1)
            do d = 1, 3
              n = ceiling(10.d0*dsqrt(ex1)*cell(d)/(2.d0*acos(-1.d0))) + 2*lq(d,mu)
              if (n > NG_CAP) npair_cap = npair_cap + 1
              NG = max(NG, min(n, NG_CAP))
            end do
          end do
        end do
      end do
      NR = min(NG, NR_MAX)
      ! largest neglected transform factor exp(-G^2/4p) at the cutoff (1e-11 if not capped)
      neglect = dexp(-(2.d0*acos(-1.d0)*dble(NG)/maxval(cell))**2/(4.d0*pmax))
      ! the reduction copies one grid per thread: limit to ~2 GB in total
      nthr = max(1, min(omp_get_max_threads(), int(2.d9/(16.d0*dble(2*NG+1)**3))))

      ! ---------------- electrons: rho_el(n1,n2,n3), |n_i| <= NG ---------------- !
      allocate(rho_el(-NG:NG,-NG:NG,-NG:NG))
      rho_el = (0.d0,0.d0)
      nmx = NG*maxval(Kc) + 12
      !$omp parallel num_threads(nthr) default(shared) private(nbes,mu,nu,ip,jp,d,k,n,m,Mm,pe,cen,r1,r2,pref,ex1,ex2,A,cm,dco,Iv,i1,i2,i3,ph,coef,nl)
      allocate(cm(-NG:NG,3), dco(-8:8), Iv(0:40000))
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
                coef(d) = Lv(d)*dexp(-2.d0*(ex1+ex2-pe(d))/av(d)**2)/av(d)**(lq(d,mu)+lq(d,nu))
              end do
              if (abs(pref*coef(1)*coef(2)*coef(3)) < 1.d-14) cycle
              do d = 1, 3
                Mm = lq(d,mu) + lq(d,nu)
                ! this pair's own cutoff along d (transform ~ exp(-G^2/4p))
                nl(d) = min(NG, ceiling(10.d0*dsqrt(ex1+ex2)*cell(d)/(2.d0*acos(-1.d0))) + Mm)
                A = 2.d0*pe(d)/av(d)**2
                call fourier_sin_product(lq(d,mu), lq(d,nu), av(d)*(cen(d)-r1(d)), av(d)*(cen(d)-r2(d)), dco(-Mm:Mm))
                nbes = min(40000, max(nl(d)*Kc(d)+12+Mm, Mm + 40 + int(10.d0*dsqrt(A))))   ! backward recurrence start
                call bessel_scaled_vector_fast(nbes, A, Iv(0:nbes))
                do n = -nl(d), nl(d)
                  m = n*Kc(d)
                  cm(n,d) = (0.d0,0.d0)
                  do k = -Mm, Mm
                    cm(n,d) = cm(n,d) + dco(k)*Iv(abs(m-k))
                  end do
                  cm(n,d) = coef(d)*cm(n,d)*exp(dcmplx(0.d0, -dble(m)*av(d)*cen(d)))
                end do
              end do
              do i3 = -nl(3), nl(3)
                do i2 = -nl(2), nl(2)
                  ph = pref*cm(i2,2)*cm(i3,3)
                  do i1 = -nl(1), nl(1)
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

      ! ---------------- sums over G ---------------- !
      ! |n_i| <= NG : full rho = rho_nuc - rho_el (a^2 kernel)
      ! |n_i| >  NG : nuclei only (electrons have decayed): point-charge tail, from
      !               cutoffs NGN(1..5) extrapolated c0 + c1/n + c2/n^2, minus the part <= NG
      e_nn_in = 0.d0 ; e_ne = 0.d0 ; e_ee = 0.d0
      !$omp parallel do default(shared) private(i1,i2,i3,G,dv,zn,ze,j) reduction(+:e_nn_in,e_ne,e_ee) collapse(2)
      do i3 = -NG, NG
        do i2 = -NG, NG
          do i1 = -NG, NG
            if (i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
            G = 2.d0*acos(-1.d0)*(/ dble(i1), dble(i2), dble(i3) /)/cell
            dv = dvG(G)
            zn = (0.d0,0.d0)
            do j = 1, nat_cell
              zn = zn + dble(charge(j))*exp(dcmplx(0.d0, -sum(G*geometry(j,:))))
            end do
            ze = rho_el(i1,i2,i3)
            e_nn_in = e_nn_in + dv*dble(zn*conjg(zn))
            e_ne = e_ne + dv*(-2.d0*dble(zn*conjg(ze)))
            e_ee = e_ee + dv*dble(ze*conjg(ze))
          end do
        end do
      end do
      !$omp end parallel do
      NGN = max(30, NG+10)*(/ 3, 4, 5, 6, 7 /)/3
      do k = 1, 5
        e_nn_cut(k) = nuc_sum(NGN(k))
      end do
      e_nn_tail = extrap(e_nn_cut) - e_nn_in
      e_nn = (e_nn_in + e_nn_tail)/(2.d0*Omega)
      e_ne = e_ne/(2.d0*Omega) ; e_ee = e_ee/(2.d0*Omega)

      ! [Nov 2026] higher orders of the kernel: R(G) = dv_exact(G) - dv_a2(G), |n_i| <= NR
      ! (tail |n_i| > NR of R ~ dv/n^2 neglected)
      allocate(Rg(0:NR,0:NR,0:NR))
      call es_remainder(Rg)
      r_nn = 0.d0 ; r_ne = 0.d0 ; r_ee = 0.d0
      do i3 = -NR, NR
        do i2 = -NR, NR
          do i1 = -NR, NR
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
      write(outfile,'(a,i4,a,i4)')    '  electrons: G cutoff |n_i| <=    : ', NG, '   remainder: ', NR
      if (npair_cap > 0) write(outfile,'(a,i4,a,es9.1,a)') '  cutoff capped at ', NG_CAP, &
                         ': neglected transform factor ', neglect, ' (tightest pair)'
      write(outfile,'(a,5i5)')        '  nuclei tail: cutoffs            : ', NGN
      write(outfile,'(a,5f14.8)')     '  nuc-nuc at these cutoffs (x 1/2Omega):', e_nn_cut/(2.d0*Omega)
      write(outfile,'(a,f18.10,a)')   '  dE nuc-nuc (extrapolated)       : ', e_nn, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE nuc-el                       : ', e_ne, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE el-el                        : ', e_ee, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_ES total                     : ', e_nn+e_ne+e_ee, ' Eh/cell'
      write(outfile,'(a)')            '  higher orders of the kernel (exact - a^2 term):'
      write(outfile,'(a,f18.10,a)')   '    nuc-nuc                       : ', r_nn, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '    nuc-el                        : ', r_ne, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '    el-el                         : ', r_ee, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_ES all orders                : ', e_nn+e_ne+e_ee+r_nn+r_ne+r_ee, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  E_HF per cell                   : ', EHF/dble(Kc(1)*Kc(2)*Kc(3)), ' Eh'
      write(outfile,'(a,f18.10,a)')   '  E_HF per cell - dE_ES (a^2)     : ', EHF/dble(Kc(1)*Kc(2)*Kc(3)) - (e_nn+e_ne+e_ee), ' Eh'
      write(outfile,'(a,f10.2,a)')    '  wall time                       : ', omp_get_wtime()-w0, ' s'
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      flush(outfile)
      deallocate(rho_el, Rg)

      ! [Nov 2026] exchange part (localized orbitals, Gaussian model)
      call fs_exchange_3D(nBas, AOn, lq, P, Lv, av, Kc, EHF, e_nn+e_ne+e_ee, e_nn+e_ne+e_ee+r_nn+r_ne+r_ee)
      deallocate(AOn)

contains

      ! [Nov 2026] R(G) = dv_exact(G) - dv_a2(G) for G = 2 pi n / cell, 0 <= n_i <= NR
      !   dv_exact = (2/sqrt pi) int dt [ prod_d L_d I~_{m_d}(z_d) - P0 ],  P0 = pi^3/2 t^-3 exp(-G^2/4t^2)
      !   dv_a2    = (2/sqrt pi) int dt  P0 sum_d eps_d                     (analytic: dvG)
      !   eps_d    = a_d^2/(48 t^2) (G_d^4/(4 t^4) - 3 G_d^2/t^2 + 3),  m_d = n_d K_d,  z_d = 2 t^2/a_d^2
      ! t = exp(s), trapezoid in s.
      subroutine es_remainder(R)
        double precision, intent(out) :: R(0:NR,0:NR,0:NR)
        double precision, allocatable :: Ivt(:), fac(:,:)
        double precision :: gmin, gmax, smin, smax, ds, t, wt, z, P0, eps, gv(3), g2, pi
        integer :: ns, is, dd, nb, j1, j2, j3
        pi = acos(-1.d0)
        gmin = 2.d0*pi/maxval(cell)
        gmax = 2.d0*pi*dble(NR)*dsqrt(sum(1.d0/cell**2))
        smin = dlog(gmin/40.d0) ; smax = dlog(20.d0*gmax)
        ds = 0.01d0 ; ns = int((smax-smin)/ds) + 1
        R = 0.d0
        !$omp parallel default(shared) private(is,t,wt,dd,z,nb,Ivt,fac,j1,j2,j3,gv,g2,P0,eps)
        allocate(fac(0:NR,3))
        !$omp do schedule(dynamic,4) reduction(+:R)
        do is = 0, ns
          t = dexp(smin + dble(is)*ds)
          wt = t*ds ; if (is == 0 .or. is == ns) wt = 0.5d0*wt
          do dd = 1, 3
            z = 2.d0*t*t/av(dd)**2
            nb = NR*Kc(dd) + 40 + int(10.d0*dsqrt(z))
            allocate(Ivt(0:nb))
            call bessel_scaled_vector_fast(nb, z, Ivt)
            do j1 = 0, NR
              fac(j1,dd) = Lv(dd)*Ivt(j1*Kc(dd))
            end do
            deallocate(Ivt)
          end do
          do j3 = 0, NR
            do j2 = 0, NR
              do j1 = 0, NR
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

      double precision function nuc_sum(nc)       ! sum_{0<|n_i|<=nc} |rho_nuc(G)|^2 dv(G)
        integer, intent(in) :: nc
        integer :: q1, q2, q3, jj
        double precision :: Gq(3)
        complex(8) :: zq
        nuc_sum = 0.d0
        !$omp parallel do default(shared) private(q1,q2,q3,Gq,zq,jj) reduction(+:nuc_sum) collapse(2)
        do q3 = -nc, nc
          do q2 = -nc, nc
            do q1 = -nc, nc
              if (q1 == 0 .and. q2 == 0 .and. q3 == 0) cycle
              Gq = 2.d0*acos(-1.d0)*(/ dble(q1), dble(q2), dble(q3) /)/cell
              zq = (0.d0,0.d0)
              do jj = 1, nat_cell
                zq = zq + dble(charge(jj))*exp(dcmplx(0.d0, -sum(Gq*geometry(jj,:))))
              end do
              nuc_sum = nuc_sum + dvG(Gq)*dble(zq*conjg(zq))
            end do
          end do
        end do
        !$omp end parallel do
      end function nuc_sum

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
      use keywords,           only : c_fs_exact
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
                          B, Cc, th, cs, sn, gain, gsum, Om0, Om1, dEx, s2(3), dEi, sgm, ocmin, ocmax, dExa, ncell, &
                          dXp, d2p(3), sp2p(3)
      integer :: npair
      double precision, external :: fs_pair_model
      double precision, allocatable :: cci(:), ccj(:)
      double precision :: XB_self, XB_pair, XB_lat, XB_part(3), XB_chk, XB_kap(3), XB_out, w0B
      integer :: XB_nhome, XB_npair, XB_M(3)
      double precision :: XB_gcut, XB_spart(3), XB_tthr, XB_nterm
      integer, allocatable :: nlist(:), alist(:,:)         ! [Nov 2026] option B: AO lists of the localized orbitals
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
      deallocate(S, Zi, T1, T2)            ! [Nov 2026] C kept: localized orbitals for option B

      ! ---------- Jacobi localization: maximize sum_d sum_i |z_ii^d|^2 ---------- !
      Om0 = omega_val()
      allocate(zci(nocc), zcj(nocc), cci(nBas), ccj(nBas))
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
            cci = C(:,i) ; ccj = C(:,j)                 ! [Nov 2026] same rotation on the orbitals
            C(:,i) =  cs*cci + sn*ccj
            C(:,j) = -sn*cci + cs*ccj
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

      ! [Nov 2026] option A: pair terms i /= j, dipolar-Gaussian model
      !   |d_k|^2 = |z^k_ij|^2 / a_k^2 ,  pair width sp2 = 2 s_i s_j/(s_i + s_j)
      dXp = 0.d0 ; npair = 0
      !$omp parallel do default(shared) private(i,j,d,d2p,sp2p) reduction(+:dXp,npair) schedule(dynamic,4)
      do i = 1, nocc-1
        do j = i+1, nocc
          do d = 1, 3
            d2p(d)  = abs(z(i,j,d))**2/av(d)**2
            sp2p(d) = 2.d0*sig2(d,i)*sig2(d,j)/(sig2(d,i) + sig2(d,j))
          end do
          if (sum(d2p) < 1.d-10) cycle
          npair = npair + 1
          dXp = dXp - 2.d0*fs_pair_model(sp2p, d2p, av, 3)      ! ordered pairs ij and ji
        end do
      end do
      !$omp end parallel do

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
      write(outfile,'(a)') ''
      write(outfile,'(a)') '  Exchange refinement A: pair terms i /= j (dipolar-Gaussian model)'
      write(outfile,'(a,i10)')        '  pairs with |d_ij|^2 > 1e-10     : ', npair
      write(outfile,'(a,f18.10,a)')   '  dE_X pairs (model A)            : ', dXp/ncell, ' Eh/cell'
      write(outfile,'(a,f18.10)')     '  E_HF - dE_ES - dE_X - dE_Xpair  : ', (EHF - dExa - dXp)/ncell - dE_ES_all

      ! [Nov 2026] option B: exact first-order exchange (all orders in the kernel), keyword fs_correction_exact
      if (c_fs_exact) then
      w0B = omp_get_wtime()
      call exchange_B()
      write(outfile,'(a)') ''
      write(outfile,'(a)') '  Exchange refinement B: exact pair terms + self-term shape (all orders)'
      write(outfile,'(a,3i5)')        '  supercell g grid |m_d| <=       : ', XB_M
      write(outfile,'(a,3f12.6)')     '  lattice dipole factors kappa    : ', XB_kap
      write(outfile,'(a,2i8)')        '  reference-cell orbitals, pairs  : ', XB_nhome, XB_npair
      write(outfile,'(a,es12.3)')     '  check max ||rho_ii(a e_d)|-|z_ii||: ', XB_chk
      write(outfile,'(a,es10.1,f12.0)') '  term threshold, mean terms/pair : ', XB_tthr, XB_nterm/dble(max(XB_npair+XB_nhome,1))
      write(outfile,'(a,f8.2)')       '  g cutoff (bohr^-1)              : ', XB_gcut
      write(outfile,'(a,3f14.10)')    '  self shape at 0.6/0.8/1.0 g_cut : ', XB_spart
      write(outfile,'(a,3f14.10)')    '  pair kernel at 0.6/0.8/1.0 g_cut: ', XB_part
      write(outfile,'(a,f18.10,a)')   '  dE_X self-term shape (extrap 1/g)   : ', XB_self, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_X pairs, kernel (extrap 1/g^2)   : ', XB_pair, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_X pairs, image (dipole) part : ', XB_lat, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_X pairs from outer 25% of R  : ', XB_out, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_X (B) total                  : ', dExa/ncell + XB_self + XB_pair + XB_lat, ' Eh/cell'
      write(outfile,'(a,f18.10)')     '  E_HF - dE_ES - dE_X(B)          : ', &
           EHF/ncell - dE_ES_all - (dExa/ncell + XB_self + XB_pair + XB_lat)
      write(outfile,'(a,f10.2,a)')    '  wall time (B)                   : ', omp_get_wtime()-w0B, ' s'
      else
        write(outfile,'(a)') ''
        write(outfile,'(a)') '  (exact exchange, option B: keyword fs_correction_exact)'
      end if
      write(outfile,'(a,f10.2,a)')    '  wall time                       : ', omp_get_wtime()-w0, ' s'
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      flush(outfile)
      deallocate(z, zci, zcj, sig2, Rx, C, cci, ccj)

contains

      !=================================================================!
      ! [Nov 2026] option B: exact first-order exchange correction      !
      !   dE_X = - sum_i [ d_Gauss(i) + S_i ] - sum_{i/=j} [ P_ij + L_ij ] !
      !   S_i  = (1/V) sum'_g ( |rho_ii(g)|^2 - prod_d |z_ii^d|^(2 m_d^2) ) D(g)
      !   P_ij = (1/V) sum'_g |rho_ij(g)|^2 D(g)                        !
      !   L_ij = - (4 pi/V) sum_k kappa_k |z_ij^k|^2 / a_k^2            !
      !   D(g) = v_L(g) - 4 pi/g^2  (exact torus coefficient, t-integral)
      !   g = (m_d a_d) supercell lattice, |g| <= g_cut                 !
      !   rho_ij(g) = sum_{mu nu} C_mu,i C_nu,j (mu nu)(g)  (Clifford pair FT)
      !   Only orbitals i centred in one reference cell (translation    !
      !   copies); per cell = (n_occ/ncell)/n_home * sum_{i home}.      !
      !=================================================================!
      subroutine exchange_B()
        double precision :: GCUT
        double precision :: CTHR
        character(len=32) :: envs
        integer :: ln, ist
        double precision, allocatable :: Dt(:,:,:), cen_o(:,:), sig_o(:)
        integer, allocatable :: home(:)
        double precision :: Vs, cellv(3), dd(3), dist, rcut, wfac, gc(3), g2
        double precision :: Sacc(3), Pacc, Lacc, Pp(3), chkv, Oacc
        integer :: Mv(3), ii, jj, kk, ih, nh, np_loc, mmax, Wc
        double precision :: wn(3)
        complex(8), allocatable :: grid(:,:,:)


        Vs = Lv(1)*Lv(2)*Lv(3)
        cellv = Lv/dble(Kc)
        GCUT = 12.d0
        call get_environment_variable('FS_GCUT', envs, ln, ist)          ! test hook
        if (ist == 0 .and. ln > 0) read(envs,*) GCUT
        XB_gcut = GCUT
        XB_tthr = 1.d-9
        CTHR = 1.d-3
        call get_environment_variable('FS_CTHR', envs, ln, ist)          ! test hook
        if (ist == 0 .and. ln > 0) read(envs,*) CTHR
        call get_environment_variable('FS_TTHR', envs, ln, ist)          ! test hook
        if (ist == 0 .and. ln > 0) read(envs,*) XB_tthr
        XB_nterm = 0.d0
        do kk = 1, 3
          Mv(kk) = ceiling(GCUT/av(kk))
        end do
        XB_M = Mv
        allocate(Dt(0:Mv(1),0:Mv(2),0:Mv(3)))
        call kernel_D(Mv, Dt)
        call lattice_kappa(XB_kap)

        ! centres, sizes, reference-cell orbitals
        allocate(cen_o(3,nocc), sig_o(nocc), home(nocc))
        nh = 0
        do ii = 1, nocc
          do kk = 1, 3
            cen_o(kk,ii) = modulo(-atan2(aimag(z(ii,ii,kk)), dble(z(ii,ii,kk)))/av(kk), Lv(kk))
          end do
          sig_o(ii) = dsqrt(sum(sig2(:,ii)))
        end do
        ! [Nov 2026] reference window: W cells along x (one cell along y, z), smallest W holding an orbital centre
        do Wc = 1, Kc(1)
          nh = 0
          do ii = 1, nocc
            wn = cellv ; wn(1) = dble(Wc)*cellv(1)
            if (all(modulo(cen_o(:,ii) - 0.1234d0*cellv, Lv) < wn)) then
              nh = nh + 1 ; home(nh) = ii
            end if
          end do
          if (nh > 0) exit
        end do
        XB_nhome = nh
        wfac = dble(nocc)/(dble(ncell)*dble(nh))

        ! AO lists per orbital (|C| > CTHR * max|C|)
        allocate(nlist(nocc), alist(nBas,nocc))
        do ii = 1, nocc
          nlist(ii) = 0
          do kk = 1, nBas
            if (abs(C(kk,ii)) > CTHR*maxval(abs(C(:,ii)))) then
              nlist(ii) = nlist(ii) + 1 ; alist(nlist(ii),ii) = kk
            end if
          end do
        end do

        Sacc = 0.d0 ; Pacc = 0.d0 ; Lacc = 0.d0 ; Pp = 0.d0 ; chkv = 0.d0 ; Oacc = 0.d0 ; np_loc = 0
        do ih = 1, nh
          call home_block(home(ih), Mv, Dt, cen_o, sig_o, Vs, Sacc, Pacc, Pp, Lacc, Oacc, chkv, np_loc)
        end do

        ! dE_X = - sum [...] : signs and per-cell weights
        XB_spart = -wfac*Sacc
        XB_part = -wfac*Pp
        ! tail extrapolation from the partial sums at 0.8 and 1.0 g_cut:
        !   self terms ~ 1/g (tight cores), pair terms ~ 1/g^2
        XB_self = (1.0d0*XB_spart(3) - 0.8d0*XB_spart(2))/(1.0d0 - 0.8d0)
        XB_pair = (1.0d0*XB_part(3) - 0.64d0*XB_part(2))/(1.0d0 - 0.64d0)
        XB_lat  = -wfac*Lacc
        XB_out  = -wfac*Oacc
        XB_chk  = chkv
        XB_npair = np_loc
        deallocate(Dt, cen_o, sig_o, home, nlist, alist)
      end subroutine exchange_B


      !=================================================================!
      ! [Nov 2026] option B, one reference orbital i and all its pairs j !
      !   H_nu(g)  = sum_mu C_mu,i (mu nu)(g)      (half transform)      !
      !   rho_ij(g) = sum_nu C_nu,j H_nu(g)         (zgemm, blocks of nu)!
      ! The basis-function pair transforms are computed once per i.     !
      ! Grid: |m_1|,|m_2| <= M, 0 <= m_3 <= M (rho(-g) = rho(g)*).        !
      !=================================================================!
      subroutine home_block(io, Mv, Dt, cen_o, sig_o, Vs, Sacc, Pacc, Pp, Lacc, Oacc, chkv, np_loc)
        integer, intent(in) :: io, Mv(3)
        double precision, intent(in) :: Dt(0:Mv(1),0:Mv(2),0:Mv(3)), cen_o(3,nocc), sig_o(nocc), Vs
        double precision, intent(inout) :: Sacc(3), Pacc, Pp(3), Lacc, Oacc, chkv
        integer, intent(inout) :: np_loc
        integer, allocatable :: jl(:), nub(:)
        double precision, allocatable :: gj(:), rc(:), dj(:)
        complex(8), allocatable :: R(:,:), Hb(:,:), Cb(:,:)
        integer :: nj, jj, kk, Mi(3), NG1, NG2, NG3, Gt, ax1(3), ax2(3), nbk, b0, b1, nbl, ib, nu, jx, m1, m2, m3, idx
        double precision :: dd(3), dist, rcut, gmax, g2, wm, r2, gau, e_s(3), e_p(3), e_l, y1, y2, c0, f1, f2, d2k
        complex(8) :: v1, v2

        ! pairs j of this i
        allocate(jl(nocc), gj(nocc), rc(nocc), dj(nocc))
        nj = 0
        do jj = 1, nocc
          dd = cen_o(:,jj) - cen_o(:,io)
          dd = dd - Lv*dnint(dd/Lv)
          dist = dsqrt(sum(dd**2))
          rcut = 4.d0*(sig_o(io) + sig_o(jj)) + 2.d0
          if (dist > rcut .and. jj /= io) cycle
          nj = nj + 1 ; jl(nj) = jj ; rc(nj) = rcut ; dj(nj) = dist
          gj(nj) = min(XB_gcut, max(4.d0, 5.d0/dsqrt(min(minval(sig2(:,io)), minval(sig2(:,jj))))))
        end do
        gmax = maxval(gj(1:nj))
        do kk = 1, 3
          Mi(kk) = max(2, ceiling(gmax/av(kk)))
        end do
        NG1 = 2*Mi(1)+1 ; NG2 = 2*Mi(2)+1 ; NG3 = Mi(3)+1 ; Gt = NG1*NG2*NG3
        do kk = 1, 3                                       ! flat indices of e_k and 2 e_k
          ax1(kk) = flat_idx(Mi, kk, 1) ; ax2(kk) = flat_idx(Mi, kk, 2)
        end do
        allocate(R(Gt,nj)) ; R = (0.d0,0.d0)

        ! blocks of nu: H for a block, then R += H * C(nu,j)
        nbk = max(4, min(128, int(3.d8/(16.d0*dble(Gt)))))
        allocate(nub(nBas))
        nbl = 0
        do nu = 1, nBas                                   ! nu relevant to at least one j
          if (maxval(abs(C(nu,jl(1:nj)))) < 1.d-6) cycle
          nbl = nbl + 1 ; nub(nbl) = nu
        end do
        allocate(Hb(Gt,nbk), Cb(nbk,nj))
        do b0 = 1, nbl, nbk
          b1 = min(nbl, b0+nbk-1)
          !$omp parallel do default(shared) private(ib) schedule(dynamic,1)
          do ib = b0, b1
            call half_transform(io, nub(ib), Mi, Hb(:,ib-b0+1))
          end do
          !$omp end parallel do
          do ib = b0, b1
            Cb(ib-b0+1,:) = dcmplx(C(nub(ib),jl(1:nj)), 0.d0)
          end do
          call zgemm('N', 'N', Gt, nj, b1-b0+1, (1.d0,0.d0), Hb, Gt, Cb, nbk, (1.d0,0.d0), R, Gt)
        end do

        ! energies per pair
        !$omp parallel do default(shared) private(jx,jj,e_s,e_p,e_l,m1,m2,m3,idx,g2,wm,r2,gau,kk,y1,y2,c0,f1,f2,d2k,v1,v2) &
        !$omp reduction(+:Sacc,Pacc,Pp,Lacc,Oacc,np_loc) reduction(max:chkv) schedule(dynamic,1)
        do jx = 1, nj
          jj = jl(jx)
          e_s = 0.d0 ; e_p = 0.d0
          do m3 = 0, Mi(3)
            do m2 = -Mi(2), Mi(2)
              do m1 = -Mi(1), Mi(1)
                if (m1 == 0 .and. m2 == 0 .and. m3 == 0) cycle
                if (m3 == 0 .and. (m2 < 0 .or. (m2 == 0 .and. m1 < 0))) cycle
                wm = 2.d0
                g2 = (m1*av(1))**2 + (m2*av(2))**2 + (m3*av(3))**2
                if (g2 > gj(jx)**2) cycle
                idx = (m1+Mi(1)+1) + NG1*(m2+Mi(2)) + NG1*NG2*m3
                if (jj == io) then
                  gau = (abs(z(io,io,1))**2)**(m1*m1)*(abs(z(io,io,2))**2)**(m2*m2)*(abs(z(io,io,3))**2)**(m3*m3)
                  r2 = wm*(abs(R(idx,jx))**2 - gau)*Dt(abs(m1),abs(m2),abs(m3))
                  if (g2 <= (0.6d0*gj(jx))**2) e_s(1) = e_s(1) + r2
                  if (g2 <= (0.8d0*gj(jx))**2) e_s(2) = e_s(2) + r2
                  e_s(3) = e_s(3) + r2
                else
                  r2 = wm*abs(R(idx,jx))**2*Dt(abs(m1),abs(m2),abs(m3))
                  if (g2 <= (0.6d0*gj(jx))**2) e_p(1) = e_p(1) + r2
                  if (g2 <= (0.8d0*gj(jx))**2) e_p(2) = e_p(2) + r2
                  e_p(3) = e_p(3) + r2
                end if
              end do
            end do
          end do
          e_s = e_s/Vs ; e_p = e_p/Vs
          if (jj == io) then
            do kk = 1, 3
              v1 = R(ax1(kk),jx) ; v2 = R(ax2(kk),jx)
              chkv = max(chkv, abs(abs(v1) - abs(z(io,io,kk))))
              y1 = dlog(max(abs(v1)**2, 1.d-300)) ; y2 = dlog(max(abs(v2)**2, 1.d-300))
              c0 = -y1/av(kk)**2 + (16.d0*y1 - y2)/(12.d0*av(kk)**2)
              e_s = e_s - 4.d0*acos(-1.d0)/Vs*XB_kap(kk)*c0
            end do
            Sacc = Sacc + e_s
          else
            np_loc = np_loc + 1
            Pacc = Pacc + e_p(3) ; Pp = Pp + e_p
            e_l = 0.d0
            do kk = 1, 3
              f1 = abs(R(ax1(kk),jx))**2/av(kk)**2
              f2 = abs(R(ax2(kk),jx))**2/(4.d0*av(kk)**2)
              if (f1 > 1.d-300 .and. f2 > 1.d-300) then
                d2k = dexp((4.d0*dlog(f1) - dlog(f2))/3.d0)
              else
                d2k = f1
              end if
              e_l = e_l - 4.d0*acos(-1.d0)/Vs*XB_kap(kk)*d2k
            end do
            Lacc = Lacc + e_l
            if (dj(jx) > 0.75d0*rc(jx)) Oacc = Oacc + e_p(3) + e_l
          end if
        end do
        !$omp end parallel do
        deallocate(jl, gj, rc, dj, R, Hb, Cb, nub)

      end subroutine home_block

      integer function flat_idx(Mi, k, m)              ! flat index of m e_k on the half grid
        integer, intent(in) :: Mi(3), k, m
        integer :: q(3)
        q = 0 ; q(k) = m
        flat_idx = (q(1)+Mi(1)+1) + (2*Mi(1)+1)*(q(2)+Mi(2)) + (2*Mi(1)+1)*(2*Mi(2)+1)*q(3)
      end function flat_idx

      ! H_nu(g) = sum_{mu in list(io)} C_mu,io (mu nu)(g), flat grid (m1 fastest, m3 >= 0)
      subroutine half_transform(io, nu, Mi, H)
        integer, intent(in) :: io, nu, Mi(3)
        complex(8), intent(out) :: H(2*Mi(1)+1, 2*Mi(2)+1, Mi(3)+1)
        complex(8), allocatable :: F1(:,:), F2(:,:), F3(:,:), B2(:,:)
        double precision, allocatable :: Iv(:)
        complex(8) :: dco(-8:8), cv
        double precision :: r1(3), r2(3), pe(3), cen(3), e1, e2, A, coef(3), pref
        integer :: a1, mu, ip, jp, dd, Mm, k, m, nT, nTmax, nb, m3
        nTmax = 0
        do a1 = 1, nlist(io)
          nTmax = nTmax + size(AOn(alist(a1,io))%exponent)*size(AOn(nu)%exponent)
        end do
        allocate(F1(-Mi(1):Mi(1),nTmax), F2(-Mi(2):Mi(2),nTmax), F3(0:Mi(3),nTmax), Iv(0:40000))
        nT = 0
        r2 = (/ AOn(nu)%x, AOn(nu)%y, AOn(nu)%z /)
        do a1 = 1, nlist(io)
          mu = alist(a1,io)
          r1 = (/ AOn(mu)%x, AOn(mu)%y, AOn(mu)%z /)
          do ip = 1, size(AOn(mu)%exponent)
            e1 = AOn(mu)%exponent(ip)
            do jp = 1, size(AOn(nu)%exponent)
              e2 = AOn(nu)%exponent(jp)
              call clifford_pair(3, e1, e2, r1, r2, pe, cen)
              pref = C(mu,io)*AOn(mu)%coefficient(ip)*AOn(nu)%coefficient(jp)
              do dd = 1, 3
                coef(dd) = Lv(dd)*dexp(-2.d0*(e1+e2-pe(dd))/av(dd)**2)/av(dd)**(lq(dd,mu)+lq(dd,nu))
              end do
              if (abs(pref*coef(1)*coef(2)*coef(3)) < XB_tthr) cycle
              nT = nT + 1
              do dd = 1, 3
                Mm = lq(dd,mu) + lq(dd,nu)
                A = 2.d0*pe(dd)/av(dd)**2
                call fourier_sin_product(lq(dd,mu), lq(dd,nu), av(dd)*(cen(dd)-r1(dd)), av(dd)*(cen(dd)-r2(dd)), dco(-Mm:Mm))
                nb = min(40000, max(Mi(dd)+Mm+2, Mm + 40 + int(10.d0*dsqrt(A))))
                call bessel_scaled_vector_fast(nb, A, Iv(0:nb))
                do m = merge(0, -Mi(dd), dd == 3), Mi(dd)
                  cv = (0.d0,0.d0)
                  do k = -Mm, Mm
                    cv = cv + dco(k)*Iv(abs(m-k))
                  end do
                  cv = coef(dd)*cv*exp(dcmplx(0.d0, -dble(m)*av(dd)*cen(dd)))
                  if (dd == 1) F1(m,nT) = pref*cv
                  if (dd == 2) F2(m,nT) = cv
                  if (dd == 3) F3(m,nT) = cv
                end do
              end do
            end do
          end do
        end do
        if (nT == 0) then
          H = (0.d0,0.d0)
        else
          allocate(B2(-Mi(2):Mi(2),nT))
          do m3 = 0, Mi(3)
            do k = 1, nT
              B2(:,k) = F2(:,k)*F3(m3,k)
            end do
            call zgemm('N', 'T', 2*Mi(1)+1, 2*Mi(2)+1, nT, (1.d0,0.d0), F1, 2*Mi(1)+1, B2, 2*Mi(2)+1, &
                       (0.d0,0.d0), H(:,:,m3+1), 2*Mi(1)+1)
          end do
          deallocate(B2)
        end if
        !$omp atomic
        XB_nterm = XB_nterm + dble(nT)
        deallocate(F1, F2, F3, Iv)
      end subroutine half_transform

      complex(8) function axis_val(grid, Mv, k, m)       ! grid value at m e_k
        integer, intent(in) :: Mv(3), k, m
        complex(8), intent(in) :: grid(-Mv(1):Mv(1),-Mv(2):Mv(2),0:Mv(3))
        if (k == 1) axis_val = grid(m,0,0)
        if (k == 2) axis_val = grid(0,m,0)
        if (k == 3) axis_val = grid(0,0,m)
      end function axis_val

      ! rho_ij(g) on |m_d| <= M_d ; separable sum filled with zgemm per m3 slice
      subroutine rho_pair(io, jo, Mv, grid)
        integer, intent(in) :: io, jo, Mv(3)
        complex(8), intent(out) :: grid(-Mv(1):Mv(1),-Mv(2):Mv(2),0:Mv(3))    ! m3 >= 0 (rho(-g) = rho(g)*)
        complex(8), allocatable :: F1(:,:), F2(:,:), F3(:,:), B2(:,:)
        double precision, allocatable :: Iv(:)
        complex(8) :: dco(-8:8), cv
        double precision :: r1(3), r2(3), pe(3), cen(3), e1, e2, A, coef(3), pref, w
        integer :: a1, a2, mu, nu, ip, jp, dd, Mm, k, m, nT, nTmax, nb, m3
        nTmax = 0
        do a1 = 1, nlist(io)
          do a2 = 1, nlist(jo)
            nTmax = nTmax + size(AOn(alist(a1,io))%exponent)*size(AOn(alist(a2,jo))%exponent)
          end do
        end do
        allocate(F1(-Mv(1):Mv(1),nTmax), F2(-Mv(2):Mv(2),nTmax), F3(-Mv(3):Mv(3),nTmax), Iv(0:40000))
        nT = 0
        do a1 = 1, nlist(io)
          mu = alist(a1,io)
          r1 = (/ AOn(mu)%x, AOn(mu)%y, AOn(mu)%z /)
          do a2 = 1, nlist(jo)
            nu = alist(a2,jo)
            w = C(mu,io)*C(nu,jo)
            r2 = (/ AOn(nu)%x, AOn(nu)%y, AOn(nu)%z /)
            do ip = 1, size(AOn(mu)%exponent)
              e1 = AOn(mu)%exponent(ip)
              do jp = 1, size(AOn(nu)%exponent)
                e2 = AOn(nu)%exponent(jp)
                call clifford_pair(3, e1, e2, r1, r2, pe, cen)
                pref = w*AOn(mu)%coefficient(ip)*AOn(nu)%coefficient(jp)
                do dd = 1, 3
                  coef(dd) = Lv(dd)*dexp(-2.d0*(e1+e2-pe(dd))/av(dd)**2)/av(dd)**(lq(dd,mu)+lq(dd,nu))
                end do
                if (abs(pref*coef(1)*coef(2)*coef(3)) < XB_tthr) cycle     ! term screening (|rho| ~ 1)
                nT = nT + 1
                do dd = 1, 3
                  Mm = lq(dd,mu) + lq(dd,nu)
                  A = 2.d0*pe(dd)/av(dd)**2
                  call fourier_sin_product(lq(dd,mu), lq(dd,nu), av(dd)*(cen(dd)-r1(dd)), av(dd)*(cen(dd)-r2(dd)), dco(-Mm:Mm))
                  nb = min(40000, max(Mv(dd)+Mm+2, Mm + 40 + int(10.d0*dsqrt(A))))
                  call bessel_scaled_vector_fast(nb, A, Iv(0:nb))
                  do m = -Mv(dd), Mv(dd)
                    cv = (0.d0,0.d0)
                    do k = -Mm, Mm
                      cv = cv + dco(k)*Iv(abs(m-k))
                    end do
                    cv = coef(dd)*cv*exp(dcmplx(0.d0, -dble(m)*av(dd)*cen(dd)))
                    if (dd == 1) F1(m,nT) = pref*cv
                    if (dd == 2) F2(m,nT) = cv
                    if (dd == 3) F3(m,nT) = cv
                  end do
                end do
              end do
            end do
          end do
        end do
        allocate(B2(-Mv(2):Mv(2),max(nT,1)))
        do m3 = 0, Mv(3)
          do k = 1, nT
            B2(:,k) = F2(:,k)*F3(m3,k)
          end do
          if (nT > 0) then
            call zgemm('N', 'T', 2*Mv(1)+1, 2*Mv(2)+1, nT, (1.d0,0.d0), F1, 2*Mv(1)+1, B2, 2*Mv(2)+1, &
                       (0.d0,0.d0), grid(:,:,m3), 2*Mv(1)+1)
          else
            grid(:,:,m3) = (0.d0,0.d0)
          end if
        end do
        deallocate(F1, F2, F3, B2, Iv)
        !$omp atomic
        XB_nterm = XB_nterm + dble(nT)
      end subroutine rho_pair

      ! e_self (i = j): (1/V) sum' (|rho|^2 - Gaussian) D ; e_pair(1:3): (1/V) sum' |rho|^2 D, g_cut 3, 4, 5
      subroutine pair_energy(io, jo, Mp, gij, grid, Mv, Dt, e_self, e_pair)
        integer, intent(in) :: io, jo, Mp(3), Mv(3)
        double precision, intent(in) :: gij
        complex(8), intent(in) :: grid(-Mp(1):Mp(1),-Mp(2):Mp(2),0:Mp(3))
        double precision, intent(in) :: Dt(0:Mv(1),0:Mv(2),0:Mv(3))
        double precision, intent(out) :: e_self(3), e_pair(3)
        integer :: m1, m2, m3
        double precision :: g2, r2, gau, Vs, wm
        Vs = Lv(1)*Lv(2)*Lv(3)
        e_self = 0.d0 ; e_pair = 0.d0
        do m3 = 0, Mp(3)
          wm = 2.d0 ; if (m3 == 0) wm = 1.d0          ! m3 < 0 by symmetry
          do m2 = -Mp(2), Mp(2)
            do m1 = -Mp(1), Mp(1)
              if (m1 == 0 .and. m2 == 0 .and. m3 == 0) cycle
              if (m3 == 0 .and. (m2 < 0 .or. (m2 == 0 .and. m1 < 0))) cycle   ! half of the m3 = 0 plane, weight 2
              if (m3 == 0) wm = 2.d0
              g2 = (m1*av(1))**2 + (m2*av(2))**2 + (m3*av(3))**2
              if (g2 > gij**2) cycle
              if (io == jo) then
                gau = (abs(z(io,io,1))**2)**(m1*m1)*(abs(z(io,io,2))**2)**(m2*m2)*(abs(z(io,io,3))**2)**(m3*m3)
                r2 = wm*(abs(grid(m1,m2,m3))**2 - gau)*Dt(abs(m1),abs(m2),abs(m3))
                if (g2 <= (0.6d0*gij)**2) e_self(1) = e_self(1) + r2
                if (g2 <= (0.8d0*gij)**2) e_self(2) = e_self(2) + r2
                e_self(3) = e_self(3) + r2
              else
                r2 = wm*abs(grid(m1,m2,m3))**2*Dt(abs(m1),abs(m2),abs(m3))
                if (g2 <= (0.6d0*gij)**2) e_pair(1) = e_pair(1) + r2
                if (g2 <= (0.8d0*gij)**2) e_pair(2) = e_pair(2) + r2
                e_pair(3) = e_pair(3) + r2
              end if
            end do
          end do
        end do
        e_self = e_self/Vs ; e_pair = e_pair/Vs
      end subroutine pair_energy

      ! D(g) = v_L(g) - 4 pi/g^2 = (2/sqrt pi) int dt [ prod_d L_d I~_{m_d}(z_d) - pi^3/2 t^-3 e^{-g^2/4t^2} ]
      subroutine kernel_D(Mv, Dt)
        integer, intent(in) :: Mv(3)
        double precision, intent(out) :: Dt(0:Mv(1),0:Mv(2),0:Mv(3))
        double precision, allocatable :: Ivt(:), fac(:,:), Dl(:,:,:)   ! [Nov 2026] Dl: per-thread sum on the heap
        double precision :: smin, smax, dss, t, wt, zz, g2, P0
        integer :: ns, is, dd, nb, m1, m2, m3
        smin = dlog(1.d-2*minval(av)) ; smax = dlog(40.d0*XB_gcut)
        dss = 0.01d0 ; ns = int((smax-smin)/dss) + 1
        Dt = 0.d0
        !$omp parallel default(shared) private(is,t,wt,dd,zz,nb,Ivt,fac,m1,m2,m3,g2,P0,Dl)
        allocate(fac(0:maxval(Mv),3), Dl(0:Mv(1),0:Mv(2),0:Mv(3))) ; Dl = 0.d0
        ! [Nov 2026] was reduction(+:Dt): the private copy lived on the thread stack (SIGBUS risk on macOS)
        !$omp do schedule(dynamic,4)
        do is = 0, ns
          t = dexp(smin + dble(is)*dss)
          wt = t*dss ; if (is == 0 .or. is == ns) wt = 0.5d0*wt
          do dd = 1, 3
            zz = 2.d0*t*t/av(dd)**2
            nb = Mv(dd) + 40 + int(10.d0*dsqrt(zz))
            allocate(Ivt(0:nb))
            call bessel_scaled_vector_fast(nb, zz, Ivt)
            fac(0:Mv(dd),dd) = Lv(dd)*Ivt(0:Mv(dd))
            deallocate(Ivt)
          end do
          do m3 = 0, Mv(3)
            do m2 = 0, Mv(2)
              do m1 = 0, Mv(1)
                if (m1 == 0 .and. m2 == 0 .and. m3 == 0) cycle
                g2 = (m1*av(1))**2 + (m2*av(2))**2 + (m3*av(3))**2
                if (g2 > XB_gcut**2*1.01d0) cycle
                P0 = pi**1.5d0/t**3*dexp(-g2/(4.d0*t*t))
                Dl(m1,m2,m3) = Dl(m1,m2,m3) + wt*(fac(m1,1)*fac(m2,2)*fac(m3,3) - P0)
              end do
            end do
          end do
        end do
        !$omp end do
        !$omp critical
        Dt = Dt + Dl
        !$omp end critical
        deallocate(fac, Dl)
        !$omp end parallel
        Dt = 2.d0/dsqrt(pi)*Dt
      end subroutine kernel_D

      ! kappa_k = - lim [ sum'_g (g_k/g)^2 e^{-eps g^2} - V/(8 (pi eps)^1.5)/3 ] (sum_k kappa_k = 1; cubic: 1/3)
      subroutine lattice_kappa(kap)
        double precision, intent(out) :: kap(3)
        double precision :: eps, Vs, g(3), g2, w, integ
        integer :: Mk(3), m1, m2, m3
        Vs = Lv(1)*Lv(2)*Lv(3)
        eps = 0.025d0/maxval(av)**2
        do m1 = 1, 3
          Mk(m1) = int(7.d0/(dsqrt(eps)*av(m1))) + 1
        end do
        kap = 0.d0
        !$omp parallel do default(shared) private(m1,m2,m3,g,g2,w) reduction(+:kap) collapse(2)
        do m3 = -Mk(3), Mk(3)
          do m2 = -Mk(2), Mk(2)
            do m1 = -Mk(1), Mk(1)
              if (m1 == 0 .and. m2 == 0 .and. m3 == 0) cycle
              g = (/ m1*av(1), m2*av(2), m3*av(3) /)
              g2 = sum(g**2)
              w = dexp(-eps*g2)
              kap = kap + w*g**2/g2
            end do
          end do
        end do
        !$omp end parallel do
        integ = Vs/(8.d0*(pi*eps)**1.5d0)/3.d0
        kap = -(kap - integ)
      end subroutine lattice_kappa

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
        allocate(Iv(0:40000))
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
                  nbes = min(40000, Mm + 40 + int(10.d0*dsqrt(A)))  ! backward recurrence: start well above sqrt(A)
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



!=====================================================================!
!  [Nov 2026] exchange refinement, option A: model for the pair       !
!  terms delta(ij|ji), i /= j (neglected before).                     !
!                                                                     !
!  rho_ij = phi_i phi_j is neutral; modelled as a dipolar Gaussian    !
!     rho_ij(r) = - d . grad g(r - c),  g normalized Gaussian with    !
!     variances sp2(k) (pair width), d = <i| r |j> (transition dipole)!
!  Then (integration by parts)                                        !
!     delta(ij|ji) = int int rho rho dv = - d^T H d ,                 !
!     H = Hessian at R = 0 of W(R) = < dv(r + R) >, r ~ N(0, 2 sp2)   !
!  With diagonal variances H is diagonal:                             !
!     delta(ij|ji) = - sum_k d_k^2 H_kk                               !
!  dv = sum_{per d} a_d^2 x_d^4/(24 r^3)                              !
!     = sum_d (a_d^2/24)(4/sqrt pi) int t^2 x_d^4 exp(-t^2 r^2) dt    !
!  Gaussian average of (x_d + R_d)^4 exp(-t^2 |r+R|^2), per direction !
!  k (s2 = 2 sp2(k), q = 1 + 2 t^2 s2):                               !
!     B_k = q^-1/2 exp(-t^2 R_k^2 / q)                                !
!     M4  = mu^4 + 6 mu^2 v + 3 v^2,  mu = R_d/q, v = s2/q            !
!  H_kk by central differences (h = 1e-2 sqrt(min s2)).               !
!  Sign: H > 0 (W ~ <|r+R|> is convex) -> delta < 0 -> the pair terms !
!  RAISE dE_X (dE_X = - sum_i delta_ii - sum_{i/=j} delta_ij).        !
!=====================================================================!
double precision function fs_pair_model(sp2, d2, av, nper)

      implicit none
      double precision, intent(in) :: sp2(3), d2(3), av(3)
      integer, intent(in)          :: nper
      double precision :: h, W0, Wp, Wm, R(3)
      integer :: k

      h = 1.d-2*dsqrt(2.d0*minval(sp2))
      R = 0.d0
      W0 = Wfun(R)
      fs_pair_model = 0.d0
      do k = 1, 3
        if (d2(k) <= 0.d0) cycle
        R = 0.d0 ; R(k) =  h ; Wp = Wfun(R)
        R(k) = -h            ; Wm = Wfun(R)
        fs_pair_model = fs_pair_model - d2(k)*(Wp + Wm - 2.d0*W0)/(h*h)
      end do

contains

      double precision function Wfun(Rv)
        double precision, intent(in) :: Rv(3)
        integer, parameter :: NQ = 2000
        double precision, parameter :: pi = 3.14159265358979323846d0
        double precision :: s2(3), tau, u, t, f, wq, q(3), Bp, mu, v, ssum, fd
        integer :: iq, d, kk
        s2 = 2.d0*sp2
        tau = 1.d0/dsqrt(sum(s2)/3.d0)
        ssum = 0.d0
        do iq = 1, NQ-1
          u = dble(iq)/dble(NQ)
          t = tau*u/(1.d0-u)
          q = 1.d0 + 2.d0*t*t*s2
          Bp = 1.d0
          do kk = 1, 3
            Bp = Bp*dexp(-t*t*Rv(kk)**2/q(kk))/dsqrt(q(kk))
          end do
          f = 0.d0
          do d = 1, nper
            mu = Rv(d)/q(d) ; v = s2(d)/q(d)
            f = f + av(d)**2/24.d0*(mu**4 + 6.d0*mu*mu*v + 3.d0*v*v)
          end do
          fd = 4.d0/dsqrt(pi)*t*t*Bp*f*tau/(1.d0-u)**2
          wq = 2.d0 ; if (mod(iq,2) == 1) wq = 4.d0
          ssum = ssum + wq*fd
        end do
        Wfun = ssum/(3.d0*dble(NQ))
      end function Wfun

end function fs_pair_model
