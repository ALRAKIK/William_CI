!=====================================================================!
!                                                                     !
!   [Nov 2026] keyword fs_correction for Tori1D and Tori2D (SAO RHF)  !
!                                                                     !
!   Same theory as fs_correction_3D (doc/fs_correction_theory.pdf),   !
!   with ndp periodic directions (1 = x ; 2 = x,y) and no = 3 - ndp   !
!   open directions (ordinary Gaussians).                             !
!                                                                     !
!   Kernel: d^2 = sum_{periodic} (2/a^2)(1-cos a x) + sum_{open} x^2  !
!     1/d - 1/r = sum_{periodic i} a_i^2 x_i^4/(24 r^3) + O(a^4)      !
!                                                                     !
!   Electrostatics, mixed representation (G periodic, k open, g^2 =   !
!   G^2 + k^2), per unit cell (Omega_p = length / area of the cell):  !
!     dE_ES = 1/(2 Omega_p) sum_G int d^no k/(2pi)^no |rho(G,k)|^2 w  !
!     w = dv(G,k) = pi sum_{i periodic} a_i^2 [8g_i^4/g^8 - 8g_i^2/g^6 + 1/g^4]
!   (G = 0 included: rho(0,k) ~ k^2 for a neutral cell without a      !
!   transverse dipole, so |rho|^2 dv stays finite at k -> 0).         !
!   All orders: w -> v_exact(G,k) - 4pi/g^2 with                      !
!     v_exact = (2/sqrt pi) int dt prod_per L I~_m(z) prod_open (sqrt(pi)/t) e^{-k^2/4t^2}
!   computed as remainder R = v_exact - 4pi/g^2 - dv for |n| <= NR.   !
!   Open-direction pair FT (u = x - P, p = al + be):                  !
!     int (x-A)^l1 (x-B)^l2 e^{-al(x-A)^2 - be(x-B)^2} e^{-ikx} dx    !
!       = e^{-mu (A-B)^2} e^{-ikP} sum_j c_j M_j(k)                   !
!     M_0 = sqrt(pi/p) e^{-k^2/4p}, M_1 = -ik/(2p) M_0,               !
!     M_{j+1} = (j M_{j-1} - ik M_j)/(2p)                             !
!   1D: k is 2D (open y,z): radial log grid x angular grid; nuclei    !
!   angular averages analytic (J0). 2D: k = +-kappa on a log grid.    !
!   Nuclei beyond the electron cutoffs: analytic angle, separate      !
!   extended grid; G tail extrapolated (c0 + c1/n + c2/n^2).          !
!                                                                     !
!   Exchange: localized orbitals (Jacobi), functional                 !
!     sum_{per} |<e^{-i a x}>_ii|^2 / a^2 + sum_{open} <x>_ii^2       !
!   spreads: periodic -ln|z|^2/a^2, open <x^2> - <x>^2 ; Gaussian     !
!   holes, a^2 term and all orders (open directions identical in      !
!   torus and bulk).                                                  !
!=====================================================================!
subroutine fs_correction_lowd(ndp, nBas, AO, n_atoms, geometry, charge, P, EHF)

      use files,              only : outfile
      use classification_ERI, only : ERI_function
      use unitcell_module,    only : nfuc, nx, ny
      use torus_init,         only : Lx, Ly, ax, ay
      use clifford_tools,     only : init_clifford, orbital_to_lxyz, bessel_scaled_vector_fast, &
                                     fourier_sin_product, clifford_pair
      use omp_lib

      implicit none

      integer           , intent(in) :: ndp, nBas, n_atoms
      type(ERI_function), intent(in) :: AO(nBas)
      double precision  , intent(in) :: geometry(n_atoms,3)
      integer           , intent(in) :: charge(n_atoms)
      double precision  , intent(in) :: P(nBas,nBas)
      double precision  , intent(in) :: EHF

      double precision, parameter :: pi = 3.14159265358979323846d0
      ! [Nov 2026] NG_CAP = 120 shells was too few for long cells (tubes): electrons cut at G = 120*2pi/L
      ! (Be tube, L = 378 bohr: G = 2 bohr^-1, error +10 uEh). Now the cap is the larger of 120 shells and
      ! a physical cutoff G_ES (default 20 bohr^-1, env FS_GES).
      integer, parameter :: NG_CAP = 120, NR_MAX = 16
      double precision :: G_ES, e_in_h, dES_h
      integer :: NG_h, lne, iste
      character(len=32) :: enve
      type(ERI_function), allocatable :: AOn(:)
      integer :: no, Kc(2), lq(3,nBas), ncell, nat_cell, NG, NR, nGt, Nk, Nang, Nkv, mu, nu, ip, jp, d, k, j, l
      integer :: iG, n1, n2, q, ja, kv, qlim, nl(2), Mm, m, nbes, nthr, ncut(5)
      double precision :: Lv(2), av(2), cell(2), Omega_p, shift(3), pmax, kel, kmin, ds, Rmax, w0
      double precision, allocatable :: kap(:), wk(:), kvec(:,:), Wkv(:), Gv(:,:), Rg(:,:)
      integer, allocatable :: nidx(:,:)
      complex(8), allocatable :: rho(:,:), cm(:,:), fo(:)
      type(ERI_function) :: dum
      double precision :: r1(3), r2(3), pe(2), cen(2), ex1, ex2, pref, coef(2), A, g2, dv, e_in, e_ne, e_ee, e_nn_in
      double precision :: e_nn_out, e_nn_tail, e_nn_cut(5), r_nn, r_ne, r_ee, nel0, ztot, scal
      double precision, allocatable :: Iv(:)
      complex(8) :: dco(-8:8), zn, ze, ph
      double precision :: dES_a2, dES_all, r_tot
      integer :: nGr
      integer, allocatable :: rIdx(:)

      w0 = omp_get_wtime()
      no = 3 - ndp
      Kc = 1 ; Kc(1) = max(nx,1) ; if (ndp == 2) Kc(2) = max(ny,1)
      Lv = (/ Lx, Ly /) ; av = (/ ax, ay /)
      cell = Lv/dble(Kc)
      Omega_p = cell(1) ; if (ndp == 2) Omega_p = cell(1)*cell(2)
      ncell = Kc(1) ; if (ndp == 2) ncell = Kc(1)*Kc(2)
      nat_cell = n_atoms/ncell

      allocate(AOn(nBas))
      AOn = AO
      if (ndp == 1) then
        call normalize_AO_toroidal_1D(nBas, AOn, .false.)
      else
        call normalize_AO_toroidal_2D(nBas, AOn, .false.)
      end if
      call init_clifford()
      do mu = 1, nBas
        call orbital_to_lxyz(AOn(mu)%orbital, lq(:,mu))
      end do

      ! transverse origin: centroid of the cell's nuclei (phases only)
      shift = 0.d0
      do d = ndp+1, 3
        shift(d) = sum(geometry(1:nat_cell,d))/dble(nat_cell)
      end do
      Rmax = 0.d0
      do j = 1, nat_cell
        Rmax = max(Rmax, dsqrt(sum((geometry(j,ndp+1:3)-shift(ndp+1:3))**2)))
      end do
      ztot = dble(sum(charge(1:nat_cell)))

      ! ---------------- cutoffs and k grid ---------------- !
      pmax = 0.d0
      do mu = 1, nfuc
        do ip = 1, size(AOn(mu)%exponent)
          pmax = max(pmax, 2.d0*AOn(mu)%exponent(ip))
        end do
      end do
      G_ES = 20.d0                                                        ! [Nov 2026]
      call get_environment_variable('FS_GES', enve, lne, iste)
      if (iste == 0 .and. lne > 0) read(enve,*) G_ES
      NG = 2
      do d = 1, ndp
        NG = max(NG, min(max(NG_CAP, ceiling(G_ES*cell(d)/(2.d0*pi))), &
                         ceiling(10.d0*dsqrt(pmax)*cell(d)/(2.d0*pi)) + 4))   ! [Nov 2026]
      end do
      NG_h = NG/2                                                         ! [Nov 2026] convergence check
      NR = min(NG, NR_MAX)
      kel  = 10.d0*dsqrt(pmax) + 2.d0          ! electrons: |k| <= kel
      kmin = 1.d-2*min(1.d0, minval(av(1:ndp)))
      ds   = 0.04d0
      Nk   = int(dlog(kel/kmin)/ds) + 1
      allocate(kap(Nk), wk(Nk))
      do q = 1, Nk                                         ! log grid, trapezoid in s = ln kappa
        kap(q) = kmin*dexp(dble(q-1)*ds)
        wk(q)  = ds*kap(q) ; if (q == 1 .or. q == Nk) wk(q) = 0.5d0*wk(q)
      end do
      if (no == 2) then
        Nang = min(256, max(16, 2*ceiling(kel*Rmax) + 16))
      else
        Nang = 2                                           ! k = +kappa, -kappa
      end if
      Nkv = Nk*Nang
      allocate(kvec(no,Nkv), Wkv(Nkv))
      do q = 1, Nk
        do ja = 1, Nang
          kv = (q-1)*Nang + ja
          if (no == 2) then
            kvec(1,kv) = kap(q)*cos(2.d0*pi*dble(ja-1)/dble(Nang))
            kvec(2,kv) = kap(q)*sin(2.d0*pi*dble(ja-1)/dble(Nang))
            Wkv(kv) = wk(q)*kap(q)*(2.d0*pi/dble(Nang))/(4.d0*pi*pi)     ! d^2k/(2pi)^2
          else
            kvec(1,kv) = kap(q)*dble(3-2*ja)                              ! +kappa, -kappa
            Wkv(kv) = wk(q)/(2.d0*pi)                                     ! dk/(2pi)
          end if
        end do
      end do

      ! G list (periodic), flattened
      nGt = 2*NG+1 ; if (ndp == 2) nGt = nGt*nGt
      allocate(nidx(2,nGt), Gv(2,nGt))
      iG = 0
      do n2 = merge(-NG, 0, ndp == 2), merge(NG, 0, ndp == 2)
        do n1 = -NG, NG
          iG = iG + 1
          nidx(:,iG) = (/ n1, n2 /)
          Gv(:,iG) = 0.d0
          Gv(1,iG) = 2.d0*pi*dble(n1)/cell(1)
          if (ndp == 2) Gv(2,iG) = 2.d0*pi*dble(n2)/cell(2)
        end do
      end do

      ! ---------------- electrons rho(G,k) ---------------- !
      allocate(rho(nGt,Nkv)) ; rho = (0.d0,0.d0)
      nel0 = 0.d0
      nthr = max(1, min(omp_get_max_threads(), int(2.d9/(16.d0*dble(nGt)*dble(Nkv)))))
      !$omp parallel num_threads(nthr) default(shared) private(mu,nu,ip,jp,d,k,m,Mm,r1,r2,ex1,ex2,pe,cen,pref,coef,A,dco,Iv,cm,fo,nl,qlim,iG,kv,q,ph,nbes)
      allocate(Iv(0:40000), cm(-NG:NG,2), fo(Nkv))
      !$omp do schedule(dynamic,4) reduction(+:rho,nel0)
      do mu = 1, nfuc
        do nu = 1, nBas
          if (abs(P(mu,nu)) < 1.d-12) cycle
          r1 = (/ AOn(mu)%x, AOn(mu)%y, AOn(mu)%z /) - shift
          r2 = (/ AOn(nu)%x, AOn(nu)%y, AOn(nu)%z /) - shift
          do ip = 1, size(AOn(mu)%exponent)
            ex1 = AOn(mu)%exponent(ip)
            do jp = 1, size(AOn(nu)%exponent)
              ex2 = AOn(nu)%exponent(jp)
              call clifford_pair(ndp, ex1, ex2, r1, r2, pe(1:ndp), cen(1:ndp))
              pref = P(mu,nu)*AOn(mu)%coefficient(ip)*AOn(nu)%coefficient(jp)
              do d = 1, ndp
                coef(d) = Lv(d)*dexp(-2.d0*(ex1+ex2-pe(d))/av(d)**2)/av(d)**(lq(d,mu)+lq(d,nu))
              end do
              if (abs(pref*product(coef(1:ndp))*open_pref(ex1,ex2,r1,r2)) < 1.d-14) cycle
              nl = 0
              do d = 1, ndp
                Mm = lq(d,mu) + lq(d,nu)
                nl(d) = min(NG, ceiling(10.d0*dsqrt(ex1+ex2)*cell(d)/(2.d0*pi)) + Mm)
                A = 2.d0*pe(d)/av(d)**2
                call fourier_sin_product(lq(d,mu), lq(d,nu), av(d)*(cen(d)-r1(d)), av(d)*(cen(d)-r2(d)), dco(-Mm:Mm))
                nbes = min(40000, max(nl(d)*Kc(d)+12+Mm, Mm + 40 + int(10.d0*dsqrt(A))))
                call bessel_scaled_vector_fast(nbes, A, Iv(0:nbes))
                do iG = -nl(d), nl(d)
                  m = iG*Kc(d)
                  cm(iG,d) = (0.d0,0.d0)
                  do k = -Mm, Mm
                    cm(iG,d) = cm(iG,d) + dco(k)*Iv(abs(m-k))
                  end do
                  cm(iG,d) = coef(d)*cm(iG,d)*exp(dcmplx(0.d0, -dble(m)*av(d)*cen(d)))
                end do
              end do
              ! open directions on the k nodes, up to this pair's own cutoff
              qlim = 0
              do q = 1, Nk
                if (kap(q) <= 10.d0*dsqrt(ex1+ex2) + 2.d0) qlim = q
              end do
              do kv = 1, qlim*Nang
                fo(kv) = open_ft(mu, nu, ex1, ex2, r1, r2, kvec(:,kv))
              end do
              ph = pref*cm(0,1) ; if (ndp == 2) ph = ph*cm(0,2)
              nel0 = nel0 + dble(ph*open_ft(mu, nu, ex1, ex2, r1, r2, (/ 0.d0, 0.d0 /)))
              do iG = 1, nGt
                if (abs(nidx(1,iG)) > nl(1)) cycle
                if (ndp == 2) then
                  if (abs(nidx(2,iG)) > nl(2)) cycle
                  ph = pref*cm(nidx(1,iG),1)*cm(nidx(2,iG),2)
                else
                  ph = pref*cm(nidx(1,iG),1)
                end if
                do kv = 1, qlim*Nang
                  rho(iG,kv) = rho(iG,kv) + ph*fo(kv)
                end do
              end do
            end do
          end do
        end do
      end do
      !$omp end do
      deallocate(Iv, cm, fo)
      !$omp end parallel

      ! electrons per cell (k = 0, G = 0) and exact neutrality
      scal = ztot/nel0
      rho = rho*scal

      ! ---------------- ES sums ---------------- !
      ! inside the electron cutoffs: total rho = nuc - el (nuclei on the same k nodes)
      e_ne = 0.d0 ; e_ee = 0.d0 ; e_nn_in = 0.d0 ; e_in = 0.d0 ; e_in_h = 0.d0
      !$omp parallel do default(shared) private(iG,kv,q,g2,dv,zn,ze) reduction(+:e_ne,e_ee,e_nn_in,e_in,e_in_h)
      do iG = 1, nGt
        do kv = 1, Nkv
          q = (kv-1)/Nang + 1
          g2 = sum(Gv(1:ndp,iG)**2) + kap(q)**2
          dv = dv_a2(Gv(:,iG), g2)
          zn = nuc_rho(Gv(:,iG), kvec(:,kv))
          ze = rho(iG,kv)
          e_nn_in = e_nn_in + Wkv(kv)*dv*dble(zn*conjg(zn))
          e_ne = e_ne + Wkv(kv)*dv*(-2.d0*dble(zn*conjg(ze)))
          e_ee = e_ee + Wkv(kv)*dv*dble(ze*conjg(ze))
          e_in = e_in + Wkv(kv)*dv*abs(zn-ze)**2        ! total charge: no cancellation at G = 0, k -> 0
          if (abs(nidx(1,iG)) <= NG_h .and. abs(nidx(2,iG)) <= NG_h) &
            e_in_h = e_in_h + Wkv(kv)*dv*abs(zn-ze)**2  ! [Nov 2026] same with half the G cutoff
        end do
      end do
      !$omp end parallel do
      ! nuclei: kappa > kel for |n| <= NG, and all kappa for |n| > NG (G tail extrapolated)
      e_nn_out = nuc_shell(NG, kel)
      ncut = max(NG+8, 20)*(/ 3, 4, 5, 6, 7 /)/3
      do k = 1, 5
        e_nn_cut(k) = nuc_full(ncut(k))
      end do
      e_nn_tail = extrap5(e_nn_cut, ncut) - nuc_full(NG)
      dES_a2 = (e_in + e_nn_out + e_nn_tail)/(2.d0*Omega_p)
      ! [Nov 2026] the same with |n| <= NG/2 (nuclei beyond NG/2 from the same tail): convergence indicator
      dES_h = (e_in_h + nuc_shell(NG_h, kel) + extrap5(e_nn_cut, ncut) - nuc_full(NG_h))/(2.d0*Omega_p)

      ! exact-kernel remainder for |n| <= NR, kappa <= kel (nuclei + electrons)
      nGr = 0
      do iG = 1, nGt
        if (abs(nidx(1,iG)) <= NR .and. abs(nidx(2,iG)) <= NR) nGr = nGr + 1
      end do
      allocate(Rg(nGr,Nk), rIdx(nGr)) ; Rg = 0.d0
      nGr = 0
      do iG = 1, nGt
        if (abs(nidx(1,iG)) <= NR .and. abs(nidx(2,iG)) <= NR) then
          nGr = nGr + 1 ; rIdx(nGr) = iG
        end if
      end do
      call lowd_remainder()
      r_tot = 0.d0
      do j = 1, nGr
        iG = rIdx(j)
        do kv = 1, Nkv
          q = (kv-1)/Nang + 1
          zn = nuc_rho(Gv(:,iG), kvec(:,kv))
          r_tot = r_tot + Wkv(kv)*Rg(j,q)*abs(zn-rho(iG,kv))**2
        end do
      end do
      r_tot = r_tot/(2.d0*Omega_p)
      dES_all = dES_a2 + r_tot

      write(outfile,'(a)') ''
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      write(outfile,'(a,i1,a)') '  Finite-size correction (keyword fs_correction), electrostatics, ', ndp, 'D'
      write(outfile,'(a)') '  dE_ES = 1/(2 Omega) sum_G int dk |rho(G,k)|^2 dv(G,k)'
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      write(outfile,'(a,2i5)')        '  cells per periodic direction    : ', Kc(1:ndp)
      write(outfile,'(a,f14.8,a,es9.1)') '  electrons per cell              : ', nel0, '   (rescaled by ', scal-1.d0
      write(outfile,'(a,i5,a,i5,a,f8.2,a,i6,a,i4)') '  cutoffs: |n| <=', NG, '  remainder |n| <=', NR, '  |k| <=', kel, &
                                      '  k nodes', Nk, ' x', Nang
      write(outfile,'(a,f18.10)')     '  inside cutoffs (nuc + el)       : ', e_in/(2.d0*Omega_p)
      write(outfile,'(a,f18.10)')     '  nuclei outside (+ G tail)       : ', (e_nn_out+e_nn_tail)/(2.d0*Omega_p)
      write(outfile,'(a,f8.2,a,f18.10)') '  G cutoff check: G_ES =', G_ES, ' ; a^2 term with |n| <= NG/2 :', dES_h
      write(outfile,'(a,f18.10,a)')   '  dE_ES (a^2 term)                : ', dES_a2, ' Eh/cell'
      write(outfile,'(a,f18.10)')     '  higher orders of the kernel     : ', r_tot
      write(outfile,'(a,f18.10,a)')   '  dE_ES (all orders)              : ', dES_all, ' Eh/cell'
      write(outfile,'(a,f10.2,a)')    '  wall time                       : ', omp_get_wtime()-w0, ' s'
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      flush(outfile)
      deallocate(rho, Rg, rIdx, kap, wk, kvec, Wkv, nidx, Gv)

      call fs_exchange_lowd(ndp, nBas, AOn, lq, P, Lv, av, ncell, EHF, dES_a2, dES_all)
      deallocate(AOn)

contains

      ! e^{-mu (A-B)^2} over the open directions (screening)
      double precision function open_pref(e1, e2, ra, rb)
        double precision, intent(in) :: e1, e2, ra(3), rb(3)
        open_pref = dexp(-e1*e2/(e1+e2)*sum((ra(ndp+1:3)-rb(ndp+1:3))**2))
      end function open_pref

      ! product over the open directions of the 1D pair FT at k
      complex(8) function open_ft(i1, i2, e1, e2, ra, rb, kk)
        integer, intent(in) :: i1, i2
        double precision, intent(in) :: e1, e2, ra(3), rb(3), kk(:)
        integer :: dd, o
        open_ft = (1.d0,0.d0)
        do o = 1, no
          dd = ndp + o
          open_ft = open_ft*ft_1d(lq(dd,i1), lq(dd,i2), e1, e2, ra(dd), rb(dd), kk(o))
        end do
      end function open_ft

      ! int (x-A)^l1 (x-B)^l2 exp(-e1(x-A)^2 - e2(x-B)^2) exp(-ikx) dx
      complex(8) function ft_1d(l1, l2, e1, e2, xa, xb, kk)
        integer, intent(in) :: l1, l2
        double precision, intent(in) :: e1, e2, xa, xb, kk
        double precision :: pp, xp, c(0:8)
        complex(8) :: Mj(0:9)
        integer :: jj
        pp = e1 + e2 ; xp = (e1*xa + e2*xb)/pp
        call poly_shift(l1, l2, xp-xa, xp-xb, 0, 0.d0, c)
        Mj(0) = dsqrt(pi/pp)*dexp(-kk*kk/(4.d0*pp))
        Mj(1) = dcmplx(0.d0, -kk/(2.d0*pp))*Mj(0)
        do jj = 1, l1+l2-1
          Mj(jj+1) = (dble(jj)*Mj(jj-1) + dcmplx(0.d0,-kk)*Mj(jj))/(2.d0*pp)
        end do
        ft_1d = (0.d0,0.d0)
        do jj = 0, l1+l2
          ft_1d = ft_1d + c(jj)*Mj(jj)
        end do
        ft_1d = ft_1d*dexp(-e1*e2/pp*(xa-xb)**2)*exp(dcmplx(0.d0, -kk*xp))
      end function ft_1d



      complex(8) function nuc_rho(G, kk)
        double precision, intent(in) :: G(2), kk(:)
        integer :: jj
        double precision :: ph1
        nuc_rho = (0.d0,0.d0)
        do jj = 1, nat_cell
          ph1 = sum(G(1:ndp)*geometry(jj,1:ndp)) + sum(kk(1:no)*(geometry(jj,ndp+1:3)-shift(ndp+1:3)))
          nuc_rho = nuc_rho + dble(charge(jj))*exp(dcmplx(0.d0, -ph1))
        end do
      end function nuc_rho

      ! angle-averaged |rho_nuc(G,k)|^2 at |k| = kk (1D: J0, 2D: cos)
      double precision function nuc_rho2_avg(G, kk)
        double precision, intent(in) :: G(2), kk
        integer :: j1, j2
        double precision :: dper, dtr
        nuc_rho2_avg = 0.d0
        do j1 = 1, nat_cell
          do j2 = 1, nat_cell
            dper = sum(G(1:ndp)*(geometry(j1,1:ndp)-geometry(j2,1:ndp)))
            if (no == 2) then
              dtr = dsqrt(sum((geometry(j1,2:3)-geometry(j2,2:3))**2))
              nuc_rho2_avg = nuc_rho2_avg + dble(charge(j1)*charge(j2))*cos(dper)*bessel_j0(kk*dtr)
            else
              dtr = geometry(j1,3) - geometry(j2,3)
              nuc_rho2_avg = nuc_rho2_avg + dble(charge(j1)*charge(j2))*cos(dper)*cos(kk*dtr)
            end if
          end do
        end do
      end function nuc_rho2_avg

      ! int over all k with |k| measure, weight dv, of <|rho_nuc|^2>, for kappa in (k0, inf)
      double precision function nuc_kint(G, k0)
        double precision, intent(in) :: G(2), k0
        integer :: qq, nq
        double precision :: s0, s1, dss, kk, w, gg2, fac
        s0 = dlog(max(k0, 1.d-6)) ; s1 = dlog(max(k0,1.d0)*1.d4) ; dss = 0.02d0
        nq = int((s1-s0)/dss) + 1
        nuc_kint = 0.d0
        do qq = 0, nq
          kk = dexp(s0 + dble(qq)*dss)
          w = dss*kk ; if (qq == 0 .or. qq == nq) w = 0.5d0*w
          if (no == 2) then
            fac = kk/(2.d0*pi)              ! d^2k/(2pi)^2 = kappa dkappa / (2pi)
          else
            fac = 2.d0/(2.d0*pi)            ! +-k
          end if
          gg2 = sum(G(1:ndp)**2) + kk*kk
          nuc_kint = nuc_kint + w*fac*dv_a2(G, gg2)*nuc_rho2_avg(G, kk)
        end do
      end function nuc_kint

      ! nuclei for |n| <= nc (G /= 0 or kappa > k0 for G = 0) over kappa > kel ... (shell)
      double precision function nuc_shell(nc, k0)
        integer, intent(in) :: nc
        double precision, intent(in) :: k0
        integer :: ii
        nuc_shell = 0.d0
        do ii = 1, nGt
          if (abs(nidx(1,ii)) > nc .or. abs(nidx(2,ii)) > nc) cycle
          nuc_shell = nuc_shell + nuc_kint(Gv(:,ii), k0)
        end do
      end function nuc_shell

      ! nuclei, all kappa, 0 < |n| <= nc (G = 0 excluded: it is in the combined grid)
      double precision function nuc_full(nc)
        integer, intent(in) :: nc
        integer :: a1, a2
        double precision :: G(2)
        nuc_full = 0.d0
        !$omp parallel do default(shared) private(a1,a2,G) reduction(+:nuc_full) collapse(2)
        do a2 = merge(-nc, 0, ndp == 2), merge(nc, 0, ndp == 2)
          do a1 = -nc, nc
            if (a1 == 0 .and. a2 == 0) cycle
            G = (/ 2.d0*pi*dble(a1)/cell(1), 0.d0 /)
            if (ndp == 2) G(2) = 2.d0*pi*dble(a2)/cell(2)
            nuc_full = nuc_full + nuc_kint(G, 1.d-6)
          end do
        end do
        !$omp end parallel do
      end function nuc_full

      double precision function dv_a2(G, gg2)
        double precision, intent(in) :: G(2), gg2
        integer :: ii
        dv_a2 = 0.d0
        do ii = 1, ndp
          dv_a2 = dv_a2 + av(ii)**2*(8.d0*G(ii)**4/gg2**4 - 8.d0*G(ii)**2/gg2**3 + 1.d0/gg2**2)
        end do
        dv_a2 = pi*dv_a2
      end function dv_a2

      double precision function extrap5(v, nn)
        double precision, intent(in) :: v(5)
        integer, intent(in) :: nn(5)
        double precision :: Am(5,3), At(3,3), rhs(3), x(3), f
        integer :: r, i, jj
        do r = 1, 5
          Am(r,:) = (/ 1.d0, 1.d0/dble(nn(r)), 1.d0/dble(nn(r))**2 /)
        end do
        At = matmul(transpose(Am), Am) ; rhs = matmul(transpose(Am), v)
        do i = 1, 3
          do jj = i+1, 3
            f = At(jj,i)/At(i,i) ; At(jj,:) = At(jj,:) - f*At(i,:) ; rhs(jj) = rhs(jj) - f*rhs(i)
          end do
        end do
        do i = 3, 1, -1
          x(i) = (rhs(i) - sum(At(i,i+1:3)*x(i+1:3)))/At(i,i)
        end do
        extrap5 = x(1)
      end function extrap5

      ! R(G,kappa) = (2/sqrt pi) int dt (sqrt(pi)/t)^no e^{-k^2/4t^2}
      !              [ prod_per L I~_m(z) - prod_per (sqrt(pi)/t) e^{-G_d^2/4t^2} (1 + sum eps_d) ]
      subroutine lowd_remainder()
        double precision, allocatable :: Ivt(:), fac(:,:)
        double precision :: smin, smax, dss, t, wt, z, P0, eps, op, gmx
        integer :: ns, is, dd, nb, ii, qq, jr
        gmx = dsqrt(sum((2.d0*pi*dble(NR)/cell(1:ndp))**2) + kel**2)
        smin = dlog(min(kmin, 2.d0*pi/maxval(cell(1:ndp)))/40.d0) ; smax = dlog(20.d0*gmx)
        dss = 0.01d0 ; ns = int((smax-smin)/dss) + 1
        !$omp parallel default(shared) private(is,t,wt,dd,z,nb,Ivt,fac,ii,qq,P0,eps,op,jr)
        allocate(fac(0:NR,2))
        !$omp do schedule(dynamic,4) reduction(+:Rg)
        do is = 0, ns
          t = dexp(smin + dble(is)*dss)
          wt = t*dss ; if (is == 0 .or. is == ns) wt = 0.5d0*wt
          fac = 1.d0
          do dd = 1, ndp
            z = 2.d0*t*t/av(dd)**2
            nb = NR*Kc(dd) + 40 + int(10.d0*dsqrt(z))
            allocate(Ivt(0:nb))
            call bessel_scaled_vector_fast(nb, z, Ivt)
            do ii = 0, NR
              fac(ii,dd) = Lv(dd)*Ivt(ii*Kc(dd))
            end do
            deallocate(Ivt)
          end do
          do jr = 1, nGr
            ii = rIdx(jr)
            P0 = 1.d0 ; eps = 0.d0
            do dd = 1, ndp
              P0 = P0*dsqrt(pi)/t*dexp(-Gv(dd,ii)**2/(4.d0*t*t))
              eps = eps + av(dd)**2/(48.d0*t*t)*(Gv(dd,ii)**4/(4.d0*t**4) - 3.d0*Gv(dd,ii)**2/(t*t) + 3.d0)
            end do
            do qq = 1, Nk
              op = (dsqrt(pi)/t)**no*dexp(-kap(qq)**2/(4.d0*t*t))
              if (op < 1.d-300) cycle
              if (ndp == 2) then
                Rg(jr,qq) = Rg(jr,qq) + wt*op*(fac(abs(nidx(1,ii)),1)*fac(abs(nidx(2,ii)),2) - P0*(1.d0+eps))
              else
                Rg(jr,qq) = Rg(jr,qq) + wt*op*(fac(abs(nidx(1,ii)),1) - P0*(1.d0+eps))
              end if
            end do
          end do
        end do
        !$omp end do
        deallocate(fac)
        !$omp end parallel
        Rg = 2.d0/dsqrt(pi)*Rg
      end subroutine lowd_remainder

end subroutine fs_correction_lowd


! coefficients c_j of (u + da)^l1 (u + db)^l2 (u + dc)^l3 in powers of u (degree <= 8)
subroutine poly_shift(l1, l2, da, db, l3, dc, c)
      implicit none
      integer, intent(in) :: l1, l2, l3
      double precision, intent(in) :: da, db, dc
      double precision, intent(out) :: c(0:8)
      double precision :: t(0:8)
      integer :: i, j
      c = 0.d0 ; c(0) = 1.d0
      do i = 1, l1
        t = 0.d0
        do j = 0, 7
          t(j+1) = t(j+1) + c(j) ; t(j) = t(j) + da*c(j)
        end do
        c = t
      end do
      do i = 1, l2
        t = 0.d0
        do j = 0, 7
          t(j+1) = t(j+1) + c(j) ; t(j) = t(j) + db*c(j)
        end do
        c = t
      end do
      do i = 1, l3
        t = 0.d0
        do j = 0, 7
          t(j+1) = t(j+1) + c(j) ; t(j) = t(j) + dc*c(j)
        end do
        c = t
      end do
end subroutine poly_shift


!=====================================================================!
!  [Nov 2026] exchange part for Tori1D / Tori2D (see fs_exchange_3D)  !
!=====================================================================!
subroutine fs_exchange_lowd(ndp, nBas, AOn, lq, P, Lv, av, ncell, EHF, dE_ES, dE_ES_all)

      use files,              only : outfile
      use classification_ERI, only : ERI_function
      use clifford_tools,     only : bessel_scaled_vector_fast, fourier_sin_product, clifford_pair
      use keywords,           only : c_fs_exact
      use omp_lib

      implicit none

      integer           , intent(in) :: ndp, nBas, lq(3,nBas), ncell
      type(ERI_function), intent(in) :: AOn(nBas)
      double precision  , intent(in) :: P(nBas,nBas), Lv(2), av(2), EHF, dE_ES, dE_ES_all

      double precision, parameter :: pi = 3.14159265358979323846d0
      integer, parameter :: MAXSWEEP = 500
      double precision, allocatable :: S(:,:), Zi(:,:), ev(:), Sh(:,:), Smh(:,:), C(:,:), T1(:,:), T2(:,:), &
                                       sig2(:,:), work(:), Rx(:), x2(:,:,:)
      complex(8), allocatable :: z(:,:,:), zci(:), zcj(:)
      double precision :: wgt(3), ntot, ocmin, ocmax, Om0, Om1, B, Cc, th, cs, sn, gain, gsum, dEx, dExa, s2(3), dEi, w0
      double precision :: dXp, d2p(3), sp2p(3), avp(3)
      double precision, allocatable :: cci(:), ccj(:)
      double precision :: XB_self, XB_pair, XB_lat, XB_spart(3), XB_part(3), XB_chk, XB_kap(3), XB_out, XB_gcut, XB_tthr, XB_nterm, w0B
      integer :: XB_nhome, XB_npair, XB_M(2), XB_Nk, XB_Nang
      integer, allocatable :: nlist(:), alist(:,:)
      integer :: npair
      double precision, external :: fs_pair_model
      complex(8) :: vv, w
      integer :: nocc, lwork, info, i, j, d, isweep, no, k

      w0 = omp_get_wtime()
      no = 3 - ndp

      allocate(S(nBas,nBas))
      call pair_mat(0, 0, S)
      ntot = sum(P*S)
      nocc = nint(0.5d0*ntot)
      allocate(ev(nBas))
      lwork = -1 ; allocate(work(1))
      call dsyev('V', 'U', nBas, S, nBas, ev, work, lwork, info)
      lwork = int(work(1)) ; deallocate(work) ; allocate(work(lwork))
      call dsyev('V', 'U', nBas, S, nBas, ev, work, lwork, info)
      allocate(T1(nBas,nBas), Sh(nBas,nBas), Smh(nBas,nBas))
      do i = 1, nBas
        T1(:,i) = S(:,i)*dsqrt(ev(i))
      end do
      call dgemm('N', 'T', nBas, nBas, nBas, 1.d0, T1, nBas, S, nBas, 0.d0, Sh, nBas)
      do i = 1, nBas
        T1(:,i) = S(:,i)/dsqrt(ev(i))
      end do
      call dgemm('N', 'T', nBas, nBas, nBas, 1.d0, T1, nBas, S, nBas, 0.d0, Smh, nBas)
      call dgemm('N', 'N', nBas, nBas, nBas, 0.5d0, P, nBas, Sh, nBas, 0.d0, T1, nBas)
      call dgemm('N', 'N', nBas, nBas, nBas, 1.d0, Sh, nBas, T1, nBas, 0.d0, S, nBas)
      deallocate(Sh, T1)
      deallocate(work) ; lwork = -1 ; allocate(work(1))
      call dsyev('V', 'U', nBas, S, nBas, ev, work, lwork, info)
      lwork = int(work(1)) ; deallocate(work) ; allocate(work(lwork))
      call dsyev('V', 'U', nBas, S, nBas, ev, work, lwork, info)
      ocmin = ev(nBas-nocc+1)
      ocmax = 0.d0 ; if (nBas > nocc) ocmax = ev(nBas-nocc)
      allocate(C(nBas,nocc))
      call dgemm('N', 'N', nBas, nocc, nBas, 1.d0, Smh, nBas, S(:,nBas-nocc+1:nBas), nBas, 0.d0, C, nBas)
      deallocate(work, ev, Smh)

      ! z(:,:,d): periodic d -> <e^{-i a x}>, weight 1/a^2 ; open d -> <x>, weight 1
      ! x2(:,:,d): open d -> <x^2> (rotated along, not in the functional)
      allocate(z(nocc,nocc,3), x2(nocc,nocc,3), T1(nBas,nocc), T2(nocc,nocc), Zi(nBas,nBas))
      x2 = 0.d0
      do d = 1, 3
        if (d <= ndp) then
          call pair_mat(d, 1, S, Zi)
          wgt(d) = 1.d0/av(d)**2
        else
          call pair_mat(d, 1, S)
          Zi = 0.d0
          wgt(d) = 1.d0
        end if
        call dgemm('N', 'N', nBas, nocc, nBas, 1.d0, S, nBas, C, nBas, 0.d0, T1, nBas)
        call dgemm('T', 'N', nocc, nocc, nBas, 1.d0, C, nBas, T1, nBas, 0.d0, T2, nocc)
        z(:,:,d) = dcmplx(T2, 0.d0)
        call dgemm('N', 'N', nBas, nocc, nBas, 1.d0, Zi, nBas, C, nBas, 0.d0, T1, nBas)
        call dgemm('T', 'N', nocc, nocc, nBas, 1.d0, C, nBas, T1, nBas, 0.d0, T2, nocc)
        z(:,:,d) = z(:,:,d) + dcmplx(0.d0, T2)
        if (d > ndp) then
          call pair_mat(d, 2, S)
          call dgemm('N', 'N', nBas, nocc, nBas, 1.d0, S, nBas, C, nBas, 0.d0, T1, nBas)
          call dgemm('T', 'N', nocc, nocc, nBas, 1.d0, C, nBas, T1, nBas, 0.d0, x2(:,:,d), nocc)
        end if
      end do
      deallocate(S, Zi, T1, T2)            ! [Nov 2026] C kept: localized orbitals for option B

      Om0 = omega_val()
      allocate(zci(nocc), zcj(nocc), cci(nBas), ccj(nBas))
      do isweep = 1, MAXSWEEP
        gsum = 0.d0
        do i = 1, nocc-1
          do j = i+1, nocc
            B = 0.d0 ; Cc = 0.d0
            do d = 1, 3
              vv = 0.5d0*(z(i,i,d) - z(j,j,d)) ; w = z(i,j,d)
              B  = B  + wgt(d)*(abs(vv)**2 - abs(w)**2)
              Cc = Cc + wgt(d)*2.d0*dble(vv*conjg(w))
            end do
            gain = dsqrt(B*B + Cc*Cc) - B
            if (gain < 1.d-14*maxval(wgt)) cycle
            gsum = gsum + gain
            th = 0.25d0*atan2(Cc, B)
            cs = cos(th) ; sn = sin(th)
            do d = 1, 3
              zci = z(:,i,d) ; zcj = z(:,j,d)
              z(:,i,d) =  cs*zci + sn*zcj ; z(:,j,d) = -sn*zci + cs*zcj
              zci = z(i,:,d) ; zcj = z(j,:,d)
              z(i,:,d) =  cs*zci + sn*zcj ; z(j,:,d) = -sn*zci + cs*zcj
              if (d > ndp) then
                zci = dcmplx(x2(:,i,d),0.d0) ; zcj = dcmplx(x2(:,j,d),0.d0)
                x2(:,i,d) = dble( cs*zci + sn*zcj) ; x2(:,j,d) = dble(-sn*zci + cs*zcj)
                zci = dcmplx(x2(i,:,d),0.d0) ; zcj = dcmplx(x2(j,:,d),0.d0)
                x2(i,:,d) = dble( cs*zci + sn*zcj) ; x2(j,:,d) = dble(-sn*zci + cs*zcj)
              end if
            end do
            cci = C(:,i) ; ccj = C(:,j)                 ! [Nov 2026] same rotation on the orbitals
            C(:,i) =  cs*cci + sn*ccj
            C(:,j) = -sn*cci + cs*ccj
          end do
        end do
        if (gsum < 1.d-11*dble(nocc)*maxval(wgt)) exit
      end do
      Om1 = omega_val()

      allocate(sig2(3,nocc), Rx(nocc))
      dEx = 0.d0
      do i = 1, nocc
        do d = 1, 3
          if (d <= ndp) then
            sig2(d,i) = -dlog(abs(z(i,i,d))**2)/av(d)**2
          else
            sig2(d,i) = x2(i,i,d) - dble(z(i,i,d))**2
          end if
        end do
        s2 = 2.d0*sig2(:,i)
        dEi = 0.d0
        do d = 1, ndp
          dEi = dEi + av(d)**2/24.d0*x4r3(s2, d)
        end do
        dEx = dEx - dEi
      end do
      call x_rem(Rx)
      dExa = dEx - sum(Rx)

      ! [Nov 2026] option A: pair terms i /= j (see fs_pair_model in fs_correction.f90)
      !   periodic: |d_k|^2 = |z^k_ij|^2/a_k^2 ; open: d_k = <i|x_k|j>
      dXp = 0.d0 ; npair = 0
      do d = 1, 3
        avp(d) = 0.d0 ; if (d <= ndp) avp(d) = av(d)
      end do
      !$omp parallel do default(shared) private(i,j,d,d2p,sp2p) reduction(+:dXp,npair) schedule(dynamic,4)
      do i = 1, nocc-1
        do j = i+1, nocc
          do d = 1, 3
            if (d <= ndp) then
              d2p(d) = abs(z(i,j,d))**2/av(d)**2
            else
              d2p(d) = dble(z(i,j,d))**2
            end if
            sp2p(d) = 2.d0*sig2(d,i)*sig2(d,j)/(sig2(d,i) + sig2(d,j))
          end do
          if (sum(d2p) < 1.d-10) cycle
          npair = npair + 1
          dXp = dXp - 2.d0*fs_pair_model(sp2p, d2p, avp, ndp)
        end do
      end do
      !$omp end parallel do

      write(outfile,'(a)') ''
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      write(outfile,'(a,i1,a)') '  Finite-size correction (keyword fs_correction), exchange, ', ndp, 'D'
      write(outfile,'(a)') '  localized orbitals, Gaussian holes, self terms'
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      write(outfile,'(a,f14.6,i8)')   '  electrons in supercell, n_occ   : ', ntot, nocc
      write(outfile,'(a,2es12.3)')    '  Loewdin check (occ min, virt max): ', ocmin, ocmax
      write(outfile,'(a,i6)')         '  Jacobi sweeps                   : ', min(isweep, MAXSWEEP)
      write(outfile,'(a,2es16.6)')    '  functional (start/end)          : ', Om0, Om1
      write(outfile,'(a,3f12.6)')     '  sigma_d rms over orbitals (x,y,z):', (dsqrt(sum(sig2(d,:))/nocc), d=1,3)
      write(outfile,'(a,2f12.6)')     '  sigma (x) min / max             : ', dsqrt(minval(sig2(1,:))), dsqrt(maxval(sig2(1,:)))
      write(outfile,'(a,f18.10,a)')   '  dE_X  (a^2 term)                : ', dEx/dble(ncell), ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_X  (all orders)              : ', dExa/dble(ncell), ' Eh/cell'
      write(outfile,'(a)') ''
      write(outfile,'(a)') '  Summary (Eh per cell)                 a^2 term          all orders'
      write(outfile,'(a,2f18.10)')    '  dE_ES                           : ', dE_ES, dE_ES_all
      write(outfile,'(a,2f18.10)')    '  dE_X                            : ', dEx/dble(ncell), dExa/dble(ncell)
      write(outfile,'(a,2f18.10)')    '  dE_ES + dE_X                    : ', dE_ES + dEx/dble(ncell), dE_ES_all + dExa/dble(ncell)
      write(outfile,'(a,f18.10)')     '  E_HF per cell                   : ', EHF/dble(ncell)
      write(outfile,'(a,2f18.10)')    '  E_HF per cell - dE_ES - dE_X    : ', &
           (EHF - dEx)/dble(ncell) - dE_ES, (EHF - dExa)/dble(ncell) - dE_ES_all
      write(outfile,'(a)') ''
      write(outfile,'(a)') '  Exchange refinement A: pair terms i /= j (dipolar-Gaussian model)'
      write(outfile,'(a,i10)')        '  pairs with |d_ij|^2 > 1e-10     : ', npair
      write(outfile,'(a,f18.10,a)')   '  dE_X pairs (model A)            : ', dXp/dble(ncell), ' Eh/cell'
      write(outfile,'(a,f18.10)')     '  E_HF - dE_ES - dE_X - dE_Xpair  : ', (EHF - dExa - dXp)/dble(ncell) - dE_ES_all

      ! [Nov 2026] option B (exact first-order exchange), mixed representation, keyword fs_correction_exact
      if (c_fs_exact) then
      w0B = omp_get_wtime()
      call exchange_B_lowd()
      write(outfile,'(a)') ''
      write(outfile,'(a)') '  Exchange refinement B: exact pair terms + self-term shape (all orders)'
      write(outfile,'(a,2i6)')        '  periodic grid |m| <=            : ', XB_M
      write(outfile,'(a,2i6)')        '  k nodes (radial x angular)      : ', XB_Nk, XB_Nang
      write(outfile,'(a,3f12.6)')     '  lattice dipole factors kappa    : ', XB_kap
      write(outfile,'(a,2i8)')        '  reference-cell orbitals, pairs  : ', XB_nhome, XB_npair
      write(outfile,'(a,es12.3)')     '  check max | |rho_ii(0)|^2 - 1 |  : ', XB_chk
      write(outfile,'(a,es10.1,f12.0)') '  term threshold, mean terms/nu   : ', XB_tthr, XB_nterm
      write(outfile,'(a,f8.2)')       '  g cutoff (bohr^-1)              : ', XB_gcut
      write(outfile,'(a,3f14.10)')    '  self shape at 0.6/0.8/1.0 g_cut : ', XB_spart
      write(outfile,'(a,3f14.10)')    '  pair kernel at 0.6/0.8/1.0 g_cut: ', XB_part
      write(outfile,'(a,f18.10,a)')   '  dE_X self-term shape (extrap 1/g)   : ', XB_self, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_X pairs, kernel (extrap 1/g^2)   : ', XB_pair, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_X pairs, image (dipole) part : ', XB_lat, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_X pairs from outer 25% of R  : ', XB_out, ' Eh/cell'
      write(outfile,'(a,f18.10,a)')   '  dE_X (B) total                  : ', dExa/dble(ncell) + XB_self + XB_pair + XB_lat, ' Eh/cell'
      write(outfile,'(a,f18.10)')     '  E_HF - dE_ES - dE_X(B)          : ', &
           EHF/dble(ncell) - dE_ES_all - (dExa/dble(ncell) + XB_self + XB_pair + XB_lat)
      write(outfile,'(a,f10.2,a)')    '  wall time (B)                   : ', omp_get_wtime()-w0B, ' s'
      else
        write(outfile,'(a)') ''
        write(outfile,'(a)') '  (exact exchange, option B: keyword fs_correction_exact)'
      end if
      write(outfile,'(a,f10.2,a)')    '  wall time                       : ', omp_get_wtime()-w0, ' s'
      write(outfile,'(a)') ' ------------------------------------------------------------------'
      flush(outfile)
      deallocate(z, x2, zci, zcj, sig2, Rx, C, cci, ccj)

contains

      !=================================================================!
      ! [Nov 2026] option B for Tori1D / Tori2D (mixed representation)   !
      !   g = (g_p, k): g_p on the periodic supercell lattice (m_d a_d),  !
      !   k continuous in the open directions (log radial grid x angles   !
      !   in 1D, +-kappa in 2D, plus a k = 0 node for the small-q data).  !
      !   dE_X = -sum_i [d_Gauss(i) + S_i] - sum_{i/=j} [P_ij + L_ij]     !
      !   P_ij = (1/V_p) sum_gp int dk/(2pi)^no |rho_ij|^2 D(g_p,|k|)     !
      !   S_i  = same with |rho_ii|^2 - Gaussian (Resta along periodic,   !
      !          exact moments along open directions) - image term       !
      !   L_ij = -(4 pi/V_p) sum_c kappa_c d_c^2  (c over all 3 dims)     !
      !   kappa for the mixed lattice: sum_c kappa_c = 0;                 !
      !   1D: kappa ~ 1/L^2, 2D: kappa ~ 1/L (image line / sheet of dipoles)
      !   D = v_L - 4pi/g^2, v_L = (2/sqrt pi) int dt prod_per L I~_m(z)  !
      !       prod_open (sqrt pi/t) e^{-k^2/4t^2}                          !
      !=================================================================!
      subroutine exchange_B_lowd()
        double precision :: GCUT, CTHR, Vp, cellv(2), wfac, dd(3), dist, rcut
        double precision, allocatable :: kap(:), wk(:), kvec(:,:), Wkv(:), Dt(:,:,:), cen_o(:,:), sig_o(:)
        integer, allocatable :: home(:)
        integer :: Nk, Nang, Nkv, Mv(2), ii, kk, nh, q, ja, kv, ih, np_loc, ln, ist, Wc
        double precision :: wn(2)
        double precision :: kmin, kmax, dss, Rmax, Sacc(3), Pacc, Pp(3), Lacc, Oacc, chkv
        character(len=32) :: envs

        Vp = Lv(1) ; if (ndp == 2) Vp = Lv(1)*Lv(2)
        cellv = 1.d0
        do kk = 1, ndp
          cellv(kk) = Lv(kk)/dble(ncellv(kk))
        end do
        GCUT = 12.d0 ; XB_tthr = 1.d-9 ; CTHR = 1.d-4
        call get_environment_variable('FS_GCUT', envs, ln, ist) ; if (ist == 0 .and. ln > 0) read(envs,*) GCUT
        call get_environment_variable('FS_TTHR', envs, ln, ist) ; if (ist == 0 .and. ln > 0) read(envs,*) XB_tthr
        call get_environment_variable('FS_CTHR', envs, ln, ist) ; if (ist == 0 .and. ln > 0) read(envs,*) CTHR
        XB_gcut = GCUT ; XB_nterm = 0.d0
        Mv = 0
        do kk = 1, ndp
          Mv(kk) = max(2, ceiling(GCUT/av(kk)))
        end do
        XB_M = Mv

        ! centres (periodic: Berry phase; open: <x>), widths, reference cell (periodic dims)
        allocate(cen_o(3,nocc), sig_o(nocc), home(nocc))
        nh = 0 ; Rmax = 0.d0
        do ii = 1, nocc
          do kk = 1, 3
            if (kk <= ndp) then
              cen_o(kk,ii) = modulo(-atan2(aimag(z(ii,ii,kk)), dble(z(ii,ii,kk)))/av(kk), Lv(kk))
            else
              cen_o(kk,ii) = dble(z(ii,ii,kk))
            end if
          end do
          sig_o(ii) = dsqrt(sum(sig2(:,ii)))
        end do
        ! [Nov 2026] reference window: W cells along direction 1 (one cell along the others), the smallest W
        ! holding at least one orbital centre (e.g. 1 orbital per 2 cells in an equidistant H chain)
        do Wc = 1, ncellv(1)
          nh = 0
          do ii = 1, nocc
            wn = cellv ; wn(1) = dble(Wc)*cellv(1)
            if (all(modulo(cen_o(1:ndp,ii) - 0.1234d0*cellv(1:ndp), Lv(1:ndp)) < wn(1:ndp))) then
              nh = nh + 1 ; home(nh) = ii
            end if
          end do
          if (nh > 0) exit
        end do
        ! [Nov 2026] angular grid: the phase exp(-i k.c) needs ~2 k_max |c_perp| points; the orbital shapes
        ! themselves only need low harmonics (16 points). Before: + 3 sigma, i.e. 128 points for atoms on the axis.
        do ii = 1, nocc
          Rmax = max(Rmax, dsqrt(sum((cen_o(ndp+1:3,ii) - sum(cen_o(ndp+1:3,:),2)/dble(nocc))**2)))
        end do
        XB_nhome = nh
        wfac = dble(nocc)/(dble(ncell)*dble(max(nh,1)))

        ! k grid: node 1 = k = 0 (weight 0, small-q data), then log radial x angles
        kmin = 1.d-2*min(1.d0, minval(av(1:ndp))) ; kmax = GCUT ; dss = 0.05d0
        Nk = int(dlog(kmax/kmin)/dss) + 1
        allocate(kap(0:Nk), wk(0:Nk))
        kap(0) = 0.d0 ; wk(0) = 0.d0
        do q = 1, Nk
          kap(q) = kmin*dexp(dble(q-1)*dss)
          wk(q) = dss*kap(q) ; if (q == 1 .or. q == Nk) wk(q) = 0.5d0*wk(q)
        end do
        if (no == 2) then
          Nang = min(128, max(16, 2*ceiling(kmax*Rmax) + 16))
        else
          Nang = 2
        end if
        XB_Nk = Nk ; XB_Nang = Nang
        Nkv = 1 + Nk*Nang
        allocate(kvec(2,Nkv), Wkv(Nkv))
        kvec = 0.d0 ; Wkv = 0.d0
        do q = 1, Nk
          do ja = 1, Nang
            kv = 1 + (q-1)*Nang + ja
            if (no == 2) then
              kvec(1,kv) = kap(q)*cos(2.d0*pi*dble(ja-1)/dble(Nang))
              kvec(2,kv) = kap(q)*sin(2.d0*pi*dble(ja-1)/dble(Nang))
              Wkv(kv) = wk(q)*kap(q)*(2.d0*pi/dble(Nang))/(4.d0*pi*pi)
            else
              kvec(1,kv) = kap(q)*dble(3-2*ja)
              Wkv(kv) = wk(q)/(2.d0*pi)
            end if
          end do
        end do

        allocate(Dt(0:Mv(1),0:Mv(2),0:Nk))
        call kernel_D_mixed(Mv, Nk, kap, Dt)
        call kappa_mixed(XB_kap)

        ! AO lists of the localized orbitals
        allocate(nlist(nocc), alist(nBas,nocc))
        do ii = 1, nocc
          nlist(ii) = 0
          do kk = 1, nBas
            if (abs(C(kk,ii)) > CTHR*maxval(abs(C(:,ii)))) then
              nlist(ii) = nlist(ii) + 1 ; alist(nlist(ii),ii) = kk
            end if
          end do
        end do

        Sacc = 0.d0 ; Pacc = 0.d0 ; Pp = 0.d0 ; Lacc = 0.d0 ; Oacc = 0.d0 ; chkv = 0.d0 ; np_loc = 0
        do ih = 1, nh
          call home_block_mixed(home(ih), Mv, Nk, Nkv, kap, kvec, Wkv, Dt, cen_o, sig_o, Vp, &
                                Sacc, Pacc, Pp, Lacc, Oacc, chkv, np_loc)
        end do
        XB_spart = -wfac*Sacc
        XB_part  = -wfac*Pp
        XB_self  = (1.0d0*XB_spart(3) - 0.8d0*XB_spart(2))/(1.0d0 - 0.8d0)
        XB_pair  = (1.0d0*XB_part(3) - 0.64d0*XB_part(2))/(1.0d0 - 0.64d0)
        XB_lat   = -wfac*Lacc
        XB_out   = -wfac*Oacc
        XB_chk   = chkv
        XB_npair = np_loc
        XB_nterm = XB_nterm/dble(max(1,nh))
        deallocate(kap, wk, kvec, Wkv, Dt, cen_o, sig_o, home, nlist, alist)
      end subroutine exchange_B_lowd

      integer function ncellv(k)                          ! cells along periodic direction k
        use unitcell_module, only : nx, ny
        integer, intent(in) :: k
        if (k == 1) ncellv = max(nx,1)
        if (k == 2) ncellv = max(ny,1)
      end function ncellv

      ! D(m_p, kappa) = v_L - 4 pi/g^2 for |m_d| <= M_d, kappa nodes 0..Nk
      subroutine kernel_D_mixed(Mv, Nk, kap, Dt)
        integer, intent(in) :: Mv(2), Nk
        double precision, intent(in) :: kap(0:Nk)
        double precision, intent(out) :: Dt(0:Mv(1),0:Mv(2),0:Nk)
        double precision, allocatable :: Ivt(:), fac(:,:), Dl(:,:,:)   ! [Nov 2026] Dl: per-thread sum on the heap
        double precision :: smin, smax, dss, t, wt, zz, g2, P0, op, gp2
        integer :: ns, is, dd, nb, m1, m2, q
        smin = dlog(min(kap(1), minval(av(1:ndp)))/40.d0) ; smax = dlog(40.d0*XB_gcut)
        dss = 0.01d0 ; ns = int((smax-smin)/dss) + 1
        Dt = 0.d0
        !$omp parallel default(shared) private(is,t,wt,dd,zz,nb,Ivt,fac,m1,m2,q,g2,P0,op,gp2,Dl)
        allocate(fac(0:maxval(Mv),2), Dl(0:Mv(1),0:Mv(2),0:Nk)) ; Dl = 0.d0
        ! [Nov 2026] was reduction(+:Dt): the private copy lived on the thread stack -> SIGBUS for N = 100 on macOS
        !$omp do schedule(dynamic,4)
        do is = 0, ns
          t = dexp(smin + dble(is)*dss)
          wt = t*dss ; if (is == 0 .or. is == ns) wt = 0.5d0*wt
          fac = 1.d0
          do dd = 1, ndp
            zz = 2.d0*t*t/av(dd)**2
            nb = Mv(dd) + 40 + int(10.d0*dsqrt(zz))
            allocate(Ivt(0:nb))
            call bessel_scaled_vector_fast(nb, zz, Ivt)
            fac(0:Mv(dd),dd) = Lv(dd)*Ivt(0:Mv(dd))
            deallocate(Ivt)
          end do
          do q = 0, Nk
            op = (dsqrt(pi)/t)**no*dexp(-kap(q)**2/(4.d0*t*t))
            if (op < 1.d-300) cycle
            do m2 = 0, Mv(2)
              do m1 = 0, Mv(1)
                gp2 = (m1*av(1))**2 ; if (ndp == 2) gp2 = gp2 + (m2*av(2))**2
                g2 = gp2 + kap(q)**2
                if (g2 <= 0.d0 .or. g2 > XB_gcut**2*1.01d0) cycle
                P0 = pi**1.5d0/t**3*dexp(-g2/(4.d0*t*t))
                if (ndp == 2) then
                  Dl(m1,m2,q) = Dl(m1,m2,q) + wt*(fac(m1,1)*fac(m2,2)*op - P0)
                else
                  Dl(m1,m2,q) = Dl(m1,m2,q) + wt*(fac(m1,1)*op - P0)
                end if
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
      end subroutine kernel_D_mixed

      ! kappa_c for the mixed lattice (Gaussian regulator; open k integrals analytic)
      !   1D: periodic g^2 E1(eps g^2) ; open (1/2)[e^{-eps g^2}/eps - g^2 E1]
      !   2D: periodic 2 pi g_c^2/g erfc(sqrt(eps) g) ; open 2[sqrt(pi/eps) e^{-eps g^2} - pi g erfc]
      subroutine kappa_mixed(kp)
        double precision, intent(out) :: kp(3)
        double precision :: eps, cont, sper(2), sopn, g(2), gg, Vp
        integer :: Mk(2), m1, m2
        eps = 0.02d0
        Vp = Lv(1) ; if (ndp == 2) Vp = Lv(1)*Lv(2)
        Mk = 0
        do m1 = 1, ndp
          Mk(m1) = int(8.d0/(dsqrt(eps)*av(m1))) + 1
        end do
        cont = (4.d0*pi/3.d0)*4.d0*pi/(2.d0*pi)**3*dsqrt(pi)/(4.d0*eps**1.5d0)
        sper = 0.d0 ; sopn = 0.d0
        do m2 = -Mk(2), Mk(2)
          do m1 = -Mk(1), Mk(1)
            g = 0.d0 ; g(1) = m1*av(1) ; if (ndp == 2) g(2) = m2*av(2)
            gg = dsqrt(sum(g**2))
            if (ndp == 1) then
              if (gg > 0.d0) then
                sper(1) = sper(1) + gg*gg*e1x(eps*gg*gg)
                sopn = sopn + 0.5d0*(dexp(-eps*gg*gg)/eps - gg*gg*e1x(eps*gg*gg))
              else
                sopn = sopn + 0.5d0/eps
              end if
            else
              if (gg > 0.d0) then
                sper = sper + 2.d0*pi*g**2/gg*erfc(dsqrt(eps)*gg)
                sopn = sopn + 2.d0*(dsqrt(pi/eps)*dexp(-eps*gg*gg) - pi*gg*erfc(dsqrt(eps)*gg))
              else
                sopn = sopn + 2.d0*dsqrt(pi/eps)
              end if
            end if
          end do
        end do
        kp = 0.d0
        kp(1) = -(Vp/(4.d0*pi))*(sper(1)/Vp - cont)
        if (ndp == 2) kp(2) = -(Vp/(4.d0*pi))*(sper(2)/Vp - cont)
        do m1 = ndp+1, 3
          kp(m1) = -(Vp/(4.d0*pi))*(sopn/Vp - cont)
        end do
      end subroutine kappa_mixed

      double precision function e1x(x)                    ! exponential integral E1(x), x > 0
        double precision, intent(in) :: x
        double precision :: s, tm, b, c, d, h, del
        integer :: k
        if (x < 1.d0) then
          s = 0.d0 ; tm = 1.d0
          do k = 1, 60
            tm = -tm*x/dble(k)
            s = s + tm/dble(k)
            if (abs(tm/dble(k)) < 1.d-17) exit
          end do
          e1x = -0.5772156649015329d0 - dlog(x) - s
        else
          b = x + 1.d0 ; c = 1.d300 ; d = 1.d0/b ; h = d
          do k = 1, 300
            tm = -dble(k)*dble(k)
            b = b + 2.d0
            d = 1.d0/(tm*d + b) ; c = b + tm/c
            del = c*d ; h = h*del
            if (abs(del-1.d0) < 1.d-15) exit
          end do
          e1x = h*dexp(-x)
        end if
      end function e1x

      ! one reference orbital i: half transforms H_nu, all pairs j by zgemm, energies
      subroutine home_block_mixed(io, Mv, Nk, Nkv, kap, kvec, Wkv, Dt, cen_o, sig_o, Vp, &
                                  Sacc, Pacc, Pp, Lacc, Oacc, chkv, np_loc)
        integer, intent(in) :: io, Mv(2), Nk, Nkv
        double precision, intent(in) :: kap(0:Nk), kvec(2,Nkv), Wkv(Nkv), Dt(0:Mv(1),0:Mv(2),0:Nk), cen_o(3,nocc), sig_o(nocc), Vp
        double precision, intent(inout) :: Sacc(3), Pacc, Pp(3), Lacc, Oacc, chkv
        integer, intent(inout) :: np_loc
        integer, allocatable :: jl(:), nub(:)
        double precision, allocatable :: gj(:), rc(:), dj(:)
        complex(8), allocatable :: R(:,:), Hb(:,:), Cb(:,:)
        integer :: nj, jj, kk, Mi(2), NP1, NP2, Npt, Gt, nbk, b0, b1, nbl, ib, nu, jx, m1, m2, kv, q, idx, ip0
        double precision :: dd(3), dist, rcut, gmax, g2, gp2, r2, gau, e_s(3), e_p(3), e_l, y1, y2, c0, f1, f2, d2k, ws, N0, zr(2)
        complex(8) :: v1, v2

        allocate(jl(nocc), gj(nocc), rc(nocc), dj(nocc))
        nj = 0
        do jj = 1, nocc
          dd = cen_o(:,jj) - cen_o(:,io)
          dd(1:ndp) = dd(1:ndp) - Lv(1:ndp)*dnint(dd(1:ndp)/Lv(1:ndp))
          dist = dsqrt(sum(dd**2))
          rcut = 4.d0*(sig_o(io) + sig_o(jj)) + 2.d0
          if (dist > rcut .and. jj /= io) cycle
          nj = nj + 1 ; jl(nj) = jj ; rc(nj) = rcut ; dj(nj) = dist
          gj(nj) = min(XB_gcut, max(4.d0, 5.d0/dsqrt(min(minval(sig2(:,io)), minval(sig2(:,jj))))))
        end do
        gmax = maxval(gj(1:nj))
        Mi = 0
        do kk = 1, ndp
          Mi(kk) = min(Mv(kk), max(2, ceiling(gmax/av(kk))))
        end do
        NP1 = 2*Mi(1)+1 ; NP2 = 2*Mi(2)+1 ; Npt = NP1*NP2 ; Gt = Npt*Nkv
        ip0 = (Mi(1)+1) + NP1*Mi(2)                       ! flat periodic index of m = 0
        allocate(R(Gt,nj)) ; R = (0.d0,0.d0)

        nbk = max(2, min(64, int(4.d8/(16.d0*dble(Gt)))))
        allocate(nub(nBas))
        nbl = 0
        do nu = 1, nBas
          if (maxval(abs(C(nu,jl(1:nj)))) < 1.d-6) cycle
          nbl = nbl + 1 ; nub(nbl) = nu
        end do
        allocate(Hb(Gt,nbk), Cb(nbk,nj))
        do b0 = 1, nbl, nbk
          b1 = min(nbl, b0+nbk-1)
          !$omp parallel do default(shared) private(ib) schedule(dynamic,1)
          do ib = b0, b1
            call half_transform_mixed(io, nub(ib), Mi, Nkv, kvec, Hb(:,ib-b0+1))
          end do
          !$omp end parallel do
          do ib = b0, b1
            Cb(ib-b0+1,:) = dcmplx(C(nub(ib),jl(1:nj)), 0.d0)
          end do
          call zgemm('N', 'N', Gt, nj, b1-b0+1, (1.d0,0.d0), Hb, Gt, Cb, nbk, (1.d0,0.d0), R, Gt)
        end do

        !$omp parallel do default(shared) private(jx,jj,e_s,e_p,e_l,m1,m2,kv,q,idx,g2,gp2,r2,gau,kk,y1,y2,c0,f1,f2,d2k,v1,v2,ws,N0,zr) &
        !$omp reduction(+:Sacc,Pacc,Pp,Lacc,Oacc,np_loc) reduction(max:chkv) schedule(dynamic,1)
        do jx = 1, nj
          jj = jl(jx)
          e_s = 0.d0 ; e_p = 0.d0
          ! reference Gaussian matched to the computed rho_ii itself: norm N0 = |rho(0,0)|^2 and
          ! periodic decay |rho(e_d,0)|^2/N0, so that rho - Gaussian -> 0 as k -> 0 (D ~ -4pi/k^2 there)
          if (jj == io) then
            N0 = abs(R(pflat(Mi,1,0),jx))**2
            do kk = 1, ndp
              zr(kk) = abs(R(pflat(Mi,kk,1),jx))**2/N0
            end do
            chkv = max(chkv, abs(N0 - 1.d0))
          end if
          do kv = 2, Nkv
            q = (kv-2)/((Nkv-1)/Nk) + 1           ! radial node of kv (kv = 1 is k = 0)
            do m2 = -Mi(2), Mi(2)
              do m1 = -Mi(1), Mi(1)
                gp2 = (m1*av(1))**2 ; if (ndp == 2) gp2 = gp2 + (m2*av(2))**2
                g2 = gp2 + kap(q)**2
                if (g2 > gj(jx)**2) cycle
                idx = (m1+Mi(1)+1) + NP1*(m2+Mi(2)) + Npt*(kv-1)
                if (jj == io) then
                  gau = N0*zr(1)**(m1*m1)
                  if (ndp == 2) gau = gau*zr(2)**(m2*m2)
                  do kk = 1, no
                    gau = gau*dexp(-sig2(ndp+kk,io)*kvec(kk,kv)**2)
                  end do
                  r2 = Wkv(kv)*(abs(R(idx,jx))**2 - gau)*Dt(abs(m1),abs(m2),q)
                  if (g2 <= (0.6d0*gj(jx))**2) e_s(1) = e_s(1) + r2
                  if (g2 <= (0.8d0*gj(jx))**2) e_s(2) = e_s(2) + r2
                  e_s(3) = e_s(3) + r2
                else
                  r2 = Wkv(kv)*abs(R(idx,jx))**2*Dt(abs(m1),abs(m2),q)
                  if (g2 <= (0.6d0*gj(jx))**2) e_p(1) = e_p(1) + r2
                  if (g2 <= (0.8d0*gj(jx))**2) e_p(2) = e_p(2) + r2
                  e_p(3) = e_p(3) + r2
                end if
              end do
            end do
          end do
          e_s = e_s/Vp ; e_p = e_p/Vp
          if (jj == io) then
            do kk = 1, ndp                                  ! periodic directions: image of the curvature mismatch
              v1 = R(pflat(Mi,kk,1),jx) ; v2 = R(pflat(Mi,kk,2),jx)
              y1 = dlog(max(abs(v1)**2, 1.d-300)) ; y2 = dlog(max(abs(v2)**2, 1.d-300))
              c0 = -y1/av(kk)**2 + (16.d0*y1 - y2)/(12.d0*av(kk)**2)
              e_s = e_s - 4.d0*pi/Vp*XB_kap(kk)*c0
            end do
            Sacc = Sacc + e_s
          else
            np_loc = np_loc + 1
            Pacc = Pacc + e_p(3) ; Pp = Pp + e_p
            e_l = 0.d0
            do kk = 1, 3
              if (kk <= ndp) then
                f1 = abs(R(pflat(Mi,kk,1),jx))**2/av(kk)**2
                f2 = abs(R(pflat(Mi,kk,2),jx))**2/(4.d0*av(kk)**2)
                if (f1 > 1.d-300 .and. f2 > 1.d-300) then
                  d2k = dexp((4.d0*dlog(f1) - dlog(f2))/3.d0)
                else
                  d2k = f1
                end if
              else
                d2k = dble(z(io,jj,kk))**2                     ! open: exact <i|x|j>
              end if
              e_l = e_l - 4.d0*pi/Vp*XB_kap(kk)*d2k
            end do
            Lacc = Lacc + e_l
            if (dj(jx) > 0.75d0*rc(jx)) Oacc = Oacc + e_p(3) + e_l
          end if
        end do
        !$omp end parallel do
        deallocate(jl, gj, rc, dj, R, Hb, Cb, nub)

      end subroutine home_block_mixed

      integer function pflat(Mi, k, m)                    ! flat index of m e_k (periodic k) at the k = 0 node
        integer, intent(in) :: Mi(2), k, m
        integer :: q(2)
        q = 0 ; q(k) = m
        pflat = (q(1)+Mi(1)+1) + (2*Mi(1)+1)*(q(2)+Mi(2))
      end function pflat

      ! H_nu(m_p, k) = sum_{mu in list(io)} C_mu,io (mu nu)(g_p, k); flat (periodic fastest, then k node)
      subroutine half_transform_mixed(io, nu, Mi, Nkv, kvec, H)
        integer, intent(in) :: io, nu, Mi(2), Nkv
        double precision, intent(in) :: kvec(2,Nkv)
        complex(8), intent(out) :: H((2*Mi(1)+1)*(2*Mi(2)+1), Nkv)
        complex(8), allocatable :: Fp(:,:), Fo(:,:), f1(:), f2(:)
        double precision, allocatable :: Iv(:)
        complex(8) :: dco(-8:8), cv, oo
        double precision :: r1(3), r2(3), pe(2), cen(2), e1, e2, A, coef(2), pref, opf
        integer :: a1, mu, ip, jp, dd, Mm, k, m, nT, nTmax, nb, kv, m1, m2, Npt
        Npt = (2*Mi(1)+1)*(2*Mi(2)+1)
        nTmax = 0
        do a1 = 1, nlist(io)
          nTmax = nTmax + size(AOn(alist(a1,io))%exponent)*size(AOn(nu)%exponent)
        end do
        allocate(Fp(Npt,nTmax), Fo(Nkv,nTmax), f1(-Mi(1):Mi(1)), f2(-Mi(2):Mi(2)), Iv(0:40000))
        nT = 0
        r2 = (/ AOn(nu)%x, AOn(nu)%y, AOn(nu)%z /)
        do a1 = 1, nlist(io)
          mu = alist(a1,io)
          r1 = (/ AOn(mu)%x, AOn(mu)%y, AOn(mu)%z /)
          do ip = 1, size(AOn(mu)%exponent)
            e1 = AOn(mu)%exponent(ip)
            do jp = 1, size(AOn(nu)%exponent)
              e2 = AOn(nu)%exponent(jp)
              call clifford_pair(ndp, e1, e2, r1, r2, pe(1:ndp), cen(1:ndp))
              pref = C(mu,io)*AOn(mu)%coefficient(ip)*AOn(nu)%coefficient(jp)
              opf = dexp(-e1*e2/(e1+e2)*sum((r1(ndp+1:3)-r2(ndp+1:3))**2))
              coef = 1.d0
              do dd = 1, ndp
                coef(dd) = Lv(dd)*dexp(-2.d0*(e1+e2-pe(dd))/av(dd)**2)/av(dd)**(lq(dd,mu)+lq(dd,nu))
              end do
              if (abs(pref*coef(1)*coef(2)*opf) < XB_tthr) cycle
              nT = nT + 1
              f2 = (1.d0,0.d0)
              do dd = 1, ndp
                Mm = lq(dd,mu) + lq(dd,nu)
                A = 2.d0*pe(dd)/av(dd)**2
                call fourier_sin_product(lq(dd,mu), lq(dd,nu), av(dd)*(cen(dd)-r1(dd)), av(dd)*(cen(dd)-r2(dd)), dco(-Mm:Mm))
                nb = min(40000, max(Mi(dd)+Mm+2, Mm + 40 + int(10.d0*dsqrt(A))))
                call bessel_scaled_vector_fast(nb, A, Iv(0:nb))
                do m = -Mi(dd), Mi(dd)
                  cv = (0.d0,0.d0)
                  do k = -Mm, Mm
                    cv = cv + dco(k)*Iv(abs(m-k))
                  end do
                  cv = coef(dd)*cv*exp(dcmplx(0.d0, -dble(m)*av(dd)*cen(dd)))
                  if (dd == 1) f1(m) = cv
                  if (dd == 2) f2(m) = cv
                end do
              end do
              do m2 = -Mi(2), Mi(2)
                do m1 = -Mi(1), Mi(1)
                  Fp((m1+Mi(1)+1) + (2*Mi(1)+1)*(m2+Mi(2)), nT) = pref*f1(m1)*f2(m2)
                end do
              end do
              do kv = 1, Nkv
                oo = (1.d0,0.d0)
                do dd = 1, no
                  oo = oo*ftx(lq(ndp+dd,mu), lq(ndp+dd,nu), e1, e2, r1(ndp+dd), r2(ndp+dd), kvec(dd,kv))
                end do
                Fo(kv,nT) = oo
              end do
            end do
          end do
        end do
        if (nT == 0) then
          H = (0.d0,0.d0)
        else
          call zgemm('N', 'T', Npt, Nkv, nT, (1.d0,0.d0), Fp, Npt, Fo, Nkv, (0.d0,0.d0), H, Npt)
        end if
        !$omp atomic
        XB_nterm = XB_nterm + dble(nT)
        deallocate(Fp, Fo, f1, f2, Iv)
      end subroutine half_transform_mixed

      ! int (x-A)^l1 (x-B)^l2 exp(-e1(x-A)^2 - e2(x-B)^2) exp(-ikx) dx  (open direction)
      complex(8) function ftx(l1, l2, e1, e2, xa, xb, kk)
        integer, intent(in) :: l1, l2
        double precision, intent(in) :: e1, e2, xa, xb, kk
        double precision :: pp, xp, c(0:8)
        complex(8) :: Mj(0:9)
        integer :: jj
        pp = e1 + e2 ; xp = (e1*xa + e2*xb)/pp
        call poly_shift(l1, l2, xp-xa, xp-xb, 0, 0.d0, c)
        Mj(0) = dsqrt(pi/pp)*dexp(-kk*kk/(4.d0*pp))
        Mj(1) = dcmplx(0.d0, -kk/(2.d0*pp))*Mj(0)
        do jj = 1, l1+l2-1
          Mj(jj+1) = (dble(jj)*Mj(jj-1) + dcmplx(0.d0,-kk)*Mj(jj))/(2.d0*pp)
        end do
        ftx = (0.d0,0.d0)
        do jj = 0, l1+l2
          ftx = ftx + c(jj)*Mj(jj)
        end do
        ftx = ftx*dexp(-e1*e2/pp*(xa-xb)**2)*exp(dcmplx(0.d0, -kk*xp))
      end function ftx

      double precision function omega_val()
        integer :: ii, dd
        omega_val = 0.d0
        do dd = 1, 3
          do ii = 1, nocc
            omega_val = omega_val + wgt(dd)*abs(z(ii,ii,dd))**2
          end do
        end do
      end function omega_val

      ! supercell AO matrices:
      !   dir = 0           : S
      !   dir = d periodic  : (Ar, Ai) = <mu| e^{-i a_d x_d} |nu>       (pw = 1)
      !   dir = d open      : Ar = <mu| x_d^pw |nu>                      (pw = 1, 2)
      subroutine pair_mat(dir, pw, Ar, Ai)
        integer, intent(in) :: dir, pw
        double precision, intent(out) :: Ar(nBas,nBas)
        double precision, intent(out), optional :: Ai(nBas,nBas)
        double precision, allocatable :: Iv(:)
        complex(8) :: dco(-8:8), cc, val
        double precision :: r1(3), r2(3), pe(2), cen(2), e1, e2, A, coef, pref, oo, pp, xp, c(0:8), Mj
        integer :: mu, nu, ip, jp, dd, kk, Mm, nbes, mdir, jj, l3
        Ar = 0.d0
        if (present(Ai)) Ai = 0.d0
        !$omp parallel default(shared) private(mu,nu,ip,jp,dd,kk,Mm,nbes,mdir,r1,r2,pe,cen,e1,e2,A,coef,pref,dco,cc,val,Iv,oo,pp,xp,c,Mj,jj,l3)
        allocate(Iv(0:40000))
        !$omp do schedule(dynamic,4)
        do mu = 1, nBas
          r1 = (/ AOn(mu)%x, AOn(mu)%y, AOn(mu)%z /)
          do nu = mu, nBas
            r2 = (/ AOn(nu)%x, AOn(nu)%y, AOn(nu)%z /)
            val = (0.d0,0.d0)
            do ip = 1, size(AOn(mu)%exponent)
              e1 = AOn(mu)%exponent(ip)
              do jp = 1, size(AOn(nu)%exponent)
                e2 = AOn(nu)%exponent(jp)
                pref = AOn(mu)%coefficient(ip)*AOn(nu)%coefficient(jp)
                call clifford_pair(ndp, e1, e2, r1, r2, pe(1:ndp), cen(1:ndp))
                pref = pref*dexp(-e1*e2/(e1+e2)*sum((r1(ndp+1:3)-r2(ndp+1:3))**2))
                do dd = 1, ndp
                  pref = pref*Lv(dd)*dexp(-2.d0*(e1+e2-pe(dd))/av(dd)**2)/av(dd)**(lq(dd,mu)+lq(dd,nu))
                end do
                if (abs(pref) < 1.d-14) cycle
                cc = dcmplx(pref, 0.d0)
                do dd = 1, ndp
                  Mm = lq(dd,mu) + lq(dd,nu)
                  mdir = 0 ; if (dd == dir) mdir = 1
                  A = 2.d0*pe(dd)/av(dd)**2
                  call fourier_sin_product(lq(dd,mu), lq(dd,nu), av(dd)*(cen(dd)-r1(dd)), av(dd)*(cen(dd)-r2(dd)), dco(-Mm:Mm))
                  nbes = min(40000, Mm + 40 + int(10.d0*dsqrt(A)))
                  call bessel_scaled_vector_fast(nbes, A, Iv(0:nbes))
                  val = (0.d0,0.d0)
                  do kk = -Mm, Mm
                    val = val + dco(kk)*Iv(abs(mdir-kk))
                  end do
                  if (mdir == 1) val = val*exp(dcmplx(0.d0, -av(dd)*cen(dd)))
                  cc = cc*val
                end do
                ! open directions: moments int (x-A)^l1 (x-B)^l2 x^l3 e^{-p(x-P)^2}
                do dd = ndp+1, 3
                  pp = e1 + e2 ; xp = (e1*r1(dd) + e2*r2(dd))/pp
                  l3 = 0 ; if (dd == dir) l3 = pw
                  call poly_shift(lq(dd,mu), lq(dd,nu), xp-r1(dd), xp-r2(dd), l3, xp, c)
                  oo = 0.d0
                  do jj = 0, lq(dd,mu)+lq(dd,nu)+l3, 2
                    Mj = gamma(0.5d0*dble(jj+1))/pp**(0.5d0*dble(jj+1))
                    oo = oo + c(jj)*Mj
                  end do
                  cc = cc*oo
                end do
                Ar(mu,nu) = Ar(mu,nu) + dble(cc)
                if (present(Ai)) Ai(mu,nu) = Ai(mu,nu) + aimag(cc)
              end do
            end do
            Ar(nu,mu) = Ar(mu,nu)
            if (present(Ai)) Ai(nu,mu) = Ai(mu,nu)
          end do
        end do
        !$omp end do
        deallocate(Iv)
        !$omp end parallel
      end subroutine pair_mat

      ! <x_d^4/r^3> for r Gaussian with variances s2v (d periodic)
      double precision function x4r3(s2v, dd)
        double precision, intent(in) :: s2v(3)
        integer, intent(in) :: dd
        integer, parameter :: NQ = 4000
        integer :: q, kk
        double precision :: tau, uq, t, f, wq, ssum
        tau = 1.d0/dsqrt(sum(s2v)/3.d0)
        ssum = 0.d0
        do q = 1, NQ-1
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
        x4r3 = 4.d0/dsqrt(pi)*ssum/(3.d0*dble(NQ))
      end function x4r3

      ! all orders minus a^2 (Gaussian holes); open directions identical in torus and bulk
      subroutine x_rem(R)
        double precision, intent(out) :: R(nocc)
        double precision, allocatable :: Ivt(:,:), s2a(:,:)
        double precision :: smin, smax, dss, t, wt, F(3), Bq(3), a2t, smn
        integer :: ns, is, dd, nb(2), io, kk, kmx, nbm
        allocate(s2a(3,nocc))
        s2a = 2.d0*max(sig2, 1.d-12)
        smn = dsqrt(minval(s2a))
        smin = dlog(1.d-3*minval(av(1:ndp))) ; smax = dlog(200.d0/smn)
        dss = 0.01d0 ; ns = int((smax-smin)/dss) + 1
        R = 0.d0
        !$omp parallel default(shared) private(is,t,wt,dd,nb,nbm,Ivt,io,F,Bq,a2t,kk,kmx)
        !$omp do schedule(dynamic,4) reduction(+:R)
        do is = 0, ns
          t = dexp(smin + dble(is)*dss)
          wt = t*dss ; if (is == 0 .or. is == ns) wt = 0.5d0*wt
          nb = 0
          do dd = 1, ndp
            nb(dd) = int(12.d0/(av(dd)*smn)) + 40 + int(10.d0*dsqrt(2.d0*t*t/av(dd)**2))
          end do
          nbm = maxval(nb)
          allocate(Ivt(0:nbm,2))
          do dd = 1, ndp
            call bessel_scaled_vector_fast(nb(dd), 2.d0*t*t/av(dd)**2, Ivt(0:nb(dd),dd))
          end do
          do io = 1, nocc
            do dd = 1, 3
              Bq(dd) = 1.d0/dsqrt(1.d0 + 2.d0*t*t*s2a(dd,io))
              F(dd) = Bq(dd)
              if (dd <= ndp) then
                kmx = min(nb(dd), int(12.d0/(av(dd)*dsqrt(s2a(dd,io)))) + 30)
                F(dd) = Ivt(0,dd)
                do kk = 1, kmx
                  F(dd) = F(dd) + 2.d0*Ivt(kk,dd)*dexp(-0.5d0*dble(kk*kk)*av(dd)**2*s2a(dd,io))
                end do
              end if
            end do
            a2t = 0.d0
            do dd = 1, ndp
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
      end subroutine x_rem

end subroutine fs_exchange_lowd
