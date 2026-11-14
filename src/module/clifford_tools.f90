!=====================================================================!
!                                                                     !
!   module clifford_tools : shared tools for the Clifford-Gaussian    !
!   integrals in 1D, 2D and 3D                                        !
!                                                                     !
!   ndp = number of periodic directions:                              !
!     Tori1D : ndp = 1  (x periodic,   y,z real)                      !
!     Tori2D : ndp = 2  (x,y periodic, z real)                        !
!     Tori3D : ndp = 3  (x,y,z periodic)                              !
!                                                                     !
!   Contents                                                          !
!   --------                                                          !
!   init_clifford             binomials, double factorials and the    !
!                             derivative table T (call once)          !
!   clifford_bessel_ratios    R(e,nu) = exp(-A) I_nu(A) / A^e         !
!   Icliff_from_R             exp(-A) Int sin^n cos^m exp(A cos) dx   !
!   overlap_1D_clifford(_ext) 1D overlap factor, periodic direction   !
!   kinetic_1D_clifford       1D Laplacian factor, periodic direction !
!   overlap_1D_real           1D overlap factor, real direction       !
!   kinetic_1D_real           1D Laplacian factor, real direction     !
!   bessel_scaled_vector      exp(-A) I_n(A), n = 0..N (Miller)       !
!   bessel_scaled_vector_fast same, start at N, no allocation (ERI)   !
!   bessel_nmax_estimate      truncation order of exp(-A) I_n(A)      !
!   fourier_sin_product       Fourier coefficients of sin^i sin^l     !
!   setup_t_tables(ncap,ndp)  t nodes/weights + Bessel tables (NA,ERI)!
!   bary_dir                  barycentre exponent/centre along d      !
!   clifford_pair             bary_dir for d = 1..ndp                 !
!   hermite_E_1D              McMurchie-Davidson E coefficients       !
!   orbital_to_lxyz           'dxy' -> (1,1,0)                        !
!                                                                     !
!   t-tables (filled by setup_t_tables):                              !
!     nq_t            number of t nodes (64 or 128)                   !
!     t_node(k)       nodes   of  Int_0^inf f(t) dt                   !
!     t_weight(k)     weights                                         !
!     t_bessel(n,k,d) exp(-C) I_n(C),  C = 2 t_k^2 / a_d^2             !
!     t_nmax(k,d)     last order kept (1e-16 relative)                !
!     d = 1..ndp  (a_1 = ax, a_2 = ay, a_3 = az)                      !
!                                                                     !
!=====================================================================!

module clifford_tools

      implicit none

      private
      public :: NMAX_CLIFF, init_clifford
      public :: clifford_bessel_ratios, Icliff_from_R
      public :: overlap_1D_clifford, overlap_1D_clifford_ext, overlap_1D_real
      public :: kinetic_1D_clifford, kinetic_1D_real
      public :: bessel_scaled_vector, bessel_scaled_vector_fast, bessel_nmax_estimate
      public :: fourier_sin_product, hermite_E_1D, orbital_to_lxyz
      public :: setup_t_tables, nq_t, t_node, t_weight, t_bessel, t_nmax
      public :: bary_dir, clifford_pair

      ! ---- t-quadrature tables (nuclear attraction and ERIs) ---- !
      integer,          save              :: nq_t = 0
      double precision, save, allocatable :: t_node(:), t_weight(:), t_bessel(:,:,:)
      integer,          save, allocatable :: t_nmax(:,:)

      ! highest total power n+m of sin/cos handled (8 = g-g overlap,
      ! 10 = g-g kinetic). Raise if you need h, i functions.
      integer, parameter :: NMAX_CLIFF = 16
      integer, parameter :: NM = NMAX_CLIFF

      double precision, parameter :: pi = 3.14159265358979323846d0

      ! T(e,nu,s,m): coefficients of d^m/dA^m [A^-s I_s] = sum T A^-e I_nu
      double precision, save :: T(0:NM, 0:NM, 0:NM/2, 0:NM)
      double precision, save :: binom(0:NM, 0:NM)
      double precision, save :: dfact(-1:2*NM)
      logical,          save :: initialized = .false.

contains

!---------------------------------------------------------------------!
subroutine init_clifford()

      integer :: s, m, e, nu, n, k
      double precision :: Tcur(0:NM,0:NM), Tnew(0:NM,0:NM)

      if (initialized) return

      ! binomial coefficients
      binom = 0.d0
      binom(0,0) = 1.d0
      do n = 1, NM
        binom(n,0) = 1.d0
        do k = 1, n
          binom(n,k) = binom(n-1,k-1) + binom(n-1,k)   ! binom(n-1,n) = 0
        end do
      end do

      ! double factorials, (-1)!! = 0!! = 1
      dfact(-1) = 1.d0 ; dfact(0) = 1.d0
      do n = 1, 2*NM
        dfact(n) = dble(n) * dfact(max(n-2,-1))
      end do

      ! derivative table
      T = 0.d0
      do s = 0, NM/2
        Tcur = 0.d0
        Tcur(s,s) = 1.d0
        T(:,:,s,0) = Tcur
        do m = 1, NM - 2*s
          Tnew = 0.d0
          do e = 0, NM
            do nu = e, NM - 1
              if (Tcur(e,nu) == 0.d0) cycle
              Tnew(e,nu+1) = Tnew(e,nu+1) + Tcur(e,nu)
              if (nu > e) Tnew(e+1,nu) = Tnew(e+1,nu) + dble(nu-e)*Tcur(e,nu)
            end do
          end do
          Tcur = Tnew
          T(:,:,s,m) = Tcur
        end do
      end do

      initialized = .true.

end subroutine init_clifford

!---------------------------------------------------------------------!
! R(e,nu) = exp(-A) I_nu(A) / A^e   for 0 <= e <= nu <= nmax           !
! Call ONCE per primitive pair, then reuse R for every (n,m).          !
!---------------------------------------------------------------------!
subroutine clifford_bessel_ratios(nmax, A, R)

      use bessel_functions, only : iv_scaled

      integer,          intent(in)  :: nmax
      double precision, intent(in)  :: A
      double precision, intent(out) :: R(0:NM,0:NM)

      double precision, parameter   :: A_switch = 1.d0
      integer,          parameter   :: kmax = 14
      integer          :: e, nu, k
      double precision :: ive(0:NM), ser, term, q, eA, Ae

      R = 0.d0

      if (A < A_switch) then
        ! power series: e^-A A^(nu-e) sum_k (A^2/4)^k / (2^nu k! (k+nu)!)
        eA = dexp(-A)
        q  = 0.25d0 * A * A
        do nu = 0, nmax
          term = 1.d0
          do k = 1, nu
            term = term / (2.d0 * dble(k))          ! 1/(2^nu nu!)
          end do
          ser = term
          do k = 1, kmax
            term = term * q / (dble(k) * dble(k+nu))
            ser  = ser + term
          end do
          do e = 0, nu
            if (nu == e) then
              R(e,nu) = eA * ser
            else
              R(e,nu) = eA * ser * A**(nu-e)
            end if
          end do
        end do
      else
        do nu = 0, nmax
          ive(nu) = iv_scaled(nu, A)
        end do
        do e = 0, nmax
          Ae = A**e
          do nu = e, nmax
            R(e,nu) = ive(nu) / Ae
          end do
        end do
      end if

end subroutine clifford_bessel_ratios

!---------------------------------------------------------------------!
! Icliff(n,m) = exp(-A) Int_0^L sin^n cos^m exp(A cos(a x)) dx         !
! (n = power of sin, m = power of cos: same convention as Icliff in    !
!  your module functions). R from clifford_bessel_ratios.             !
!---------------------------------------------------------------------!
double precision function Icliff_from_R(n, m, R, L)

      integer,          intent(in) :: n, m
      double precision, intent(in) :: R(0:NM,0:NM), L
      integer :: s, e, nu
      double precision :: acc

      Icliff_from_R = 0.d0
      if (mod(n,2) /= 0) return

      s   = n / 2
      acc = 0.d0
      do e = s, s + m
        do nu = e, s + m
          if (T(e,nu,s,m) /= 0.d0) acc = acc + T(e,nu,s,m) * R(e,nu)
        end do
      end do
      Icliff_from_R = L * dfact(n-1) * acc

end function Icliff_from_R


!---------------------------------------------------------------------!
! 1D Clifford overlap factor for powers i (centre A) and l (centre B) !
! WITHOUT the prefactor kx = exp(-2(alpha+beta-px)/ax^2).             !
!---------------------------------------------------------------------!
double precision function overlap_1D_clifford(i, l, spa, cpa, spb, cpb, inv_ax, R, Lcell)

      integer,          intent(in) :: i, l
      double precision, intent(in) :: spa, cpa, spb, cpb, inv_ax, Lcell
      double precision, intent(in) :: R(0:NM,0:NM)

      overlap_1D_clifford = overlap_1D_clifford_ext(i, l, 0, 0, spa, cpa, spb, cpb, inv_ax, R, Lcell)

end function overlap_1D_clifford

!---------------------------------------------------------------------!
! Same, with an extra factor S^ds C^dc in the integrand,              !
!   S = sin(a(x-xP)), C = cos(a(x-xP))                                !
!   a^-(i+l) Int [sin a(x-xA)]^i [sin a(x-xB)]^l S^ds C^dc e^(A C)    !
!---------------------------------------------------------------------!
double precision function overlap_1D_clifford_ext(i, l, ds, dc, spa, cpa, spb, cpb, inv_ax, R, Lcell)

      integer,          intent(in) :: i, l, ds, dc
      double precision, intent(in) :: spa, cpa, spb, cpb, inv_ax, Lcell
      double precision, intent(in) :: R(0:NM,0:NM)
      integer :: k, q, ns, nc
      double precision :: acc, fa, fb

      overlap_1D_clifford_ext = 0.d0
      if (i < 0 .or. l < 0) return

      acc = 0.d0
      do k = 0, i
        fa = binom(i,k) * cpa**k * spa**(i-k)
        if (fa == 0.d0) cycle
        do q = 0, l
          ns = k + q + ds                 ! power of sin(a(x-xP))
          if (mod(ns,2) /= 0) cycle       ! odd sin power integrates to 0
          nc = (i - k) + (l - q) + dc     ! power of cos(a(x-xP))
          fb = binom(l,q) * cpb**q * spb**(l-q)
          acc = acc + fa * fb * Icliff_from_R(ns, nc, R, Lcell)
        end do
      end do
      overlap_1D_clifford_ext = acc * inv_ax**(i+l)

end function overlap_1D_clifford_ext

!---------------------------------------------------------------------!
! 1D Clifford Laplacian factor  < g_A^i | d2/dx2 | g_B^l >            !
! WITHOUT kx. With u = a(x-xB), g^l = (sin u / a)^l exp(B(cos u - 1)),!
! B = 2 beta/a^2:                                                     !
!   d2 g^l/dx2 = l(l-1) g^(l-2) - l^2 a^2 g^l                         !
!                - 2 beta (2l+1) cos(u) g^l + 4 beta^2 g^(l+2)        !
! and cos(u) = C cos(a(xP-xB)) - S sin(a(xP-xB)).                     !
! (a -> 0 gives the usual GTO formula.)                               !
! R must hold Bessel ratios up to order i+l+2.                        !
!---------------------------------------------------------------------!
double precision function kinetic_1D_clifford(i, l, beta, ax, spa, cpa, spb, cpb, R, Lcell)

      integer,          intent(in) :: i, l
      double precision, intent(in) :: beta, ax, spa, cpa, spb, cpb, Lcell
      double precision, intent(in) :: R(0:NM,0:NM)
      double precision :: inv_ax, acc, cosu

      inv_ax = 1.d0 / ax
      acc = 0.d0

      if (l >= 2) acc = acc + dble(l*(l-1)) * &
                    overlap_1D_clifford_ext(i, l-2, 0, 0, spa, cpa, spb, cpb, inv_ax, R, Lcell)

      if (l >= 1) acc = acc - dble(l*l) * ax * ax * &
                    overlap_1D_clifford_ext(i, l,   0, 0, spa, cpa, spb, cpb, inv_ax, R, Lcell)

      cosu = cpb * overlap_1D_clifford_ext(i, l, 0, 1, spa, cpa, spb, cpb, inv_ax, R, Lcell) &
           - spb * overlap_1D_clifford_ext(i, l, 1, 0, spa, cpa, spb, cpb, inv_ax, R, Lcell)
      acc = acc - 2.d0 * beta * dble(2*l+1) * cosu

      acc = acc + 4.d0 * beta * beta * &
                    overlap_1D_clifford_ext(i, l+2, 0, 0, spa, cpa, spb, cpb, inv_ax, R, Lcell)

      kinetic_1D_clifford = acc

end function kinetic_1D_clifford

!---------------------------------------------------------------------!
! 1D real-Gaussian overlap factor WITHOUT the prefactor exp(-mu R^2). !
! p = alpha+beta, pa = P-A, pb = P-B                                  !
!---------------------------------------------------------------------!
double precision function overlap_1D_real(j, m, p, pa, pb)

      integer,          intent(in) :: j, m
      double precision, intent(in) :: p, pa, pb
      integer :: k, q, n
      double precision :: acc, base, mom(0:NM)

      base = dsqrt(pi / p)
      do n = 0, j + m
        if (mod(n,2) == 0) then
          mom(n) = dfact(n-1) / (2.d0*p)**(n/2) * base
        else
          mom(n) = 0.d0
        end if
      end do

      acc = 0.d0
      do k = 0, j
        do q = 0, m
          if (mod(k+q,2) /= 0) cycle
          acc = acc + binom(j,k) * binom(m,q) * pa**(j-k) * pb**(m-q) * mom(k+q)
        end do
      end do
      overlap_1D_real = acc

end function overlap_1D_real

!---------------------------------------------------------------------!
! 1D real-Gaussian Laplacian factor < g_A^j | d2/dy2 | g_B^m >        !
! WITHOUT exp(-mu R^2):                                               !
!   m(m-1) S(j,m-2) - 2 beta (2m+1) S(j,m) + 4 beta^2 S(j,m+2)        !
!---------------------------------------------------------------------!
double precision function kinetic_1D_real(j, m, beta, p, pa, pb)

      integer,          intent(in) :: j, m
      double precision, intent(in) :: beta, p, pa, pb
      double precision :: acc

      acc = - 2.d0 * beta * dble(2*m+1) * overlap_1D_real(j, m, p, pa, pb) &
            + 4.d0 * beta * beta * overlap_1D_real(j, m+2, p, pa, pb)
      if (m >= 2) acc = acc + dble(m*(m-1)) * overlap_1D_real(j, m-2, p, pa, pb)
      kinetic_1D_real = acc

end function kinetic_1D_real

!---------------------------------------------------------------------!
! Scaled modified Bessel functions  Iv(n) = exp(-A) I_n(A), n=0..N     !
! Miller backward recurrence with rescaling (no overflow for small A), !
! normalised with iv_scaled(0,A).                                      !
!---------------------------------------------------------------------!
subroutine bessel_scaled_vector(N, A, Iv)

      use bessel_functions, only : iv_scaled

      integer,          intent(in)  :: N
      double precision, intent(in)  :: A
      double precision, intent(out) :: Iv(0:N)

      integer          :: n0, k
      double precision :: ip1, ik, im1, norm
      double precision, allocatable :: w(:)

      if (A < 1.d-12) then
        Iv    = 0.d0
        Iv(0) = 1.d0
        return
      end if

      n0 = N + 30 + int(dsqrt(80.d0*A))
      allocate(w(0:n0))
      ip1 = 0.d0 ; ik = 1.d-300
      w(n0) = ik
      do k = n0, 1, -1
        im1 = ip1 + (2.d0*dble(k)/A) * ik
        ip1 = ik ; ik = im1
        w(k-1) = ik
        if (abs(ik) > 1.d250) then          ! rescale everything computed so far
          w(k-1:n0) = w(k-1:n0) * 1.d-250
          ip1 = ip1 * 1.d-250 ; ik = ik * 1.d-250
        end if
      end do
      norm = iv_scaled(0, A) / w(0)
      Iv(0:N) = w(0:N) * norm
      deallocate(w)

end subroutine bessel_scaled_vector

!---------------------------------------------------------------------!
! Fast variant for the ERIs: N is already the truncation order (the    !
! values at N are negligible), so the backward recurrence starts at N  !
! exactly as bessel_I_scaled_backward of the original code; writes     !
! into the caller's array (no allocation), rescales against overflow.  !
!---------------------------------------------------------------------!
subroutine bessel_scaled_vector_fast(N, A, Iv)

      use bessel_functions, only : iv_scaled

      integer,          intent(in)  :: N
      double precision, intent(in)  :: A
      double precision, intent(out) :: Iv(0:N)

      integer          :: k
      double precision :: inv_A2

      if (A < 1.d-10) then
        Iv    = 0.d0
        Iv(0) = 1.d0
        return
      end if
      if (N < 2) then
        do k = 0, N
          Iv(k) = iv_scaled(k, A)
        end do
        return
      end if

      inv_A2 = 2.d0 / A
      Iv(N)   = 1.d-280
      Iv(N-1) = Iv(N) * dble(N) * inv_A2          ! I_{N+1} = 0 seed
      do k = N-1, 1, -1
        Iv(k-1) = Iv(k+1) + dble(k) * inv_A2 * Iv(k)
        if (Iv(k-1) > 1.d250) Iv(k-1:N) = Iv(k-1:N) * 1.d-250
      end do
      Iv(0:N) = Iv(0:N) * (iv_scaled(0, A) / Iv(0))

end subroutine bessel_scaled_vector_fast

!---------------------------------------------------------------------!
! Order beyond which exp(-A) I_n(A) / exp(-A) I_0(A) < ~1e-17          !
! ( ~ exp(-n^2/2A) for large A, (A/2)^n/n! for small A )               !
!---------------------------------------------------------------------!
integer function bessel_nmax_estimate(A)

      double precision, intent(in) :: A

      bessel_nmax_estimate = 20 + ceiling(9.d0 * dsqrt(max(A,0.d0)))

end function bessel_nmax_estimate

!---------------------------------------------------------------------!
! Fourier coefficients of  sin^i(u+phA) sin^l(u+phB) = sum_k d(k) e^{iku}
! k = -(i+l) .. (i+l),  d(-k) = conj(d(k))                             !
!---------------------------------------------------------------------!
subroutine fourier_sin_product(i, l, phA, phB, d)

      integer,          intent(in)  :: i, l
      double precision, intent(in)  :: phA, phB
      complex(8),       intent(out) :: d(-(i+l):(i+l))

      integer    :: m, k, deg
      complex(8) :: cp, cm, tmp(-(i+l)-1:(i+l)+1)

      d   = (0.d0,0.d0)
      d(0) = (1.d0,0.d0)
      deg = 0
      do m = 1, i + l
        if (m <= i) then
          cp =  exp(dcmplx(0.d0, phA)) / dcmplx(0.d0, 2.d0)   !  e^{i ph}/(2i)
          cm = -exp(dcmplx(0.d0,-phA)) / dcmplx(0.d0, 2.d0)   ! -e^{-i ph}/(2i)
        else
          cp =  exp(dcmplx(0.d0, phB)) / dcmplx(0.d0, 2.d0)
          cm = -exp(dcmplx(0.d0,-phB)) / dcmplx(0.d0, 2.d0)
        end if
        tmp = (0.d0,0.d0)
        do k = -deg, deg
          tmp(k+1) = tmp(k+1) + cp * d(k)
          tmp(k-1) = tmp(k-1) + cm * d(k)
        end do
        deg = deg + 1
        d(-deg:deg) = tmp(-deg:deg)
      end do

end subroutine fourier_sin_product

!---------------------------------------------------------------------!
! t-grid + Bessel tables for the ndp periodic directions             !
!   exp(-t^2 |x_d-C_d|_E^2) = sum_n Ĩ_n(C) e^{i n a_d (x_d-C_d)},      !
!   C = 2 t^2 / a_d^2                                                  !
! Same nodes/weights as the old ERI code (gauss_legendre_points,       !
! nquad_point = 64 or 128); t_val/w_val of precomputed_bessel are left !
! as they were found. Orders 0..ncap only (the pair vectors vanish     !
! beyond); per-t truncation t_nmax(k,d) at 1e-16 relative (+5).        !
! Called by the nuclear attraction and by the ERIs; rebuilt only if   !
! a larger ncap or more directions are requested.                     !
!---------------------------------------------------------------------!
subroutine setup_t_tables(ncap, ndp)

      use torus_init,         only : ax, ay, az
      use precomputed_bessel, only : nquad_point, t_val, w_val, gauss_legendre_points

      integer, intent(in) :: ncap , ndp

      logical          :: had_tw
      integer          :: k, nq, n, d
      double precision :: a(3)

      nq = nquad_point
      if (nq /= 64 .and. nq /= 128) nq = 64

      if (nq_t == nq .and. allocated(t_bessel)) then
        if (ubound(t_bessel,1) >= ncap .and. size(t_bessel,3) >= ndp) return
      end if
      if (allocated(t_node))   deallocate(t_node)
      if (allocated(t_weight)) deallocate(t_weight)
      if (allocated(t_bessel)) deallocate(t_bessel)
      if (allocated(t_nmax))   deallocate(t_nmax)

      had_tw = allocated(t_val)
      if (.not. had_tw) then
        nquad_point = nq
        call gauss_legendre_points()
      end if
      nq = size(t_val)
      allocate(t_node(nq), t_weight(nq), t_nmax(nq,ndp), t_bessel(0:ncap, nq, ndp))
      t_node   = t_val
      t_weight = w_val
      if (.not. had_tw) deallocate(t_val, w_val)    ! restore state

      a = (/ ax, ay, az /)
      !$omp parallel do default(shared) private(k,d) schedule(dynamic,1) collapse(2)
      do d = 1, ndp
        do k = 1, nq
          call bessel_scaled_vector(ncap, 2.d0 * t_node(k)**2 / (a(d)*a(d)), t_bessel(:,k,d))
        end do
      end do
      !$omp end parallel do

      do d = 1, ndp
        do k = 1, nq
          t_nmax(k,d) = ncap
          do n = 0, ncap
            if (t_bessel(n,k,d) < 1.d-16 * t_bessel(0,k,d)) then
              t_nmax(k,d) = min(ncap, n + 5)
              exit
            end if
          end do
        end do
      end do
      nq_t = nq

end subroutine setup_t_tables

!---------------------------------------------------------------------!
! Product of two Clifford primitives along direction d (1=x,2=y,3=z): !
!   pe  : barycentric exponent,  cen : barycentre                     !
!---------------------------------------------------------------------!
subroutine bary_dir(d, e1, e2, x1, x2, pe, cen)

      integer,          intent(in)  :: d
      double precision, intent(in)  :: e1, e2, x1, x2
      double precision, intent(out) :: pe, cen

      select case (d)
      case (1)
        call bary_exponent_x       (e1, e2, x1-x2, pe)
        call bary_center_toroidal_x(e1, e2, x1, x2, cen)
      case (2)
        call bary_exponent_y       (e1, e2, x1-x2, pe)
        call bary_center_toroidal_y(e1, e2, x1, x2, cen)
      case default
        call bary_exponent_z       (e1, e2, x1-x2, pe)
        call bary_center_toroidal_z(e1, e2, x1, x2, cen)
      end select

end subroutine bary_dir

!---------------------------------------------------------------------!
! bary_dir for the periodic directions d = 1..ndp                     !
!---------------------------------------------------------------------!
subroutine clifford_pair(ndp, e1, e2, r1, r2, pe, cen)

      integer,          intent(in)  :: ndp
      double precision, intent(in)  :: e1, e2, r1(3), r2(3)
      double precision, intent(out) :: pe(ndp), cen(ndp)
      integer :: d

      do d = 1, ndp
        call bary_dir(d, e1, e2, r1(d), r2(d), pe(d), cen(d))
      end do

end subroutine clifford_pair

!---------------------------------------------------------------------!
! McMurchie-Davidson Hermite coefficients (real Gaussian, 1D):         !
!   (y-A)^i (y-B)^j exp(-p (y-P)^2) = sum_t E(t) d^t/dP^t exp(-p(y-P)^2)
!   PA = P - A, PB = P - B ; the prefactor exp(-mu AB^2) is NOT included
!   E(t) for t = 0 .. i+j                                              !
!---------------------------------------------------------------------!
subroutine hermite_E_1D(i, j, p, PA, PB, E)

      integer,          intent(in)  :: i, j
      double precision, intent(in)  :: p, PA, PB
      double precision, intent(out) :: E(0:NM)

      double precision :: Enew(0:NM+1), Eold(-1:NM+1), inv2p
      integer :: m, t, deg

      inv2p = 0.5d0 / p
      Eold  = 0.d0
      Eold(0) = 1.d0
      deg = 0
      do m = 1, i + j
        do t = 0, deg + 1
          if (m <= i) then
            Enew(t) = inv2p*Eold(t-1) + PA*Eold(t) + dble(t+1)*Eold(t+1)
          else
            Enew(t) = inv2p*Eold(t-1) + PB*Eold(t) + dble(t+1)*Eold(t+1)
          end if
        end do
        deg = deg + 1
        Eold(0:deg) = Enew(0:deg)
        Eold(deg+1) = 0.d0
      end do
      E = 0.d0
      E(0:deg) = Eold(0:deg)

end subroutine hermite_E_1D

!---------------------------------------------------------------------!
! Cartesian label -> (lx,ly,lz) by counting x, y, z letters:          !
!   's' -> 000, 'px' -> 100, 'dxy' -> 110, 'dzz' -> 002,              !
!   'fxyz' -> 111, 'gxxyy' -> 220 ...                                 !
! ADAPT if your AO%orbital uses another convention.                   !
!---------------------------------------------------------------------!
subroutine orbital_to_lxyz(label, lxyz)

      character(len=*), intent(in)  :: label
      integer,          intent(out) :: lxyz(3)
      integer :: c
      character :: ch

      lxyz = 0
      do c = 2, len_trim(label)              ! skip shell letter s,p,d,f,g
        ch = label(c:c)
        select case (ch)
        case ('x','X') ; lxyz(1) = lxyz(1) + 1
        case ('y','Y') ; lxyz(2) = lxyz(2) + 1
        case ('z','Z') ; lxyz(3) = lxyz(3) + 1
        end select
      end do

end subroutine orbital_to_lxyz

end module clifford_tools
