      !===================================================================!
      !                                                                   !
      ! module clifford_tools : shared tools for the Clifford-Gaussian    !
      ! integrals in 1D, 2D and 3D                                        !
      !                                                                   !
      ! ndp = number of periodic directions:                              !
      !   Tori1D : ndp = 1  (x periodic,   y,z real)                      !
      !   Tori2D : ndp = 2  (x,y periodic, z real)                        !
      !   Tori3D : ndp = 3  (x,y,z periodic)                              !
      !                                                                   !
      ! Contents                                                          !
      ! --------                                                          !
      ! init_clifford             binomials, double factorials and the    !
      !                           derivative table T (call once)          !
      ! clifford_bessel_ratios    R(e,nu) = exp(-A) I_nu(A) / A^e         !
      ! Icliff_from_R             exp(-A) Int sin^n cos^m exp(A cos) dx   !
      ! overlap_1D_clifford(_ext) 1D overlap factor, periodic direction   !
      ! kinetic_1D_clifford       1D Laplacian factor, periodic direction !
      ! overlap_1D_real           1D overlap factor, real direction       !
      ! kinetic_1D_real           1D Laplacian factor, real direction     !
      ! bessel_scaled_vector      exp(-A) I_n(A), n = 0..N (Miller)       !
      ! bessel_scaled_vector_fast same, start at N, no allocation (ERI)   !
      ! bessel_nmax_estimate      truncation order of exp(-A) I_n(A)      !
      ! fourier_sin_product       Fourier coefficients of sin^i sin^l     !
      ! setup_t_tables(ncap,ndp)  t nodes/weights + Bessel tables (NA,ERI)!
      ! bary_dir                  barycentre exponent/centre along d      !
      ! clifford_pair             bary_dir for d = 1..ndp                 !
      ! hermite_E_1D              McMurchie-Davidson E coefficients       !
      ! orbital_to_lxyz           'dxy' -> (1,1,0)                        !
      !                                                                   !
      ! t-tables (filled by setup_t_tables):                              !
      !   nq_t            number of t nodes (64 or 128)                   !
      !   t_node(k)       nodes   of  Int_0^inf f(t) dt                   !
      !   t_weight(k)     weights                                         !
      !   t_bessel(n,k,d) exp(-C) I_n(C),  C = 2 t_k^2 / a_d^2            !
      !   t_nmax(k,d)     last order kept (1e-16 relative)                !
      !   d = 1..ndp  (a_1 = ax, a_2 = ay, a_3 = az)                      !
      !                                                                   !
      !===================================================================!

module clifford_tools

      implicit none

      private
      public :: NMAX_CLIFF, init_clifford
      public :: clifford_bessel_ratios, Icliff_from_R
      public :: overlap_1D_clifford, overlap_1D_clifford_ext, overlap_1D_real
      public :: kinetic_1D_clifford, kinetic_1D_real
      public :: bessel_scaled_vector, bessel_scaled_vector_fast, bessel_nmax_estimate
      public :: pair_vector_cached, pair_cache_free
      public :: fourier_sin_product, hermite_E_1D, orbital_to_lxyz
      public :: setup_t_tables, nq_t, t_node, t_weight, t_bessel, t_nmax, t_bt, t_nmaxall, t_kfirst, t_nsplit
      public :: bary_dir, clifford_pair
      public :: ps_on, ps_noscreen, ps_tol, ps_qp, ps_bra, ps_ket
      public :: rys_req, rys_on, rys_ratio, rys_nodes, rys_build, rys_rule, bessel_column_at, bessel_columns_multi
      public :: rys_nq_used, rys_nq_total, rys_nhi

      ! ---- t-quadrature tables (nuclear attraction and ERIs) ---- !
      integer,          save              :: nq_t = 0
      double precision, save, allocatable :: t_node(:), t_weight(:), t_bessel(:,:,:)
      integer,          save, allocatable :: t_nmax(:,:)
      ! [Nov 2026] transposed table for the ERI kernel: t_bt(k,n,d) = c_n t_bessel(n,k,d)
      ! (c_0 = 1, c_n = 2) for n <= t_nmax(k,d), 0 beyond; t_nmaxall(d) = max_k t_nmax(k,d)
      double precision, save, allocatable :: t_bt(:,:,:)
      integer,          save, allocatable :: t_nmaxall(:)
      integer,          save, allocatable :: t_kfirst(:,:)   ! t_bt(k,n,d) = 0 for k < t_kfirst(n,d)
      integer,          save, allocatable :: t_nsplit(:)     ! last n with >= 8 nodes active

      ! ---- [Nov 2026] keyword prim_screen: primitive-pair Schwarz screening ---- !
      ! ps_qp(i,j,SI,SJ) = max_comp sqrt|(ij|ij)| for primitives i,j of shells   !
      ! SI (unit cell) and SJ; a primitive quartet is skipped in the ERI kernel !
      ! when ps_bra(i,j)*ps_ket(k,l) < ps_tol. Off by default.                  !
      logical,          save :: ps_on = .false., ps_noscreen = .false.
      double precision, save :: ps_tol = 1.d-14
      double precision, save, allocatable, target :: ps_qp(:,:,:,:)
      double precision, pointer :: ps_bra(:,:) => null(), ps_ket(:,:) => null()

      ! ---- [Nov 2026] keyword rys_quad: Rys-type t quadrature per primitive quartet ---- !
      ! u^2 = t^2/(rho+t^2), Gauss rule of weight exp(-T u^2) on [0,1],                   !
      ! T = rho * (chord distance P-Q)^2.  Used only when rho/a_d^2 >= rys_ratio for every !
      ! periodic direction (localized quartets); otherwise the global t tables are used.   !
      ! Rules are tabulated on v = sqrt(T) = 0, 0.1, ..., 12; above, the scaled v = 12 rule.!
      logical,          save :: rys_req = .false.     ! keyword given
      logical,          save :: rys_on  = .false.     ! tables built and self-check passed
      double precision, parameter :: rys_ratio = 10.d0
      integer,          parameter :: rys_nv = 120
      double precision, parameter :: rys_dv = 0.1d0
      integer,          save :: rys_nlo = 0, rys_nhi = -1
      double precision, save, allocatable :: rys_u(:,:,:), rys_w(:,:,:)   ! (node, 0:rys_nv, n)
      integer(8),       save :: rys_nq_used = 0, rys_nq_total = 0          ! statistics (primitive quartets)
      !$omp threadprivate(ps_bra, ps_ket)

      ! highest total power n+m of sin/cos handled (8 = g-g overlap,
      ! 10 = g-g kinetic). Raise if you need h, i functions.
      integer, parameter :: NMAX_CLIFF = 16
      integer, parameter :: NM = NMAX_CLIFF

      double precision, parameter :: pi = 3.14159265358979323846d0

      ! ---- [Nov 2026] memo of the pair Fourier vectors of the ERI kernel ---- !
      !   V_n = sum_k d_k exp(-A) I_|n - s k|(A),  n = 0..N,  s = +1 (ket) / -1 (bra)
      ! d_k = Fourier coefficients of sin^i(u+phA) sin^l(u+phB). Keyed on the  !
      ! bit patterns of all inputs, so a hit returns exactly what would be     !
      ! computed. One table shared by all threads: a slot is filled once, under!
      ! a critical section, and published last (pv_off), so readers need no   !
      ! lock. The pool is reserved once (virtual memory, pages are used only   !
      ! when filled). When it is full, vectors are just computed.             !
      integer,    parameter :: pv_slots = 262144                   ! hash slots (power of 2)
      integer(8), parameter :: pv_pool_max = 16000000_8            ! complex numbers (256 MB max)
      integer(8),       save, allocatable :: pv_key(:,:), pv_off(:)
      complex(8),       save, allocatable :: pv_pool(:)
      integer(8),       save :: pv_used = 0
      logical,          save :: pv_ready = .false.

      ! T(e,nu,s,m): coefficients of d^m/dA^m [A^-s I_s] = sum T A^-e I_nu
      double precision, save :: T(0:NM, 0:NM, 0:NM/2, 0:NM)
      double precision, save :: binom(0:NM, 0:NM)
      double precision, save :: dfact(-1:2*NM)
      logical,          save :: initialized = .false.

contains

!---------------------------------------------------------------------!
subroutine init_clifford()

      integer          :: s, m, e, nu, n, k
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
! Pair Fourier vector of the ERI kernel, memoized (shared by threads). !
!---------------------------------------------------------------------!
subroutine pair_vector_cached(N, A, i, l, phA, phB, sgn, V)

      integer,          intent(in)  :: N, i, l, sgn
      double precision, intent(in)  :: A, phA, phB
      complex(8),       intent(out) :: V(0:N)

      integer(8) :: k1, k2, k3, k4, hh, o
      integer    :: h0, h, probe, m, kk
      double precision :: Iv(0:N+i+l)
      complex(8)       :: dco(-(i+l):(i+l)), acc
      logical          :: ready

      !$omp atomic read
      ready = pv_ready
      if (.not. ready) then
        !$omp critical (pv_cache)
        if (.not. pv_ready) then
          allocate(pv_key(4,0:pv_slots-1), pv_off(0:pv_slots-1), pv_pool(pv_pool_max))
          pv_off = -1 ; pv_used = 0
          !$omp flush
          !$omp atomic write
          pv_ready = .true.
        end if
        !$omp end critical (pv_cache)
      end if

      k1 = transfer(A, k1) ; k2 = transfer(phA, k2) ; k3 = transfer(phB, k3)
      k4 = int(N,8) + 1000000_8*(int(i,8) + 64_8*(int(l,8) + 64_8*int(sgn+1,8)))
      hh = k1*31_8 + k2
      hh = hh*31_8 + k3
      hh = hh*31_8 + k4
      hh = ieor(hh, ishft(hh,-31))
      h0 = int(iand(hh, int(pv_slots-1,8)))

      ! ---- lookup (lock free) ---- !
      h = h0
      do probe = 1 , 64
        !$omp atomic read
        o = pv_off(h)
        if (o < 0) exit
        !$omp flush
        if (pv_key(1,h) == k1 .and. pv_key(2,h) == k2 .and. pv_key(3,h) == k3 .and. pv_key(4,h) == k4) then
          V(0:N) = pv_pool(o:o+N)
          return
        end if
        h = iand(h+1, pv_slots-1)
      end do

      ! ---- compute ---- !
      call bessel_scaled_vector_fast(N+i+l, A, Iv)
      if (i + l == 0) then
        V(0:N) = cmplx(Iv(0:N), 0.d0, kind=8)
      else
        call fourier_sin_product(i, l, phA, phB, dco)
        do m = 0 , N
          acc = (0.d0,0.d0)
          do kk = -(i+l) , i+l
            acc = acc + dco(kk) * Iv(abs(m - sgn*kk))
          end do
          V(m) = acc
        end do
      end if

      ! ---- store (one writer at a time; the slot is published last) ---- !
      !$omp critical (pv_cache)
      h = h0
      do probe = 1 , 64
        if (pv_off(h) < 0) exit
        if (pv_key(1,h) == k1 .and. pv_key(2,h) == k2 .and. pv_key(3,h) == k3 .and. pv_key(4,h) == k4) exit
        h = iand(h+1, pv_slots-1)
      end do
      if (pv_off(h) < 0 .and. pv_used + N + 1 <= pv_pool_max) then
        o = pv_used + 1
        pv_pool(o:o+N) = V(0:N)
        pv_used = pv_used + N + 1
        pv_key(:,h) = (/ k1, k2, k3, k4 /)
        !$omp flush
        !$omp atomic write
        pv_off(h) = o
      end if
      !$omp end critical (pv_cache)

end subroutine pair_vector_cached

!---------------------------------------------------------------------!
! Release the pair-vector memo (call outside parallel regions).        !
!---------------------------------------------------------------------!
subroutine pair_cache_free()
      if (allocated(pv_off))  deallocate(pv_off)
      if (allocated(pv_key))  deallocate(pv_key)
      if (allocated(pv_pool)) deallocate(pv_pool)
      pv_used = 0
      pv_ready = .false.
end subroutine pair_cache_free

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
      if (nq < 1) nq = 64          ! 64/128 (default tables) or 32/40 (keyword fast_quad)

      if (nq_t == nq .and. allocated(t_bessel)) then
        if (ubound(t_bessel,1) >= ncap .and. size(t_bessel,3) >= ndp) return
      end if
      if (allocated(t_node))   deallocate(t_node)
      if (allocated(t_weight)) deallocate(t_weight)
      if (allocated(t_bessel)) deallocate(t_bessel)
      if (allocated(t_nmax))   deallocate(t_nmax)
      if (allocated(t_bt))     deallocate(t_bt)
      if (allocated(t_nmaxall)) deallocate(t_nmaxall)
      if (allocated(t_kfirst)) deallocate(t_kfirst)
      if (allocated(t_nsplit)) deallocate(t_nsplit)

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
      allocate(t_bt(nq, 0:ncap, ndp), t_nmaxall(ndp))
      t_bt = 0.d0
      do d = 1, ndp
        t_nmaxall(d) = maxval(t_nmax(:,d))
        do k = 1, nq
          t_bt(k,0,d) = t_bessel(0,k,d)
          do n = 1, t_nmax(k,d)
            t_bt(k,n,d) = 2.d0*t_bessel(n,k,d)
          end do
        end do
      end do
      allocate(t_kfirst(0:ncap, ndp))
      do d = 1, ndp
        do n = 0, ncap
          t_kfirst(n,d) = nq + 1
          do k = 1, nq
            if (any(t_bt(k:nq,n,d) /= 0.d0)) then
              t_kfirst(n,d) = k ; exit
            end if
          end do
          ! keep only leading zeros: rows k >= t_kfirst are summed as they are
        end do
      end do
      allocate(t_nsplit(ndp))
      do d = 1, ndp
        t_nsplit(d) = 0
        do n = 0, ncap
          if (count(t_nmax(:,d) >= n) >= min(8,nq)) t_nsplit(d) = n
        end do
      end do
      nq_t = nq

      ! [Nov 2026] keyword rys_quad: Rys-type rules (up to (gg|gg), 24 nodes)
      if (rys_req) then
        call rys_build(rys_nodes(16))
        rys_on = .true.
      end if

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
      do c = 2, len_trim(label)
        ch = label(c:c)
        select case (ch)
        case ('x','X') ; lxyz(1) = lxyz(1) + 1
        case ('y','Y') ; lxyz(2) = lxyz(2) + 1
        case ('z','Z') ; lxyz(3) = lxyz(3) + 1
        end select
      end do

end subroutine orbital_to_lxyz


!---------------------------------------------------------------------!
! [Nov 2026] keyword rys_quad                                         !
! number of Rys nodes for a shell quartet of total angular momentum NN !
! (12 + 3/4 NN: 12 for s, 18 for (dd|dd), 24 for (gg|gg))            !
!---------------------------------------------------------------------!
integer function rys_nodes(NN)
      integer, intent(in) :: NN
      rys_nodes = 12 + (3*NN + 3)/4
end function rys_nodes

!---------------------------------------------------------------------!
! Gauss rules of weight exp(-T u^2) on [0,1] for n = 10 .. nhi and     !
! T = (j*rys_dv)^2, j = 0..rys_nv (discretized Stieltjes on 1000 GL    !
! points of [0, min(1, sqrt(60/T))], then the Jacobi matrix, dstev).   !
! Stored weights are w_i * exp(T u_i^2) (they multiply the integrand). !
!---------------------------------------------------------------------!
subroutine rys_build(nhi)
      integer, intent(in) :: nhi
      integer, parameter  :: M = 1000
      double precision    :: xg(M), wg(M), X(M), W(M), p0(M), p1(M), p2(M)
      double precision    :: al(nhi), bt(nhi), dd(nhi), ee(nhi), Z(nhi,nhi), work(max(1,2*nhi-2))
      double precision    :: T, ub, nrm, nrm0, w0
      integer             :: n, j, k, info, i
      double precision, parameter :: pi = 3.14159265358979323846d0
      if (rys_nhi >= nhi .and. allocated(rys_u)) return
      if (allocated(rys_u)) deallocate(rys_u, rys_w)
      rys_nlo = 10 ; rys_nhi = max(nhi, rys_nlo)
      allocate(rys_u(rys_nhi, 0:rys_nv, rys_nlo:rys_nhi), rys_w(rys_nhi, 0:rys_nv, rys_nlo:rys_nhi))
      rys_u = 0.d0 ; rys_w = 0.d0
      call gl_unit(M, xg, wg)
      !$omp parallel do default(shared) private(j,T,ub,X,W,n,al,bt,p0,p1,p2,k,nrm,nrm0,dd,ee,Z,work,info,w0,i) schedule(dynamic,1)
      do j = 0 , rys_nv
        T  = (j*rys_dv)**2
        ub = 1.d0
        if (T > 60.d0) ub = dsqrt(60.d0/T)
        X  = ub*xg
        W  = ub*wg*dexp(-T*X*X)
        w0 = sum(W)
        do n = rys_nlo , rys_nhi
          p0 = 0.d0 ; p1 = 1.d0 ; nrm0 = 1.d0
          do k = 1 , n
            nrm   = sum(W*p1*p1)
            al(k) = sum(W*X*p1*p1)/nrm
            if (k > 1) then
              bt(k) = nrm/nrm0
            else
              bt(k) = 0.d0
            end if
            p2 = (X - al(k))*p1 - bt(k)*p0
            p0 = p1 ; p1 = p2 ; nrm0 = nrm
          end do
          dd(1:n) = al(1:n)
          do k = 1 , n-1
            ee(k) = dsqrt(bt(k+1))
          end do
          call dstev('V', n, dd, ee, Z, nhi, work, info)
          if (info /= 0) error stop "rys_build: dstev failed"
          do i = 1 , n
            rys_u(i,j,n) = dd(i)
            rys_w(i,j,n) = w0*Z(1,i)**2*dexp(T*dd(i)*dd(i))
          end do
        end do
      end do
      !$omp end parallel do
end subroutine rys_build

!---- Gauss-Legendre nodes/weights on [0,1] ----!
subroutine gl_unit(n, x, w)
      integer, intent(in) :: n
      double precision, intent(out) :: x(n), w(n)
      double precision, parameter :: pi = 3.14159265358979323846d0
      integer :: i, j, it
      double precision :: z, q0, q1, q2, dp, dz
      do i = 1, n
        z = cos(pi*(dble(i)-0.25d0)/(dble(n)+0.5d0))
        do it = 1, 100
          q0 = 1.d0 ; q1 = z
          do j = 2, n
            q2 = ((2.d0*j-1.d0)*z*q1 - (j-1.d0)*q0)/dble(j) ; q0 = q1 ; q1 = q2
          end do
          dp = dble(n)*(z*q1 - q0)/(z*z - 1.d0)
          dz = q1/dp ; z = z - dz
          if (abs(dz) < 1.d-15) exit
        end do
        q0 = 1.d0 ; q1 = z
        do j = 2, n
          q2 = ((2.d0*j-1.d0)*z*q1 - (j-1.d0)*q0)/dble(j) ; q0 = q1 ; q1 = q2
        end do
        dp = dble(n)*(z*q1 - q0)/(z*z - 1.d0)
        x(i) = 0.5d0*(1.d0 - z)
        w(i) = 1.d0/((1.d0 - z*z)*dp*dp)
      end do
end subroutine gl_unit

!---------------------------------------------------------------------!
! n-point rule for exp(-T u^2): nearest tabulated v = sqrt(T); above   !
! v = 12 the v = 12 rule scaled (u -> u*12/v, w -> w*12/v).            !
!---------------------------------------------------------------------!
subroutine rys_rule(n, T, u, w)
      integer, intent(in) :: n
      double precision, intent(in)  :: T
      double precision, intent(out) :: u(n), w(n)
      double precision :: v
      integer :: j
      v = dsqrt(max(T,0.d0))
      j = nint(v/rys_dv)
      if (j <= rys_nv) then
        u = rys_u(1:n,j,n) ; w = rys_w(1:n,j,n)
      else
        u = rys_u(1:n,rys_nv,n)*(rys_nv*rys_dv/v)
        w = rys_w(1:n,rys_nv,n)*(rys_nv*rys_dv/v)
      end if
end subroutine rys_rule

!---------------------------------------------------------------------!
! exp(-C) I_n(C), n = 0..ncol, ~1e-16 relative to n = 0 (no calls to  !
! iv_scaled, which is slow):                                          !
!   C <  30 : Miller backward from max(ncol, 9 sqrt C)+10, normalized  !
!             by the sum rule  I~_0 + 2 sum_n I~_n = 1                 !
!   C >= 30 : I~_0, I~_1 from the large-C asymptotic series (exact to  !
!             1e-16 there); forward recurrence when ncol <= 2 sqrt C   !
!             (stable), else Miller from max(ncol, 6.5 sqrt C)+10      !
!             normalized by I~_0                                       !
!---------------------------------------------------------------------!
subroutine bessel_column_at(C, ncol, Iv)
      double precision, intent(in)  :: C
      integer,          intent(in)  :: ncol
      double precision, intent(out) :: Iv(0:ncol)
      double precision, parameter   :: pi = 3.14159265358979323846d0
      integer :: k, ns
      double precision :: i2c, i0, i1, t0, t1, s0, s1, sm
      if (C < 1.d-12) then
        Iv = 0.d0 ; Iv(0) = 1.d0 ; return
      end if
      i2c = 2.d0/C
      if (C < 30.d0) then
        ns = max(ncol, ceiling(9.d0*dsqrt(C))) + 10
        call miller(ns)
        return
      end if
      ! asymptotic series for I~_0 and I~_1
      s0 = 1.d0 ; s1 = 1.d0 ; t0 = 1.d0 ; t1 = 1.d0
      do k = 1 , 40
        t0 = t0*(dble((2*k-1)**2))/(8.d0*k*C)
        t1 = t1*(dble((2*k-1)**2) - 4.d0)/(8.d0*k*C)
        s0 = s0 + t0 ; s1 = s1 + t1
        if (dabs(t0) < 1.d-17 .and. dabs(t1) < 1.d-17) exit
      end do
      i0 = s0/dsqrt(2.d0*pi*C) ; i1 = s1/dsqrt(2.d0*pi*C)
      Iv(0) = i0
      if (ncol == 0) return
      if (dble(ncol) <= 2.d0*dsqrt(C)) then
        Iv(1) = i1
        do k = 1 , ncol-1
          Iv(k+1) = Iv(k-1) - dble(k)*i2c*Iv(k)
        end do
      else
        ns = max(ncol, ceiling(6.5d0*dsqrt(C))) + 10
        call miller(ns)
      end if
contains
      subroutine miller(ns)
        integer, intent(in) :: ns
        double precision    :: wk(0:ns+1)
        integer             :: kk
        wk(ns+1) = 0.d0 ; wk(ns) = 1.d-280
        do kk = ns , 1 , -1
          wk(kk-1) = wk(kk+1) + dble(kk)*i2c*wk(kk)
          if (wk(kk-1) > 1.d250) wk(kk-1:ns) = wk(kk-1:ns)*1.d-250
        end do
        if (C < 30.d0) then
          sm = wk(0) + 2.d0*sum(wk(1:ns))
          Iv(0:ncol) = wk(0:ncol)/sm
        else
          Iv(0:ncol) = wk(0:ncol)*(i0/wk(0))
        end if
      end subroutine miller
end subroutine bessel_column_at


!---------------------------------------------------------------------!
! Same as bessel_column_at for nr arguments at once (node j: C(j),     !
! orders 0..ncol(j) into col(0:ncol(j),j)). The recurrences of the     !
! different nodes run side by side (vectorized over j), which hides    !
! the latency of the single recurrence.                                !
!---------------------------------------------------------------------!
subroutine bessel_columns_multi(nr, C, ncol, ld, col)
      integer,          intent(in)  :: nr, ld
      double precision, intent(in)  :: C(nr)
      integer,          intent(in)  :: ncol(nr)
      double precision, intent(out) :: col(0:ld, nr)
      double precision, parameter   :: pi = 3.14159265358979323846d0
      integer          :: j, k, nf, nmil, kmax, jj
      integer          :: jf(nr), jm(nr), ns(nr), nc(nr)
      double precision :: i2c(nr), i0(nr), i1(nr), t0, t1, s0, s1, r, sm
      double precision, allocatable :: W(:,:)
      nf = 0 ; nmil = 0
      do j = 1 , nr
        if (C(j) < 1.d-12) then
          col(0:ncol(j),j) = 0.d0 ; col(0,j) = 1.d0
          cycle
        end if
        i2c(j) = 2.d0/C(j)
        if (C(j) >= 30.d0) then
          s0 = 1.d0 ; s1 = 1.d0 ; t0 = 1.d0 ; t1 = 1.d0 ; r = 1.d0/(8.d0*C(j))
          do k = 1 , 40
            t0 = t0*(dble((2*k-1)**2))*r/dble(k)
            t1 = t1*(dble((2*k-1)**2) - 4.d0)*r/dble(k)
            s0 = s0 + t0 ; s1 = s1 + t1
            if (dabs(t0) < 1.d-17 .and. dabs(t1) < 1.d-17) exit
          end do
          i0(j) = s0/dsqrt(2.d0*pi*C(j)) ; i1(j) = s1/dsqrt(2.d0*pi*C(j))
          if (ncol(j) == 0) then
            col(0,j) = i0(j)
          else if (dble(ncol(j)) <= 2.d0*dsqrt(C(j))) then
            nf = nf + 1 ; jf(nf) = j
          else
            nmil = nmil + 1 ; jm(nmil) = j ; ns(nmil) = max(ncol(j), ceiling(6.5d0*dsqrt(C(j)))) + 10
          end if
        else
          nmil = nmil + 1 ; jm(nmil) = j ; ns(nmil) = max(ncol(j), ceiling(9.d0*dsqrt(C(j)))) + 10
        end if
      end do
      ! ---- forward group ---- !
      if (nf > 0) then
        kmax = 0
        do jj = 1 , nf
          nc(jj) = ncol(jf(jj)) ; kmax = max(kmax, nc(jj))
        end do
        allocate(W(nf, 0:kmax))
        do jj = 1 , nf
          W(jj,0) = i0(jf(jj)) ; W(jj,1) = i1(jf(jj)) ; r = i2c(jf(jj)) ; i2c(jf(jj)) = r
        end do
        block
          double precision :: q(nf)
          do jj = 1 , nf
            q(jj) = i2c(jf(jj))
          end do
          do k = 1 , kmax-1
            do jj = 1 , nf
              W(jj,k+1) = W(jj,k-1) - dble(k)*q(jj)*W(jj,k)
            end do
          end do
        end block
        do jj = 1 , nf
          col(0:nc(jj),jf(jj)) = W(jj,0:nc(jj))
        end do
        deallocate(W)
      end if
      ! ---- Miller group ---- !
      if (nmil > 0) then
        kmax = maxval(ns(1:nmil))
        allocate(W(nmil, 0:kmax+1))
        W = 0.d0
        block
          double precision :: q(nmil), v
          do jj = 1 , nmil
            q(jj) = i2c(jm(jj)) ; W(jj,ns(jj)) = 1.d-280
          end do
          do k = kmax , 1 , -1
            do jj = 1 , nmil
              v = W(jj,k+1) + dble(k)*q(jj)*W(jj,k)
              W(jj,k-1) = merge(v, W(jj,k-1), k <= ns(jj))
            end do
            if (mod(k,8) == 0 .or. k == 1) then
              do jj = 1 , nmil
                if (W(jj,k-1) > 1.d200) W(jj,k-1:ns(jj)) = W(jj,k-1:ns(jj))*1.d-200
              end do
            end if
          end do
          do jj = 1 , nmil
            j = jm(jj)
            if (C(j) < 30.d0) then
              sm = W(jj,0) + 2.d0*sum(W(jj,1:ns(jj)))
              col(0:ncol(j),j) = W(jj,0:ncol(j))/sm
            else
              col(0:ncol(j),j) = W(jj,0:ncol(j))*(i0(j)/W(jj,0))
            end if
          end do
        end block
        deallocate(W)
      end if
end subroutine bessel_columns_multi

end module clifford_tools
