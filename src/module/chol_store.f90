!=====================================================================!
!                                                                     !
!   module chol_store  [Nov 2026]  keyword cholesky                   !
!                                                                     !
!   Cholesky decomposition of the k-space two-electron integrals for  !
!   the fast k-CCSD(T) (fast_ccsdt / fast_ccsdt_real).                !
!                                                                     !
!   g(p,q,r,s; k1,k2,k3) = (p k1 q k2 | r k3 s k4) , k4 = k1-k2+k3    !
!                        = sum_Q X^Q_pq(k1,k2) conj( X^Q_sr(k4,k3) )   !
!                                                                     !
!   The pair (p k1, q k2) carries the momentum K = k2 - k1. For every !
!   K the matrix V(K)_{(pq),(sr)} = (pq|rs) is Hermitian and positive !
!   semidefinite (Coulomb self-repulsion of the pair densities), so a !
!   pivoted Cholesky decomposition V(K) = L L^H exists. It is done in !
!   the Bloch AO pair basis (a, b, k2) from the real-space SAO        !
!   integrals ERI_AO(a, j2, j3, j4) (a in the unit cell), with the    !
!   same phases and 1/N as transform_ERI_AO_to_kMO_half, and then     !
!   transformed to the k-MO basis:                                    !
!                                                                     !
!     X^Q_pq(k1,k2) = sum_ab conj(C_ap(k1)) C_bq(k2) L_{(a,b,k2),Q}    !
!                                                                     !
!   Stored: LK(K)%x(Q, p, q, k2)  (k1 = k2 - K),  Q = 1..LK(K)%m      !
!                                                                     !
!   Error: every integral of V(K) is reproduced to <= chol_tol        !
!   (largest remaining diagonal; Cauchy-Schwarz on the residual).     !
!                                                                     !
!   Use in the CC:                                                    !
!     chol_build_g   : the compact g with an OCCUPIED first index only !
!                      (g(1:o_hi,:,:,:,slot)), all other CC terms     !
!     chol_vvvv_gm   : Wvvvv(t1) blocks for the particle-particle     !
!                      ladder, built from X on the fly (never stored) !
!=====================================================================!

module chol_store

      use constants_module, only: dp, dpc, pi
      implicit none

      private
      public :: chol_on, chol_tol, chol_nb, LK, kdiff, chol_rank_total
      public :: chol_build, chol_reduce_active, chol_build_g, chol_bytes
      public :: chol_vvvv_gm, chol_vvvv_gm_r, chol_free

      type lk_t
        integer :: m = 0
        complex(dpc), allocatable :: x(:,:,:,:)      ! x(Q, p, q, k2)
      end type lk_t

      logical , save :: chol_on  = .false.
      real(dp), save :: chol_tol = 1.0e-6_dp
      integer , save :: chol_nb  = 0                  ! bands stored in x (nfuc or active)
      type(lk_t), allocatable, save :: LK(:)
      integer   , allocatable, save :: kdiff(:,:)     ! kdiff(k2,k1) = index of k2 - k1

contains

      !-----------------------------------------------------------------!
      ! cell coordinates (same order as the transform) and k2 - k1       !
      !-----------------------------------------------------------------!
      subroutine set_coords(coord)
        use unitcell_module, only: nx, ny, nz, N_cell
        integer, intent(out) :: coord(3,N_cell)
        integer :: ix, iy, iz, idx, k1, k2
        do iz = 0, nz-1
          do iy = 0, ny-1
            do ix = 0, nx-1
              idx = iz*nx*ny + iy*nx + ix + 1
              coord(:,idx) = [ix,iy,iz]
            end do
          end do
        end do
        if (allocated(kdiff)) deallocate(kdiff)
        allocate(kdiff(N_cell,N_cell))
        do k1 = 1, N_cell
          do k2 = 1, N_cell
            kdiff(k2,k1) = modulo(coord(3,k2)-coord(3,k1),nz)*nx*ny &
                         + modulo(coord(2,k2)-coord(2,k1),ny)*nx    &
                         + modulo(coord(1,k2)-coord(1,k1),nx) + 1
          end do
        end do
      end subroutine set_coords

      !-----------------------------------------------------------------!
      ! integral-matrix Cholesky per momentum transfer K, then AO -> MO  !
      !-----------------------------------------------------------------!
      subroutine chol_build(nBas, ERI_AO, coeff_SAO, tol, iout)

        use unitcell_module, only: nfuc, nx, ny, nz, N_cell
        use omp_lib

        integer     , intent(in) :: nBas, iout
        real(dp)    , intent(in) :: ERI_AO(nfuc,nBas,nBas,nBas)
        complex(dpc), intent(in) :: coeff_SAO(nBas,nBas)
        real(dp)    , intent(in) :: tol

        integer :: f, ncl, n, iK, k1, k2, k3, k4, a, b, c, d, r2, r3, r4, j2, j3, j4, P, m, cap, Q, off1, off2
        integer :: coord(3,N_cell), piv, mK
        real(dp) :: asym
        complex(dpc), allocatable :: Xs(:,:,:,:), Xd(:,:,:,:)
        real(dp) :: theta, invN, dmax0, dmin_seen, t0
        complex(dpc) :: ph(N_cell,N_cell), w
        real(dp)    , allocatable :: diag(:,:), dg(:)
        complex(dpc), allocatable :: T3(:,:,:), G2(:,:), col(:), Z(:,:), Lao(:,:), tmpL(:,:), Lm(:,:), tmp(:,:), Xm(:,:)
        complex(dpc), parameter :: one=(1.0_dp,0.0_dp), zero=(0.0_dp,0.0_dp)

        f = nfuc ; ncl = N_cell ; n = f*f*ncl
        chol_tol = tol
        chol_nb  = f
        invN = 1.0_dp/real(ncl,dp)
        t0 = omp_get_wtime()

        call set_coords(coord)
        do k1 = 1, ncl
          do r2 = 1, ncl
            theta = 2.0_dp*pi * ( &
                real(coord(1,k1),dp)*real(coord(1,r2),dp)/real(nx,dp) &
              + real(coord(2,k1),dp)*real(coord(2,r2),dp)/real(ny,dp) &
              + real(coord(3,k1),dp)*real(coord(3,r2),dp)/real(nz,dp))
            ph(k1,r2) = exp(cmplx(0.0_dp,theta,kind=dpc))        ! e^{+i k.R}
          end do
        end do

        ! ---- diagonals of every V(K):  V_{(a,b,k2),(a,b,k2)} ---- !
        !   = 1/N sum_{R2,R3,R4} ERI(a,bR2,bR3,aR4) e^{i k2 R2} e^{-i k2 R3} e^{i k1 R4} , k1 = k2-K
        allocate(diag(n,ncl))
        !$omp parallel default(shared) private(a,b,r2,r3,r4,k2,iK,k1,T3,G2,w,P)
        allocate(T3(ncl,ncl,ncl), G2(ncl,ncl))
        !$omp do collapse(2) schedule(dynamic)
        do b = 1, f
          do a = 1, f
            do r4 = 1, ncl
              do r3 = 1, ncl
                do r2 = 1, ncl
                  T3(r2,r3,r4) = ERI_AO(a, b+(r2-1)*f, b+(r3-1)*f, a+(r4-1)*f)
                end do
              end do
            end do
            ! G2(k2,R4) = sum_{R2,R3} T3 e^{i k2 (R2-R3)}
            do r4 = 1, ncl
              do k2 = 1, ncl
                w = zero
                do r3 = 1, ncl
                  do r2 = 1, ncl
                    w = w + T3(r2,r3,r4)*ph(k2,r2)*conjg(ph(k2,r3))
                  end do
                end do
                G2(k2,r4) = w
              end do
            end do
            do iK = 1, ncl
              do k2 = 1, ncl
                k1 = kdiff(k2,iK)                       ! k1 = k2 - K
                w = zero
                do r4 = 1, ncl
                  w = w + G2(k2,r4)*ph(k1,r4)
                end do
                P = a + (b-1)*f + (k2-1)*f*f
                diag(P,iK) = real(w,dp)*invN
              end do
            end do
          end do
        end do
        !$omp end do
        deallocate(T3, G2)
        !$omp end parallel

        if (allocated(LK)) call chol_free()
        allocate(LK(ncl))
        allocate(col(n), Z(f,nBas), dg(n), Lm(f,f), tmp(f,f), Xm(f,f))
        dmin_seen = minval(diag)

        do iK = 1, ncl

          ! exchange symmetry (pq|rs) = (rs|pq): V(-K) follows from V(K), only
          ! one K of every pair {K,-K} is decomposed (the other is set below)
          if (kdiff(1,iK) < iK) cycle

          ! ---------------- pivoted Cholesky of V(K) ---------------- !
          dg = diag(:,iK)
          dmax0 = maxval(dg)
          cap = max(16, min(n, 8*f))
          allocate(Lao(n,cap))
          m = 0
          do
            piv = maxloc(dg,1)
            if (dg(piv) <= tol .or. m >= n) exit
            ! pivot pair (d, c, k3): first band d at k4 = k3 - K, second band c at k3
            d  = mod(piv-1, f) + 1
            c  = mod((piv-1)/f, f) + 1
            k3 = (piv-1)/(f*f) + 1
            k4 = kdiff(k3,iK)
            ! Z(a, j2) = sum_{R3,R4} ERI(a, j2, c R3, d R4) e^{-i k3 R3} e^{+i k4 R4}
            !$omp parallel do schedule(static) private(j2,r3,r4,j3,j4,w)
            do j2 = 1, nBas
              Z(:,j2) = zero
              do r4 = 1, ncl
                j4 = d + (r4-1)*f
                do r3 = 1, ncl
                  j3 = c + (r3-1)*f
                  w  = conjg(ph(k3,r3))*ph(k4,r4)
                  Z(:,j2) = Z(:,j2) + ERI_AO(:,j2,j3,j4)*w
                end do
              end do
            end do
            !$omp end parallel do
            ! col(a,b,k2) = 1/N sum_R2 Z(a, b R2) e^{i k2 R2}
            !$omp parallel do collapse(2) schedule(static) private(k2,b,r2,P)
            do k2 = 1, ncl
              do b = 1, f
                P = (b-1)*f + (k2-1)*f*f
                col(P+1:P+f) = zero
                do r2 = 1, ncl
                  col(P+1:P+f) = col(P+1:P+f) + Z(:, b+(r2-1)*f)*ph(k2,r2)
                end do
                col(P+1:P+f) = col(P+1:P+f)*invN
              end do
            end do
            !$omp end parallel do
            ! col -= L(:,1:m) conj(L(piv,1:m))
            if (m > 0) then
              allocate(tmpL(m,1))
              tmpL(:,1) = conjg(Lao(piv,1:m))
              call zgemv('N', n, m, -one, Lao, n, tmpL(:,1), 1, one, col, 1)
              deallocate(tmpL)
            end if
            if (m == cap) call grow(Lao, cap)
            m = m + 1
            Lao(:,m) = col / sqrt(dg(piv))
            dg = dg - real(Lao(:,m)*conjg(Lao(:,m)),dp)
            dg(piv) = 0.0_dp
          end do

          ! ---------------- AO pairs -> MO pairs ---------------- !
          LK(iK)%m = m
          allocate(LK(iK)%x(max(m,1), f, f, ncl))
          LK(iK)%x = zero
          do k2 = 1, ncl
            k1   = kdiff(k2,iK)
            off1 = (k1-1)*f
            off2 = (k2-1)*f
            do Q = 1, m
              P = (k2-1)*f*f
              Lm = reshape(Lao(P+1:P+f*f,Q), [f,f])            ! Lm(a,b)
              call zgemm('N','N', f, f, f, one, Lm, f, coeff_SAO(off2+1,off2+1), nBas, zero, tmp, f)
              call zgemm('C','N', f, f, f, one, coeff_SAO(off1+1,off1+1), nBas, tmp, f, zero, Xm, f)
              LK(iK)%x(Q,:,:,k2) = Xm
            end do
          end do
          deallocate(Lao)

          write(iout,'(a,i5,a,i7,a,i9,a,es10.2,a,es10.2)') '  [cholesky] K =', iK, '   rank =', m, &
               ' of ', n, '   max diag =', dmax0, '   left =', maxval(dg)
        end do

        deallocate(col, Z, dg, Lm, tmp, Xm, diag)

        ! ---- -K from K:  X^Q_rs(kr,ks; -K) = conj X^Q_sr(ks,kr; K) ---- !
        ! (this keeps (pq|rs) = (rs|pq) exact); a self-conjugate K (K = -K)
        ! whose vectors break this symmetry is symmetrized with the doubled set
        ! [X, conj-swap(X)]/sqrt(2) (V = (L L^H + L' L'^H)/2, same accuracy)
        do iK = 1, ncl
          mK = kdiff(1,iK)
          if (mK < iK) then
            m = LK(mK)%m
            LK(iK)%m = m
            allocate(LK(iK)%x(max(m,1), f, f, ncl))
            call conj_swap(LK(mK)%x, LK(iK)%x, iK)
            write(iout,'(a,i5,a,i7,a,i5)') '  [cholesky] K =', iK, '   rank =', m, '   from -K =', mK
          else if (mK == iK .and. LK(iK)%m > 0) then
            m = LK(iK)%m
            allocate(Xs(m, f, f, ncl))
            call conj_swap(LK(iK)%x, Xs, iK)
            asym = maxval(abs(Xs - LK(iK)%x(1:m,:,:,:)))
            if (asym > 1.0e-10_dp) then
              allocate(Xd(2*m, f, f, ncl))
              Xd(1:m,:,:,:)     = LK(iK)%x(1:m,:,:,:) / sqrt(2.0_dp)
              Xd(m+1:2*m,:,:,:) = Xs / sqrt(2.0_dp)
              call move_alloc(Xd, LK(iK)%x)
              LK(iK)%m = 2*m
              write(iout,'(a,i5,a,es9.2,a,i7)') '  [cholesky] K =', iK, '   self-conjugate, asymmetry', asym, &
                   ' -> symmetrized rank', 2*m
            end if
            deallocate(Xs)
          end if
        end do

        write(iout,'(a,es9.2,a,i9,a,f8.2,a,f10.3,a)') '  [cholesky] tol =', tol, '   total rank =', &
             chol_rank_total(), '   vectors:', chol_bytes()/1024.d0**3, ' GB   time =', omp_get_wtime()-t0, ' s'
        if (dmin_seen < -10.0_dp*tol) write(iout,'(a,es10.2)') &
             '  [cholesky] WARNING: negative diagonal integral ', dmin_seen
        flush(iout)

      contains

        ! Y(Q,r,s,ks) = conj X(Q,s,r,kr) , kr = ks - K'(of Y) , K' = target K index
        subroutine conj_swap(X, Y, iKy)
          complex(dpc), intent(in)  :: X(:,:,:,:)
          complex(dpc), intent(out) :: Y(:,:,:,:)
          integer, intent(in) :: iKy
          integer :: ks, kr, r, s, mm
          mm = min(size(X,1), size(Y,1))
          Y = zero
          do ks = 1, ncl
            kr = kdiff(ks,iKy)
            do s = 1, f
              do r = 1, f
                Y(1:mm,r,s,ks) = conjg(X(1:mm,s,r,kr))
              end do
            end do
          end do
        end subroutine conj_swap

        subroutine grow(A, cp)
          complex(dpc), allocatable, intent(inout) :: A(:,:)
          integer, intent(inout) :: cp
          complex(dpc), allocatable :: B(:,:)
          integer :: nc
          nc = min(size(A,1), 2*cp)
          allocate(B(size(A,1), nc))
          B(:,1:cp) = A(:,1:cp)
          call move_alloc(B, A)
          cp = nc
        end subroutine grow

      end subroutine chol_build

      integer function chol_rank_total()
        integer :: i
        chol_rank_total = 0
        if (.not. allocated(LK)) return
        do i = 1, size(LK)
          chol_rank_total = chol_rank_total + LK(i)%m
        end do
      end function chol_rank_total

      double precision function chol_bytes()
        integer :: i
        chol_bytes = 0.d0
        if (.not. allocated(LK)) return
        do i = 1, size(LK)
          if (allocated(LK(i)%x)) chol_bytes = chol_bytes + 16.d0*dble(size(LK(i)%x,kind=8))
        end do
      end function chol_bytes

      subroutine chol_free()
        integer :: i
        if (.not. allocated(LK)) return
        do i = 1, size(LK)
          if (allocated(LK(i)%x)) deallocate(LK(i)%x)
        end do
        deallocate(LK)
      end subroutine chol_free

      !-----------------------------------------------------------------!
      ! frozen core: drop the nfc lowest bands of every k               !
      !-----------------------------------------------------------------!
      subroutine chol_reduce_active(nfc)
        integer, intent(in) :: nfc
        integer :: i, nb, m, C
        complex(dpc), allocatable :: y(:,:,:,:)
        if (nfc <= 0) return
        nb = chol_nb - nfc
        do i = 1, size(LK)
          m = size(LK(i)%x,1) ; C = size(LK(i)%x,4)
          allocate(y(m, nb, nb, C))
          y = LK(i)%x(:, nfc+1:chol_nb, nfc+1:chol_nb, :)
          call move_alloc(y, LK(i)%x)
        end do
        chol_nb = nb
      end subroutine chol_reduce_active

      !-----------------------------------------------------------------!
      ! compact g with an occupied first index:                         !
      !   g(i,q,r,s,slot) = sum_Q X^Q_iq(k1,k2) conj(X^Q_sr(k4,k3))      !
      !-----------------------------------------------------------------!
      subroutine chol_build_g(n1, g)
        use unitcell_module, only: N_cell
        use cc_dims        , only: gidx, nslot
        integer     , intent(in)  :: n1
        complex(dpc), intent(out) :: g(n1,chol_nb,chol_nb,chol_nb,nslot)
        integer :: k1, k2, k3, k4, iK, ic, m, nb, r, s
        integer :: fourth_k
        complex(dpc), allocatable :: A(:,:,:), B(:,:,:)
        complex(dpc), parameter :: one=(1.0_dp,0.0_dp), zero=(0.0_dp,0.0_dp)
        nb = chol_nb
        do k3 = 1, N_cell
          do k2 = 1, N_cell
            do k1 = 1, N_cell
              ic = gidx(k1,k2,k3)
              if (ic <= 0) cycle
              k4 = fourth_k(k1,k3,k2)              ! k4 = k1 - k2 + k3 (fourth_k(a,b,c) = a+b-c)
              iK = kdiff(k2,k1)
              if (kdiff(k3,k4) /= iK) error stop 'chol_build_g: momentum mismatch'
              m  = LK(iK)%m
              if (m == 0) then
                g(:,:,:,:,ic) = zero
                cycle
              end if
              allocate(A(m,n1,nb), B(m,nb,nb))
              A = LK(iK)%x(1:m, 1:n1, :, k2)
              do s = 1, nb
                do r = 1, nb
                  B(:,r,s) = conjg(LK(iK)%x(1:m, s, r, k3))
                end do
              end do
              call zgemm('T','N', n1*nb, nb*nb, m, one, A, m, B, m, zero, g(1,1,1,1,ic), n1*nb)
              deallocate(A, B)
            end do
          end do
        end do
      end subroutine chol_build_g

      !-----------------------------------------------------------------!
      ! Wvvvv(t1) for the ladder (same as g after update_vvvv_in_g):     !
      !   W(a,c,b,d; ka,kc,kb) = (ac|bd) - sum_l t1(l,b) (ac|ld)         !
      !                                  - sum_l t1(l,a) (lc|bd)         !
      !    = sum_Q Xt_ac Y_bd - X_ac dY_bd                               !
      !      Xt_ac = X_ac - sum_l t1(l,a;ka) X_lc     (pair ka,kc)       !
      !      Y_bd  = conj X_db(kd,kb), dY_bd = sum_l t1(l,b;kb) conj X_dl !
      ! Output T(a,c,b,d) for b = b0..b1 (band ranges v_lo..chol_nb)     !
      !-----------------------------------------------------------------!
      subroutine vvvv_core(ka, kc, kb, b0, b1, no, v_lo, Dm, T)
        integer     , intent(in)  :: ka, kc, kb, b0, b1, no, v_lo
        complex(dpc), intent(in)  :: Dm(no, v_lo:chol_nb, *)       ! masked t1(l,a;k)
        complex(dpc), intent(out) :: T(v_lo:chol_nb, v_lo:chol_nb, b0:b1, v_lo:chol_nb)
        integer :: iK, m, nv, nbb, c, d, nf
        complex(dpc), allocatable :: A2(:,:,:), B2(:,:,:), Xo(:,:), Yo(:,:)
        complex(dpc), parameter :: one=(1.0_dp,0.0_dp), zero=(0.0_dp,0.0_dp)
        integer :: fourth_k, kd

        nf  = chol_nb
        nv  = nf - v_lo + 1
        nbb = b1 - b0 + 1
        iK  = kdiff(kc,ka)
        kd  = fourth_k(ka,kb,kc)              ! kd = ka - kc + kb
        m   = LK(iK)%m
        if (m == 0) then
          T = zero
          return
        end if
        allocate(A2(2*m, v_lo:nf, v_lo:nf), B2(2*m, b0:b1, v_lo:nf))
        ! A2(1:m,a,c) = Xt_ac , A2(m+1:2m,a,c) = -X_ac
        A2(1:m, :, :)     =  LK(iK)%x(1:m, v_lo:nf, v_lo:nf, kc)
        A2(m+1:2*m, :, :) = -LK(iK)%x(1:m, v_lo:nf, v_lo:nf, kc)
        if (no > 0) then
          allocate(Xo(m,no))
          do c = v_lo, nf
            Xo = LK(iK)%x(1:m, 1:no, c, kc)
            ! A2(1:m, a, c) -= sum_l Xo(Q,l) Dm(l,a;ka)
            call zgemm('N','N', m, nv, no, -one, Xo, m, Dm(1,v_lo,ka), no, one, A2(1,v_lo,c), 2*m)
          end do
          deallocate(Xo)
        end if
        ! B2(1:m,b,d) = conj X_db(kd,kb) , B2(m+1:2m,b,d) = sum_l Dm(l,b;kb) conj X_dl(kd,kb)
        allocate(Yo(m,no))
        do d = v_lo, nf
          B2(1:m, :, d) = conjg(LK(iK)%x(1:m, d, b0:b1, kb))
          if (no > 0) then
            Yo = conjg(LK(iK)%x(1:m, d, 1:no, kb))
            call zgemm('N','N', m, nbb, no, one, Yo, m, Dm(1,b0,kb), no, zero, B2(m+1,b0,d), 2*m)
          else
            B2(m+1:2*m, :, d) = zero
          end if
        end do
        deallocate(Yo)
        call zgemm('T','N', nv*nv, nbb*nv, 2*m, one, A2, 2*m, B2, 2*m, zero, T, nv*nv)
        deallocate(A2, B2)
      end subroutine vvvv_core

      ! complex: Gm(c,d,a,b-b0+1) = W(a,c,b,d; ka,kc,kb)
      subroutine chol_vvvv_gm(ka, kc, kb, b0, b1, no, v_lo, occ, t1, Gm)
        use unitcell_module, only: nfuc, N_cell
        integer     , intent(in)  :: ka, kc, kb, b0, b1, no, v_lo
        logical     , intent(in)  :: occ(nfuc,N_cell)
        complex(dpc), intent(in)  :: t1(nfuc,nfuc,N_cell)
        complex(dpc), intent(inout) :: Gm(v_lo:nfuc,v_lo:nfuc,v_lo:nfuc,*)
        complex(dpc), allocatable :: D(:,:,:), T(:,:,:,:)
        integer :: a, b, c, k, l
        if (chol_nb /= nfuc) error stop 'chol_vvvv_gm: band count mismatch'
        allocate(D(no, v_lo:nfuc, N_cell))
        do k = 1, N_cell
          do a = v_lo, nfuc
            do l = 1, no
              D(l,a,k) = (0.0_dp,0.0_dp)
              if (occ(l,k) .and. .not.occ(a,k)) D(l,a,k) = t1(l,a,k)
            end do
          end do
        end do
        allocate(T(v_lo:nfuc, v_lo:nfuc, b0:b1, v_lo:nfuc))
        call vvvv_core(ka, kc, kb, b0, b1, no, v_lo, D, T)
        !$omp parallel do collapse(2) private(b,a,c)
        do b = b0, b1
          do a = v_lo, nfuc
            do c = v_lo, nfuc
              Gm(c,:,a,b-b0+1) = T(a,c,b,:)
            end do
          end do
        end do
        !$omp end parallel do
        deallocate(D, T)
      end subroutine chol_vvvv_gm

      ! real (fast_ccsdt_real: every k = -k, real gauge): real part
      subroutine chol_vvvv_gm_r(ka, kc, kb, b0, b1, no, v_lo, occ, t1, Gm)
        use unitcell_module, only: nfuc, N_cell
        integer , intent(in)  :: ka, kc, kb, b0, b1, no, v_lo
        logical , intent(in)  :: occ(nfuc,N_cell)
        real(dp), intent(in)  :: t1(nfuc,nfuc,N_cell)
        real(dp), intent(inout) :: Gm(v_lo:nfuc,v_lo:nfuc,v_lo:nfuc,*)
        complex(dpc), allocatable :: D(:,:,:), T(:,:,:,:)
        integer :: a, b, c, k, l
        if (chol_nb /= nfuc) error stop 'chol_vvvv_gm_r: band count mismatch'
        allocate(D(no, v_lo:nfuc, N_cell))
        do k = 1, N_cell
          do a = v_lo, nfuc
            do l = 1, no
              D(l,a,k) = (0.0_dp,0.0_dp)
              if (occ(l,k) .and. .not.occ(a,k)) D(l,a,k) = cmplx(t1(l,a,k), 0.0_dp, kind=dpc)
            end do
          end do
        end do
        allocate(T(v_lo:nfuc, v_lo:nfuc, b0:b1, v_lo:nfuc))
        call vvvv_core(ka, kc, kb, b0, b1, no, v_lo, D, T)
        !$omp parallel do collapse(2) private(b,a,c)
        do b = b0, b1
          do a = v_lo, nfuc
            do c = v_lo, nfuc
              Gm(c,:,a,b-b0+1) = real(T(a,c,b,:), dp)
            end do
          end do
        end do
        !$omp end parallel do
        deallocate(D, T)
      end subroutine chol_vvvv_gm_r

end module chol_store
