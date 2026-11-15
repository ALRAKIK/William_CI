!---------------------------------------------------------------------!
! [Nov 2026] FAST version of build_ccsd_residuals (get_R.f90)         !
! used only with the keyword  fast_ccsdt  (the default is get_R.f90)  !
!                                                                     !
! Same equations as get_R.f90 (same result up to round-off), but      !
! organised in k blocks without per-element momentum tests:          !
!                                                                     !
!  R1 : explicit loops over the momentum-conserving k blocks only     !
!                                                                     !
!  R2(i,j,a,b;K) = conj g(i,a,j,b) + P(i,j,a,b;K) + P(j,i,b,a;K')     !
!                  - (e_i+e_j-e_a-e_b) t2(i,j,a,b)                    !
!                  + sum_{m,l} Woooo(m,l,i,j) tau(m,l,a,b)  (zgemm)   !
!                  + particle-particle ladder        (ladder_pp_gemm) !
!                  + ring terms W*t2                 (ring_gemm)      !
!  with the "pair" intermediate (one k slab at a time, zgemm):        !
!   P(i,j,a,b) =  sum_c Lvv(a,c) t2(i,j,c,b)  - e_a t2(i,j,a,b)       !
!               - sum_m Loo(m,i) t2(m,j,a,b)  + e_i t2(i,j,a,b)       !
!               + sum_{c in kj} z1(i,a,c,b) t1(j,c)                   !
!               - sum_{m in kb} z2(i,a,j,m) t1(m,b)                   !
!   z1(i,a,c,b) = conj g(i,a,c,b) - sum_{m in ka} g(m,i,b,c) t1(m,a)  !
!   z2(i,a,j,m) = conj g(i,a,j,m) + sum_{c in kj} g(m,c,a,i) t1(j,c)  !
!                                                                     !
! k <-> -k symmetry (use_ksym): P and R2 are computed only for the    !
! canonical k slabs, the others are complex conjugates.               !
! g is the COMPACT k-space integral array (module cc_dims, gq/gv1).   !
!---------------------------------------------------------------------!
subroutine build_ccsd_residuals_fast(blk,g,ktab,t1,t2,Foo,Fvv,Fov,Loo,Lvv, &
                                     Woooo,Wvoov,Wvovo,R1,R2)
      use constants_module, only: dp, dpc
      use unitcell_module , only: nfuc, N_cell
      use cc_dims         , only: o_hi, v_lo, use_ksym, gidx, nslot
      use k_blocks        , only: K_block
      use omp_lib                                   ! [Nov 2026]
      implicit none
      ! input !
      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(nfuc,nfuc,nfuc,nfuc,nslot)
      complex(dpc) , intent(in) :: t1(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: t2(o_hi,o_hi,v_lo:nfuc,v_lo:nfuc,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: Foo(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: Fvv(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: Fov(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: Loo(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: Lvv(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: Woooo(o_hi,o_hi,o_hi,o_hi,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: Wvoov(v_lo:nfuc,o_hi,o_hi,v_lo:nfuc,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: Wvovo(v_lo:nfuc,o_hi,v_lo:nfuc,o_hi,N_cell,N_cell,N_cell)
      integer      , intent(in) :: ktab(N_cell,N_cell,N_cell)
      ! output !
      complex(dpc), intent(out) :: R1(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: R2(o_hi,o_hi,v_lo:nfuc,v_lo:nfuc,N_cell,N_cell,N_cell)
      ! local !
      integer      :: no, nv, k, p, ki, kj, ka, kc, kd, km, kl, i, a, c, d, m, l
      integer      :: ii, aa, ip, nR1, ns, is
      logical      :: isO(nfuc,N_cell), isV(nfuc,N_cell)
      integer      :: kOcc(nfuc,N_cell), nOcc(N_cell), kVir(nfuc,N_cell), nVir(N_cell)
      integer      :: r1i(nfuc*nfuc*N_cell), r1a(nfuc*nfuc*N_cell), r1k(nfuc*nfuc*N_cell)
      integer      :: slab(3,N_cell**3), rslab(3,N_cell**3), nrs   ! [Nov 2026] rslab: pair-orbit representatives
      real(dp)     :: eps(nfuc,N_cell)
      complex(dpc) :: v
      complex(dpc), allocatable :: Pr(:,:,:,:,:,:,:)   ! pair intermediate
      complex(dpc), allocatable :: R1p(:,:,:,:)
      complex(dpc), parameter :: czero=(0.0_dp,0.0_dp), one=(1.0_dp,0.0_dp), mone=(-1.0_dp,0.0_dp)

      no = o_hi ; nv = nfuc - v_lo + 1
      nOcc=0; nVir=0; isO=.false.; isV=.false.
      do k=1,N_cell
        do p=1,nfuc
          eps(p,k) = blk(k)%energy(p)
          if (blk(k)%occupation(p)==2) then
            nOcc(k)=nOcc(k)+1; kOcc(nOcc(k),k)=p; isO(p,k)=.true.
          else if (blk(k)%occupation(p)==0) then
            nVir(k)=nVir(k)+1; kVir(nVir(k),k)=p; isV(p,k)=.true.
          else
            error stop 'CC residual: occupation must be 0 or 2'
          end if
        end do
      end do

      !=================================================================!
      ! R1 : one (k,i,a) per task, momentum-conserving blocks only       !
      !=================================================================!
      nR1=0
      do k=1,N_cell
        do ii=1,nOcc(k)
          do aa=1,nVir(k)
            nR1=nR1+1; r1i(nR1)=kOcc(ii,k); r1a(nR1)=kVir(aa,k); r1k(nR1)=k
          end do
        end do
      end do
      R1 = czero
      !$omp parallel do schedule(dynamic) default(shared) &
      !$omp private(ip,i,a,ki,km,kc,kd,kl,m,c,d,l,v)
      do ip=1,nR1
        i=r1i(ip); a=r1a(ip); ki=r1k(ip)
        v=czero
        ! Fvv t1 - Foo t1   (k_c = k_m = k_i)
        do c=1,nfuc
          if (isV(c,ki)) v=v+Fvv(a,c,ki)*t1(i,c,ki)
          if (isO(c,ki)) v=v-Foo(c,i,ki)*t1(c,a,ki)
        end do
        do km=1,N_cell
          ! Fov (2 t2(m,i,c,a) - t2(i,m,c,a) + t1 t1)  and  (2(mc|ai)-(mi|ac)) t1(m,c) : k_c = k_m
          do m=1,nfuc
            if (.not.isO(m,km)) cycle
            do c=1,nfuc
              if (.not.isV(c,km)) cycle
              v=v+Fov(m,c,km)*(2.0_dp*t2(m,i,c,a,km,ki,km)-t2(i,m,c,a,ki,km,km))
              if (km==ki) v=v+Fov(m,c,km)*t1(i,c,ki)*t1(m,a,ki)
              v=v+(2.0_dp*gq(m,c,a,i,km,km,ki)-gq(m,i,a,c,km,ki,ki))*t1(m,c,km)
            end do
          end do
        end do
        R1(i,a,ki)=v
      end do
      !$omp end parallel do

      ! ---- the two g*t2 terms of R1 with zgemm, one (ki,km) pair per task ----
      !   E: sum_{kc} sum_{m,c,d} (2(md|ac)-(mc|ad)) (t2(i,m,c,d) + t1(m,d)t1(i,c))
      !   F: sum_{kl} sum_{m,l,c} (-2(lc|mi)+(mc|li)) (t2(m,l,a,c) + t1(l,c)t1(m,a))
      allocate(R1p(no,v_lo:nfuc,N_cell,N_cell))
      ! [Nov 2026] few k blocks: serial loop, threaded BLAS inside
      !$omp parallel do if(N_cell*N_cell >= omp_get_max_threads()) schedule(dynamic) default(shared) private(ip,ki,km)
      do ip=1,N_cell*N_cell
        ki=(ip-1)/N_cell+1; km=mod(ip-1,N_cell)+1
        call r1_gemm(ki,km)
      end do
      !$omp end parallel do
      do ki=1,N_cell
        do km=1,N_cell
          do a=v_lo,nfuc; do i=1,no
            if (isO(i,ki).and.isV(a,ki)) R1(i,a,ki)=R1(i,a,ki)+R1p(i,a,ki,km)
          end do; end do
        end do
      end do
      deallocate(R1p)

      !=================================================================!
      ! R2                                                                !
      !=================================================================!
      ns=0
      do ka=1,N_cell; do kj=1,N_cell; do ki=1,N_cell
        if (use_ksym) then
          if (.not.canon3(ki,kj,ka)) cycle
        end if
        ns=ns+1; slab(:,ns)=[ki,kj,ka]
      end do; end do; end do

      ! ---- pair intermediate P on the canonical slabs ----
      allocate(Pr(o_hi,o_hi,v_lo:nfuc,v_lo:nfuc,N_cell,N_cell,N_cell))
      ! [Nov 2026] few k blocks: serial loop, threaded BLAS inside
      !$omp parallel do if(ns >= omp_get_max_threads()) schedule(dynamic) default(shared) private(is)
      do is=1,ns
        call pair_slab(slab(1,is),slab(2,is),slab(3,is))
      end do
      !$omp end parallel do
      if (use_ksym) call fill_conj(Pr)

      ! ---- R2 on the orbit representatives ----
      ! [Nov 2026] pair symmetry R2(i,j,a,b;ki,kj,ka) = R2(j,i,b,a;kj,ki,kb) (and
      ! k -> -k with use_ksym): only one slab of every orbit is computed, r2_slab
      ! stores it into all members of the orbit
      nrs=0
      do ka=1,N_cell; do kj=1,N_cell; do ki=1,N_cell
        if (.not.rep3(ki,kj,ka)) cycle
        nrs=nrs+1; rslab(:,nrs)=[ki,kj,ka]
      end do; end do; end do
      R2 = czero
      !$omp parallel do if(nrs >= omp_get_max_threads()) schedule(dynamic) default(shared) private(is)
      do is=1,nrs
        call r2_slab(rslab(1,is),rslab(2,is),rslab(3,is))
      end do
      !$omp end parallel do
      deallocate(Pr)

      call ladder_pp_gemm(blk,g,ktab,t1,t2,R2)
      call ring_gemm(blk,ktab,t2,Wvoov,Wvovo,R2)

contains

      !-----------------------------------------------------------------!
      ! R1p(:,:,ki,km): terms E and F of R1 for one (ki,km) pair         !
      !-----------------------------------------------------------------!
      subroutine r1_gemm(ki_,km_)
        integer, intent(in) :: ki_, km_
        integer :: kc_, kd_, kl_, i_, m_, l_, c_, d_, a_
        complex(dpc), allocatable :: V(:,:,:,:), TE(:,:,:,:), U(:,:,:,:), TF(:,:,:,:)
        allocate(V(no,v_lo:nfuc,v_lo:nfuc,v_lo:nfuc), TE(no,no,v_lo:nfuc,v_lo:nfuc))
        allocate(U(no,no,v_lo:nfuc,no), TF(no,no,v_lo:nfuc,v_lo:nfuc))
        R1p(:,:,ki_,km_) = czero
        ! ---- E ----
        do kc_=1,N_cell
          kd_ = ktab(ki_,km_,kc_)
          ! V(m,c,d,a) = 2 g(m,d,a,c; km,kd,ki) - g(m,c,a,d; km,kc,ki)
          do a_=v_lo,nfuc; do d_=v_lo,nfuc; do c_=v_lo,nfuc
            if (isV(c_,kc_).and.isV(d_,kd_)) then
              V(:,c_,d_,a_) = 2.0_dp*gv1(1,no,d_,a_,c_,km_,kd_,ki_) - gv1(1,no,c_,a_,d_,km_,kc_,ki_)
              do m_=1,no
                if (.not.isO(m_,km_)) V(m_,c_,d_,a_) = czero
              end do
            else
              V(:,c_,d_,a_) = czero
            end if
          end do; end do; end do
          ! TE(i,m,c,d) = t2(i,m,c,d; ki,km,kc) + t1(m,d) t1(i,c) [kc==ki]
          TE = t2(:,:,:,:,ki_,km_,kc_)
          if (kc_==ki_) then
            do d_=v_lo,nfuc; do c_=v_lo,nfuc; do m_=1,no; do i_=1,no
              TE(i_,m_,c_,d_) = TE(i_,m_,c_,d_) + t1(m_,d_,km_)*t1(i_,c_,ki_)
            end do; end do; end do; end do
          end if
          call zgemm('N','N', no, nv, no*nv*nv, one, TE, no, V, no*nv*nv, one, R1p(1,v_lo,ki_,km_), no)
        end do
        ! ---- F ----
        do kl_=1,N_cell
          kc_ = ktab(km_,kl_,ki_)
          ! U(m,l,c,i) = -2 g(l,c,m,i; kl,kc,km) + g(m,c,l,i; km,kc,kl)
          do i_=1,no; do c_=v_lo,nfuc; do m_=1,no
            U(m_,:,c_,i_) = -2.0_dp*gv1(1,no,c_,m_,i_,kl_,kc_,km_)
          end do; end do; end do
          do i_=1,no; do c_=v_lo,nfuc; do l_=1,no
            U(:,l_,c_,i_) = U(:,l_,c_,i_) + gv1(1,no,c_,l_,i_,km_,kc_,kl_)
          end do; end do; end do
          do i_=1,no; do c_=v_lo,nfuc; do l_=1,no; do m_=1,no
            if (.not.(isO(m_,km_).and.isO(l_,kl_).and.isV(c_,kc_))) U(m_,l_,c_,i_) = czero
          end do; end do; end do; end do
          ! TF(m,l,c,a) = t2(m,l,a,c; km,kl,ki) + t1(l,c) t1(m,a) [kc==kl, km==ki]
          do a_=v_lo,nfuc; do c_=v_lo,nfuc
            TF(:,:,c_,a_) = t2(:,:,a_,c_,km_,kl_,ki_)
            if (kc_==kl_ .and. km_==ki_) then
              do l_=1,no; do m_=1,no
                TF(m_,l_,c_,a_) = TF(m_,l_,c_,a_) + t1(l_,c_,kl_)*t1(m_,a_,km_)
              end do; end do
            end if
          end do; end do
          call zgemm('T','N', no, nv, no*no*nv, one, U, no*no*nv, TF, no*no*nv, one, R1p(1,v_lo,ki_,km_), no)
        end do
        deallocate(V, TE, U, TF)
      end subroutine r1_gemm

      !-----------------------------------------------------------------!
      ! P(:,:,:,:,ki,kj,ka) (see header), dense band ranges, masked     !
      !-----------------------------------------------------------------!
      subroutine pair_slab(ki_,kj_,ka_)
        integer, intent(in) :: ki_, kj_, ka_
        integer :: kb_, i_, j_, a_, b_, c_, m_
        complex(dpc), allocatable :: Lv(:,:), Lo(:,:), t1j(:,:), t1a(:,:), t1b(:,:)
        complex(dpc), allocatable :: Pb(:,:,:,:), Z1(:,:,:,:), G(:,:,:,:), X(:,:,:,:), Q(:,:,:)
        complex(dpc), allocatable :: Y(:,:,:,:), H(:,:,:,:), Y2(:,:,:,:), Q2(:,:,:,:)

        kb_ = ktab(ki_,kj_,ka_)
        allocate(Lv(v_lo:nfuc,v_lo:nfuc), Lo(no,no), t1j(no,v_lo:nfuc), t1a(no,v_lo:nfuc), t1b(no,v_lo:nfuc))
        allocate(Pb(no,no,v_lo:nfuc,v_lo:nfuc))
        ! masked one-body pieces
        do c_=v_lo,nfuc; do a_=v_lo,nfuc
          Lv(a_,c_)=czero; if (isV(a_,ka_).and.isV(c_,ka_)) Lv(a_,c_)=Lvv(a_,c_,ka_)
        end do; end do
        do i_=1,no; do m_=1,no
          Lo(m_,i_)=czero; if (isO(m_,ki_).and.isO(i_,ki_)) Lo(m_,i_)=Loo(m_,i_,ki_)
        end do; end do
        do c_=v_lo,nfuc; do j_=1,no
          t1j(j_,c_)=czero; if (isO(j_,kj_).and.isV(c_,kj_)) t1j(j_,c_)=t1(j_,c_,kj_)
          t1a(j_,c_)=czero; if (isO(j_,ka_).and.isV(c_,ka_)) t1a(j_,c_)=t1(j_,c_,ka_)
          t1b(j_,c_)=czero; if (isO(j_,kb_).and.isV(c_,kb_)) t1b(j_,c_)=t1(j_,c_,kb_)
        end do; end do

        ! ---- Lvv t2 - e_a t2 ;  - Loo t2 + e_i t2 ----
        Pb = czero
        do b_=v_lo,nfuc
          call zgemm('N','T', no*no, nv, nv, one, t2(1,1,v_lo,b_,ki_,kj_,ka_), no*no, Lv, nv, &
                     one, Pb(1,1,v_lo,b_), no*no)
        end do
        call zgemm('T','N', no, no*nv*nv, no, mone, Lo, no, t2(1,1,v_lo,v_lo,ki_,kj_,ka_), no, one, Pb, no)
        do b_=v_lo,nfuc; do a_=v_lo,nfuc; do j_=1,no; do i_=1,no
          Pb(i_,j_,a_,b_)=Pb(i_,j_,a_,b_)+(eps(i_,ki_)-eps(a_,ka_))*t2(i_,j_,a_,b_,ki_,kj_,ka_)
        end do; end do; end do; end do

        ! ---- + sum_{c in kj} z1(i,a,c,b) t1(j,c) ----
        allocate(Z1(no,v_lo:nfuc,v_lo:nfuc,v_lo:nfuc), G(no,no,v_lo:nfuc,v_lo:nfuc), &
                 X(no,v_lo:nfuc,v_lo:nfuc,v_lo:nfuc), Q(no,v_lo:nfuc,no))
        do b_=v_lo,nfuc; do c_=v_lo,nfuc; do a_=v_lo,nfuc
          Z1(:,a_,c_,b_) = conjg(gv1(1,no,a_,c_,b_,ki_,ka_,kj_))
        end do; end do; end do
        ! G(m,i,b,c) = g(m,i,b,c; ka,ki,kb)  ->  X(i,b,c,a) = sum_m G(m,i,b,c) t1(m,a)
        do c_=v_lo,nfuc; do b_=v_lo,nfuc; do i_=1,no
          G(:,i_,b_,c_) = gv1(1,no,i_,b_,c_,ka_,ki_,kb_)
        end do; end do; end do
        call zgemm('T','N', no*nv*nv, nv, no, one, G, no, t1a, no, czero, X, no*nv*nv)
        do b_=v_lo,nfuc; do c_=v_lo,nfuc; do a_=v_lo,nfuc
          Z1(:,a_,c_,b_) = Z1(:,a_,c_,b_) - X(:,b_,c_,a_)
        end do; end do; end do
        do b_=v_lo,nfuc
          call zgemm('N','T', no*nv, no, nv, one, Z1(1,v_lo,v_lo,b_), no*nv, t1j, no, czero, Q, no*nv)
          do a_=v_lo,nfuc; do j_=1,no
            Pb(:,j_,a_,b_) = Pb(:,j_,a_,b_) + Q(:,a_,j_)
          end do; end do
        end do
        deallocate(Z1, G, X, Q)

        ! ---- - sum_{m in kb} z2(i,a,j,m) t1(m,b) ----
        allocate(Y(no,v_lo:nfuc,no,no), H(v_lo:nfuc,no,v_lo:nfuc,no), Y2(no,v_lo:nfuc,no,no), Q2(no,v_lo:nfuc,no,v_lo:nfuc))
        do m_=1,no; do j_=1,no; do a_=v_lo,nfuc
          Y(:,a_,j_,m_) = conjg(gv1(1,no,a_,j_,m_,ki_,ka_,kj_))
        end do; end do; end do
        ! H(c,m,a,i) = g(m,c,a,i; kb,kj,ka)  ->  Y2(m,a,i,j) = sum_c H(c,m,a,i) t1(j,c)
        do i_=1,no; do a_=v_lo,nfuc; do c_=v_lo,nfuc
          H(c_,:,a_,i_) = gv1(1,no,c_,a_,i_,kb_,kj_,ka_)
        end do; end do; end do
        call zgemm('T','T', no*nv*no, no, nv, one, H, nv, t1j, no, czero, Y2, no*nv*no)
        do m_=1,no; do j_=1,no; do a_=v_lo,nfuc; do i_=1,no
          Y(i_,a_,j_,m_) = Y(i_,a_,j_,m_) + Y2(m_,a_,i_,j_)
        end do; end do; end do; end do
        call zgemm('N','N', no*nv*no, nv, no, mone, Y, no*nv*no, t1b, no, czero, Q2, no*nv*no)
        do b_=v_lo,nfuc; do a_=v_lo,nfuc; do j_=1,no
          Pb(:,j_,a_,b_) = Pb(:,j_,a_,b_) + Q2(:,a_,j_,b_)
        end do; end do; end do
        deallocate(Y, H, Y2, Q2)

        ! ---- store (physical entries only) ----
        do b_=v_lo,nfuc; do a_=v_lo,nfuc; do j_=1,no; do i_=1,no
          if (isO(i_,ki_).and.isO(j_,kj_).and.isV(a_,ka_).and.isV(b_,kb_)) then
            Pr(i_,j_,a_,b_,ki_,kj_,ka_) = Pb(i_,j_,a_,b_)
          else
            Pr(i_,j_,a_,b_,ki_,kj_,ka_) = czero
          end if
        end do; end do; end do; end do
        deallocate(Lv, Lo, t1j, t1a, t1b, Pb)
      end subroutine pair_slab

      !-----------------------------------------------------------------!
      ! R2(:,:,:,:,ki,kj,ka) without the pp ladder and ring terms        !
      !-----------------------------------------------------------------!
      subroutine r2_slab(ki_,kj_,ka_)
        integer, intent(in) :: ki_, kj_, ka_
        integer :: kb_, km_, kl_, i_, j_, a_, b_, m_, l_
        integer :: mem(3,4), nm_, e_, q_, k1_, k2_, k3_
        logical :: msw(4), mcj(4)
        complex(dpc), allocatable :: R(:,:,:,:), TAU(:,:,:,:)

        kb_ = ktab(ki_,kj_,ka_)
        allocate(R(no,no,v_lo:nfuc,v_lo:nfuc), TAU(no,no,v_lo:nfuc,v_lo:nfuc))
        do b_=v_lo,nfuc; do a_=v_lo,nfuc; do j_=1,no
          R(:,j_,a_,b_) = conjg(gv1(1,no,a_,j_,b_,ki_,ka_,kj_))
        end do; end do; end do
        do b_=v_lo,nfuc; do a_=v_lo,nfuc; do j_=1,no; do i_=1,no
          R(i_,j_,a_,b_) = R(i_,j_,a_,b_) + Pr(i_,j_,a_,b_,ki_,kj_,ka_) + Pr(j_,i_,b_,a_,kj_,ki_,kb_) &
                         - (eps(i_,ki_)+eps(j_,kj_)-eps(a_,ka_)-eps(b_,kb_))*t2(i_,j_,a_,b_,ki_,kj_,ka_)
        end do; end do; end do; end do
        ! hole-hole ladder : R(ij,ab) += sum_{km} Woooo(ml,ij)^T tau(ml,ab) , kl = ki+kj-km
        do km_=1,N_cell
          kl_ = ktab(ki_,kj_,km_)
          TAU = t2(:,:,:,:,km_,kl_,ka_)
          if (km_==ka_) then
            do b_=v_lo,nfuc; do a_=v_lo,nfuc
              if (.not.(isV(a_,ka_).and.isV(b_,kb_))) cycle
              do l_=1,no
                if (.not.isO(l_,kl_)) cycle
                do m_=1,no
                  if (isO(m_,km_)) TAU(m_,l_,a_,b_) = TAU(m_,l_,a_,b_) + t1(m_,a_,km_)*t1(l_,b_,kl_)
                end do
              end do
            end do; end do
          end if
          call zgemm('T','N', no*no, nv*nv, no*no, one, Woooo(1,1,1,1,km_,kl_,ki_), no*no, &
                     TAU, no*no, one, R, no*no)
        end do
        ! store physical entries into every distinct member of the orbit
        mem(:,1)=[ki_,kj_,ka_]; msw(1)=.false.; mcj(1)=.false.
        mem(:,2)=[kj_,ki_,kb_]; msw(2)=.true. ; mcj(2)=.false.
        nm_=2
        if (use_ksym) then
          nm_=4
          mem(:,3)=[ktab(1,1,ki_),ktab(1,1,kj_),ktab(1,1,ka_)]; msw(3)=.false.; mcj(3)=.true.
          mem(:,4)=[ktab(1,1,kj_),ktab(1,1,ki_),ktab(1,1,kb_)]; msw(4)=.true. ; mcj(4)=.true.
        end if
        do e_=1,nm_
          if (any([(all(mem(:,q_)==mem(:,e_)), q_=1,e_-1)])) cycle
          k1_=mem(1,e_); k2_=mem(2,e_); k3_=mem(3,e_)
          if (.not.msw(e_)) then
            ! same index order: contiguous copy (conjugated for the -k member)
            do b_=v_lo,nfuc; if (.not.isV(b_,kb_)) cycle
            do a_=v_lo,nfuc; if (.not.isV(a_,ka_)) cycle
            do j_=1,no;      if (.not.isO(j_,kj_)) cycle
              if (mcj(e_)) then
                do i_=1,no
                  if (isO(i_,ki_)) R2(i_,j_,a_,b_,k1_,k2_,k3_) = conjg(R(i_,j_,a_,b_))
                end do
              else
                do i_=1,no
                  if (isO(i_,ki_)) R2(i_,j_,a_,b_,k1_,k2_,k3_) = R(i_,j_,a_,b_)
                end do
              end if
            end do; end do; end do
          else
            ! swapped member R2(j,i,b,a): write contiguously along j
            do a_=v_lo,nfuc; if (.not.isV(a_,ka_)) cycle
            do b_=v_lo,nfuc; if (.not.isV(b_,kb_)) cycle
            do i_=1,no;      if (.not.isO(i_,ki_)) cycle
              if (mcj(e_)) then
                do j_=1,no
                  if (isO(j_,kj_)) R2(j_,i_,b_,a_,k1_,k2_,k3_) = conjg(R(i_,j_,a_,b_))
                end do
              else
                do j_=1,no
                  if (isO(j_,kj_)) R2(j_,i_,b_,a_,k1_,k2_,k3_) = R(i_,j_,a_,b_)
                end do
              end if
            end do; end do; end do
          end if
        end do
        deallocate(R, TAU)
      end subroutine r2_slab

      ! X(-ki,-kj,-ka) = conj X(ki,kj,ka) for the non-canonical slabs
      subroutine fill_conj(Xa)
        complex(dpc), intent(inout) :: Xa(o_hi,o_hi,v_lo:nfuc,v_lo:nfuc,N_cell,N_cell,N_cell)
        integer :: k1, k2, k3
        !$omp parallel do collapse(3) default(shared) private(k1,k2,k3)
        do k1=1,N_cell
          do k2=1,N_cell
            do k3=1,N_cell
              if (.not.canon3(k1,k2,k3)) Xa(:,:,:,:,k1,k2,k3) = &
                 conjg(Xa(:,:,:,:,ktab(1,1,k1),ktab(1,1,k2),ktab(1,1,k3)))
            end do
          end do
        end do
        !$omp end parallel do
      end subroutine fill_conj

      ! (k1,k2,k3) <= (-k1,-k2,-k3) in lexicographic order
      ! [Nov 2026] (k1,k2,k3) is the lexicographically smallest member of its
      ! orbit {(k1,k2,k3), (k2,k1,k4)} (plus the negatives with use_ksym)
      logical function rep3(k1,k2,k3)
        integer, intent(in) :: k1,k2,k3
        integer :: m(3,4), nm, r
        m(:,1)=[k1,k2,k3]
        m(:,2)=[k2,k1,ktab(k1,k2,k3)]
        nm=2
        if (use_ksym) then
          nm=4
          m(:,3)=[ktab(1,1,k1),ktab(1,1,k2),ktab(1,1,k3)]
          m(:,4)=[ktab(1,1,m(1,2)),ktab(1,1,m(2,2)),ktab(1,1,m(3,2))]
        end if
        rep3=.true.
        do r=2,nm
          if (lexlt(m(:,r),m(:,1))) rep3=.false.
        end do
      end function rep3
      logical function lexlt(x,y)
        integer, intent(in) :: x(3), y(3)
        lexlt = (x(1)<y(1)) .or. (x(1)==y(1) .and. (x(2)<y(2) .or. (x(2)==y(2) .and. x(3)<y(3))))
      end function lexlt

      logical function canon3(k1,k2,k3)
        integer, intent(in) :: k1,k2,k3
        integer :: m1,m2,m3
        m1=ktab(1,1,k1); m2=ktab(1,1,k2); m3=ktab(1,1,k3)
        canon3 = (k1<m1) .or. (k1==m1 .and. (k2<m2 .or. (k2==m2 .and. k3<=m3)))
      end function canon3

      ! compact g: block (k1,k2,k3) = g(:,:,:,:,|gidx|), conjugated if gidx < 0
      complex(dpc) function gq(p_,q_,r_,s_,k1_,k2_,k3_)
        integer, intent(in) :: p_,q_,r_,s_,k1_,k2_,k3_
        integer :: ic_
        ic_ = gidx(k1_,k2_,k3_)
        if (ic_ > 0) then
          gq = g(p_,q_,r_,s_,ic_)
        else
          gq = conjg(g(p_,q_,r_,s_,-ic_))
        end if
      end function gq
      ! vector g(lo:hi,q,r,s,k1,k2,k3) (first index range)
      function gv1(lo_,hi_,q_,r_,s_,k1_,k2_,k3_) result(v_)
        integer, intent(in) :: lo_,hi_,q_,r_,s_,k1_,k2_,k3_
        complex(dpc) :: v_(hi_-lo_+1)
        integer :: ic_
        ic_ = gidx(k1_,k2_,k3_)
        if (ic_ > 0) then
          v_ = g(lo_:hi_,q_,r_,s_,ic_)
        else
          v_ = conjg(g(lo_:hi_,q_,r_,s_,-ic_))
        end if
      end function gv1

end subroutine build_ccsd_residuals_fast
