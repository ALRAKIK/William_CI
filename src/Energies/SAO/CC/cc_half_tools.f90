!---------------------------------------------------------------------!
! [Nov 2026] fast_ccsdt: copies of get_t_term, build_ccsd_f,          !
! build_ccsd_l (get_t.f90, get_F.f90, get_L.f90) and update_vvvv_in_g  !
! (get_W.f90) for the COMPACT k-space integrals                        !
!   g(p,q,r,s,islot), block (k1,k2,k3) -> slot gidx(k1,k2,k3)          !
!   (conjugated when gidx < 0, k <-> -k symmetry; module cc_dims).     !
! Same formulas; every g(...,k1,k2,k3) read is gq(...,k1,k2,k3).       !
!---------------------------------------------------------------------!

subroutine get_t_term_h(blk,g,ktab,t2)

      use constants_module, only: dp, dpc
      use cc_dims         , only: gidx, nslot, g_n1      ! [Nov 2026] compact g
      use unitcell_module , only: nfuc, N_cell
      use cc_dims       , only: o_hi, v_lo
      use k_blocks        , only: K_block

      implicit none 

      ! input ! 

      type(K_block), intent(in)  :: blk(N_cell)
      complex(dpc) , intent(in)  :: g(g_n1,nfuc,nfuc,nfuc,nslot)   ! [Nov 2026] compact g (cc_dims gidx)
      integer      , intent(in)  :: ktab(N_cell,N_cell,N_cell)

      ! local ! 

      integer                    :: ki,kj,ka,kb
      integer                    :: i,j,a,b
      integer                    :: k,nocc
      real(dp)                   :: gap
      complex(dpc)               :: gd

      ! output ! 

      complex(dpc), intent(out)  :: t2(o_hi,o_hi,v_lo:nfuc,v_lo:nfuc,N_cell,N_cell,N_cell)

      nocc = 0
      do k=1,N_cell
        nocc = nocc + count(blk(k)%occupation == 2)
      end do

      t2 = cmplx(0.0_dp,0.0_dp,kind=dpc)

      do ka=1,N_cell
        do kj=1,N_cell
          do ki=1,N_cell
            kb = ktab(ki,kj,ka)
            do b=1,nfuc
              if (blk(kb)%occupation(b) /= 0) cycle
              do a=1,nfuc
                if (blk(ka)%occupation(a) /= 0) cycle
                do j=1,nfuc
                  if (blk(kj)%occupation(j) /= 2) cycle
                  do i=1,nfuc
                    if (blk(ki)%occupation(i) /= 2) cycle
                    gap = blk(ka)%energy(a)+blk(kb)%energy(b) &
                        - blk(ki)%energy(i)-blk(kj)%energy(j)
                    if (gap <= 1.0e-12_dp) error stop 'CC baseline: nonpositive excitation gap'
                    gd = gq(i,a,j,b,ki,ka,kj)
                    t2(i,j,a,b,ki,kj,ka) = -conjg(gd)/gap
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      contains
      ! [Nov 2026] compact g: block (k1,k2,k3) = g(:,:,:,:,|gidx|), conjugated if gidx < 0
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
      
end subroutine get_t_term_h

subroutine build_ccsd_f_h(blk,g,ktab,t1,t2,Foo,Fvv,Fov)
      ! [Nov 2026] F intermediates with zgemm (same formulas as get_F.f90)
      !   Foo(m,i) = f_mi + sum_{n,e,f} V(m;n,e,f) tau(i;n,e,f)
      !   Fvv(a,e) = f_ae - sum_{m,n,f} tau(m,n,a,f) V(m,e,n,f)
      !   Fov(m,e) = sum_{n,f} V(m,e,n,f) t1(n,f)
      !   V(m,e,n,f) = 2(me|nf) - (mf|ne) ,  tau = t2 + t1 t1
      ! partial sums per k pair (parallel), then summed in k order.
      use constants_module, only: dp, dpc
      use unitcell_module , only: nfuc, N_cell
      use cc_dims         , only: o_hi, v_lo, gidx, nslot, g_n1
      use k_blocks        , only: K_block
      use omp_lib                                   ! [Nov 2026]
      implicit none
      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(g_n1,nfuc,nfuc,nfuc,nslot)
      complex(dpc) , intent(in) :: t1(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: t2(o_hi,o_hi,v_lo:nfuc,v_lo:nfuc,N_cell,N_cell,N_cell)
      integer      , intent(in) :: ktab(N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: Foo(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: Fvv(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: Fov(nfuc,nfuc,N_cell)
      integer      :: no, nv, k, p, k2, ip, kn, ke, kf, km, m, n, e, f, i, a
      logical      :: isO(nfuc,N_cell), isV(nfuc,N_cell)
      complex(dpc), allocatable :: Fo(:,:,:,:), Fv(:,:,:,:), Vb(:,:,:,:), Tb(:,:,:,:)
      complex(dpc) :: v
      complex(dpc), parameter :: czero=(0.0_dp,0.0_dp), one=(1.0_dp,0.0_dp)

      no = o_hi ; nv = nfuc - v_lo + 1
      do k=1,N_cell
        do p=1,nfuc
          isO(p,k) = (blk(k)%occupation(p)==2)
          isV(p,k) = (blk(k)%occupation(p)==0)
        end do
      end do
      allocate(Fo(no,no,N_cell,N_cell), Fv(v_lo:nfuc,v_lo:nfuc,N_cell,N_cell))

      ! ---- Foo : pair (k, kn), sum over ke inside ----
      ! [Nov 2026] few k blocks: serial loop, threaded BLAS inside
      !$omp parallel do if(N_cell*N_cell >= omp_get_max_threads()) schedule(dynamic) default(shared) private(ip,k,kn,ke,kf,n,e,f,i,Vb,Tb)
      do ip=1,N_cell*N_cell
        k  = (ip-1)/N_cell + 1
        kn = mod(ip-1,N_cell) + 1
        allocate(Vb(no,no,v_lo:nfuc,v_lo:nfuc), Tb(no,no,v_lo:nfuc,v_lo:nfuc))
        Fo(:,:,k,kn) = czero
        do ke=1,N_cell
          kf = ktab(k,kn,ke)
          do f=v_lo,nfuc; do e=v_lo,nfuc; do n=1,no
            if (isO(n,kn).and.isV(e,ke).and.isV(f,kf)) then
              Vb(:,n,e,f) = 2.0_dp*gv1(1,no,e,n,f,k,ke,kn) - gv1(1,no,f,n,e,k,kf,kn)
            else
              Vb(:,n,e,f) = czero
            end if
          end do; end do; end do
          Tb = t2(:,:,:,:,k,kn,ke)
          if (ke==k .and. kf==kn) then
            do f=v_lo,nfuc; do e=v_lo,nfuc; do n=1,no; do i=1,no
              Tb(i,n,e,f) = Tb(i,n,e,f) + t1(i,e,k)*t1(n,f,kn)
            end do; end do; end do; end do
          end if
          call zgemm('N','T', no, no, no*nv*nv, one, Vb, no, Tb, no, one, Fo(1,1,k,kn), no)
        end do
        deallocate(Vb, Tb)
      end do
      !$omp end parallel do

      ! ---- Fvv : pair (k, km), sum over kn inside ----
      !$omp parallel do if(N_cell*N_cell >= omp_get_max_threads()) schedule(dynamic) default(shared) private(ip,k,km,kn,kf,m,n,e,f,a,Vb,Tb)
      do ip=1,N_cell*N_cell
        k  = (ip-1)/N_cell + 1
        km = mod(ip-1,N_cell) + 1
        allocate(Vb(no,no,v_lo:nfuc,v_lo:nfuc), Tb(no,no,v_lo:nfuc,v_lo:nfuc))
        Fv(:,:,k,km) = czero
        do kn=1,N_cell
          kf = ktab(km,kn,k)
          ! Vb(m,n,f,e) = 2(me|nf) - (mf|ne)  (k_m, k_e = k, k_n, k_f)
          do e=v_lo,nfuc; do f=v_lo,nfuc; do n=1,no
            if (isO(n,kn).and.isV(f,kf).and.isV(e,k)) then
              Vb(:,n,f,e) = 2.0_dp*gv1(1,no,e,n,f,km,k,kn) - gv1(1,no,f,n,e,km,kf,kn)
              do m=1,no
                if (.not.isO(m,km)) Vb(m,n,f,e) = czero
              end do
            else
              Vb(:,n,f,e) = czero
            end if
          end do; end do; end do
          ! Tb(m,n,f,a) = tau(m,n,a,f)
          do a=v_lo,nfuc; do f=v_lo,nfuc
            Tb(:,:,f,a) = t2(:,:,a,f,km,kn,k)
            if (km==k .and. kf==kn) then
              do n=1,no; do m=1,no
                Tb(m,n,f,a) = Tb(m,n,f,a) + t1(m,a,km)*t1(n,f,kn)
              end do; end do
            end if
          end do; end do
          call zgemm('T','N', nv, nv, no*no*nv, -one, Tb, no*no*nv, Vb, no*no*nv, one, Fv(v_lo,v_lo,k,km), nv)
        end do
        deallocate(Vb, Tb)
      end do
      !$omp end parallel do

      Foo = czero ; Fvv = czero ; Fov = czero
      do k=1,N_cell
        do i=1,nfuc
          if (isO(i,k)) Foo(i,i,k) = cmplx(blk(k)%energy(i),0.0_dp,kind=dpc)
          if (isV(i,k)) Fvv(i,i,k) = cmplx(blk(k)%energy(i),0.0_dp,kind=dpc)
        end do
        do k2=1,N_cell
          do i=1,no; do m=1,no
            if (isO(m,k).and.isO(i,k)) Foo(m,i,k) = Foo(m,i,k) + Fo(m,i,k,k2)
          end do; end do
          do e=v_lo,nfuc; do a=v_lo,nfuc
            if (isV(a,k).and.isV(e,k)) Fvv(a,e,k) = Fvv(a,e,k) + Fv(a,e,k,k2)
          end do; end do
        end do
      end do
      deallocate(Fo, Fv)

      ! ---- Fov(m,e,k) = sum_{kn,n,f} V(m,e,n,f) t1(n,f,kn) ----
      !$omp parallel do collapse(2) default(shared) private(k,m,e,kn,n,f,v)
      do k=1,N_cell
        do m=1,nfuc
          if (.not.isO(m,k)) cycle
          do e=1,nfuc
            if (.not.isV(e,k)) cycle
            v = czero
            do kn=1,N_cell
              do n=1,nfuc
                if (.not.isO(n,kn)) cycle
                do f=1,nfuc
                  if (.not.isV(f,kn)) cycle
                  v = v + (2.0_dp*gq(m,e,n,f,k,k,kn) - gq(m,f,n,e,k,kn,kn))*t1(n,f,kn)
                end do
              end do
            end do
            Fov(m,e,k) = v
          end do
        end do
      end do
      !$omp end parallel do

contains
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
end subroutine build_ccsd_f_h

subroutine build_ccsd_l_h(blk,g,t1,Foo,Fvv,Loo,Lvv)
      
      use constants_module, only: dp, dpc
      use cc_dims         , only: gidx, nslot, g_n1      ! [Nov 2026] compact g
      use unitcell_module , only: nfuc, N_cell
      use k_blocks        , only: K_block

      implicit none 

      ! input ! 

      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(g_n1,nfuc,nfuc,nfuc,nslot)   ! [Nov 2026] compact g (cc_dims gidx)
      complex(dpc) , intent(in) :: t1(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: Foo(nfuc,nfuc,N_cell),Fvv(nfuc,nfuc,N_cell)

      ! local ! 

      integer                   :: k,kl,kd,m,i,l,c,a,e,d
      integer                   :: kOcc(nfuc,N_cell), nOcc(N_cell)
      integer                   :: kVir(nfuc,N_cell), nVir(N_cell)
      integer                   :: ii , mm , ee, aa , ll , cc , mm2 , dd 
      integer                   :: p 

      ! output ! 

      complex(dpc), intent(out) :: Loo(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: Lvv(nfuc,nfuc,N_cell)
      
      

      ! Canonical HF assumption: bare f_ov=0, matching build_ccsd_f.
      ! Loo(m,i)=Foo(m,i)+sum_lc [2(lc|mi)-(mc|li)]*t1(l,c)
      ! Lvv(a,e)=Fvv(a,e)+sum_md [2(md|ae)-(me|ad)]*t1(m,d)
      ! The bare f_ov here must NOT be replaced by dressed Fov.

      nOcc=0; nVir=0
      do k=1,N_cell
        do p=1,nfuc
          if (blk(k)%occupation(p)==2) then
            nOcc(k)=nOcc(k)+1; kOcc(nOcc(k),k)=p
          else if (blk(k)%occupation(p)==0) then
            nVir(k)=nVir(k)+1; kVir(nVir(k),k)=p
          end if
        end do
      end do


      Loo=cmplx(0.0_dp,0.0_dp,kind=dpc)
      Lvv=cmplx(0.0_dp,0.0_dp,kind=dpc)

            ! ---- Loo: m,i occupied at k ; l occupied, c virtual at kl ----
      do k=1,N_cell
        do mm=1,nOcc(k)
          m=kOcc(mm,k)
          do ii=1,nOcc(k)
            i=kOcc(ii,k)
            Loo(m,i,k)=Foo(m,i,k)
            do kl=1,N_cell
              do ll=1,nOcc(kl)
                l=kOcc(ll,kl)
                do cc=1,nVir(kl)
                  c=kVir(cc,kl)
                  Loo(m,i,k)=Loo(m,i,k)+(2.0_dp*gq(l,c,m,i,kl,kl,k) &
                    -gq(m,c,l,i,k,kl,kl))*t1(l,c,kl)
                end do
              end do
            end do
          end do
        end do
      end do

      ! ---- Lvv: a,e virtual at k ; m occupied, d virtual at kd ----
      do k=1,N_cell
        do aa=1,nVir(k)
          a=kVir(aa,k)
          do ee=1,nVir(k)
            e=kVir(ee,k)
            Lvv(a,e,k)=Fvv(a,e,k)
            do kd=1,N_cell
              do mm2=1,nOcc(kd)
                m=kOcc(mm2,kd)
                do dd=1,nVir(kd)
                  d=kVir(dd,kd)
                  Lvv(a,e,k)=Lvv(a,e,k)+(2.0_dp*gq(m,d,a,e,kd,kd,k) &
                    -gq(m,e,a,d,kd,k,k))*t1(m,d,kd)
                end do
              end do
            end do
          end do
        end do
      end do
      contains
      ! [Nov 2026] compact g: block (k1,k2,k3) = g(:,:,:,:,|gidx|), conjugated if gidx < 0
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

end subroutine build_ccsd_l_h

subroutine update_vvvv_in_g_h(blk,g,ktab,dt1)

      ! [Nov 2026] zgemm version. For every stored block (kp,kq,kr) and every s:
      !   W(p,r,q) = sum_l g(l,s,p,r; kq,ks,kp) dt1(l,q;kq)      (gemm 1)
      !            + sum_l dt1(l,p;kp) g(l,r,q,s; kp,kr,kq)      (gemm 2)
      !   g(p,r,q,s; kp,kr,kq) -= W            (physical virtual entries only)
      ! Same result as the former scalar O(o v^4) loop, which ran every iteration.

      use constants_module, only: dp, dpc
      use cc_dims         , only: gidx, nslot, g_n1, o_hi, v_lo
      use unitcell_module , only: nfuc, N_cell
      use k_blocks        , only: K_block
      use chol_store      , only: chol_on          ! [Nov 2026] keyword cholesky
      use cc_disk         , only: disk_on          ! [Nov 2026] keyword CC_to_disk

      implicit none

      type(K_block), intent(in)    :: blk(N_cell)
      complex(dpc) , intent(inout) :: g(g_n1,nfuc,nfuc,nfuc,nslot)
      integer      , intent(in)    :: ktab(N_cell,N_cell,N_cell)
      complex(dpc) , intent(in)    :: dt1(nfuc,nfuc,N_cell)

      integer      :: k, p, q, r, kp, kq, kr, ks, ic, i1, i2, no, nv, it, ntk, bs
      logical      :: isV(nfuc,N_cell), isO(nfuc,N_cell), allphys
      integer, allocatable :: tk(:,:)
      complex(dpc), allocatable :: D(:,:,:), A1(:,:,:), A2(:,:,:), W(:,:,:)
      complex(dpc), parameter :: zero=(0.0_dp,0.0_dp), one=(1.0_dp,0.0_dp)

      ! [Nov 2026] keyword cholesky: Wvvvv is built from the Cholesky vectors
      ! in the ladder (chol_vvvv_gm); g holds no vvvv block
      ! (CC_to_disk: g holds no vvvv block either, Wvvvv is formed in the ladder)
      if (chol_on .or. disk_on) return
      no = o_hi ; nv = nfuc - v_lo + 1
      if (no <= 0 .or. nv <= 0) return
      do k=1,N_cell
        isO(:,k) = blk(k)%occupation(:)==2
        isV(:,k) = blk(k)%occupation(:)==0
      end do
      allphys = all(isV(v_lo:nfuc,:))
      ! masked dt1 block D(l,a,k): occupied l, virtual a only
      allocate(D(no,v_lo:nfuc,N_cell))
      do k=1,N_cell
        do q=v_lo,nfuc; do p=1,no
          D(p,q,k) = zero
          if (isO(p,k) .and. isV(q,k)) D(p,q,k) = dt1(p,q,k)
        end do; end do
      end do
      ! task list: (stored block, s)
      ntk = 0
      do kr=1,N_cell; do kq=1,N_cell; do kp=1,N_cell
        if (gidx(kp,kr,kq) > 0) ntk = ntk + nv
      end do; end do; end do
      allocate(tk(4,ntk)); ntk = 0
      do kr=1,N_cell; do kq=1,N_cell; do kp=1,N_cell
        if (gidx(kp,kr,kq) <= 0) cycle
        do bs=v_lo,nfuc
          ntk = ntk+1; tk(:,ntk) = [kp,kq,kr,bs]
        end do
      end do; end do; end do

      !$omp parallel default(shared) private(it,kp,kq,kr,ks,bs,ic,i1,i2,p,q,r,A1,A2,W)
      allocate(A1(no,v_lo:nfuc,v_lo:nfuc), A2(no,v_lo:nfuc,v_lo:nfuc), W(v_lo:nfuc,v_lo:nfuc,v_lo:nfuc))
      !$omp do schedule(dynamic,4)
      do it=1,ntk
        kp=tk(1,it); kq=tk(2,it); kr=tk(3,it); bs=tk(4,it)
        ks=ktab(kp,kq,kr)
        if (.not.isV(bs,ks)) cycle
        ic=gidx(kp,kr,kq)
        i1=gidx(kq,ks,kp)          ! A1(l,p,r) = g(l,s,p,r; kq,ks,kp)
        i2=gidx(kp,kr,kq)          ! A2(l,r,q) = g(l,r,q,s; kp,kr,kq)  (= stored block ic)
        if (i1 > 0) then
          A1 = g(1:no,bs,v_lo:nfuc,v_lo:nfuc,i1)
        else
          A1 = conjg(g(1:no,bs,v_lo:nfuc,v_lo:nfuc,-i1))
        end if
        A2 = g(1:no,v_lo:nfuc,v_lo:nfuc,bs,i2)
        ! W(p,r,q) = A1(l,[p,r])^T D(l,q;kq)
        call zgemm('T','N', nv*nv, nv, no, one, A1, no, D(1,v_lo,kq), no, zero, W, nv*nv)
        ! W(p,[r,q]) += D(l,p;kp)^T A2(l,[r,q])
        call zgemm('T','N', nv, nv*nv, no, one, D(1,v_lo,kp), no, A2, no, one, W, nv)
        if (allphys) then
          do q=v_lo,nfuc; do r=v_lo,nfuc
            g(v_lo:nfuc,r,q,bs,ic) = g(v_lo:nfuc,r,q,bs,ic) - W(:,r,q)
          end do; end do
        else
          do q=v_lo,nfuc; if (.not.isV(q,kq)) cycle
          do r=v_lo,nfuc; if (.not.isV(r,kr)) cycle
          do p=v_lo,nfuc; if (.not.isV(p,kp)) cycle
            g(p,r,q,bs,ic) = g(p,r,q,bs,ic) - W(p,r,q)
          end do; end do; end do
        end if
      end do
      !$omp end do
      deallocate(A1, A2, W)
      !$omp end parallel
      deallocate(D, tk)

end subroutine update_vvvv_in_g_h

