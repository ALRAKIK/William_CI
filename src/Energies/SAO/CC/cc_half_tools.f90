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
      use cc_dims         , only: gidx, nslot      ! [Nov 2026] compact g
      use unitcell_module , only: nfuc, N_cell
      use cc_dims       , only: o_hi, v_lo
      use k_blocks        , only: K_block

      implicit none 

      ! input ! 

      type(K_block), intent(in)  :: blk(N_cell)
      complex(dpc) , intent(in)  :: g(nfuc,nfuc,nfuc,nfuc,nslot)   ! [Nov 2026] compact g (cc_dims gidx)
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

      use constants_module, only: dp, dpc
      use cc_dims         , only: gidx, nslot      ! [Nov 2026] compact g
      use unitcell_module , only: nfuc, N_cell
      use cc_dims       , only: o_hi, v_lo
      use k_blocks        , only: K_block

      implicit none

      ! input ! 

      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(nfuc,nfuc,nfuc,nfuc,nslot)   ! [Nov 2026] compact g (cc_dims gidx)
      complex(dpc) , intent(in) :: t1(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: t2(o_hi,o_hi,v_lo:nfuc,v_lo:nfuc,N_cell,N_cell,N_cell)
      integer      , intent(in) :: ktab(N_cell,N_cell,N_cell)
      ! local ! 

      integer      :: k,km,kn,ke,kf,m,n,i,a,e,f
      integer      :: fourth_k
      complex(dpc) :: v,tau

      integer :: kOcc(nfuc,N_cell), nOcc(N_cell)
      integer :: kVir(nfuc,N_cell), nVir(N_cell)
      
      integer :: ii , mm , nn , ee , ff, aa
      integer :: p 
      integer :: kp , kq , kr

      ! output !

      complex(dpc), intent(out) :: Foo(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: Fvv(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: Fov(nfuc,nfuc,N_cell)



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

      Foo=cmplx(0.0_dp,0.0_dp,kind=dpc)
      Fvv=cmplx(0.0_dp,0.0_dp,kind=dpc)
      Fov=cmplx(0.0_dp,0.0_dp,kind=dpc)

      ! Foo(m,i,k) = f_mi + sum_nef (2(me|nf)-(mf|ne))*tau_in^ef
      do k=1,N_cell
        do ii=1,Nocc(k)
            i=kOcc(ii,k)
          do mm=1,nOcc(k)
              m=kOcc(mm,k)
              if (m==i) Foo(m,i,k)=cmplx(blk(k)%energy(i),0.0_dp,kind=dpc)
              do kn=1,N_cell
                do ke=1,N_cell
                  kf=ktab(k,kn,ke)
                  do nn=1,nOcc(kn)
                    n=kOcc(nn,kn)
                    do ee=1,nVir(ke)
                      e=kVir(ee,ke)
                      do ff=1,nVir(kf)
                        f=kVir(ff,kf)
                        tau=t2(i,n,e,f,k,kn,ke)
                        if (ke==k .and. kf==kn) tau=tau+t1(i,e,k)*t1(n,f,kn)
                        v=2.0_dp*gq(m,e,n,f,k,ke,kn)-gq(m,f,n,e,k,kf,kn)
                        Foo(m,i,k)=Foo(m,i,k)+v*tau
                      end do
                    end do
                  end do
                end do
              end do
          end do
        end do
      end do

      ! Fvv(a,e,k) = f_ae - sum_mnf (2(me|nf)-(mf|ne))*tau_mn^af

      do k=1,N_cell
        do ee=1,nVir(k)
            e = kVir(ee,k)
          do aa=1,nVir(k)
            a = kVir(aa,k)
            if (a==e) Fvv(a,e,k)=cmplx(blk(k)%energy(a),0.0_dp,kind=dpc)
            do km=1,N_cell
              do kn=1,N_cell
                kf=ktab(km,kn,k)
                do mm=1,nOcc(km)
                  m = kOcc(mm,km)
                  do nn=1,nOcc(kn)
                    n = kOcc(nn,kn)
                    do ff=1,nVir(kf)
                      f = kVir(ff,kf)
                      tau=t2(m,n,a,f,km,kn,k)
                      if (km==k .and. kf==kn) tau=tau+t1(m,a,km)*t1(n,f,kn)
                      v=2.0_dp*gq(m,e,n,f,km,k,kn)-gq(m,f,n,e,km,kf,kn)
                      Fvv(a,e,k)=Fvv(a,e,k)-v*tau
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

      ! Fov(m,e,k) = sum_nf (2(me|nf)-(mf|ne))*t1(n,f,kn).
      ! Translational symmetry makes ALL three F blocks diagonal in k.

      do k=1,N_cell
        do ee=1,nVir(k)
            e = kVir(ee,k)
          do mm=1,nOcc(k)
            m = kOcc(mm,k)
            do kn=1,N_cell
              do nn=1,nOcc(kn)
                n = kOcc(nn,kn)
                do ff=1,nVir(kn)
                  f = kVir(ff,kn)
                  v=2.0_dp*gq(m,e,n,f,k,k,kn)-gq(m,f,n,e,k,kn,kn)
                  Fov(m,e,k)=Fov(m,e,k)+v*t1(n,f,kn)
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

end subroutine build_ccsd_f_h

subroutine build_ccsd_l_h(blk,g,t1,Foo,Fvv,Loo,Lvv)
      
      use constants_module, only: dp, dpc
      use cc_dims         , only: gidx, nslot      ! [Nov 2026] compact g
      use unitcell_module , only: nfuc, N_cell
      use k_blocks        , only: K_block

      implicit none 

      ! input ! 

      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(nfuc,nfuc,nfuc,nfuc,nslot)   ! [Nov 2026] compact g (cc_dims gidx)
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

      use constants_module, only: dp, dpc
      use cc_dims         , only: gidx, nslot      ! [Nov 2026] compact g
      use unitcell_module , only: nfuc, N_cell
      use k_blocks        , only: K_block

      implicit none

      type(K_block), intent(in)    :: blk(N_cell)
      complex(dpc) , intent(inout) :: g(nfuc,nfuc,nfuc,nfuc,nslot)   ! [Nov 2026] compact g (cc_dims gidx)
      integer      , intent(in)    :: ktab(N_cell,N_cell,N_cell)
      complex(dpc) , intent(in)    :: dt1(nfuc,nfuc,N_cell)

      integer      :: kOcc(nfuc,N_cell), nOcc(N_cell), kVir(nfuc,N_cell), nVir(N_cell)
      integer      :: k, p, kp, kq, kr, ks, ii, jj, mm, nn, ll, bp, bq, br, bs, bl
      complex(dpc) :: w
      integer      :: ic
      complex(dpc), parameter :: zero=(0.0_dp,0.0_dp)

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

      !$omp parallel do collapse(3) schedule(dynamic) default(shared) &
      !$omp private(kp,kq,kr,ks,ii,jj,mm,nn,ll,bp,bq,br,bs,bl,w,ic)
      do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
        ks=ktab(kp,kq,kr)
        ic=gidx(kp,kr,kq)                     ! [Nov 2026] stored block of Wvvvv(kp,kq,kr)
        if (ic <= 0) cycle                    ! conjugate of a stored block: nothing to do
        do ii=1,nVir(kp); bp=kVir(ii,kp)
        do jj=1,nVir(kq); bq=kVir(jj,kq)
        do mm=1,nVir(kr); br=kVir(mm,kr)
        do nn=1,nVir(ks); bs=kVir(nn,ks)
          w=zero
          do ll=1,nOcc(kq); bl=kOcc(ll,kq)
            w=w+gq(bl,bs,bp,br,kq,ks,kp)*dt1(bl,bq,kq)
          end do
          do ll=1,nOcc(kp); bl=kOcc(ll,kp)
            w=w+gq(bl,br,bq,bs,kp,kr,kq)*dt1(bl,bp,kp)
          end do
          g(bp,br,bq,bs,ic)=g(bp,br,bq,bs,ic)-w
        end do; end do; end do; end do
      end do; end do; end do
      !$omp end parallel do
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

end subroutine update_vvvv_in_g_h

