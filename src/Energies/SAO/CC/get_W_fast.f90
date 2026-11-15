subroutine build_ccsd_w_fast(blk,g,ktab,t1,t2,Woooo,Wvoov,Wvovo)

      use constants_module, only: dp, dpc
      use cc_dims         , only: gidx, nslot      ! [Nov 2026] compact g
      use unitcell_module , only: nfuc, N_cell
      use cc_dims         , only: o_hi, v_lo, use_ksym   ! [Nov 2026] use_ksym
      use k_blocks        , only: K_block
      use omp_lib

      implicit none

      ! input !
      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(nfuc,nfuc,nfuc,nfuc,nslot)   ! [Nov 2026] compact g (cc_dims gidx)
      complex(dpc) , intent(in) :: t1(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: t2(o_hi,o_hi,v_lo:nfuc,v_lo:nfuc,N_cell,N_cell,N_cell)
      integer      , intent(in) :: ktab(N_cell,N_cell,N_cell)

      ! output !
      complex(dpc), intent(out) :: Woooo(o_hi,o_hi,o_hi,o_hi,N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: Wvoov(v_lo:nfuc,o_hi,o_hi,v_lo:nfuc,N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: Wvovo(v_lo:nfuc,o_hi,v_lo:nfuc,o_hi,N_cell,N_cell,N_cell)

      ! local !
      integer      :: kOcc(nfuc,N_cell), nOcc(N_cell), kVir(nfuc,N_cell), nVir(N_cell)
      logical      :: isO(nfuc,N_cell), isV(nfuc,N_cell)
      integer      :: k, p, no, nv, nov, nthr
      integer      :: kp,kq,kr,ks,ke,kd,kl,kpr
      integer      :: ii,jj,mm,nn,ll
      integer      :: bp,bq,br,bs,bl,bd,a,i,l,d,m,e,r,s,q
      complex(dpc) :: w
      logical      :: outer, partner
      integer      :: mkp, mkq, mkr
      complex(dpc), parameter :: zero=(0.0_dp,0.0_dp), one=(1.0_dp,0.0_dp)
      complex(dpc), allocatable :: Aoo(:,:,:,:,:), Too(:,:,:,:,:,:), Coo(:,:,:,:,:)
      complex(dpc), allocatable :: TT(:,:,:,:,:), XX(:,:,:,:,:,:), CC(:,:,:,:,:)

      nOcc=0; nVir=0; isO=.false.; isV=.false.
      do k=1,N_cell
        do p=1,nfuc
          if (blk(k)%occupation(p)==2) then
            nOcc(k)=nOcc(k)+1; kOcc(nOcc(k),k)=p; isO(p,k)=.true.
          else if (blk(k)%occupation(p)==0) then
            nVir(k)=nVir(k)+1; kVir(nVir(k),k)=p; isV(p,k)=.true.
          end if
        end do
      end do
      no  = o_hi
      nv  = nfuc - v_lo + 1
      nov = no*nv
      nthr  = omp_get_max_threads()
      outer = (N_cell*N_cell >= nthr)     ! k-pairs in parallel, else threaded BLAS

      !=================================================================!
      ! 1) cheap part: same loops as get_W.f90 without the t2 terms     !
      !=================================================================!
      !$omp parallel default(shared) &
      !$omp private(kp,kq,kr,ks,ii,jj,mm,nn,ll,bp,bq,br,bs,bl,bd,w)

      !$omp do collapse(3) schedule(dynamic)
      do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
        ks=ktab(kp,kq,kr)
        if (use_ksym) then                     ! [Nov 2026] filled below from (-k)
          if (.not.canon2(kp,kq)) cycle
        end if
        Woooo(:,:,:,:,kp,kq,kr)=zero
        do ii=1,nOcc(kp); bp=kOcc(ii,kp)
        do jj=1,nOcc(kq); bq=kOcc(jj,kq)
        do mm=1,nOcc(kr); br=kOcc(mm,kr)
        do nn=1,nOcc(ks); bs=kOcc(nn,ks)
          w=gq(bp,br,bq,bs,kp,kr,kq)
          do ll=1,nVir(ks); bd=kVir(ll,ks)
            w=w+gq(bq,bd,bp,br,kq,ks,kp)*t1(bs,bd,ks)
          end do
          do ll=1,nVir(kr); bd=kVir(ll,kr)
            w=w+gq(bp,bd,bq,bs,kp,kr,kq)*t1(br,bd,kr)
          end do
          Woooo(bp,bq,br,bs,kp,kq,kr)=w
        end do; end do; end do; end do
      end do; end do; end do
      !$omp end do nowait

      !$omp do collapse(3) schedule(dynamic)
      do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
        ks=ktab(kp,kq,kr)
        if (use_ksym) then                     ! [Nov 2026] filled below from (-k)
          if (.not.canon2(kp,kr)) cycle
        end if
        Wvoov(:,:,:,:,kp,kq,kr)=zero
        do ii=1,nVir(kp); bp=kVir(ii,kp)
        do jj=1,nOcc(kq); bq=kOcc(jj,kq)
        do mm=1,nOcc(kr); br=kOcc(mm,kr)
        do nn=1,nVir(ks); bs=kVir(nn,ks)
          w=gq(bq,bs,bp,br,kq,ks,kp)
          do ll=1,nVir(kr); bd=kVir(ll,kr)
            w=w+gq(bq,bs,bp,bd,kq,ks,kp)*t1(br,bd,kr)
          end do
          do ll=1,nOcc(kp); bl=kOcc(ll,kp)
            w=w-gq(bq,bs,bl,br,kq,ks,kp)*t1(bl,bp,kp)
          end do
          Wvoov(bp,bq,br,bs,kp,kq,kr)=w
        end do; end do; end do; end do
      end do; end do; end do
      !$omp end do nowait

      !$omp do collapse(3) schedule(dynamic)
      do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
        ks=ktab(kp,kq,kr)
        if (use_ksym) then                     ! [Nov 2026] filled below from (-k)
          if (.not.canon2(kp,ks)) cycle
        end if
        Wvovo(:,:,:,:,kp,kq,kr)=zero
        do ii=1,nVir(kp); bp=kVir(ii,kp)
        do jj=1,nOcc(kq); bq=kOcc(jj,kq)
        do mm=1,nVir(kr); br=kVir(mm,kr)
        do nn=1,nOcc(ks); bs=kOcc(nn,ks)
          w=gq(bq,bs,bp,br,kq,ks,kp)
          do ll=1,nVir(ks); bd=kVir(ll,ks)
            w=w+gq(bq,bd,bp,br,kq,ks,kp)*t1(bs,bd,ks)
          end do
          do ll=1,nOcc(kp); bl=kOcc(ll,kp)
            w=w-gq(bl,br,bq,bs,kp,kr,kq)*t1(bl,bp,kp)
          end do
          Wvovo(bp,bq,br,bs,kp,kq,kr)=w
        end do; end do; end do; end do
      end do; end do; end do
      !$omp end do

      !$omp end parallel

      ! [Nov 2026] k <-> -k symmetry: X(-kp,-kq,-kr) = X(kp,kq,kr)^*
      if (use_ksym) then
        !$omp parallel do collapse(3) default(shared) private(kp,kq,kr,ks)
        do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
          ks=ktab(kp,kq,kr)
          if (.not.canon2(kp,kq)) Woooo(:,:,:,:,kp,kq,kr) = conjg(Woooo(:,:,:,:,ktab(1,1,kp),ktab(1,1,kq),ktab(1,1,kr)))
          if (.not.canon2(kp,kr)) Wvoov(:,:,:,:,kp,kq,kr) = conjg(Wvoov(:,:,:,:,ktab(1,1,kp),ktab(1,1,kq),ktab(1,1,kr)))
          if (.not.canon2(kp,ks)) Wvovo(:,:,:,:,kp,kq,kr) = conjg(Wvovo(:,:,:,:,ktab(1,1,kp),ktab(1,1,kq),ktab(1,1,kr)))
        end do; end do; end do
        !$omp end parallel do
      end if

      if (no<=0 .or. nv<=0) return

      !=================================================================!
      ! 2) Woooo t2 term, fixed (kp,kq):                                 !
      !    W(p,q,r,s) += sum_{ke,e,d} g(p,e,q,d;kp,ke,kq) tau(r,s,e,d)   !
      !    tau = t2(r,s,e,d;kr,ks,ke) + t1(r,e,kr) t1(s,d,ks) [ke==kr]   !
      !=================================================================!
      !$omp parallel do if(outer) collapse(2) schedule(dynamic) default(shared) &
      !$omp private(kp,kq,kr,ks,ke,kd,p,q,r,s,e,d,Aoo,Too,Coo,partner,mkp,mkq,mkr)
      do kp=1,N_cell
        do kq=1,N_cell
          if (use_ksym) then                   ! [Nov 2026] block (-kp,-kq) added as conjugate
            if (.not.canon2(kp,kq)) cycle
          end if
          partner = use_ksym .and. .not.(ktab(1,1,kp)==kp .and. ktab(1,1,kq)==kq)
          mkp=ktab(1,1,kp); mkq=ktab(1,1,kq)
          allocate(Aoo(no,no,v_lo:nfuc,v_lo:nfuc,N_cell), Too(no,no,N_cell,v_lo:nfuc,v_lo:nfuc,N_cell))
          allocate(Coo(no,no,no,no,N_cell))
          do ke=1,N_cell
            kd=ktab(kp,kq,ke)
            do d=v_lo,nfuc; do e=v_lo,nfuc
              if (isV(e,ke) .and. isV(d,kd)) then
                do q=1,no
                  Aoo(:,q,e,d,ke)=gv1(1,no,e,q,d,kp,ke,kq)
                end do
              else
                Aoo(:,:,e,d,ke)=zero
              end if
            end do; end do
            do kr=1,N_cell
              ks=ktab(kp,kq,kr)
              do d=v_lo,nfuc; do e=v_lo,nfuc
                Too(:,:,kr,e,d,ke)=t2(:,:,e,d,kr,ks,ke)
                if (ke==kr) then
                  do s=1,no; do r=1,no
                    Too(r,s,kr,e,d,ke)=Too(r,s,kr,e,d,ke)+t1(r,e,kr)*t1(s,d,ks)
                  end do; end do
                end if
              end do; end do
            end do
          end do
          call zgemm('N','T', no*no, no*no*N_cell, nv*nv*N_cell, one, Aoo, no*no, &
                     Too, no*no*N_cell, zero, Coo, no*no)
          do kr=1,N_cell
            ks=ktab(kp,kq,kr)
            do s=1,no; if (.not.isO(s,ks)) cycle
            do r=1,no; if (.not.isO(r,kr)) cycle
            do q=1,no; if (.not.isO(q,kq)) cycle
            do p=1,no; if (.not.isO(p,kp)) cycle
              Woooo(p,q,r,s,kp,kq,kr)=Woooo(p,q,r,s,kp,kq,kr)+Coo(p,q,r,s,kr)
              if (partner) Woooo(p,q,r,s,mkp,mkq,ktab(1,1,kr)) = &
                           Woooo(p,q,r,s,mkp,mkq,ktab(1,1,kr)) + conjg(Coo(p,q,r,s,kr))
            end do; end do; end do; end do
          end do
          deallocate(Aoo,Too,Coo)
        end do
      end do
      !$omp end parallel do

      !=================================================================!
      ! 3) Wvoov t2 terms, fixed (kp,kr,kl), all kq at once:             !
      !    W(a,m,i,e) += sum_{l,d} T1(a,i,l,d) g(l,d,m,e;kl,kd,kq)       !
      !                          + T2(a,i,l,d) g(l,e,m,d;kl,ks,kq)       !
      !    T1 = -1/2 t2(i,l,d,a;kr,kl,kd) + t2(i,l,a,d;kr,kl,kp)         !
      !         - t1(i,d,kr) t1(l,a,kp) [kl==kp]                         !
      !    T2 = -1/2 t2(i,l,a,d;kr,kl,kp)                                !
      !=================================================================!
      !$omp parallel do if(outer) collapse(2) schedule(dynamic) default(shared) &
      !$omp private(kp,kr,kq,ks,kl,kd,a,i,l,d,m,e,TT,XX,CC,partner,mkp,mkq,mkr)
      do kp=1,N_cell
        do kr=1,N_cell
          if (use_ksym) then                   ! [Nov 2026] block (-kp,-kr) added as conjugate
            if (.not.canon2(kp,kr)) cycle
          end if
          partner = use_ksym .and. .not.(ktab(1,1,kp)==kp .and. ktab(1,1,kr)==kr)
          mkp=ktab(1,1,kp); mkr=ktab(1,1,kr)
          allocate(TT(v_lo:nfuc,no,no,v_lo:nfuc,2), XX(no,v_lo:nfuc,2,no,v_lo:nfuc,N_cell))
          allocate(CC(v_lo:nfuc,no,no,v_lo:nfuc,N_cell))
          CC=zero
          do kl=1,N_cell
            kd=ktab(kr,kl,kp)
            do d=v_lo,nfuc; do l=1,no
              if (isO(l,kl) .and. isV(d,kd)) then
                do i=1,no; do a=v_lo,nfuc
                  TT(a,i,l,d,1)=-0.5_dp*t2(i,l,d,a,kr,kl,kd)+t2(i,l,a,d,kr,kl,kp)
                  TT(a,i,l,d,2)=-0.5_dp*t2(i,l,a,d,kr,kl,kp)
                end do; end do
                if (kl==kp) then
                  do i=1,no; do a=v_lo,nfuc
                    TT(a,i,l,d,1)=TT(a,i,l,d,1)-t1(i,d,kr)*t1(l,a,kp)
                  end do; end do
                end if
              else
                TT(:,:,l,d,:)=zero
              end if
            end do; end do
            do kq=1,N_cell
              ks=ktab(kp,kq,kr)
              do e=v_lo,nfuc; do m=1,no
                do d=v_lo,nfuc
                  XX(:,d,1,m,e,kq)=gv1(1,no,d,m,e,kl,kd,kq)
                  XX(:,d,2,m,e,kq)=gv1(1,no,e,m,d,kl,ks,kq)
                end do
              end do; end do
            end do
            ! TT(:,:,l,d,1:2) as rows (a,i), cols (l,d,t) ; XX rows (l,d,t)
            ! the row order of XX is (l,d,t) and the col order of TT is (l,d,t)
            call zgemm('N','N', nov, nov*N_cell, 2*nov, one, TT, nov, XX, 2*nov, one, CC, nov)
          end do
          do kq=1,N_cell
            ks=ktab(kp,kq,kr)
            do e=v_lo,nfuc; if (.not.isV(e,ks)) cycle
            do m=1,no;      if (.not.isO(m,kq)) cycle
            do i=1,no;      if (.not.isO(i,kr)) cycle
            do a=v_lo,nfuc; if (.not.isV(a,kp)) cycle
              Wvoov(a,m,i,e,kp,kq,kr)=Wvoov(a,m,i,e,kp,kq,kr)+CC(a,i,m,e,kq)
              if (partner) Wvoov(a,m,i,e,mkp,ktab(1,1,kq),mkr) = &
                           Wvoov(a,m,i,e,mkp,ktab(1,1,kq),mkr) + conjg(CC(a,i,m,e,kq))
            end do; end do; end do; end do
          end do
          deallocate(TT,XX,CC)
        end do
      end do
      !$omp end parallel do

      !=================================================================!
      ! 4) Wvovo t2 terms, fixed (kp,ks,kl), all kq at once (kr=kp+kq-ks)!
      !    W(a,m,e,i) += sum_{l,d} T3(a,i,l,d) g(l,e,m,d;kl,kr,kq)       !
      !    T3 = -1/2 t2(i,l,d,a;ks,kl,kd) - t1(i,d,ks) t1(l,a,kp) [kl==kp]!
      !=================================================================!
      !$omp parallel do if(outer) collapse(2) schedule(dynamic) default(shared) &
      !$omp private(kp,ks,kq,kr,kl,kd,a,i,l,d,m,e,TT,XX,CC,partner,mkp,mkq,mkr)
      do kp=1,N_cell
        do ks=1,N_cell
          if (use_ksym) then                   ! [Nov 2026] block (-kp,-ks) added as conjugate
            if (.not.canon2(kp,ks)) cycle
          end if
          partner = use_ksym .and. .not.(ktab(1,1,kp)==kp .and. ktab(1,1,ks)==ks)
          mkp=ktab(1,1,kp)
          allocate(TT(v_lo:nfuc,no,no,v_lo:nfuc,1), XX(no,v_lo:nfuc,1,no,v_lo:nfuc,N_cell))
          allocate(CC(v_lo:nfuc,no,no,v_lo:nfuc,N_cell))
          CC=zero
          do kl=1,N_cell
            kd=ktab(ks,kl,kp)
            do d=v_lo,nfuc; do l=1,no
              if (isO(l,kl) .and. isV(d,kd)) then
                do i=1,no; do a=v_lo,nfuc
                  TT(a,i,l,d,1)=-0.5_dp*t2(i,l,d,a,ks,kl,kd)
                end do; end do
                if (kl==kp) then
                  do i=1,no; do a=v_lo,nfuc
                    TT(a,i,l,d,1)=TT(a,i,l,d,1)-t1(i,d,ks)*t1(l,a,kp)
                  end do; end do
                end if
              else
                TT(:,:,l,d,1)=zero
              end if
            end do; end do
            do kq=1,N_cell
              kr=ktab(kp,kq,ks)
              do e=v_lo,nfuc; do m=1,no
                do d=v_lo,nfuc
                  XX(:,d,1,m,e,kq)=gv1(1,no,e,m,d,kl,kr,kq)
                end do
              end do; end do
            end do
            call zgemm('N','N', nov, nov*N_cell, nov, one, TT, nov, XX, nov, one, CC, nov)
          end do
          do kq=1,N_cell
            kr=ktab(kp,kq,ks)
            do i=1,no;      if (.not.isO(i,ks)) cycle
            do e=v_lo,nfuc; if (.not.isV(e,kr)) cycle
            do m=1,no;      if (.not.isO(m,kq)) cycle
            do a=v_lo,nfuc; if (.not.isV(a,kp)) cycle
              Wvovo(a,m,e,i,kp,kq,kr)=Wvovo(a,m,e,i,kp,kq,kr)+CC(a,i,m,e,kq)
              if (partner) Wvovo(a,m,e,i,mkp,ktab(1,1,kq),ktab(1,1,kr)) = &
                           Wvovo(a,m,e,i,mkp,ktab(1,1,kq),ktab(1,1,kr)) + conjg(CC(a,i,m,e,kq))
            end do; end do; end do; end do
          end do
          deallocate(TT,XX,CC)
        end do
      end do
      !$omp end parallel do

contains
      ! [Nov 2026] (k1,k2) <= (-k1,-k2) in lexicographic order
      logical function canon2(k1,k2)
        integer, intent(in) :: k1,k2
        canon2 = (k1<ktab(1,1,k1)) .or. (k1==ktab(1,1,k1) .and. k2<=ktab(1,1,k2))
      end function
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
end subroutine build_ccsd_w_fast
