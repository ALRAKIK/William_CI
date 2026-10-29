subroutine build_ccsd_w(blk,g,ktab,t1,t2,Woooo,Wvvvv,Wvoov,Wvovo)

      use constants_module, only: dp, dpc
      use unitcell_module , only: nfuc, N_cell
      use k_blocks        , only: K_block

      implicit none

      ! input !
      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: t1(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      integer      , intent(in) :: ktab(N_cell,N_cell,N_cell)

      ! output !

      complex(dpc), intent(out) :: Woooo(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: Wvvvv(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: Wvoov(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: Wvovo(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      ! local !

      integer      :: kOcc(nfuc,N_cell), nOcc(N_cell), kVir(nfuc,N_cell), nVir(N_cell)
      integer      :: k, p
      integer      :: kp,kq,kr,ks,ke,kd,kl
      integer      :: ii,jj,mm,nn,ee,dd,ll
      integer      :: bp,bq,br,bs,be,bd,bl
      complex(dpc) :: w, gld, tps
      complex(dpc), parameter :: zero=(0.0_dp,0.0_dp)

      ! Layouts (unchanged):
      ! Woooo(m,n,i,j,km,kn,ki), Wvvvv(a,b,e,f,ka,kb,ke),
      ! Wvoov(a,m,i,e,ka,km,ki), Wvovo(a,m,e,i,ka,km,ke), ks=ktab(kp,kq,kr)
      ! g(b1,b2,b3,b4,k1,k2,k3), t1(bi,ba,k), t2(bi,bj,ba,bb,ki,kj,ka)

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

      !$omp parallel default(shared) &
      !$omp private(kp,kq,kr,ks,ke,kd,kl,ii,jj,mm,nn,ee,dd,ll, &
      !$omp         bp,bq,br,bs,be,bd,bl,w,gld,tps)

      ! ---- Woooo: all occupied ----
      !$omp do collapse(3) schedule(dynamic)
      do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
        ks=ktab(kp,kq,kr)
        Woooo(:,:,:,:,kp,kq,kr)=zero
        do ii=1,nOcc(kp); bp=kOcc(ii,kp)
        do jj=1,nOcc(kq); bq=kOcc(jj,kq)
        do mm=1,nOcc(kr); br=kOcc(mm,kr)
        do nn=1,nOcc(ks); bs=kOcc(nn,ks)
          w=g(bp,br,bq,bs,kp,kr,kq)
          do ee=1,nVir(ks); be=kVir(ee,ks)
            w=w+g(bq,be,bp,br,kq,ks,kp)*t1(bs,be,ks)
          end do
          do ee=1,nVir(kr); be=kVir(ee,kr)
            w=w+g(bp,be,bq,bs,kp,kr,kq)*t1(br,be,kr)
          end do
          do ke=1,N_cell
            kd=ktab(kr,ks,ke)
            do ee=1,nVir(ke); be=kVir(ee,ke)
              do dd=1,nVir(kd); bd=kVir(dd,kd)
                tps=t2(br,bs,be,bd,kr,ks,ke)
                if (ke==kr) tps=tps+t1(br,be,kr)*t1(bs,bd,ks)
                w=w+g(bp,be,bq,bd,kp,ke,kq)*tps
              end do
            end do
          end do
          Woooo(bp,bq,br,bs,kp,kq,kr)=w
        end do; end do; end do; end do
      end do; end do; end do
      !$omp end do nowait

      ! ---- Wvvvv: all virtual ----
      !$omp do collapse(3) schedule(dynamic)
      do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
        ks=ktab(kp,kq,kr)
        Wvvvv(:,:,:,:,kp,kq,kr)=zero
        do ii=1,nVir(kp); bp=kVir(ii,kp)
        do jj=1,nVir(kq); bq=kVir(jj,kq)
        do mm=1,nVir(kr); br=kVir(mm,kr)
        do nn=1,nVir(ks); bs=kVir(nn,ks)
          w=g(bp,br,bq,bs,kp,kr,kq)
          do ll=1,nOcc(kq); bl=kOcc(ll,kq)
            w=w-g(bl,bs,bp,br,kq,ks,kp)*t1(bl,bq,kq)
          end do
          do ll=1,nOcc(kp); bl=kOcc(ll,kp)
            w=w-g(bl,br,bq,bs,kp,kr,kq)*t1(bl,bp,kp)
          end do
          Wvvvv(bp,bq,br,bs,kp,kq,kr)=w
        end do; end do; end do; end do
      end do; end do; end do
      !$omp end do nowait

      ! ---- Wvoov: p,s virtual ; q,r occupied ----
      !$omp do collapse(3) schedule(dynamic)
      do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
        ks=ktab(kp,kq,kr)
        Wvoov(:,:,:,:,kp,kq,kr)=zero
        do ii=1,nVir(kp); bp=kVir(ii,kp)
        do jj=1,nOcc(kq); bq=kOcc(jj,kq)
        do mm=1,nOcc(kr); br=kOcc(mm,kr)
        do nn=1,nVir(ks); bs=kVir(nn,ks)
          w=g(bq,bs,bp,br,kq,ks,kp)
          do dd=1,nVir(kr); bd=kVir(dd,kr)
            w=w+g(bq,bs,bp,bd,kq,ks,kp)*t1(br,bd,kr)
          end do
          do ll=1,nOcc(kp); bl=kOcc(ll,kp)
            w=w-g(bq,bs,bl,br,kq,ks,kp)*t1(bl,bp,kp)
          end do
          do kl=1,N_cell
            kd=ktab(kr,kl,kp)
            do ll=1,nOcc(kl); bl=kOcc(ll,kl)
              do dd=1,nVir(kd); bd=kVir(dd,kd)
                gld=g(bl,bd,bq,bs,kl,kd,kq)
                w=w-0.5_dp*gld*t2(br,bl,bd,bp,kr,kl,kd)                    &
                   -0.5_dp*g(bl,bs,bq,bd,kl,ks,kq)*t2(br,bl,bp,bd,kr,kl,kp) &
                   +gld*t2(br,bl,bp,bd,kr,kl,kp)
                if (kl==kp) w=w-gld*t1(br,bd,kr)*t1(bl,bp,kp)
              end do
            end do
          end do
          Wvoov(bp,bq,br,bs,kp,kq,kr)=w
        end do; end do; end do; end do
      end do; end do; end do
      !$omp end do nowait

      ! ---- Wvovo: p,r virtual ; q,s occupied ----
      !$omp do collapse(3) schedule(dynamic)
      do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
        ks=ktab(kp,kq,kr)
        Wvovo(:,:,:,:,kp,kq,kr)=zero
        do ii=1,nVir(kp); bp=kVir(ii,kp)
        do jj=1,nOcc(kq); bq=kOcc(jj,kq)
        do mm=1,nVir(kr); br=kVir(mm,kr)
        do nn=1,nOcc(ks); bs=kOcc(nn,ks)
          w=g(bq,bs,bp,br,kq,ks,kp)
          do dd=1,nVir(ks); bd=kVir(dd,ks)
            w=w+g(bq,bd,bp,br,kq,ks,kp)*t1(bs,bd,ks)
          end do
          do ll=1,nOcc(kp); bl=kOcc(ll,kp)
            w=w-g(bl,br,bq,bs,kp,kr,kq)*t1(bl,bp,kp)
          end do
          do kl=1,N_cell
            kd=ktab(ks,kl,kp)
            do ll=1,nOcc(kl); bl=kOcc(ll,kl)
              do dd=1,nVir(kd); bd=kVir(dd,kd)
                gld=g(bl,br,bq,bd,kl,kr,kq)
                w=w-0.5_dp*gld*t2(bs,bl,bd,bp,ks,kl,kd)
                if (kl==kp) w=w-gld*t1(bs,bd,ks)*t1(bl,bp,kp)
              end do
            end do
          end do
          Wvovo(bp,bq,br,bs,kp,kq,kr)=w
        end do; end do; end do; end do
      end do; end do; end do
      !$omp end do

      !$omp end parallel
end subroutine build_ccsd_w
