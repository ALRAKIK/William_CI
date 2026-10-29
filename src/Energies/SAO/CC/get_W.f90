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
      ! local !

      integer      :: occ(nfuc*N_cell), vir(nfuc*N_cell), no,nv,x
      integer      :: p,q,r,s,kp,kq,kr,ks,bp,m,e,l,d
      integer      :: fourth_k
      
      complex(dpc) :: w, tau

      integer      :: kOcc(nfuc,N_cell), nOcc(N_cell)
      integer      :: kVir(nfuc,N_cell), nVir(N_cell)
      
      integer      :: k , ii , jj , mm , nn , ke , kd , kl , ee , dd ,ll 


      ! output !

      
      complex(dpc), intent(out) :: Woooo(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: Wvvvv(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: Wvoov(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: Wvovo(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      

      ! fourth_k MUST return k1+k2-k3, as in your build_ccsd_f.
      ! All W arrays use their literal orbital order and three first momenta:
      ! Woooo(m,n,i,j,km,kn,ki), Wvvvv(a,b,e,f,ka,kb,ke),
      ! Wvoov(a,m,i,e,ka,km,ki), Wvovo(a,m,e,i,ka,km,ke).
      ! In EVERY case ks=fourth_k(kp,kq,kr).

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

      no=0; nv=0

      do kp=1,N_cell
        do bp=1,nfuc
          x=(kp-1)*nfuc+bp
          if (blk(kp)%occupation(bp)==2) then
            no=no+1; occ(no)=x
          else if (blk(kp)%occupation(bp)==0) then
            nv=nv+1; vir(nv)=x
          end if
        end do
      end do

      Woooo=cmplx(0.0_dp,0.0_dp,kind=dpc)
      Wvvvv=cmplx(0.0_dp,0.0_dp,kind=dpc)
      Wvoov=cmplx(0.0_dp,0.0_dp,kind=dpc)
      Wvovo=cmplx(0.0_dp,0.0_dp,kind=dpc)

      ! ---- Woooo: all occupied ----
      do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
        ks=ktab(kp,kq,kr)
        do ii=1,nOcc(kp); do jj=1,nOcc(kq); do mm=1,nOcc(kr); do nn=1,nOcc(ks)
          p=(kp-1)*nfuc+kOcc(ii,kp); q=(kq-1)*nfuc+kOcc(jj,kq)
          r=(kr-1)*nfuc+kOcc(mm,kr); s=(ks-1)*nfuc+kOcc(nn,ks)
          w=eri(p,r,q,s)
          do ee=1,nVir(ks)
            e=(ks-1)*nfuc+kVir(ee,ks)
            w=w+eri(q,e,p,r)*single(s,e)
          end do
          do ee=1,nVir(kr)
            e=(kr-1)*nfuc+kVir(ee,kr)
            w=w+eri(p,e,q,s)*single(r,e)
          end do
          do ke=1,N_cell
            kd=ktab(kr,ks,ke)
            do ee=1,nVir(ke)
              e=(ke-1)*nfuc+kVir(ee,ke)
              do dd=1,nVir(kd)
                d=(kd-1)*nfuc+kVir(dd,kd)
                tau=double(r,s,e,d)+single(r,e)*single(s,d)
                w=w+eri(p,e,q,d)*tau
              end do
            end do
          end do
          Woooo(band(p),band(q),band(r),band(s),kp,kq,kr)=w
        end do; end do; end do; end do
      end do; end do; end do

      ! ---- Wvvvv: all virtual ----
      do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
        ks=ktab(kp,kq,kr)
        do ii=1,nVir(kp); do jj=1,nVir(kq); do mm=1,nVir(kr); do nn=1,nVir(ks)
          p=(kp-1)*nfuc+kVir(ii,kp); q=(kq-1)*nfuc+kVir(jj,kq)
          r=(kr-1)*nfuc+kVir(mm,kr); s=(ks-1)*nfuc+kVir(nn,ks)
          w=eri(p,r,q,s)
          do ee=1,nOcc(kq)
            m=(kq-1)*nfuc+kOcc(ee,kq)
            w=w-eri(m,s,p,r)*single(m,q)
          end do
          do ee=1,nOcc(kp)
            m=(kp-1)*nfuc+kOcc(ee,kp)
            w=w-eri(m,r,q,s)*single(m,p)
          end do
          Wvvvv(band(p),band(q),band(r),band(s),kp,kq,kr)=w
        end do; end do; end do; end do
      end do; end do; end do

      ! ---- Wvoov: p,s virtual ; q,r occupied ----
      do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
        ks=ktab(kp,kq,kr)
        do ii=1,nVir(kp); do jj=1,nOcc(kq); do mm=1,nOcc(kr); do nn=1,nVir(ks)
          p=(kp-1)*nfuc+kVir(ii,kp); q=(kq-1)*nfuc+kOcc(jj,kq)
          r=(kr-1)*nfuc+kOcc(mm,kr); s=(ks-1)*nfuc+kVir(nn,ks)
          w=eri(q,s,p,r)
          do dd=1,nVir(kr)
            d=(kr-1)*nfuc+kVir(dd,kr)
            w=w+eri(q,s,p,d)*single(r,d)
          end do
          do ll=1,nOcc(kp)
            l=(kp-1)*nfuc+kOcc(ll,kp)
            w=w-eri(q,s,l,r)*single(l,p)
          end do
          do kl=1,N_cell
            kd=ktab(kr,kl,kp)
            do ll=1,nOcc(kl)
              l=(kl-1)*nfuc+kOcc(ll,kl)
              do dd=1,nVir(kd)
                d=(kd-1)*nfuc+kVir(dd,kd)
                w=w-0.5_dp*eri(l,d,q,s)*double(r,l,d,p) &
                   -0.5_dp*eri(l,s,q,d)*double(r,l,p,d) &
                   -eri(l,d,q,s)*single(r,d)*single(l,p) &
                   +eri(l,d,q,s)*double(r,l,p,d)
              end do
            end do
          end do
          Wvoov(band(p),band(q),band(r),band(s),kp,kq,kr)=w
        end do; end do; end do; end do
      end do; end do; end do

      ! ---- Wvovo: p,r virtual ; q,s occupied ----
      do kp=1,N_cell; do kq=1,N_cell; do kr=1,N_cell
        ks=ktab(kp,kq,kr)
        do ii=1,nVir(kp); do jj=1,nOcc(kq); do mm=1,nVir(kr); do nn=1,nOcc(ks)
          p=(kp-1)*nfuc+kVir(ii,kp); q=(kq-1)*nfuc+kOcc(jj,kq)
          r=(kr-1)*nfuc+kVir(mm,kr); s=(ks-1)*nfuc+kOcc(nn,ks)
          w=eri(q,s,p,r)
          do dd=1,nVir(ks)
            d=(ks-1)*nfuc+kVir(dd,ks)
            w=w+eri(q,d,p,r)*single(s,d)
          end do
          do ll=1,nOcc(kp)
            l=(kp-1)*nfuc+kOcc(ll,kp)
            w=w-eri(l,r,q,s)*single(l,p)
          end do
          do kl=1,N_cell
            kd=ktab(ks,kl,kp)
            do ll=1,nOcc(kl)
              l=(kl-1)*nfuc+kOcc(ll,kl)
              do dd=1,nVir(kd)
                d=(kd-1)*nfuc+kVir(dd,kd)
                w=w-0.5_dp*eri(l,r,q,d)*double(s,l,d,p) &
                   -eri(l,r,q,d)*single(s,d)*single(l,p)
              end do
            end do
          end do
          Wvovo(band(p),band(q),band(r),band(s),kp,kq,kr)=w
        end do; end do; end do; end do
      end do; end do; end do

      contains

      integer function mom(a)
        integer, intent(in) :: a
        mom=(a-1)/nfuc+1
      end function mom

      integer function band(a)
        integer, intent(in) :: a
        band=modulo(a-1,nfuc)+1
      end function band

      logical function isocc(a)
        integer, intent(in) :: a
        isocc=blk(mom(a))%occupation(band(a))==2
      end function isocc

      logical function isvir(a)
        integer, intent(in) :: a
        isvir=blk(mom(a))%occupation(band(a))==0
      end function isvir

      function eri(a,b,c,k) result(v)
        integer, intent(in) :: a,b,c,k
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        ! Chemist momentum: kd=ka-kb+kc.
        if (mom(k)/=fourth_k(mom(a),mom(c),mom(b))) return
        v=g(band(a),band(b),band(c),band(k),mom(a),mom(b),mom(c))
      end function eri

      function single(i,a) result(v)
        integer, intent(in) :: i,a
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (mom(i)/=mom(a)) return
        v=t1(band(i),band(a),mom(i))
      end function single

      function double(i,j,a,b) result(v)
        integer, intent(in) :: i,j,a,b
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (mom(b)/=fourth_k(mom(i),mom(j),mom(a))) return
        v=t2(band(i),band(j),band(a),band(b),mom(i),mom(j),mom(a))
      end function double

end subroutine build_ccsd_w