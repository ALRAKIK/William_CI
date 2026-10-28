subroutine build_ccsd_w(blk,g,t1,t2,Woooo,Wvvvv,Wvoov,Wvovo)

      use constants_module, only: dp, dpc
      use unitcell_module , only: nfuc, N_cell
      use k_blocks        , only: K_block

      implicit none 

      ! input ! 

      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: t1(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      ! local !

      integer      :: occ(nfuc*N_cell), vir(nfuc*N_cell), no,nv,x,y
      integer      :: p,q,r,s,kp,kq,kr,ks,bp,bq,br,bs,m,e,l,d
      integer      :: fourth_k
      
      complex(dpc) :: w, tau

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

      do kp=1,N_cell
        do kq=1,N_cell
          do kr=1,N_cell
            ks=fourth_k(kp,kq,kr)
            do bp=1,nfuc
              p=(kp-1)*nfuc+bp
              do bq=1,nfuc
                q=(kq-1)*nfuc+bq
                do br=1,nfuc
                  r=(kr-1)*nfuc+br
                  do bs=1,nfuc
                    s=(ks-1)*nfuc+bs
                    if (isocc(p).and.isocc(q).and.isocc(r).and.isocc(s)) then
                      ! W(p,q,r,s): bare (pr|qs).
                      w=eri(p,r,q,s)
                      do x=1,nv
                        e=vir(x)
                        w=w+eri(q,e,p,r)*single(s,e) &
                           +eri(p,e,q,s)*single(r,e)
                        do y=1,nv
                          d=vir(y)
                          tau=double(r,s,e,d)+single(r,e)*single(s,d)
                          w=w+eri(p,e,q,d)*tau
                        end do
                      end do
                      Woooo(bp,bq,br,bs,kp,kq,kr)=w
                    end if
                    if (isvir(p).and.isvir(q).and.isvir(r).and.isvir(s)) then
                      ! PySCF cc_Wvvvv: no T2 or T1*T1 dressing here.
                      w=eri(p,r,q,s)
                      do x=1,no
                        m=occ(x)
                        w=w-eri(m,s,p,r)*single(m,q) &
                           -eri(m,r,q,s)*single(m,p)
                      end do
                      Wvvvv(bp,bq,br,bs,kp,kq,kr)=w
                    end if
                    if (isvir(p).and.isocc(q).and.isocc(r).and.isvir(s)) then
                      ! Wvoov(a,m,i,e): bare (me|ai).
                      w=eri(q,s,p,r)
                      do x=1,nv
                        d=vir(x)
                        w=w+eri(q,s,p,d)*single(r,d)
                      end do
                      do y=1,no
                        l=occ(y)
                        w=w-eri(q,s,l,r)*single(l,p)
                        do x=1,nv
                          d=vir(x)
                          w=w-0.5_dp*eri(l,d,q,s)*double(r,l,d,p) &
                             -0.5_dp*eri(l,s,q,d)*double(r,l,p,d) &
                             -eri(l,d,q,s)*single(r,d)*single(l,p) &
                             +eri(l,d,q,s)*double(r,l,p,d)
                        end do
                      end do
                      Wvoov(bp,bq,br,bs,kp,kq,kr)=w
                    end if
                    if (isvir(p).and.isocc(q).and.isvir(r).and.isocc(s)) then
                      ! Wvovo(a,m,e,i): bare (mi|ae).
                      w=eri(q,s,p,r)
                      do x=1,nv
                        d=vir(x)
                        w=w+eri(q,d,p,r)*single(s,d)
                      end do
                      do y=1,no
                        l=occ(y)
                        w=w-eri(l,r,q,s)*single(l,p)
                        do x=1,nv
                          d=vir(x)
                          w=w-0.5_dp*eri(l,r,q,d)*double(s,l,d,p) &
                             -eri(l,r,q,d)*single(s,d)*single(l,p)
                        end do
                      end do
                      Wvovo(bp,bq,br,bs,kp,kq,kr)=w
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

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