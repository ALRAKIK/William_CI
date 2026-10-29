subroutine build_ccsd_residuals(blk,g,ktab,t1,t2,Foo,Fvv,Fov,Loo,Lvv, &
                                 Woooo,Wvvvv,Wvoov,Wvovo,R1,R2)

      use constants_module, only: dp, dpc
      use unitcell_module , only: nfuc, N_cell
      use k_blocks        , only: K_block

      implicit none

      ! input !
      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: t1(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: Foo(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: Fvv(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: Fov(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: Loo(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: Lvv(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: Woooo(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: Wvvvv(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: Wvoov(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: Wvovo(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      integer      , intent(in) :: ktab(N_cell,N_cell,N_cell)
      ! local !
      integer      :: occ(nfuc*N_cell), vir(nfuc*N_cell)
      integer      :: no,nv,p,k,i,j,a,b,c,d,l,m
      integer      :: ii,jj,aa,bb,cc,dd,ll,mm
      complex(dpc) :: v
      integer      :: kOcc(nfuc,N_cell), nOcc(N_cell)
      integer      :: kVir(nfuc,N_cell), nVir(N_cell)

      integer      :: momA(nfuc*N_cell), bandA(nfuc*N_cell)
      integer      :: ki,kj,ka,kb, kc ,kd ,km , kl
      ! output !
      complex(dpc), intent(out) :: R1(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: R2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      ! ---- setup: lists, lookup tables, momentum table ----
      
      nOcc=0; nVir=0
      no=0;  nv=0
      do k=1,N_cell
        do p=1,nfuc
          momA((k-1)*nfuc+p)=k
          bandA((k-1)*nfuc+p)=p
          if (blk(k)%occupation(p)==2) then
            nOcc(k)=nOcc(k)+1; kOcc(nOcc(k),k)=p
            no=no+1;  occ(no)=(k-1)*nfuc+p
          else if (blk(k)%occupation(p)==0) then
            nVir(k)=nVir(k)+1; kVir(nVir(k),k)=p
            nv=nv+1;  vir(nv)=(k-1)*nfuc+p
          else
            error stop 'CC residual: occupation must be 0 or 2'
          end if
        end do
      end do

      R1=cmplx(0.0_dp,0.0_dp,kind=dpc)
      R2=cmplx(0.0_dp,0.0_dp,kind=dpc)

      ! ---- R1: k-diagonal by construction ----
      do ki=1,N_cell
        do ii=1,nOcc(ki)
          i=(ki-1)*nfuc+kOcc(ii,ki)
          do aa=1,nVir(ki)
            a=(ki-1)*nfuc+kVir(aa,ki)
            v=cmplx(0.0_dp,0.0_dp,kind=dpc)
            do cc=1,nv
              c=vir(cc); v=v+f2(Fvv,a,c)*s1(i,c)
            end do
            do mm=1,no
              m=occ(mm); v=v-f2(Foo,m,i)*s1(m,a)
              do cc=1,nv
                c=vir(cc)
                v=v+f2(Fov,m,c)*(2.0_dp*s2(m,i,c,a)-s2(i,m,c,a)+s1(i,c)*s1(m,a)) &
                  +(2.0_dp*eri(m,c,a,i)-eri(m,i,a,c))*s1(m,c)
                do dd=1,nv
                  d=vir(dd)
                  v=v+(2.0_dp*eri(m,d,a,c)-eri(m,c,a,d))*(s2(i,m,c,d)+s1(m,d)*s1(i,c))
                end do
                do ll=1,no
                  l=occ(ll)
                  v=v+(-2.0_dp*eri(l,c,m,i)+eri(m,c,l,i))*(s2(m,l,a,c)+s1(l,c)*s1(m,a))
                end do
              end do
            end do
            R1(bandA(i),bandA(a),ki)=v
          end do
        end do
      end do

      ! ---- R2: kb fixed by momentum conservation ----
      do ki=1,N_cell
        do kj=1,N_cell
          do ka=1,N_cell
            kb=ktab(ki,kj,ka)
            do ii=1,nOcc(ki)
              i=(ki-1)*nfuc+kOcc(ii,ki)
              do jj=1,nOcc(kj)
                j=(kj-1)*nfuc+kOcc(jj,kj)
                do aa=1,nVir(ka)
                  a=(ka-1)*nfuc+kVir(aa,ka)
                  do bb=1,nVir(kb)
                    b=(kb-1)*nfuc+kVir(bb,kb)
                    v=conjg(eri(i,a,j,b))+pair(i,j,a,b)+pair(j,i,b,a) &
                      -(blk(ki)%energy(bandA(i))+blk(kj)%energy(bandA(j)) &
                       -blk(ka)%energy(bandA(a))-blk(kb)%energy(bandA(b)))*s2(i,j,a,b)
                    do km=1,N_cell
                      kl=ktab(ki,kj,km)
                      do mm=1,nOcc(km)
                        m=(km-1)*nfuc+kOcc(mm,km)
                        do ll=1,nOcc(kl)
                          l=(kl-1)*nfuc+kOcc(ll,kl)
                          v=v+w4(Woooo,m,l,i,j)*(s2(m,l,a,b)+s1(m,a)*s1(l,b))
                        end do
                      end do
                    end do
                    do kc=1,N_cell
                      kd=ktab(ka,kb,kc)
                      do cc=1,nVir(kc)
                        c=(kc-1)*nfuc+kVir(cc,kc)
                        do dd=1,nVir(kd)
                          d=(kd-1)*nfuc+kVir(dd,kd)
                          v=v+w4(Wvvvv,a,b,c,d)*(s2(i,j,c,d)+s1(i,c)*s1(j,d))
                        end do
                      end do
                    end do
                    R2(bandA(i),bandA(j),bandA(a),bandA(b),ki,kj,ka)=v
                  end do
                end do
              end do
            end do
          end do
        end do
      end do


      contains

      integer function mom(x)
        integer, intent(in) :: x
        mom=momA(x)
      end function

      integer function band(x)
        integer, intent(in) :: x
        band=bandA(x)
      end function

      function eri(p,q,r,s) result(v)
        integer, intent(in) :: p,q,r,s
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (momA(s)/=ktab(momA(p),momA(r),momA(q))) return
        v=g(bandA(p),bandA(q),bandA(r),bandA(s),momA(p),momA(q),momA(r))
      end function

      function s1(i,a) result(v)
        integer, intent(in) :: i,a
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (momA(i)==momA(a)) v=t1(bandA(i),bandA(a),momA(i))
      end function

      function s2(i,j,a,b) result(v)
        integer, intent(in) :: i,j,a,b
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (momA(b)/=ktab(momA(i),momA(j),momA(a))) return
        v=t2(bandA(i),bandA(j),bandA(a),bandA(b),momA(i),momA(j),momA(a))
      end function

      function f2(f,p,q) result(v)
        complex(dpc), intent(in) :: f(nfuc,nfuc,N_cell)
        integer, intent(in) :: p,q
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (momA(p)==momA(q)) v=f(bandA(p),bandA(q),momA(p))
      end function

      function w4(w,p,q,r,s) result(v)
        complex(dpc), intent(in) :: w(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
        integer, intent(in) :: p,q,r,s
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (momA(s)/=ktab(momA(p),momA(q),momA(r))) return
        v=w(bandA(p),bandA(q),bandA(r),bandA(s),momA(p),momA(q),momA(r))
      end function

      function pair_ref(i,j,a,b) result(v)
        integer, intent(in) :: i,j,a,b
        integer :: mm,ll,cc,m,l,c
        complex(dpc) :: v,z
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        do cc=1,nv
          c=vir(cc)
          z=conjg(eri(i,a,c,b))
          do mm=1,no
            m=occ(mm); z=z-eri(m,i,b,c)*s1(m,a)
          end do
          v=v+z*s1(j,c)+f2(Lvv,a,c)*s2(i,j,c,b)
          if (a==c) v=v-blk(momA(a))%energy(bandA(a))*s2(i,j,c,b)
        end do
        do mm=1,no
          m=occ(mm)
          z=conjg(eri(i,a,j,m))
          do cc=1,nv
            c=vir(cc); z=z+eri(m,c,a,i)*s1(j,c)
          end do
          v=v-z*s1(m,b)-f2(Loo,m,i)*s2(m,j,a,b)
          if (m==i) v=v+blk(momA(i))%energy(bandA(i))*s2(m,j,a,b)
          do cc=1,nv
            c=vir(cc)
            v=v+(2.0_dp*w4(Wvoov,a,m,i,c)-w4(Wvovo,a,m,c,i))*s2(m,j,c,b) &
              -w4(Wvoov,a,m,i,c)*s2(m,j,b,c) &
              -w4(Wvovo,b,m,c,i)*s2(m,j,a,c)
          end do
        end do
      end function pair_ref


        function pair(i,j,a,b) result(v)
        integer, intent(in) :: i,j,a,b
        integer :: mm,cc,kc,im,ic,m,c
        complex(dpc) :: v,z

        ! Unsymmetrized part; caller adds pair_fast(j,i,b,a).
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)

        ! ---- terms 1-6: reference loops, verbatim ----
        do cc=1,nv
          c=vir(cc)
          z=conjg(eri(i,a,c,b))
          do mm=1,no
            m=occ(mm); z=z-eri(m,i,b,c)*s1(m,a)
          end do
          v=v+z*s1(j,c)+f2(Lvv,a,c)*s2(i,j,c,b)
          if (a==c) v=v-blk(mom(a))%energy(band(a))*s2(i,j,c,b)
        end do

        do mm=1,no
          m=occ(mm)
          z=conjg(eri(i,a,j,m))
          do cc=1,nv
            c=vir(cc); z=z+eri(m,c,a,i)*s1(j,c)
          end do
          v=v-z*s1(m,b)-f2(Loo,m,i)*s2(m,j,a,b)
          if (m==i) v=v+blk(mom(i))%energy(band(i))*s2(m,j,a,b)

          ! ---- W dressing, for this m ----
          ! kc = km + kj - kb, with km = mom(m)
          kc=ktab(mom(m),mom(j),mom(b))
          do ic=1,nVir(kc)
            c=(kc-1)*nfuc+kVir(ic,kc)
            v=v+(2.0_dp*w4(Wvoov,a,m,i,c)-w4(Wvovo,a,m,c,i))*s2(m,j,c,b) &
                 -w4(Wvoov,a,m,i,c)*s2(m,j,b,c)
          end do

          ! kc = km + kj - ka, with km = mom(m)
          kc=ktab(mom(m),mom(j),mom(a))
          do ic=1,nVir(kc)
            c=(kc-1)*nfuc+kVir(ic,kc)
            v=v-w4(Wvovo,b,m,c,i)*s2(m,j,a,c)
          end do

        end do
      end function pair



end subroutine build_ccsd_residuals
