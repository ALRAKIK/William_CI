subroutine build_ccsd_residuals(blk,g,t1,t2,Foo,Fvv,Fov,Loo,Lvv, &
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

      ! local ! 

      integer      :: occ(nfuc*N_cell)
      integer      :: vir(nfuc*N_cell)
      integer      :: no,nv,p,k,i,j,a,b,c,d,l,m
      integer      :: ii,jj,aa,bb,cc,dd,ll,mm
      integer      :: fourth_k
      complex(dpc) :: v

      ! output ! 
      
      complex(dpc), intent(out) :: R1(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: R2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    
  
      ! Canonical HF only: bare f_ov=0; no level shift.
      ! All supplied F/L blocks MUST retain their orbital-energy diagonals.
      ! Residual R=N-D*T; update later via T_new=T+R/D.

      no=0; nv=0

      do k=1,N_cell; do p=1,nfuc
        if (blk(k)%occupation(p)==2) then
          no=no+1; occ(no)=(k-1)*nfuc+p
        else if (blk(k)%occupation(p)==0) then
          nv=nv+1; vir(nv)=(k-1)*nfuc+p
        else
          error stop 'CC residual: occupation must be 0 or 2'
        end if
      end do; end do

      R1=cmplx(0.0_dp,0.0_dp,kind=dpc)
      R2=cmplx(0.0_dp,0.0_dp,kind=dpc)

      do ii=1,no 
        do aa=1,nv
          i=occ(ii); a=vir(aa)
          if (mom(i)/=mom(a)) cycle
            v=cmplx(0.0_dp,0.0_dp,kind=dpc)
            do cc=1,nv
              c=vir(cc); v=v+f2(Fvv,a,c)*s1(i,c)
            end do
            do mm=1,no
              m=occ(mm); v=v-f2(Foo,m,i)*s1(m,a)
              do cc=1,nv
                c=vir(cc)
                v=v+f2(Fov,m,c)*(2.0_dp*s2(m,i,c,a)-s2(i,m,c,a)+s1(i,c)*s1(m,a))+(2.0_dp*eri(m,c,a,i)-eri(m,i,a,c))*s1(m,c)
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
          R1(band(i),band(a),mom(i))=v
        end do
      end do


      do ii=1,no
        do jj=1,no
          do aa=1,nv
            do bb=1,nv
              i=occ(ii); j=occ(jj); a=vir(aa); b=vir(bb)
              if (mom(b)/=fourth_k(mom(i),mom(j),mom(a))) cycle
              v=conjg(eri(i,a,j,b))+pair(i,j,a,b)+pair(j,i,b,a)-(blk(mom(i))%energy(band(i))+blk(mom(j))%energy(band(j))-blk(mom(a))%energy(band(a))-blk(mom(b))%energy(band(b)))*s2(i,j,a,b)
              do mm=1,no
                do ll=1,no
                  m=occ(mm)
                  l=occ(ll)
                  v=v+w4(Woooo,m,l,i,j)*(s2(m,l,a,b)+s1(m,a)*s1(l,b))
                end do
              end do
              do cc=1,nv
                do dd=1,nv
                  c=vir(cc)
                  d=vir(dd)
                  v=v+w4(Wvvvv,a,b,c,d)*(s2(i,j,c,d)+s1(i,c)*s1(j,d))
                end do
              end do
              R2(band(i),band(j),band(a),band(b),mom(i),mom(j),mom(a))=v
            end do
          end do
        end do
      end do

      contains

      integer function mom(x)
        integer, intent(in) :: x
        mom=(x-1)/nfuc+1
      end function

      integer function band(x)
        integer, intent(in) :: x
        band=modulo(x-1,nfuc)+1
      end function

      function eri(p,q,r,s) result(v)
        integer, intent(in) :: p,q,r,s
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (mom(s)/=fourth_k(mom(p),mom(r),mom(q))) return
        v=g(band(p),band(q),band(r),band(s),mom(p),mom(q),mom(r))
      end function

      function s1(i,a) result(v)
        integer, intent(in) :: i,a
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (mom(i)==mom(a)) v=t1(band(i),band(a),mom(i))
      end function

      function s2(i,j,a,b) result(v)
        integer, intent(in) :: i,j,a,b
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (mom(b)/=fourth_k(mom(i),mom(j),mom(a))) return
        v=t2(band(i),band(j),band(a),band(b),mom(i),mom(j),mom(a))
      end function

      function f2(f,p,q) result(v)
        complex(dpc), intent(in) :: f(nfuc,nfuc,N_cell)
        integer, intent(in) :: p,q
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (mom(p)==mom(q)) v=f(band(p),band(q),mom(p))
      end function

      function w4(w,p,q,r,s) result(v)
        complex(dpc), intent(in) :: w(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
        integer, intent(in) :: p,q,r,s
        complex(dpc) :: v
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (mom(s)/=fourth_k(mom(p),mom(q),mom(r))) return
        v=w(band(p),band(q),band(r),band(s),mom(p),mom(q),mom(r))
      end function

      function pair(i,j,a,b) result(v)
        integer, intent(in) :: i,j,a,b
        integer :: mm,ll,cc,m,l,c
        complex(dpc) :: v,z
        ! Unsymmetrized part; caller adds simultaneous i<->j,a<->b.
        v=cmplx(0.0_dp,0.0_dp,kind=dpc)
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
          do cc=1,nv
            c=vir(cc)
            v=v+(2.0_dp*w4(Wvoov,a,m,i,c)-w4(Wvovo,a,m,c,i))*s2(m,j,c,b) &
              -w4(Wvoov,a,m,i,c)*s2(m,j,b,c) &
              -w4(Wvovo,b,m,c,i)*s2(m,j,a,c)
          end do
        end do
      end function pair

end subroutine build_ccsd_residuals