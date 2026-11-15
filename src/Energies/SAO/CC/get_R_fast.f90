subroutine build_ccsd_residuals_fast(blk,g,ktab,t1,t2,Foo,Fvv,Fov,Loo,Lvv, &
                                 Woooo,Wvoov,Wvovo,R1,R2)
      use constants_module, only: dp, dpc
      use cc_dims         , only: gidx, nslot      ! [Nov 2026] compact g
      use unitcell_module , only: nfuc, N_cell
      use cc_dims       , only: o_hi, v_lo, use_ksym   ! [Nov 2026] use_ksym
      use k_blocks        , only: K_block
      implicit none
      ! input !
      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(nfuc,nfuc,nfuc,nfuc,nslot)   ! [Nov 2026] compact g (cc_dims gidx)
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
      ! local !
      integer      :: occ(nfuc*N_cell), vir(nfuc*N_cell)
      integer      :: no,nv,p,k,i,j,a,b,c,d,l,m
      integer      :: ii,jj,aa,bb,cc,dd,ll,mm
      integer      :: oi,oj,va,vb,om,ol,vc,vd
      complex(dpc) :: v,tau
      integer      :: kOcc(nfuc,N_cell), nOcc(N_cell)
      integer      :: kVir(nfuc,N_cell), nVir(N_cell)
      integer      :: momA(nfuc*N_cell), bandA(nfuc*N_cell)
      integer      :: ki,kj,ka,kb,kc,kd,km,kl
      integer      :: nR1, ip
      integer      :: r1i(nfuc*nfuc*N_cell), r1a(nfuc*nfuc*N_cell)
      complex(dpc), parameter :: czero=(0.0_dp,0.0_dp)
      ! output !
      complex(dpc), intent(out) :: R1(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: R2(o_hi,o_hi,v_lo:nfuc,v_lo:nfuc,N_cell,N_cell,N_cell)

      ! ---- setup: lists, lookup tables ----
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

      ! ---- flat list of k-diagonal (i,a) pairs for R1 ----
      nR1=0
      do ki=1,N_cell
        do ii=1,nOcc(ki)
          do aa=1,nVir(ki)
            nR1=nR1+1
            r1i(nR1)=(ki-1)*nfuc+kOcc(ii,ki)
            r1a(nR1)=(ki-1)*nfuc+kVir(aa,ki)
          end do
        end do
      end do

      R1=czero   ! tiny, serial

      !$omp parallel default(shared) &
      !$omp private(ip,ki,kj,ka,kb,kc,kd,km,kl,ii,jj,aa,bb,cc,dd,ll,mm, &
      !$omp         i,j,a,b,c,d,l,m,oi,oj,va,vb,om,ol,vc,vd,v,tau)

      ! ---- R1: one (i,a) pair per task ----
      !$omp do schedule(dynamic)
      do ip=1,nR1
        i=r1i(ip); a=r1a(ip); ki=momA(i)
        v=czero
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
      !$omp end do nowait

      ! ---- R2: one (ki,kj,ka) slab per task, kb fixed ----
      !$omp do collapse(3) schedule(dynamic)
      do ki=1,N_cell
        do kj=1,N_cell
          do ka=1,N_cell
            ! [Nov 2026] k <-> -k symmetry: only canonical slabs, the others
            ! are filled below with R2(-ki,-kj,-ka) = R2(ki,kj,ka)^*
            if (use_ksym) then
              if (.not.canon3(ki,kj,ka)) cycle
            end if
            R2(:,:,:,:,ki,kj,ka)=czero
            kb=ktab(ki,kj,ka)
            do ii=1,nOcc(ki)
              oi=kOcc(ii,ki); i=(ki-1)*nfuc+oi
              do jj=1,nOcc(kj)
                oj=kOcc(jj,kj); j=(kj-1)*nfuc+oj
                do aa=1,nVir(ka)
                  va=kVir(aa,ka); a=(ka-1)*nfuc+va
                  do bb=1,nVir(kb)
                    vb=kVir(bb,kb); b=(kb-1)*nfuc+vb

                    v=conjg(gq(oi,va,oj,vb,ki,ka,kj)) &
                      +pair(i,j,a,b)+pair(j,i,b,a) &
                      -(blk(ki)%energy(oi)+blk(kj)%energy(oj) &
                       -blk(ka)%energy(va)-blk(kb)%energy(vb)) &
                       *t2(oi,oj,va,vb,ki,kj,ka)

                    ! hole-hole ladder: km+kl = ki+kj
                    do km=1,N_cell
                      kl=ktab(ki,kj,km)
                      do mm=1,nOcc(km)
                        om=kOcc(mm,km)
                        do ll=1,nOcc(kl)
                          ol=kOcc(ll,kl)
                          tau=t2(om,ol,va,vb,km,kl,ka)
                          if (km==ka) tau=tau+t1(om,va,km)*t1(ol,vb,kl)
                          v=v+Woooo(om,ol,oi,oj,km,kl,ki)*tau
                        end do
                      end do
                    end do

                    ! particle-particle ladder: done with zgemm after this loop

                    R2(oi,oj,va,vb,ki,kj,ka)=v
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      !$omp end do nowait
      !$omp end parallel

      ! [Nov 2026] k <-> -k symmetry: fill the non-canonical slabs
      if (use_ksym) then
        !$omp parallel do collapse(3) default(shared) private(ki,kj,ka)
        do ki=1,N_cell
          do kj=1,N_cell
            do ka=1,N_cell
              if (.not.canon3(ki,kj,ka)) R2(:,:,:,:,ki,kj,ka) = &
                 conjg(R2(:,:,:,:,ktab(1,1,ki),ktab(1,1,kj),ktab(1,1,ka)))
            end do
          end do
        end do
        !$omp end parallel do
      end if

      call ladder_pp_gemm(blk,g,ktab,t1,t2,R2)
      call ring_gemm(blk,ktab,t2,Wvoov,Wvovo,R2)

      contains
      ! [Nov 2026] (k1,k2,k3) <= (-k1,-k2,-k3) in lexicographic order
      logical function canon3(k1,k2,k3)
        integer, intent(in) :: k1,k2,k3
        integer :: m1,m2,m3
        m1=ktab(1,1,k1); m2=ktab(1,1,k2); m3=ktab(1,1,k3)
        canon3 = (k1<m1) .or. (k1==m1 .and. (k2<m2 .or. (k2==m2 .and. k3<=m3)))
      end function
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
        v=czero
        if (momA(s)/=ktab(momA(p),momA(r),momA(q))) return
        v=gq(bandA(p),bandA(q),bandA(r),bandA(s),momA(p),momA(q),momA(r))
      end function
      function s1(i,a) result(v)
        integer, intent(in) :: i,a
        complex(dpc) :: v
        v=czero
        if (momA(i)==momA(a)) v=t1(bandA(i),bandA(a),momA(i))
      end function
      function s2(i,j,a,b) result(v)
        integer, intent(in) :: i,j,a,b
        complex(dpc) :: v
        v=czero
        if (momA(b)/=ktab(momA(i),momA(j),momA(a))) return
        v=t2(bandA(i),bandA(j),bandA(a),bandA(b),momA(i),momA(j),momA(a))
      end function
      function f2(f,p,q) result(v)
        complex(dpc), intent(in) :: f(nfuc,nfuc,N_cell)
        integer, intent(in) :: p,q
        complex(dpc) :: v
        v=czero
        if (momA(p)==momA(q)) v=f(bandA(p),bandA(q),momA(p))
      end function
      function w4_voov(w,p,q,r,s) result(v)
        complex(dpc), intent(in) :: w(v_lo:nfuc,o_hi,o_hi,v_lo:nfuc,N_cell,N_cell,N_cell)
        integer, intent(in) :: p,q,r,s
        complex(dpc) :: v
        v=czero
        if (momA(s)/=ktab(momA(p),momA(q),momA(r))) return
        v=w(bandA(p),bandA(q),bandA(r),bandA(s),momA(p),momA(q),momA(r))
      end function
      function w4_vovo(w,p,q,r,s) result(v)
        complex(dpc), intent(in) :: w(v_lo:nfuc,o_hi,v_lo:nfuc,o_hi,N_cell,N_cell,N_cell)
        integer, intent(in) :: p,q,r,s
        complex(dpc) :: v
        v=czero
        if (momA(s)/=ktab(momA(p),momA(q),momA(r))) return
        v=w(bandA(p),bandA(q),bandA(r),bandA(s),momA(p),momA(q),momA(r))
      end function
      function pair(i,j,a,b) result(v)
        ! momentum-restricted: only the k blocks where the terms are non-zero
        ! (the skipped terms are exact zeros; same summation order)
        integer, intent(in) :: i,j,a,b
        integer :: mm,cc,kc,ic,m,c,kj_,ka_,kb_,kblk(2),nb,ib,km_
        complex(dpc) :: v,z
        v=czero
        kj_=mom(j); ka_=mom(a); kb_=mom(b)
        ! c runs over the virtuals of k_a (Lvv, energy term) and k_j (s1(j,c)), ascending k
        nb=1; kblk(1)=min(ka_,kj_)
        if (ka_/=kj_) then; nb=2; kblk(2)=max(ka_,kj_); end if
        do ib=1,nb
          kc=kblk(ib)
          do cc=1,nVir(kc)
            c=(kc-1)*nfuc+kVir(cc,kc)
            z=czero
            if (kc==kj_) then
              z=conjg(eri(i,a,c,b))
              do mm=1,nOcc(ka_)                     ! s1(m,a) needs k_m = k_a
                m=(ka_-1)*nfuc+kOcc(mm,ka_); z=z-eri(m,i,b,c)*s1(m,a)
              end do
            end if
            v=v+z*s1(j,c)+f2(Lvv,a,c)*s2(i,j,c,b)
            if (a==c) v=v-blk(mom(a))%energy(band(a))*s2(i,j,c,b)
          end do
        end do
        do mm=1,no
          m=occ(mm); km_=mom(m)
          z=czero
          if (km_==kb_) then                        ! s1(m,b) needs k_m = k_b
            z=conjg(eri(i,a,j,m))
            do cc=1,nVir(kj_)                       ! s1(j,c) needs k_c = k_j
              c=(kj_-1)*nfuc+kVir(cc,kj_); z=z+eri(m,c,a,i)*s1(j,c)
            end do
          end if
          v=v-z*s1(m,b)-f2(Loo,m,i)*s2(m,j,a,b)
          if (m==i) v=v+blk(mom(i))%energy(band(i))*s2(m,j,a,b)
          ! ring (W t2) terms: done with zgemm (ring_gemm)
        end do
      end function pair
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
end subroutine build_ccsd_residuals_fast
