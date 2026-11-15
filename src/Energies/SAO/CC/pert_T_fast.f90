subroutine build_ccsd_t_fast(blk,g,ktab,t1,t2,et)

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
      real(dp)     , intent(out):: et
      ! local !
      integer, parameter :: sg(3,6) = reshape([1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1],[3,6])
      integer(8), parameter :: budget = 256_8*1024_8*1024_8
      logical  :: isO(nfuc,N_cell), isV(nfuc,N_cell)
      integer  :: no, nv, na, nchunk, ntask, it, nthr, k, p
      integer  :: ka, kb, kc, ki, kj, kk, Kv, wk, cmpk
      integer, allocatable :: task(:,:)
      complex(dpc) :: energy
      logical  :: outer

      isO=.false.; isV=.false.
      do k=1,N_cell
        do p=1,nfuc
          isO(p,k) = (blk(k)%occupation(p)==2)
          isV(p,k) = (blk(k)%occupation(p)==0)
          if (.not.(isO(p,k).or.isV(p,k))) error stop '(T): occupation must be 0 or 2'
        end do
      end do
      no = o_hi
      nv = nfuc - v_lo + 1
      et = 0.0_dp
      if (no<=0 .or. nv<=0) return

      ! a-chunk: W and H (6 orderings each) of size na*nv*nv*no^3
      na = int(max(1_8, min(int(nv,8), budget/(12_8*16_8*int(nv,8)**2*int(no,8)**3))))
      nchunk = (nv+na-1)/na

      ! ---- task list: (ka<=kb<=kc) x (ki<=kj<=kk, same total k) x a-chunk ----
      ntask = 0
      do it=1,2
        if (it==2) allocate(task(8,ntask))
        ntask = 0
        do ka=1,N_cell; do kb=ka,N_cell; do kc=kb,N_cell
          Kv = ktab(ktab(ka,kb,1),kc,1)
          do ki=1,N_cell; do kj=ki,N_cell
            kk = ktab(Kv,1,ktab(ki,kj,1))
            if (kk < kj) cycle
            ! [Nov 2026] k <-> -k symmetry: the task with all k negated gives the
            ! complex-conjugate energy -> keep the smaller one with weight 2
            wk = 1
            if (use_ksym) then
              cmpk = cmp_neg([ka,kb,kc,ki,kj,kk])
              if (cmpk > 0) cycle
              if (cmpk < 0) wk = 2
            end if
            do k=1,nchunk
              ntask = ntask+1
              if (it==2) task(:,ntask) = [ka,kb,kc,ki,kj,kk,k,wk]
            end do
          end do; end do
        end do; end do; end do
      end do

      nthr  = omp_get_max_threads()
      outer = (ntask >= 2*nthr)
      energy = (0.0_dp,0.0_dp)

      !$omp parallel do if(outer) schedule(dynamic,1) default(shared) private(it) reduction(+:energy)
      do it=1,ntask
        call t_task(task(:,it), energy)
      end do
      !$omp end parallel do

      energy = (2.0_dp/6.0_dp)*energy
      et = real(energy,dp)
      write(*,*)
      write(*,'(A)')          '*-----------------------------------*'
      write(*,'(A,f16.12,A)') '*-*   T-ccsd = ', et, '   *-*'
      write(*,'(A)')          '*-----------------------------------*'
      write(*,*)
      flush(6)
      deallocate(task)

contains

      !-----------------------------------------------------------------!
      ! one task: virtual k-triple (sorted), occupied multiset, a-chunk !
      !-----------------------------------------------------------------!
      subroutine t_task(tk, eacc)
        integer     , intent(in)    :: tk(8)
        complex(dpc), intent(inout) :: eacc
        integer :: kv3(3), ko(3,6), nord, io, s, q, a0, a1, nx
        integer :: o_jki, o_kij, o_kji, o_ikj, o_jik, i, j, kq, a, b, c
        integer :: wv, nkey, ikey, kkey(9), keys(9,36), use_key(6,6)
        real(dp) :: dijk, d
        complex(dpc) :: y, eloc
        complex(dpc), allocatable :: Wf(:,:,:,:,:,:,:), Hf(:,:,:,:,:,:,:)
        real(dp) :: eO(no,N_cell), eV(v_lo:nfuc,N_cell)
        integer :: kk1(3)

        kv3 = tk(1:3)
        ! weight of the virtual k-multiset (distinct orderings)
        if (kv3(1)==kv3(2) .and. kv3(2)==kv3(3)) then
          wv = 1
        else if (kv3(1)==kv3(2) .or. kv3(2)==kv3(3)) then
          wv = 3
        else
          wv = 6
        end if
        ! distinct orderings of the occupied multiset
        nord = 0
        do s=1,6
          kk1 = [tk(3+sg(1,s)), tk(3+sg(2,s)), tk(3+sg(3,s))]
          if (find_ord(ko,nord,kk1) == 0) then
            nord = nord+1; ko(:,nord) = kk1
          end if
        end do
        a0 = v_lo + (tk(7)-1)*na
        a1 = min(nfuc, a0+na-1)
        nx = a1-a0+1

        do kq=1,N_cell
          do i=1,no
            eO(i,kq) = 0.0_dp
            if (isO(i,kq)) eO(i,kq) = blk(kq)%energy(i)
          end do
          do a=v_lo,nfuc
            eV(a,kq) = 0.0_dp
            if (isV(a,kq)) eV(a,kq) = blk(kq)%energy(a)
          end do
        end do

        allocate(Wf(nx,nv,nv,no,no,no,nord), Hf(nx,nv,nv,no,no,no,nord))
        Wf = (0.0_dp,0.0_dp); Hf = (0.0_dp,0.0_dp)

        ! each distinct raw block (slot momenta + slot band ranges) is computed
        ! once and added to every (ordering, permutation) where it appears
        nkey = 0
        do io=1,nord
          do s=1,6
            kkey(1:3) = [kv3(sg(1,s)), kv3(sg(2,s)), kv3(sg(3,s))]
            kkey(4:6) = [ko(sg(1,s),io), ko(sg(2,s),io), ko(sg(3,s),io)]
            do q=1,3
              kkey(6+q) = 0                             ! one chunk: all slots full range
              if (nchunk > 1 .and. sg(q,s)==1) kkey(6+q) = 1   ! slot that holds the a-chunk
            end do
            ikey = 0
            do q=1,nkey
              if (all(keys(:,q)==kkey)) then
                ikey = q; exit
              end if
            end do
            if (ikey==0) then
              nkey = nkey+1; keys(:,nkey) = kkey; ikey = nkey
            end if
            use_key(s,io) = ikey
          end do
        end do
        do q=1,nkey
          call raw_block(keys(:,q), q, a0, a1, Wf, Hf, nord, use_key)
        end do
        Hf = Hf + Wf

        ! ---- energy ----
        eloc = (0.0_dp,0.0_dp)
        do io=1,nord
          o_jki = find_ord(ko,nord,[ko(2,io),ko(3,io),ko(1,io)])
          o_kij = find_ord(ko,nord,[ko(3,io),ko(1,io),ko(2,io)])
          o_kji = find_ord(ko,nord,[ko(3,io),ko(2,io),ko(1,io)])
          o_ikj = find_ord(ko,nord,[ko(1,io),ko(3,io),ko(2,io)])
          o_jik = find_ord(ko,nord,[ko(2,io),ko(1,io),ko(3,io)])
          do q=1,no;   if (.not.isO(q,ko(3,io))) cycle
          do j=1,no;   if (.not.isO(j,ko(2,io))) cycle
          do i=1,no;   if (.not.isO(i,ko(1,io))) cycle
            dijk = eO(i,ko(1,io)) + eO(j,ko(2,io)) + eO(q,ko(3,io))
            do c=1,nv;   if (.not.isV(v_lo-1+c,kv3(3))) cycle
            do b=1,nv;   if (.not.isV(v_lo-1+b,kv3(2))) cycle
            do a=1,nx;   if (.not.isV(a0-1+a,kv3(1))) cycle
              d = dijk - eV(a0-1+a,kv3(1)) - eV(v_lo-1+b,kv3(2)) - eV(v_lo-1+c,kv3(3))
              y =  4.0_dp*Hf(a,b,c,i,j,q,io) + Hf(a,b,c,j,q,i,o_jki) + Hf(a,b,c,q,i,j,o_kij) &
                 - 2.0_dp*Hf(a,b,c,q,j,i,o_kji) - 2.0_dp*Hf(a,b,c,i,q,j,o_ikj)              &
                 - 2.0_dp*Hf(a,b,c,j,i,q,o_jik)
              eloc = eloc + Wf(a,b,c,i,j,q,io)*conjg(y)/d
            end do; end do; end do
          end do; end do; end do
        end do
        eacc = eacc + real(wv*tk(8),dp)*eloc        ! [Nov 2026] tk(8) = k <-> -k weight
        deallocate(Wf,Hf)
      end subroutine t_task

      ! [Nov 2026] compare (sorted v-triple, sorted o-triple) with the same
      ! for all k negated: -1 smaller, 0 equal, +1 larger
      integer function cmp_neg(t6)
        integer, intent(in) :: t6(6)
        integer :: n6(6), r
        n6(1:3) = sort3([ktab(1,1,t6(1)), ktab(1,1,t6(2)), ktab(1,1,t6(3))])
        n6(4:6) = sort3([ktab(1,1,t6(4)), ktab(1,1,t6(5)), ktab(1,1,t6(6))])
        cmp_neg = 0
        do r=1,6
          if (t6(r) < n6(r)) then
            cmp_neg = -1; return
          else if (t6(r) > n6(r)) then
            cmp_neg = 1; return
          end if
        end do
      end function cmp_neg

      function sort3(x) result(y)
        integer, intent(in) :: x(3)
        integer :: y(3), tmp
        y = x
        if (y(1) > y(2)) then; tmp=y(1); y(1)=y(2); y(2)=tmp; end if
        if (y(2) > y(3)) then; tmp=y(2); y(2)=y(3); y(3)=tmp; end if
        if (y(1) > y(2)) then; tmp=y(1); y(1)=y(2); y(2)=tmp; end if
      end function sort3

      integer function find_ord(ko, nord, kx)
        integer, intent(in) :: ko(3,6), nord, kx(3)
        integer :: r
        find_ord = 0
        do r=1,nord
          if (all(ko(:,r)==kx)) then
            find_ord = r; return
          end if
        end do
      end function find_ord

      !-----------------------------------------------------------------!
      ! one raw block, key = (kx,ky,kz, kl,km,kn, chunk flags of slots)  !
      !   raw_xyz(l,m,n) = sum_f  g*(l,x,f,y) t2(n,m,z,f)                !
      !                  - sum_m' g*(l,x,m,m') t2(m',n,y,z)              !
      !   v_xyz(l,m,n)   = g*(l,x,m,y) t1(n,z)          (only kn == kz)  !
      ! added to W(:,io) (+raw) and H(:,io) (+v/2) for every            !
      ! (permutation s, ordering io) that uses this key: slot p goes to  !
      ! position s(p) of (a,b,c) and of (i,j,k)                          !
      !-----------------------------------------------------------------!
      subroutine raw_block(key, ikey, a0, a1, W, H, nord, use_key)
        integer, intent(in) :: key(9), ikey, a0, a1, nord, use_key(6,6)
        complex(dpc), intent(inout) :: W(a1-a0+1,nv,nv,no,no,no,nord), H(a1-a0+1,nv,nv,no,no,no,nord)
        integer :: kx,ky,kz,kl,km,kn,kf,kmp
        integer :: lo(3), hi(3), nn(3), x0,x1,y0,y1,z0,z1,nx,ny,nz
        integer :: x,y,z,m,f,mp,p
        integer :: ext(6), sc(6), nw(6)
        complex(dpc), allocatable :: Am1(:,:,:,:), Cm1(:,:,:,:,:,:), Am2(:,:,:,:), Bm2(:,:,:,:), Cm2(:,:,:,:,:,:)
        complex(dpc), allocatable :: Gv(:,:,:,:), Tv(:,:), Vb(:,:,:,:,:,:)
        complex(dpc), parameter :: one=(1.0_dp,0.0_dp), zero=(0.0_dp,0.0_dp)

        kx=key(1); ky=key(2); kz=key(3); kl=key(4); km=key(5); kn=key(6)
        do p=1,3
          if (key(6+p)==1) then
            lo(p)=a0; hi(p)=a1
          else
            lo(p)=v_lo; hi(p)=nfuc
          end if
          nn(p)=hi(p)-lo(p)+1
        end do
        x0=lo(1); x1=hi(1); y0=lo(2); y1=hi(2); z0=lo(3); z1=hi(3)
        nx=nn(1); ny=nn(2); nz=nn(3)
        nw = [a1-a0+1, nv, nv, no, no, no]

        ! ---------------- term 1 : Cm1(l,x,y,n,m,z) ----------------
        kf = ktab(kx,ky,kl)
        allocate(Am1(no,x0:x1,y0:y1,v_lo:nfuc), Cm1(no,x0:x1,y0:y1,no,no,z0:z1))
        do f=v_lo,nfuc
          if (isV(f,kf)) then
            do y=y0,y1; do x=x0,x1
              Am1(:,x,y,f) = conjg(gv1(1,no,x,f,y,kl,kx,kf))
            end do; end do
          else
            Am1(:,:,:,f) = zero
          end if
        end do
        call zgemm('N','T', no*nx*ny, no*no*nz, nv, one, Am1, no*nx*ny, &
                   t2(1,1,z0,v_lo,kn,km,kz), no*no*nv, zero, Cm1, no*nx*ny)
        deallocate(Am1)
        ! index order of Cm1: l x y n m z  -> slot of each: occ1 vir1 vir2 occ3 occ2 vir3
        ext = [no, nx, ny, no, no, nz]
        sc  = [1, no, no*nx, no*nx*ny, no*nx*ny*no, no*nx*ny*no*no]
        call scatter(Cm1, ext, sc, [4,1,2,6,5,3], one, W, ikey, use_key, nord, nw)
        deallocate(Cm1)

        ! ---------------- term 2 : Cm2(l,x,m,n,y,z) ----------------
        kmp = ktab(kl,km,kx)
        allocate(Am2(no,x0:x1,no,no), Bm2(no,no,y0:y1,z0:z1), Cm2(no,x0:x1,no,no,y0:y1,z0:z1))
        do mp=1,no
          if (isO(mp,kmp)) then
            do m=1,no; do x=x0,x1
              Am2(:,x,m,mp) = conjg(gv1(1,no,x,m,mp,kl,kx,km))
            end do; end do
          else
            Am2(:,:,:,mp) = zero
          end if
        end do
        do z=z0,z1; do y=y0,y1
          Bm2(:,:,y,z) = t2(:,:,y,z,kmp,kn,ky)
        end do; end do
        call zgemm('N','N', no*nx*no, no*ny*nz, no, one, Am2, no*nx*no, &
                   Bm2, no, zero, Cm2, no*nx*no)
        deallocate(Am2, Bm2)
        ext = [no, nx, no, no, ny, nz]
        sc  = [1, no, no*nx, no*nx*no, no*nx*no*no, no*nx*no*no*ny]
        call scatter(Cm2, ext, sc, [4,1,5,6,2,3], -one, W, ikey, use_key, nord, nw)
        deallocate(Cm2)

        ! ---------------- v term : Vb(l,x,m,y,n,z), only kn == kz ----------------
        if (kn==kz) then
          allocate(Gv(no,x0:x1,no,y0:y1), Tv(no,z0:z1), Vb(no,x0:x1,no,y0:y1,no,z0:z1))
          do y=y0,y1; do m=1,no; do x=x0,x1
            Gv(:,x,m,y) = conjg(gv1(1,no,x,m,y,kl,kx,km))
          end do; end do; end do
          do z=z0,z1
            do m=1,no
              Tv(m,z) = zero
              if (isO(m,kn) .and. isV(z,kz)) Tv(m,z) = t1(m,z,kz)
            end do
          end do
          call zgemm('N','T', no*nx*no*ny, no*nz, 1, one, Gv, no*nx*no*ny, Tv, no*nz, zero, Vb, no*nx*no*ny)
          ext = [no, nx, no, ny, no, nz]
          sc  = [1, no, no*nx, no*nx*no, no*nx*no*ny, no*nx*no*ny*no]
          call scatter(Vb, ext, sc, [4,1,5,2,6,3], (0.5_dp,0.0_dp), H, ikey, use_key, nord, nw)
          deallocate(Gv, Tv, Vb)
        end if

      end subroutine raw_block

      !-----------------------------------------------------------------!
      ! add alpha*C to W(:,io) for every (s,io) using block ikey.        !
      ! code(c) of C index c: 1..3 virtual slot, 4..6 occupied slot;    !
      ! slot p goes to position s(p) (virtual) / 3+s(p) (occupied)      !
      !-----------------------------------------------------------------!
      subroutine scatter(Cin, ext, sc, code, alpha, W, ikey, use_key, nord, nw)
        complex(dpc), intent(in)    :: Cin(*)
        integer     , intent(in)    :: ext(6), sc(6), code(6), ikey, nord, use_key(6,6), nw(6)
        complex(dpc), intent(in)    :: alpha
        complex(dpc), intent(inout) :: W(*)
        integer :: s, io, ic, q, e6(6), sc6(6), sw6(6), qc, ordq(6), nq
        integer(8) :: off
        sw6(1) = 1
        do q=2,6
          sw6(q) = sw6(q-1)*nw(q-1)
        end do
        do io=1,nord
          do s=1,6
            if (use_key(s,io) /= ikey) cycle
            do ic=1,6
              if (code(ic) <= 3) then
                q = sg(code(ic),s)
              else
                q = 3 + sg(code(ic)-3,s)
              end if
              e6(q)  = ext(ic)
              sc6(q) = sc(ic)
            end do
            off = int(io-1,8)*int(sw6(6),8)*int(nw(6),8)
            ! loop order: W-contiguous index innermost, C-contiguous index next
            ! (cache-friendly 2D transpose), then the other W positions
            qc = 1
            do q=1,6
              if (sc6(q)==1) qc = q
            end do
            ordq(1) = 1; nq = 1
            if (qc /= 1) then
              nq = 2; ordq(2) = qc
            end if
            do q=2,6
              if (q /= qc) then
                nq = nq+1; ordq(nq) = q
              end if
            end do
            call padd(W, off, Cin, e6(ordq), sw6(ordq), sc6(ordq), alpha)
          end do
        end do
      end subroutine scatter

      ! W(off + sum_q i_q*sw(q)) += alpha*C(sum_q i_q*sc(q)), loop 1 innermost
      subroutine padd(W, off, Cin, e, sw, sc, alpha)
        complex(dpc), intent(inout) :: W(*)
        complex(dpc), intent(in)    :: Cin(*), alpha
        integer(8)  , intent(in)    :: off
        integer     , intent(in)    :: e(6), sw(6), sc(6)
        integer :: i1,i2,i3,i4,i5,i6
        integer(8) :: w2,w3,w4,w5,w6
        integer :: c2,c3,c4,c5,c6
        do i6=0,e(6)-1
          w6 = off + int(i6,8)*sw(6); c6 = i6*sc(6)
          do i5=0,e(5)-1
            w5 = w6 + i5*sw(5); c5 = c6 + i5*sc(5)
            do i4=0,e(4)-1
              w4 = w5 + i4*sw(4); c4 = c5 + i4*sc(4)
              do i3=0,e(3)-1
                w3 = w4 + i3*sw(3); c3 = c4 + i3*sc(3)
                do i2=0,e(2)-1
                  w2 = w3 + i2*sw(2); c2 = c3 + i2*sc(2)
                  do i1=0,e(1)-1
                    W(w2+i1*sw(1)+1) = W(w2+i1*sw(1)+1) + alpha*Cin(c2+i1*sc(1)+1)
                  end do
                end do
              end do
            end do
          end do
        end do
      end subroutine padd
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

end subroutine build_ccsd_t_fast
