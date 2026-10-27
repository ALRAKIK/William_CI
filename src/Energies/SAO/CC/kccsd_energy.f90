
! This is the verified MP2 starting point, NOT a CCSD solver.
! Spatial restricted amplitudes: t1(i,a,ki), t2(i,j,a,b,ki,kj,ka).
! Bare chemist-order ERIs. Do not antisymmetrize the driving term.
! All energies below are total supercell correlation energies.

subroutine initialize_ccsd_baseline(blk,g,t1,t2,emp2,eimag,mp2_defect)

      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, N_cell, nx, ny, nz
      use k_blocks, only: K_block
      implicit none 


      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: t1(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      real(dp), intent(out) :: emp2,eimag,mp2_defect
      integer :: ki,kj,ka,kb,i,j,a,b,k,nocc
      integer :: fourth_k
      integer :: cc_integral_k4
      real(dp) :: gap,eden,direct_mp2
      complex(dpc) :: gd,gx,energy,res

      if (min(nx,ny,nz) < 1) error stop 'CC baseline: invalid grid'
      if (N_cell /= nx*ny*nz) error stop 'CC baseline: grid size mismatch'
      nocc = 0
      do k=1,N_cell
        if (blk(k)%dim /= nfuc) error stop 'CC baseline: block dimension mismatch'
        if (.not. allocated(blk(k)%occupation)) error stop 'CC baseline: missing occupations'
        if (.not. allocated(blk(k)%energy)) error stop 'CC baseline: missing energies'
        if (size(blk(k)%occupation) /= nfuc) error stop 'CC baseline: occupation size mismatch'
        if (size(blk(k)%energy) /= nfuc) error stop 'CC baseline: energy size mismatch'
        if (any(blk(k)%occupation /= 0 .and. blk(k)%occupation /= 2)) &
          error stop 'CC baseline: occupations must be 0 or 2'
        nocc = nocc + count(blk(k)%occupation == 2)
      end do

      t1 = cmplx(0.0_dp,0.0_dp,kind=dpc)
      t2 = cmplx(0.0_dp,0.0_dp,kind=dpc)
      energy = cmplx(0.0_dp,0.0_dp,kind=dpc)
      direct_mp2 = 0.0_dp
      mp2_defect = 0.0_dp

    do ka=1,N_cell
      do kj=1,N_cell
        do ki=1,N_cell
          kb = fourth_k(ki,kj,ka)
          if (kb /= cc_integral_k4(ki,ka,kj)) &
            error stop 'CC baseline: inconsistent momentum mapping'
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
                  gd = g(i,a,j,b,ki,ka,kj)
                  gx = g(i,b,j,a,ki,kb,kj)
                  t2(i,j,a,b,ki,kj,ka) = -conjg(gd)/gap
                  energy = energy + (2.0_dp*gd-gx)*t2(i,j,a,b,ki,kj,ka)
                  eden = real(2.0_dp*gd*conjg(gd)-gd*conjg(gx),dp)
                  direct_mp2 = direct_mp2 - eden/gap
                  ! No-coupling residual must vanish at the MP2 guess.
                  res = conjg(gd) + gap*t2(i,j,a,b,ki,kj,ka)
                  mp2_defect = max(mp2_defect,abs(res))
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    emp2 = real(energy,dp)
    eimag = aimag(energy)
    write(*,'(A,I6)') 'CC baseline occupied spatial orbitals: ', nocc
    write(*,'(A,F24.16)') 'CC baseline MP2 energy:              ', emp2
    write(*,'(A,ES14.6)') 'Amplitude/direct MP2 energy defect:  ', abs(emp2-direct_mp2)
    write(*,'(A,ES14.6)') 'No-coupling MP2 residual:            ', mp2_defect
    write(*,'(A,ES14.6)') 'MP2 imaginary energy:               ', eimag
end subroutine initialize_ccsd_baseline


  
  subroutine build_ccsd_f(blk,g,t1,t2,Foo,Fvv,Fov)
      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, N_cell
      use k_blocks, only: K_block

      implicit none


    type(K_block), intent(in) :: blk(N_cell)
    complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(in) :: t1(nfuc,nfuc,N_cell)
    complex(dpc), intent(in) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(out) :: Foo(nfuc,nfuc,N_cell)
    complex(dpc), intent(out) :: Fvv(nfuc,nfuc,N_cell)
    complex(dpc), intent(out) :: Fov(nfuc,nfuc,N_cell)
    integer :: k,km,kn,ke,kf,m,n,i,a,e,f
    integer :: fourth_k
    complex(dpc) :: v,tau
    Foo=cmplx(0.0_dp,0.0_dp,kind=dpc)
    Fvv=cmplx(0.0_dp,0.0_dp,kind=dpc)
    Fov=cmplx(0.0_dp,0.0_dp,kind=dpc)

    ! Foo(m,i,k) = f_mi + sum_nef (2(me|nf)-(mf|ne))*tau_in^ef.
    do k=1,N_cell
      do i=1,nfuc
        if (blk(k)%occupation(i)/=2) cycle
        do m=1,nfuc
          if (blk(k)%occupation(m)/=2) cycle
          if (m==i) Foo(m,i,k)=cmplx(blk(k)%energy(i),0.0_dp,kind=dpc)
          do kn=1,N_cell
            do ke=1,N_cell
              kf=fourth_k(k,kn,ke)
              do n=1,nfuc
                if (blk(kn)%occupation(n)/=2) cycle
                do e=1,nfuc
                  if (blk(ke)%occupation(e)/=0) cycle
                  do f=1,nfuc
                    if (blk(kf)%occupation(f)/=0) cycle
                    tau=t2(i,n,e,f,k,kn,ke)
                    if (ke==k .and. kf==kn) tau=tau+t1(i,e,k)*t1(n,f,kn)
                    v=2.0_dp*g(m,e,n,f,k,ke,kn)-g(m,f,n,e,k,kf,kn)
                    Foo(m,i,k)=Foo(m,i,k)+v*tau
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    ! Fvv(a,e,k) = f_ae - sum_mnf (2(me|nf)-(mf|ne))*tau_mn^af.
    do k=1,N_cell
      do e=1,nfuc
        if (blk(k)%occupation(e)/=0) cycle
        do a=1,nfuc
          if (blk(k)%occupation(a)/=0) cycle
          if (a==e) Fvv(a,e,k)=cmplx(blk(k)%energy(a),0.0_dp,kind=dpc)
          do km=1,N_cell
            do kn=1,N_cell
              kf=fourth_k(km,kn,k)
              do m=1,nfuc
                if (blk(km)%occupation(m)/=2) cycle
                do n=1,nfuc
                  if (blk(kn)%occupation(n)/=2) cycle
                  do f=1,nfuc
                    if (blk(kf)%occupation(f)/=0) cycle
                    tau=t2(m,n,a,f,km,kn,k)
                    if (km==k .and. kf==kn) tau=tau+t1(m,a,km)*t1(n,f,kn)
                    v=2.0_dp*g(m,e,n,f,km,k,kn)-g(m,f,n,e,km,kf,kn)
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
      do e=1,nfuc
        if (blk(k)%occupation(e)/=0) cycle
        do m=1,nfuc
          if (blk(k)%occupation(m)/=2) cycle
          do kn=1,N_cell
            do n=1,nfuc
              if (blk(kn)%occupation(n)/=2) cycle
              do f=1,nfuc
                if (blk(kn)%occupation(f)/=0) cycle
                v=2.0_dp*g(m,e,n,f,k,k,kn)-g(m,f,n,e,k,kn,kn)
                Fov(m,e,k)=Fov(m,e,k)+v*t1(n,f,kn)
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine build_ccsd_f

  subroutine test_ccsd_f_zero(blk,g)


      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, N_cell
      use k_blocks, only: K_block

      implicit none


    type(K_block), intent(in) :: blk(N_cell)
    complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), allocatable :: z1(:,:,:),z2(:,:,:,:,:,:,:)
    complex(dpc) :: fo(nfuc,nfuc,N_cell),fv(nfuc,nfuc,N_cell), &
                    fx(nfuc,nfuc,N_cell)
    integer :: k,p
    allocate(z1(nfuc,nfuc,N_cell))
    allocate(z2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell))
    z1=cmplx(0.0_dp,0.0_dp,kind=dpc)
    z2=cmplx(0.0_dp,0.0_dp,kind=dpc)
    call build_ccsd_f(blk,g,z1,z2,fo,fv,fx)
    do k=1,N_cell
      do p=1,nfuc
        if (blk(k)%occupation(p)==2) fo(p,p,k)=fo(p,p,k)-blk(k)%energy(p)
        if (blk(k)%occupation(p)==0) fv(p,p,k)=fv(p,p,k)-blk(k)%energy(p)
      end do
    end do
    write(*,'(A,ES14.6)') 'F zero-limit Foo defect: ',maxval(abs(fo))
    write(*,'(A,ES14.6)') 'F zero-limit Fvv defect: ',maxval(abs(fv))
    write(*,'(A,ES14.6)') 'F zero-limit Fov defect: ',maxval(abs(fx))
    deallocate(z1,z2)
  end subroutine test_ccsd_f_zero


  subroutine build_ccsd_w(blk,g,t1,t2,Woooo,Wvvvv,Wvoov,Wvovo)

      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, N_cell
      use k_blocks, only: K_block
      implicit none 

    type(K_block), intent(in) :: blk(N_cell)
    complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(in) :: t1(nfuc,nfuc,N_cell)
    complex(dpc), intent(in) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(out) :: Woooo(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(out) :: Wvvvv(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(out) :: Wvoov(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(out) :: Wvovo(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    integer :: occ(nfuc*N_cell), vir(nfuc*N_cell), no,nv,x,y,z,h
    integer :: p,q,r,s,kp,kq,kr,ks,bp,bq,br,bs,m,e,l,d
    integer, external :: fourth_k
    complex(dpc) :: w, tau
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
    function eri(a,b,c,d) result(v)
      integer, intent(in) :: a,b,c,d
      complex(dpc) :: v
      v=cmplx(0.0_dp,0.0_dp,kind=dpc)
      ! Chemist momentum: kd=ka-kb+kc.
      if (mom(d)/=fourth_k(mom(a),mom(c),mom(b))) return
      v=g(band(a),band(b),band(c),band(d),mom(a),mom(b),mom(c))
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

  subroutine test_ccsd_w_zero(blk,g)

      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, N_cell
      use k_blocks, only: K_block

    type(K_block), intent(in) :: blk(N_cell)
    complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), allocatable :: z1(:,:,:),z2(:,:,:,:,:,:,:)
    complex(dpc), allocatable :: wo(:,:,:,:,:,:,:),wv(:,:,:,:,:,:,:)
    complex(dpc), allocatable :: wx(:,:,:,:,:,:,:),wy(:,:,:,:,:,:,:)
    complex(dpc) :: ref
    real(dp) :: defect(4)
    integer :: kp,kq,kr,ks,p,q,r,s,op,oq,orr,os
    integer, external :: fourth_k
    allocate(z1(nfuc,nfuc,N_cell))
    allocate(z2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell))
    allocate(wo, mold=z2); allocate(wv, mold=z2)
    allocate(wx, mold=z2); allocate(wy, mold=z2)
    z1=cmplx(0.0_dp,0.0_dp,kind=dpc)
    z2=cmplx(0.0_dp,0.0_dp,kind=dpc)
    call build_ccsd_w(blk,g,z1,z2,wo,wv,wx,wy)
    defect=0.0_dp
    do kp=1,N_cell
      do kq=1,N_cell
        do kr=1,N_cell
          ks=fourth_k(kp,kq,kr)
          do p=1,nfuc
            op=blk(kp)%occupation(p)
            do q=1,nfuc
              oq=blk(kq)%occupation(q)
              do r=1,nfuc
                orr=blk(kr)%occupation(r)
                do s=1,nfuc
                  os=blk(ks)%occupation(s)
                  ref=cmplx(0.0_dp,0.0_dp,kind=dpc)
                  if (op==2.and.oq==2.and.orr==2.and.os==2) ref=g(p,r,q,s,kp,kr,kq)
                  defect(1)=max(defect(1),abs(wo(p,q,r,s,kp,kq,kr)-ref))
                  ref=cmplx(0.0_dp,0.0_dp,kind=dpc)
                  if (op==0.and.oq==0.and.orr==0.and.os==0) ref=g(p,r,q,s,kp,kr,kq)
                  defect(2)=max(defect(2),abs(wv(p,q,r,s,kp,kq,kr)-ref))
                  ref=cmplx(0.0_dp,0.0_dp,kind=dpc)
                  if (op==0.and.oq==2.and.orr==2.and.os==0) ref=g(q,s,p,r,kq,ks,kp)
                  defect(3)=max(defect(3),abs(wx(p,q,r,s,kp,kq,kr)-ref))
                  ref=cmplx(0.0_dp,0.0_dp,kind=dpc)
                  if (op==0.and.oq==2.and.orr==0.and.os==2) ref=g(q,s,p,r,kq,ks,kp)
                  defect(4)=max(defect(4),abs(wy(p,q,r,s,kp,kq,kr)-ref))
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    write(*,'(A,ES14.6)') 'W zero-limit Woooo defect: ',defect(1)
    write(*,'(A,ES14.6)') 'W zero-limit Wvvvv defect: ',defect(2)
    write(*,'(A,ES14.6)') 'W zero-limit Wvoov defect: ',defect(3)
    write(*,'(A,ES14.6)') 'W zero-limit Wvovo defect: ',defect(4)
    deallocate(z1,z2,wo,wv,wx,wy)
  end subroutine test_ccsd_w_zero



  subroutine test_ccsd_w_dense(blk,g,t1,t2)

      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, N_cell
      use k_blocks, only: K_block

    type(K_block), intent(in) :: blk(N_cell)
    complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(in) :: t1(nfuc,nfuc,N_cell)
    complex(dpc), intent(in) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), allocatable :: s1(:,:,:),s2(:,:,:,:,:,:,:)
    complex(dpc), allocatable :: wo(:,:,:,:,:,:,:),wv(:,:,:,:,:,:,:)
    complex(dpc), allocatable :: wx(:,:,:,:,:,:,:),wy(:,:,:,:,:,:,:)
    complex(dpc), allocatable :: v(:,:,:,:),u(:,:),t(:,:,:,:)
    integer :: occ(nfuc*N_cell),vir(nfuc*N_cell),km(nfuc*N_cell),bn(nfuc*N_cell)
    integer :: ng,no,nv,p,q,r,s,i,j,a,b,k1,k2,k3,k4,mode,l,d,x,y
    integer, external :: fourth_k
    real(dp) :: err(4),scale(4),tol,phase
    complex(dpc) :: ref(4),got(4)
    character(len=6) :: names(4)
    names=[character(len=6) :: 'Woooo','Wvvvv','Wvoov','Wvovo']
    ng=nfuc*N_cell; no=0; nv=0
    do p=1,ng
      km(p)=(p-1)/nfuc+1; bn(p)=modulo(p-1,nfuc)+1
      if (blk(km(p))%occupation(bn(p))==2) then
        no=no+1; occ(no)=p
      else if (blk(km(p))%occupation(bn(p))==0) then
        nv=nv+1; vir(nv)=p
      else
        error stop 'W validation: occupation must be 0 or 2'
      end if
    end do
    ! For nfuc*N_cell=60, the three dense rank-four arrays below
    ! would be expensive. Only two are allocated: v and t (~415 MB
    ! total with 16-byte complex), in addition to compressed arrays.
    allocate(v(ng,ng,ng,ng),t(ng,ng,ng,ng),u(ng,ng))
    allocate(s1,source=t1); allocate(s2,source=t2)
    allocate(wo,mold=t2); allocate(wv,mold=t2)
    allocate(wx,mold=t2); allocate(wy,mold=t2)
    v=cmplx(0.0_dp,0.0_dp,kind=dpc)
    ! Expand chemist ERIs directly; no builder accessors are used.
    do p=1,ng; do q=1,ng; do r=1,ng; do s=1,ng
      if (km(s)/=fourth_k(km(p),km(r),km(q))) cycle
      v(p,q,r,s)=g(bn(p),bn(q),bn(r),bn(s),km(p),km(q),km(r))
    end do; end do; end do; end do
    do mode=1,2
      if (mode==1) then
        write(*,'(A)') 'W dense validation: supplied amplitudes'
      else
        write(*,'(A)') 'W dense validation: synthetic complex T1/T2'
        s1=cmplx(0.0_dp,0.0_dp,kind=dpc)
        s2=cmplx(0.0_dp,0.0_dp,kind=dpc)
        do x=1,no; do y=1,nv
          i=occ(x); a=vir(y)
          if (km(i)/=km(a)) cycle
          phase=real(3*i+5*a,dp)
          s1(bn(i),bn(a),km(i))=1.0e-3_dp*cmplx(sin(phase),cos(phase),kind=dpc)
        end do; end do
        do x=1,no; do y=1,no
          i=occ(x); j=occ(y)
          do p=1,nv; do q=1,nv
            a=vir(p); b=vir(q)
            if (km(b)/=fourth_k(km(i),km(j),km(a))) cycle
            ! Symmetric under simultaneous (i,j)<->(j,i),(a,b)<->(b,a).
            phase=real(i+j+3*(a+b)+i*a+j*b,dp)
            s2(bn(i),bn(j),bn(a),bn(b),km(i),km(j),km(a))= &
              1.0e-3_dp*cmplx(sin(phase),cos(phase),kind=dpc)
          end do; end do
        end do; end do
      end if
      u=cmplx(0.0_dp,0.0_dp,kind=dpc)
      t=cmplx(0.0_dp,0.0_dp,kind=dpc)
      do x=1,no; do y=1,nv
        i=occ(x); a=vir(y)
        if (km(i)==km(a)) u(i,a)=s1(bn(i),bn(a),km(i))
      end do; end do
      do x=1,no; do y=1,no
        i=occ(x); j=occ(y)
        do p=1,nv; do q=1,nv
          a=vir(p); b=vir(q)
          if (km(b)/=fourth_k(km(i),km(j),km(a))) cycle
          t(i,j,a,b)=s2(bn(i),bn(j),bn(a),bn(b),km(i),km(j),km(a))
        end do; end do
      end do; end do
      call build_ccsd_w(blk,g,s1,s2,wo,wv,wx,wy)
      err=0.0_dp; scale=0.0_dp
      ! Dense reference: loop over ALL external global indices.
      ! Thus forbidden fourth momenta are also tested against zero.
      do p=1,ng; do q=1,ng; do r=1,ng; do s=1,ng
        ref=cmplx(0.0_dp,0.0_dp,kind=dpc)
        got=cmplx(0.0_dp,0.0_dp,kind=dpc)
        k1=km(p); k2=km(q); k3=km(r); k4=km(s)
        if (blk(k1)%occupation(bn(p))==2.and. &
            blk(k2)%occupation(bn(q))==2.and. &
            blk(k3)%occupation(bn(r))==2.and. &
            blk(k4)%occupation(bn(s))==2) then
          ref(1)=v(p,r,q,s)
          do a=1,nv
            d=vir(a)
            ref(1)=ref(1)+v(q,d,p,r)*u(s,d)+v(p,d,q,s)*u(r,d)
            do b=1,nv
              l=vir(b)
              ref(1)=ref(1)+v(p,d,q,l)*(t(r,s,d,l)+u(r,d)*u(s,l))
            end do
          end do
        end if
        if (blk(k1)%occupation(bn(p))==0.and. &
            blk(k2)%occupation(bn(q))==0.and. &
            blk(k3)%occupation(bn(r))==0.and. &
            blk(k4)%occupation(bn(s))==0) then
          ref(2)=v(p,r,q,s)
          do i=1,no
            l=occ(i)
            ref(2)=ref(2)-v(l,s,p,r)*u(l,q)-v(l,r,q,s)*u(l,p)
          end do
        end if
        if (blk(k1)%occupation(bn(p))==0.and.blk(k2)%occupation(bn(q))==2) then
          if (blk(k3)%occupation(bn(r))==2.and.blk(k4)%occupation(bn(s))==0) then
            ref(3)=v(q,s,p,r)
            do a=1,nv
              d=vir(a); ref(3)=ref(3)+v(q,s,p,d)*u(r,d)
            end do
            do i=1,no
              l=occ(i); ref(3)=ref(3)-v(q,s,l,r)*u(l,p)
              do a=1,nv
                d=vir(a)
                ref(3)=ref(3)+v(l,d,q,s)*(t(r,l,p,d)-0.5_dp*t(r,l,d,p)-u(r,d)*u(l,p)) &
                  -0.5_dp*v(l,s,q,d)*t(r,l,p,d)
              end do
            end do
          end if
          if (blk(k3)%occupation(bn(r))==0.and.blk(k4)%occupation(bn(s))==2) then
            ref(4)=v(q,s,p,r)
            do a=1,nv
              d=vir(a); ref(4)=ref(4)+v(q,d,p,r)*u(s,d)
            end do
            do i=1,no
              l=occ(i); ref(4)=ref(4)-v(l,r,q,s)*u(l,p)
              do a=1,nv
                d=vir(a)
                ref(4)=ref(4)-v(l,r,q,d)*(0.5_dp*t(s,l,d,p)+u(s,d)*u(l,p))
              end do
            end do
          end if
        end if
        if (k4==fourth_k(k1,k2,k3)) then
          got(1)=wo(bn(p),bn(q),bn(r),bn(s),k1,k2,k3)
          got(2)=wv(bn(p),bn(q),bn(r),bn(s),k1,k2,k3)
          got(3)=wx(bn(p),bn(q),bn(r),bn(s),k1,k2,k3)
          got(4)=wy(bn(p),bn(q),bn(r),bn(s),k1,k2,k3)
        end if
        err=max(err,abs(got-ref)); scale=max(scale,abs(ref))
      end do; end do; end do; end do
      do x=1,4
        tol=1.0e-11_dp*max(1.0_dp,scale(x))
        write(*,'(A,A,A,ES14.6,A,ES14.6)') 'Dense ',names(x),' defect: ',err(x),' scale: ',scale(x)
        if (.not.(err(x)<=tol)) error stop 'W dense-reference validation failed'
      end do
      write(*,'(A)') 'W dense-reference validation PASS'
    end do
    deallocate(s1,s2,wo,wv,wx,wy,v,t,u)
  end subroutine test_ccsd_w_dense




  subroutine build_ccsd_l(blk,g,t1,Foo,Fvv,Loo,Lvv)
      
      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, N_cell
      use k_blocks, only: K_block

    type(K_block), intent(in) :: blk(N_cell)
    complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(in) :: t1(nfuc,nfuc,N_cell)
    complex(dpc), intent(in) :: Foo(nfuc,nfuc,N_cell),Fvv(nfuc,nfuc,N_cell)
    complex(dpc), intent(out) :: Loo(nfuc,nfuc,N_cell),Lvv(nfuc,nfuc,N_cell)
    integer :: k,kl,kd,m,i,l,c,a,e,d
    ! Canonical HF assumption: bare f_ov=0, matching build_ccsd_f.
    ! Loo(m,i)=Foo(m,i)+sum_lc [2(lc|mi)-(mc|li)]*t1(l,c)
    ! Lvv(a,e)=Fvv(a,e)+sum_md [2(md|ae)-(me|ad)]*t1(m,d)
    ! The bare f_ov here must NOT be replaced by dressed Fov.
    Loo=cmplx(0.0_dp,0.0_dp,kind=dpc)
    Lvv=cmplx(0.0_dp,0.0_dp,kind=dpc)
    do k=1,N_cell
      do m=1,nfuc
        if (blk(k)%occupation(m)/=2) cycle
        do i=1,nfuc
          if (blk(k)%occupation(i)/=2) cycle
          Loo(m,i,k)=Foo(m,i,k)
          do kl=1,N_cell
            do l=1,nfuc
              if (blk(kl)%occupation(l)/=2) cycle
              do c=1,nfuc
                if (blk(kl)%occupation(c)/=0) cycle
                Loo(m,i,k)=Loo(m,i,k)+(2.0_dp*g(l,c,m,i,kl,kl,k) &
                  -g(m,c,l,i,k,kl,kl))*t1(l,c,kl)
              end do
            end do
          end do
        end do
      end do
      do a=1,nfuc
        if (blk(k)%occupation(a)/=0) cycle
        do e=1,nfuc
          if (blk(k)%occupation(e)/=0) cycle
          Lvv(a,e,k)=Fvv(a,e,k)
          do kd=1,N_cell
            do m=1,nfuc
              if (blk(kd)%occupation(m)/=2) cycle
              do d=1,nfuc
                if (blk(kd)%occupation(d)/=0) cycle
                Lvv(a,e,k)=Lvv(a,e,k)+(2.0_dp*g(m,d,a,e,kd,kd,k) &
                  -g(m,e,a,d,kd,k,k))*t1(m,d,kd)
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine build_ccsd_l

  subroutine test_ccsd_l(blk,g,Foo,Fvv)

      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, N_cell
      use k_blocks, only: K_block

    type(K_block), intent(in) :: blk(N_cell)
    complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(in) :: Foo(nfuc,nfuc,N_cell),Fvv(nfuc,nfuc,N_cell)
    complex(dpc) :: singles(nfuc,nfuc,N_cell)
    complex(dpc) :: lo(nfuc,nfuc,N_cell),lv(nfuc,nfuc,N_cell)
    complex(dpc), allocatable :: v(:,:,:,:),u(:,:)
    complex(dpc) :: ref,got
    real(dp) :: defect(2),scale(2),tol,phase
    integer :: ng,p,q,r,s,k,m,i,a,e,l,c,mode,x
    integer :: km(nfuc*N_cell),bn(nfuc*N_cell),oc(nfuc*N_cell)
    integer, external :: fourth_k
    ! Tests corrections relative to supplied F, not F itself.
    ! Thus callers can pass their already validated MP2-built Foo/Fvv.
    singles=cmplx(0.0_dp,0.0_dp,kind=dpc)
    call build_ccsd_l(blk,g,singles,Foo,Fvv,lo,lv)
    defect=0.0_dp
    do k=1,N_cell; do p=1,nfuc; do q=1,nfuc
      if (blk(k)%occupation(p)==2.and.blk(k)%occupation(q)==2) &
        defect(1)=max(defect(1),abs(lo(p,q,k)-Foo(p,q,k)))
      if (blk(k)%occupation(p)==0.and.blk(k)%occupation(q)==0) &
        defect(2)=max(defect(2),abs(lv(p,q,k)-Fvv(p,q,k)))
    end do; end do; end do
    write(*,'(A,ES14.6)') 'L zero-singles Loo-Foo defect: ',defect(1)
    write(*,'(A,ES14.6)') 'L zero-singles Lvv-Fvv defect: ',defect(2)
    if (.not.all(defect<=1.0e-12_dp)) error stop 'L zero-singles validation failed'
    ng=nfuc*N_cell
    allocate(v(ng,ng,ng,ng),u(ng,ng))
    do p=1,ng
      km(p)=(p-1)/nfuc+1; bn(p)=modulo(p-1,nfuc)+1
      oc(p)=blk(km(p))%occupation(bn(p))
    end do
    v=cmplx(0.0_dp,0.0_dp,kind=dpc)
    do p=1,ng; do q=1,ng; do r=1,ng; do s=1,ng
      if (km(s)/=fourth_k(km(p),km(r),km(q))) cycle
      v(p,q,r,s)=g(bn(p),bn(q),bn(r),bn(s),km(p),km(q),km(r))
    end do; end do; end do; end do
    u=cmplx(0.0_dp,0.0_dp,kind=dpc)
    do i=1,ng; do a=1,ng
      if (oc(i)/=2.or.oc(a)/=0.or.km(i)/=km(a)) cycle
      phase=real(3*i+5*a,dp)
      u(i,a)=1.0e-3_dp*cmplx(sin(phase),cos(phase),kind=dpc)
      singles(bn(i),bn(a),km(i))=u(i,a)
    end do; end do
    call build_ccsd_l(blk,g,singles,Foo,Fvv,lo,lv)
    defect=0.0_dp; scale=0.0_dp
    ! Compare only L-F corrections, including forbidden k off-diagonals.
    do m=1,ng; do i=1,ng
      if (oc(m)/=2.or.oc(i)/=2) cycle
      ref=cmplx(0.0_dp,0.0_dp,kind=dpc)
      do l=1,ng; do c=1,ng
        if (oc(l)/=2.or.oc(c)/=0) cycle
        ref=ref+(2.0_dp*v(l,c,m,i)-v(m,c,l,i))*u(l,c)
      end do; end do
      got=cmplx(0.0_dp,0.0_dp,kind=dpc)
      if (km(m)==km(i)) got=lo(bn(m),bn(i),km(m))-Foo(bn(m),bn(i),km(m))
      defect(1)=max(defect(1),abs(got-ref)); scale(1)=max(scale(1),abs(ref))
    end do; end do
    do a=1,ng; do e=1,ng
      if (oc(a)/=0.or.oc(e)/=0) cycle
      ref=cmplx(0.0_dp,0.0_dp,kind=dpc)
      do m=1,ng; do c=1,ng
        if (oc(m)/=2.or.oc(c)/=0) cycle
        ref=ref+(2.0_dp*v(m,c,a,e)-v(m,e,a,c))*u(m,c)
      end do; end do
      got=cmplx(0.0_dp,0.0_dp,kind=dpc)
      if (km(a)==km(e)) got=lv(bn(a),bn(e),km(a))-Fvv(bn(a),bn(e),km(a))
      defect(2)=max(defect(2),abs(got-ref)); scale(2)=max(scale(2),abs(ref))
    end do; end do
    write(*,'(A,ES14.6,A,ES14.6)') 'Dense Loo correction defect: ',defect(1),' scale: ',scale(1)
    write(*,'(A,ES14.6,A,ES14.6)') 'Dense Lvv correction defect: ',defect(2),' scale: ',scale(2)
    do x=1,2
      tol=1.0e-11_dp*max(1.0_dp,scale(x))
      if (.not.(defect(x)<=tol)) error stop 'L dense-reference validation failed'
    end do
    write(*,'(A)') 'L dense-reference validation PASS'
    deallocate(v,u)
  end subroutine test_ccsd_l


subroutine build_ccsd_residuals(blk,g,t1,t2,Foo,Fvv,Fov,Loo,Lvv, &
                                 Woooo,Wvvvv,Wvoov,Wvovo,R1,R2)

      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, N_cell
      use k_blocks, only: K_block

    type(K_block), intent(in) :: blk(N_cell)
    complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(in) :: t1(nfuc,nfuc,N_cell)
    complex(dpc), intent(in) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(in) :: Foo(nfuc,nfuc,N_cell),Fvv(nfuc,nfuc,N_cell)
    complex(dpc), intent(in) :: Fov(nfuc,nfuc,N_cell)
    complex(dpc), intent(in) :: Loo(nfuc,nfuc,N_cell),Lvv(nfuc,nfuc,N_cell)
    complex(dpc), intent(in) :: Woooo(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(in) :: Wvvvv(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(in) :: Wvoov(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(in) :: Wvovo(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(out) :: R1(nfuc,nfuc,N_cell)
    complex(dpc), intent(out) :: R2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    integer :: occ(nfuc*N_cell),vir(nfuc*N_cell),no,nv,p,k,i,j,a,b,c,d,l,m
    integer :: ii,jj,aa,bb,cc,dd,ll,mm
    integer, external :: fourth_k
    complex(dpc) :: v
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
    do ii=1,no; do aa=1,nv
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
          v=v+f2(Fov,m,c)*(2.0_dp*s2(m,i,c,a)-s2(i,m,c,a)+s1(i,c)*s1(m,a)) &
            +(2.0_dp*eri(m,c,a,i)-eri(m,i,a,c))*s1(m,c)
          do dd=1,nv
            d=vir(dd)
            v=v+(2.0_dp*eri(m,d,a,c)-eri(m,c,a,d))* &
              (s2(i,m,c,d)+s1(m,d)*s1(i,c))
          end do
          do ll=1,no
            l=occ(ll)
            v=v+(-2.0_dp*eri(l,c,m,i)+eri(m,c,l,i))* &
              (s2(m,l,a,c)+s1(l,c)*s1(m,a))
          end do
        end do
      end do
      R1(band(i),band(a),mom(i))=v
    end do; end do
    do ii=1,no; do jj=1,no; do aa=1,nv; do bb=1,nv
      i=occ(ii); j=occ(jj); a=vir(aa); b=vir(bb)
      if (mom(b)/=fourth_k(mom(i),mom(j),mom(a))) cycle
      v=conjg(eri(i,a,j,b))+pair(i,j,a,b)+pair(j,i,b,a) &
        -(blk(mom(i))%energy(band(i))+blk(mom(j))%energy(band(j)) &
          -blk(mom(a))%energy(band(a))-blk(mom(b))%energy(band(b)))*s2(i,j,a,b)
      do mm=1,no; do ll=1,no
        m=occ(mm); l=occ(ll)
        v=v+w4(Woooo,m,l,i,j)*(s2(m,l,a,b)+s1(m,a)*s1(l,b))
      end do; end do
      do cc=1,nv; do dd=1,nv
        c=vir(cc); d=vir(dd)
        v=v+w4(Wvvvv,a,b,c,d)*(s2(i,j,c,d)+s1(i,c)*s1(j,d))
      end do; end do
      R2(band(i),band(j),band(a),band(b),mom(i),mom(j),mom(a))=v
    end do; end do; end do; end do
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

  subroutine test_ccsd_residual_zero(blk,g)
    use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, N_cell
      use k_blocks, only: K_block
    type(K_block), intent(in) :: blk(N_cell)
    complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), allocatable :: z1(:,:,:),z2(:,:,:,:,:,:,:),r1(:,:,:),r2(:,:,:,:,:,:,:)
    complex(dpc), allocatable :: wo(:,:,:,:,:,:,:),wv(:,:,:,:,:,:,:)
    complex(dpc), allocatable :: wx(:,:,:,:,:,:,:),wy(:,:,:,:,:,:,:)
    complex(dpc) :: fo(nfuc,nfuc,N_cell),fv(nfuc,nfuc,N_cell),fx(nfuc,nfuc,N_cell)
    complex(dpc) :: lo(nfuc,nfuc,N_cell),lv(nfuc,nfuc,N_cell),ref
    real(dp) :: defect
    integer :: k,ki,kj,ka,kb,i,j,a,b
    integer, external :: fourth_k
    allocate(z1(nfuc,nfuc,N_cell),r1(nfuc,nfuc,N_cell))
    allocate(z2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell))
    allocate(r2,mold=z2); allocate(wo,mold=z2); allocate(wv,mold=z2)
    allocate(wx,mold=z2); allocate(wy,mold=z2)
    z1=cmplx(0.0_dp,0.0_dp,kind=dpc); z2=cmplx(0.0_dp,0.0_dp,kind=dpc)
    fo=cmplx(0.0_dp,0.0_dp,kind=dpc); fv=fo; fx=fo
    do k=1,N_cell; do i=1,nfuc
      if (blk(k)%occupation(i)==2) fo(i,i,k)=blk(k)%energy(i)
      if (blk(k)%occupation(i)==0) fv(i,i,k)=blk(k)%energy(i)
    end do; end do
    call build_ccsd_l(blk,g,z1,fo,fv,lo,lv)
    call build_ccsd_w(blk,g,z1,z2,wo,wv,wx,wy)
    call build_ccsd_residuals(blk,g,z1,z2,fo,fv,fx,lo,lv,wo,wv,wx,wy,r1,r2)
    defect=0.0_dp
    do ki=1,N_cell; do kj=1,N_cell; do ka=1,N_cell
      kb=fourth_k(ki,kj,ka)
      do i=1,nfuc; do j=1,nfuc; do a=1,nfuc; do b=1,nfuc
        ref=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (blk(ki)%occupation(i)==2.and.blk(kj)%occupation(j)==2.and. &
            blk(ka)%occupation(a)==0.and.blk(kb)%occupation(b)==0) &
          ref=conjg(g(i,a,j,b,ki,ka,kj))
        defect=max(defect,abs(r2(i,j,a,b,ki,kj,ka)-ref))
      end do; end do; end do; end do
    end do; end do; end do
    write(*,'(A,ES14.6)') 'CCSD zero-amplitude R1 defect: ',maxval(abs(r1))
    write(*,'(A,ES14.6)') 'CCSD zero-amplitude R2 driving defect: ',defect
    if (.not.(maxval(abs(r1))<=1.0e-12_dp.and.defect<=1.0e-12_dp)) &
      error stop 'CCSD zero-amplitude residual validation failed'
    deallocate(z1,z2,r1,r2,wo,wv,wx,wy)
  end subroutine test_ccsd_residual_zero




  subroutine export_ccsd_residual_tests(blk,g,t1,t2,build_f)
      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, N_cell
      use k_blocks, only: K_block
    use iso_fortran_env, only: int32,real64
    type(K_block), intent(in) :: blk(N_cell)
    complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(in) :: t1(nfuc,nfuc,N_cell)
    complex(dpc), intent(in) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    interface
      subroutine build_f(blk,g,t1,t2,Foo,Fvv,Fov)
        import :: K_block,nfuc,N_cell,dpc
        type(K_block), intent(in) :: blk(N_cell)
        complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
        complex(dpc), intent(in) :: t1(nfuc,nfuc,N_cell)
        complex(dpc), intent(in) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
        complex(dpc), intent(out) :: Foo(nfuc,nfuc,N_cell),Fvv(nfuc,nfuc,N_cell)
        complex(dpc), intent(out) :: Fov(nfuc,nfuc,N_cell)
      end subroutine
    end interface
    complex(dpc), allocatable :: s1(:,:,:),s2(:,:,:,:,:,:,:),r1(:,:,:),r2(:,:,:,:,:,:,:)
    complex(dpc), allocatable :: wo(:,:,:,:,:,:,:),wv(:,:,:,:,:,:,:)
    complex(dpc), allocatable :: wx(:,:,:,:,:,:,:),wy(:,:,:,:,:,:,:)
    complex(dpc) :: fo(nfuc,nfuc,N_cell),fv(nfuc,nfuc,N_cell),fx(nfuc,nfuc,N_cell)
    complex(dpc) :: lo(nfuc,nfuc,N_cell),lv(nfuc,nfuc,N_cell)
    integer(int32) :: occupations(nfuc,N_cell)
    real(real64) :: energies(nfuc,N_cell)
    integer :: mode,k,ki,kj,ka,kb,i,j,a,b,gi,gj,ga,gb,unit
    integer, external :: fourth_k
    real(dp) :: phase
    character(len=32) :: filename
    if (storage_size(cmplx(0.0_dp,0.0_dp,kind=dpc))/=128) &
      error stop 'Snapshot export requires 128-bit complex'
    allocate(s1,source=t1); allocate(s2,source=t2)
    allocate(r1,mold=t1); allocate(r2,mold=t2)
    allocate(wo,mold=t2); allocate(wv,mold=t2)
    allocate(wx,mold=t2); allocate(wy,mold=t2)
    do k=1,N_cell
      occupations(:,k)=int(blk(k)%occupation,int32)
      energies(:,k)=real(blk(k)%energy,real64)
    end do
    do mode=1,2
      if (mode==1) then
        filename='ccsd_mp2.bin'
      else
        filename='ccsd_synthetic.bin'
        s1=cmplx(0.0_dp,0.0_dp,kind=dpc); s2=cmplx(0.0_dp,0.0_dp,kind=dpc)
        do ki=1,N_cell; do i=1,nfuc; do a=1,nfuc
          if (blk(ki)%occupation(i)/=2.or.blk(ki)%occupation(a)/=0) cycle
          gi=(ki-1)*nfuc+i; ga=(ki-1)*nfuc+a
          phase=real(3*gi+5*ga,dp)
          s1(i,a,ki)=1.0e-3_dp*cmplx(sin(phase),cos(phase),kind=dpc)
        end do; end do; end do
        do ki=1,N_cell; do kj=1,N_cell; do ka=1,N_cell
          kb=fourth_k(ki,kj,ka)
          do i=1,nfuc; do j=1,nfuc; do a=1,nfuc; do b=1,nfuc
            if (blk(ki)%occupation(i)/=2.or.blk(kj)%occupation(j)/=2) cycle
            if (blk(ka)%occupation(a)/=0.or.blk(kb)%occupation(b)/=0) cycle
            gi=(ki-1)*nfuc+i; gj=(kj-1)*nfuc+j
            ga=(ka-1)*nfuc+a; gb=(kb-1)*nfuc+b
            phase=real(gi+gj+3*(ga+gb)+gi*ga+gj*gb,dp)
            s2(i,j,a,b,ki,kj,ka)=1.0e-3_dp*cmplx(sin(phase),cos(phase),kind=dpc)
          end do; end do; end do; end do
        end do; end do; end do
      end if
      call build_f(blk,g,s1,s2,fo,fv,fx)
      call build_ccsd_l(blk,g,s1,fo,fv,lo,lv)
      call build_ccsd_w(blk,g,s1,s2,wo,wv,wx,wy)
      call build_ccsd_residuals(blk,g,s1,s2,fo,fv,fx,lo,lv,wo,wv,wx,wy,r1,r2)
      ! Native-endian unformatted STREAM: no record markers.
      ! Explicit fourth-k table prevents assuming the k-point enumeration.
      open(newunit=unit,file=trim(filename),access='stream',form='unformatted',status='replace')
      write(unit) int(nfuc,int32),int(N_cell,int32)
      write(unit) occupations,energies
      do ka=1,N_cell; do kj=1,N_cell; do ki=1,N_cell
        write(unit) int(fourth_k(ki,kj,ka),int32)
      end do; end do; end do
      write(unit) g,s1,s2,r1,r2
      close(unit)
      write(*,'(A,A)') 'Residual snapshot written: ',trim(filename)
    end do
    deallocate(s1,s2,r1,r2,wo,wv,wx,wy)
  end subroutine



  subroutine solve_kccsd(blk,g,t1,t2,build_f,ecorr,eimag,converged, &
                        max_iter,damping,energy_tol,residual_tol,ndiis)
      use constants_module, only: dp, dpc
      use unitcell_module, only: nfuc, N_cell
      use k_blocks, only: K_block

      interface
        subroutine diis_coefficients(errors,c,ok)
          use constants_module, only: dp,dpc
          implicit none
          complex(dpc), intent(in) :: errors(:,:)
          real(dp), intent(out) :: c(:)
          logical, intent(out) :: ok
        end subroutine diis_coefficients
      end interface

    type(K_block), intent(in) :: blk(N_cell)
    complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    complex(dpc), intent(inout) :: t1(nfuc,nfuc,N_cell)
    complex(dpc), intent(inout) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
    interface
      subroutine build_f(blk,g,t1,t2,Foo,Fvv,Fov)
        import :: K_block,nfuc,N_cell,dpc
        type(K_block), intent(in) :: blk(N_cell)
        complex(dpc), intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
        complex(dpc), intent(in) :: t1(nfuc,nfuc,N_cell)
        complex(dpc), intent(in) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
        complex(dpc), intent(out) :: Foo(nfuc,nfuc,N_cell),Fvv(nfuc,nfuc,N_cell)
        complex(dpc), intent(out) :: Fov(nfuc,nfuc,N_cell)
      end subroutine
    end interface
    real(dp), intent(out) :: ecorr,eimag
    logical, intent(out) :: converged
    integer, optional, intent(in) :: max_iter,ndiis
    real(dp), optional, intent(in) :: damping,energy_tol,residual_tol
    complex(dpc), allocatable :: wo(:,:,:,:,:,:,:),wv(:,:,:,:,:,:,:)
    complex(dpc), allocatable :: wx(:,:,:,:,:,:,:),wy(:,:,:,:,:,:,:)
    complex(dpc), allocatable :: r1(:,:,:),r2(:,:,:,:,:,:,:),x(:),err(:),trial(:)
    complex(dpc), allocatable :: hx(:,:),he(:,:)
    real(dp), allocatable :: denom(:),coef(:)
    integer, allocatable :: ix(:,:),jx(:,:)
    complex(dpc) :: fo(nfuc,nfuc,N_cell),fv(nfuc,nfuc,N_cell),fx(nfuc,nfuc,N_cell)
    complex(dpc) :: lo(nfuc,nfuc,N_cell),lv(nfuc,nfuc,N_cell),energy,tau
    integer :: n1,n2,n,p,q,k,ki,kj,ka,kb,i,j,a,b,iter,limit,nd,h,nh
    integer, external :: fourth_k
    real(dp) :: alpha,etol,rtol,previous,de,norm1,norm2,max1,max2
    logical :: ok
    limit=100; nd=6; alpha=0.7_dp; etol=1.0e-9_dp; rtol=1.0e-7_dp
    if (present(max_iter)) limit=max_iter
    if (present(ndiis)) nd=ndiis
    if (present(damping)) alpha=damping
    if (present(energy_tol)) etol=energy_tol
    if (present(residual_tol)) rtol=residual_tol
    if (limit<1.or.nd<0.or.alpha<=0.0_dp.or.alpha>1.0_dp.or.etol<=0.or.rtol<=0) &
      error stop 'Invalid CCSD solver parameters'
    n1=0; n2=0
    do k=1,N_cell; do i=1,nfuc; do a=1,nfuc
      if (blk(k)%occupation(i)==2.and.blk(k)%occupation(a)==0) n1=n1+1
    end do; end do; end do
    do ki=1,N_cell; do kj=1,N_cell; do ka=1,N_cell
      kb=fourth_k(ki,kj,ka)
      do i=1,nfuc; do j=1,nfuc; do a=1,nfuc; do b=1,nfuc
        if (blk(ki)%occupation(i)/=2.or.blk(kj)%occupation(j)/=2) cycle
        if (blk(ka)%occupation(a)/=0.or.blk(kb)%occupation(b)/=0) cycle
        n2=n2+1
      end do; end do; end do; end do
    end do; end do; end do
    n=n1+n2
    if (n1==0.or.n2==0) error stop 'CCSD solver requires occupied and virtual orbitals'
    allocate(ix(3,n1),jx(7,n2),denom(n),x(n),err(n),trial(n))
    p=0
    do k=1,N_cell; do i=1,nfuc; do a=1,nfuc
      if (blk(k)%occupation(i)/=2.or.blk(k)%occupation(a)/=0) cycle
      p=p+1; ix(:,p)=[i,a,k]
      denom(p)=blk(k)%energy(i)-blk(k)%energy(a)
    end do; end do; end do
    q=0
    do ki=1,N_cell; do kj=1,N_cell; do ka=1,N_cell
      kb=fourth_k(ki,kj,ka)
      do i=1,nfuc; do j=1,nfuc; do a=1,nfuc; do b=1,nfuc
        if (blk(ki)%occupation(i)/=2.or.blk(kj)%occupation(j)/=2) cycle
        if (blk(ka)%occupation(a)/=0.or.blk(kb)%occupation(b)/=0) cycle
        q=q+1; jx(:,q)=[i,j,a,b,ki,kj,ka]
        denom(n1+q)=blk(ki)%energy(i)+blk(kj)%energy(j) &
                      -blk(ka)%energy(a)-blk(kb)%energy(b)
      end do; end do; end do; end do
    end do; end do; end do
    if (minval(abs(denom))<1.0e-10_dp) error stop 'CCSD near-zero denominator'
    allocate(r1,mold=t1); allocate(r2,mold=t2)
    allocate(wo,mold=t2); allocate(wv,mold=t2)
    allocate(wx,mold=t2); allocate(wy,mold=t2)
    allocate(hx(n,max(1,nd)),he(n,max(1,nd)),coef(max(1,nd)))
    nh=0; converged=.false.; previous=0.0_dp
    ! Pack and clear inactive entries, preserving only physical amplitudes.
    call pack(t1,t2,x)
    call unpack(x)
    write(*,'(A)') 'CCSD: iter Ecorr dE maxR1 maxR2 RMS_R1 RMS_R2 Im(E)'
    do iter=0,limit
      call build_f(blk,g,t1,t2,fo,fv,fx)
      call build_ccsd_l(blk,g,t1,fo,fv,lo,lv)
      call build_ccsd_w(blk,g,t1,t2,wo,wv,wx,wy)
      call build_ccsd_residuals(blk,g,t1,t2,fo,fv,fx,lo,lv,wo,wv,wx,wy,r1,r2)
      call pack(r1,r2,err)
      max1=maxval(abs(err(:n1))); max2=maxval(abs(err(n1+1:)))
      norm1=sqrt(sum(abs(err(:n1))**2)/real(n1,dp))
      norm2=sqrt(sum(abs(err(n1+1:))**2)/real(n2,dp))
      energy=cmplx(0.0_dp,0.0_dp,kind=dpc)
      do q=1,n2
        i=jx(1,q); j=jx(2,q); a=jx(3,q); b=jx(4,q)
        ki=jx(5,q); kj=jx(6,q); ka=jx(7,q); kb=fourth_k(ki,kj,ka)
        tau=t2(i,j,a,b,ki,kj,ka)
        if (ki==ka.and.kj==kb) tau=tau+t1(i,a,ki)*t1(j,b,kj)
        energy=energy+(2.0_dp*g(i,a,j,b,ki,ka,kj)-g(i,b,j,a,ki,kb,kj))*tau
      end do
      ecorr=real(energy,dp); eimag=aimag(energy)
      de=0.0_dp
      if (iter>0) de=ecorr-previous
      write(*,'(I4,1X,F20.12,6(1X,ES12.4))') iter,ecorr,de,max1,max2,norm1,norm2,eimag
      if (.not.(max1<huge(1.0_dp).and.max2<huge(1.0_dp).and.abs(ecorr)<huge(1.0_dp))) &
        error stop 'CCSD nonfinite residual or energy'
      if (max(max1,max2)<rtol.and.(iter==0.or.abs(de)<etol)) then
        converged=.true.; exit
      end if
      if (iter==limit) exit
      call pack(t1,t2,x)
      ! DIIS pairs CURRENT amplitudes with their CURRENT raw residuals.
      ! Extrapolate a Jacobi step from the mixed pair. With DIIS disabled
      ! this is exactly x_new=x+alpha*R/D.
      if (nd>=2) then
        if (nh==nd) then
          hx(:,1:nd-1)=hx(:,2:nd); he(:,1:nd-1)=he(:,2:nd)
          nh=nd-1
        end if
        nh=nh+1; hx(:,nh)=x; he(:,nh)=err
        if (nh>=2.and.iter>=2) then
          call diis_coefficients(he(:,1:nh),coef(1:nh),ok)
          if (ok) then
            trial=cmplx(0.0_dp,0.0_dp,kind=dpc)
            do h=1,nh
              trial=trial+coef(h)*(hx(:,h)+alpha*he(:,h)/denom)
            end do
          else
            trial=x+alpha*err/denom
            nh=1; hx(:,1)=x; he(:,1)=err
            write(*,'(A)') 'CCSD DIIS reset: ill-conditioned subspace'
          end if
        else
          trial=x+alpha*err/denom
        end if
      else
        trial=x+alpha*err/denom
      end if
      previous=ecorr
      call unpack(trial)
    end do
    if (converged) then
      write(*,'(A)') 'CCSD CONVERGED'
    else
      write(*,'(A)') 'CCSD NOT CONVERGED: iteration limit reached'
    end if
    deallocate(ix,jx,denom,x,err,trial,hx,he,coef,r1,r2,wo,wv,wx,wy)
  contains
    subroutine pack(a1,a2,v)
      complex(dpc), intent(in) :: a1(nfuc,nfuc,N_cell)
      complex(dpc), intent(in) :: a2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: v(n)
      integer :: p,q
      do p=1,n1
        v(p)=a1(ix(1,p),ix(2,p),ix(3,p))
      end do
      do q=1,n2
        v(n1+q)=a2(jx(1,q),jx(2,q),jx(3,q),jx(4,q),jx(5,q),jx(6,q),jx(7,q))
      end do
    end subroutine
    subroutine unpack(v)
      complex(dpc), intent(in) :: v(n)
      integer :: p,q
      t1=cmplx(0.0_dp,0.0_dp,kind=dpc); t2=cmplx(0.0_dp,0.0_dp,kind=dpc)
      do p=1,n1
        t1(ix(1,p),ix(2,p),ix(3,p))=v(p)
      end do
      do q=1,n2
        t2(jx(1,q),jx(2,q),jx(3,q),jx(4,q),jx(5,q),jx(6,q),jx(7,q))=v(n1+q)
      end do
    end subroutine
  end subroutine solve_kccsd

  subroutine diis_coefficients(errors,c,ok)
          use constants_module, only: dp, dpc
    complex(dpc), intent(in) :: errors(:,:)
    real(dp), intent(out) :: c(:)
    logical, intent(out) :: ok
    real(dp) :: mat(size(c)+1,size(c)+1),rhs(size(c)+1),row(size(c)+1)
    real(dp) :: factor,pivot,scale,tmp
    integer :: m,i,j,k,p
    m=size(c); ok=.false.; c=0.0_dp
    mat=0.0_dp; rhs=0.0_dp
    do i=1,m; do j=1,m
      mat(i,j)=real(sum(conjg(errors(:,i))*errors(:,j)),dp)
    end do; end do
    scale=maxval(abs(mat(:m,:m)))
    if (scale<=tiny(1.0_dp)) return
    mat(:m,:m)=mat(:m,:m)/scale
    do i=1,m
      mat(i,i)=mat(i,i)+1.0e-12_dp
    end do
    mat(:m,m+1)=-1.0_dp; mat(m+1,:m)=-1.0_dp; rhs(m+1)=-1.0_dp
    do k=1,m+1
      p=k-1+maxloc(abs(mat(k:m+1,k)),dim=1)
      if (abs(mat(p,k))<1.0e-14_dp) return
      row=mat(k,:); mat(k,:)=mat(p,:); mat(p,:)=row
      tmp=rhs(k); rhs(k)=rhs(p); rhs(p)=tmp
      do i=k+1,m+1
        factor=mat(i,k)/mat(k,k)
        mat(i,k:m+1)=mat(i,k:m+1)-factor*mat(k,k:m+1)
        rhs(i)=rhs(i)-factor*rhs(k)
      end do
    end do
    do i=m+1,1,-1
      rhs(i)=(rhs(i)-sum(mat(i,i+1:m+1)*rhs(i+1:m+1)))/mat(i,i)
    end do
    c=rhs(:m)
    ok=all(abs(c)<100.0_dp).and.abs(sum(c)-1.0_dp)<1.0e-8_dp
  end subroutine diis_coefficients