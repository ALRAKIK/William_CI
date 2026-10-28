subroutine build_ccsd_t(blk,g,t1,t2,et)

      use constants_module, only: dp,dpc
      use unitcell_module , only: nfuc,N_cell
      use k_blocks        , only: K_block

      implicit none

      ! input ! 

      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: t1(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)

      ! local ! 

      integer      :: occ(nfuc*N_cell),vir(nfuc*N_cell),no,nv,p,q,k,x
      integer      :: a,b,c,aa,bb,cc,i,j,l,ii,jj,kk,f,m
      integer      :: vp(3,6),oi(3,6),kf,km,kl
      integer      :: fourth_k
      real(dp)     :: deg,denom
      

      integer, parameter :: wt(6,6)=reshape([ &
      1,2,3,4,5,6, 2,1,5,6,3,4, 3,4,1,2,6,5, &
      4,3,6,5,1,2, 5,6,2,1,4,3, 6,5,4,3,2,1 ],[6,6])

      complex(dpc), allocatable :: w(:,:,:,:)
      complex(dpc), allocatable :: v(:,:,:,:) 
      complex(dpc), allocatable :: z(:,:,:,:)

      complex(dpc)              :: energy

      ! output ! 

      real(dp)    , intent(out) :: et

      ! Rows: occupied W layouts ijk,ikj,jik,jki,kij,kji.
      ! Columns: virtual Z layouts abc,acb,bac,bca,cab,cba.
      
      no=0
      nv=0

      do k=1,N_cell
        do x=1,nfuc
          if (blk(k)%occupation(x)==2) then
            no=no+1
            occ(no)=(k-1)*nfuc+x
          else if (blk(k)%occupation(x)==0) then
            nv=nv+1
            vir(nv)=(k-1)*nfuc+x
          else
            error stop '(T): occupation must be 0 or 2'
          end if
        end do
      end do


      et=0.0_dp

      if (no==0.or.nv==0) return
      allocate(w(no,no,no,6),v(no,no,no,6),z(no,no,no,6))
      energy=cmplx(0.0_dp,0.0_dp,kind=dpc)
      ! Identical global virtual ordering to the snapshot reference.
      do a=1,nv
        do b=1,a
          do c=1,b
            vp(:,1)=[vir(a),vir(b),vir(c)]
            vp(:,2)=[vir(a),vir(c),vir(b)]
            vp(:,3)=[vir(b),vir(a),vir(c)]
            vp(:,4)=[vir(b),vir(c),vir(a)]
            vp(:,5)=[vir(c),vir(a),vir(b)]
            vp(:,6)=[vir(c),vir(b),vir(a)]
            w=cmplx(0.0_dp,0.0_dp,kind=dpc)
            v=cmplx(0.0_dp,0.0_dp,kind=dpc)
            do p=1,6
              aa=vp(1,p)
              bb=vp(2,p)
              cc=vp(3,p)
              do ii=1,no
                do jj=1,no
                  do kk=1,no
                    i=occ(ii)
                    j=occ(jj)
                    l=occ(kk)
                    ! Total momentum: ki+kj+kl=ka+kb+kc.
                    kl=fourth_k(mom(aa),mom(bb),mom(i))
                    kl=fourth_k(kl,mom(cc),mom(j))
                    if (mom(l)/=kl) cycle
                    ! conj(ovvv(i,aa,f,bb))*t2(l,j,cc,f).
                    ! Its integral conservation fixes kf=ka+kb-ki.
                    kf=fourth_k(mom(aa),mom(bb),mom(i))
                    do x=1,nfuc
                      if (blk(kf)%occupation(x)/=0) cycle
                      f=(kf-1)*nfuc+x
                      w(ii,jj,kk,p)=w(ii,jj,kk,p)+conjg(eri(i,aa,f,bb))*double(l,j,cc,f)
                    end do
                    ! -conj(ovoo(i,aa,j,m))*t2(m,l,bb,cc).
                    km=fourth_k(mom(i),mom(j),mom(aa))
                    do x=1,nfuc
                      if (blk(km)%occupation(x)/=2) cycle
                      m=(km-1)*nfuc+x
                      w(ii,jj,kk,p)=w(ii,jj,kk,p)-conjg(eri(i,aa,j,m))*double(m,l,bb,cc)
                    end do
                    ! Canonical HF: bare virtual-occupied Fock block is zero.
                    v(ii,jj,kk,p)=conjg(eri(i,aa,j,bb))*single(l,cc)
                  end do
                end do
              end do
            end do
            deg=1.0_dp
            if (a==c) then
              deg=6.0_dp
            else if (a==b.or.b==c) then
              deg=2.0_dp
            end if
            do p=1,6
              do ii=1,no
                do jj=1,no
                  do kk=1,no
                    i=occ(ii); j=occ(jj); l=occ(kk)
                    denom=(eps(i)+eps(j)+eps(l)-eps(vir(a))-eps(vir(b))-eps(vir(c)))*deg
                    if (abs(denom)<1.0e-12_dp) error stop '(T): near-zero denominator'
                      z(ii,jj,kk,p)=(4.0_dp*h(ii,jj,kk,p)+h(jj,kk,ii,p)+h(kk,ii,jj,p) &
                        -2.0_dp*h(kk,jj,ii,p)-2.0_dp*h(ii,kk,jj,p)-2.0_dp*h(jj,ii,kk,p))/denom
                  end do
                end do
              end do
            end do
            do ii=1,no
              do jj=1,no
                do kk=1,no
                  oi(:,1)=[ii,jj,kk]; oi(:,2)=[ii,kk,jj]; oi(:,3)=[jj,ii,kk]
                  oi(:,4)=[jj,kk,ii]; oi(:,5)=[kk,ii,jj]; oi(:,6)=[kk,jj,ii]
                  do p=1,6
                    do q=1,6
                      energy=energy+w(oi(1,q),oi(2,q),oi(3,q),wt(q,p))*conjg(z(ii,jj,kk,p))
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      
      energy=2.0_dp*energy
      et=real(energy,dp)

      write(*,*)
      write(*,'(A)')          '*-----------------------------------*'
      write(*,'(A,f16.12,A)') '*-*   T-ccsd = ', et, '    *-*'
      write(*,'(A)')          '*-----------------------------------*'
      write(*,*)
      Flush(6)

      deallocate(w,v,z)

      contains
      
      integer function mom(x)
        integer, intent(in) :: x
        mom=(x-1)/nfuc+1
      end function

      integer function band(x)
        integer, intent(in) :: x
        band=modulo(x-1,nfuc)+1
      end function

      real(dp) function eps(x)
        integer, intent(in) :: x
        eps=blk(mom(x))%energy(band(x))
      end function

      function eri(p,q,r,s) result(value)
        integer, intent(in) :: p,q,r,s
        complex(dpc) :: value
        value=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (mom(s)/=fourth_k(mom(p),mom(r),mom(q))) return
        value=g(band(p),band(q),band(r),band(s),mom(p),mom(q),mom(r))
      end function

      function single(i,a) result(value)
        integer, intent(in) :: i,a
        complex(dpc) :: value
        value=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (mom(i)==mom(a)) value=t1(band(i),band(a),mom(i))
      end function

      function double(i,j,a,b) result(value)
        integer, intent(in) :: i,j,a,b
        complex(dpc) :: value
        value=cmplx(0.0_dp,0.0_dp,kind=dpc)
        if (mom(b)/=fourth_k(mom(i),mom(j),mom(a))) return
        value=t2(band(i),band(j),band(a),band(b),mom(i),mom(j),mom(a))
      end function

      function h(i,j,k,p) result(value)
        integer, intent(in) :: i,j,k,p
        complex(dpc) :: value
        value=w(i,j,k,p)+0.5_dp*v(i,j,k,p)
      end function

end subroutine build_ccsd_t