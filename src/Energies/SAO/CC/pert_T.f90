subroutine build_ccsd_t(blk,g,ktab,t1,t2,et)

      use constants_module, only: dp,dpc
      use unitcell_module , only: nfuc,N_cell
      use k_blocks        , only: K_block
      
      implicit none
      ! input !
      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: t1(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      integer      , intent(in) :: ktab(N_cell,N_cell,N_cell)
      ! output !
      real(dp)     , intent(out):: et
      ! local !
      integer, parameter :: wt(6,6)=reshape([ &
      1,2,3,4,5,6, 2,1,5,6,3,4, 3,4,1,2,6,5, &
      4,3,6,5,1,2, 5,6,2,1,4,3, 6,5,4,3,2,1 ],[6,6])
      integer  :: no,nv,k,x,n,nt,p,q,Ktot,a,b,c,ii,jj,kk
      integer  :: ka,kb,kc,ki,kj,kl,kf,km,ba,bb,bc,bi,bj,bl
      integer  :: li,lj,ll,lf,lm
      integer  :: pidx(3,6),vk(3,6),vbd(3,6),oi(3,6)
      integer  :: nOcc(N_cell),nVir(N_cell)
      integer  :: oBand(nfuc,N_cell),oPos(nfuc,N_cell),vBand(nfuc,N_cell)
      integer  :: oK(nfuc*N_cell),vKa(nfuc*N_cell),vBa(nfuc*N_cell)
      real(dp) :: epsO(nfuc*N_cell),epsV(nfuc*N_cell)
      real(dp) :: deg,eV,denom
      integer , allocatable :: nTri(:),tri(:,:,:)
      real(dp), allocatable :: tD(:,:)
      complex(dpc) :: energy,acc,vij,gF(nfuc),gM(nfuc)
      complex(dpc), allocatable :: w(:,:,:,:),v(:,:,:,:),z(:,:,:,:)

      ! ---- orbital lists (local band / position / energy) ----
      no=0; nv=0; nOcc=0; nVir=0
      do k=1,N_cell
        do x=1,nfuc
          if (blk(k)%occupation(x)==2) then
            no=no+1
            nOcc(k)=nOcc(k)+1
            oBand(nOcc(k),k)=x
            oPos(nOcc(k),k)=no
            oK(no)=k
            epsO(no)=blk(k)%energy(x)
          else if (blk(k)%occupation(x)==0) then
            nv=nv+1
            nVir(k)=nVir(k)+1
            vBand(nVir(k),k)=x
            vKa(nv)=k
            vBa(nv)=x
            epsV(nv)=blk(k)%energy(x)
          else
            error stop '(T): occupation must be 0 or 2'
          end if
        end do
      end do
      et=0.0_dp
      if (no==0.or.nv==0) return

      ! ---- conserving occupied triples, keyed by total momentum ----
      ! key S(x,y,z) = x+y+z-2*k1 (k1 = index 1), equal keys <=> conservation
      allocate(nTri(N_cell),tri(3,no**3,N_cell),tD(no**3,N_cell))
      nTri=0
      do ii=1,no
        do jj=1,no
          do kk=1,no
            Ktot=ktab(ktab(oK(ii),oK(jj),1),oK(kk),1)
            nTri(Ktot)=nTri(Ktot)+1
            n=nTri(Ktot)
            tri(:,n,Ktot)=[ii,jj,kk]
            tD(n,Ktot)=epsO(ii)+epsO(jj)+epsO(kk)
          end do
        end do
      end do

      energy=cmplx(0.0_dp,0.0_dp,kind=dpc)

      !$omp parallel default(shared) &
      !$omp private(b,c,n,nt,p,q,x,Ktot,ii,jj,kk,ka,kb,kc,ki,kj,kl,kf,km, &
      !$omp         ba,bb,bc,bi,bj,bl,li,lj,ll,lf,lm,pidx,vk,vbd,oi, &
      !$omp         deg,eV,denom,acc,vij,gF,gM,w,v,z)
      allocate(w(no,no,no,6),v(no,no,no,6),z(no,no,no,6))
      !$omp do schedule(dynamic,1) reduction(+:energy)
      do a=1,nv
        do b=1,a
          do c=1,b
            Ktot=ktab(ktab(vKa(a),vKa(b),1),vKa(c),1)
            nt=nTri(Ktot)
            if (nt==0) cycle

            pidx(:,1)=[a,b,c]; pidx(:,2)=[a,c,b]; pidx(:,3)=[b,a,c]
            pidx(:,4)=[b,c,a]; pidx(:,5)=[c,a,b]; pidx(:,6)=[c,b,a]
            do p=1,6
              do x=1,3
                vk(x,p)=vKa(pidx(x,p))
                vbd(x,p)=vBa(pidx(x,p))
              end do
            end do

            ! ---- W and V on conserving slots only (no zeroing needed) ----
            do p=1,6
              ka=vk(1,p); kb=vk(2,p); kc=vk(3,p)
              ba=vbd(1,p); bb=vbd(2,p); bc=vbd(3,p)
              do ki=1,N_cell
                if (nOcc(ki)==0) cycle
                kf=ktab(ka,kb,ki)
                do kj=1,N_cell
                  if (nOcc(kj)==0) cycle
                  kl=ktab(kf,kc,kj)
                  if (nOcc(kl)==0) cycle
                  km=ktab(ki,kj,ka)
                  do li=1,nOcc(ki)
                    bi=oBand(li,ki); ii=oPos(li,ki)
                    ! eri(i,aa,f,bb): k_f = k_a + k_b - k_i
                    do lf=1,nVir(kf)
                      gF(lf)=conjg(g(bi,ba,vBand(lf,kf),bb,ki,ka,kf))
                    end do
                    do lj=1,nOcc(kj)
                      bj=oBand(lj,kj); jj=oPos(lj,kj)
                      ! eri(i,aa,j,m): k_m = k_i + k_j - k_a
                      do lm=1,nOcc(km)
                        gM(lm)=conjg(g(bi,ba,bj,oBand(lm,km),ki,ka,kj))
                      end do
                      ! eri(i,aa,j,bb)*t1(l,cc) survives only if k_l = k_c
                      if (kl==kc) then
                        vij=conjg(g(bi,ba,bj,bb,ki,ka,kj))
                      else
                        vij=cmplx(0.0_dp,0.0_dp,kind=dpc)
                      end if
                      do ll=1,nOcc(kl)
                        bl=oBand(ll,kl); kk=oPos(ll,kl)
                        acc=cmplx(0.0_dp,0.0_dp,kind=dpc)
                        ! double(l,j,cc,f)
                        do lf=1,nVir(kf)
                          acc=acc+gF(lf)*t2(bl,bj,bc,vBand(lf,kf),kl,kj,kc)
                        end do
                        ! double(m,l,bb,cc)
                        do lm=1,nOcc(km)
                          acc=acc-gM(lm)*t2(oBand(lm,km),bl,bb,bc,km,kl,kb)
                        end do
                        w(ii,jj,kk,p)=acc
                        if (kl==kc) then
                          v(ii,jj,kk,p)=vij*t1(bl,bc,kc)
                        else
                          v(ii,jj,kk,p)=cmplx(0.0_dp,0.0_dp,kind=dpc)
                        end if
                      end do
                    end do
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
            eV=epsV(a)+epsV(b)+epsV(c)

            ! ---- v <- h = w + 0.5 v  (in place, conserving slots) ----
            do n=1,nt
              ii=tri(1,n,Ktot); jj=tri(2,n,Ktot); kk=tri(3,n,Ktot)
              do p=1,6
                v(ii,jj,kk,p)=w(ii,jj,kk,p)+0.5_dp*v(ii,jj,kk,p)
              end do
            end do

            ! ---- Z ----
            do n=1,nt
              ii=tri(1,n,Ktot); jj=tri(2,n,Ktot); kk=tri(3,n,Ktot)
              denom=(tD(n,Ktot)-eV)*deg
              if (abs(denom)<1.0e-12_dp) error stop '(T): near-zero denominator'
              do p=1,6
                z(ii,jj,kk,p)=(4.0_dp*v(ii,jj,kk,p)+v(jj,kk,ii,p)+v(kk,ii,jj,p) &
                  -2.0_dp*v(kk,jj,ii,p)-2.0_dp*v(ii,kk,jj,p)-2.0_dp*v(jj,ii,kk,p))/denom
              end do
            end do

            ! ---- energy contraction ----
            do n=1,nt
              ii=tri(1,n,Ktot); jj=tri(2,n,Ktot); kk=tri(3,n,Ktot)
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
      !$omp end do
      deallocate(w,v,z)
      !$omp end parallel

      energy=2.0_dp*energy
      et=real(energy,dp)
      write(*,*)
      write(*,'(A)')          '*-----------------------------------*'
      write(*,'(A,f16.12,A)') '*-*   T-ccsd = ', et, '    *-*'
      write(*,'(A)')          '*-----------------------------------*'
      write(*,*)
      Flush(6)
      deallocate(nTri,tri,tD)
      
end subroutine build_ccsd_t
