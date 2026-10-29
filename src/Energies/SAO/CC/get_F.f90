subroutine build_ccsd_f(blk,g,ktab,t1,t2,Foo,Fvv,Fov)

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

      integer      :: k,km,kn,ke,kf,m,n,i,a,e,f
      integer      :: fourth_k
      complex(dpc) :: v,tau

      integer :: kOcc(nfuc,N_cell), nOcc(N_cell)
      integer :: kVir(nfuc,N_cell), nVir(N_cell)
      
      integer :: ii , mm , nn , ee , ff, aa
      integer :: p 
      integer :: kp , kq , kr

      ! output !

      complex(dpc), intent(out) :: Foo(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: Fvv(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: Fov(nfuc,nfuc,N_cell)



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

      Foo=cmplx(0.0_dp,0.0_dp,kind=dpc)
      Fvv=cmplx(0.0_dp,0.0_dp,kind=dpc)
      Fov=cmplx(0.0_dp,0.0_dp,kind=dpc)

      ! Foo(m,i,k) = f_mi + sum_nef (2(me|nf)-(mf|ne))*tau_in^ef
      do k=1,N_cell
        do ii=1,Nocc(k)
            i=kOcc(ii,k)
          do mm=1,nOcc(k)
              m=kOcc(mm,k)
              if (m==i) Foo(m,i,k)=cmplx(blk(k)%energy(i),0.0_dp,kind=dpc)
              do kn=1,N_cell
                do ke=1,N_cell
                  kf=ktab(k,kn,ke)
                  do nn=1,nOcc(kn)
                    n=kOcc(nn,kn)
                    do ee=1,nVir(ke)
                      e=kVir(ee,ke)
                      do ff=1,nVir(kf)
                        f=kVir(ff,kf)
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

      ! Fvv(a,e,k) = f_ae - sum_mnf (2(me|nf)-(mf|ne))*tau_mn^af

      do k=1,N_cell
        do ee=1,nVir(k)
            e = kVir(ee,k)
          do aa=1,nVir(k)
            a = kVir(aa,k)
            if (a==e) Fvv(a,e,k)=cmplx(blk(k)%energy(a),0.0_dp,kind=dpc)
            do km=1,N_cell
              do kn=1,N_cell
                kf=ktab(km,kn,k)
                do mm=1,nOcc(km)
                  m = kOcc(mm,km)
                  do nn=1,nOcc(kn)
                    n = kOcc(nn,kn)
                    do ff=1,nVir(kf)
                      f = kVir(ff,kf)
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
        do ee=1,nVir(k)
            e = kVir(ee,k)
          do mm=1,nOcc(k)
            m = kOcc(mm,k)
            do kn=1,N_cell
              do nn=1,nOcc(kn)
                n = kOcc(nn,kn)
                do ff=1,nVir(kn)
                  f = kVir(ff,kn)
                  v=2.0_dp*g(m,e,n,f,k,k,kn)-g(m,f,n,e,k,kn,kn)
                  Fov(m,e,k)=Fov(m,e,k)+v*t1(n,f,kn)
                end do
              end do
            end do
          end do
        end do
      end do

end subroutine build_ccsd_f