subroutine build_ccsd_f(blk,g,t1,t2,Foo,Fvv,Fov)

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

      integer      :: k,km,kn,ke,kf,m,n,i,a,e,f
      integer      :: fourth_k
      complex(dpc) :: v,tau

      ! output !

      complex(dpc), intent(out) :: Foo(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: Fvv(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: Fov(nfuc,nfuc,N_cell)

      Foo=cmplx(0.0_dp,0.0_dp,kind=dpc)
      Fvv=cmplx(0.0_dp,0.0_dp,kind=dpc)
      Fov=cmplx(0.0_dp,0.0_dp,kind=dpc)

      ! Foo(m,i,k) = f_mi + sum_nef (2(me|nf)-(mf|ne))*tau_in^ef

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

      ! Fvv(a,e,k) = f_ae - sum_mnf (2(me|nf)-(mf|ne))*tau_mn^af

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