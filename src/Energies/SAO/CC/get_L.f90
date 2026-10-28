subroutine build_ccsd_l(blk,g,t1,Foo,Fvv,Loo,Lvv)
      
      use constants_module, only: dp, dpc
      use unitcell_module , only: nfuc, N_cell
      use k_blocks        , only: K_block

      implicit none 

      ! input ! 

      type(K_block), intent(in) :: blk(N_cell)
      complex(dpc) , intent(in) :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc) , intent(in) :: t1(nfuc,nfuc,N_cell)
      complex(dpc) , intent(in) :: Foo(nfuc,nfuc,N_cell),Fvv(nfuc,nfuc,N_cell)

      ! local ! 

      integer                   :: k,kl,kd,m,i,l,c,a,e,d

      ! output ! 

      complex(dpc), intent(out) :: Loo(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: Lvv(nfuc,nfuc,N_cell)
      
      

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