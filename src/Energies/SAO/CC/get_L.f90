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
      integer                   :: kOcc(nfuc,N_cell), nOcc(N_cell)
      integer                   :: kVir(nfuc,N_cell), nVir(N_cell)
      integer                   :: ii , mm , ee, aa , ll , cc , mm2 , dd 
      integer                   :: p 

      ! output ! 

      complex(dpc), intent(out) :: Loo(nfuc,nfuc,N_cell)
      complex(dpc), intent(out) :: Lvv(nfuc,nfuc,N_cell)
      
      

      ! Canonical HF assumption: bare f_ov=0, matching build_ccsd_f.
      ! Loo(m,i)=Foo(m,i)+sum_lc [2(lc|mi)-(mc|li)]*t1(l,c)
      ! Lvv(a,e)=Fvv(a,e)+sum_md [2(md|ae)-(me|ad)]*t1(m,d)
      ! The bare f_ov here must NOT be replaced by dressed Fov.

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


      Loo=cmplx(0.0_dp,0.0_dp,kind=dpc)
      Lvv=cmplx(0.0_dp,0.0_dp,kind=dpc)

            ! ---- Loo: m,i occupied at k ; l occupied, c virtual at kl ----
      do k=1,N_cell
        do mm=1,nOcc(k)
          m=kOcc(mm,k)
          do ii=1,nOcc(k)
            i=kOcc(ii,k)
            Loo(m,i,k)=Foo(m,i,k)
            do kl=1,N_cell
              do ll=1,nOcc(kl)
                l=kOcc(ll,kl)
                do cc=1,nVir(kl)
                  c=kVir(cc,kl)
                  Loo(m,i,k)=Loo(m,i,k)+(2.0_dp*g(l,c,m,i,kl,kl,k) &
                    -g(m,c,l,i,k,kl,kl))*t1(l,c,kl)
                end do
              end do
            end do
          end do
        end do
      end do

      ! ---- Lvv: a,e virtual at k ; m occupied, d virtual at kd ----
      do k=1,N_cell
        do aa=1,nVir(k)
          a=kVir(aa,k)
          do ee=1,nVir(k)
            e=kVir(ee,k)
            Lvv(a,e,k)=Fvv(a,e,k)
            do kd=1,N_cell
              do mm2=1,nOcc(kd)
                m=kOcc(mm2,kd)
                do dd=1,nVir(kd)
                  d=kVir(dd,kd)
                  Lvv(a,e,k)=Lvv(a,e,k)+(2.0_dp*g(m,d,a,e,kd,kd,k) &
                    -g(m,e,a,d,kd,k,k))*t1(m,d,kd)
                end do
              end do
            end do
          end do
        end do
      end do

end subroutine build_ccsd_l