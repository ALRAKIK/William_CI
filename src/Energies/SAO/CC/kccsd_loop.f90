subroutine solve_kccsd(blk,g,t1,t2,converged, max_iter,damping,energy_tol,residual_tol,ndiis,ecorr)

      use files
      use constants_module, only: dp, dpc
      use unitcell_module , only: nfuc, N_cell
      use k_blocks        , only: K_block

      implicit none 

      ! input ! 

      type(K_block), intent(in)      :: blk(N_cell)
      complex(dpc) , intent(in)      :: g(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      integer      , intent(in)      :: max_iter,ndiis
      real(dp)     , intent(in)      :: damping,energy_tol,residual_tol

      ! output ! 

      complex(dpc) , intent(inout)   :: t1(nfuc,nfuc,N_cell)
      complex(dpc) , intent(inout)   :: t2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      real(dp)     , intent(out)     :: ecorr

      ! local ! 
      
      real(dp)                       :: eimag
      logical, intent(out)           :: converged
      
      complex(dpc), allocatable      :: wo(:,:,:,:,:,:,:)
      complex(dpc), allocatable      :: wv(:,:,:,:,:,:,:)
      complex(dpc), allocatable      :: wx(:,:,:,:,:,:,:)
      complex(dpc), allocatable      :: wy(:,:,:,:,:,:,:)
      complex(dpc), allocatable      :: r1(:,:,:)
      complex(dpc), allocatable      :: r2(:,:,:,:,:,:,:)
      complex(dpc), allocatable      :: err(:)
      complex(dpc), allocatable      :: x(:)   ,trial(:)
      complex(dpc), allocatable      :: hx(:,:),he(:,:)

      real(dp)    , allocatable      :: denom(:),coef(:)
      integer     , allocatable      :: ix(:,:),jx(:,:)

      complex(dpc)                   :: fo(nfuc,nfuc,N_cell),fv(nfuc,nfuc,N_cell),fx(nfuc,nfuc,N_cell)
      complex(dpc)                   :: lo(nfuc,nfuc,N_cell),lv(nfuc,nfuc,N_cell)
      complex(dpc)                   :: energy,tau
      integer                        :: n1,n2,n,p,q,iter,h,nh
      integer                        :: i,j,a,b
      integer                        :: k,ki,kj,ka,kb

      real(dp)                       :: previous,de,norm1,norm2,max1,max2
      logical                        :: ok

      integer, external              :: fourth_k

      interface
        subroutine diis_coefficients(errors,c,done)
          use constants_module, only: dp,dpc
          implicit none
          complex(dpc), intent(in) :: errors(:,:)
          real(dp), intent(out)    :: c(:)
          logical, intent(out)     :: done
        end subroutine diis_coefficients
      end interface

    
      if (max_iter<1.or.ndiis<0.or.damping<=0.0_dp.or.damping>1.0_dp.or.energy_tol<=0.or.residual_tol<=0) &
        error stop 'Invalid CCSD solver parameters'
      
      n1=0 
      do k=1,N_cell
        do i=1,nfuc
          do a=1,nfuc
            if (blk(k)%occupation(i)==2.and.blk(k)%occupation(a)==0) n1=n1+1
          end do
        end do
      end do

      n2=0

      do ki=1,N_cell
        do kj=1,N_cell
          do ka=1,N_cell
            kb=fourth_k(ki,kj,ka)
            do i=1,nfuc
                do j=1,nfuc
                  do a=1,nfuc
                    do b=1,nfuc
                      if (blk(ki)%occupation(i)/=2.or.blk(kj)%occupation(j)/=2) cycle
                      if (blk(ka)%occupation(a)/=0.or.blk(kb)%occupation(b)/=0) cycle
                      n2=n2+1
                    end do
                  end do
                end do
            end do
          end do
        end do
      end do

      n=n1+n2

      if (n1==0.or.n2==0) error stop 'CCSD solver requires occupied and virtual orbitals'

      allocate(ix(3,n1),jx(7,n2),denom(n),x(n),err(n),trial(n))

      p=0
      do k=1,N_cell
        do i=1,nfuc
          do a=1,nfuc
            if (blk(k)%occupation(i)/=2.or.blk(k)%occupation(a)/=0) cycle
            p=p+1
            ix(:,p)=[i,a,k]
            denom(p)=blk(k)%energy(i)-blk(k)%energy(a)
          end do
        end do
      end do

      q=0
      do ki=1,N_cell
        do kj=1,N_cell
          do ka=1,N_cell
            kb=fourth_k(ki,kj,ka)
            do i=1,nfuc
              do j=1,nfuc
                do a=1,nfuc
                  do b=1,nfuc
                    if (blk(ki)%occupation(i)/=2.or.blk(kj)%occupation(j)/=2) cycle
                    if (blk(ka)%occupation(a)/=0.or.blk(kb)%occupation(b)/=0) cycle
                    q=q+1
                    jx(:,q)=[i,j,a,b,ki,kj,ka]
                    denom(n1+q)=blk(ki)%energy(i)+blk(kj)%energy(j)-blk(ka)%energy(a)-blk(kb)%energy(b)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

      if (minval(abs(denom))<1.0e-10_dp) error stop 'CCSD near-zero denominator'

      allocate(r1,mold=t1); allocate(r2,mold=t2)
      allocate(wo,mold=t2); allocate(wv,mold=t2)
      allocate(wx,mold=t2); allocate(wy,mold=t2)
      allocate(hx(n,max(1,ndiis)),he(n,max(1,ndiis)),coef(max(1,ndiis)))

      nh=0
      converged=.false.
      previous=0.0_dp

      ! Pack and clear inactive entries, preserving only physical amplitudes.

      call pack(t1,t2,x)
      call unpack(x)


      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 

      write(outfile,*)
      write(outfile,*) repeat('-', 94)
      write(outfile,*) "|",repeat(' ', 38),"CCSD calculation",repeat(' ', 38),"|"
      write(outfile,*) repeat('-', 94)

      write(outfile,*) repeat('-', 94)
      write(outfile,'(1x,a1,a3,1x,a1,a,a1,a,a1,a,a1,a,a1,a,a1,a,a1)') &
      "|","  #","|","   Correlation   ","|", "      dE     ","|","    Max R1   ","|", "    Max R2   ",&
      "|","    RMS R1   ","|","    RMS R2   ","|"
      write(outfile,*) repeat('-', 94)
      flush(outfile)

      write(*,*) repeat('-', 94)
      write(*,'(1x,a1,a3,1x,a1,a,a1,a,a1,a,a1,a,a1,a,a1,a,a1)') &
      "|","  #","|","   Correlation   ","|", "      dE     ","|","    Max R1   ","|", "    Max R2   ",&
      "|","    RMS R1   ","|","    RMS R2   ","|"
      write(*,*) repeat('-', 94)
      flush(6)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !






      do iter=0,max_iter

        call build_ccsd_f(blk,g,t1,t2,fo,fv,fx)
        call build_ccsd_l(blk,g,t1,fo,fv,lo,lv)
        call build_ccsd_w(blk,g,t1,t2,wo,wv,wx,wy)
        call build_ccsd_residuals(blk,g,t1,t2,fo,fv,fx,lo,lv,wo,wv,wx,wy,r1,r2)
        call pack(r1,r2,err)

        max1=maxval(abs(err(:n1)))
        max2=maxval(abs(err(n1+1:)))

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

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
        
        write(outfile,"(1x,a1,i3,1x,a1,f16.12,1x,a1,f12.8,1x,a1,f12.8,1x,a1,f12.8,1x,a1,f12.8,1x,a1,f12.8,1x,a1)")      & 
        "|",iter,"|",ecorr,"|",dE,"|",max1,"|",max2,"|",norm1,"|",norm2,"|"
        write(outfile,*) repeat('-', 94)
        flush(outfile)

        write(*,"(1x,a1,i3,1x,a1,f16.12,1x,a1,f12.8,1x,a1,f12.8,1x,a1,f12.8,1x,a1,f12.8,1x,a1,f12.8,1x,a1)")      & 
        "|",iter,"|",ecorr,"|",dE,"|",max1,"|",max2,"|",norm1,"|",norm2,"|"
        write(*,*) repeat('-', 94)
        flush(6)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


        if (.not.(max1<huge(1.0_dp).and.max2<huge(1.0_dp).and.abs(ecorr)<huge(1.0_dp))) &
          error stop 'CCSD nonfinite residual or energy'

        if (max(max1,max2)<residual_tol.and.(iter==0.or.abs(de)<energy_tol)) then
          converged=.true.; exit
        end if

      if (iter==max_iter) exit

      call pack(t1,t2,x)
      ! DIIS pairs CURRENT amplitudes with their CURRENT raw residuals.
      ! Extrapolate a Jacobi step from the mixed pair. With DIIS disabled
      ! this is exactly x_new=x+damping*R/D.
      if (ndiis>=2) then
        if (nh==ndiis) then
          hx(:,1:ndiis-1)=hx(:,2:ndiis); he(:,1:ndiis-1)=he(:,2:ndiis)
          nh=ndiis-1
        end if
        nh=nh+1; hx(:,nh)=x; he(:,nh)=err
        if (nh>=2.and.iter>=2) then
          call diis_coefficients(he(:,1:nh),coef(1:nh),ok)
          if (ok) then
            trial=cmplx(0.0_dp,0.0_dp,kind=dpc)
            do h=1,nh
              trial=trial+coef(h)*(hx(:,h)+damping*he(:,h)/denom)
            end do
          else
            trial=x+damping*err/denom
            nh=1; hx(:,1)=x; he(:,1)=err
            write(*,'(A)') 'CCSD DIIS reset: ill-conditioned subspace'
          end if
        else
          trial=x+damping*err/denom
        end if
      else
        trial=x+damping*err/denom
      end if

      previous=ecorr
      call unpack(trial)

      end do

      if (converged) then
        write(*,'(A)')
        write(*,'(A)') '*---------------------*'
        write(*,'(A)') '*-*  CCSD CONVERGED *-*'
        write(*,'(A)') '*---------------------*'
        write(*,'(A)')
      else
        write(*,'(A)')
        write(*,'(A)') '*----------------------------------------------*'
        write(*,'(A)') '   CCSD NOT CONVERGED: iteration limit reached  '
        write(*,'(A)') '*----------------------------------------------*'
        write(*,'(A)')
      end if

      deallocate(ix,jx,denom,x,err,trial,hx,he,coef,r1,r2,wo,wv,wx,wy)
    
contains

subroutine pack(a1,a2,v)

      complex(dpc), intent(in)  :: a1(nfuc,nfuc,N_cell)
      complex(dpc), intent(in)  :: a2(nfuc,nfuc,nfuc,nfuc,N_cell,N_cell,N_cell)
      complex(dpc), intent(out) :: v(n)

      do p=1,n1
        v(p)=a1(ix(1,p),ix(2,p),ix(3,p))
      end do
      do q=1,n2
        v(n1+q)=a2(jx(1,q),jx(2,q),jx(3,q),jx(4,q),jx(5,q),jx(6,q),jx(7,q))
      end do

end subroutine


subroutine unpack(v)

      complex(dpc), intent(in) :: v(n)

      t1=cmplx(0.0_dp,0.0_dp,kind=dpc)
      t2=cmplx(0.0_dp,0.0_dp,kind=dpc)

      do p=1,n1
        t1(ix(1,p),ix(2,p),ix(3,p))=v(p)
      end do
      do q=1,n2
        t2(jx(1,q),jx(2,q),jx(3,q),jx(4,q),jx(5,q),jx(6,q),jx(7,q))=v(n1+q)
      end do

end subroutine

end subroutine solve_kccsd