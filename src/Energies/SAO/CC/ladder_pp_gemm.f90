subroutine ladder_pp_gemm(blk,g,ktab,t1,t2,R2)
  use constants_module, only: dp, dpc
      use cc_dims         , only: gidx, nslot      ! [Nov 2026] compact g
  use unitcell_module , only: nfuc, N_cell
  use cc_dims         , only: o_hi, v_lo, use_ksym   ! [Nov 2026] use_ksym
  use k_blocks        , only: K_block
  implicit none
  type(K_block), intent(in)    :: blk(N_cell)
  complex(dpc) , intent(in)    :: g(nfuc,nfuc,nfuc,nfuc,nslot)   ! [Nov 2026] compact g (cc_dims gidx)
  integer      , intent(in)    :: ktab(N_cell,N_cell,N_cell)
  complex(dpc) , intent(in)    :: t1(nfuc,nfuc,N_cell)
  complex(dpc) , intent(in)    :: t2(o_hi,o_hi,v_lo:nfuc,v_lo:nfuc,N_cell,N_cell,N_cell)
  complex(dpc) , intent(inout) :: R2(o_hi,o_hi,v_lo:nfuc,v_lo:nfuc,N_cell,N_cell,N_cell)
  integer :: nv, no2, ka, kb, kc, kd, ki, kj, a, b, i, j, c, d, k, p, nbc, b0, b1, nb
  integer(8) :: per_b
  integer    :: ic
  complex(dpc), allocatable :: Gm(:,:,:,:), TAUm(:,:,:,:,:), RWm(:,:,:,:,:)
  logical :: occ(nfuc,N_cell), partner
  complex(dpc), parameter :: one=(1.0_dp,0.0_dp), zero=(0.0_dp,0.0_dp)
  integer(8), parameter :: max_bytes = 256_8*1024_8*1024_8

  do k=1,N_cell; do p=1,nfuc; occ(p,k)=(blk(k)%occupation(p)==2); end do; end do
  nv  = nfuc - v_lo + 1
  no2 = o_hi*o_hi
  if (nv <= 0 .or. o_hi <= 0) return

  ! chunk of b so that Gm(nv,nv,nv,nbc) <= max_bytes
  per_b = 16_8*int(nv,8)**3
  nbc   = int(max(1_8, min(int(nv,8), max_bytes/per_b)))

  allocate(Gm(v_lo:nfuc,v_lo:nfuc,v_lo:nfuc,nbc))
  allocate(TAUm(o_hi,o_hi,N_cell,v_lo:nfuc,v_lo:nfuc), RWm(o_hi,o_hi,N_cell,v_lo:nfuc,v_lo:nfuc))

  do ka=1,N_cell
    do kb=1,N_cell
      ! [Nov 2026] k <-> -k symmetry: block (-ka,-kb) = conjugate of block (ka,kb)
      if (use_ksym) then
        if (.not.canon2(ka,kb)) cycle
      end if
      RWm = zero
      do kc=1,N_cell
        ! TAUm(i,j,ki,c,d) = t2(i,j,c,d; ki,kj,kc) (+ t1 t1 if kc == ki)
        !$omp parallel do private(ki,kj,kd,c,d,i,j)
        do ki=1,N_cell
          kj = ktab(ka,kb,ki)
          kd = ktab(ka,kb,kc)
          do d=v_lo,nfuc
            do c=v_lo,nfuc
              ! [Nov 2026] only physical (virtual at kc / kd) c,d contribute
              if (occ(c,kc) .or. occ(d,kd)) then
                TAUm(:,:,ki,c,d) = zero
                cycle
              end if
              TAUm(:,:,ki,c,d) = t2(:,:,c,d,ki,kj,kc)
              if (kc==ki) then
                do j=1,o_hi
                  do i=1,o_hi
                    TAUm(i,j,ki,c,d) = TAUm(i,j,ki,c,d) + t1(i,c,ki)*t1(j,d,kj)
                  end do
                end do
              end if
            end do
          end do
        end do
        !$omp end parallel do
        do b0=v_lo,nfuc,nbc
          b1 = min(nfuc, b0+nbc-1)
          nb = b1 - b0 + 1
          ! Gm(c,d,a,b) = g(a,c,b,d; ka,kc,kb)
          ic = gidx(ka,kc,kb)                   ! [Nov 2026] compact g
          !$omp parallel do collapse(2) private(a,b)
          do b=b0,b1
            do a=v_lo,nfuc
              if (ic > 0) then
                Gm(:,:,a,b-b0+1) = g(a,v_lo:nfuc,b,v_lo:nfuc,ic)
              else
                Gm(:,:,a,b-b0+1) = conjg(g(a,v_lo:nfuc,b,v_lo:nfuc,-ic))
              end if
            end do
          end do
          !$omp end parallel do
          call zgemm('N','N', no2*N_cell, nv*nb, nv*nv, one, TAUm, no2*N_cell, &
                     Gm, nv*nv, one, RWm(1,1,1,v_lo,b0), no2*N_cell)
        end do
      end do
      ! scatter into R2 (physical entries only)
      partner = use_ksym .and. .not.(ktab(1,1,ka)==ka .and. ktab(1,1,kb)==kb)
      do ki=1,N_cell
        kj = ktab(ka,kb,ki)
        do b=v_lo,nfuc; if (occ(b,kb)) cycle
          do a=v_lo,nfuc; if (occ(a,ka)) cycle
            do j=1,o_hi; if (.not.occ(j,kj)) cycle
              do i=1,o_hi; if (.not.occ(i,ki)) cycle
                R2(i,j,a,b,ki,kj,ka) = R2(i,j,a,b,ki,kj,ka) + RWm(i,j,ki,a,b)
                if (partner) R2(i,j,a,b,ktab(1,1,ki),ktab(1,1,kj),ktab(1,1,ka)) = &
                             R2(i,j,a,b,ktab(1,1,ki),ktab(1,1,kj),ktab(1,1,ka)) + conjg(RWm(i,j,ki,a,b))
              end do
            end do
          end do
        end do
      end do
    end do
  end do
  deallocate(Gm, TAUm, RWm)
contains
  ! [Nov 2026] (k1,k2) <= (-k1,-k2) in lexicographic order
  logical function canon2(k1,k2)
    integer, intent(in) :: k1,k2
    canon2 = (k1<ktab(1,1,k1)) .or. (k1==ktab(1,1,k1) .and. k2<=ktab(1,1,k2))
  end function
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
end subroutine ladder_pp_gemm
