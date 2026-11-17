!=====================================================================!
!                                                                     !
!   module cc_disk  [Nov 2026]  keyword CC_to_disk                    !
!                                                                     !
!   The vvvv k-blocks of the compact g are written to a direct-access !
!   scratch file (one record per stored slot, exact integrals) and    !
!   read back block by block in the particle-particle ladder. In      !
!   memory only g(i,q,r,s,slot) with an OCCUPIED first index is kept  !
!   (o*f^3 instead of f^4 per slot): every other CC term reads only   !
!   those (including the ovvv integrals (ia|bc)).                     !
!                                                                     !
!   Wvvvv(t1) is formed when the block is read (same formula as       !
!   update_vvvv_in_g, which is not used in this mode):                !
!     W(a,c,b,d) = (ac|bd) - sum_l t1(l,b;kb) (ac|ld)                 !
!                          - sum_l t1(l,a;ka) (lc|bd)                 !
!     (lc|bd) = g(l,c,b,d; ka,kc,kb)  ,  (ac|ld) = g(l,d,a,c; kb,kd,ka) !
!                                                                     !
!   File: ./cc_vvvv_scratch.bin (deleted at the end of the CC).       !
!=====================================================================!

module cc_disk

      use constants_module, only: dp, dpc
      implicit none

      private
      public :: disk_on, disk_open, disk_write_slot, disk_close, disk_bytes
      public :: disk_vvvv_gm, disk_vvvv_gm_r, disk_nv

      logical , save :: disk_on = .false.
      integer , save :: disk_unit = -1
      integer , save :: disk_nv   = 0           ! virtual band count of the records
      integer , save :: disk_vlo  = 0           ! first virtual band when written
      integer(8), save :: nrec_written = 0
      character(len=*), parameter :: disk_file = 'cc_vvvv_scratch.bin'

contains

      subroutine disk_open(v_lo, nfuc)
        integer, intent(in) :: v_lo, nfuc
        integer :: rl
        complex(dpc) :: probe
        disk_vlo = v_lo
        disk_nv  = nfuc - v_lo + 1
        inquire(iolength=rl) probe
        open(newunit=disk_unit, file=disk_file, access='direct', form='unformatted', &
             recl=rl*disk_nv**4, status='replace', action='readwrite')
        disk_on = .true.
        nrec_written = 0
      end subroutine disk_open

      ! one stored slot: blk(1:nfuc,1:nfuc,1:nfuc,1:nfuc) -> record islot (vvvv part)
      subroutine disk_write_slot(islot, nfuc, blk)
        integer     , intent(in) :: islot, nfuc
        complex(dpc), intent(in) :: blk(nfuc,nfuc,nfuc,nfuc)
        !$omp critical (cc_disk_io)
        write(disk_unit, rec=islot) blk(disk_vlo:nfuc, disk_vlo:nfuc, disk_vlo:nfuc, disk_vlo:nfuc)
        nrec_written = nrec_written + 1
        !$omp end critical (cc_disk_io)
      end subroutine disk_write_slot

      double precision function disk_bytes()
        disk_bytes = 16.d0*dble(disk_nv)**4*dble(nrec_written)
      end function disk_bytes

      subroutine disk_close()
        if (disk_on) close(disk_unit, status='delete')   ! newunit numbers are negative
        disk_unit = -1
        disk_on = .false.
      end subroutine disk_close

      !-----------------------------------------------------------------!
      ! W(a,c,b,d) for b in b0..b1 from the record of slot (ka,kc,kb)    !
      ! G1(l,c,b,d) = g(l,c,b,d; ka,kc,kb) , H(l,d,a,c) = g(l,d,a,c; kb,kd,ka)
      ! (both already conjugated when the slot is a conjugate block)     !
      !-----------------------------------------------------------------!
      subroutine w_core(ic, nfuc, v_lo, b0, b1, no, Da, Db, G1, H, T)
        integer     , intent(in)  :: ic, nfuc, v_lo, b0, b1, no
        complex(dpc), intent(in)  :: Da(no, v_lo:nfuc), Db(no, v_lo:nfuc)
        complex(dpc), intent(in)  :: G1(no, v_lo:nfuc, b0:b1, v_lo:nfuc)
        complex(dpc), intent(in)  :: H(no, v_lo:nfuc, v_lo:nfuc, v_lo:nfuc)     ! H(l,d,a,c)
        complex(dpc), intent(out) :: T(v_lo:nfuc, v_lo:nfuc, b0:b1, v_lo:nfuc)  ! T(a,c,b,d)
        complex(dpc), allocatable :: R(:,:,:,:), Hd(:,:,:)
        integer :: nv, nbb, d
        complex(dpc), parameter :: one=(1.0_dp,0.0_dp)
        nv  = nfuc - v_lo + 1
        nbb = b1 - b0 + 1
        if (nv /= disk_nv) error stop 'CC_to_disk: virtual band count changed since the integrals were written'
        allocate(R(nv,nv,nv,nv))
        !$omp critical (cc_disk_io)
        read(disk_unit, rec=abs(ic)) R
        !$omp end critical (cc_disk_io)
        if (ic < 0) R = conjg(R)
        ! R(a,c,b,d) = (ac|bd) on the full virtual range -> b chunk
        T = R(:, :, b0-v_lo+1:b1-v_lo+1, :)
        deallocate(R)
        if (no <= 0) return
        ! T(a,(c,b,d)) -= sum_l Da(l,a) G1(l,(c,b,d))
        call zgemm('T','N', nv, nv*nbb*nv, no, -one, Da, no, G1, no, one, T, nv)
        ! T(a,c,b,d) -= sum_l H(l,d,a,c) Db(l,b)   (for every d)
        allocate(Hd(no, nv, nv))
        do d = v_lo, nfuc
          Hd = H(:, d, :, :)                                  ! Hd(l,a,c)
          call zgemm('T','N', nv*nv, nbb, no, -one, Hd, no, Db(1,b0), no, one, T(v_lo,v_lo,b0,d), nv*nv)
        end do
        deallocate(Hd)
      end subroutine w_core

      ! complex: Gm(c,d,a,b-b0+1) = W(a,c,b,d; ka,kc,kb)
      subroutine disk_vvvv_gm(ka, kc, kb, b0, b1, no, v_lo, occ, t1, g, Gm)
        use unitcell_module, only: nfuc, N_cell
        use cc_dims        , only: gidx, nslot, g_n1
        integer     , intent(in)  :: ka, kc, kb, b0, b1, no, v_lo
        logical     , intent(in)  :: occ(nfuc,N_cell)
        complex(dpc), intent(in)  :: t1(nfuc,nfuc,N_cell)
        complex(dpc), intent(in)  :: g(g_n1,nfuc,nfuc,nfuc,nslot)
        complex(dpc), intent(inout) :: Gm(v_lo:nfuc,v_lo:nfuc,v_lo:nfuc,*)
        complex(dpc), allocatable :: Da(:,:), Db(:,:), G1(:,:,:,:), H(:,:,:,:), T(:,:,:,:)
        integer :: kd, i1, i2, a, b, c, l
        integer :: fourth_k
        kd = fourth_k(ka,kb,kc)                          ! kd = ka - kc + kb
        allocate(Da(no,v_lo:nfuc), Db(no,v_lo:nfuc))
        do a = v_lo, nfuc
          do l = 1, no
            Da(l,a) = (0.0_dp,0.0_dp) ; Db(l,a) = (0.0_dp,0.0_dp)
            if (occ(l,ka) .and. .not.occ(a,ka)) Da(l,a) = t1(l,a,ka)
            if (occ(l,kb) .and. .not.occ(a,kb)) Db(l,a) = t1(l,a,kb)
          end do
        end do
        allocate(G1(no,v_lo:nfuc,b0:b1,v_lo:nfuc), H(no,v_lo:nfuc,v_lo:nfuc,v_lo:nfuc))
        i1 = gidx(ka,kc,kb) ; i2 = gidx(kb,kd,ka)
        if (i1 > 0) then
          G1 = g(1:no, v_lo:nfuc, b0:b1, v_lo:nfuc, i1)
        else
          G1 = conjg(g(1:no, v_lo:nfuc, b0:b1, v_lo:nfuc, -i1))
        end if
        if (i2 > 0) then
          H = g(1:no, v_lo:nfuc, v_lo:nfuc, v_lo:nfuc, i2)
        else
          H = conjg(g(1:no, v_lo:nfuc, v_lo:nfuc, v_lo:nfuc, -i2))
        end if
        allocate(T(v_lo:nfuc, v_lo:nfuc, b0:b1, v_lo:nfuc))
        call w_core(i1, nfuc, v_lo, b0, b1, no, Da, Db, G1, H, T)
        !$omp parallel do collapse(2) private(b,a,c)
        do b = b0, b1
          do a = v_lo, nfuc
            do c = v_lo, nfuc
              Gm(c,:,a,b-b0+1) = T(a,c,b,:)
            end do
          end do
        end do
        !$omp end parallel do
        deallocate(Da, Db, G1, H, T)
      end subroutine disk_vvvv_gm

      ! real (fast_ccsdt_real): records are complex with zero imaginary part
      subroutine disk_vvvv_gm_r(ka, kc, kb, b0, b1, no, v_lo, occ, t1, g, Gm)
        use unitcell_module, only: nfuc, N_cell
        use cc_dims        , only: gidx, nslot, g_n1
        integer , intent(in)  :: ka, kc, kb, b0, b1, no, v_lo
        logical , intent(in)  :: occ(nfuc,N_cell)
        real(dp), intent(in)  :: t1(nfuc,nfuc,N_cell)
        real(dp), intent(in)  :: g(g_n1,nfuc,nfuc,nfuc,nslot)
        real(dp), intent(inout) :: Gm(v_lo:nfuc,v_lo:nfuc,v_lo:nfuc,*)
        complex(dpc), allocatable :: Da(:,:), Db(:,:), G1(:,:,:,:), H(:,:,:,:), T(:,:,:,:)
        integer :: kd, i1, i2, a, b, c, l
        integer :: fourth_k
        kd = fourth_k(ka,kb,kc)
        allocate(Da(no,v_lo:nfuc), Db(no,v_lo:nfuc))
        do a = v_lo, nfuc
          do l = 1, no
            Da(l,a) = (0.0_dp,0.0_dp) ; Db(l,a) = (0.0_dp,0.0_dp)
            if (occ(l,ka) .and. .not.occ(a,ka)) Da(l,a) = cmplx(t1(l,a,ka),0.0_dp,kind=dpc)
            if (occ(l,kb) .and. .not.occ(a,kb)) Db(l,a) = cmplx(t1(l,a,kb),0.0_dp,kind=dpc)
          end do
        end do
        allocate(G1(no,v_lo:nfuc,b0:b1,v_lo:nfuc), H(no,v_lo:nfuc,v_lo:nfuc,v_lo:nfuc))
        i1 = gidx(ka,kc,kb) ; i2 = gidx(kb,kd,ka)
        G1 = cmplx(g(1:no, v_lo:nfuc, b0:b1, v_lo:nfuc, abs(i1)), 0.0_dp, kind=dpc)
        H  = cmplx(g(1:no, v_lo:nfuc, v_lo:nfuc, v_lo:nfuc, abs(i2)), 0.0_dp, kind=dpc)
        allocate(T(v_lo:nfuc, v_lo:nfuc, b0:b1, v_lo:nfuc))
        call w_core(i1, nfuc, v_lo, b0, b1, no, Da, Db, G1, H, T)
        !$omp parallel do collapse(2) private(b,a,c)
        do b = b0, b1
          do a = v_lo, nfuc
            do c = v_lo, nfuc
              Gm(c,:,a,b-b0+1) = real(T(a,c,b,:), dp)
            end do
          end do
        end do
        !$omp end parallel do
        deallocate(Da, Db, G1, H, T)
      end subroutine disk_vvvv_gm_r

end module cc_disk
