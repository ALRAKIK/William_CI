!---------------------------------------------------------------------!
! [Nov 2026] keyword frozen_core (fast_ccsdt / fast_ccsdt_real)       !
!---------------------------------------------------------------------!

! number of core orbitals (doubly occupied) of an atom of charge Z:
! H-He 0, Li-Ne 1 (1s), Na-Ar 5 (1s2s2p), K-Kr 9, Rb-Xe 18, Cs-Rn 27, beyond 43
integer function n_core_orbitals(Z)
      implicit none
      integer, intent(in) :: Z
      if (Z <= 2) then
        n_core_orbitals = 0
      else if (Z <= 10) then
        n_core_orbitals = 1
      else if (Z <= 18) then
        n_core_orbitals = 5
      else if (Z <= 36) then
        n_core_orbitals = 9
      else if (Z <= 54) then
        n_core_orbitals = 18
      else if (Z <= 86) then
        n_core_orbitals = 27
      else
        n_core_orbitals = 43
      end if
end function n_core_orbitals

! compact k-space integrals g(p,q,r,s,islot) -> active bands only
! (bands nfc+1..nf of every index); g is replaced by the smaller array
subroutine reduce_kmo_to_active(nf, nfc, ns, g)
      use constants_module, only: dpc
      implicit none
      integer     , intent(in)                 :: nf, nfc, ns
      complex(dpc), allocatable, intent(inout) :: g(:,:,:,:,:)
      complex(dpc), allocatable :: ga(:,:,:,:,:)
      integer :: na, is, s, it
      na = nf - nfc
      allocate(ga(na,na,na,na,ns))
      !$omp parallel do private(is,s)
      do it = 0, ns*na-1
        is = it/na + 1
        s  = mod(it,na) + 1
        ga(:,:,:,s,is) = g(nfc+1:nf, nfc+1:nf, nfc+1:nf, nfc+s, is)
      end do
      !$omp end parallel do
      deallocate(g)
      call move_alloc(ga, g)
end subroutine reduce_kmo_to_active

!---------------------------------------------------------------------!
! [Nov 2026] keyword CC_to_disk: frozen core for the in-memory part    !
! g(1:n1,nf,nf,nf,ns) (occupied first index): keep bands nfc+1.. .     !
!---------------------------------------------------------------------!
subroutine reduce_go_active(n1, nf, nfc, ns, g)
      use constants_module, only: dpc
      implicit none
      integer, intent(in) :: n1, nf, nfc, ns
      complex(dpc), allocatable, intent(inout) :: g(:,:,:,:,:)
      complex(dpc), allocatable :: y(:,:,:,:,:)
      integer :: is
      allocate(y(n1-nfc, nf-nfc, nf-nfc, nf-nfc, ns))
      do is = 1, ns
        y(:,:,:,:,is) = g(nfc+1:n1, nfc+1:nf, nfc+1:nf, nfc+1:nf, is)
      end do
      call move_alloc(y, g)
end subroutine reduce_go_active
