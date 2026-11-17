!=====================================================================!
!                                                                     !
!   module cc_dims : compact storage of the k-CCSD(T) arrays          !
!                                                                     !
!   In every k block the occupied bands are the lowest ones           !
!   (energies sorted, occupation = energy <= Fermi level):            !
!       occupied bands  1     .. nOcc(k)   ,  nOcc(k) <= o_hi          !
!       virtual  bands  nOcc(k)+1 .. nfuc  ,  nOcc(k)+1 >= v_lo        !
!   The amplitudes and intermediates are stored only on these ranges  !
!   with custom bounds, e.g.                                          !
!       t2(1:o_hi, 1:o_hi, v_lo:nfuc, v_lo:nfuc, k, k, k)             !
!   so every access t2(i,j,a,b,...) keeps its BAND indices.           !
!                                                                     !
!   o_hi = max_k nOcc(k) ,  v_lo = min_k nOcc(k) + 1                  !
!=====================================================================!

module cc_dims

      implicit none

      integer, save :: o_hi = 0
      integer, save :: v_lo = 1
      ! [Nov 2026] fast_ccsdt: .true. when the orbitals are in the time-reversal
      ! gauge (C(-k) = C(k)^*, see ksym_gauge.f90), so that every k-space array
      ! satisfies X(-k1,-k2,-k3) = X(k1,k2,k3)^* and half of the k blocks can be skipped
      logical, save :: use_ksym = .false.

      ! [Nov 2026] fast_ccsdt: compact storage of the k-space integrals
      !   g(p,q,r,s,k1,k2,k3) is stored as gh(p,q,r,s,islot) with
      !   gidx(k1,k2,k3) =  islot : block stored in gh(:,:,:,:,islot)
      !                  = -islot : block = conjg(gh(:,:,:,:,islot))   (k <-> -k)
      !   With use_ksym only one (k1,k2,k3) of every pair {K,-K} is stored
      !   (the lexicographically smaller one): ~half the memory of g.
      !   Without use_ksym every block has its own slot (nslot = N_cell^3).
      integer, save              :: nslot = 0
      ! [Nov 2026] keyword cholesky: first extent of the compact g. nfuc normally;
      ! o_hi with cholesky (only g(i,q,r,s) with i occupied is stored: every
      ! CC access has an occupied first index except Wvvvv, built from L)
      integer, save              :: g_n1 = 0
      integer, save, allocatable :: gidx(:,:,:)

contains

      ! [Nov 2026] slots of the compact g (see gidx above); fourth_k(1,1,k) = -k
      subroutine set_gslots(N_cell)
        integer, intent(in) :: N_cell
        integer :: k1, k2, k3, m1, m2, m3
        integer :: fourth_k
        if (allocated(gidx)) deallocate(gidx)
        allocate(gidx(N_cell,N_cell,N_cell))
        gidx  = 0
        nslot = 0
        do k3 = 1, N_cell
          do k2 = 1, N_cell
            do k1 = 1, N_cell
              if (gidx(k1,k2,k3) /= 0) cycle
              nslot = nslot + 1
              gidx(k1,k2,k3) = nslot
              if (use_ksym) then
                m1 = fourth_k(1,1,k1); m2 = fourth_k(1,1,k2); m3 = fourth_k(1,1,k3)
                if (gidx(m1,m2,m3) == 0) gidx(m1,m2,m3) = -nslot
              end if
            end do
          end do
        end do
      end subroutine set_gslots

      subroutine set_cc_dims(nfuc, N_cell, occupation)

        integer, intent(in) :: nfuc, N_cell
        integer, intent(in) :: occupation(nfuc, N_cell)   ! 2 occupied, 0 virtual
        integer :: k, p, no_k, omin

        o_hi = 0 ; omin = nfuc
        do k = 1, N_cell
          no_k = count(occupation(:,k) == 2)
          do p = 1, nfuc
            if (p <= no_k .and. occupation(p,k) /= 2) &
              error stop 'cc_dims: occupied bands are not the lowest ones in a k block'
            if (p >  no_k .and. occupation(p,k) /= 0) &
              error stop 'cc_dims: occupation must be 0 or 2'
          end do
          o_hi = max(o_hi, no_k)
          omin = min(omin, no_k)
        end do
        v_lo = omin + 1

      end subroutine set_cc_dims

end module cc_dims
