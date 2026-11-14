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

contains

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
