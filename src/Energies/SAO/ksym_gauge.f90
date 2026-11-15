!---------------------------------------------------------------------!
! [Nov 2026] used only with the keyword  fast_ccsdt                   !
!                                                                     !
! Time-reversal gauge of the SAO orbitals:                            !
!   the AOs are real, so S(-k) = S(k)^*, F(-k) = F(k)^*, e(-k) = e(k) !
!   but the diagonalizer gives every k its own (random) phases, and   !
!   inside degenerate bands even a different mixing.                  !
!   Here the orbitals are re-chosen (no change of any energy):        !
!     * self-conjugate k (k = -k: Gamma, zone edge): real orbitals    !
!       (real generalized eigenproblem F c = S c e, LAPACK dsygv)     !
!     * pairs (k,-k):  C(-k) = C(k)^*  and  e(-k) = e(k)              !
!   Then the k-space integrals and the CC amplitudes satisfy          !
!     X(-k1,-k2,-k3) = X(k1,k2,k3)^*                                  !
!   which the fast CCSD(T) routines use to skip half of the k blocks. !
!                                                                     !
! The density matrix and the HF energy do not change (only phases /   !
! rotations inside degenerate subspaces).                             !
!---------------------------------------------------------------------!
subroutine ksym_gauge(nBas, nfuc, N_cell, S_SAO, F_SAO, coeff_SAO, e, ok)

      use files
      use constants_module

      implicit none

      integer     , intent(in)    :: nBas, nfuc, N_cell
      complex(dpc), intent(in)    :: S_SAO(nBas,nBas), F_SAO(nBas,nBas)
      complex(dpc), intent(inout) :: coeff_SAO(nBas,nBas)
      double precision, intent(inout) :: e(nBas)
      logical     , intent(out)   :: ok

      integer          :: k, km, b1, b2, m1, m2, info, lwork
      integer          :: fourth_k
      double precision :: Fr(nfuc,nfuc), Sr(nfuc,nfuc), w(nfuc), de, dimF
      double precision, allocatable :: work(:)

      ok = .true.
      lwork = max(1, 3*nfuc*nfuc)
      allocate(work(lwork))

      do k = 1, N_cell
        km = fourth_k(1,1,k)                  ! index of -k
        b1 = (k-1)*nfuc + 1  ; b2 = k*nfuc
        m1 = (km-1)*nfuc + 1 ; m2 = km*nfuc

        if (km == k) then
          ! ---- self-conjugate k: F and S are real, take real orbitals ----
          dimF = max(maxval(abs(aimag(F_SAO(b1:b2,b1:b2)))), maxval(abs(aimag(S_SAO(b1:b2,b1:b2)))))
          Fr = real(F_SAO(b1:b2,b1:b2))
          Sr = real(S_SAO(b1:b2,b1:b2))
          call dsygv(1, 'V', 'U', nfuc, Fr, nfuc, Sr, nfuc, w, work, lwork, info)
          de = maxval(abs(w - e(b1:b2)))
          if (info /= 0 .or. dimF > 1.0d-8 .or. de > 1.0d-8) then
            ok = .false.
            write(outfile,'(a,i5,a,i6,2es10.2)') '  ksym_gauge: real gauge failed at k =', k, &
                                                  '  info, max|Im F,S|, max|de| =', info, dimF, de
            cycle
          end if
          coeff_SAO(b1:b2,b1:b2) = cmplx(Fr, 0.0d0, kind=dpc)
          e(b1:b2) = w

        else if (k < km) then
          ! ---- pair (k,-k): copy the conjugate ----
          de = maxval(abs(e(m1:m2) - e(b1:b2)))
          if (de > 1.0d-8) then
            ok = .false.
            write(outfile,'(a,2i5,a,es10.2)') '  ksym_gauge: e(-k) /= e(k) for k,-k =', k, km, '  max|de| =', de
            cycle
          end if
          coeff_SAO(m1:m2,m1:m2) = conjg(coeff_SAO(b1:b2,b1:b2))
          e(m1:m2) = e(b1:b2)
        end if
      end do

      deallocate(work)

      if (ok) then
        write(outfile,'(a)') '  [fast_ccsdt] time-reversal gauge set: C(-k) = C(k)^*, real C at k = -k'
      else
        write(outfile,'(a)') '  [fast_ccsdt] WARNING: time-reversal gauge NOT set, k <-> -k symmetry not used'
      end if
      flush(outfile)

end subroutine ksym_gauge
