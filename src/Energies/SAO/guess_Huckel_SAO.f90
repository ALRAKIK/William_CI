subroutine guess_Huckel_RHF_SAO(nBas,c_details,nO,HC,X,ENuc,S,T,V,P)

      use files
      use unitcell_module

      implicit none 

      ! --------------------------------------------------------------- !
      integer         , intent(in)  ::                nBas , nO 
      double precision, intent(in)  ::             HC(nBas,nBas)
      double precision, intent(in)  ::              X(nBas,nBas)
      double precision, intent(in)  ::              T(nBas,nBas)
      double precision, intent(in)  ::              V(nBas,nBas)
      double precision, intent(in)  ::              S(nBas,nBas)
      double precision, intent(in)  ::                     ENuc
      logical         , intent(in)  ::                c_details
      double precision,external     ::             trace_matrix
      ! --------------------------------------------------------------- !
      double precision,allocatable  :: cp(:,:), c(:,:) , e(:) , ct(:,:)
      integer                       :: i, o, j 
      double precision              :: ET, EV
      ! --------------------------------------------------------------- !
      double precision, intent(out) :: P(nbas,nbas)
      ! --------------------------------------------------------------- !
      ! extended Huckel parameter
      double precision,parameter    :: alpha = 0.875d0 
      ! --------------------------------------------------------------- !

      allocate(cp(nbas,nbas),c(nbas,nbas),ct(nbas,nbas), e(nbas))

      cp(:,:) = 0.d0 

      do i = 1 , nbas 
        cp(i,i) = Hc(i,i)
        do j = 1 + i , nbas
          cp(i,j) = alpha * S(i,j) * ( Hc(i,i) + Hc(j,j) )
          cp(j,i) = cp(i,j)
        end do 
      end do

      ! [Nov 2026] BUG FIX: the Huckel matrix must be diagonalized in the orthonormal basis,
      ! H' = X^T H X, then C = X C'. Before, the AO matrix itself was diagonalized and its
      ! eigenvectors multiplied by X: wrong orbital energies and order (e.g. a "degenerate level
      ! at the cut" at -41 Eh for the Be chain), worst for near-linearly dependent bases.
      ! Removed vectors (zero columns of the canonical X) are put at the top (+1e4).
      c(:,:)  = matmul(cp(:,:), X(:,:))
      cp(:,:) = matmul(transpose(X(:,:)), c(:,:))
      do i = 1 , nbas
        if (all(X(:,i) == 0.d0)) cp(i,i) = 1.d4
      end do
      call diagonalize_matrix(nbas, cp, e)
      c(:,:)  = matmul(X(:,:), cp(:,:))

      ct(:,:) = transpose(c)

      ! [Nov 2026] occupation with equal fractional filling of a degenerate
      !            level at the HOMO/LUMO cut (keeps P translation invariant)
      call guess_occupy_cut(nBas, nO, e, c, P)

      ! --------------------------------------------------------------- !
      !                          Print                                  !
      ! --------------------------------------------------------------- !

      ET = trace_matrix(nBas,matmul(P,T))
      EV = trace_matrix(nBas,matmul(P,V))

      if (c_details) then 

      write(HFfile,'(3a)') "!" ,repeat('-',44), "!"
      write(HFfile,'(a)') "                 Iter  =   0"
      write(HFfile,'(3a)') "!" ,repeat('-',44), "!"
      write(HFfile,'(a)') ""
      
      call header_HF("Guess F Matrix, F = HC ", -1)

      write(HFfile,'(17x,1000(i3,15x))') (i,i=1,nfuc)
      do i = 1 , size(HC,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (HC(o,i),o=1,nfuc)
      end do 
      write(HFfile,'(a)') ""

      call header_HF("Guess X Matrix, X = S^(-1/2)", -1)
      write(HFfile,'(17x,1000(i3,15x))') (i,i=1,nfuc)
      do i = 1 , size(X,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (X(o,i),o=1,nfuc)
      end do 
      write(HFfile,'(a)') ""

      call header_HF("Guess Density Matrix, P = 2 C C^t", -1)
      write(HFfile,'(17x,1000(i3,15x))') (i,i=1,nfuc)
      do i = 1 , size(P,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (P(o,i),o=1,nfuc)
      end do 
      write(HFfile,'(a)') ""
      
      call header_HF("Guess Orbital Energies", -1)

      do i = 1 , size(e)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  e(i)
      end do 
      
      write(HFfile,'(a)') ""

      write(HFfile,'(a)') ""
      write(HFfile,'(a,f24.16)')   "      The Kinetic   Energy    = ", ET
      write(HFfile,'(a,f24.16)')   "      The Potential Energy    = ", EV
      write(HFfile,'(a,f24.16,a)') "      The Nuclear   Energy    = ", ENuc
      write(HFfile,'(a,f24.16)')   "-----------------------------------"
      write(HFfile,'(a,f24.16)') " The Energy of the guess  = ", ET+EV+ENuc
      write(HFfile,'(a)') ""

      write(HFfile,'(2a)')  repeat('*_',36) , "*"
      write(HFfile,'(a)')   repeat('_',73)

      end if 

      write(outfile,'(a)') ""
      write(outfile,'(a,f24.16)')   "      The Kinetic   Energy    = ", ET
      write(outfile,'(a,f24.16)')   "      The Potential Energy    = ", EV
      write(outfile,'(a,f24.16,a)') "      The Nuclear   Energy    = ", ENuc
      write(outfile,'(a,f24.16)')   "      -----------------------------------"
      write(outfile,'(a,f24.16)')   "      The guess Energy        = ", ET+EV+ENuc
      write(outfile,'(a)') ""


      ! --------------------------------------------------------------- !
      deallocate(cp,c,ct,e)
      ! --------------------------------------------------------------- !


end subroutine guess_Huckel_RHF_SAO

! [Nov 2026] ---------------------------------------------------------- !
! P = sum_i occ_i c_i c_i^T  (eigenvalues e ascending).                  !
! If the level at the cut (index nO) is degenerate with levels above it, !
! the whole degenerate set [lo,hi] gets the same occupation              !
!   occ = 2 (nO-lo+1)/(hi-lo+1),  levels below lo get 2.                 !
! An eigenspace of a translation-invariant matrix is itself invariant,   !
! so P is periodic; picking nO vectors out of it is arbitrary and not.   !
! --------------------------------------------------------------------- !
subroutine guess_occupy_cut(nBas, nO, e, c, P)

      use files

      implicit none

      integer         , intent(in)  :: nBas, nO
      double precision, intent(in)  :: e(nBas), c(nBas,nBas)
      double precision, intent(out) :: P(nBas,nBas)

      double precision, parameter   :: deg_tol = 1.d-8
      double precision              :: tol, occ
      double precision, allocatable :: cw(:,:)
      integer                       :: lo, hi, i

      P(:,:) = 0.d0
      if (nO <= 0) return

      tol = deg_tol * max(1.d0, abs(e(nO)))
      lo = nO
      do while (lo > 1)
        if (abs(e(lo-1) - e(nO)) >= tol) exit
        lo = lo - 1
      end do
      hi = nO
      do while (hi < nBas)
        if (abs(e(hi+1) - e(nO)) >= tol) exit
        hi = hi + 1
      end do

      if (hi == nO) then
        P(:,:) = 2d0 * matmul(c(:,1:nO), transpose(c(:,1:nO)))
        return
      end if

      occ = 2.d0 * dble(nO - lo + 1) / dble(hi - lo + 1)

      allocate(cw(nBas,hi))
      cw(:,1:lo-1) = c(:,1:lo-1) * sqrt(2.d0)
      cw(:,lo:hi)  = c(:,lo:hi)  * sqrt(occ)
      P(:,:) = matmul(cw, transpose(cw))
      deallocate(cw)

      write(outfile,'(a)') ""
      write(outfile,'(a,i0,a,i0,a,f20.12)') &
        "      Guess: degenerate level at the cut, orbitals ", lo, " - ", hi, &
        "   e = ", e(nO)
      write(outfile,'(a,i0,a,i0,a,f12.8)') &
        "      Guess: ", 2*(nO-lo+1), " electrons spread over ", hi-lo+1, &
        " orbitals, occupation = ", occ

end subroutine guess_occupy_cut
