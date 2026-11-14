!---------------------------------------------------------------------!
! [Nov 2026] the number of basis functions of each atom is taken from  !
! nao_atom(:) (module atom_basis, filled in main from the basis set).  !
! Before: nBas/natom functions were assumed for every atom, which is   !
! wrong as soon as the atoms are different (e.g. H2O).                 !
! The functions of atom iat are consecutive (same order as the AO list)!
!---------------------------------------------------------------------!
subroutine print_mulliken_charges(nBas, natom, S, P, Zat)

      use files, only: outfile
      use atom_basis, only: nao_atom

      implicit none
      integer         , intent(in) :: nBas, natom
      double precision, intent(in) :: S(nBas,nBas), P(nBas,nBas)
      integer         , intent(in) :: Zat(natom)   

      integer          :: iat, mu, nu, imu, first
      double precision :: pop(natom), q

      pop(:) = 0d0
      first = 0
      do iat = 1, natom
         do imu = 1, nao_atom(iat)
            mu = first + imu
            do nu = 1, nBas
               pop(iat) = pop(iat) + P(mu,nu) * S(nu,mu)
            end do
         end do
         first = first + nao_atom(iat)
      end do

      write(outfile,'(6x,a)') ""
      write(outfile,'(10x,a)') " Mulliken population analysis"
      write(outfile,'(10x,a)') " ----------------------------"
      write(outfile,'(6x,a)') " Atom       Population        Charge"
      do iat = 1, natom
         q = Zat(iat) - pop(iat)
         write(outfile,'(6x,i3,2x,2f16.6)') iat, pop(iat), q
      end do
      write(outfile,'(6x,a)') ""
      write(outfile,'(6x,a,f16.6)') "Total population : ", sum(pop)
      write(outfile,'(6x,a,f16.6)') "Total charge     : ", sum(Zat) - sum(pop)
      write(outfile,'(6x,a)') ""

end subroutine