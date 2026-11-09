subroutine print_mulliken_charges(nBas, natom, nao_atom, S, P, Zat)

      use files, only: outfile

      implicit none
      integer         , intent(in) :: nBas, natom, nao_atom
      double precision, intent(in) :: S(nBas,nBas), P(nBas,nBas)
      integer         , intent(in) :: Zat(natom)   

      integer          :: iat, mu, nu, imu
      double precision :: pop(natom), q

      pop(:) = 0d0
      do iat = 1, natom
         do imu = 1, nao_atom
            mu = (iat-1)*nao_atom + imu
            do nu = 1, nBas
               pop(iat) = pop(iat) + P(mu,nu) * S(nu,mu)
            end do
         end do
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