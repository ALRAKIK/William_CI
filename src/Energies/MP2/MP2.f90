subroutine MP2(nBas,nO,e,ERI,ENuc,EHF)

      use files 
      implicit none

      ! Input variables

      integer,intent(in)            ::                     nBas
      integer,intent(in)            ::                       nO
      double precision,intent(in)   ::                     ENuc
      double precision,intent(in)   ::                      EHF
      double precision,intent(in)   ::                  e(nBas)
      double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

      ! Local variables

      integer                       ::                      i,j
      integer                       ::                      a,b
      double precision              ::                    Dijab
      double precision              ::                      E2d
      double precision              ::                      E2x
      double precision              ::                    EcMP2

      ! Output variables

      call header_under("MP2 correlation energy",-1)

      !--------------------------------!
      ! Compute MP2 correlation energy !
      !--------------------------------!

      E2d   = 0d0
      E2x   = 0d0

      do i=1,nO
        do j=1,nO
          do a=nO+1,nBas
            do b=nO+1,nBas
            
              Dijab = e(a) + e(b) - e(i) - e(j)
            
              ! Second-order ring diagram
            
              E2d   = E2d   - ERI(i,a,j,b)*ERI(i,a,j,b)/Dijab
            
              ! Second-order exchange diagram
            
              E2x   = E2x   - ERI(i,a,j,b)*ERI(i,b,j,a)/Dijab
            
            enddo
          enddo
        enddo
      enddo

      EcMP2   = 2d0*E2d   - E2x

      !------------!
      ! MP2 energy !
      !------------!

      write(outfile,*)
      write(outfile,'(A32)')           '--------------------------'
      write(outfile,'(A32)')           ' MP2 calculation          '
      write(outfile,'(A32)')           '--------------------------'
      write(outfile,'(A32,1X,F24.16)') ' MP2 correlation energy = ',EcMP2
      write(outfile,'(A32,1X,F24.16)') ' Direct part            = ',2d0*E2d
      write(outfile,'(A32,1X,F24.16)') ' Exchange part          = ',-E2x
      write(outfile,'(A32)')           '--------------------------'
      write(outfile,'(A32,1X,F24.16)') ' Hartree Fock    energy = ',EHF + ENuc
      write(outfile,'(A32,1X,F24.16)') ' MP2 electronic  energy = ',EcMP2
      write(outfile,'(A32,1X,F24.16)') ' MP2 total       energy = ',ENuc + EHF + EcMP2
      write(outfile,'(A32)')           '--------------------------'
      write(outfile,*)

end subroutine MP2