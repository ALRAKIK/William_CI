subroutine print_RHF_SAO(nBas,nO,e,c,ENuc,ET,EV,EJ,Ex,EHF,blk)

      ! Print one- and two-electron energies and other stuff for RHF calculation

      use files 
      use keywords
      use unitcell_module
      use constants_module
      use k_blocks

      implicit none

      integer,intent(in)                 :: nBas,nO
      double precision,intent(in)        :: e(nBas),ENuc,ET,EV,EJ,Ex,EHF
      complex(dpc),intent(in)            :: c(nBas,nBas)

      double precision                   :: e_sorted(nBas)
      integer                            :: HOMO,LUMO , i , j , b
      double precision                   :: Gap
      double precision,parameter         :: autoev = 27.21138602d0

      type(K_block)                      :: blk(N_cell)

      ! HOMO and LUMO

      HOMO = nO
      LUMO = HOMO + 1
      
      if (nBas > nO) then
        e_sorted = e
        call sort_array(e_sorted, nBas)
        Gap = e_sorted(nO+1) - e_sorted(nO)
      else
        Gap = 0.d0
      endif


      ! Dump results

      write(outfile,*)
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A32)')           ' Summary              '
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,1X,F24.16)')      ' One-electron energy      ',ET + EV
      write(outfile,'(A36,1X,F24.16)')      ' Kinetic      energy      ',ET
      write(outfile,'(A36,1X,F24.16)')      ' Potential    energy      ',EV
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,1X,F24.16)')      ' Two-electron energy      ',EJ + Ex
      write(outfile,'(A36,1X,F24.16)')      ' Coulomb      energy      ',EJ
      write(outfile,'(A36,1X,F24.16)')      ' Exchange     energy      ',Ex
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,1X,F24.16)')      ' Electronic   energy      ',EHF
      write(outfile,'(A36,1X,F24.16)')      ' Nuclear   repulsion      ',ENuc
      write(outfile,'(A36,1X,F24.16)')      ' Hartree-Fock energy      ',EHF + ENuc
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,1X,F24.16)')      ' Potential energy        :',(EV+EJ+Ex+ENuc)
      write(outfile,'(A36,1X,F24.16)')      ' Kinetic   energy        :', ET
      write(outfile,'(A36,1X,F24.16)')      ' Virial    theorem (-V/T):',-(EV+EJ+Ex+ENuc)/ET      
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,1X,F13.6)')       ' HF HOMO      energy (au):',e_sorted(HOMO)
      write(outfile,'(A36,1X,F13.6)')       ' HF LUMO      energy (au):',e_sorted(LUMO)
      write(outfile,'(A36,1X,F13.6)')       ' HF HOMO-LUMO gap    (au):',Gap
      write(outfile,'(A36,1X,F13.6)')       ' HF HOMO-LUMO gap    (eV):',Gap*autoev
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,*)

      ! Print results

      write(outfile,*)
      write(outfile,'(A50)') '---------------------------------------'
      write(outfile,'(A50)') ' Hartree-Fock orbital energies         '
      write(outfile,'(A50)') '---------------------------------------'
      call matout(nBas,1,e)
      write(outfile,*)


      if (c_Orbitals) then 

      write(orbfile,'(A)') "|-----------------------------------------------------------------------------------------|"
      write(orbfile,'(A)') "|                                 The molecular orbitals                                  |"
      write(orbfile,'(A)') "|-----------------------------------------------------------------------------------------|"

      write(orbfile,'(a)') ""
      write(orbfile,'(A)') "|-------------------------------------  Orbital coeff c ----------------------------------|"
      write(orbfile,'(a)') ""

      
      do b = 1 , N_cell 
        write(orbfile,'(A)') ""
        write(orbfile,'(A,I0)') "k-point  ", b
        write(orbfile,'(A)') ""
        write(orbfile,'((10x,1000(i4,20x)))') (j, j=1, blk(b)%dim)
        do i = 1, blk(b)%dim
          write(orbfile,'(1000(F12.7,F12.7))') ( dble(blk(b)%Coeff(i,j)), dimag(blk(b)%Coeff(i,j)), j=1,blk(b)%dim )
        end do 
        write(orbfile,'((7x,1000(f12.7,12x)))') blk(b)%energy 
        write(orbfile,'((14x,1000(I1,23x)))') blk(b)%occupation
      end do

      Flush(orbfile)
    
      end if 

end subroutine print_RHF_SAO