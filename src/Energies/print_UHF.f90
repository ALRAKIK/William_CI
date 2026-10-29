subroutine print_UHF(nBas,n_alpha,n_beta,e_alpha,e_beta,C_alpha,C_beta,ENuc,ET,EV,EJ,Ex,EHF)

      ! Print one- and two-electron energies and other stuff for RHF calculation

      use files 
      implicit none

      integer,intent(in)                 :: nBas,n_alpha , n_beta
      double precision,intent(in)        :: e_alpha(nBas),c_alpha(nBas,nBas)
      double precision,intent(in)        :: e_beta(nBas),c_beta(nBas,nBas)
      double precision,intent(in)        :: ENuc,ET,EV,EJ,Ex,EHF

      double precision                   :: HOMO,LUMO
      double precision                   :: Gap , Gap_alpha , Gap_beta
      double precision,parameter         :: autoev = 27.21138602d0


      Gap_alpha = 0.d0 
      Gap_beta  = 0.d0

      ! HOMO and LUMO
      if (n_alpha > 0 ) then 
        if(nBas > n_alpha) then
          HOMO = e_alpha(n_alpha+1)
          LUMO = e_alpha(n_alpha)
          Gap_alpha = HOMO - LUMO
        else
          Gap_alpha = 0d0
        endif
      end if 

      if (n_beta > 0 ) then 
        if(nBas > n_beta) then
          HOMO = e_beta(n_beta+1)
          LUMO = e_beta(n_beta)
          Gap_beta = HOMO - LUMO
        else
          Gap_beta = 0d0
        endif
      end if 
      
      if (n_alpha > 0 .and. n_beta > 0) then 
        HOMO =  min(e_beta(n_beta+1), e_alpha(n_alpha+1))
        LUMO =  min(e_beta(n_beta), e_alpha(n_alpha))
        Gap = min(Gap_alpha, Gap_beta)
      elseif (n_alpha <= 0 ) then
        HOMO = e_beta(n_beta+1)
        LUMO = e_beta(n_beta)
        Gap = Gap_beta
      elseif (n_beta <= 0 ) then
        HOMO = e_alpha(n_alpha+1)
        LUMO = e_alpha(n_alpha)
        Gap = Gap_alpha
      end if 



      ! Dump results

      write(outfile,*)
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A32)')           ' Summary              '
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A32,1X,F16.10)') ' One-electron energy  ',ET + EV
      write(outfile,'(A32,1X,F16.10)') ' Kinetic      energy  ',ET
      write(outfile,'(A32,1X,F16.10)') ' Potential    energy  ',EV
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A32,1X,F16.10)') ' Two-electron energy  ',EJ + Ex
      write(outfile,'(A32,1X,F16.10)') ' Coulomb      energy  ',EJ
      write(outfile,'(A32,1X,F16.10)') ' Exchange     energy  ',Ex
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A32,1X,F16.10)') ' Electronic   energy  ',EHF
      write(outfile,'(A32,1X,F16.10)') ' Nuclear   repulsion  ',ENuc
      write(outfile,'(A32,1X,F16.10)') ' Hartree-Fock energy  ',EHF + ENuc
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,F13.6)')     ' HF HOMO      energy (au):',HOMO
      write(outfile,'(A36,F13.6)')     ' HF LUMO      energy (au):',LUMO
      write(outfile,'(A36,F13.6)')     ' HF HOMO-LUMO gap    (au):',Gap
      write(outfile,'(A36,F13.6)')     ' HF HOMO-LUMO gap    (eV):',Gap*autoev
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,F13.6)')     ' Potential energy        :',(EV+EJ+Ex+ENuc)
      write(outfile,'(A36,F13.6)')     ' Kinetic   energy        :', ET
      write(outfile,'(A36,F13.6)')     ' Virial    theorem (-V/T):',-(EV+EJ+Ex+ENuc)/ET
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,*)

      ! Print results

      write(outfile,'(A50)') '--------------------------------------------'
      write(outfile,'(A50)') 'Hartree-Fock orbital coefficients  (alpha)  '
      write(outfile,'(A50)') '--------------------------------------------'
      call matout(nBas,nBas,C_alpha)
      write(outfile,*)
      write(outfile,'(A50)') '--------------------------------------------'
      write(outfile,'(A50)') ' Hartree-Fock orbital energies     (alpha)  '
      write(outfile,'(A50)') '--------------------------------------------'
      call matout(nBas,1,e_alpha)
      write(outfile,*)

      write(outfile,'(A50)') '--------------------------------------------'
      write(outfile,'(A50)') 'Hartree-Fock orbital coefficients  (beta)  '
      write(outfile,'(A50)') '--------------------------------------------'
      call matout(nBas,nBas,C_beta)
      write(outfile,*)
      write(outfile,'(A50)') '--------------------------------------------'
      write(outfile,'(A50)') ' Hartree-Fock orbital energies     (beta)  '
      write(outfile,'(A50)') '--------------------------------------------'
      call matout(nBas,1,e_beta)
      write(outfile,*)

end subroutine print_UHF