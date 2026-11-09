subroutine sym_call(calculation_type,n_atoms,nBas,atoms,AO,geometry,S,T,V,ERI,ERI_size)

      use files 
      use atom_basis
      use classification_ERI

      implicit none

      character(len=10)  ,intent(in)    ::               calculation_type
      integer            ,intent(in)    ::                        n_atoms
      integer            ,intent(in)    ::                           nBas
      double precision   ,intent(in)    ::            geometry(n_atoms,3)
      type(atom)         ,intent(in)    ::                 atoms(n_atoms)
      type(ERI_function) ,intent(in)    ::                       AO(nBas)

      integer            , intent(in)   ::                       ERI_size

      double precision   , intent(out)  ::                   S(nBas,nBas)
      double precision   , intent(out)  ::                   T(nBas,nBas)
      double precision   , intent(out)  ::                   V(nBas,nBas)
      double precision   , intent(out)  ::   ERI(ERI_size,nBas,nBas,nBas)


      write(   *   ,'(A)')  ''
      write(   *   ,'(A)')  '$ ----------------------------------------------- $'
      write(   *   ,'(A)')  '          Using The inversion symmetry'
      write(   *   ,'(A)')  '$ ----------------------------------------------- $'
      write(   *   ,'(A)')  ''

      write(outfile,'(A)')  ''
      write(outfile,'(A)')  '$ ----------------------------------------------- $'
      write(outfile,'(A)')  '           Using The inversion symmetry'
      write(outfile,'(A)')  '$ ----------------------------------------------- $'
      write(outfile,'(A)')  ''


      select case (trim(calculation_type))
        case ("Tori1D")
          call Tori1D_sym(n_atoms,nBas,atoms,AO,geometry,S,T,V,ERI)           ! Toroidal 1D Gaussian TRR
        case ("Tori2D")
          call Tori2D(n_atoms,nBas,atoms,AO,geometry)                         ! Real Toroidal 2D Gaussian
        case ("Tori3D")
          call Tori3D_sym(n_atoms,nBas,atoms,AO,geometry,S,T,V,ERI)           ! Real Toroidal 3D Gaussian
        case default
          write(outfile,'(A)') 'Unknown calculation type: ',          &
          &                     trim(calculation_type)
          stop
      end select







end subroutine sym_call