subroutine calculation_type_select(calculation_type,n_atoms,nBas,atoms,geometry,AO,S,T,V,ERI,ERI_size)

      use atom_basis
      use classification_ERI
      use files
      use keywords
      use precomputed_bessel

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


      if (c_nquad_128) then 
        write(*,'(A)') '$ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $'
        write(*,'(A)') 'The Numerical integrals will be calculated using 128 point quadrature'
        write(*,'(A)') '$ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $'
        nquad_point = 128 
      else 
        write(*,'(A)') '$ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $'
        write(*,'(A)') 'The Numerical integrals will be calculated using 64 point quadrature'
        write(*,'(A)') '$ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $ - $'
        nquad_point = 64
      end if

      if (c_inv_sym) then
        call sym_call(calculation_type,n_atoms,nBas,atoms,AO,geometry,S,T,V,ERI)
        return 
      end if 
      

      select case (trim(calculation_type))
        case ("OBC", "Ring", "OBC2D")
          call molecule(n_atoms,nBas,atoms,geometry,S,T,V,ERI)            ! Molecule 
        case ("Torus")
          call Torus_PBC(n_atoms,nBas,atoms,AO,geometry)                  ! Torus with PBC 
        case ("Tori1D")
          call Tori1D(n_atoms,nBas,atoms,AO,geometry,S,T,V,ERI)           ! Toroidal 1D Gaussian TRR
        case ("Tori2D")
          call Tori2D(n_atoms,nBas,atoms,AO,geometry)                     ! Real Toroidal 2D Gaussian
        case ("Tori3D")
          call Tori3D(n_atoms,nBas,atoms,AO,geometry,S,T,V,ERI)           ! Real Toroidal 3D Gaussian
        case default
          write(outfile,'(A)') 'Unknown calculation type: ',          &
          &                     trim(calculation_type)
          stop
      end select   


end subroutine calculation_type_select