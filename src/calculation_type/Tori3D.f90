subroutine Tori3D(n_atoms,number_of_functions,atoms,AO,geometry,OV,K,NA,ERI)

      use files
      use torus_init
      use atom_basis
      use classification_ERI
      use keywords

      implicit none

      ! - input - !
       
      integer           ,intent(in)  :: n_atoms , number_of_functions
      type(atom)        ,intent(in)  :: atoms(n_atoms)
      double precision  ,intent(in)  :: geometry(n_atoms,3)
      type(ERI_function),intent(in)  :: AO (number_of_functions)

      ! - local - !

      double precision               :: start,end,time

      ! - output - ! 

      double precision,intent(out)   ::  OV(number_of_functions,number_of_functions)
      double precision,intent(out)   ::   K(number_of_functions,number_of_functions)
      double precision,intent(out)   ::  NA(number_of_functions,number_of_functions)
      double precision,intent(out)   :: ERI(number_of_functions,number_of_functions,number_of_functions,number_of_functions)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


      call header("Toroidal  3D Gaussian",-1)
      call header_under("Calculate the integerals",-1)
      call Torus_def()
      write(outfile,*) ""
      write(outfile,'(a,f16.8,a,f16.8,a,f16.8,a)') "The length of the box:  Lx , Ly , Lz  = (", Lx ," , ", Ly ," , ",Lz, " )"  
      write(outfile,*) ""

      call cpu_time(start)

      call overlap_matrix_toroidal_3D(n_atoms,number_of_functions,atoms,AO,OV)
      call check_the_overlap(number_of_functions,OV)
      call kinetic_matrix_toroidal_3D(n_atoms,number_of_functions,atoms,AO,K)
      call nuclear_attraction_matrix_toroidal_3D(n_atoms,number_of_functions,geometry,atoms,AO,NA)
      if (.not. c_One) call ERI_integral_toroidal_3D(n_atoms,geometry,number_of_functions,atoms,ERI)

      call cpu_time(end)

      write(outfile,*) ""

      time = end - start

      write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for integrals = ',time,' seconds'
      write(outfile,*)

      call system("rm torus_parameters.inp")

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


end subroutine Tori3D