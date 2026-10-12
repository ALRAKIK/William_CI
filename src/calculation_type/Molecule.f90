subroutine molecule(n_atoms,number_of_functions,atoms,geometry,OV,K,NA,ERI)

      use files
      use torus_init
      use atom_basis
      use keywords

      implicit none

      ! - input - !
       
      integer          ,intent(in) :: n_atoms
      integer          ,intent(in) :: number_of_functions
      double precision ,intent(in) :: geometry(n_atoms,3)
      type(atom)       ,intent(in) :: atoms(n_atoms)

      ! - local - !

      double precision             :: start,end,time

      ! - output - ! 

      double precision,intent(out)   ::  OV(number_of_functions,number_of_functions)
      double precision,intent(out)   ::   K(number_of_functions,number_of_functions)
      double precision,intent(out)   ::  NA(number_of_functions,number_of_functions)
      double precision,intent(out)   :: ERI(number_of_functions,number_of_functions,number_of_functions,number_of_functions)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


      call header(" OBC Calculation ",-1)
      call header_under(" Calculate the integerals ",-1)

      call cpu_time(start)

        call overlap_matrix             (n_atoms,number_of_functions,geometry,atoms,OV)
        call check_the_overlap          (number_of_functions,OV)
        call kinetic_matrix             (n_atoms,number_of_functions,geometry,atoms,K)
        call nuclear_attraction_matrix  (n_atoms,number_of_functions,geometry,atoms,NA)  
        if (.not. c_One) call ERI_integral               (n_atoms,number_of_functions,geometry,atoms,ERI)

      call cpu_time(end)


      time = end - start

      write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for integrals = ',time,' seconds'
      write(outfile,*)

      call system("rm torus_parameters.inp")

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


end subroutine molecule 