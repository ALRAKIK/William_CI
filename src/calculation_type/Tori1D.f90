subroutine Tori1D(n_atoms,number_of_functions,atoms,AO,geometry,OV,K,NA,ERI)

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

      ! ----------------------    Time     ---------------------------- !
      integer                        :: days, hours, minutes, seconds , time 
      double precision               :: start,end
      !-----------------------------------------------------------------!

      ! - output - ! 

      double precision,intent(out)   ::  OV(number_of_functions,number_of_functions)
      double precision,intent(out)   ::   K(number_of_functions,number_of_functions)
      double precision,intent(out)   ::  NA(number_of_functions,number_of_functions)
      double precision,intent(out)   :: ERI(number_of_functions,number_of_functions,number_of_functions,number_of_functions)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


      call header("Toroidal Gaussian",-1)
      call header_under("Calculate the integerals",-1)
      call Torus_def()
      write(outfile,*) ""
      write(outfile,'(a,f12.8)') "The length of the box:  Lx = ", Lx
      write(outfile,*) ""

      call cpu_time(start)

      call overlap_matrix_toroidal(n_atoms,number_of_functions,atoms,AO,OV)
      call check_the_overlap(number_of_functions,OV)
      call kinetic_matrix_toroidal(n_atoms,number_of_functions,atoms,AO,K)
      call nuclear_attraction_matrix_toroidal_1D(n_atoms,number_of_functions,geometry,atoms,AO,NA)
      if (c_One) then 
        ERI = 0.d0
      else 
        call ERI_integral_toroidal_1D(n_atoms,geometry,number_of_functions,atoms,ERI) 
      end if
        
      call cpu_time(end)

      write(outfile,*) ""

      time = int(end - start)
      days = (time/86400)
      hours=mod(time,86400)/3600
      minutes=mod(mod(time,86400),3600)/60
      seconds=mod(mod(mod(time,86400),3600),60)

      write(outfile,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') 'Total CPU time for integrals = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec"
      write(outfile,*)
      flush(outfile)

      call system("rm torus_parameters.inp")

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


end subroutine Tori1D