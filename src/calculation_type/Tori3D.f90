subroutine Tori3D(n_atoms,nBas,atoms,AO,geometry,OV,K,NA,ERI)

      use files
      use torus_init
      use atom_basis
      use classification_ERI
      use keywords

      implicit none

      ! - input - !
       
      integer           ,intent(in)  :: n_atoms , nBas
      type(atom)        ,intent(in)  :: atoms(n_atoms)
      double precision  ,intent(in)  :: geometry(n_atoms,3)
      type(ERI_function),intent(in)  :: AO (nBas)

      ! - local - !
  
      ! ----------------------    Time     ---------------------------- !
      integer                        :: days, hours, minutes, seconds , time 
      double precision               :: start,end
      !-----------------------------------------------------------------!


      ! - output - ! 

      double precision,intent(out)   ::  OV(nBas,nBas)
      double precision,intent(out)   ::   K(nBas,nBas)
      double precision,intent(out)   ::  NA(nBas,nBas)
      double precision,intent(out)   :: ERI(nBas,nBas,nBas,nBas)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


      call header("Toroidal  3D Gaussian",-1)
      call header_under("Calculate the integerals",-1)
      call Torus_def()
      write(outfile,*) ""
      write(outfile,'(a,f16.8,a,f16.8,a,f16.8,a)') "The length of the box:  Lx , Ly , Lz  = (", Lx ," , ", Ly ," , ",Lz, " )"  
      write(outfile,*) ""

      call cpu_time(start)

      call overlap_matrix_toroidal_3D(nBas,AO,OV)
      call check_the_overlap(nBas,OV)
      call kinetic_matrix_toroidal_3D(nBas,AO,K)
      call nuclear_attraction_matrix_toroidal_3D(n_atoms,nBas,geometry,atoms,AO,NA)
      
      if (c_One) then 
            ERI = 0.d0 
      else 
            call ERI_integral_toroidal_3D(n_atoms,geometry,nBas,atoms,ERI)
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


end subroutine Tori3D