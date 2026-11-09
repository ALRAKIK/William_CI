subroutine Tori1D(n_atoms,nBas,atoms,AO,geometry,OV,K,NA,ERI)

      use files
      use torus_init
      use atom_basis
      use classification_ERI
      use keywords
      use unitcell_module

      implicit none

      ! - input - !
       
      integer           ,intent(in)  :: n_atoms , nBas
      type(atom)        ,intent(in)  :: atoms(n_atoms)
      double precision  ,intent(in)  :: geometry(n_atoms,3)
      type(ERI_function),intent(in)  :: AO (nBas)

      ! - local - !

      ! ----------------------    Time     ---------------------------- !
      integer                        :: days, hours, minutes, seconds 
      integer                        :: time 
      double precision               :: start,end
      !-----------------------------------------------------------------!

      ! - output - ! 

      double precision,intent(out)   ::                OV(nBas,nBas)
      double precision,intent(out)   ::                 K(nBas,nBas)
      double precision,intent(out)   ::                NA(nBas,nBas)
      double precision,intent(out)   :: ERI(ERI_size,nBas,nBas,nBas)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


      call header("Toroidal Gaussian",-1)
      call header_under("Calculate the integerals",-1)
      call Torus_def()
      write(outfile,*) ""
      write(outfile,'(a,f12.8)') "The length of the box:  Lx = ", Lx
      write(outfile,*) ""

      call cpu_time(start)

      call overlap_matrix_toroidal_1D(nBas,AO,OV)
      call check_the_overlap(nBas,OV)
      call kinetic_matrix_toroidal_1D(nBas,AO,K)
      call nuclear_attraction_matrix_toroidal_1D(n_atoms,nBas,geometry,atoms,AO,NA)
      if (c_One) then 
        ERI = 0.d0
      else 
        call ERI_integral_toroidal_1D(n_atoms,geometry,nBas,atoms,ERI) 
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


subroutine Tori1D_sym(n_atoms,nBas,atoms,AO,geometry,OV,K,NA,ERI)

      use files
      use torus_init
      use atom_basis
      use classification_ERI
      use keywords
      use unitcell_module

      implicit none

      ! - input - !
       
      integer           ,intent(in)  :: n_atoms , nBas
      type(atom)        ,intent(in)  :: atoms(n_atoms)
      double precision  ,intent(in)  :: geometry(n_atoms,3)
      type(ERI_function),intent(in)  :: AO (nBas)

      ! - local - !

      ! ----------------------    Time     ---------------------------- !
      integer                        :: days, hours, minutes, seconds 
      integer                        :: time 
      double precision               :: start,end
      !-----------------------------------------------------------------!

      ! - output - ! 

      double precision,intent(out)   ::                OV(nBas,nBas)
      double precision,intent(out)   ::                 K(nBas,nBas)
      double precision,intent(out)   ::                NA(nBas,nBas)
      double precision,intent(out)   :: ERI(ERI_size,nBas,nBas,nBas)

      integer                        :: phase(nfuc)

      integer                        :: orb_atom(nfuc)
      integer                        :: orb_local(nfuc)
      double precision               :: atom_pos(nfuc,3)
      logical                        :: has_inv

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


      call header("Toroidal Gaussian",-1)
      call header_under("Calculate the integerals",-1)
      call Torus_def()
      write(outfile,*) ""
      write(outfile,'(a,f12.8)') "The length of the box:  Lx = ", Lx
      write(outfile,*) ""


      !--- group AOs into atoms for the uniqe atoms in the unitcell 

      call build_orb_to_atom(nfuc, AO, orb_atom, orb_local)

      !--- auto-detect the inversion permutation

      if (allocated(inv_perm)) deallocate(inv_perm)
      allocate(inv_perm(nfuc))


      atom_pos(1:nauc,:) = geometry(1:nauc,:)

      call auto_detect_inv_perm(orb_atom, orb_local, atom_pos, has_inv)

      if (.not. has_inv) then
        write(*,'(A)')       'No representable inversion center detected'
        write(outfile,'(A)') 'No representable inversion center detected'
        stop
      end if
      
      call one_electron_inv_sym(nfuc, AO, phase)
        

      call cpu_time(start)

      call overlap_matrix_toroidal_1D_sym(nBas,AO,OV,phase)
      call check_the_overlap(nBas,OV)
      call kinetic_matrix_toroidal_1D_sym(nBas,AO,K,phase)
      call nuclear_attraction_matrix_toroidal_1D_Sym(n_atoms,nBas,geometry,atoms,AO,NA,phase)
        
      if (c_One) then 
        ERI = 0.d0
      else 
        call ERI_integral_toroidal_1D_sym(n_atoms,geometry,nBas,atoms,ERI) 
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


end subroutine Tori1D_sym