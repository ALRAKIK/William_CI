!=====================================================================!
!                                                                     !
!   Tori3D : integrals for the 3D toroidal (Clifford) Gaussians       !
!            x,y,z periodic, any l (s ... g)                          !
!                                                                     !
!   1) Torus_def                        (Lx, ax, ...)                !
!   2) normalize_AO_toroidal_3D         (local copy of AO)           !
!   3) overlap_matrix_toroidal_3D            -> OV                   !
!      kinetic_matrix_toroidal_3D            -> K                    !
!      nuclear_attraction_matrix_toroidal_3D -> NA                   !
!      ERI_integral_toroidal (ndp = 3)       -> ERI                  !
!                                                                     !
!   files: src/Toroidal/3D/3D_Overlap_toroidal.f90                    !
!          src/Toroidal/3D/3D_Kinetic_toroidal.f90                    !
!          src/Toroidal/3D/3D_Nuclear_attraction_toroidal.f90         !
!          src/Toroidal/ERI_toroidal.f90  (+ ERI_shell_block_toroidal)!
!          src/basis/normalization/normalize_3D.f90                   !
!          src/module/clifford_tools.f90        (shared tools)        !
!=====================================================================!
subroutine Tori3D(n_atoms,nBas,atoms,AO,geometry,OV,K,NA,ERI)

      use files
      use torus_init
      use atom_basis
      use classification_ERI
      use keywords
      use unitcell_module
      use debug_tools
      use omp_lib

      implicit none

      ! - input - !
       
      integer           ,intent(in)   :: n_atoms , nBas
      type(atom)        ,intent(in)   :: atoms(n_atoms)
      double precision  ,intent(in)   :: geometry(n_atoms,3)
      type(ERI_function),intent(in)   :: AO (nBas)

      ! - local - !

      ! ----------------------    Time     ---------------------------- !
      integer                         :: days, hours, minutes, seconds 
      integer                         :: time 
      double precision                :: start,end
      double precision                :: w0 , w1
      !-----------------------------------------------------------------!

      type(ERI_function), allocatable :: AOn(:)

      ! - output - ! 

      double precision,intent(out)    ::                OV(nBas,nBas)
      double precision,intent(out)    ::                 K(nBas,nBas)
      double precision,intent(out)    ::                NA(nBas,nBas)
      double precision,intent(out)    :: ERI(ERI_size,nBas,nBas,nBas)

      
      ! --------------------------------------------------------------- !
      ! --------------------------------------------------------------- !

      call header("Toroidal Gaussian 3D",-1)
      call header_under("Calculate the integerals",-1)
      call Torus_def()
      write(outfile,*) ""
      write(outfile,'(a,f16.8,a,f16.8,a,f16.8,a)') "The length of the box:  Lx , Ly , Lz  = (", Lx ," , ", Ly ," , ",Lz, " )"
      write(outfile,*) ""

      call cpu_time(start)
      w0 = omp_get_wtime()

      ! $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ !

      allocate(AOn(nBas))
      AOn = AO
      call normalize_AO_toroidal_3D(nBas, AOn, .true.)                   ! normalize the basis functions

      ! $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ !


      call overlap_matrix_toroidal_3D(nBas,AOn,OV)
      call check_the_overlap(nBas,OV)
      call kinetic_matrix_toroidal_3D(nBas,AOn,K)
      call nuclear_attraction_matrix_toroidal_3D(n_atoms,nBas,geometry,atoms,AOn,NA)

      if (c_debug) call check_symmetric('S', OV, outfile)
      if (c_debug) call check_symmetric('T', K, outfile)
      if (c_debug) call check_symmetric('V', NA, outfile)

      if (c_One) then
        ERI = 0.d0
      else
        call ERI_integral_toroidal(3,nBas,AOn,ERI)          ! ndp = 3 periodic directions
      end if

      call cpu_time(end)
      w1 = omp_get_wtime()

      write(outfile,*) ""
      time = int(end - start)
      days = (time/86400)
      hours=mod(time,86400)/3600
      minutes=mod(mod(time,86400),3600)/60
      seconds=mod(mod(mod(time,86400),3600),60)

      write(outfile,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') 'Total CPU time for integrals = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec"
      write(outfile,*)
      write(outfile,'(A65,5X,F12.3,a)') 'Wall time for integrals = ', w1 - w0, " s"
      write(outfile,*)
      flush(outfile)

      deallocate(AOn)

      call system("rm torus_parameters.inp")

end subroutine Tori3D
