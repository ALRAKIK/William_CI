subroutine trexio_conv_init(calculation_type,n_atom_unitcell,label_tmp,n_atoms)
   
      use files 
      use trexio
  
      implicit none

      ! input ! 

      character(len=10),intent(in)    :: calculation_type
      integer          ,intent(in)    :: n_atoms
      integer          ,intent(in)    :: n_atom_unitcell
      character(len=2) ,intent(in)    :: label_tmp(n_atom_unitcell)

      ! local ! 

      integer(trexio_exit_code)       :: rc       
      integer                         :: i
      character*(128)                 :: err_msg  
      character(len=100)              :: trexio_file_name
      character(len=100)              :: label


      ! --------------------------------------------------------------- !
      !                 initializing the trexio file                    !
      ! --------------------------------------------------------------- !

      write(label,'(100A)') (trim(label_tmp(i)), i = 1, n_atom_unitcell)


      write(trexio_file_name,'(A,A,A,A,I0,A)') trim(calculation_type),"_",trim(label),"_",n_atoms,".trexio"

      call system("rm -r  "// trim(trexio_file_name))

      trexio_file = trexio_open (trim(trexio_file_name), 'w', TREXIO_HDF5, rc)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
        print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

end subroutine trexio_conv_init

subroutine trexio_conv_close()
   
      use files 
      use trexio
  
      implicit none

      ! local ! 

      integer(trexio_exit_code)       :: rc       
      character*(128)                 :: err_msg  


      ! close the trexio file !

      rc = trexio_close(trexio_file)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
        print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

end subroutine trexio_conv_close


subroutine trexio_conv_global(n_atoms,label,geometry,charge,E_nuc,n_electron,&
                              number_of_functions,number_of_primitives,&
                              number_of_shells)
   
      use files 
      use trexio
  
      implicit none

      ! input ! 

      integer          ,intent(in)    :: n_atoms
      integer          ,intent(in)    :: charge(n_atoms)
      integer          ,intent(in)    :: n_electron
      integer          ,intent(in)    :: number_of_functions
      integer          ,intent(in)    :: number_of_primitives
      integer          ,intent(in)    :: number_of_shells
      character(len=2) ,intent(in)    :: label(n_atoms)
      double precision ,intent(in)    :: geometry(n_atoms,3)
      double precision ,intent(in)    :: E_nuc
      

      ! local ! 

      integer(trexio_exit_code)       :: rc
      character*(128)                 :: err_msg  
      double precision                :: Nshell(number_of_shells)
      double precision                :: NAO(number_of_functions)

      ! --------------------------------------------------------------- !
      !                 writing into the trexio file                    ! 
      ! --------------------------------------------------------------- !

      !         - Writing the global (metadata) information -           !

      rc = trexio_write_metadata_code_num(trexio_file, 1)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_metadata_code(trexio_file,"Clifford Gaussian",20)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_metadata_author_num(trexio_file,1)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_metadata_author(trexio_file,"Amer Alrakik",15)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if
 
      rc = trexio_write_metadata_description(trexio_file,'Calculate the integrals using Clifford Gaussian',49)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      !        - Writing the system (nucleus group) information -       !


      rc = trexio_write_nucleus_num (trexio_file, n_atoms)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_nucleus_charge (trexio_file, real(charge,8))
      if (rc /= TREXIO_SUCCESS) then
         call trexio_string_of_error(rc, err_msg)
         print *, 'Error: '//trim(err_msg)
         call exit(-1)
       end if

      rc = trexio_write_nucleus_coord (trexio_file, geometry)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
        print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_nucleus_label (trexio_file, label,n_atoms)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if 

    
      rc = trexio_write_nucleus_repulsion (trexio_file, E_nuc)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      !       - Writing the Electron (electron group) information -     !

      rc = trexio_write_electron_num (trexio_file, n_electron)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_electron_up_num (trexio_file, n_electron/2)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_electron_dn_num (trexio_file, n_electron/2)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      !       - Writing the Basis set (basis group) information -       !

      rc = trexio_write_basis_type(trexio_file,"Gaussian",8)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_basis_prim_num(trexio_file, number_of_primitives)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_basis_shell_num(trexio_file, number_of_shells)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      Nshell(:) = 1.d0


      rc = trexio_write_basis_shell_factor(trexio_file, Nshell)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if
      
      !          - Writing the Orbitals (AO group) information -        !


      rc = trexio_write_ao_cartesian(trexio_file, 1)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_ao_num(trexio_file, number_of_functions)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if



      NAO(:) = 1.d0


      rc = trexio_write_ao_normalization(trexio_file, NAO)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      

end subroutine trexio_conv_global 

subroutine trexio_conv_integrals(nBas,S,T,V,Hc,ERI)
   
      use files 
      use trexio
  
      implicit none

      ! input ! 
      
      integer          ,intent(in)     :: nBas
      double precision ,intent(in)     :: S(nBas,nBas)
      double precision ,intent(in)     :: T(nBas,nBas)
      double precision ,intent(in)     :: V(nBas,nBas)
      double precision ,intent(in)     :: Hc(nBas,nBas)
      double precision ,intent(in)     :: ERI(nBas,nBas,nBas,nBas)

      ! local ! 

      integer                         :: i , j , k , l

      integer(trexio_exit_code)       :: rc       
      character*(128)                 :: err_msg  
      integer(8)                      :: offset, icount
      integer(8), parameter           :: BUFSIZE=10000_8
      integer                         :: buffer_index(4,BUFSIZE)
      double precision                :: buffer_values(BUFSIZE)
      double precision                :: integral

      ! writing into the trexio file ! 

      rc = trexio_write_ao_1e_int_overlap(trexio_file, S)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_ao_1e_int_kinetic(trexio_file, T)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_ao_1e_int_potential_n_e(trexio_file, V)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_ao_1e_int_core_hamiltonian(trexio_file, Hc)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      icount = 0_8
      offset = 0_8
      do l = 1 , nBas
        do k = 1 , nBas
          do j = 1 , nBas
            do i = 1 , nBas
              if (i==j .and. k<l) cycle
              if (i<j) cycle
              integral = ERI(i,k,j,l)
              if (integral == 0.d0) cycle
              icount = icount + 1_8
              buffer_index(1,icount) = i
              buffer_index(2,icount) = j
              buffer_index(3,icount) = k
              buffer_index(4,icount) = l
              buffer_values(icount)  = integral
              if (icount == BUFSIZE) then
                rc = trexio_write_ao_2e_int_eri(trexio_file,offset,BUFSIZE,buffer_index,buffer_values)
                call trexio_assert(rc, TREXIO_SUCCESS)
                offset = offset + icount
                icount = 0_8
              end if
            end do 
          end do
        end do 
      end do

      if (icount > 0_8) then
      rc = trexio_write_ao_2e_int_eri(trexio_file, offset, BUFSIZE, buffer_index, buffer_values)
      call trexio_assert(rc, TREXIO_SUCCESS)
      end if
      
      

end subroutine trexio_conv_integrals