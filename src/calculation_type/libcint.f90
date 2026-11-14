subroutine molecule(n_atoms,number_of_functions,atoms,geometry,OV,K,NA,ERI)

      use files
      use torus_init
      use atom_basis
      use keywords
      use libcint_integrals

      implicit none

      ! - input - !
       
      integer          ,intent(in)   :: n_atoms
      integer          ,intent(in)   :: number_of_functions
      double precision ,intent(in)   :: geometry(n_atoms,3)
      type(atom)       ,intent(in)   :: atoms(n_atoms)

      ! - local - !

      double precision               :: start,end,time

      ! - output - ! 

      double precision,intent(out)   ::  OV(number_of_functions,number_of_functions)
      double precision,intent(out)   ::   K(number_of_functions,number_of_functions)
      double precision,intent(out)   ::  NA(number_of_functions,number_of_functions)
      double precision,intent(out)   :: ERI(number_of_functions,number_of_functions,number_of_functions,number_of_functions)


      integer, allocatable           :: atm(:,:), bas(:,:), offset(:)
      double precision, allocatable  :: env(:)
      integer                        :: natm, nbas, dim
      integer                        :: ia, jb, kc, ld
      double precision, allocatable  :: dn(:)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

      call header(" OBC Calculation ",-1)
      call header_under(" Calculate the integerals ",-1)

      call cpu_time(start)

      call build_libcint_tables(n_atoms, geometry, atoms, natm, nbas, dim, atm, bas, env, offset)

      if (dim /= number_of_functions) then
        write(*,*) "ERROR: libcint dim = ", dim, " vs n_func = ", number_of_functions
        stop
      end if

      call one_electron_libcint(atm, bas, env, offset, natm, nbas, dim, OV, K, NA)

      ! ---- in-memory AO normalization (replaces the old normOBC) ---- !
      ! every Cartesian AO gets <i|i> = 1, same convention as the       !
      ! Clifford branches. dn(ia) = 1/sqrt(OV(ia,ia))                       !
      allocate(dn(dim))
      do ia = 1, dim
        if (OV(ia,ia) <= 0.d0) then
          write(*,*) "ERROR: non-positive self overlap for AO ", ia, OV(ia,ia)
          stop
        end if
        dn(ia) = 1.d0 / sqrt(OV(ia,ia))
      end do
      do jb = 1, dim
        do ia = 1, dim
          OV(ia,jb) = OV(ia,jb) * dn(ia) * dn(jb)
          K (ia,jb) = K (ia,jb) * dn(ia) * dn(jb)
          NA(ia,jb) = NA(ia,jb) * dn(ia) * dn(jb)
        end do
      end do

      call check_the_overlap (number_of_functions,OV)

      if (.not. c_One) then 
        if (c_physicist) then        ! [Nov 2026] was c_MCDE: one flag for the stored order (get_kw.f90)
          call eri_libcint_physicist(atm, bas, env, offset, natm, nbas, dim, ERI)
        else 
          call eri_libcint_chemist(atm, bas, env, offset, natm, nbas, dim, ERI)
        end if 
      end if 

      if (.not. c_One) then
        !$omp parallel do private(ia,jb,kc,ld) collapse(2)
        do ld = 1, dim
          do kc = 1, dim
            do jb = 1, dim
              do ia = 1, dim
                ERI(ia,jb,kc,ld) = ERI(ia,jb,kc,ld) * dn(ia) * dn(jb) * dn(kc) * dn(ld)
              end do
            end do
          end do
        end do
        !$omp end parallel do
      end if
      deallocate(dn)

      call destroy_libcint_tables(atm, bas, env, offset)

      print*, "ERI done"

      call cpu_time(end)

      time = end - start

      write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for integrals = ',time,' seconds'
      write(outfile,*)

      call system("rm torus_parameters.inp")

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


end subroutine molecule