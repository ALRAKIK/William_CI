program CI 

      use files
      use constants_module
      use torus_init
      use atom_basis
      use classification_ERI
      use trexio
      use keywords
      use precomputed_bessel
      use unitcell_module
      use k_blocks
      use cc_dims , only : o_hi , v_lo , set_cc_dims , set_gslots , nslot , use_ksym , g_n1   ! [Nov 2026] compact g
      use chol_store , only : chol_on , chol_build , chol_build_g , chol_reduce_active , chol_free , chol_bytes   ! [Nov 2026] keyword cholesky
      use cc_disk    , only : disk_on , disk_open , disk_close , disk_bytes   ! [Nov 2026] keyword CC_to_disk
      use omp_lib , only : omp_get_wtime
      use df_toroidal, only : df_kMO , df_build_chol   ! [Nov 2026] keyword density_fitting
      use symmetry_clifford, only : sym_run   ! [Nov 2026] keyword symmetry

      implicit none

      !-----------------------------------------------------------------!
      !                         The  variables                          !
      !-----------------------------------------------------------------!


      integer                         ::                  i ,    j
      integer                         ::            n_atoms , nBAS
      integer                         ::           nO , n_electron 

      double precision                ::  geometry_tmp(max_atom,3)
      double precision                ::  hfw0                    ! [Nov 2026] wall clock at HF start
      integer                         ::      charge_tmp(max_atom)
      character*(2)                   ::       label_tmp(max_atom)
      integer                         ::           n_atom_unitcell 
      integer                         ::      number_of_primitives
      integer                         ::          number_of_shells
      character*(10)                  ::              keyword(100)

      double precision  ,allocatable  ::             geometry(:,:)
      integer           ,allocatable  ::                 charge(:)
      type(atom)        ,allocatable  ::                  atoms(:)
      type(ERI_function),allocatable  ::                    AO (:)

      double precision  ,allocatable  ::                    S(:,:)
      double precision  ,allocatable  ::                    T(:,:)
      double precision  ,allocatable  ::                    V(:,:)
      double precision  ,allocatable  ::                   Hc(:,:)
      double precision  ,allocatable  ::                    X(:,:)
      double precision  ,allocatable  ::              ERI(:,:,:,:)
      double precision  ,allocatable  ::           ERI_MO(:,:,:,:)
      double precision  ,allocatable  ::                      e(:)
      double precision  ,allocatable  ::                coeff(:,:)

      complex(dpc)      ,allocatable  ::            coeff_SAO(:,:)
      integer           ,allocatable  ::                occ_tab(:,:)
      integer                         ::     p_lo , p_hi , q_lo , q_hi
      double precision                ::          w_start , w_end
      complex(dpc)      ,allocatable  ::    ERI_kMO(:,:,:,:,:,:,:)
      complex(dpc)      ,allocatable  ::    ERI_kMO_h(:,:,:,:,:)     ! [Nov 2026] fast_ccsdt compact g
      real(dp)          ,allocatable  ::    ERI_kMO_r(:,:,:,:,:)     ! [Nov 2026] fast_ccsdt_real: real compact g
      logical                         ::    real_path                ! [Nov 2026]
      real(dp)                        ::    g_imax, g_rmax           ! [Nov 2026]
      integer                         ::    k_r, fourth_k_r          ! [Nov 2026]
      integer                         ::    fourth_k
      type(K_block)     ,allocatable  ::    blk_save(:)              ! [Nov 2026] frozen core: full bands
      integer                         ::    n_fc, n_act, nfuc_save   ! [Nov 2026] frozen core
      integer                         ::    n_core_orbitals          ! [Nov 2026]
      interface                                                      ! [Nov 2026]
        subroutine reduce_kmo_to_active(nf, nfc, ns, g)
          use constants_module, only: dpc
          integer     , intent(in)                 :: nf, nfc, ns
          complex(dpc), allocatable, intent(inout) :: g(:,:,:,:,:)
        end subroutine reduce_kmo_to_active
        subroutine reduce_go_active(n1, nf, nfc, ns, g)               ! [Nov 2026] CC_to_disk
          use constants_module, only: dpc
          integer, intent(in) :: n1, nf, nfc, ns
          complex(dpc), allocatable, intent(inout) :: g(:,:,:,:,:)
        end subroutine reduce_go_active
      end interface
      type(K_block)     ,allocatable  ::                    blk(:)

      double precision                ::               E_nuc , EHF
      double precision                ::               E_MP2
      double precision                ::            start,end,time

      character(len=10)               ::          calculation_type 
      character(len=2),allocatable    ::                  label(:)

      integer                         ::                   io_stat
      character(len=100)              ::                      line
      integer                         ::          n_alpha , n_beta


      integer                         :: start_count, end_count
      integer                         :: count_rate
      double precision                :: elapsed

      logical                         :: transformation

      !-----------------------------------------------------------------!
      !                        END variables                            !
      !-----------------------------------------------------------------!

      ! *************************************************************** !

      ! --------------------------------------------------------------- !
      !                   build the super molecule                      !
      ! --------------------------------------------------------------- !
      
      call build_super_molecule(keyword,n_atom_unitcell)

      ! --------------------------------------------------------------- !
      !                        Read key words                           !
      ! --------------------------------------------------------------- !

      call get_keywords(keyword)
                        
      ! --------------------------------------------------------------- !
      !            Read number of alpha and beta electrons              !
      ! --------------------------------------------------------------- !

      if (c_UHF) then 
        n_alpha = 0 
        n_beta  = 0 
        open(1,file="unitcell.mol")
        do
        read(1, '(A)',iostat=io_stat) line
        if (io_stat /= 0) then
          close(1)
          error stop 'keyword UHF: the line UHF was not found in unitcell.mol'   ! [Nov 2026] was a silent stop
        end if

        if (trim(line) == 'UHF') then
          read(1, *,iostat=io_stat) n_alpha , n_beta
          ! [Nov 2026] without this line the run went on with 0 electrons
          if (io_stat /= 0 .or. n_alpha < 0 .or. n_beta < 0 .or. n_alpha + n_beta == 0) then
            write(*,'(a)') ' ERROR: keyword UHF must be followed by a line "n_alpha n_beta" (e.g. "5 5")'
            error stop 'keyword UHF: missing or wrong line "n_alpha n_beta" after UHF'
          end if
          exit
        end if

        end do
        close(1)

      end if 

      ! --------------------------------------------------------------- !
      ! *************************************************************** !

      ! --------------------------------------------------------------- !
      ! Read the Geometry, Label and the Charge from the super molecule !
      ! --------------------------------------------------------------- !

      call read_geometry(n_atoms,charge_tmp,geometry_tmp,               &
      &                  calculation_type,label_tmp)                      

      ! --------------------------------------------------------------- !
      !           prepare the folder to write the output files          !
      ! --------------------------------------------------------------- !

      call initialize_ff(calculation_type,n_atom_unitcell,label_tmp,n_atoms)

      ! --------------------------------------------------------------- !
      !       If the calculation is on a torus read the parameters      !  
      ! --------------------------------------------------------------- !
       
      if (   calculation_type == "Tori1D" .or.             &                            
          &  calculation_type == "Tori2D" .or.             & 
          &  calculation_type == "Tori3D")  then 
        
          call Torus_def()

      end if 

      ! --------------------------------------------------------------- !

      allocate(geometry    (n_atoms,3))
      allocate(charge        (n_atoms))
      allocate(label         (n_atoms))
      allocate(atoms         (n_atoms))

      ! --------------------------------------------------------------- !
      !        Convert the geometry to Bohr if in Angstrom              !
      ! --------------------------------------------------------------- !

      do i = 1 , n_atoms
        label (i)       =      label_tmp(i)
        charge(i)       =     charge_tmp(i)
        if (c_Angstrom) then 
          geometry(i,1) = geometry_tmp(i,1) * Ang_par
          geometry(i,2) = geometry_tmp(i,2) * Ang_par
          geometry(i,3) = geometry_tmp(i,3) * Ang_par
        else
          geometry(i,1) = geometry_tmp(i,1)
          geometry(i,2) = geometry_tmp(i,2)
          geometry(i,3) = geometry_tmp(i,3)
        end if
      end do 

      ! --------------------------------------------------------------- !
      !                        Start Printing                           !
      ! --------------------------------------------------------------- !

      Call Title()

      call system_clock(count=start_count, count_rate=count_rate)

      Call Print_the_input_file() 

      call print_basis(n_atoms,charge)

      ! --------------------------------------------------------------- !

      ! --------------------------------------------------------------- !
      !             Print the geometry of the supermolecule             !
      ! --------------------------------------------------------------- !

      Call HEADER ('The Geometry',-1)

      do i = 1 , n_atoms
        write(outfile,"(6x,I2,3f16.8)") charge(i), (geometry(i,j),j=1,3)
      end do 

      write(outfile,*)
      write(outfile,*)
      flush(outfile)

      ! --------------------------------------------------------------- !

      ! --------------------------------------------------------------- !
      !                  build the atomic basis set                     !
      ! --------------------------------------------------------------- !

      call basis(n_atoms,charge,atoms)

      ! [Nov 2026] number of basis functions of every atom (Mulliken analysis)
      allocate(nao_atom(n_atoms))
      do i = 1 , n_atoms
        nao_atom(i) =       atoms(i)%num_s_function + 3 *atoms(i)%num_p_function &
                    + 6  *  atoms(i)%num_d_function + 10*atoms(i)%num_f_function &
                    + 15 *  atoms(i)%num_g_function
      end do

      ! --------------------------------------------------------------- !
      !                 Information from the Basis set                  !
      ! --------------------------------------------------------------- !

      n_electron = 0 
      do i = 1 , n_atoms
        n_electron = n_electron + atoms(i)%charge
      end do 

      nO = n_electron / 2

      nBAS = 0 
      do i = 1 , n_atoms
        nBAS = nBAS &
                  +      atoms(i)%num_s_function &              ! 1
                  + 3  * atoms(i)%num_p_function &              ! 3
                  + 6  * atoms(i)%num_d_function &              ! 6
                  + 10 * atoms(i)%num_f_function &              ! 10
                  + 15 * atoms(i)%num_g_function                ! 15
      end do

      number_of_primitives = 0
      do i = 1 , n_atoms
        number_of_primitives = number_of_primitives &
                  + atoms(i)%num_exponent_s &
                  + atoms(i)%num_exponent_p &
                  + atoms(i)%num_exponent_d &
                  + atoms(i)%num_exponent_f &
                  + atoms(i)%num_exponent_g
      end do

      number_of_shells = 0
      do i = 1 , n_atoms
        number_of_shells = number_of_shells &
                  + atoms(i)%num_s_function &
                  + atoms(i)%num_p_function &
                  + atoms(i)%num_d_function & 
                  + atoms(i)%num_f_function &
                  + atoms(i)%num_g_function
      end do

      ! ---- non-periodic (OBC/Ring/OBC2D): no unit cell, no SAO ---- !
      if (.not. any(calculation_type == ["Tori1D","Tori2D","Tori3D"])) then   
        if (c_SAO) then
          write(outfile,'(6x,A)') 'WARNING: SAO keyword ignored for non-periodic calculation ('//trim(calculation_type)//')'
          write(*,'(A)')          ' WARNING: SAO keyword ignored for non-periodic calculation ('//trim(calculation_type)//')'
          c_SAO = .false.
        end if
        number_of_atom_in_unitcell = n_atoms
      end if

      nfuc = 0 

      do i = 1 , number_of_atom_in_unitcell
        nfuc = nfuc  & 
                  +      atoms(i)%num_s_function &              ! 1
                  + 3  * atoms(i)%num_p_function &              ! 3
                  + 6  * atoms(i)%num_d_function &              ! 6
                  + 10 * atoms(i)%num_f_function &              ! 10
                  + 15 * atoms(i)%num_g_function                ! 15
      end do

      ! --------------------------------------------------------------- !
      !          Allocate the memory for the atomic orbitals            !
      ! --------------------------------------------------------------- !

      allocate(AO(nBas))

      ! --------------------------------------------------------------- !
      !          Classify the atomic orbitals and print them            !
      ! --------------------------------------------------------------- !

      if (any(calculation_type == ["Tori1D","Tori2D","Tori3D"])) then
        N_cell = nBas / nfuc
      else
        N_cell = 1
      end if
           
      call classification    (n_atoms,nBas,geometry,atoms,AO)
            
      ! --------------------------------------------------------------- !
      !                    print  The orbital table                     !
      ! --------------------------------------------------------------- !
      
        call print_shell_table(AO,nBas)      

      ! --------------------------------------------------------------- !
      !           print  Informations from the Basis set                !
      ! --------------------------------------------------------------- !

      call basis_info(nBas,number_of_primitives,number_of_shells,n_atom_unitcell)

      ! --------------------------------------------------------------- !
      !                    Nuclear repulsion energy                     !
      ! --------------------------------------------------------------- !

      call NRE(calculation_type,n_atoms,geometry,atoms,E_nuc)

      ! --------------------------------------------------------------- !
      !             Write the geometry to a TREXIO file                 !
      ! --------------------------------------------------------------- !

      if (c_trexio) then
        call trexio_conv_init(calculation_type,n_atom_unitcell,label_tmp,n_atoms)
        call trexio_conv_global(n_atoms,label,geometry,charge,          &
        &                       E_nuc,n_electron,nBas,   &
        &                       number_of_primitives,number_of_shells)
      end if


      ! --------------------------------------------------------------- !

      ! --------------------------------------------------------------- !
      !        Memory estimate of every step (before allocating)        !
      ! --------------------------------------------------------------- !

      call memory_report(calculation_type,nBas,nfuc,N_cell,nO)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !               ---------------------------------                 !
      !                       Allocate the memory                       !
      !               ---------------------------------                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas))
      
      allocate(X(nBas,nBas),e(nBas),coeff(nBas,nBas))

      allocate(coeff_SAO(nBas,nBas))

      if (c_SAO) then 
        if (c_direct) then
          allocate(ERI(1,1,1,1))            ! [Nov 2026] keyword direct: no ERIs stored
        else
        allocate(ERI (nfuc,nBas,nBas,nBas))
        end if
        allocate(blk(N_cell))
        ERI_size = nfuc
      else 
        allocate(ERI  (nBas,nBas,nBas,nBas))
        ERI_size = nBas 
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !               ---------------------------------                 !
      !                       Plot the gussians                         !
      !               ---------------------------------                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (c_plot) then 
        call plot(n_atoms,geometry,atoms)
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !              ---------------------------------                 !
      !                  calculate the integrals                       !
      !              ---------------------------------                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (c_read) then

        write(outfile,'(A)') 'The integrals will be read from file'
        call read_integrals_from_file(nBas,S,T,V,Hc,ERI,calculation_type)

      else 

        write(outfile,'(A)') 'The integrals will be calculated'
        call calculation_type_select(calculation_type,n_atoms,nBas,atoms,geometry,AO,S,T,V,ERI,ERI_size)

        Hc(:,:) = T(:,:) + V(:,:)

        ! [Nov 2026] keyword symmetry: detect the usable point operations and check
        ! that S, T, V and the ERIs are invariant under them (nothing uses them yet)
        if (c_symmetry) call sym_run(calculation_type, n_atoms, geometry, atoms(:)%charge, nBas, nfuc, &
                                     AO, S, T, V, ERI, ERI_size, .not. c_direct, outfile)

      end if

      ! --------------------------------------------------------------- !
      !          Stop if only the integrals calculation needed          !
      ! --------------------------------------------------------------- !      

      if (c_Integral) then
        write(outfile,'(A)') 'All of The integrals calculated'
        stop 
      end if

      ! --------------------------------------------------------------- !
      !                  write the integrals to outputfile              !
      ! --------------------------------------------------------------- !

      if (c_details .and. .not. c_direct) then   ! [Nov 2026] direct: no ERI to print
        if (c_SAO) then 
          call details_integrals_SAO(nBas,S,T,V,ERI)
        else 
          call details_integrals(nBas,S,T,V,ERI)
        end if 
      end if 

      ! --------------------------------------------------------------- !

      if (c_trexio .and. .not. c_direct) then
        call trexio_conv_integrals(nBas,S,T,V,Hc,ERI)
      end if

      ! --------------------------------------------------------------- !

      !------------------------------------------------------!
      !                                  (-1/2)         t    !
      !  orthogonalization and get  X = S       =  U s U     !
      !                                                      !
      !------------------------------------------------------!

      !CALL HEADER ('The Overlap Matrix',-1)

      !call matout(nBas,nBas,S)

      ! --------------------------------------------------------------- !


      ! --------------------------------------------------------------- ! 
      !        pure 2D toroidal case (not 2D in 3D normal case)         !
      ! --------------------------------------------------------------- !

    
      call get_X_from_overlap(nBAS,S,X)      
       
      ! ---------------------------------------------------------------- !
      !                                                                  !
      !                           HF code start                          !
      !                                                                  !
      ! ---------------------------------------------------------------- ! 
       
      hfw0 = omp_get_wtime()
      if (c_SAO) then 
        call cpu_time(start)
          call Hartree_Fock_SAO(nBas,nO,S,T,V,Hc,ERI,X,E_nuc,EHF,e,coeff,coeff_SAO,n_alpha,n_beta,blk,n_atoms,charge)
          ! [Nov 2026] restart files: the CC amplitudes are only valid with the same
          ! orbitals (phases), so the HF orbitals of the first run are saved and reused
          if (c_checkpoint) call ckpt_orbitals(nBas, nO, coeff_SAO, e, blk)
        call cpu_time(end)
      else 
        call cpu_time(start)
          call Hartree_Fock(nBas,nO,S,T,V,Hc,ERI,X,E_nuc,EHF,e,coeff,n_alpha,n_beta,n_atoms,charge)
        call cpu_time(end)
      end if 
      
      time = end - start

      write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for HF = ',time,' seconds'
      write(outfile,'(A65,1X,F9.3,A8)') 'Wall time for HF = ',omp_get_wtime()-hfw0,' seconds'   ! [Nov 2026]

      write(outfile,*)
      
      if (c_details) then 
        call system("tar -czf " // trim(output_file_name) // ".tar.gz "  // trim(tmp_file_name) )
      end if 

      ! --------------------------------------------------------------- !
      !       ____________________________________________________      !
      !      |____________________________________________________|     !
      !      | __     __   ____   ___ ||  ____    ____     _  __  |     !
      !      ||  |__ |--|_| || |_|   |||_|**|*|__|+|+||___| ||  | |     !
      !      ||==|^^||--| |=||=| |=*=||| |~~|~|  |=|=|| | |~||==| |     !
      !      ||  |##||  | | || | |JRO|||-|  | |==|+|+||-|-|~||__| |     !
      !      ||__|__||__|_|_||_|_|___|||_|__|_|__|_|_||_|_|_||__|_|     !
      !      ||_______________________||__________________________|     !
      !      | _____________________  ||      __   __  _  __    _ |     !
      !      ||=|=|=|=|=|=|=|=|=|=|=| __..\/ |  |_|  ||#||==|  / /|     !
      !      || | | | | | | | | | | |/\ \  \\|++|=|  || ||==| / / |     !
      !      ||_|_|_|_|_|_|_|_|_|_|_/_/\_.___\__|_|__||_||__|/_/__|     !
      !      |____________________ /\~()/()~//\ __________________|     !
      !      | __   __    _  _     \_  (_ .  _/ _    ___     _____|     !
      !      ||~~|_|..|__| || |_ _   \ //\\ /  |=|__|~|~|___| | | |     !
      !      ||--|+|^^|==|1||2| | |__/\ __ /\__| |==|x|x|+|+|=|=|=|     !
      !      ||__|_|__|__|_||_|_| /  \ \  / /  \_|__|_|_|_|_|_|_|_|     !
      !      |_________________ _/    \/\/\/    \_ _______________|     !
      !      | _____   _   __  |/      \../      \|  __   __   ___|     !
      !      ||_____|_| |_|##|_||   |   \/ __|   ||_|==|_|++|_|-|||     !
      !      ||______||=|#|--| |\   \   o    /   /| |  |~|  | | |||     !
      !      ||______||_|_|__|_|_\   \  o   /   /_|_|__|_|__|_|_|||     !
      !      |_________ __________\___\____/___/___________ ______|     !
      !      |__    _  /    ________     ______           /| _ _ _|     !
      !      |\ \  |=|/   //    /| //   /  /  / |        / ||%|%|%|     !
      !      | \/\ |*/  .//____//.//   /__/__/ (_)      /  ||=|=|=|     !
      !    __|  \/\|/   /(____|/ //                    /  /||~|~|~|__   !
      !      |___\_/   /________//   ________         /  / ||_|_|_|     !
      !      |___ /   (|________/   |\_______\       /  /| |______|     !
      !          /                  \|________)     /  / | |            !
      ! --------------------------------------------------------------- !
      !                            Post-HF                              !
      ! --------------------------------------------------------------- !

      !-----------------------------------------------------------------!
      !                  transformation AO --> MO                       !
      !-----------------------------------------------------------------!

      transformation = (c_MP2 .or. c_MCDE) .and. (.not. c_SAO)

      if (transformation) then 

        ! in place: the AO array becomes the MO array (no second N^4 array)
        w_start = omp_get_wtime()
          ! [Nov 2026] the MO integrals are kept in physicist order when MCDE
          ! runs (MCDE reads <pq|rs>), otherwise chemist order (pq|rs) for MP2.
          ! With physicist-ordered AO integrals, swapping indices 2,3 gives chemist.
          call AO_to_MO_ERI_inplace(nBas,coeff,ERI,c_physicist .and. .not. c_MCDE)
          call move_alloc(ERI, ERI_MO)
        w_end = omp_get_wtime()

        !call print_MO_file(nBas,ERI_MO)

        write(outfile,'(A65,1X,F9.3,A8)') 'Wall time for AO to MO transformation = ',w_end-w_start,' seconds'
        write(outfile,*)

      end if

      !-----------------------------------------------------------------!
      !                     MP2 Energy calculation                      !
      !-----------------------------------------------------------------!

      if (c_MP2) then 

        if (c_SAO) then 

          ! occupied / virtual band ranges of the k blocks (module cc_dims)
          allocate(occ_tab(nfuc,N_cell))
          do i = 1 , N_cell
            occ_tab(:,i) = blk(i)%occupation(1:nfuc)
          end do
          call set_cc_dims(nfuc, N_cell, occ_tab)
          deallocate(occ_tab)

          ! [Nov 2026] keyword fast_ccsdt: COMPACT k-space integrals
          !   g(p,q,r,s,islot), one block of every pair {K,-K} (k <-> -k symmetry,
          !   ~half the memory of ERI_kMO), then the fast k-CCSD(T).
          !   (ccsdt / MP2 alone: unchanged path below)
          if (c_fast_ccsdt) then
            call set_gslots(N_cell)
            write(outfile,'(a,i8,a,i8,a,f8.2,a,l2,a)') '  [fast_ccsdt] k-space integrals: ', nslot, &
                 ' of ', N_cell**3, ' k-blocks stored (', 16.d0*dble(nfuc)**4*dble(nslot)/1024.d0**3, &
                 ' GB, k<->-k symmetry =', use_ksym, ')'
            g_n1 = nfuc
            if (c_cholesky) then
              ! [Nov 2026] keyword cholesky: pivoted Cholesky of the integrals per
              ! momentum transfer K (from the AO integrals of the SCF), then only
              ! g(i,q,r,s) with an occupied first index is stored; Wvvvv is built
              ! from the Cholesky vectors inside the ladder.
              if (c_physicist) error stop 'keyword cholesky needs chemist-ordered integrals (no physicist / MCDE)'
              g_n1 = o_hi
              w_start = omp_get_wtime()
                if (c_df) then
                  ! [Nov 2026] density_fitting + ccsdt: DF vectors instead of Cholesky (no 4-index ERIs)
                  call df_build_chol(nBas, coeff_SAO, outfile)
                else
                call chol_build(nBas, ERI, coeff_SAO, cholesky_tol, outfile)
                end if
                deallocate(ERI)                ! the real-space ERI is not needed any more
                chol_on = .true.
                write(outfile,'(a,f8.2,a,f8.2,a)') '  [cholesky] stored k-space integrals (occupied first index): ', &
                     16.d0*dble(g_n1)*dble(nfuc)**3*dble(nslot)/1024.d0**3, ' GB instead of ', &
                     16.d0*dble(nfuc)**4*dble(nslot)/1024.d0**3, ' GB'
                allocate(ERI_kMO_h(g_n1,nfuc,nfuc,nfuc,nslot))
                call chol_build_g(g_n1, ERI_kMO_h)
              w_end = omp_get_wtime()
              write(outfile,'(A65,1X,F9.3,A8)') 'Wall time for Cholesky + k-MO integrals = ',w_end-w_start,' seconds'
            else if (c_cc_disk) then
              ! [Nov 2026] keyword CC_to_disk: exact integrals; only g(i,q,r,s) with
              ! an occupied first index in memory, the vvvv blocks on disk
              g_n1 = o_hi
              allocate(ERI_kMO_h(g_n1,nfuc,nfuc,nfuc,nslot))
              w_start = omp_get_wtime()
                call disk_open(v_lo, nfuc)
                call transform_ERI_AO_to_kMO_disk(nBas, ERI, coeff_SAO, g_n1, ERI_kMO_h)
                deallocate(ERI)                ! the real-space ERI is not needed any more
              w_end = omp_get_wtime()
              write(outfile,'(a,f8.2,a,f9.2,a,f8.2,a)') '  [CC_to_disk] in memory (occupied first index): ', &
                   16.d0*dble(g_n1)*dble(nfuc)**3*dble(nslot)/1024.d0**3, ' GB , vvvv on disk: ', &
                   disk_bytes()/1024.d0**3, ' GB  (all in memory: ', 16.d0*dble(nfuc)**4*dble(nslot)/1024.d0**3, ' GB)'
              write(outfile,'(A65,1X,F9.3,A8)') 'Wall time for AO to k-MO transformation (disk) = ',w_end-w_start,' seconds'
            else
            allocate(ERI_kMO_h(nfuc,nfuc,nfuc,nfuc,nslot))

            w_start = omp_get_wtime()
              call transform_ERI_AO_to_kMO_half(nBas, ERI, coeff_SAO, 1, nfuc, 1, nfuc, ERI_kMO_h)
              deallocate(ERI)                  ! the real-space ERI is not needed any more
            w_end = omp_get_wtime()
            write(outfile,'(A65,1X,F9.3,A8)') 'Wall time for AO to k-MO transformation = ',w_end-w_start,' seconds'
            end if
            write(outfile,*)

            call cpu_time(start)
              call kMP2_calc_half(N_cell, blk, ERI_kMO_h, E_nuc, EHF, E_MP2)
            call cpu_time(end)
            time = end - start
            write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for K-MP2 calculation = ',time,' seconds'
            write(outfile,*)

            ! [Nov 2026] keyword frozen_core: the n_fc lowest bands of every k
            ! (1s of Li-Ne, 1s2s2p of Na-Ar, ...) are not correlated. The k-space
            ! integrals and the bands are reduced to the active bands, nfuc is set
            ! to the number of active bands during the CC, and the unchanged CC
            ! code runs on the active space (Fock diagonal = HF orbital energies,
            ! which already contain the core). MP2 is repeated with frozen core.
            nfuc_save = nfuc
            if (c_frozen_core) then
              n_fc = 0
              do i = 1, n_atoms
                n_fc = n_fc + n_core_orbitals(charge(i))
              end do
              if (mod(n_fc, N_cell) /= 0) then
                write(*,'(a)') ' ERROR frozen_core: core orbitals not the same in every cell'
                stop
              end if
              n_fc = n_fc / N_cell
              n_act = nfuc - n_fc
              do i = 1, N_cell
                if (count(blk(i)%occupation == 2) <= n_fc) then
                  write(*,'(a)') ' ERROR frozen_core: no active occupied band left'
                  stop
                end if
                if (any(blk(i)%occupation(1:n_fc) /= 2)) then
                  write(*,'(a)') ' ERROR frozen_core: the lowest bands are not occupied'
                  stop
                end if
              end do
              write(outfile,'(a,i5,a,i5,a,i5,a)') '  [frozen_core] ', n_fc, ' core bands frozen per k, ', &
                   n_act, ' of ', nfuc, ' bands active'
              write(outfile,'(a,f14.6,a,f14.6)') '  [frozen_core] highest core band energy =', &
                   maxval([(blk(i)%energy(n_fc), i=1,N_cell)]), '   lowest active band energy =', &
                   minval([(blk(i)%energy(n_fc+1), i=1,N_cell)])
              if (n_fc > 0) then
                if (c_cholesky) then
                  ! [Nov 2026] cholesky: drop the core bands from the vectors, g is rebuilt below
                  call chol_reduce_active(n_fc)
                  deallocate(ERI_kMO_h)
                else if (c_cc_disk) then
                  ! [Nov 2026] CC_to_disk: drop the core bands of the in-memory part
                  ! (the virtual bands on disk are unchanged, only renumbered)
                  call reduce_go_active(g_n1, nfuc, n_fc, nslot, ERI_kMO_h)
                  g_n1 = g_n1 - n_fc
                else
                call reduce_kmo_to_active(nfuc, n_fc, nslot, ERI_kMO_h)
                end if
                allocate(blk_save(N_cell))
                do i = 1, N_cell
                  blk_save(i) = blk(i)
                end do
                deallocate(blk)
                allocate(blk(N_cell))
                do i = 1, N_cell
                  blk(i)%dim = n_act
                  allocate(blk(i)%occupation(n_act), blk(i)%energy(n_act), blk(i)%Coeff(1,1))
                  blk(i)%occupation = blk_save(i)%occupation(n_fc+1:nfuc)
                  blk(i)%energy     = blk_save(i)%energy(n_fc+1:nfuc)
                  blk(i)%Coeff      = (0.0d0,0.0d0)
                end do
                nfuc = n_act
                if (c_cholesky) then
                  g_n1 = g_n1 - n_fc
                  allocate(ERI_kMO_h(g_n1,nfuc,nfuc,nfuc,nslot))
                  call chol_build_g(g_n1, ERI_kMO_h)
                else if (.not. c_cc_disk) then
                  g_n1 = nfuc
                end if
                write(outfile,'(a)') '  [frozen_core] k-MP2 with frozen core:'
                call kMP2_calc_half(N_cell, blk, ERI_kMO_h, E_nuc, EHF, E_MP2)
              end if
              flush(outfile)
            end if

            ! [Nov 2026] keyword fast_ccsdt_real: when every k is self-conjugate
            ! (Gamma / zone edge: 1D N=1,2 ; 2x2 ; 2x2x2 ...) and the orbitals are in
            ! the real time-reversal gauge, all integrals and amplitudes are real:
            ! the CC runs in real arithmetic (dgemm, half the memory of g)
            real_path = .false.
            if (c_fast_real) then
              real_path = use_ksym
              do k_r = 1, N_cell
                fourth_k_r = fourth_k(1,1,k_r)
                if (fourth_k_r /= k_r) real_path = .false.
              end do
              if (.not. real_path) then
                write(outfile,'(a)') '  [fast_ccsdt_real] WARNING: not every k is self-conjugate (or no real gauge):'
                write(outfile,'(a)') '                    complex fast_ccsdt is used instead'
              else
                allocate(ERI_kMO_r(g_n1,nfuc,nfuc,nfuc,nslot))
                call eri_to_real(size(ERI_kMO_h,kind=8), ERI_kMO_h, ERI_kMO_r, g_imax, g_rmax)
                write(outfile,'(a,es10.2,a,es10.2)') '  [fast_ccsdt_real] max |Im g| =', g_imax, '   max |g| =', g_rmax
                if (g_imax > 1.0d-8*max(1.0d0,g_rmax)) then
                  write(outfile,'(a)') '  [fast_ccsdt_real] WARNING: integrals are not real, complex fast_ccsdt is used'
                  deallocate(ERI_kMO_r)
                  real_path = .false.
                else
                  deallocate(ERI_kMO_h)
                  write(outfile,'(a,f8.2,a)') '  [fast_ccsdt_real] real arithmetic, k-space integrals: ', &
                       8.d0*dble(g_n1)*dble(nfuc)**3*dble(nslot)/1024.d0**3, ' GB'
                end if
              end if
              flush(outfile)
            end if

            call cpu_time(start)
            if (real_path) then
              call k_ccsd_T_fast_r(blk,ERI_kMO_r,EHF,E_nuc,E_MP2)
            else
              call k_ccsd_T_fast(blk,ERI_kMO_h,EHF,E_nuc,E_MP2)
            end if
            call cpu_time(end)
            time = end - start
            write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for K-CCSD(T) calculation = ',time,' seconds'
            write(outfile,*)
            if (allocated(ERI_kMO_h)) deallocate(ERI_kMO_h)
            if (allocated(ERI_kMO_r)) deallocate(ERI_kMO_r)
            if (chol_on) then                ! [Nov 2026] keyword cholesky
              call chol_free()
              chol_on = .false.
            end if
            if (disk_on) call disk_close()   ! [Nov 2026] CC_to_disk: scratch file deleted
            ! [Nov 2026] frozen core: back to the full bands
            if (allocated(blk_save)) then
              nfuc = nfuc_save
              deallocate(blk)
              call move_alloc(blk_save, blk)
              allocate(occ_tab(nfuc,N_cell))
              do i = 1, N_cell
                occ_tab(:,i) = blk(i)%occupation(1:nfuc)
              end do
              call set_cc_dims(nfuc, N_cell, occ_tab)
              deallocate(occ_tab)
            end if

          else

          ! real-space AO -> k-space MO in one pass (no ERI_SAO stored)
          !   with CCSD(T): all (pq|rs)      ,  MP2 alone: only (ov|ov)
          if (c_ccsdt) then
            p_lo = 1    ; p_hi = nfuc ; q_lo = 1    ; q_hi = nfuc
          else
            p_lo = 1    ; p_hi = o_hi ; q_lo = v_lo ; q_hi = nfuc
          end if
          allocate(ERI_kMO(p_lo:p_hi,q_lo:q_hi,p_lo:p_hi,q_lo:q_hi,N_cell,N_cell,N_cell))

          w_start = omp_get_wtime()
          if (c_df) then
            ! [Nov 2026] keyword density_fitting: (ov|ov) = 1/N sum_Q L^Q L^Q
            call df_kMO(nBas, coeff_SAO, p_lo, p_hi, q_lo, q_hi, ERI_kMO, outfile)
            deallocate(ERI)
          else if (c_direct) then
            ! [Nov 2026] keyword direct: (ov|ov) from integrals computed shell by shell
            call transform_ERI_direct_kMO(nBas, coeff_SAO, p_lo, p_hi, q_lo, q_hi, ERI_kMO, outfile)
            deallocate(ERI)
          else
            call transform_ERI_AO_to_kMO(nBas, ERI, coeff_SAO, p_lo, p_hi, q_lo, q_hi, ERI_kMO)
            deallocate(ERI)                  ! the real-space ERI is not needed any more
          end if
          w_end = omp_get_wtime()

          write(outfile,'(A65,1X,F9.3,A8)') 'Wall time for AO to k-MO transformation = ',w_end-w_start,' seconds'
          write(outfile,*)


          call cpu_time(start)
            call kMP2_calc(N_cell, blk, p_lo, p_hi, q_lo, q_hi, ERI_kMO, E_nuc, EHF,E_MP2)
          call cpu_time(end)

          time = end - start
          write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for K-MP2 calculation = ',time,' seconds'
          write(outfile,*)


          if (c_ccsdt) then
            call cpu_time(start)
              call k_ccsd_T(blk,ERI_kMO,EHF,E_nuc,E_MP2)
            call cpu_time(end)
            time = end - start
            write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for K-CCSD(T) calculation = ',time,' seconds'
            write(outfile,*)
          end if 

          end if   ! [Nov 2026] c_fast_ccsdt

        else 

          call cpu_time(start)
            call MP2(nBas,nO,e,ERI_MO,E_nuc,EHF,c_MCDE)   ! [Nov 2026] MO order: physicist if MCDE
          call cpu_time(end)

          time = end - start
          write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for MP2 calculation = ',time,' seconds'
          write(outfile,*)

        end if 

      end if 

      !-----------------------------------------------------------------!
      !                       MCDE  calculation                         !
      !-----------------------------------------------------------------!

      if (c_MCDE) then
        call  MCDE_3_1(nBas, nO, e, ERI_MO, c_Lanczos, c_Spin_adapt, c_screen)
      end if 
    
      !-----------------------------------------------------------------!
      ! FCI Energy calculation
      !-----------------------------------------------------------------!

      !call FCI(Hc_MO,ERI_MO,nBAS,E_nuc)
        
      ! --------------------------------------------------------------- !
      ! --------------------------------------------------------------- !

      call system_clock(count=end_count)

      elapsed = real(end_count - start_count) / real(count_rate)

      write(outfile, '(6x, a, i2.2, ":", i2.2, ":", i2.2)') &
      "Elapsed time : ", &
      int(elapsed) / 3600, &
      mod(int(elapsed), 3600) / 60, &
      mod(int(elapsed), 60)
      
      close(outfile)

      call system("rm torus_parameters.inp")
      if ( .not. c_details) then 
        call system("rm -r " // trim(tmp_file_name))
      end if 


      ! --------------------------------------------------------------- !
      !            close the TREXIO file and exit                       !
      ! --------------------------------------------------------------- !

      if (c_trexio) then
        call trexio_conv_close()
      end if

      ! --------------------------------------------------------------- !

end program CI