subroutine build_super_molecule(keyword,num_atom_per_unitcell)

      use constants_module
      use unitcell_module

      implicit none 
      
      !-----------------------------------------------------------------!

      ! local ! 

      integer                        :: num_atoms
      integer                        :: io_stat, i, j, k, l, atom_index

      double precision               :: distance_between_unitcells
      double precision               :: Lx , Ly , Lz
      double precision               :: rx, theta

      double precision,allocatable   ::  geometry_unitcell(:,:)
      double precision,allocatable   ::           unitcell(:,:)
      double precision,allocatable   ::     super_geometry(:,:)

      character(len=10)              ::    type_of_calculation
      character(len=100)             ::                   line

      character(len=2),allocatable   ::          atom_names(:)
      character(len=2),allocatable   ::             a_names(:)
      character(len=2),allocatable   ::         super_atoms(:)

      logical                        ::           reading_cell


      ! output !

      character(len=10),intent(out)  ::            keyword(100)
      integer                        :: num_atom_per_unitcell

      !-----------------------------------------------------------------!

      open(1,file="unitcell.mol")

        read(1,*) type_of_calculation

        call read_2nd_line(type_of_calculation, distance_between_unitcells, Lx , Ly , Lz)

        call read_keywords(keyword)

        allocate(atom_names(max_atom))
        allocate(geometry_unitcell(max_atom,3))

        num_atoms              = 0
        geometry_unitcell(:,3) = 0.d0 

      ! --------------------------------------------------------------- !
      !               Read unitcell geometry from file                  !
      ! --------------------------------------------------------------- !

        do

          read(1, '(A)', iostat=io_stat) line
          if (trim(line) == '$$') then
            reading_cell = .false.
            exit
          end if

          num_atoms = num_atoms + 1

          read(line, *, iostat=io_stat) atom_names(num_atoms),            &
          &                             geometry_unitcell(num_atoms,1),   &
          &                             geometry_unitcell(num_atoms,2),   &
          &                             geometry_unitcell(num_atoms,3)

          if (io_stat /= 0) then
            print *, 'Error parsing atom data:', trim(line)
            num_atoms = num_atoms - 1
          end if

        end do

      close(1)

      nauc = num_atoms

      ! --------------------------------------------------------------- !
      ! --------------------------------------------------------------- !

      num_atom_per_unitcell = num_atoms

      allocate(unitcell(num_atoms,3))
      allocate(a_names(num_atoms))

      do i = 1 , num_atoms
        unitcell(i,1) = geometry_unitcell(i,1)
        unitcell(i,2) = geometry_unitcell(i,2)
        unitcell(i,3) = geometry_unitcell(i,3)
        a_names(i)    = atom_names(i)
      end do

      deallocate(atom_names)
      deallocate(geometry_unitcell)

      allocate(super_atoms(num_atoms * nx * ny * nz))
      allocate(super_geometry(num_atoms * nx * ny * nz, 3))

      atom_index = 0

      ! --------------------------------------------------------------- !

      if (   type_of_calculation == 'OBC'    .or. &
          &  type_of_calculation == 'Tori1D' .or. & 
          &  type_of_calculation == 'Tori2D' .or. & 
          &  type_of_calculation == 'Tori3D')        then 

        ! BUG FIX: cells are generated with x FASTEST (then y, then z),
        ! the order used by fill_with_symmetry (S,T,V), SAO_Transfer, AO_to_MO_SAO,
        ! MP2_SAO and fourth_k. Before, z was fastest (loops i,j,k), so for 2D/3D
        ! cells with nx /= ny (e.g. 3x2) the first-cell rows were copied to the
        ! wrong cells -> wrong S,T,V (negative overlap eigenvalue) and wrong SAO.
        ! 1D (ny = nz = 1), nx = ny 2D cells and Ring are not affected.

        do k = 0 , nz - 1
          do j = 0, ny - 1
            do i = 0, nx - 1
              do l = 1 , num_atoms
                atom_index = atom_index + 1
                super_atoms(atom_index) = a_names(l)
              ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - !
                super_geometry(atom_index, 1) = unitcell(l, 1) + i * distance_between_unitcells
                super_geometry(atom_index, 2) = unitcell(l, 2) + j * distance_between_unitcells
                super_geometry(atom_index, 3) = unitcell(l, 3) + k * distance_between_unitcells
              ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - !
              end do 
            end do 
          end do
        end do

      end if 

      if (type_of_calculation == "Ring") then 

        rx         = Lx / (2.d0*pi) 

        do i = 0, nx-1
          do j = 1, num_atoms
              atom_index = atom_index + 1
              super_atoms(atom_index) = a_names(j)
              theta = (i * distance_between_unitcells + unitcell(j, 1)) &
            & * 2.d0*pi/(nx * distance_between_unitcells)
              super_geometry(atom_index, 1) = rx * cos(theta)
              super_geometry(atom_index, 2) = rx * sin(theta)
              super_geometry(atom_index, 3) = 0.d0
          end do
        end do

      end if

      open(2,file="supermolecule.mol")
        do i = 1, num_atoms * nx * ny * nz
          write(2,*)  super_atoms(i),             & 
          &           super_geometry(i, 1),       &
          &           super_geometry(i, 2),       &
          &           super_geometry(i, 3)
        end do
      close(2)

      open(3,file="torus_parameters.inp")
        if (keyword(4) == 'Angstrom') then 
          write(3,*) Lx * Ang_par , Ly * Ang_par , Lz * Ang_par 
        else
          write(3,*) Lx , Ly , Lz
        end if
          write(3,*) num_atoms
      close(3)

      ! --------------------------------------------------------------- !

      deallocate(super_atoms)
      deallocate(super_geometry)

      open(4 , file =  "general_parameters.dat" )
        write (4,'(a10)') type_of_calculation
      close(4)

      ! --------------------------------------------------------------- !

end subroutine


subroutine read_2nd_line(type_of_calculation, distance_between_unitcells, Lx , Ly , Lz)
      
      use unitcell_module

      implicit none

      character(len=10)              :: type_of_calculation
      double precision               :: distance_between_unitcells
      double precision               :: Lx , Ly , Lz


      if (type_of_calculation == "OBC") then 

        write(*,'(a)') ""
        write(*,'(a)') "************************************"
        write(*,'(a)') "* Type of calculation: OBC         *"
        write(*,'(a)') "************************************"
        write(*,'(a)') ""

        read(1,*) nx , distance_between_unitcells, Lx

        ny = 1 
        nz = 1 
 
      else if (type_of_calculation == "Ring") then

        write(*,'(a)') ""
        write(*,'(a)') "************************************"
        write(*,'(a)') "* Type of calculation: Ring        *"
        write(*,'(a)') "************************************"
        write(*,'(a)') ""

        read(1,*) nx , distance_between_unitcells, Lx

        ny    = 1 
        nz    = 1

      else if (type_of_calculation == "Tori1D") then 

        write(*,'(a)') ""
        write(*,'(a)') "************************************"
        write(*,'(a)') "* Type of calculation: Toroidal 1D *"
        write(*,'(a)') "************************************"
        write(*,'(a)') ""

        read(1,*) nx , distance_between_unitcells, Lx

        ny   = 1
        nz   = 1

      else if (type_of_calculation == "Tori2D") then 

        write(*,'(a)') ""
        write(*,'(a)') "************************************"
        write(*,'(a)') "* Type of calculation: Toroidal 2D *"
        write(*,'(a)') "************************************"
        write(*,'(a)') ""

        read(1,*) nx , ny , distance_between_unitcells, Lx , Ly

        nz   = 1 

      else if (type_of_calculation == "Tori3D") then 

        write(*,'(a)') ""
        write(*,'(a)') "************************************"
        write(*,'(a)') "* Type of calculation: Toroidal 3D *"
        write(*,'(a)') "************************************"
        write(*,'(a)') ""

        read(1,*) nx , ny , nz , distance_between_unitcells, Lx , Ly , Lz

      else 

        write(*,'(a)') ""
        write(*,'(a)') "Error: Unknown type of calculation"
        write(*,'(a)') ""
        write(*,'(a)') "Please use either    'OBC',   'Ring',          "
        write(*,'(a)') "                  'Tori1D', 'Tori2D', 'Tori3D' "
        write(*,'(a)') ""
        stop

      end if


      Lx = dble(Lx)
      Ly = dble(Ly)
      Lz = dble(Lz)

end subroutine


subroutine read_keywords(keyword)

      implicit none

      character(len=100)             :: line
      integer                        :: io_stat

      character(len=10),intent(out)  :: keyword(100)
      logical                        :: reading_cell

      do
        read(1, '(A)',iostat=io_stat) line
        if (io_stat /= 0) then
          write(*,'(a)') ""
          write(*,'(a)')'Error: Reached end of file without finding'
          write(*,'(a)')'       opening $$ marker in your unitcell'
          write(*,'(a)') ""
          close(1)
          stop
        end if

        if (trim(line) == '$$') then
          reading_cell = .true.
          exit
        end if

        if (trim(line) == 'Read'     )   keyword(1) = 'Read'

        if (trim(line) == 'Integrals')   keyword(2) = 'Integrals'

        if (trim(line) == 'Trexio'   )   keyword(3) = 'Trexio'

        if (trim(line) == 'Angstrom' )   keyword(4) = 'Angstrom'

        if (trim(line) == 'Plot'     )   keyword(5) = 'Plot'

        if (trim(line) == 'Details'  )   keyword(6) = 'Details'

        if (trim(line) == 'MO'       )   keyword(7) = 'MO'

        if (trim(line) == 'UHF'      )   keyword(8) = 'UHF'

        if (trim(line) == 'Huckel'   )   keyword(9) = 'Huckel'

        if (trim(line) == 'MP2'      )   keyword(10) = 'MP2'

        if (trim(line) == 'OV'       )   keyword(11) = 'OV'

        if (trim(line) == 'K'        )   keyword(12) = 'K' 

        if (trim(line) == 'NA'       )   keyword(13) = 'NA' 

        if (trim(line) == 'ERI'      )   keyword(14) = 'ERI'

        if (trim(line) == 'One'      )   keyword(15) = 'One'

        if (trim(line) == 'Orbitals' )   keyword(16) = 'Orbitals'

        if (trim(line) == 'SAO'      )   keyword(17) = 'SAO'

        if (trim(line) == 'Save'     )   keyword(18) = 'Save'

        if (trim(line) == 'DIIS'     )   keyword(19) = 'DIIS'

        if (trim(line) == 'MCDE'     )   keyword(20) = 'MCDE'

        if (trim(line) == 'nquad_128' )  keyword(21) = 'nquad_128'

        if (trim(line) == 'physicist' )  keyword(22) = 'physicist'

        if (trim(line) == 'ccsdt'     )  keyword(23) = 'ccsdt'

        if (trim(line) == 'fast_ccsdt')  keyword(24) = 'fast_ccsdt'   ! [Nov 2026] zgemm k-CCSD residual
        if (trim(line) == 'fast_ccsdt_real') keyword(26) = 'fcc_real' ! [Nov 2026] real-arithmetic fast_ccsdt (tag <= 10 chars)
        if (trim(line) == 'frozen_core') keyword(27) = 'frozcore'     ! [Nov 2026] frozen core in fast k-CCSD(T)
        if (trim(line) == 'checkpoint')  keyword(28) = 'checkpoint'   ! [Nov 2026] binary restart files
        if (trim(line) == 'restart')     keyword(29) = 'restart'      ! [Nov 2026] restart from them
        if (trim(line) == 'fast_quad')   keyword(30) = 'fast_quad'    ! [Nov 2026] mapped 32/40-point t quadrature
        if (trim(line) == 'prim_screen') keyword(31) = 'primscreen'   ! [Nov 2026] primitive Schwarz screening (tag <= 10 chars)
        if (trim(line) == 'cholesky')    keyword(32) = 'cholesky'     ! [Nov 2026] Cholesky ERIs in k-CCSD(T), tol 1e-6
        if (trim(line) == 'cholesky_4')  keyword(33) = 'cholesky_4'   ! [Nov 2026] same, tol 1e-4
        if (trim(line) == 'cholesky_8')  keyword(34) = 'cholesky_8'   ! [Nov 2026] same, tol 1e-8
        if (trim(line) == 'CC_to_disk')  keyword(35) = 'cc_disk'      ! [Nov 2026] vvvv k-blocks on disk in fast k-CCSD(T)
        if (trim(line) == 'direct')      keyword(36) = 'direct'       ! [Nov 2026] integral-direct HF and MP2 (no stored ERIs)
        if (trim(line) == 'density_fitting' .or. trim(line) == 'DF' .or. trim(line) == 'df') &
                                         keyword(37) = 'df'           ! [Nov 2026] density fitting for HF and MP2 (file AuxBasis)
        if (trim(line) == 'df_etb')      keyword(38) = 'df_etb'       ! [Nov 2026] DF with the even-tempered aux basis (not Cholesky fitting)
        if (trim(line) == 'df_cd_6')     keyword(39) = 'df_cd_6'      ! [Nov 2026] Cholesky fitting tolerance 1e-6 (default 1e-8)
        if (trim(line) == 'df_cd_10')    keyword(40) = 'df_cd_10'     ! [Nov 2026] Cholesky fitting tolerance 1e-10



        if (trim(line) == 'debug'     )  keyword(25) = 'debug'

        ! MCDE keywords !

        if (trim(line) == 'Lanczos'  )   keyword(50) = 'Lanczos'

        if (trim(line) == 'Spin'     )   keyword(51) = 'Spin'

        if (trim(line) == 'screen'   )   keyword(52) = 'screen'

        if (trim(line) == 'block_lanczos') keyword(53) = 'blanczos'   ! [Nov 2026] MCDE block Lanczos spectral function

      end do

end subroutine