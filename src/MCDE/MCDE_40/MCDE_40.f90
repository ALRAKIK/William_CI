subroutine MCDE_4_0(nBas, nOcc, energies, ERI)

      ! ==================================================================================================================
      ! -This (4,0)-MCDE module contains all the routines for the creation of the (4,0)-MCDE effective Hamiltonian
      !  Physics notation is assumed for the Electron Repulsion Integrals (ERI)
      !
      !    || __   ||                 ====================================================================================
      !    ||=\_`\=||                                                       Notation
      !    || (__/ ||                 ====================================================================================
      !    ||  | | :-"""-.       
      !    ||==| \/-=-.   \           ERI (ERI_Z)        : Electron Repulsion Integrals (ERI in spin-basis)
      !    ||  |(_|o o/   |_          nBas (nBasZ)       : Amount of (spin) molecular orbitals
      !    ||   \/ "  \   ,_)         nOcc               : Amount of occupied orbitals
      !    ||====\ V  /__/            heffbas            : Dimensions effective Hamiltonian
      !    ||     ;--'  `-.           g0space, nSpace    : 2e1h and 2h1e configuration space and total size
      !    ||    /      .  \          fi, fj, fl         : Occupation numbers
      !    ||===;        \  \         heffMatrix         : Effective Hamiltonian
      !    ||   |         | |         heffbas            : Dimensions effective Hamiltonian
      !    || .-\ '     _/_/          energies           : Hartree-Fock energies in molecular orbital basis
      !    |:'  _;.    (_  \          energiesZ          : Hartree-Fock energies in spin-orbital basis
      !    /  .'  `;\   \\_/          dim1, dim2         : Dimensions of 2p1h and 2h1p rows respectively
      !   |_ /     |||  |\\
      !  /  _)=====|||  | ||          ====================================================================================
      ! /  /|      ||/  / //                                                Routines
      ! \_/||      ( `-/ ||           ====================================================================================
      !    ||======/  /  \\ .-.       g3_spin_potential  : Converts ERI to ERI_Z (integrals to spin basis)
      !    ||      \_/    \'-'/       get_g3space        : Obtains the allowed 3-particle states (2e1h/2h1e) in an array
      !    ||      ||      `"`        heff_4_0             : Constructs the (4,0)-MCDE effective Hamiltonian
      !    ||======|| Horacio         MCDE_4_0             : Combines the above routines into 1 routine to compute the 
      !    ||      ||                                      (4,0)-MCDE effective Hamiltonian
      !    ||      ||                 
      ! ==================================================================================================================


      ! Input variables
      integer,intent(in)                  :: nBas, nOcc
      double precision,intent(in)         :: energies(:)
      double precision, intent(in)        :: ERI(:,:,:,:)

      ! Local Variables
      integer                             :: nBasZ, nSpace, nSpace02, nSpace04, heffbas, nV, n, k, e3, h3
      integer, allocatable                :: g02space(:,:)
      integer, allocatable                :: g04space(:,:)
      double precision, allocatable       :: ERI_Z(:,:,:,:), heffMatrix(:,:)
      double precision, allocatable       :: heff_energies(:)
      double precision, allocatable       :: complex_energies(:)
      double precision, allocatable       :: alpha(:), beta(:)
      double precision, allocatable       :: VL(:)
      double precision, allocatable       :: VR(:)
      
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '            |  (4,0)-MCDE Routine |'
      write(*,*) '---------------------------------------------------'
      write(*,*)  'Total Orbital size       = ', nBas
      write(*,*)  'Occupied Orbital size    = ', nOcc
      write(*,*) '---------------------------------------------------'
      write(*,*)

      ! Variable initialization
      nBasZ   = nBas * 2
      nV = nBas - nOcc
      nSpace = 8*nOcc*nV + 2*(nOcc*(2*nOcc - 1)*nV*(2*nV-1))
      nSpace02 = 8*nOcc*nV
      nSpace04 = 2*(nOcc*(2*nOcc - 1)*nV*(2*nV-1))
      
      heffbas = nSpace
      allocate(heff_energies(heffbas))

      ! 2-particle space creation
      allocate(g02space(2, nSpace02))
      call get_g2space(nBas, nOcc, nSpace02, g02Space)

      ! 4-particle space creation
      allocate(g04space(4, nSpace04))
      call get_g4space(nBas, nOcc, nSpace04, g04Space)

      ! Create the ERI in spin-basis
      allocate(ERI_Z(nBasZ, nBasZ, nBasZ, nBasZ))
      call g4_spin_potential(nBas, ERI, ERI_Z)

      ! Allocate and construct the (4,0)-MCDE effective Hamiltonian
      allocate(heffMatrix(heffbas, heffbas))
      call heff_4_0(nBas, nOcc, nSpace02, nSpace04, heffbas, g04space, g02space, energies, ERI_Z, heffMatrix)

      print*, "Diagonalizing the effective Hamiltonian..."

      allocate(complex_energies(heffbas))
      ! call diagonalize_matrix_nonh(heffbas,heffMatrix,heff_energies,complex_energies)  ! Diagonalize the effective Hamiltonian
      
      ! call matout(heffbas,1,heff_energies)
      write(*,*)
      ! allocate(alpha(heffbas,1))
      ! allocate(beta(heffbas,1))
      call Lanczos_Haydock(heffMatrix, nSpace02, nSpace04, heffbas, 200, alpha, beta)
      write(*,*)'a = ',alpha
      write(*,*)'b = ',beta

      ! NOTE: One can now do an internal diagonalization routine to obtain the eigenvalues/eigenvectors

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '            |  (4,0)-MCDE Routine Finished|'
      write(*,*) '---------------------------------------------------'
      write(*,*)

end subroutine MCDE_4_0