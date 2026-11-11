subroutine MCDE_3_1(nBas, nOcc, energies, ERI, lanczos, spin, screen)


      ! ==================================================================================================================
      ! -This (3,1)-MCDE module contains all the routines for the creation of the (3,1)-MCDE effective Hamiltonian
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
      !    ||      ||      `"`        heff               : Constructs the (3,1)-MCDE effective Hamiltonian
      !    ||======|| Horacio         MCDE               : Combines the above routines into 1 routine to compute the 
      !    ||      ||                                      (3,1)-MCDE effective Hamiltonian
      !    ||      ||                 
      ! ==================================================================================================================

      use files
      use OMP_LIB

      implicit none

      ! Input variables
      integer         , intent(in)        :: nBas, nOcc
      double precision, intent(in)        :: energies(nBas)
      double precision, intent(in)        :: ERI(nBas,nBas,nBas,nBas)
      logical         , intent(in)        :: lanczos, spin, screen

      ! Local Variables
      integer                             :: nBasZ, nSpace, nSpace2, heffbas, nspaceS, heffbasS, nV, i, it
      integer, allocatable                :: g03space(:,:), g03space_spin(:,:), g02space(:,:)

      double precision                    :: eta, step, e_min, e_max
      double precision, allocatable       :: ERI_Z(:,:,:,:), heffMatrix(:,:), S_inv(:,:), VR(:,:), VL(:,:), ERI_S(:,:,:,:)
      double precision, allocatable       :: heff_energies(:)
      double precision, allocatable       :: alpha(:), beta(:)
      double precision, allocatable       :: A_w(:,:)
      double precision, allocatable       :: complex_energies(:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '            |  (3,1)-MCDE Routine |'
      write(*,*) '---------------------------------------------------'
      write(*,*)  'Total Orbital size       = ', nBas
      write(*,*)  'Occupied Orbital size    = ', nOcc
      write(*,*) '---------------------------------------------------'
      write(*,*)

      ! Variable initialization
      nBasZ    = nBas * 2
      nV       = nBas - nOcc
      nSpace   = nOcc*(nBasZ - nOcc*2)*(nBasZ - 2)
      heffbas  = nBasZ + nSpace
      eta      = 0.02d0
      nSpaceS  = nOcc * (nV * (nV+1)/2) + nV * (nOcc * (nOcc+1)/2)
      heffbasS = nBas + nSpaceS*2
      nSpace2  = 8*nOcc*nV

      ! - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / ! 
      ! read third line of unitcell.mol for e_min, e_max, step, it ! 
      open (unit=3,file="unitcell.mol",action="read")
      read(3,*)
      read(3,*)
      read(3,*)e_min, e_max, step, it
      close(3)
      ! - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / !
   
      ! 3-particle space creation     
      write(*,*)'screen = ', screen

      if (screen) then

         ! Create the ERI in spin-basis

         allocate(ERI_Z(nBasZ, nBasZ, nBasZ, nBasZ))
         allocate(ERI_S(nBasZ, nBasZ, nBasZ, nBasZ))
         call g3_spin_potential_RPA(nBas,nBasZ, ERI, ERI_Z)
         
         allocate(g02space(2, nSpace2))
         call get_g2space(nBas, nOcc, nSpace2, g02Space)
         allocate(heff_energies(nSpace2))
         allocate(heffMatrix(nSpace2,nSpace2))
         call G2RPA(nBas, nOcc, nSpace2, g02space, energies, ERI_Z, heffMatrix)
         heffMatrix = heffMatrix*27.2114d0
         allocate(complex_energies(nSpace2))
         allocate(VR(nspace2,nSpace2))
         allocate(VL(nspace2,nSpace2))
         call diagonalize_matrix_nonh(nSpace2,heffMatrix,heff_energies,complex_energies,VR,VL)
         ! call spectral_function(A_w, VR, eta, heff_energies, e_min, e_max, step, nSpace2, nSpace2)
         ! open (unit=10,file="G2RPA.txt",action="write")
         !    do i=1, size(A_W(:,1))
         !       write (10,*)  A_w(i,1) , A_w(i,2) 
         !    enddo
         ! deallocate(A_w)
         deallocate(heffMatrix)
         allocate(S_inv(nSpace2,nSpace2))
         call inv_S(nspace2,VR, S_inv)   
         ! write(*,*)S_inv
         !write(*,*)VR                                          
         call calculate_ERI_S(nspace2,S_inv, VR, heff_energies, ERI_Z, ERI_S,nOcc*2,nBas)
         do i = 1, nBas
            write(*,*)ERI_S(i,i,i,i) - ERI_Z(i,i,i,i)
         enddo
         deallocate(heff_energies)
      else
         allocate(ERI_Z(nBasZ, nBasZ, nBasZ, nBasZ))
         call g3_spin_potential(nBas,nBasz, ERI, ERI_Z)
      end if     

      if (spin) then 
         ! Spin adaptation
         allocate(g03space_spin(3, nSpaceS))
         call get_g3space_spin_adapted(nBas, nOcc, nSpaceS, g03space_spin)
         allocate(heff_energies(heffbasS))
         allocate(heffMatrix(heffbasS,heffbasS))
         call spin_adapt(nBas, nOcc, energies, ERI, heffMatrix, g03space_spin, heffbasS, nSpaceS)
      else     
            ! Allocate and construct the (3,1)-MCDE effective Hamiltonian
         allocate(g03space(3, nSpace))
         call get_g3space(nBas, nOcc, nSpace, g03Space)
         allocate(heff_energies(heffbas))
         if (screen) then
            allocate(heffMatrix(heffbas, heffbas))
            do i = 1, nBas
               write(*,*)ERI_S(i,i,i,i) - ERI_Z(i,i,i,i)
            enddo
            call heff_screened_3_1(nBas, nOcc, nSpace, heffbas, g03space, energies, ERI_Z, ERI_S, heffMatrix)  
            deallocate(heff_energies)
            allocate(heff_energies(heffbas))
         else
            allocate(heffMatrix(heffbas, heffbas))      
            call heff_3_1(nBas, nOcc, nSpace, heffbas, g03space, energies, ERI_Z, heffMatrix)
         endif
      endif
      heffMatrix = heffMatrix*27.2114d0

      ! NOTE: One can now do an internal diagonalization routine to obtain the eigenvalues/eigenvectors

      print*, "Diagonalizing the effective Hamiltonian..."
      if (lanczos) then
         call Lanczos_Haydock(heffMatrix,nBas,nSpace,heffbas,it,alpha,beta)
         ! Diagonalize the effective Hamiltonian
         call Lanczos_spectral_function(e_min, e_max, step, alpha, beta, A_w, eta)
         if (screen) then
            open (unit=10,file="Lanczos MCDE with screening.txt",action="write")
            do i=1, size(A_W(:,1))
               write (10,*)  A_w(i,1) , A_w(i,2)
            enddo
         else
            open (unit=10,file="Lanczos MCDE.txt",action="write")
            do i=1, size(A_W(:,1))
               write (10,*)  A_w(i,1) , A_w(i,2)
            enddo
         endif
         close (10)
      else
         if (spin) then
            call diagonalize_matrix(heffbasS,heffMatrix,heff_energies)
            call matout(heffbasS,1,heff_energies)
         else
            call diagonalize_matrix_d(heffbas,heffMatrix,heff_energies)
            call matout(heffbas,1,heff_energies)
            allocate(A_w(int((e_max-e_min)/step,4),2))
            call spectral_function(A_w, heffMatrix, eta, heff_energies, e_min, e_max, step, nBas, heffbas)
            if (screen) then
               open (unit=10,file="MCDE screening.txt",action="write")
               do i=1, size(A_W(:,1))
                  write (10,*)  A_w(i,1) , A_w(i,2) 
               enddo
            else
               open (unit=10,file="MCDE.txt",action="write")
               do i=1, size(A_W(:,1))
                  write (10,*)  A_w(i,1) , A_w(i,2) 
               enddo
            endif           
         endif
         close(10)
      endif

      !call diagonalize_matrix_d(heffbas,heffMatrix,heff_energies)  ! Diagonalize the effective Hamiltonian

      write(outfile,'(A50)') '---------------------------------------'
      write(outfile,'(A50)') '             MCDE energies             '
      write(outfile,'(A50)') '---------------------------------------'
      !call matout(heffbas,1,heff_energies)
      write(outfile,*)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '            |  (3,1)-MCDE Routine Finished|'
      write(*,*) '---------------------------------------------------'
      write(*,*)

end subroutine MCDE_3_1