module MCDE_module

      implicit none
      private
      save
   
      public :: MCDE

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

      contains

      subroutine g3_spin_potential(nBas, ERI, ERI_Z)

      ! Input variables
      integer, intent(in)                      :: nBas
      double precision, intent(in)             :: ERI(:, :, :, :)

      ! Output variables
      double precision, intent(out)            :: ERI_Z(:,:,:,:)

      ! Local variables
      integer                                  :: mu, nu, la, si, spinmu, mu2, nu2, la2, si2, spinnu, spinla, spinsi
      double precision                         :: spinmatch1, spinmatch2, spintotalmatch

      ! Debug: Print Header and array sizes
      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '           |  Spin Integral calculation   |'
      write(*,*) '---------------------------------------------------'
      write(*,'(A, 4(I5, 1X))') 'Size of ERI       = ',size(ERI, 1), size(ERI, 2), size(ERI, 3), size(ERI, 4)
      write(*,'(A, 4(I5, 1X))') 'Size of ERI_Z     = ',size(ERI_Z, 1), size(ERI_Z, 2), size(ERI_Z, 3), size(ERI_Z, 4)
      write(*,*) '---------------------------------------------------'
      write(*,*)

      ! Initialize the output array
      ERI_Z = 0.0d0

      ! Loop over spin and basis indices
      do spinmu = 1, 2
         do spinnu = 1, 2
            do spinla = 1, 2
               do spinsi = 1, 2
                  do mu = 1, nBas
                     do nu = 1, nBas
                        do la = 1, nBas
                           do si = 1, nBas

                              ! Set spin-related variables for matrix elements
                              mu2 = (mu) * 2 - 1
                              nu2 = (nu) * 2 - 1
                              la2 = (la) * 2 - 1
                              si2 = (si) * 2 - 1

                              ! Adjust indices for spin states
                              if (spinmu == 1) mu2 = mu2 - 1
                              if (spinnu == 1) nu2 = nu2 - 1
                              if (spinla == 1) la2 = la2 - 1
                              if (spinsi == 1) si2 = si2 - 1

                              ! Initialize spin matching variables
                              spinmatch1 = 0.0d0
                              spinmatch2 = 0.0d0
                              spintotalmatch = 0.0d0

                              ! Check for matching spins between mu, nu, la, and si
                              if (spinmu == spinnu .and. spinla == spinsi) then
                                 spinmatch1 = 1.0d0
                              end if
                              if (spinmu == spinsi .and. spinla == spinnu) then
                                 spinmatch2 = 1.0d0
                              end if

                              ! Condition for total spin matching based on the sum of spin states
                              if (MOD(spinmu + spinnu + spinla + spinsi, 2) == 0) then
                                 spintotalmatch = 1.0d0
                              end if

                              ! Update the output array with the Coulomb matrix element
                              ERI_Z(mu2 + 1, la2 + 1, si2 + 1, nu2 + 1) = spintotalmatch *& 
                                                                        (ERI(mu, la, nu, si) * spinmatch1-&
                                                                         ERI(mu, la, si, nu) * spinmatch2)

                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '      |  Spin Integral calculation finished   |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      end subroutine g3_spin_potential

      subroutine get_g3space(nBas, nOcc, nSpace, g03space)
   
      ! Input variables
      integer, intent(in)           :: nBas, nOcc,  nSpace

      ! Output 
      integer, intent(out)          :: g03space(:,:)
         
      ! Local variables
      integer                       :: i, j, l, idx, nBasZ, nOccZ
      double precision              :: fi, fj, fl
      integer, allocatable          :: tempSpace(:,:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '         |  Auxiliary Subspace Creation  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      ! Temporary workspace (maximum possible size)
      allocate(tempSpace(3, nSpace))

      nOccZ = 2 * nOcc
      nBasZ = 2 * nBas
         
      ! Selection algorithm
      idx = 1
      do i = 1, nBasZ
         fi = 1.0d0
         ! nO times 2 because we use the spin orbital basis set
         if (i > nOccZ) fi = 0.0d0
            
         do j = 1, nBasZ

            fj = 1.0d0
            if (j > nOccZ) fj = 0.0d0
                  
            if (i > j) then
               do l = 1, nBasZ
                  fl = 1.0d0
                        
                  if (l > nOccZ) fl = 0.0d0

                  ! Makes sure it's either 2e1h or 2h1e
                  if ((fi - fl) * (fj - fl) /= 0.0d0) then
                     tempSpace(1, idx) = i
                     tempSpace(2, idx) = j
                     tempSpace(3, idx) = l
                     idx = idx + 1
                  end if
               end do
            end if
         end do
      end do 

      ! Copy to output list
      g03space = tempSpace(:, 1:nSpace)

      deallocate(tempSpace)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '    |  Auxiliary Subspace Creation Finished  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      end subroutine get_g3space

      subroutine heff(nBas, nOcc, nSpace, heffbas, g03space, energies, ERI_Z, heffMatrix)
      
      ! Input variables
      integer, intent(in)                    :: nBas, nOcc, nSpace, heffbas
      double precision, intent(in)           :: energies(:)
      integer, intent(in)                    :: g03space(:,:)
      double precision, intent(in)           :: ERI_Z(:,:,:,:)

      ! Output variables     
      double precision, intent(inout)        :: heffMatrix(:,:)

      ! Local variables
      integer                                :: idx, i, p, q, a, b, c, x, y, z,nOccZ, nBasZ, dim1, dim2
      double precision                       :: fa,fb,fc, dlk, dmj, dio, doj, dmi, body_tmp
      double precision, allocatable          :: energiesZ(:), prefac(:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '        |  Effective Hamiltonian Creation  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      nOccZ = 2 * nOcc
      nBasZ = 2 * nBas
      dim1  = (nOccZ) * (nBasZ - nOccZ) * (nBasZ - nOccZ - 1) / 2
      dim2  = (nOccZ) * (nOccZ - 1) * (nBasZ - nOccZ) / 2
      
      if (nSpace /= dim1 + dim2) then
         write(*,*) "Dimensions incorrect"
         stop
      endif

      ! Writing the double basis hf energies
      allocate(energiesZ(nBasZ))
      do i = 1, nBasZ
         idx=(i-1)/2+1
         energiesZ(i) = energies(idx)
      enddo

      ! Diagonal elements
      heffMatrix = 0.0d0
      do i=1, nBasZ
         heffMatrix(i,i) = energiesZ(i)
      enddo

      do p = 1, nSpace
         a = g03space(1, p)
         b = g03space(2, p)
         c = g03space(3, p)
         heffMatrix(nBasZ + p, nBasZ + p) = energiesZ(a) + energiesZ(b) - energiesZ(c)
      enddo

      ! Right and Left Shoulder (U^I and U^II in ADC notation)
      do i = 1, nBasZ
         do p = nBasZ + 1, heffbas
            idx = p - nBasZ

            a = g03space(1, idx)
            b = g03space(2, idx)
            c = g03space(3, idx)

            heffMatrix(i,p) = heffMatrix(i,p) + ERI_Z(i,c,b,a)
            heffMatrix(p,i) = heffMatrix(p,i) + ERI_Z(i,c,b,a)
         enddo
      enddo

      ! Precompute prefactors
      allocate(prefac(nSpace))
      do p = 1, nSpace
         a = g03space(1, p) !i
         b = g03space(2, p) !j 
         c = g03space(3, p) !l

         fa = merge(0.0d0, 1.0d0, a > nOccZ)
         fb = merge(0.0d0, 1.0d0, b > nOccZ)
         fc = merge(0.0d0, 1.0d0, c > nOccZ)
         if ((fa - fc) * (fb - fc) == 0.0d0) then
            prefac(p) = 0.0d0
         else
            prefac(p) = ((1.0d0-fa)*(1.0d0-fb)*fc - fa*fb*(1.0d0-fc))
         end if
      end do

      !Off-diagonal components for the 3-body (C in ADC notation)
      do p = 1, nSpace
         do q = 1, nSpace

            a = g03space(1, p) !i
            b = g03space(2, p) !j 
            c = g03space(3, p) !l
            x = g03space(1, q) !m
            y = g03space(2, q) !o
            z = g03space(3, q) !k

            dlk = merge(1.0d0, 0.0d0, c == z)
            dmj = merge(1.0d0, 0.0d0, x == b)
            dio = merge(1.0d0, 0.0d0, a == y) 
            doj = merge(1.0d0, 0.0d0, y == b)
            dmi = merge(1.0d0, 0.0d0, a == x)

            body_tmp = prefac(p) * &
                     (dlk*ERI_Z(a,b,y,x) + dmj*ERI_Z(a,z,c,y) + &
                     dio*ERI_Z(b,z,c,x) - doj*ERI_Z(a,z,c,x) - &
                     dmi*ERI_Z(b,z,c,y))

            heffMatrix(nBasZ+p, nBasZ+q) = heffMatrix(nBasZ+p, nBasZ+q) + body_tmp
         enddo
      enddo

      ! NOTE: One can enter their own diagonalization routine here to obtain the eigenvalues and eigenvectors directly.

      deallocate(prefac, energiesZ)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '   |  Effective Hamiltonian Creation Finished |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      end subroutine heff

      subroutine MCDE(nBas, nOcc, energies, ERI)

      use files

      ! Input variables
      integer,intent(in)                  :: nBas, nOcc
      double precision,intent(in)         :: energies(:)
      double precision, intent(in)        :: ERI(:,:,:,:)

      ! Local Variables
      integer                             :: nBasZ, nSpace, heffbas
      integer, allocatable                :: g03space(:,:)
      double precision, allocatable       :: ERI_Z(:,:,:,:), heffMatrix(:,:)
      double precision, allocatable       :: heff_energies(:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '            |  (3,1)-MCDE Routine |'
      write(*,*) '---------------------------------------------------'
      write(*,*)  'Total Orbital size       = ', nBas
      write(*,*)  'Occupied Orbital size    = ', nOcc
      write(*,*) '---------------------------------------------------'
      write(*,*)

      ! Variable initialization
      nBasZ   = nBas * 2
      nSpace  = nOcc*(nBasZ - nOcc*2)*(nBasZ - 2)
      heffbas = nBasZ + nSpace

      ! 3-particle space creation
      allocate(g03space(3, nSpace))
      call get_g3space(nBas, nOcc, nSpace, g03Space)

      ! Create the ERI in spin-basis
      allocate(ERI_Z(nBasZ, nBasZ, nBasZ, nBasZ))
      call g3_spin_potential(nBas, ERI, ERI_Z)

      ! Allocate and construct the (3,1)-MCDE effective Hamiltonian
      allocate(heffMatrix(heffbas, heffbas))
      allocate(heff_energies(heffbas))

      call heff(nBas, nOcc, nSpace, heffbas, g03space, energies, ERI_Z, heffMatrix)

      ! NOTE: One can now do an internal diagonalization routine to obtain the eigenvalues/eigenvectors

      print*, "Diagonalizing the effective Hamiltonian..."

      call diagonalize_matrix(heffbas,heffMatrix,heff_energies)  ! Diagonalize the effective Hamiltonian

      write(outfile,'(A50)') '---------------------------------------'
      write(outfile,'(A50)') '             MCDE energies             '
      write(outfile,'(A50)') '---------------------------------------'
      call matout(heffbas,1,heff_energies)
      write(outfile,*)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '            |  (3,1)-MCDE Routine Finished|'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      end subroutine MCDE

endmodule MCDE_module
