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

   subroutine g3_spin_potential_RPA(nBas, ERI, ERI_Z)

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
                                                                        (ERI(mu, la, nu, si) * spinmatch1)

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

   end subroutine g3_spin_potential_RPA

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

   subroutine get_g2space(nBas, nOcc, nSpace, g02space)
   
      ! Input variables
      integer, intent(in)           :: nBas, nOcc,  nSpace

      ! Output 
      integer, intent(out)          :: g02space(:,:)
         
      ! Local variables
      integer                       :: j, l, idx, nBasZ, nOccZ
      double precision              :: fj, fl
      integer, allocatable          :: tempSpace(:,:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '         |  Auxiliary Subspace Creation  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      ! Temporary workspace (maximum possible size)
      allocate(tempSpace(2, nSpace))

      nOccZ = 2 * nOcc
      nBasZ = 2 * nBas
         
      ! Selection algorithm
      idx = 1
      do j = 1, nBasZ
         fj = 1.0d0
         ! nO times 2 because we use the spin orbital basis set
         if (j > nOccZ) fj = 0.0d0
            
         do l = 1, nBasZ

            fl = 1.0d0
            if (l > nOccZ) fl = 0.0d0

                        ! Makes sure it's 1e1h
            if (fj-fl /= 0.0d0) then
               tempSpace(1, idx) = j
               tempSpace(2, idx) = l
               idx = idx + 1
            end if
         end do
      end do

      ! Copy to output list
      g02space = tempSpace(:, 1:nSpace)

      deallocate(tempSpace)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '    |  Auxiliary Subspace Creation Finished  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)
   end subroutine get_g2space


   subroutine heff(nBas, nOcc, nSpace, heffbas, g03space, energies, ERI_Z, heffMatrix)
      
      ! Input variables
      integer, intent(in)                    :: nBas, nOcc, nSpace, heffbas
      double precision, intent(in)           :: energies(:)
      integer, intent(in)                    :: g03space(:,:)
      double precision, intent(in)           :: ERI_Z(:,:,:,:)

      ! Output variables     
      double precision, intent(inout)        :: heffMatrix(:,:)

      ! Local variables
      integer                                :: idx, i, j, p, q, a, b, c, x, y, z, k, nOccZ, nBasZ, dim1, dim2
      double precision                       :: fa,fb,fc, fp, fq, dlk, dmj, dio, doj, dmi, body_tmp
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
            if (abs(heffMatrix(i,p))<1.0d-15) then
               heffMatrix(i,p) = 0
               heffMatrix(p,i) = 0
            endif
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
            if (abs(heffMatrix(nBasZ+p, nBasZ+q))<1.0d-12) then
               heffMatrix(nBasZ+p, nBasZ+q) = 0
            endif
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



   subroutine G2RPA(nBas, nOcc, nSpace, g02space, energies, ERI_Z, heffMatrix)
      
      ! Input variables
      integer, intent(in)                    :: nBas, nOcc, nSpace
      double precision, intent(in)           :: energies(:)
      integer, intent(in)                    :: g02space(:,:)
      double precision, intent(in)           :: ERI_Z(:,:,:,:)

      ! Output variables     
      double precision, intent(inout)        :: heffMatrix(:,:)

      ! Local variables
      integer                                :: idx, i, j, p, q, a, b, c, x, y, z, k, nOccZ, nBasZ, dim1, dim2
      double precision                       :: fo, fk
      double precision, allocatable          :: energiesZ(:), prefac(:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '        |  Effective G2RPA Creation  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      nBasZ = nBas*2
      nOccZ = nOcc*2
      ! Writing the double basis hf energies
      allocate(energiesZ(nBasZ))
      do i = 1, nBasZ
         idx=(i-1)/2+1
         energiesZ(i) = energies(idx)
      enddo

      ! Diagonal elements
      heffMatrix = 0.0d0
      do i = 1, nSpace
         a = g02space(1, i)
         b = g02space(2, i)
         heffMatrix(i,i) = energiesZ(a) - energiesZ(b)
      enddo

      allocate(prefac(nSpace))
      do q = 1, nSpace
         x = g02space(1, q)
         y = g02space(2, q)

         fo = merge(0.0d0, 1.0d0, x > nOccZ)
         fk = merge(0.0d0, 1.0d0, y > nOccZ)
         if ((fo-fk)==0.0d0) then
            prefac(q) = 0.0d0
         else
            prefac(q) = fo-fk
         end if
      end do

      do p = 1, nSpace
         x = g02space(1, p)
         y = g02space(2, p)
         do q = 1, nSpace
            a = g02space(1, q)
            b = g02space(2, q)

            heffMatrix(p,q) = heffMatrix(p,q) - prefac(q)*ERI_Z(x,b,a,y)
            if (abs(heffMatrix(p,q)) < 1.0d-10) then
               heffMatrix(p,q) = 0.0d0
            end if
         enddo
      enddo
      

      ! NOTE: One can enter their own diagonalization routine here to obtain the eigenvalues and eigenvectors directly.

      deallocate(prefac, energiesZ)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '      |  Effective G2RPA Creation finished |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

   end subroutine G2RPA

   subroutine Lanczos_Haydock(heff_matrix, nSpace1, nSpace3, heffbas, it, a, b)

      ! Input variables
      integer, intent(in)                    :: nSpace1, nSpace3, heffbas
      double precision, intent(in)           :: heff_matrix(:,:)
      integer, intent(in)                    :: it

      ! Output variables     
      double precision, allocatable, intent(out)          :: a(:)
      double precision, allocatable, intent(out)          :: b(:)

      ! Local variables
      double precision, allocatable          :: s0(:,:), w0(:,:), w(:,:), Temp(:,:), alpha(:,:), beta(:,:), v_prec(:,:), v_curr(:,:)
      double precision, allocatable          :: Krylov(:,:,:)
      double precision, allocatable          :: s(:,:)
      integer                                :: i,j,f

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '             |  Lanczos algorithm  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)
      
      allocate(s0(heffbas,1))
      s0=0
      allocate(w(heffbas,1))
      allocate(alpha(heffbas,1))
      allocate(beta(heffbas,1))
      do i = 1, nSpace1
         s0(i,1) = 1.0d0
      enddo
      allocate(Krylov(heffbas,heffbas,1))
      s0=s0/norm2(s0)
      
      Krylov(1,:,1)=s0(:,1)
      w0=Transpose(MATMUL(heff_matrix,s0))

      alpha(1,1) = dot_product(s0(:,1),w0(1,:))
      beta(1,1)=0.0d0
      w=TRANSPOSE(w0)-alpha(1,1)*s0
      beta(2,1)=NORM2(w)
      Krylov(2,:,:)=w/beta(2,1)

      do i = 2, (it-1)
         v_prec=Krylov(i-1,:,:)
         v_curr=Krylov(i,:,:)
         w0=Transpose(matmul(heff_matrix,v_curr))
         alpha(i,1)=dot_product((v_curr(:,1)),w0(1,:)) 
         w=TRANSPOSE(w0)-(alpha(i,1)*v_curr+beta(i,1)*v_prec)
         do j = 1, i
            w = w - matmul(matmul(Krylov(j,:,:),Transpose(w)),Krylov(j,:,:))
         enddo
         beta(i+1,1) = norm2(w)
         Krylov(i+1,:,:) = w/(beta(i+1,1))
         if (abs(beta(i+1,1)) < 1.0d-12) then
            EXIT 
         endif
      enddo
      f = i
      s=matmul(matmul(TRANSPOSE(Krylov(f,:,:)),heff_matrix),(Krylov(f,:,:)))
      alpha(f,1)=s(1,1)
      allocate(a(f))
      allocate(b(f))
      a=alpha(:f,1)
      b=beta(2:f+1,1)     
   
   end subroutine Lanczos_Haydock



   subroutine Lanczos_spectral_function(e_min,e_max,delta,a,b,A_w,eta)

      ! Input variables
      double precision, intent(in)           :: e_min, e_max, delta
      double precision, intent(in)           :: a(:), b(:)
      double precision, intent(in)           :: eta

      ! Output variables     
      double precision, allocatable, intent(out)          :: A_w(:,:)

      ! Local variables
      integer                                             :: steps
      integer                                             :: i
      double precision, allocatable                       :: omega(:)
      complex, allocatable                                :: z(:), result(:)
      double precision, parameter                         :: PI = 4*atan(1.0d0)

      steps = int((e_max-e_min)/delta)
      allocate(omega(steps))
      allocate(z(steps))
      allocate(result(steps))
      allocate(A_w(steps,2))
      do i = 1,steps
         omega(i)=(i-1)*delta+e_min
      enddo
      z = cmplx(omega,eta)
      result = z-a(size(a))
      do i = size(b), 1, -1
         result = z - a(i) - b(i)**2 / result
      end do
      result = -PI*AIMAG(1/result)
      A_w(:,1) = omega
      A_w(:,2) = result

   end subroutine Lanczos_spectral_function

   subroutine get_g3space_spin_adapted(nBas, nOcc, nSpace, g03space)
   
      ! Input variables
      integer, intent(in)           :: nBas, nOcc,  nSpace

      ! Output 
      integer, intent(out)          :: g03space(:,:)
         
      ! Local variables
      integer                       :: i, j, l, idx
      double precision              :: fi, fj, fl
      integer, allocatable          :: tempSpace(:,:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '         |  Auxiliary Subspace Creation  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      ! Temporary workspace (maximum possible size)
      allocate(tempSpace(3, nSpace))
         
      ! Selection algorithm
      idx = 1
      do i = 1, nBas
         fi = 1.0d0

         if (i > nOcc) fi = 0.0d0
            
         do j = 1, nBas

            fj = 1.0d0
            if (j > nOcc) fj = 0.0d0
                  
            if (i >= j) then
               do l = 1, nBas
                  fl = 1.0d0
                        
                  if (l > nOcc) fl = 0.0d0

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
      write(*,*) '   |  Auxiliary Subspace Spin Adapted Finished  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)
   end subroutine get_g3space_spin_adapted

   integer function d(a,b)

   integer                                   :: a,b
   
   if (a == b) then
      d = 1.0d0
   else
      d = 0.0d0
   end if

   return
   end

   subroutine inv_S(eigve, S)

      ! Input variables
      double precision, intent(in)           :: eigve(:,:)
      double precision, intent(inout)        :: S(:,:)

      ! Local variable
      double precision, allocatable          :: Piv(:,:), LDA(:,:)
      integer                                :: i, j, length, err


      length = size(eigve(1,:))
      allocate(Piv(length,length))
      allocate(LDA(length,length))

      do i=1, length
         do j=1, length
            S(i,j)=dot_product(eigve(i,:),eigve(j,:))
         enddo
      enddo

      call dgetrf(length, length, S, length, Piv, err)

      call dgetri(length, S, length, Piv, LDA, length, err)

   end subroutine inv_S

   subroutine calculate_ERI_S(S_inv, eigve, eigva, ERI_Z, ERI_S, nOccZ, nBas)

      ! Input variables
      double precision, intent(in)           ::ERI_Z(:,:,:,:), S_inv(:,:), eigve(:,:)
      double precision, intent(inout)        ::ERI_S(:,:,:,:)
      double precision, intent(in)           ::eigva(:)
      integer, intent(in)                    ::nOccZ, nBas

      ! Local variables
      double precision, allocatable          ::P(:,:,:,:), X(:,:,:,:)
      integer                                ::length, length_eri
      integer                                ::a,b,c,d,i,j,k,l,l1,l2
      double precision                       ::fc, fd

      length = size(eigva)
      length_eri = nBas*2
      allocate(P(length_eri,length_eri,length_eri,length_eri))
      allocate(X(length_eri,length_eri,length_eri,length_eri))
      P = 0.0d0
      X = 0.0d0
      
      do a = 1, length_eri
         do b = 1, length_eri
            do c = 1, length_eri
               do d = 1, length_eri
                  do l1 = 1, length
                     do l2 = 1, length

                        fc = merge(0.0d0, 1.0d0, c > nOccZ)
                        fd = merge(0.0d0, 1.0d0, d > nOccZ)

                        P(a,b,c,d) = P(a,b,c,d) + eigve(a,b)*S_inv(l1,l2)*eigve(c,d)/eigva(l1)*(fc - fd)

                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo      

      do a = 1, length_eri
         do b = 1, length_eri
            do j = 1, length_eri
               do l = 1, length_eri
                  do c = 1, length_eri
                     do d = 1, length_eri

                        X(a,b,j,l) = X(a,b,j,l) + ERI_Z(c,j,d,l)*P(a,b,c,d)

                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      ERI_S = ERI_Z
      do i = 1, length_eri
         do j = 1, length_eri
            do k = 1, length_eri
               do l = 1, length_eri
                  do a = 1, length_eri
                     do b = 1, length_eri

                        ERI_S(i,j,k,l) = ERI_S(i,j,k,l) + ERI_Z(i,b,k,a)*X(a,b,j,l)

                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

   end subroutine calculate_ERI_S



   subroutine heff_screened(nBas, nOcc, nSpace, heffbas, g03space, energies, ERI_Z, ERI_S, heffMatrix)
      
      ! Input variables
      integer, intent(in)                    :: nBas, nOcc, nSpace, heffbas
      double precision, intent(in)           :: energies(:)
      integer, intent(in)                    :: g03space(:,:)
      double precision, intent(in)           :: ERI_Z(:,:,:,:), ERI_S(:,:,:,:)

      ! Output variables     
      double precision, intent(inout)        :: heffMatrix(:,:)

      ! Local variables
      integer                                :: idx, i, j, p, q, a, b, c, x, y, z, k, nOccZ, nBasZ, dim1, dim2
      double precision                       :: fa,fb,fc, fp, fq, dlk, dmj, dio, doj, dmi, body_tmp
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

            heffMatrix(i,p) = heffMatrix(i,p) + ERI_Z(i,c,a,b) - ERI_Z(i,c,b,a)
            heffMatrix(p,i) = heffMatrix(p,i) + ERI_Z(i,c,a,b) - ERI_Z(i,c,b,a)
            if (abs(heffMatrix(i,p))<1.0d-15) then
               heffMatrix(i,p) = 0
               heffMatrix(p,i) = 0
            endif
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
                     (dlk*(ERI_S(a,b,x,y) - ERI_S(a,b,y,x)) - dmj*(ERI_Z(a,z,c,y) - ERI_S(a,z,y,c)) - &
                     dio*(ERI_Z(b,z,c,x) - ERI_S(b,z,x,c)) + doj*(ERI_Z(a,z,c,x) - ERI_S(a,z,x,c)) + &
                     dmi*(ERI_Z(b,z,c,y) - ERI_S(b,z,y,c)))
            
            heffMatrix(nBasZ+p, nBasZ+q) = heffMatrix(nBasZ+p, nBasZ+q) - body_tmp
            if (abs(heffMatrix(nBasZ+p, nBasZ+q))<1.0d-12) then
               heffMatrix(nBasZ+p, nBasZ+q) = 0
            endif
         enddo
      enddo

      ! NOTE: One can enter their own diagonalization routine here to obtain the eigenvalues and eigenvectors directly.

      deallocate(prefac, energiesZ)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '   |  Effective Hamiltonian Creation Finished |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

   end subroutine heff_screened



   subroutine spin_adapt(nBas, nOcc, energies, ERI, heffMatrix, g03space, heffbas, nSpace)

      ! Input variables
      integer, intent(in)                    :: nBas, nOcc, nSpace, heffbas
      double precision, intent(in)           :: energies(:)
      integer, intent(in)                    :: g03space(:,:)
      double precision, intent(in)           :: ERI(:,:,:,:)

      ! Output variables     
      double precision, intent(inout)        :: heffMatrix(:,:)

      ! Local variables
      integer                                :: idx, i, j, p, q, a, b, c, x, y, z, k, dim1, dim2
      double precision                       :: fa,fb,fc, fp, fq, body_tmp
      double precision, allocatable          :: energiesZ(:), prefac(:), e3(:)
      double precision                       :: C1, C2, A1, A2, F1, F2
      integer, allocatable                   :: indices(:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '            |  Spin Adapted Matrix  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      ! Diagonal elements
      allocate(e3(nSpace))
      heffMatrix = 0.0d0
      do i=1, nBas
         heffMatrix(i,i) = energies(i)
      enddo

      do p = 0, nSpace-1
         a = g03space(1, p+1)
         b = g03space(2, p+1)
         c = g03space(3, p+1)
         e3(p+1) = energies(a) + energies(b) - energies(c)
      enddo

      ! Right and Left Shoulder (U^I and U^II in ADC notation)
      do i = 1, nBas
         do p = 0, nSpace-1

            a = g03space(1, p+1)
            b = g03space(2, p+1)
            c = g03space(3, p+1)

            C1 = SQRT(0.5d0)**d(a,b)*SQRT(0.5d0)*(ERI(a,b,i,c)+ERI(a,b,c,i))
            C2 = SQRT(3.0d0/2.0d0)*(ERI(a,b,i,c)-ERI(a,b,c,i))
            
            heffMatrix(i,nBas+2*p+1) = heffMatrix(i,nBas+2*p+1) + C1
            heffMatrix(i,nBas+2*p+2) = heffMatrix(i,nBas+2*p+2) + C2
            heffMatrix(nBas+2*p+1,i) = heffMatrix(nBas+2*p+1,i) + C1
            heffMatrix(nBas+2*p+2,i) = heffMatrix(nBas+2*p+1,i) + C2

            if (abs(heffMatrix(i,nBas+2*p+1))<1.0d-15) then
               heffMatrix(i,nBas+2*p+1) = 0.0d0
               heffMatrix(nBas+2*p+1,i) = 0.0d0
            endif
            if (abs(heffMatrix(i,nBas+2*p+2))<1.0d-15) then
               heffMatrix(i,nBas+2*p+2) = 0.0d0
               heffMatrix(nBas+2*p+2,i) = 0.0d0
            endif
         enddo
      enddo

      ! Precompute prefactors
      allocate(prefac(nSpace))
      do p = 1, nSpace
         a = g03space(1, p) !i
         b = g03space(2, p) !j 
         c = g03space(3, p) !l

         fa = merge(0.0d0, 1.0d0, a >= nOcc)
         fb = merge(0.0d0, 1.0d0, b >= nOcc)
         fc = merge(0.0d0, 1.0d0, c >= nOcc)
         if ((fa - fc) * (fb - fc) == 0.0d0) then
            prefac(p) = 0.0d0
         else
            prefac(p) = -((1-fa)*(1-fb)*fc-fa*fb*(1-fc))
         end if
      end do

      !Off-diagonal components for the 3-body (C in ADC notation)
      do p = 0, nSpace-1
         do q = 0, nSpace-1

            a = g03space(1, p+1) !i
            b = g03space(2, p+1) !j 
            c = g03space(3, p+1) !l
            x = g03space(1, q+1) !m
            y = g03space(2, q+1) !o
            z = g03space(3, q+1) !k
            
            A1 = SQRT(0.5d0)**d(x,y)*SQRT(0.5d0)**d(a,b)*(-d(c,z)*(ERI(a,b,x,y)+ERI(a,b,y,x))+d(x,b)*(ERI(a,z,y,c)-0.5d0*ERI(a,z,c,y)) & 
            +d(a,y)*(ERI(b,z,x,c)-0.5d0*ERI(b,z,c,x))+d(b,y)*(ERI(a,z,x,c)-0.5d0*ERI(a,z,c,x))+d(a,x)*(ERI(b,z,y,c)-0.5d0*ERI(b,z,c,y)))
            A2 = -d(c,z)*(ERI(a,b,x,y)-ERI(a,b,y,x))-d(x,b)*(ERI(a,z,y,c)-1.5d0*ERI(a,z,c,y))-d(a,y)*(ERI(b,z,x,c)-1.5d0*ERI(b,z,c,x)) & 
            +d(b,y)*(ERI(a,z,x,c)-1.5d0*ERI(a,z,c,x))+d(a,x)*(ERI(b,z,y,c)-1.5d0*ERI(b,z,c,y))
            F1 = SQRT(0.5d0)**d(a,b)*SQRT(3.0d0)/2.0d0*(-d(x,b)*ERI(a,z,c,y)+d(a,y)*ERI(b,z,c,x)+d(y,b)*ERI(a,z,c,x)-d(a,x)*ERI(b,z,c,y))
            F2 = SQRT(0.5d0)**d(x,y)*SQRT(3.0d0)/2.0d0*(d(x,b)*ERI(a,z,c,y)-d(a,y)*ERI(b,z,c,x)+d(y,b)*ERI(a,z,c,x)-d(a,x)*ERI(b,z,c,y))

            heffMatrix(nBas+2*p+1, nBas+2*q+1) = d(a,x)*d(b,y)*d(c,z)*e3(p+1) + prefac(p+1)*A1
            heffMatrix(nBas+2*p+2, nBas+2*q+2) = d(a,x)*d(b,y)*d(c,z)*e3(p+1) + prefac(p+1)*A2
            heffMatrix(nBas+2*p+1, nBas+2*q+2) = prefac(p+1)*F1
            heffMatrix(nBas+2*p+2, nBas+2*q+1) = prefac(p+1)*F2
            
            if (abs(heffMatrix(nBas+2*p+1, nBas+2*q+1))<1.0d-12) then
               heffMatrix(nBas+2*p+1, nBas+2*q+1) = 0.0d0
            endif
            if (abs(heffMatrix(nBas+2*p+2, nBas+2*q+1))<1.0d-12) then
               heffMatrix(nBas+2*p+2, nBas+2*q+1) = 0.0d0
            endif
            if (abs(heffMatrix(nBas+2*p+1, nBas+2*q+2))<1.0d-12) then
               heffMatrix(nBas+2*p+1, nBas+2*q+2) = 0.0d0
            endif
            if (abs(heffMatrix(nBas+2*p+2, nBas+2*q+2))<1.0d-12) then
               heffMatrix(nBas+2*p+2, nBas+2*q+2) = 0.0d0
            endif
         enddo
      enddo

      q=1
      allocate(indices(heffbas))
      do p = 1, heffbas
         if (maxval(heffMatrix(p,:))==sum(ABS(heffMatrix(p,:)))) then
            indices(q)=p
            heffMatrix(p,p)=0.0d0
            q = q+1
         endif
      enddo

      ! NOTE: One can enter their own diagonalization routine here to obtain the eigenvalues and eigenvectors directly.

      deallocate(prefac)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '   |  Effective Hamiltonian Creation Finished |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

   end subroutine spin_adapt



   subroutine MCDE(nBas, nOcc, energies, ERI, lanczos, spin, screen)

      use files

      ! Input variables
      integer,intent(in)                  :: nBas, nOcc
      double precision,intent(in)         :: energies(:)
      double precision, intent(in)        :: ERI(:,:,:,:)
      logical, intent(in)                 :: lanczos, spin, screen

      ! Local Variables
      integer                             :: nBasZ, nSpace, nSpace2, heffbas, nspaceS, heffbasS, nV, i
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
      nBasZ   = nBas * 2
      nV = nBas - nOcc
      nSpace  = nOcc*(nBasZ - nOcc*2)*(nBasZ - 2)
      heffbas = nBasZ + nSpace
      eta = 0.02d0
      nSpaceS = nOcc * (nV * (nV+1)/2) + nV * (nOcc * (nOcc+1)/2)
      heffbasS = nBas + nSpaceS*2
      nSpace2 = 8*nOcc*nV
      open (unit=3,file="unitcell.mol",action="read")
      read(3,*)
      read(3,*)
      read(3,*)e_min, e_max, step
      
      ! 3-particle space creation     
      

      if (screen) then
         ! Create the ERI in spin-basis
         allocate(ERI_Z(nBasZ, nBasZ, nBasZ, nBasZ))
         ! call g3_spin_potential(nBas, ERI, ERI_Z)
         allocate(ERI_S(nBasZ, nBasZ, nBasZ, nBasZ))
         call g3_spin_potential_RPA(nBas, ERI, ERI_Z)
         

         allocate(g02space(2, nSpace2))
         call get_g2space(nBas, nOcc, nSpace2, g02Space)
         allocate(heff_energies(nSpace2))
         allocate(heffMatrix(nSpace2,nSpace2))
         call G2RPA(nBas, nOcc, nSpace2, g02space, energies, ERI_Z, heffMatrix)
         allocate(complex_energies(nSpace2))
         allocate(VR(nspace2,nSpace2))
         allocate(VL(nspace2,nSpace2))
         call diagonalize_matrix_nonh(nSpace2,heffMatrix,heff_energies,complex_energies,VR,VL)
         deallocate(heffMatrix)

         allocate(S_inv(nSpace2,nSpace2))   
         call calculate_ERI_S(S_inv, VR, heff_energies, ERI_Z, ERI_S,nOcc*2,nBas)
         do i = 1, nBas*2
            write(*,*)ERI_Z(i,i,i,i), ERI_S(i,i,i,i)
         enddo
         allocate(heffMatrix(heffbas, heffbas))
         allocate(g03space(3, nSpace))
         call get_g3space(nBas, nOcc, nSpace, g03Space)
         call heff_screened(nBas, nOcc, nSpace, heffbas, g03space, energies, ERI_Z, ERI_S, heffMatrix)  
         deallocate(heff_energies)
         allocate(heff_energies(heffbas))

      else
         allocate(ERI_Z(nBasZ, nBasZ, nBasZ, nBasZ))
         call g3_spin_potential(nBas, ERI, ERI_Z)
      
      

         ! if (spin) then 
         !    ! Spin adaptation
         !    allocate(g03space_spin(3, nSpaceS))
         !    call get_g3space_spin_adapted(nBas, nOcc, nSpaceS, g03space_spin)
         !    allocate(heff_energies(heffbasS))
         !    allocate(heffMatrix(heffbasS,heffbasS))
         !    call spin_adapt(nBas, nOcc, energies, ERI, heffMatrix, g03space_spin, heffbasS, nSpaceS)
         ! else     
         !    ! Allocate and construct the (3,1)-MCDE effective Hamiltonian
         !    allocate(g03space(3, nSpace))
         !    call get_g3space(nBas, nOcc, nSpace, g03Space)
         !    allocate(heff_energies(heffbas))
         !    allocate(heffMatrix(heffbas, heffbas))      
         !    call heff(nBas, nOcc, nSpace, heffbas, g03space, energies, ERI_Z, heffMatrix)
         ! endif
      endif

      ! NOTE: One can now do an internal diagonalization routine to obtain the eigenvalues/eigenvectors

      print*, "Diagonalizing the effective Hamiltonian..."
      if (lanczos) then
         call Lanczos_Haydock(heffMatrix,nBas,nSpace,heffbas,300,alpha,beta)
         ! Diagonalize the effective Hamiltonian
         call Lanczos_spectral_function(e_min, e_max, step, alpha, beta, A_w, eta)
         open (unit=10,file="out.txt",action="write")
         do i=1, size(A_W(:,1))
            write (10,*)  A_w(i,1) , A_w(i,2) 
         enddo
         close (10)
      else
         if (spin) then
            call diagonalize_matrix(heffbasS,heffMatrix,heff_energies)
            call matout(heffbasS,1,heff_energies)
         else
            call diagonalize_matrix(heffbas,heffMatrix,heff_energies)
            call matout(heffbas,1,heff_energies)               
         endif
      endif


      write(outfile,'(A50)') '---------------------------------------'
      write(outfile,'(A50)') '             MCDE energies             '
      write(outfile,'(A50)') '---------------------------------------'
      write(outfile,*)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '            |  (3,1)-MCDE Routine Finished|'
      write(*,*) '---------------------------------------------------'
      write(*,*)

   end subroutine MCDE

endmodule MCDE_module
