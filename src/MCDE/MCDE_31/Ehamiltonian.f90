subroutine heff_3_1(nBas, nOcc, nSpace, heffbas, g03space, energies, ERI_Z, heffMatrix)
      
      implicit none

      ! Input variables
      integer         , intent(in)           :: nBas, nOcc, nSpace, heffbas
      double precision, intent(in)           :: energies(heffbas)
      integer         , intent(in)           :: g03space(3, nSpace)
      double precision, intent(in)           :: ERI_Z(nBas*2, nBas*2, nBas*2, nBas*2)

      ! Output variables     
      double precision, intent(inout)        :: heffMatrix(heffbas, heffbas)

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
      !$omp parallel do collapse(2) default(none) schedule(static) &
      !$omp& shared(nSpace, nBasZ, g03space, prefac, ERI_Z, heffMatrix) &
      !$omp& private(p, q, a, b, c, x, y, z, dlk, dmj, dio, doj, dmi, body_tmp)
      do q = 1, nSpace
         do p = 1, nSpace

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
      !$omp end parallel do

      ! NOTE: One can enter their own diagonalization routine here to obtain the eigenvalues and eigenvectors directly.

      deallocate(prefac, energiesZ)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '   |  Effective Hamiltonian Creation Finished |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

end subroutine heff_3_1


subroutine heff_screened_3_1(nBas, nOcc, nSpace, heffbas, g03space, energies, ERI_Z, ERI_S, heffMatrix)
      
      ! Input variables
      integer         , intent(in)           :: nBas, nOcc, nSpace, heffbas
      double precision, intent(in)           :: energies(heffbas)
      integer         , intent(in)           :: g03space(3, nSpace)
      double precision, intent(in)           :: ERI_Z(nbas*2, nbas*2, nbas*2, nbas*2), ERI_S(nbas*2, nbas*2, nbas*2, nbas*2)

      ! Output variables     
      double precision, intent(inout)        :: heffMatrix(heffbas,heffbas)

      ! Local variables
      integer                                :: idx, i, p, q, a, b, c, x, y, z, nOccZ, nBasZ, dim1, dim2
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

end subroutine heff_screened_3_1