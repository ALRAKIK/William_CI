subroutine G2RPA(nBas, nOcc, nSpace, g02space, energies, ERI_Z, heffMatrix)

      implicit none
      
      ! Input variables
      integer, intent(in)                    :: nBas, nOcc, nSpace
      double precision, intent(in)           :: energies(nSpace)
      integer, intent(in)                    :: g02space(2, nSpace)
      double precision, intent(in)           :: ERI_Z(nBas*2, nBas*2, nBas*2, nBas*2)

      ! Output variables     
      double precision, intent(inout)        :: heffMatrix(nSpace, nSpace)

      ! Local variables
      integer                                :: idx, i, p, q, a, b, x, y, nOccZ, nBasZ
      double precision                       :: fo, fk, fj, fl
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

            fo = merge(0.0d0, 1.0d0, x > nOccZ)
            fk = merge(0.0d0, 1.0d0, y > nOccZ)

            fj = merge(0.0d0, 1.0d0, a > nOccZ)
            fl = merge(0.0d0, 1.0d0, b > nOccZ)

            ! if (fj == fo .AND. fk == fl) then

            heffMatrix(p,q) = heffMatrix(p,q) - prefac(q)*ERI_Z(x,b,a,y)
            if (abs(heffMatrix(p,q)) < 1.0d-10) then
               heffMatrix(p,q) = 0.0d0
            end if
            ! end if
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