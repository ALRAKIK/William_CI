subroutine spin_adapt(nBas, nOcc, energies, ERI, heffMatrix, g03space, heffbas, nSpace)

      implicit none

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
      
      integer                                :: d

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