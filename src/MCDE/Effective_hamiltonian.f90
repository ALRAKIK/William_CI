
subroutine heff_4_0(nBas, nOcc, nSpace2, nSpace4, heffbas, g04space, g02space, energies, ERI_Z, heffMatrix)
      
      ! Input variables
      integer         , intent(in)           :: nBas, nOcc, nSpace2, nSpace4, heffbas
      double precision, intent(in)           :: energies(:)
      integer         , intent(in)           :: g02space(:,:)
      integer         , intent(in)           :: g04space(:,:)
      double precision, intent(in)           :: ERI_Z(:,:,:,:)

      ! Output variables     
      double precision, intent(out)          :: heffMatrix(:,:)

      ! Local variables
      integer                                :: idx, i, j, p, q, a, b, c, d, w, x, y, z, k, nOccZ, nBasZ, dim1, dim2
      double precision                       :: fi,fj,fl,fn, fo, fk, dlk, djm, dio, djo, dim, dnp, dlp, dnk, body_tmp
      double precision, allocatable          :: energiesZ(:), prefac2(:), prefac4(:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '        |  Effective Hamiltonian Creation  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      nOccZ = 2 * nOcc
      nBasZ = 2 * nBas
      dim1  = (nOccZ) * (nBasZ - nOccZ) * (nBasZ - nOccZ - 1) / 2
      dim2  = (nOccZ) * (nOccZ - 1) * (nBasZ - nOccZ) / 2
      

      ! Writing the double basis hf energies
      allocate(energiesZ(nBasZ))
      do i = 1, nBasZ
         idx=(i-1)/2+1
         energiesZ(i) = energies(idx)
      enddo

      ! Diagonal elements
      heffMatrix = 0.0d0
      do i = 1, nSpace2
         a = g02space(1, i)
         b = g02space(2, i)
         heffMatrix(i,i) = energiesZ(a) - energiesZ(b)
      enddo

      do p = 1, nSpace4
         a = g04space(1, p)
         b = g04space(2, p)
         c = g04space(3, p)
         d = g04space(4, p)
         heffMatrix(nSpace2 + p, nSpace2 + p) = energiesZ(a) + energiesZ(b) - energiesZ(c) - energiesZ(d)
      enddo

      ! Precompute prefactors
      allocate(prefac2(nSpace2))
      do q = 1, nSpace2
         x = g02space(1, q)
         y = g02space(2, q)

         fo = merge(0.0d0, 1.0d0, x > nOccZ)
         fk = merge(0.0d0, 1.0d0, y > nOccZ)
         if ((fo-fk)==0.0d0) then
            prefac2(q) = 0.0d0
         else
            prefac2(q) = fo-fk
         end if
      end do
      
      allocate(prefac4(nSpace4))
      do p = 1, nSpace4
         a = g04space(1, p) !i
         b = g04space(2, p) !j 
         c = g04space(3, p) !l
         d = g04space(4, p) !n

         fi = merge(0.0d0, 1.0d0, a > nOccZ)
         fj = merge(0.0d0, 1.0d0, b > nOccZ)
         fl = merge(0.0d0, 1.0d0, c > nOccZ)
         fn = merge(0.0d0, 1.0d0, d > nOccz)
         if ((fi-fn)*(fi-fl)*(fj-fn)==0.0d0) then
            prefac4(p) = 0.0d0
         else
            prefac4(p) = (fi-fn)*(fi-fl)*(fj-fn)
         end if
      end do

      ! Head (2-body part)
      do p = 1, nSpace2
         x = g02space(1, p)
         y = g02space(2, p)
         do q = 1, nSpace2
            a = g02space(1, q)
            b = g02space(2, q)

            heffMatrix(p,q) = heffMatrix(p,q) - prefac2(q)*ERI_Z(x,b,a,y)
            if (abs(heffMatrix(p,q)) < 1.0d-10) then
               heffMatrix(p,q) = 0.0d0
            end if
         enddo
      enddo
      
            

      ! Right and Left Shoulder (U^I and U^II in ADC notation)
      do i = 1, nSpace2
        x = g02space(1, i) !o
        y = g02space(2, i) !k
        do p = 1, nSpace4

           a = g04space(1, p) !i
           b = g04space(2, p) !j
           c = g04space(3, p) !l
           d = g04space(4, p) !n

           dlk = merge(1.0d0, 0.0d0, c == y)
           dio = merge(1.0d0, 0.0d0, a == x) 
           djo = merge(1.0d0, 0.0d0, x == b)
           dnk = merge(1.0d0, 0.0d0, d == y)

           heffMatrix(i,nspace2+p) = heffMatrix(i,nspace2+p) - prefac2(i)*&
                (djo*ERI_Z(a,y,c,d) + dlk*ERI_Z(a,b,d,x) + &
                dnk*ERI_Z(a,b,x,c) + dio*ERI_Z(b,y,d,c))
           heffMatrix(nspace2+p,i) = heffMatrix(nspace2+p,i) - prefac2(i)*&
                 (djo*ERI_Z(a,y,c,d) + dlk*ERI_Z(a,b,d,x) + &
                 dnk*ERI_Z(a,b,x,c) + dio*ERI_Z(b,y,d,c))
           if (abs(heffMatrix(nspace2+p,i)) < 1.0d-10) then
               heffMatrix(nspace2+p,i) = 0.0d0
               heffMatrix(i,nspace2+p) = 0.0d0
           end if
         enddo
      enddo

      ! Off-diagonal components for the 3-body (C in ADC notation)
      do p = 1, nSpace4
         do q = 1, nSpace4

            a = g04space(1, p) !i
            b = g04space(2, p) !j 
            c = g04space(3, p) !l
            d = g04space(4, p) !n
            x = g04space(1, q) !m
            y = g04space(2, q) !o
            z = g04space(3, q) !k
            w = g04space(4, q) !p

            dlk = merge(1.0d0, 0.0d0, c == z)
            djm = merge(1.0d0, 0.0d0, x == b)
            dio = merge(1.0d0, 0.0d0, a == y) 
            djo = merge(1.0d0, 0.0d0, y == b)
            dim = merge(1.0d0, 0.0d0, a == x)
            dnp = merge(1.0d0, 0.0d0, d == w)
            dnk = merge(1.0d0, 0.0d0, d == z)
            dlp = merge(1.0d0, 0.0d0, c == w)

            body_tmp = prefac4(p) * &
                     (dim*djo*ERI_Z(w,z,d,c) + dlk*dnp*ERI_Z(a,b,x,y) + &
                     dlp*dnk*ERI_Z(a,b,y,x) + dio*djm*ERI_Z(w,z,d,c) + &
                     dim*dnp*ERI_Z(b,z,c,y) + dim*dlk*ERI_Z(b,w,d,y) + &
                     djo*dnp*ERI_Z(a,z,c,x) + djo*dlk*ERI_Z(a,w,d,x) - &
                     dim*dnk*ERI_Z(b,w,c,y) - dim*dlp*ERI_Z(b,z,d,y) - &
                     djo*dnk*ERI_Z(a,w,c,x) - djo*dlp*ERI_Z(a,z,d,x) - &
                     dio*dnp*ERI_Z(b,z,c,x) - djm*dnp*ERI_Z(a,z,c,y) - &
                     dio*dlk*ERI_Z(b,w,d,x) - djm*dlk*ERI_Z(a,w,d,y) + &
                     dio*dlp*ERI_Z(b,z,d,x) + dio*dnk*ERI_Z(b,w,c,x) + &
                     djm*dlp*ERI_Z(a,z,d,y) + djm*dnk*ERI_Z(a,w,c,y))

            heffMatrix(nSpace2+p, nSpace2+q) = heffMatrix(nSpace2+p, nSpace2+q) - body_tmp
            if (abs(heffMatrix(nspace2+p,nspace2+q)) < 1.0d-10) then
               heffMatrix(nspace2+p,nspace2+q) = 0.0d0
            end if
         enddo
      enddo

      ! NOTE: One can enter their own diagonalization routine here to obtain the eigenvalues and eigenvectors directly.

      deallocate(prefac2, prefac4, energiesZ)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '   |  Effective Hamiltonian Creation Finished |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

end subroutine heff_4_0



subroutine heff_4_2(nBas, nOcc, nSpace2, nSpace4, heffbas, g04space, g02space, energies, ERI_Z, heffMatrix)
      
      ! Input variables
      integer, intent(in)                    :: nBas, nOcc, nSpace2, nSpace4, heffbas
      double precision, intent(in)           :: energies(:)
      integer, intent(in)                    :: g02space(:,:)
      integer, intent(in)                    :: g04space(:,:)
      double precision, intent(in)           :: ERI_Z(:,:,:,:)

      ! Output variables     
      double precision, intent(out)          :: heffMatrix(:,:)

      ! Local variables
      integer                                :: idx, i, j, p, q, a, b, c, d, w, x, y, z, k, nOccZ, nBasZ, dim1, dim2
      double precision                       :: fi,fj,fl,fn, fo, fk, dik, djk, dlk, djm, dio, djo, dlo, dim, dpn, dlp, dlm, dnk, body_tmp
      double precision, allocatable          :: energiesZ(:), prefac2(:), prefac4(:)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '        |  Effective Hamiltonian Creation  |'
      write(*,*) '---------------------------------------------------'
      write(*,*)

      nOccZ = 2 * nOcc
      nBasZ = 2 * nBas
      dim1  = (nOccZ) * (nBasZ - nOccZ) * (nBasZ - nOccZ - 1) / 2
      dim2  = (nOccZ) * (nOccZ - 1) * (nBasZ - nOccZ) / 2
      
      !if (nSpace /= dim1 + dim2) then
      !   write(*,*) "Dimensions incorrect"
      !   stop
      !endif

      ! Writing the double basis hf energies
      allocate(energiesZ(nBasZ))
      do i = 1, nBasZ
         idx=(i-1)/2+1
         energiesZ(i) = energies(idx)
      enddo

      ! Diagonal elements
      heffMatrix = 0.0d0
      do i = 1, nSpace2
         a = g02space(1, i)
         b = g02space(2, i)
         heffMatrix(i,i) = energiesZ(a) + energiesZ(b)
      enddo

      do p = 1, nSpace4
         a = g04space(1, p)
         b = g04space(2, p)
         c = g04space(3, p)
         d = g04space(4, p)
         heffMatrix(nSpace2 + p, nSpace2 + p) = energiesZ(a) + energiesZ(b) + energiesZ(c) - energiesZ(d)
      enddo

      ! Precompute prefactors
      allocate(prefac2(nSpace2))
      do q = 1, nSpace2
         x = g02space(1, q)
         y = g02space(2, q)

         fo = merge(0.0d0, 1.0d0, x > nOccZ)
         fk = merge(0.0d0, 1.0d0, y > nOccZ)
         if ((fo+fk-1)==0.0d0) then
            prefac2(q) = 0.0d0
         else
            prefac2(q) = fo+fk-1
         end if
      end do
      
      allocate(prefac4(nSpace4))
      do p = 1, nSpace4
         a = g04space(1, p) !i
         b = g04space(2, p) !j 
         c = g04space(3, p) !l
         d = g04space(4, p) !n

         fi = merge(0.0d0, 1.0d0, a > nOccZ)
         fj = merge(0.0d0, 1.0d0, b > nOccZ)
         fl = merge(0.0d0, 1.0d0, c > nOccZ)
         fn = merge(0.0d0, 1.0d0, d > nOccz)
         if ((fi-fn)*(fj-fn)*(fl-fn)==0.0d0) then
            prefac4(p) = 0.0d0
         else
            prefac4(p) = (fi-fn)*(fj-fn)*(fl-fn)
         end if
      end do

      ! Head (2-body part)
      do p = 1, nSpace2
         x = g02space(1, p)
         y = g02space(2, p)
         do q = 1, nSpace2
            a = g02space(1, q)
            b = g02space(2, q)

            heffMatrix(p,q) = heffMatrix(p,q) - prefac2(q)*ERI_Z(x,y,b,a)
            if (abs(heffMatrix(p,q))<1.0d-10) then
                 heffMatrix(p,q) = 0.0d0
            end if
         enddo
      enddo
      
            

      ! Right and Left Shoulder (U^I and U^II in ADC notation)
      do i = 1, nSpace2
        x = g02space(1, i) !o
        y = g02space(2, i) !k
        
        do p = 1, nSpace4

            a = g04space(1, p) !i
            b = g04space(2, p) !j
            c = g04space(3, p) !l
            d = g04space(4, p) !n

            dik = merge(1.0d0, 0.0d0, a == y)
            djk = merge(1.0d0, 0.0d0, b == y)
            dlk = merge(1.0d0, 0.0d0, c == y)
            dio = merge(1.0d0, 0.0d0, a == x) 
            djo = merge(1.0d0, 0.0d0, x == b)
            dlo = merge(1.0d0, 0.0d0, c == x)

           heffMatrix(i,nspace2+p) = heffMatrix(i,nspace2+p) - prefac2(i)*&
                 (dik*ERI_Z(c,b,x,d) + dlk*ERI_Z(b,a,x,d) + &
                 djk*ERI_Z(a,c,x,d) + dio*ERI_Z(b,c,y,d) + &
                 djo*ERI_Z(c,a,y,d) + dlo*ERI_Z(a,b,y,d))
           heffMatrix(nspace2+p,i) = heffMatrix(nspace2+p,i) - prefac2(i)*&
                 (dik*ERI_Z(c,b,x,d) + dlk*ERI_Z(b,a,x,d) + &
                 djk*ERI_Z(a,c,x,d) + dio*ERI_Z(b,c,y,d) + &
                 djo*ERI_Z(c,a,y,d) + dlo*ERI_Z(a,b,y,d))
            if (abs(heffMatrix(nSpace2+p,i))<1.0d-10) then
                 heffMatrix(nSpace2+p,i) = 0.0d0
                 heffMatrix(i,nSpace2+p) = 0.0d0
            end if
         enddo
      enddo

      ! Off-diagonal components for the 4-body (C in ADC notation)
      do p = 1, nSpace4
         do q = 1, nSpace4

            a = g04space(1, p) !i
            b = g04space(2, p) !j 
            c = g04space(3, p) !l
            d = g04space(4, p) !n

            x = g04space(1, q) !m
            y = g04space(2, q) !o
            z = g04space(3, q) !k
            w = g04space(4, q) !p

            dlk = merge(1.0d0, 0.0d0, c == z)
            djm = merge(1.0d0, 0.0d0, x == b)
            dio = merge(1.0d0, 0.0d0, a == y) 
            djo = merge(1.0d0, 0.0d0, y == b)
            dim = merge(1.0d0, 0.0d0, a == x)
            dpn = merge(1.0d0, 0.0d0, d == w)
            dlm = merge(1.0d0, 0.0d0, c == x)
            dlo = merge(1.0d0, 0.0d0, c == y)
            dik = merge(1.0d0, 0.0d0, a == z)
            djk = merge(1.0d0, 0.0d0, b == z)

            body_tmp = prefac4(p) * &
                     (- dlk*dpn*ERI_Z(a,b,x,y) - djo*dpn*ERI_Z(a,c,x,z) - &
                     dim*dpn*ERI_Z(b,c,y,z) - dlm*dpn*ERI_Z(a,b,y,z) + &
                     dlo*dpn*ERI_Z(a,b,x,z) + djk*dpn*ERI_Z(a,c,x,y) + &
                     djm*dpn*ERI_Z(a,c,y,z) - dik*dpn*ERI_Z(b,c,x,y) + &
                     dio*dpn*ERI_Z(b,c,x,z) - dim*djo*ERI_Z(c,w,d,z) + &
                     dim*djk*ERI_Z(c,w,d,y) - dio*djk*ERI_Z(c,w,d,x) + &
                     dim*dlo*ERI_Z(b,w,d,z) + dio*dlk*ERI_Z(b,w,d,x) - &
                     djo*dlk*ERI_Z(a,w,d,x) + djm*dlk*ERI_Z(a,w,d,y) - &
                     djm*dlo*ERI_Z(a,w,d,z) - dim*dlk*ERI_Z(b,w,d,y))

            heffMatrix(nSpace2+p, nSpace2+q) = heffMatrix(nSpace2+p, nSpace2+q) - body_tmp
            if (abs(heffMatrix(nSpace2+p,nSpace2+q))<1.0d-10) then
                 heffMatrix(nSpace2+p,nSpace2+q) = 0.0d0
            end if
         enddo
      enddo

      ! NOTE: One can enter their own diagonalization routine here to obtain the eigenvalues and eigenvectors directly.

      deallocate(prefac2, prefac4, energiesZ)

      write(*,*)
      write(*,*) '---------------------------------------------------'
      write(*,*) '   |  Effective Hamiltonian Creation Finished |'
      write(*,*) '---------------------------------------------------'
      write(*,*)
   end subroutine heff_4_2