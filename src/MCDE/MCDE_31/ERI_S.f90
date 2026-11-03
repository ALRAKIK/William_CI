subroutine calculate_ERI_S(nspace,S_inv, eigve, eigva, ERI_Z, ERI_S, nOccZ, nBas)

      implicit none
      
      ! Input variables
      integer         , intent(in)           :: nspace
      double precision, intent(in)           :: ERI_Z(nBas*2, nBas*2, nBas*2, nBas*2), S_inv(nspace, nspace), eigve(nspace, nspace)
      double precision, intent(inout)        :: ERI_S(nBas*2, nBas*2, nBas*2, nBas*2)
      double precision, intent(in)           :: eigva(nBas*2)
      integer         , intent(in)           :: nOccZ, nBas

      ! Local variables
      double precision, allocatable          ::P(:,:,:,:), X(:,:,:,:)
      integer                                ::length, length_eri
      integer                                ::a,b,c,d1,i,j,k,l,l1,l2
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
               do d1 = 1, length_eri
                  do l1 = 1, length
                     do l2 = 1, length

                        fc = merge(0.0d0, 1.0d0, c > nOccZ)
                        fd = merge(0.0d0, 1.0d0, d1 > nOccZ)

                        P(a,b,c,d1) = P(a,b,c,d1) + eigve(a,b)*S_inv(l1,l2)*eigve(c,d1)/eigva(l1)*(fc - fd)
                        ! 

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
                     do d1 = 1, length_eri

                        X(a,b,j,l) = X(a,b,j,l) + P(a,b,c,d1)*ERI_Z(c,j,d1,l)

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
      do i = 1, nBas*2
         write(*,*)"P = ",P(i,i,i,i)
         write(*,*)"X = ",X(i,i,i,i)
      enddo
      deallocate(P,X)

end subroutine calculate_ERI_S