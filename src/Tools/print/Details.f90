subroutine details_integrals(nBas,S,T,V,ERI)

      use files 
      implicit none 

      integer , intent(in)         ::                    nBas 
      double precision, intent(in) ::             S(nBas,nBas)
      double precision, intent(in) ::             T(nBas,nBas)
      double precision, intent(in) ::             V(nBas,nBas)
      double precision, intent(in) :: ERI(nBas,nBas,nBas,nBas)


      ! local ! 

      integer                      :: i, j, k, l 

      open(1,file=trim(tmp_file_name)//"/OV.dat ")
        do i = 1 , size(S,1)
          do j = i , size(S,1)
            if (abs(S(i,j)) > 1e-15 ) write(1,*) i , j , S(i,j)
          end do 
        end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/OV_matrix.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(S,1))
      do i = 1 , size(S,1)
        write(1,'(i3,6x,1000(f16.12,2x))') i ,  (S(i,j),j=1,size(S,1))
      end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/KI.dat ")
      do i = 1 , size(T,1)
        do j = i , size(T,1)
          if (abs(T(i,j)) > 1e-15 ) write(1,*) i, j , T(i,j)
        end do 
      end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/KI_matrix.dat")
        write(1,'(15x,1000(i3,15x))') (i,i=1,size(T,1))
        do i = 1 , size(T,1)
          write(1,'(i3,6x,1000(f16.12,2x))') i ,   (T(i,j),j=1,size(T,1))
        end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/NA.dat ")
        do i = 1 , size(V,1)
          do j = i , size(V,1)
            if (abs(V(i,j)) > 1e-15 ) write(1,*) i , j ,  V(i,j)
          end do 
        end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/NA_matrix.dat ")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(V,1))
      do i = 1 , size(V,1)
        write(1,'(i3,6x,1000(f16.12,2x))') i , (V(i,j),j=1,size(V,1))
      end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/ERI.dat")

        do i = 1, nBas
          do j = 1 , nBas
            do k = 1 , nBas
              do l = 1 , nBas
                if (abs(ERI(i,j,k,l)) > 1.d-12) write(1,*) i , j , k , l , ERI(i,j,k,l)
              end do 
            end do 
          end do 
        end do 
      close(1)

end subroutine details_integrals


subroutine details_integrals_SAO(nBas,S,T,V,ERI)

      use files 
      use unitcell_module
      implicit none 

      integer , intent(in)         ::                    nBas 
      double precision, intent(in) ::             S(nBas,nBas)
      double precision, intent(in) ::             T(nBas,nBas)
      double precision, intent(in) ::             V(nBas,nBas)
      double precision, intent(in) :: ERI(nfuc,nBas,nBas,nBas)


      ! local ! 

      integer                      :: i, j, k, l 

      open(1,file=trim(tmp_file_name)//"/OV.dat ")
        do i = 1 , size(S,1)
          do j = i , size(S,1)
            if (abs(S(i,j)) > 1e-15 ) write(1,*) i , j , S(i,j)
          end do 
        end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/OV_matrix.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(S,1))
      do i = 1 , size(S,1)
        write(1,'(i3,6x,1000(f16.12,2x))') i ,  (S(i,j),j=1,size(S,1))
      end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/KI.dat ")
      do i = 1 , size(T,1)
        do j = i , size(T,1)
          if (abs(T(i,j)) > 1e-15 ) write(1,*) i, j , T(i,j)
        end do 
      end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/KI_matrix.dat")
        write(1,'(15x,1000(i3,15x))') (i,i=1,size(T,1))
        do i = 1 , size(T,1)
          write(1,'(i3,6x,1000(f16.12,2x))') i ,   (T(i,j),j=1,size(T,1))
        end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/NA.dat ")
        do i = 1 , size(V,1)
          do j = i , size(V,1)
            if (abs(V(i,j)) > 1e-15 ) write(1,*) i , j ,  V(i,j)
          end do 
        end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/NA_matrix.dat ")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(V,1))
      do i = 1 , size(V,1)
        write(1,'(i3,6x,1000(f16.12,2x))') i , (V(i,j),j=1,size(V,1))
      end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/ERI.dat")

        do i = 1, nfuc
          do j = 1 , nBas
            do k = 1 , nBas
              do l = 1 , nBas
                if (abs(ERI(i,j,k,l)) > 1.d-15) write(1,*) i , j , k , l , ERI(i,j,k,l)
              end do 
            end do 
          end do 
        end do 
      close(1)

end subroutine details_integrals_SAO


subroutine details_print_HF(nBas,nSCF,Hartree,Exchange,Fock,type)

      use files
      use unitcell_module
      implicit none

      integer,intent(in)           :: nBas
      integer,intent(in)           :: nscf 
      integer,intent(in)           :: type
      double precision, intent(in) :: Hartree (nBas,nBas)
      double precision, intent(in) :: Exchange(nBas,nBas)
      double precision, intent(in) :: Fock    (nBas,nBas)

      ! local !

      integer :: i , j 
      integer :: n_size

      if (type == 0) then 
        n_size = nBas
      else 
        n_size = nfuc
      end if 


      write(HFfile,'(a)') ""
      write(HFfile,'(3a)') "!" ,repeat('-',44), "!"
      write(HFfile,'(a,I3)') "                 Iter  = ", nSCF
      write(HFfile,'(3a)') "!" ,repeat('-',44), "!"
      write(HFfile,'(a)') ""

      call header_HF("Hartree potential J = P(la,si)*ERI(mu,nu,la,si)",-1)
      write(HFfile,'(17x,*(i3,15x))') (i,i=1,n_size)
      do i = 1 , size(Hartree,1)
        write(HFfile,'(i3,6x,*(f16.10,2x))') i ,  (Hartree(j,i),j=1,n_size)
      end do 
      write(HFfile,'(a)') ""

      call header_HF("Exchange potential K = -0.5 P(la,si)*ERI(mu,la,si,nu)",-1)
      write(HFfile,'(17x,*(i3,15x))') (i,i=1,n_size)
      do i = 1 , size(Exchange,1)
        write(HFfile,'(i3,6x,*(f16.10,2x))') i ,  (Exchange(j,i),j=1,n_size)
      end do 
      write(HFfile,'(a)') ""

      call header_HF("Fock matrix F = Hc + J + K", -1)
      write(HFfile,'(17x,*(i3,15x))') (i,i=1,n_size)
      do i = 1 , size(Fock,1)
        write(HFfile,'(i3,6x,*(f16.10,2x))') i ,  (Fock(j,i),j=1,n_size)
      end do 
      write(HFfile,'(a)') ""

end subroutine details_print_HF