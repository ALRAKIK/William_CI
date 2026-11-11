subroutine normalize_basis()

      use files 
      
      implicit none 

      ! local ! 

      integer                       :: n_gaussian , n_contraction
      integer                       :: i , j , n_type
      character(len=100)            :: lines

      double precision,allocatable  :: exponent(:) , contraction(:,:) , contractionN(:,:) 

      open(1,file=trim(tmp_file_name)//"/Basis_scratch")
      open(2,file=trim(tmp_file_name)//"/Basis_normalized")

      do
      
        read(1,'(A)',end=3) lines 
      
        write(2,'(A)') trim(lines) 
      
        if (lines == "$ S-TYPE FUNCTIONS") then  
          read(1,*) n_gaussian , n_contraction
          allocate(exponent(n_gaussian),contraction(n_gaussian,n_contraction),contractionN(n_gaussian,n_contraction))
          do i = 1 , n_gaussian
            read(1,*) exponent(i) , (contraction(i,j),j=1,n_contraction)
          end do 
          n_type = 1 
          call norm_orb(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,*) exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN)
        end if
      
        ! P function ! 
      
        if (lines == "$ P-TYPE FUNCTIONS") then 
          read(1,*) n_gaussian , n_contraction
          allocate(exponent(n_gaussian),contraction(n_gaussian,n_contraction),contractionN(n_gaussian,n_contraction))
          do i = 1 , n_gaussian
            read(1,*) exponent(i) , (contraction(i,j),j=1,n_contraction)
          end do 
          n_type = 2 
          call norm_orb(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,*) exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN) 
        end if

        ! d function ! 

        if (lines == "$ D-TYPE FUNCTIONS") then 
          read(1,*) n_gaussian , n_contraction
          allocate(exponent(n_gaussian),contraction(n_gaussian,n_contraction),contractionN(n_gaussian,n_contraction))
          do i = 1 , n_gaussian
            read(1,*) exponent(i) , (contraction(i,j),j=1,n_contraction)
          end do 
          n_type = 3 
          call norm_orb(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,*) exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN) 
        end if

        ! F function ! 

        if (lines == "$ F-TYPE FUNCTIONS") then 
          read(1,*) n_gaussian , n_contraction
          allocate(exponent(n_gaussian),contraction(n_gaussian,n_contraction),contractionN(n_gaussian,n_contraction))
          do i = 1 , n_gaussian
            read(1,*) exponent(i) , (contraction(i,j),j=1,n_contraction)
          end do 
          n_type = 4 
          call norm_orb(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,*) exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN) 
        end if

        ! G function ! 

        if (lines == "$ G-TYPE FUNCTIONS") then 
          read(1,*) n_gaussian , n_contraction
          allocate(exponent(n_gaussian),contraction(n_gaussian,n_contraction),contractionN(n_gaussian,n_contraction))
          do i = 1 , n_gaussian
            read(1,*) exponent(i) , (contraction(i,j),j=1,n_contraction)
          end do 
          n_type = 5 
          call norm_orb(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,*) exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN) 
        end if
      
      end do 

3     close(1)

      close(2)

end subroutine normalize_basis



subroutine norm_orb(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN)

      implicit none 

      integer          ,intent(in)  :: n_gaussian , n_contraction
      integer          ,intent(in)  :: n_type 
      double precision ,intent(in)  :: exponent(n_gaussian) , contraction(n_gaussian,n_contraction)
      double precision ,intent(out) :: contractionN(n_gaussian,n_contraction)

      ! local ! 

      integer          :: n , i , j

      double precision :: D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0   ,&
                          & DP25 = 0.25D0, DP5 = 0.5D0, DP75 = 0.75D0      ,&
                          & THRMIN = 1.D-17

      double precision :: PIPPI , sum , T

      double precision , parameter :: pi = 3.14159265358979323846D00

      ! code !

      PIPPI = (DP5/PI)**DP75

      do n = 1 , n_contraction
        sum = D0 
        do i = 1 , n_gaussian
          do j = 1 , n_gaussian
            T = D2*SQRT(exponent(i)*exponent(j))/(exponent(i)+exponent(j))
            sum  = sum + contraction(i,n)*contraction(j,n)*(T**(n_type + DP5))
          end do 
        end do 
        IF (SQRT(sum) .LT. THRMIN) GOTO 10
        sum=D1/SQRT(sum)
        do j = 1, n_gaussian
          contractionN(j,n)= contraction(j,n)*sum*(D4*exponent(j))**(DP5*n_type+DP25)*PIPPI
        end do 
        
      end do 

      RETURN
10    continue 
      
      write(*,"(a)") " You Have Zero Norm" 

end subroutine