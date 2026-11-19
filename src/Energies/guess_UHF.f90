subroutine guess_UHF(nBas,c_details,nO,HC,X,ENuc,T,V,P,type_calc,ET1,EV1,ET2,EV2)

      use files

      implicit none 

      ! --------------------------------------------------------------- !
      integer         , intent(in)  :: nBas , nO 
      double precision, intent(in)  :: HC(nBas,nBas), X(nBas,nBas)
      double precision, intent(in)  ::  T(nBas,nBas), V(nBas,nBas)
      double precision, intent(in)  :: ENuc
      logical         , intent(in)  :: c_details
      double precision,external     :: trace_matrix
      ! --------------------------------------------------------------- !
      double precision,allocatable  :: cp(:,:), c(:,:) , e(:) , ct(:,:)
      integer                       :: i      , o 
      double precision              :: ET     , EV
      ! --------------------------------------------------------------- !
      double precision, intent(out) :: P(nbas,nbas)
      double precision, intent(out) :: ET1, EV1, ET2, EV2
      character(len=10), intent(in) :: type_calc
      ! --------------------------------------------------------------- !


      allocate(cp(nbas,nbas),c(nbas,nbas),ct(nbas,nbas), e(nbas))

      cp(:,:) = matmul(transpose(X(:,:)), matmul(HC(:,:), X(:,:)))
      ! [Nov 2026] removed (linearly dependent) vectors -> top of the spectrum
      call shift_dropped_real(nBas,X,cp)

      call diagonalize_matrix(nbas, cp, e)

      c(:,:)  = matmul(X(:,:), cp(:,:))

      if (trim(type_calc) == 'beta' .and. nO < nBas) then
      ! Mix HOMO with LUMO (small rotation)
        do i = 1, nO
          c(:,i) = 0.95d0 * c(:,i) + 0.05d0 * c(:,i+1)
        end do
      end if

      ct(:,:) = transpose(c)

      P(:,:)  = matmul(c(:,1:nO), transpose(c(:,1:nO)))


      ! --------------------------------------------------------------- !
      !                          Print                                  !
      ! --------------------------------------------------------------- !


      ET = trace_matrix(nBas,matmul(P,T))
      EV = trace_matrix(nBas,matmul(P,V))

      if (trim(type_calc) == 'alpha' ) then
        ET1 = ET
        EV1 = EV
      write(outfile,'(a)') ""
      write(outfile,'(a,f16.10)')   "      The Kinetic   Energy (alpha)    = ", ET1
      write(outfile,'(a,f16.10)')   "      The Potential Energy (alpha)    = ", EV1
      write(outfile,'(a,f16.10,a)') "      The Nuclear   Energy (alpha)    = ", ENuc , "  +"
      write(outfile,'(a,f16.10)')   "      -----------------------------------"
      write(outfile,'(a,f16.10)')   "      The guess Energy     (alpha)    = ", ET1+EV1+ENuc
      write(outfile,'(a)') ""
      end if

      if (trim(type_calc) == 'beta' ) then
        ET2 = ET
        EV2 = EV
      write(outfile,'(a)') ""
      write(outfile,'(a,f16.10)')   "      The Kinetic   Energy (beta)     = ", ET2
      write(outfile,'(a,f16.10)')   "      The Potential Energy (beta)     = ", EV2
      write(outfile,'(a,f16.10,a)') "      The Nuclear   Energy (beta)     = ", ENuc , "  +"
      write(outfile,'(a,f16.10)')   "      -----------------------------------"
      write(outfile,'(a,f16.10)')   "      The guess Energy     (beta)     = ", ET2+EV2+ENuc
      write(outfile,'(a)') ""
      end if


      ! --------------------------------------------------------------- !
      deallocate(cp,c,ct,e)
      ! --------------------------------------------------------------- !


end subroutine guess_UHF