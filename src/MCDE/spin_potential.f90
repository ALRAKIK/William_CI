subroutine g4_spin_potential(nBas, ERI, ERI_Z)

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

end subroutine g4_spin_potential