!=====================================================================!
!                                                                     !
!          1D toroidal : NORMALIZATION OF THE ATOMIC ORBITALS         !
!                                                                     !
!   Every Cartesian component (s, px, dxy, fxyz, gxxyy, ...) is       !
!   normalized separately with the 1D toroidal overlap               !
!   (x = Clifford, y,z = real):                                       !
!                                                                     !
!     1) primitive  :  c_i <- c_i / sqrt( <g_i|g_i> )                 !
!     2) contracted :  c   <- c   / sqrt( <phi|phi> )                 !
!                                                                     !
!   Input : raw library coefficients (read_atom_shells_gen +          !
!           classification_tor). Call ONCE (not idempotent: step 1    !
!           assumes raw coefficients). Needs ax, Lx (after Torus_def).!
!=====================================================================!

subroutine normalize_AO_toroidal_1D(nBas, AO, report)

      use files
      use torus_init
      use classification_ERI
      use clifford_tools, only : init_clifford, orbital_to_lxyz

      implicit none

      ! - input / output - !

      integer           ,intent(in)    :: nBas
      logical           ,intent(in)    :: report            ! print the check
      type(ERI_function),intent(inout) :: AO(nBas)

      ! - local - !

      type(ERI_function)               :: prim
      integer                          :: j , i , lxyz(3)
      double precision                 :: r(3) , S , dev

      !-----------------------------------------------------------------!

      call init_clifford()

      dev = 0.d0

      do j = 1 , nBas

        call orbital_to_lxyz(AO(j)%orbital, lxyz)
        r = (/ AO(j)%x , AO(j)%y , AO(j)%z /)

        ! 1) primitives
        prim%x = r(1) ; prim%y = r(2) ; prim%z = r(3) ; prim%orbital = AO(j)%orbital
        do i = 1 , size(AO(j)%exponent)
          prim%exponent    = (/ AO(j)%exponent(i) /)
          prim%coefficient = (/ 1.d0 /)
          call overlap_integral_toroidal_1D(lxyz, lxyz, r, r, prim, prim, S)
          AO(j)%coefficient(i) = AO(j)%coefficient(i) / dsqrt(S)
        end do

        ! 2) contraction
        call overlap_integral_toroidal_1D(lxyz, lxyz, r, r, AO(j), AO(j), S)
        AO(j)%coefficient = AO(j)%coefficient / dsqrt(S)

        ! check
        call overlap_integral_toroidal_1D(lxyz, lxyz, r, r, AO(j), AO(j), S)
        dev = max(dev, abs(S - 1.d0))

      end do

      if (report) then
        write(outfile,'(a)')          ""
        write(outfile,'(a)')          "------------------------------------------------------------"
        write(outfile,'(a)')          " Normalization of the AO (each Cartesian component)       "
        write(outfile,'(a)')          "------------------------------------------------------------"
        write(outfile,'(a,i10)')      " functions normalized             : ", nBas
        write(outfile,'(a,es12.4)')   " max |<phi|phi> - 1|              : ", dev
        write(outfile,'(a)')          "------------------------------------------------------------"
        flush(outfile)
      end if

end subroutine normalize_AO_toroidal_1D
