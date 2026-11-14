subroutine classification(number_of_atoms, number_of_functions, geometry, atoms, AO)

      use atom_basis
      use classification_ERI
      use filter_module

      implicit none

      integer           ,intent(in)   :: number_of_functions
      integer           ,intent(in)   :: number_of_atoms
      double precision  ,intent(in)   :: geometry(number_of_atoms,3)
      type(atom)        ,intent(in)   :: atoms(number_of_atoms)
      type(ERI_function),intent(out)  :: AO(number_of_functions)

      character(len=5), parameter     :: Lp(3)  = [character(len=5) :: 'px','py','pz']
      character(len=5), parameter     :: Ld(6)  = [character(len=5) :: 'dxx','dxy','dxz','dyy','dyz','dzz']
      character(len=5), parameter     :: Lf(10) = [character(len=5) :: 'fxxx','fxxy','fxxz','fxyy','fxyz', &
                                                                       'fxzz','fyyy','fyyz','fyzz','fzzz']
      character(len=5), parameter     :: Lg(15) = [character(len=5) :: 'gxxxx','gxxxy','gxxxz','gxxyy','gxxyz', &
                                                                       'gxxzz','gxyyy','gxyyz','gxyzz','gxzzz', &
                                                                       'gyyyy','gyyyz','gyyzz','gyzzz','gzzzz']
      integer :: i , j , k , c

      j = 0
      do k = 1 , number_of_atoms
        do i = 1 , atoms(k)%num_s_function
          call add("s", atoms(k)%coefficient_s(:,i), atoms(k)%exponent_s, geometry(k,:))
        end do
        do i = 1 , atoms(k)%num_p_function
          do c = 1 , 3
            call add(Lp(c), atoms(k)%coefficient_p(:,i), atoms(k)%exponent_p, geometry(k,:))
          end do
        end do
        do i = 1 , atoms(k)%num_d_function
          do c = 1 , 6
            call add(Ld(c), atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
          end do
        end do
        do i = 1 , atoms(k)%num_f_function
          do c = 1 , 10
            call add(Lf(c), atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
          end do
        end do
        do i = 1 , atoms(k)%num_g_function
          do c = 1 , 15
            call add(Lg(c), atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
          end do
        end do
      end do

      if (j /= number_of_functions) then
        write(*,'(A,I0,A,I0)') ' classification_tor_gen: built ', j, ' functions, expected ', number_of_functions
        stop
      end if

contains

      subroutine add(label, coef_col, expo_shell, geom)
        character(*),     intent(in) :: label
        double precision, intent(in) :: coef_col(:), expo_shell(:), geom(3)
        double precision, allocatable :: tmp_coef(:), tmp_expo(:)
        j = j + 1
        AO(j)%x = geom(1) ; AO(j)%y = geom(2) ; AO(j)%z = geom(3)
        AO(j)%orbital = label
        call filter(coef_col, expo_shell, tmp_coef, tmp_expo)       ! drop zero coefficients
        AO(j)%exponent    = tmp_expo
        AO(j)%coefficient = tmp_coef
        if (allocated(tmp_coef)) deallocate(tmp_coef)
        if (allocated(tmp_expo)) deallocate(tmp_expo)
      end subroutine add

end subroutine classification