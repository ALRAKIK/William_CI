subroutine classification(number_of_atoms,number_of_functions,geometry,atoms,ERI)
    
      use atom_basis
      use classification_ERI
      use filter_module

      implicit none 

      integer                       :: i , j , k
      integer                       :: number_of_functions
      integer                       :: number_of_atoms 
      double precision              :: geometry(number_of_atoms,3)
      type(atom)                    :: atoms(number_of_atoms)
      type(ERI_function)            :: ERI  (number_of_functions)
      double precision, allocatable :: f_array_coef(:), f_array_expo(:)
      
      j = 0 
      
      do k = 1 , number_of_atoms

        do i = 1 , atoms(k)%num_s_function
              
          j = j + 1

          ERI(j)%x = geometry(k,1) 
          ERI(j)%y = geometry(k,2) 
          ERI(j)%z = geometry(k,3)

          ERI(j)%orbital     = "s"

          call filter(atoms(k)%coefficient_s(:,i),atoms(k)%exponent_s,f_array_coef,f_array_expo)

          ERI(j)%exponent    = f_array_expo
          ERI(j)%coefficient = f_array_coef

          if (allocated(f_array_coef)) deallocate(f_array_coef)
          if (allocated(f_array_expo)) deallocate(f_array_expo)

        end do 

        do i = 1 , atoms(k)%num_p_function

          j = j + 1

          ERI(j)%x = geometry(k,1) 
          ERI(j)%y = geometry(k,2) 
          ERI(j)%z = geometry(k,3)

          ERI(j)%orbital     = "px"

          call filter(atoms(k)%coefficient_p(:,i),atoms(k)%exponent_p,f_array_coef,f_array_expo)

          ERI(j)%exponent    = f_array_expo
          ERI(j)%coefficient = f_array_coef

          if (allocated(f_array_coef)) deallocate(f_array_coef)
          if (allocated(f_array_expo)) deallocate(f_array_expo)

          j = j + 1

          ERI(j)%x = geometry(k,1) 
          ERI(j)%y = geometry(k,2) 
          ERI(j)%z = geometry(k,3)

          ERI(j)%orbital     = "py"

          call filter(atoms(k)%coefficient_p(:,i),atoms(k)%exponent_p,f_array_coef,f_array_expo)

          ERI(j)%exponent    = f_array_expo
          ERI(j)%coefficient = f_array_coef

          if (allocated(f_array_coef)) deallocate(f_array_coef)
          if (allocated(f_array_expo)) deallocate(f_array_expo)

          j = j + 1

          ERI(j)%x = geometry(k,1) 
          ERI(j)%y = geometry(k,2) 
          ERI(j)%z = geometry(k,3)

          ERI(j)%orbital     = "pz"

          call filter(atoms(k)%coefficient_p(:,i),atoms(k)%exponent_p,f_array_coef,f_array_expo)

          ERI(j)%exponent    = f_array_expo
          ERI(j)%coefficient = f_array_coef

          if (allocated(f_array_coef)) deallocate(f_array_coef)
          if (allocated(f_array_expo)) deallocate(f_array_expo)

          end do 

      end do 


end subroutine classification

subroutine classification_tor(number_of_atoms, number_of_functions, &
                               geometry, atoms, ERI)
      use atom_basis
      use classification_ERI
      use filter_module

      implicit none

      integer                       :: i, j, k
      integer                       :: number_of_functions
      integer                       :: number_of_atoms
      double precision              :: geometry(number_of_atoms,3)
      type(atom)                    :: atoms(number_of_atoms)
      type(ERI_function)            :: ERI(number_of_functions)
      double precision, allocatable :: f_array_coef(:), f_array_expo(:)

      j = 0
      do k = 1, number_of_atoms

        ! --- S ---
        do i = 1, atoms(k)%num_s_function
          j = j + 1
          ERI(j)%x       = geometry(k,1)
          ERI(j)%y       = geometry(k,2)
          ERI(j)%z       = geometry(k,3)
          ERI(j)%orbital = "s"
          call filter(atoms(k)%coefficient_s(:,i), atoms(k)%exponent_s, &
                      f_array_coef, f_array_expo)
          ERI(j)%exponent    = f_array_expo
          ERI(j)%coefficient = f_array_coef
          if (allocated(f_array_coef)) deallocate(f_array_coef)
          if (allocated(f_array_expo)) deallocate(f_array_expo)
        end do

        ! --- PX / PY / PZ ---
        do i = 1, atoms(k)%num_p_function
          
          j = j + 1
          ERI(j)%x       = geometry(k,1)
          ERI(j)%y       = geometry(k,2)
          ERI(j)%z       = geometry(k,3)
          ERI(j)%orbital = "px"
          call filter(atoms(k)%coefficient_px(:,i), atoms(k)%exponent_p, &
                      f_array_coef, f_array_expo)
          ERI(j)%exponent    = f_array_expo
          ERI(j)%coefficient = f_array_coef
          if (allocated(f_array_coef)) deallocate(f_array_coef)
          if (allocated(f_array_expo)) deallocate(f_array_expo)

          j = j + 1
          ERI(j)%x       = geometry(k,1)
          ERI(j)%y       = geometry(k,2)
          ERI(j)%z       = geometry(k,3)
          ERI(j)%orbital = "py"
          call filter(atoms(k)%coefficient_py(:,i), atoms(k)%exponent_p, &
                      f_array_coef, f_array_expo)
          ERI(j)%exponent    = f_array_expo
          ERI(j)%coefficient = f_array_coef
          if (allocated(f_array_coef)) deallocate(f_array_coef)
          if (allocated(f_array_expo)) deallocate(f_array_expo)

          j = j + 1
          ERI(j)%x       = geometry(k,1)
          ERI(j)%y       = geometry(k,2)
          ERI(j)%z       = geometry(k,3)
          ERI(j)%orbital = "pz"
          call filter(atoms(k)%coefficient_pz(:,i), atoms(k)%exponent_p, &
                      f_array_coef, f_array_expo)
          ERI(j)%exponent    = f_array_expo
          ERI(j)%coefficient = f_array_coef
          if (allocated(f_array_coef)) deallocate(f_array_coef)
          if (allocated(f_array_expo)) deallocate(f_array_expo)
        end do

      end do

end subroutine classification_tor


subroutine classification_orbital(number_of_atoms,number_of_functions,geometry,atoms,AO)
    
      use atom_basis
      use classification_ERI
      use filter_module
  
      implicit none 
  

      integer           ,intent(in)   :: number_of_functions
      integer           ,intent(in)   :: number_of_atoms 
      double precision  ,intent(in)   :: geometry(number_of_atoms,3)
      type(atom)        ,intent(in)   :: atoms(number_of_atoms)
      type(ERI_function),intent(out)  :: AO  (number_of_functions)

      integer                         :: i , j , k
        
      j = 0

      do k = 1, number_of_atoms
      
         ! ---- s ---- (1 Cartesian components)
         do i = 1, atoms(k)%num_s_function
            call fill_p("s", atoms(k)%coefficient_s(:,i), atoms(k)%exponent_s, geometry(k,:))
         end do
       
         ! ---- p ---- (3 Cartesian components)
         do i = 1, atoms(k)%num_p_function
            call fill_p("px", atoms(k)%coefficient_p(:,i), atoms(k)%exponent_p, geometry(k,:))
            call fill_p("py", atoms(k)%coefficient_p(:,i), atoms(k)%exponent_p, geometry(k,:))
            call fill_p("pz", atoms(k)%coefficient_p(:,i), atoms(k)%exponent_p, geometry(k,:))
         end do
       
         ! ---- d ----  (6 Cartesian components)
         do i = 1, atoms(k)%num_d_function
            call fill_p("dxx", atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
            call fill_p("dxy", atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
            call fill_p("dxz", atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
            call fill_p("dyy", atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
            call fill_p("dyz", atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
            call fill_p("dzz", atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
         end do
       
         ! ---- f ----  (10 Cartesian components)
         do i = 1, atoms(k)%num_f_function
            call fill_p("fxxx", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fxxy", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fxxz", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fxyy", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fxyz", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fxzz", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fyyy", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fyyz", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fyzz", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fzzz", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
         end do
       
         ! ---- g ----  (15 Cartesian components)
         do i = 1, atoms(k)%num_g_function
            call fill_p("gxxxx", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxxxy", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxxxz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxxyy", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxxyz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxxzz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxyyy", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxyyz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxyzz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxzzz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gyyyy", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gyyyz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gyyzz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gyzzz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gzzzz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
         end do
       
      end do

      contains

      ! Helper: advance j, set position, set label, filter and store
      subroutine fill_p(label, coef_col, expo_shell, geom)

        character(*),     intent(in) :: label
        double precision, intent(in) :: coef_col(:), expo_shell(:), geom(3)

        double precision, allocatable :: tmp_coef(:), tmp_expo(:)

        j = j + 1
        AO(j)%x = geom(1)
        AO(j)%y = geom(2)
        AO(j)%z = geom(3)
        AO(j)%orbital = label
        call filter(coef_col, expo_shell, tmp_coef, tmp_expo)
        AO(j)%exponent    = tmp_expo
        AO(j)%coefficient = tmp_coef
        if (allocated(tmp_coef)) deallocate(tmp_coef)
        if (allocated(tmp_expo)) deallocate(tmp_expo)

      end subroutine fill_p
  
end subroutine classification_orbital

subroutine classification_orbital_tor(number_of_atoms,number_of_functions,geometry,atoms,AO)
    
      use atom_basis
      use classification_ERI
      use filter_module
  
      implicit none 

      integer           ,intent(in)   :: number_of_functions
      integer           ,intent(in)   :: number_of_atoms 
      double precision  ,intent(in)   :: geometry(number_of_atoms,3)
      type(atom)        ,intent(in)   :: atoms(number_of_atoms)
      type(ERI_function),intent(out)  :: AO  (number_of_functions)

      integer                         :: i , j , k
        
      j = 0

      do k = 1, number_of_atoms
      
         ! ---- s ---- (1 Cartesian components)
         do i = 1, atoms(k)%num_s_function
            call fill_p("s", atoms(k)%coefficient_s(:,i), atoms(k)%exponent_s, geometry(k,:))
         end do
       
         ! ---- p ---- (3 Cartesian components)
         do i = 1, atoms(k)%num_p_function
            call fill_p("px", atoms(k)%coefficient_px(:,i), atoms(k)%exponent_p, geometry(k,:))
            call fill_p("py", atoms(k)%coefficient_py(:,i), atoms(k)%exponent_p, geometry(k,:))
            call fill_p("pz", atoms(k)%coefficient_pz(:,i), atoms(k)%exponent_p, geometry(k,:))
         end do
       
         ! ---- d ----  (6 Cartesian components)
         do i = 1, atoms(k)%num_d_function
            call fill_p("dxx", atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
            call fill_p("dxy", atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
            call fill_p("dxz", atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
            call fill_p("dyy", atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
            call fill_p("dyz", atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
            call fill_p("dzz", atoms(k)%coefficient_d(:,i), atoms(k)%exponent_d, geometry(k,:))
         end do
       
         ! ---- f ----  (10 Cartesian components)
         do i = 1, atoms(k)%num_f_function
            call fill_p("fxxx", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fxxy", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fxxz", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fxyy", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fxyz", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fxzz", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fyyy", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fyyz", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fyzz", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
            call fill_p("fzzz", atoms(k)%coefficient_f(:,i), atoms(k)%exponent_f, geometry(k,:))
         end do
       
         ! ---- g ----  (15 Cartesian components)
         do i = 1, atoms(k)%num_g_function
            call fill_p("gxxxx", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxxxy", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxxxz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxxyy", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxxyz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxxzz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxyyy", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxyyz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxyzz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gxzzz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gyyyy", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gyyyz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gyyzz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gyzzz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
            call fill_p("gzzzz", atoms(k)%coefficient_g(:,i), atoms(k)%exponent_g, geometry(k,:))
         end do
       
      end do

        contains

      ! Helper: advance j, set position, set label, filter and store
      subroutine fill_p(label, coef_col, expo_shell, geom)

        character(*),     intent(in) :: label
        double precision, intent(in) :: coef_col(:), expo_shell(:), geom(3)

        double precision, allocatable :: tmp_coef(:), tmp_expo(:)

        j = j + 1
        AO(j)%x = geom(1)
        AO(j)%y = geom(2)
        AO(j)%z = geom(3)
        AO(j)%orbital = label
        call filter(coef_col, expo_shell, tmp_coef, tmp_expo)
        AO(j)%exponent    = tmp_expo
        AO(j)%coefficient = tmp_coef
        if (allocated(tmp_coef)) deallocate(tmp_coef)
        if (allocated(tmp_expo)) deallocate(tmp_expo)

      end subroutine fill_p
  
  
end subroutine classification_orbital_tor 

