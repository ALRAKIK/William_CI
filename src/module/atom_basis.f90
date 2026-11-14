module atom_basis
      
      implicit none

      type :: atom
        integer                         :: charge 
        ! s !
        integer                         :: num_s_function 
        integer                         :: num_exponent_s
        double precision  , allocatable :: exponent_s(:)
        double precision  , allocatable :: coefficient_s(:,:)
        ! p !
        integer                         :: num_p_function
        integer                         :: num_exponent_p
        double precision  , allocatable :: exponent_p(:)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        double precision  , allocatable :: coefficient_p(:,:)   ! for normal Gaussian
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        double precision  , allocatable :: coefficient_px(:,:)
        double precision  , allocatable :: coefficient_py(:,:)
        double precision  , allocatable :: coefficient_pz(:,:)
        ! d !
        integer                         :: num_d_function
        integer                         :: num_exponent_d 
        double precision  , allocatable :: exponent_d(:)        
        double precision  , allocatable :: coefficient_d(:,:)
        ! f ! 
        integer                         :: num_f_function
        integer                         :: num_exponent_f
        double precision  , allocatable :: exponent_f(:)        
        double precision  , allocatable :: coefficient_f(:,:)
        ! g ! 
        integer                         :: num_g_function
        integer                         :: num_exponent_g
        double precision  , allocatable :: exponent_g(:)        
        double precision  , allocatable :: coefficient_g(:,:)

        
      end type atom

      ! [Nov 2026] number of Cartesian basis functions of each atom of the
      ! (super)molecule, filled in main after the basis is read; used by the
      ! Mulliken analysis (the functions of an atom are consecutive)
      integer, allocatable :: nao_atom(:)

end module 