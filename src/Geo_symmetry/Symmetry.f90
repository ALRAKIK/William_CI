subroutine one_electron_inv_sym(nfuc, AO, phase)

      use torus_init
      use classification_ERI

      implicit none
      integer,             intent(in)  :: nfuc
      type(ERI_function),  intent(in)  :: AO(*)
      integer,             intent(out) :: phase(nfuc)

      integer                          :: i

      do i = 1, nfuc
         phase(i) = orbital_phase(AO(i)%orbital)
      end do

      contains

      integer function orbital_phase(str)
        
        implicit none 

        character(len=*), intent(in) :: str


        select case(trim(str))

        case("s")  ; orbital_phase =  1
        case("px") ; orbital_phase = -1
        case("py") ; orbital_phase = -1
        case("pz") ; orbital_phase = -1

        case default
            print*, "Orbital not implemented: ", trim(str)
            stop
        end select

      end function orbital_phase

end subroutine one_electron_inv_sym


subroutine setup_cell_inversion(nfuc, nx, ny, nz, cell_in_C)
   integer, intent(in)  :: nfuc, nx, ny, nz
   logical, intent(out), allocatable :: cell_in_C(:)
   integer :: ncell, c, cx, cy, cz, icx, icy, icz, inv_c
   ncell = nx*ny*nz
   allocate(cell_in_C(0:ncell-1))
   do c = 0, ncell-1
      cx = mod(c,       nx)
      cy = mod(c/nx,    ny)
      cz =     c/(nx*ny)
      icx = modulo(-cx, nx); icy = modulo(-cy, ny); icz = modulo(-cz, nz)
      inv_c = icx + nx*icy + nx*ny*icz
      cell_in_C(c) = (c <= inv_c)
   end do
end subroutine

integer function inv_index(j, nfuc, nx, ny, nz)
    use unitcell_module, only : inv_perm
   implicit none
   integer, intent(in) :: j, nfuc, nx, ny, nz
   integer :: c, orb, cx, cy, cz, icx, icy, icz, inv_c, inv_orb
   c   = (j-1)/nfuc
   orb = mod(j-1, nfuc) + 1
   cx = mod(c,    nx); cy = mod(c/nx, ny); cz = c/(nx*ny)
   icx = modulo(-cx, nx); icy = modulo(-cy, ny); icz = modulo(-cz, nz)
   inv_c = icx + nx*icy + nx*ny*icz
   if (allocated(inv_perm)) then
     inv_orb = inv_perm(orb)
   else
     inv_orb = orb
   end if
   inv_index = inv_c*nfuc + inv_orb
end function

logical function valid_phase1(a, b, c, d, nfuc, cell_in_C)
   integer, intent(in) :: a, b, c, d, nfuc
   logical, intent(in) :: cell_in_C(0:*)
   valid_phase1 = .false.
   if (a > nfuc) return
   if (b < a)    return
   if (d < c)    return
   if (.not. (a < c .or. (a == c .and. b <= d))) return
   if (.not. cell_in_C((b-1)/nfuc)) return
   valid_phase1 = .true.
end function

!-----------------------------------------------------------------------!
! Auto-detect the inversion permutation (and the inversion center) from !
! unit cell geometry.                                                   !
!                                                                       !
! Tries every plausible inversion center x0 = (frac(1)+frac(b)+n)/2     !
! with n in {0,1}^3, verifies that:                                     !
!   (a) every atom reflects onto some atom,                             !
!   (b) for the resulting orbital map, the CELL SHIFT is 0 for every    !
!       atom (so that  (c, a) -> (-c, perm(a))  holds).                 !
!                                                                       !
! Works for SC, BCC primitive, FCC primitive, diamond primitive, and    !
! the 1D case 2 (2 atoms at 0 and L/2). Requires that an achievable     !
! inversion center exists with the cell map  c -> -c.                   !
!-----------------------------------------------------------------------!

subroutine auto_detect_inv_perm(orb_atom, orb_local,atom_pos, has_inv)
      
      use torus_init
      use unitcell_module

      implicit none

      integer,             intent(in)  :: orb_atom(nfuc), orb_local(nfuc)
      double precision,    intent(in)  :: atom_pos(nauc,3)
      logical,             intent(out) :: has_inv

      double precision, parameter      :: tol = 1.d-10
      double precision                 :: L(3), frac(nauc,3), x0(3), img(3), d(3), cs
      integer                          :: inv_atom(nauc)
      integer                          :: a, atom_pair, orb, orb2, n1, n2, n3, k, b_img
      logical                          :: all_ok

      L = (/ Lx, Ly, Lz /)

      do a = 1, nauc
        do k = 1, 3
          if (L(k) < 1.d-12) then
            frac(a,k) = 0.d0
          else
            frac(a,k) = modulo(atom_pos(a,k), L(k)) / L(k)
          end if
        end do
      end do

      has_inv  = .false.
      inv_perm = 0

      do atom_pair = 1, nauc
        do n1 = 0, 1
          do n2 = 0, 1
            do n3 = 0, 1
            
              x0(1) = 0.5d0*(frac(1,1) + frac(atom_pair,1) + n1)
              x0(2) = 0.5d0*(frac(1,2) + frac(atom_pair,2) + n2)
              x0(3) = 0.5d0*(frac(1,3) + frac(atom_pair,3) + n3)
            
              ! (a) match each atom to its image
              all_ok   = .true.
              inv_atom = 0
              do a = 1, nauc
                do k = 1, 3
                  img(k) = modulo(2.d0*x0(k) - frac(a,k), 1.d0)
                end do
                inv_atom(a) = 0
                do orb2 = 1, nauc
                  do k = 1, 3
                    d(k) = frac(orb2,k) - img(k)
                    d(k) = d(k) - nint(d(k))
                  end do
                  if (all(abs(d) < tol)) then
                    inv_atom(a) = orb2
                    exit
                  end if
                end do
                if (inv_atom(a) == 0) then
                  all_ok = .false.
                  exit
                end if
              end do
              if (.not. all_ok) cycle
            
              ! (b) cell shift must be zero
              do a = 1, nauc
                do k = 1, 3
                  if (L(k) < 1.d-12) cycle
                  cs = 2.d0*x0(k) - frac(a,k) - frac(inv_atom(a),k)
                  if (abs(cs) > tol) then
                    all_ok = .false.
                    exit
                  end if
                end do
                if (.not. all_ok) exit
              end do
              if (.not. all_ok) cycle
            
              ! (c) build inv_perm
              inv_perm = 0
              do orb = 1, nfuc
                a     = orb_atom(orb)
                b_img = inv_atom(a)
                do orb2 = 1, nfuc
                  if (orb_atom(orb2) == b_img .and. &
                      orb_local(orb2) == orb_local(orb)) then
                    inv_perm(orb) = orb2
                    exit
                  end if
                end do
                if (inv_perm(orb) == 0) then
                  all_ok = .false.
                  exit
                end if
              end do
              if (all_ok) then
                has_inv = .true.
                return
              end if
            
            end do
          end do
        end do
      end do

end subroutine auto_detect_inv_perm




subroutine build_orb_to_atom(nfuc, AO, orb_atom, orb_local)


      use classification_ERI
      implicit none
      integer,             intent(in)  :: nfuc
      type(ERI_function),  intent(in)  :: AO(*)
      integer,             intent(out) :: orb_atom(nfuc), orb_local(nfuc)
      

      double precision, parameter      :: tol = 1.d-10
      double precision                 :: pos(nfuc,3), atom_ref(nfuc,3)
      integer                          :: i, a
      logical                          :: new_atom
      integer                          :: n_atoms_unit

      do i = 1, nfuc
        pos(i,1) = AO(i)%x ; pos(i,2) = AO(i)%y ; pos(i,3) = AO(i)%z
      end do

      n_atoms_unit = 0
      orb_local    = 0

      do i = 1, nfuc
        new_atom = .true.
        do a = 1, n_atoms_unit
          if (maxval(abs(pos(i,:) - atom_ref(a,:))) < tol) then
            new_atom = .false.
            orb_atom(i) = a
            exit
          end if
        end do
        if (new_atom) then
          n_atoms_unit = n_atoms_unit + 1
          atom_ref(n_atoms_unit,:) = pos(i,:)
          orb_atom(i) = n_atoms_unit
        end if
        orb_local(i) = count(orb_atom(1:i-1) == orb_atom(i)) + 1
      end do

end subroutine