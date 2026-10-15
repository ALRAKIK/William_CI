integer function encode_orbital_pattern(o1, o2, o3, o4)

      implicit none
      character(len=*), intent(in) :: o1, o2, o3, o4
      integer                      :: c1, c2, c3, c4

      ! Direct comparison – Fortran pads "s" with spaces to match the length of o1
      
      c1 = orbital_code(o1)
      c2 = orbital_code(o2)
      c3 = orbital_code(o3)
      c4 = orbital_code(o4)

      encode_orbital_pattern = c1*1000 + c2*100 + c3*10 + c4

contains

      integer function orbital_code(str)
        character(len=*), intent(in) :: str

    
      select case(str)
        case("s")    ; orbital_code = 0
        case("px")   ; orbital_code = 1
        case("py")   ; orbital_code = 2
        case("pz")   ; orbital_code = 3
        case("dxy")  ; orbital_code = 4
        case("dxz")  ; orbital_code = 5
        case("dyz")  ; orbital_code = 6
        case("dx2")  ; orbital_code = 7
        case("dy2")  ; orbital_code = 8
        case("dz2")  ; orbital_code = 9
        case default
        print*, "Orbital not implemented: ", trim(str)
        stop
      end select

      end function orbital_code

end function


integer function encode_orbital_pattern_AO(o1, o2)
      implicit none
      character(len=*), intent(in) :: o1, o2
      integer                      :: c1, c2

      c1 = orbital_code(o1)
      c2 = orbital_code(o2)

      encode_orbital_pattern_AO = c1*10 + c2

      contains

      integer function orbital_code(str)
        character(len=*), intent(in) :: str

    
      select case(str)
        case("s")    ; orbital_code = 0
        case("px")   ; orbital_code = 1
        case("py")   ; orbital_code = 2
        case("pz")   ; orbital_code = 3
        case("dxy")  ; orbital_code = 4
        case("dxz")  ; orbital_code = 5
        case("dyz")  ; orbital_code = 6
        case("dx2")  ; orbital_code = 7
        case("dy2")  ; orbital_code = 8
        case("dz2")  ; orbital_code = 9
        case default
        print*, "Orbital not implemented: ", trim(str)
        stop
      end select

      end function orbital_code

end function