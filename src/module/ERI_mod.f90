module classification_ERI
      
      implicit none

      type :: ERI_function
    
        double precision               :: x , y , z 
        double precision , allocatable :: coefficient(:)
        double precision , allocatable :: exponent(:)
        character(LEN=5)               :: orbital 
    
      end type ERI_function

end module