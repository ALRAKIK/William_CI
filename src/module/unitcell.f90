module unitcell_module
      integer              :: nx , ny , nz 
      integer              :: nfuc , N_cell, nauc
      integer              :: ERI_size
      integer, allocatable :: inv_perm(:)
end module unitcell_module