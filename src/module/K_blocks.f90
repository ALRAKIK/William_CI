module k_blocks
  
      implicit none

      private
      public :: K_block , build_k_blocks

      type :: K_block

        integer                      :: dim
        integer        , allocatable :: occupation(:)
        complex(kind=8), allocatable :: Coeff(:,:)
        real(kind=8)   , allocatable :: energy(:)

      end type K_block

contains


      subroutine build_k_blocks(nBas,nO, nBlocks, Cfull, efull, blk)

      use constants_module
      use unitcell_module
      
      implicit none 

      integer, intent(in)                     :: nBas
      integer, intent(in)                     :: nO
      integer, intent(in)                     :: nBlocks
      complex(kind=8) , intent(in)            :: Cfull(nBas,nBas)
      double precision, intent(in)            :: efull(nBas)
      type(K_block),    intent(out)           :: blk(N_cell)
      
      integer                                 :: b, s, e, d

      integer                                 :: i

      double precision                        :: e_sorted(nBas)

      ! --------------------------------------------------------------- !  

        e_sorted = efull

        call sort_array(e_sorted, nBas)
      
        do b = 1, nBlocks
          
          s = (b-1)*nBas/nBlocks + 1
          e = b*nBas/nBlocks
          d = e - s + 1
        
          blk(b)%dim    = d
        
          allocate(blk(b)%Coeff(d,d), blk(b)%energy(d),blk(b)%occupation(d))

          blk(b)%Coeff  = Cfull(s:e, s:e)
          blk(b)%energy = efull(s:e)

          do i = 1 , d 
            if ( blk(b)%energy(i) <= e_sorted(nO) ) blk(b)%occupation(i) = 2
          end do 

        end do



      end subroutine build_k_blocks


end module k_blocks