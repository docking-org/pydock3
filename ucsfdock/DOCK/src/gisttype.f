c this sets up a derived type of the 4 important pieces of a phimap, 
c the grid data notwithstanding
      module gisttype
 
      implicit none
      type gistparm
        real :: gistspace !angstroms per grid point
        real, dimension(3) :: orgin !front,lower,left corner
        !integer :: nsize !how big the grid is in one dimension
        integer :: xnsize !how big the grid is in x dimension
        integer :: ynsize !how big the grid is in y dimension
        integer :: znsize !how big the grid is in z dimension
      end type gistparm

      end module gisttype
