c this sets up a derived type of the 4 important pieces of a phimap, 
c the grid data notwithstanding
      module phimaptype
 
      implicit none
      type phimap
        real :: phiscale !angstroms per grid point
        real, dimension(3) :: oldmid !real space center 
        real :: goff !distance from center to lower grid edges (cubic)
        integer :: nsize !how big the grid is in one dimension, usually 193
      end type phimap

      end module phimaptype
