c this sets up a derived type of the 4 important pieces of a solvmap, 
c the grid data notwithstanding
      module solvmaptype
 
      implicit none
      type solvmap
        integer :: sperang !grid points per angstrom
        integer, dimension(3) :: scadif !grid dimensions in each direction, rect
        real, dimension(3) :: smax !max grid realspace values
      end type solvmap

      end module solvmaptype
