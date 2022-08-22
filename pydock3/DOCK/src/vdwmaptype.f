c this sets up a derived type of the important pieces of a vdw map
c the grid data notwithstanding
      module vdwmaptype
 
      implicit none
      type vdwmap
        integer :: ngrd !number of total grid points
        real :: grddiv !grid points per angstroms
        real :: invgrid !angstroms per grid point
        real, dimension(3) :: offset ! box xmin, ymin, zmin
        integer, dimension(3) :: grdpts ! number of grid points along box dim
      end type vdwmap

      end module vdwmaptype
