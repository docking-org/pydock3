c-----------------------------------------------------------------------
c     functions for CHEMGRID                           ECMeng  4/91
c-----------------------------------------------------------------------
      real function dist2(c1, c2)
c
c  --calculates the square of the distance between two points
c
      real c1(3), c2(3)
      dist2 = (c1(1)-c2(1))**2 + (c1(2)-c2(2))**2 +
     &      (c1(3)-c2(3))**2
      return
      end
c-----------------------------------------------------------------------
      integer function indx1(i,j,k,grdpts)
c
c  --converts the 3-dimensional (virtual) indices of a grid point to the
c    actual index in a 1-dimensional array
c
      integer i, j, k
      integer grdpts(3)
c
      indx1 = grdpts(1)*grdpts(2)*(k-1) + grdpts(1)*(j-1) + i
      return
      end
c-----------------------------------------------------------------------
