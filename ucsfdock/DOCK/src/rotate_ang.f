c     a helper subroutine that produces a rotation matrix to 
c     rotate any vector, by the given angle, around a given axis.
c     NL 01/2012

      subroutine rotate_ang(ang, axis, R)

      implicit none

c     input
      real ang         !rotation angle in radians
      real axis(3)     !the axis around which to rotate the input

c     output
      real R(3,3)        !R is the rotation matrix 
     
c     Local variables
      real axis_len    !norm of axis
      real na(3)       !normalized rotation axis

c     main code
c     =========

c     calculate the lengths of the given vectors
      axis_len   = sqrt(axis(1)**2   + axis(2)**2   + axis(3)**2);

c     normalize the axis you want to rotate around
      na(1) = axis(1)/axis_len
      na(2) = axis(2)/axis_len
      na(3) = axis(3)/axis_len

      R(1,1) = cos(ang) + (na(1)**2)*(1-cos(ang)) 
      R(1,2) = na(1)*na(2)*(1-cos(ang)) - (na(3)*sin(ang))
      R(1,3) = na(1)*na(3)*(1-cos(ang)) + (na(2)*sin(ang))

      R(2,1) = na(2)*na(1)*(1-cos(ang)) + (na(3)*sin(ang))
      R(2,2) = cos(ang) + (na(2)**2)*(1-cos(ang))
      R(2,3) = na(2)*na(3)*(1-cos(ang)) - (na(1)*sin(ang))

      R(3,1) = na(3)*na(1)*(1-cos(ang)) - (na(2)*sin(ang))
      R(3,2) = na(3)*na(2)*(1-cos(ang)) + (na(1)*sin(ang))
      R(3,3) = cos(ang) + (na(3)**2)*(1-cos(ang))

      return
      end 
