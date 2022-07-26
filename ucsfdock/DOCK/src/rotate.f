c     a helper subroutine that produces a rotation matrix to rotate 
c     one input vector onto a second input vector, around a given 
c     axis.
c     NL 01/2012

      subroutine rotate(invec, target, axis, R)

      implicit none

c     input
      real invec(3)    !input vector
      real target(3)   !target vector to align the input to
      real axis(3)     !the axis around which to rotate the input

c     output
      real R(3,3)        !R is the rotation matrix 
     
c     Local variables
      real invec_len   !norm of invec
      real target_len  !norm of target
      real axis_len    !norm of axis
      real i_norm(3), t_norm(3), a_norm(3) !normalized input
      real dot_invec_target !dot product of invec and target
      real ang         !the angle in RADIANS btw invec and target

      !TEST
      real cross(3), c_len, c_norm(3), a_dot_c


c     main code
c     =========

c     calculate the lengths of the given vectors
      invec_len  = sqrt(invec(1)**2  + invec(2)**2  + invec(3)**2);
      target_len = sqrt(target(1)**2 + target(2)**2 + target(3)**2);
      axis_len   = sqrt(axis(1)**2   + axis(2)**2   + axis(3)**2)

c     normalize input vectors
      i_norm = invec/invec_len
      t_norm = target/target_len
      a_norm = axis/axis_len

c     calc the dot product of the normalized vectors
      dot_invec_target = i_norm(1)*t_norm(1) + 
     &                   i_norm(2)*t_norm(2) +
     &                   i_norm(3)*t_norm(3)

c     calc the cross product of normalized invec and target
      cross(1)=i_norm(2)*t_norm(3)-i_norm(3)*t_norm(2)
      cross(2)=i_norm(3)*t_norm(1)-i_norm(1)*t_norm(3)
      cross(3)=i_norm(1)*t_norm(2)-i_norm(2)*t_norm(1)

c     calc the dot product of the axis with the input vector normal
c     this is actualy were the sign of the angle comes from
      a_dot_c = a_norm(1)*cross(1)+
     &          a_norm(2)*cross(2)+
     &          a_norm(3)*cross(3)

c     calculate the SIGNED angle between invec and target
      ang=atan2(a_dot_c,dot_invec_target)

c     produce a rotation matrix around axis by this ang
      call rotate_ang(ang,axis,R)

      return
      end 
