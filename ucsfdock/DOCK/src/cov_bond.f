c     Main 'machinary' for sampling atom locations based on 
c     dihedrals, bond lengths and angles. Based on fwd-kinematics
c     Cite: Kavraki 1992. 
c     NL 01/2012

      subroutine cov_bond(p1,p2,p3,a,fi,d,out)

      implicit none

c                             out
c      * p1                    * ---
c       \                   t /  /
c        \ u             ___ /  d 
c         \             /   /  /
c          \    /\  v  |  a/  / 
c           *--|----------* ---
c           p2  \/          p3
c               fi

c     input variables
      real p1(3), p2(3), p3(3)  !three points representing 
      ! the covalent atom (p3) and two prev. atoms in chain
      real fi, a, d    !covalent bond parameters (dihedral,angle,len)

c     output
      real out(3)      !output coordinates

c     Local variables
      real Tr(4,4)     !translation matrix
      real v(3)        !vector representing the last bond in the
                       !receptor prior to the covalent bond
      real u(3)        !vector representing the prior to last bond
                       !in the receptor 
      real t(3)        !vector representing the covalent bond
      real new_t(3)    !t after rotation by fi deg
      real tmp_t(4)    !tmp result of 4x4 X 4x1 matmul
      real ang_vu      !the bond angle btw v and u
      real v_len       !v bond length
      real u_len       !u bond length
      real dot_vu      !dot product of v and u
      integer tr_vec(4)!null translation vector
      real T0(4,4),  T1(4,4) !DH matrices
      real E1s(3,3), E2s(3,3) !rotation matrices
      real E1(4,4),  E2(4,4)  !rotation+transformation matrices
      real rot_fi(3,3) !rotation matrix to sample the dihedral
      real z(3), x(3)  ! z and x axis'
      real cross_zu(3), cross_uv(3)! intermediate cross products
      real rot_z(3), rot_x(3) ! axis after rotation
      integer i,j      !loop variables

c     main code
c     =========

c     initialize axis
      z(1)=0; z(2)=0; z(3)=1            !z axis
      x(1)=1; x(2)=0; x(3)=0            !x axis

c     initialize Tr to translate to the coordinates of the anchor
c     atom p2, and tr_vec as the final term in the multiplication chain
      Tr(1,1)=1; Tr(1,2)=0; Tr(1,3)=0; Tr(1,4)=p2(1); 
      Tr(2,1)=0; Tr(2,2)=1; Tr(2,3)=0; Tr(2,4)=p2(2);
      Tr(3,1)=0; Tr(3,2)=0; Tr(3,3)=1; Tr(3,4)=p2(3);
      Tr(4,1)=0; Tr(4,2)=0; Tr(4,3)=0; Tr(4,4)=1;

      tr_vec(1)=0; tr_vec(2)=0; tr_vec(3)=0; tr_vec(4)=1;

c     initialize the vectors that represent the two bonds before the
c     covalent bond.
      v(1)=p3(1)-p2(1); v(2)=p3(2)-p2(2); v(3)=p3(3)-p2(3);
      u(1)=p2(1)-p1(1); u(2)=p2(2)-p1(2); u(3)=p2(3)-p1(3);

c     calculate the bond angle btw v and u (in RADIANS)
      v_len = sqrt(v(1)**2+v(2)**2+v(3)**2)
      u_len = sqrt(u(1)**2+u(2)**2+u(3)**2)
      dot_vu=(v(1)*u(1))+(v(2)*u(2))+(v(3)*u(3))
      ang_vu=acos(dot_vu/(v_len*u_len))

c     calculate the rotation matrices to align the local frame at the
c     anchor atom with the global reference frame.
c     z-axis pointing to the previous bond and x-axis normal to the 
c     plane created by u-v:

      cross_zu(1)=z(2)*u(3)-z(3)*u(2)   !cross product of 
      cross_zu(2)=z(3)*u(1)-z(1)*u(3)   !z axis and u vector
      cross_zu(3)=z(1)*u(2)-z(2)*u(1)
      
      cross_uv(1)=u(2)*v(3)-u(3)*v(2)   !cross product of
      cross_uv(2)=u(3)*v(1)-u(1)*v(3)   !u and v bonds
      cross_uv(3)=u(1)*v(2)-u(2)*v(1)

c     roatate the z-axis to match the vector u (prev bond)
      call rotate (z,u,cross_zu,E1s)
      rot_z = matmul(E1s,z)
      rot_x = matmul(E1s,x)

c     now rotate (the new) x axis to match the normal to the
c     u-v bonds plane. rotation is around the (new) z axis.
      call rotate (rot_x,cross_uv,rot_z,E2s)

c     fill the larger transformation matrices
      do i=1,3 
        E1(i,4)=0; E2(i,4)=0
        E1(4,i)=0; E2(4,i)=0
      enddo
      
      do i=1,3 
        do j=1,3 
          E1(i,j)=E1s(i,j); E2(i,j)=E2s(i,j)
        enddo
      enddo
      
      E1(4,4)=1; E2(4,4)=1

c     calculate the 'initial' Denavit-Hartenberg matrix to roatate
c     around the bond angle v-u
c     notice the dihedral angle is set to 0 is there are only 3 atoms 
c     this is left in (and not zeroed) for clarity
      call dh_mat(ang_vu,0,v_len,T0)

c     calculate the second (T1) Denavit-Hartenberg matrix to roatate
c     around the bond angle v-covalent bond
c     notice the dihedral angle only matters if we continue the chain
      call dh_mat(a,fi,d,T1)

c     chain the matrices to yield the final position 
      tmp_t = matmul(Tr, matmul(E2, matmul(E1, matmul(T0, 
     & matmul(T1,tr_vec)))))

c     t is the vector btw p3 to the new atom
      t(1) = tmp_t(1)-p3(1) 
      t(2) = tmp_t(2)-p3(2) 
      t(3) = tmp_t(3)-p3(3)

c     rotate this vector to the correct fi angle
      call rotate_ang(fi,v,rot_fi)
      new_t = matmul(rot_fi,t)
      
c     output the translated coordinates
      out(1) = new_t(1)+p3(1);
      out(2) = new_t(2)+p3(2);
      out(3) = new_t(3)+p3(3);

      return
      end
