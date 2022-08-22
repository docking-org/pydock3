c     Given a dihedral (fi), bond angle (a) and bond length (d)
c     calculate the Denavit-Hartenberg transformation matrix
c     NL 01/2012

      subroutine dh_mat(a,fi,d,T)

      implicit none

c     input
      real a,fi,d   !dihedral (fi), bond angle (a) and length (d)
c     output
      real T(4,4)

      T(1,1)=cos(fi)       
      T(1,2)=-sin(fi)
      T(1,3)=0
      T(1,4)=0
      T(2,1)=sin(fi)*cos(a) 
      T(2,2)=cos(fi)*cos(a)
      T(2,3)=-1*sin(a)
      T(2,4)=-1*sin(a)*d
      T(3,1)=sin(fi)*sin(a)
      T(3,2)=cos(fi)*sin(a)
      T(3,3)=cos(a)
      T(3,4)=cos(a)*d 
      T(4,1)=0
      T(4,2)=0
      T(4,3)=0
      T(4,4)=1

      return
      end
