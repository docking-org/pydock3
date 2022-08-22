c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c     converts Euler angles to rotation matrix
      subroutine euler2rot(v, match, matchnum)


      use matchtype


      implicit none

      type(matcht), intent(inout) :: match
      integer, intent(in):: matchnum !for getting rotation matrices
      real, dimension(6), intent(in) :: v

c     local variables
      integer :: i,j
      real :: conv
      real :: rephi, rethet, repsi
      real :: coss, cosp, cost, sins, sinp, sint

      conv = (2.0d0 * asin(1.0d0)) / 180.0d0      ! = pi/180
      rephi = v(1) * conv
      rethet = v(2) * conv
      repsi = v(3) * conv

      coss = cos(repsi)
      cosp = cos(rephi)
      cost = cos(rethet)
      sins = sin(repsi)
      sinp = sin(rephi)
      sint = sin(rethet)


C**** Calculate rotation matrix from B to C
      match%rotbc(1, 1, matchnum) =  coss * cosp - cost * sinp * sins
      match%rotbc(1, 2, matchnum) =  coss * sinp + cost * cosp * sins
      match%rotbc(1, 3, matchnum) =  sins * sint
      match%rotbc(2, 1, matchnum) = -sins * cosp - cost * sinp * coss
      match%rotbc(2, 2, matchnum) = -sins * sinp + cost * cosp * coss
      match%rotbc(2, 3, matchnum) =  coss * sint
      match%rotbc(3, 1, matchnum) =  sint * sinp
      match%rotbc(3, 2, matchnum) = -sint * cosp
      match%rotbc(3, 3, matchnum) =  cost

C***  Multiply BC x AB to get AC matrix, add to get new COMR
      do i = 1, 3
        do j = 1, 3
           match%rot(i, j, matchnum) =
     &      match%rotbc(i, 1, matchnum) * match%rotab(1, j, matchnum) +
     &      match%rotbc(i, 2, matchnum) * match%rotab(2, j, matchnum) +
     &      match%rotbc(i, 3, matchnum) * match%rotab(3, j, matchnum)
        enddo
        match%comr(i, matchnum) = match%comrab(i, matchnum) + v(3 + i)
      enddo

      return
      end
