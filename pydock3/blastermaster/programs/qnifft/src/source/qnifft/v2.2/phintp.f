      subroutine phintp(gp,phi)
c
c interpolates the potential at any point inside
c a cubical volume using the potential values at the
c 8 vertices by means of a trilinear function:
C	
C W = A1.XYZ + A2.XY + A3.XZ + A4.YZ + A5.X +
C		A6.Y + A7.Z + A8
C	
C where Ai coefficients are linear combinations of Wi at
C the cube corners
c
	include 'qdiffpar.h'
c---------------------------------------------------

      dimension gp(3)
c
c return 0.0 if outside grid
c
	rgrid = igrid
c	print *,'in phintp: gp ',gp ! debug
	do 9000 k = 1,3
	  if((gp(k).le.1).or.(gp(k).ge.rgrid)) then ! kas 1-feb-01 fix limit cases
	    phi = 0.0
	    return
	  end if
c	end do
9000	continue
c
c find lower left bottom grid point
c
	nx = int(gp(1))
      ny = int(gp(2))
      nz = int(gp(3))
c
c calculate cube coordinates of point
c
	xg = nx
      yg = ny
	zg = nz
	xgr = gp(1) - xg
	ygr = gp(2) - yg
	zgr = gp(3) - zg
c
c calculate coefficients of trilinear function
c
	a8 = phimap(nx,ny,nz)
	a7 = phimap(nx,ny,nz+1) - a8
	a6 = phimap(nx,ny+1,nz) - a8
	a5 = phimap(nx+1,ny,nz) - a8
	a4 = phimap(nx,ny+1,nz+1) - a8 - a7 - a6
	a3 = phimap(nx+1,ny,nz+1) - a8 - a7 - a5
	a2 = phimap(nx+1,ny+1,nz) - a8 - a6 - a5
	a1 = phimap(nx+1,ny+1,nz+1) - a8 - a7 - a6
     &   - a5 - a4 - a3 - a2
c
c determine value of phi
c
	phi = a1*xgr*ygr*zgr + a2*xgr*ygr + a3*xgr*zgr
     &     + a4*ygr*zgr + a5*xgr + a6*ygr + a7*zgr + a8

      return
      end
