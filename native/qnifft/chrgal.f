      subroutine chrgal(xo,gchrgv,xrad,nout)
c
c     assign fractional charges to grid points using 
c     anti-aliasing using sphere to approximate subcube
c     and determining common volume in lieu of counting
c     subcubes within atom
c     fractional assignment added for atoms which overlap 
c     grid edge 4/16/96
c
      include 'qdiffpar.h'
      dimension xo(3)
      logical outgr
      data pi / 3.1415929 /
c---------------------------------------------------

c     set radius of sphere which is to be centered
c     around grid point and calculate volume of sphere
      sprad=0.62
      spvol=(4./3.)*pi*(sprad**3)

c     find grid coordinate range that lies within radius
c     plus one grid point of atom
      ixup=nint(xo(1)+xrad+1)
      ixlow=nint(xo(1)-xrad-1)
      iyup=nint(xo(2)+xrad+1)
      iylow=nint(xo(2)-xrad-1)
      izup=nint(xo(3)+xrad+1)
      izlow=nint(xo(3)-xrad-1)

c     initialize variables needed for section dealing
c     with atoms falling outside of the grid
      total=0.
      outgr=.false.
      vingr=1.0
      pin=0.
      ptotal=0.

c     check to make sure range falls within grid
c     if in grid, proceed; if not, count atom out
c     of grid, reset range to 1 or ngrid, and determine
c     the fraction of in-range grid points which fall
c     within grid
      nxlow=ixlow
      nxup=ixup
      nylow=iylow
      nyup=iyup
      nzlow=izlow
      nzup=izup

      if (ixup.gt.ngrid) then
       nxup=ngrid
       outgr=.true.
      else if (ixup.lt.1) then
       nxup=1
       outgr=.true.
      end if
      if (iyup.gt.ngrid) then
       nyup=ngrid
       outgr=.true.
      else if (iyup.lt.1) then
       nyup=1
       outgr=.true.
      end if
      if (izup.gt.ngrid) then
       nzup=ngrid
       outgr=.true.
      else if (izup.lt.1) then
       nzup=1
       outgr=.true.
      end if
      if (ixlow.lt.1) then
       nxlow=1
       outgr=.true.
      else if (ixlow.gt.ngrid) then
       nxlow=ngrid
       outgr=.true.
      end if
      if (iylow.lt.1) then
       nylow=1
       outgr=.true.
      else if (iylow.gt.ngrid) then
       nylow=ngrid
       outgr=.true.
      end if
      if (izlow.lt.1) then
       nzlow=1
       outgr=.true.
      else if (izlow.gt.ngrid) then
       nzlow=ngrid
       outgr=.true.
      end if
     
      if (outgr) then
       nout=nout+1
       do i=ixlow,ixup
        do j=iylow,iyup
         do k=izlow,izup
          if (i.lt.1.or.i.gt.ngrid) then
          else if (j.lt.1.or.j.gt.ngrid) then
          else if (k.lt.1.or.k.gt.ngrid) then
          else   
           pin=pin+1
          end if
          ptotal=ptotal+1
         end do
        end do
       end do
       vingr=(pin/ptotal)
      end if

c     scale atom charge by fraction of grid points
c     assigned to atom which fall in the grid
      gchrgv=gchrgv*vingr

c     intialize volume sum for entire atom to zero
       vfsum=0.
      
c     clear phimap array to store in-atom volumes for each grid point
        do iphi=ixlow,ixup
         do jphi=iylow,iyup
          do kphi=izlow,izup
           phimap(iphi,jphi,kphi)=0.
          end do
         end do
        end do

c     loop over x, y, and z grid coordinates, calculating common
c     volume between sphere surrounding grid point and atom
      do i=nxlow, nxup
       do j=nylow, nyup
        do k=nzlow, nzup

c     calculate distance between grid point and atom center
      dsx=((xo(1)-i)**2+(xo(2)-j)**2+(xo(3)-k)**2)**0.5
      rsum=xrad+sprad
      rdiff=abs(xrad-sprad)

      if (dsx.gt.rsum) then
       phimap(i,j,k)=0.
      else if (dsx.le.rdiff) then
       phimap(i,j,k)=1.0
      else

c     determine distance from center of sphere and center of
c     atom to center of the circle of intersection(COI)
      dx=((xrad**2)-(sprad**2)+(dsx**2))/(2.*dsx)
      ds=((sprad**2)-(xrad**2)+(dsx**2))/(2.*dsx)

c     determine height from center of COI to buried surfaces
      hx=xrad-dx
      hs=sprad-ds

c     calculate volume of each hemisphere with COI as base, sum
c     volumes, and store fraction of volume that lies within atom 
c     in phimap array
      vx=pi*(((hx**2)*xrad)-((hx**3)/3.))
      vs=pi*(((hs**2)*sprad)-((hs**3)/3.))
      phimap(i,j,k)=(vx+vs)/spvol
  
      end if

c     keep running sum of volume fractions for entire atom
       vfsum=vfsum+phimap(i,j,k)

c     end loops
        end do
       end do
      end do

c     loop through grid points and assign fractional charges
c     store fractions in phimap and charges in qmap
      do i=ixlow, ixup
       do j=iylow, iyup
        do k=izlow, izup
         if (vfsum.ne.0) then
           phimap(i,j,k)=phimap(i,j,k)/vfsum
           qmap(i,j,k)=qmap(i,j,k)+phimap(i,j,k)*gchrgv
         end if
        end do
       end do
      end do

      return
      end
