	subroutine chrgit(charge,chrgv)
c
c assign charge at each corner of box to a grid volume
c and assign fractional charge to grid points of this volume
c to preserve as many electrical moments as possible
c using simple linear weighting scheme
c
	include 'qdiffpar.h'
c---------------------------------------------------

c
      dimension kb(3),kt(3),dg(3),cg(3)      
	dimension charge(3)
c
c place charge if inside grid
c
      iout = 0
	rgrid = igrid
      do 9000 j = 1,3
        if((charge(j).le.1.).or.(charge(j).ge.rgrid)) iout = 1
c      end do
9000	continue
      if(iout.eq.0) then
        do 9001 i = 1,3
c
c trucate to nearest grid point
c
          kb(i) = charge(i)
c
c find upper grid point of box
c
          kt(i) = kb(i) + 1
c
c find position of charge in box in terms
c of fractional distance along edges
c
          dg(i) = charge(i) - kb(i)
          cg(i) = 1-dg(i)
c
c keep within grid- 
c
          if(kb(i).lt.1) kb(i) = 1
          if(kt(i).gt.igrid) kt(i) = igrid
c        end do
9001    continue
c
c put fractional charges at 8 corners
c
        qmap(kb(1),kb(2),kb(3)) = qmap(kb(1),kb(2),kb(3)) + 
     &  cg(1)*cg(2)*cg(3)*chrgv
        qmap(kt(1),kb(2),kb(3)) = qmap(kt(1),kb(2),kb(3)) + 
     &  dg(1)*cg(2)*cg(3)*chrgv
        qmap(kb(1),kt(2),kb(3)) = qmap(kb(1),kt(2),kb(3)) +
     &  cg(1)*dg(2)*cg(3)*chrgv
        qmap(kt(1),kt(2),kb(3)) = qmap(kt(1),kt(2),kb(3)) +
     &  dg(1)*dg(2)*cg(3)*chrgv
        qmap(kb(1),kb(2),kt(3)) = qmap(kb(1),kb(2),kt(3)) +
     &  cg(1)*cg(2)*dg(3)*chrgv
        qmap(kt(1),kb(2),kt(3)) = qmap(kt(1),kb(2),kt(3)) +
     &  dg(1)*cg(2)*dg(3)*chrgv
        qmap(kb(1),kt(2),kt(3)) = qmap(kb(1),kt(2),kt(3)) +
     &  cg(1)*dg(2)*dg(3)*chrgv
        qmap(kt(1),kt(2),kt(3)) = qmap(kt(1),kt(2),kt(3)) +
     &  dg(1)*dg(2)*dg(3)*chrgv
	end if	!if inside grid

	return
	end
c--------------------------------------------------
	subroutine chrgit1(charge,chrgv)
c
c assign charge to grid points using mike gilsons
c scheme which preserves spherical symmetry better
c
	include 'qdiffpar.h'
c--------------------------------------------------
c
c parameters controlling size and resolution of sphere within which
c grid points assign charge will lie- determined empirically by mkg to
c give good accuracy for interactions of point charges in spheres
c
	parameter (crit = .0001)
	data radmin,radmax,radstep / 1.2, 2.0, .005 /
c--------------------------------------------------
	dimension ddist2(500), rsave(500)
	dimension clist(4,50,500), nlist(500)
	logical flag
      dimension kb(3),dg(3)      
	dimension charge(3)
c--------------------------------------------------
c
c place charge if inside grid
c
      iout = 0
	rgrid = igrid
      do 9000 j = 1,3
        if((charge(j).le.1.).or.(charge(j).ge.ngrid)) iout = 1
c      end do
9000	continue
      if(iout.eq.0) then
        do 9001 i = 1,3
c
c trucate to nearest grid point
c
          kb(i) = charge(i)
c
c find position of charge in box in terms
c of fractional distance along edges
c
          dg(i) = charge(i) - kb(i)
c        end do
9001	  continue
c
c loop over increasing radius to find suitable one
c
	  irmax = 0
	  do 9002 rmax = radmin, radmax, radstep
	    rmax2 = rmax*rmax
	    irmax = irmax + 1
	    rsave (irmax) = rmax
	    ilist = 0
	    ir = rmax + dg(1) + 1
c
c loop over all grid points within range of rmax:
c place base point of grid box containing the charge, at the origin, for now.
c
	    do 9003 i = -ir, ir
	      xdist = i - dg(1)
	      do 9004 j = -ir, ir
	        ydist = j - dg(2)
	        do 9005 k = -ir, ir
	          zdist = k - dg(3)
c
c calculate distance from current grid point to charge:
c
	          dist2 =  xdist*xdist + ydist*ydist + zdist*zdist
	          dist = sqrt(dist2)
c
c if grid point is closer than R, store the point and its distance:
c
		    if (dist2.le. Rmax2) then
	     		ilist = ilist + 1
	     		clist(1,ilist, irmax)  = i
	     		clist(2,ilist, irmax) = j
	     		clist(3,ilist, irmax) = k
	     		clist(4,ilist, irmax) = dist
	          endif
c	        end do
c	      end do
c	    end do
9005		  continue
9004        continue
9003      continue
	    nlist(irmax) = ilist
c
c generate weighting factors for this rmax:
c sum normalizes the weighting to one:
c
	    sum = 0
	    do 9006 il = 1, nlist(irmax)
	      clist(4,il, irmax) = Rmax - clist(4,il, irmax)
	      sum = sum + clist(4,il, irmax)
c	    end do
9006	    continue
c
c normalize the weighting factors to sum to one:
c
	    do 9007 il = 1, nlist(irmax)
	      clist(4,il, irmax) = clist(4,il, irmax)/sum
c	    end do
9007	    continue
c
c calculate center of charge for this rmax:
c
	    xsum = 0
	    ysum = 0
	    zsum = 0
	    do 9008 il = 1, nlist(irmax)
	      xsum = xsum + clist(4,il, irmax) * clist(1,il, irmax)
	      ysum = ysum + clist(4,il, irmax) * clist(2,il, irmax)
	      zsum = zsum + clist(4,il, irmax) * clist(3,il, irmax)
c	    end do
9008	    continue
c	  print *, 'center of charge in map:' , xsum,ysum, zsum
c	  print *, 'actual location:        ', dg(1), dg(2), dg(3)
c
c check whether criterion is satisfied, and if so, exit rmax loop.
c
	    ddx = dg(1) - xsum
	    ddy = dg(2) - ysum
	    ddz = dg(3) - zsum
	    ddist2(irmax) = ddx * ddx + ddy*ddy + ddz * ddz
	    if (ddist2(irmax) .le.crit) goto 1000
c
c otherwise, try another cutoff radius:		     
c
c	  end do
9002    continue
c
c if loop gets finished without a radius yielding good enough
c results, print warning,and use the best cutoff radius found:
c	print *, 'Criterion not satisfied for charge #, ', iq
c	print *, ' Using best cutoff radius found...'
	  call minimum (ddist2, irmax, ichoice, dmin)
	  dmin = sqrt (dmin)
	  rmax = rsave(ichoice)
	  goto 1500

1000	  continue
	  ichoice = irmax
	  dmin = sqrt(ddist2(ichoice))

1500	  continue

c
c now we know what set of grids we're distributing the charge over
c (clist(1-3, 1-ilist(ichoice), ichoice) for ichoice), and the weighting
c (clist(4,1-ilist, ichoice)...
c now, distribute the charge:
c
	  do 9009 ilist = 1, nlist(ichoice)
c
c get grid point by adding offset to it
c
	    kx = kb(1) + clist(1,ilist,ichoice)
	    ky = kb(2) + clist(2,ilist,ichoice)
	    kz = kb(3) + clist(3,ilist,ichoice)
c
c make sure grid point is within the big box (should be a problem only
c in cases where box edge cuts through
c or very near the protein):
c
	    flag = .false.
   	    if (kx.lt.1 ) then
	      kx = 1
	      flag = .true.
	    endif
   	    if (ky .lt.1) then
	      ky = 1
	      flag = .true.
	    endif
   	    if (kz .lt.1) then
	      kz = 1
	      flag = .true.
	    endif
   	    if (kx.gt.ngrid) then
	      kx = ngrid
	      flag = .true.
	    endif
   	    if (ky.gt.ngrid)then
	      ky = ngrid
	      flag = .true.
	    endif
   	    if (kz.gt.ngrid) then
	      kz = ngrid
	      flag = .true.
	    endif
  
	    if(flag)then
	      write(6,*)' problem for charge at', (charge(j),j=1,3)
	    end if

	    qmap(kx,ky,kz)= qmap(kx,ky,kz) + clist(4,ilist,ichoice)*chrgv

c	  end do
9009	  continue
	end if

	return
	end
c---------------------------------------------------------
	subroutine minimum (ddist2, irmax, ichoice, dmin)
	dimension ddist2(500)

	dmin = 100000
	do 9000 i = 1, irmax
	  if (ddist2(i).lt.dmin) then
	    dmin = ddist2(i)
	    ichoice = i
	  endif
c	end do
9000	continue
	return
	end

c---------------------------------------------------
