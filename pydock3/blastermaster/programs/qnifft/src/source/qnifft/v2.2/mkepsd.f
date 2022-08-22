	subroutine mkepsd(natom,exrad,radprb)
c
c implement anthony nicholls algorithm for making eps, debye map
c as in qdiffxs
c cleaned up version kas 14-feb-01
c
c------------------------------------------------------------------
	include 'qdiffpar.h'
	parameter (nsp=500000)
	integer ibgrd(3,nsp)
	integer iepsmp2(ngrid,ngrid,ngrid,3)
	dimension xo(3),xn(3)
	logical   it1(0:6),ihs
	integer iline(65)
	integer inear(6)
	character*80 filnam
c------------------------------------------------------------------
c
c set eps and deb map values to their exterior values
c prior to resetting interior values using atom coords and radii
c
	start = cputime(0.)
	write(6,*)'   '
	write(6,*)'initializing dielectric and debye maps...'
	write(6,*)'   '
	do k = 1,igrid
	  do j = 1,igrid
	    do i = 1,igrid
		iepsmp(i,j,k,1) = 0
		iepsmp(i,j,k,2) = 0
		iepsmp(i,j,k,3) = 0
		iepsmp2(i,j,k,1) = 0
		iepsmp2(i,j,k,2) = 0
		iepsmp2(i,j,k,3) = 0
		debmap(i,j,k)  = 1 ! kas 8-feb-01
		iatmap(i,j,k)  = 0 ! kas 14-feb-01
		iatmapd(i,j,k)  = 0 ! kas 7-may-01
	    end do
	  end do
	end do
c
c set in all midline (eps) & grid points within expanded (solvent accessible) surface
c
	call setout(atcrd,atrad,natom,iepsmp2,iatmapd,igrid,ngrid,
     &exrad,radprb,scale)
c
c set iepsmap to inside (1) wherever iespmp2 is set to any atom number,
c i.e. is inside any atom
c
	finish=cputime(start)
	write(6,*) 'time to turn everything out is',finish 
	do iz=1,igrid
	  do iy=1,igrid
	    do ix=1,igrid
	      do i=1,3
	        if(iepsmp2(ix,iy,iz,i).ne.0) iepsmp(ix,iy,iz,i)=1
	      end do
	      if(iatmapd(ix,iy,iz).ne.0) debmap(ix,iy,iz)=0
	    end do
	  end do
	end do
c
c check for boundary gridpoints (any that are neighbored by both in and out midlines)
c and store indices
c
	it1(0)= .false.
	it1(6)= .false.
	do ix=1,5
	  it1(ix)= .true.
	end do
	ibnum=0
	do k=2,igrid-1
	  do j=2,igrid-1
	    do i=2,igrid-1
            ieps=   iepsmp(i,j,k,1) + iepsmp(i,j,k,2) +
     &      iepsmp(i,j,k,3) + iepsmp(i-1,j,k,1)+ 
     &      iepsmp(i,j-1,k,2) + iepsmp(i,j,k-1,3) 
	      if(it1(ieps)) then
	        ibnum=ibnum+1
	        if(ibnum.gt.nsp)then
		    print *,'ibnum in mkmaps: ',ibnum
          print *,'max # eps bndry pts exceeded, increase nsp: ',nsp
	          stop
	        end if
	        ibgrd(1,ibnum)=i
	        ibgrd(2,ibnum)=j
	        ibgrd(3,ibnum)=k
	      end if
	    end do
	  end do
	end do
	write(6,*) "number of grid points on expanded surface= ",ibnum
c
c finish off dielectric map, producing  molecular surface
c volume defined by prob radius
c by resetting to outside any midlines that lie within solvent probe radius
c placing probe on each boundary gridpoint in turn. iepsmp2 contains atom #
c for each midline, so boundary gridpoints can be moved to exact atom surface
c before use
c
	call fleps(radprb,ibnum,ibgrd,nsp,atrad,atcrd,natom,igrid,ngrid,
     &scale,iepsmp,iepsmp2)
c
c have used iepsmp2, now reassign iepsmp2 so it contains atom # for
c inside midlines that are within vdw volume.  These atom numbers can be
c used to reposition surface gridpoints exactly on convex (vdw) part of
c molecular surface later in phierg
c
	finish=cputime(start)
	write(6,*) 'redo insides at ',finish 
c
	call setin(atcrd,atrad,natom,iepsmp2,igrid,ngrid,scale)
c
	ibnum = 0
	do k = 2,igrid-1
	  do j = 2,igrid-1
	    do i = 2,igrid-1
            ieps =    iepsmp(i,j,k,1) + iepsmp(i,j,k,2) +
     &      iepsmp(i,j,k,3) + iepsmp(i-1,j,k,1)+ 
     &      iepsmp(i,j-1,k,2) + iepsmp(i,j,k-1,3) 
c
c iepsmp2 contains indices of atoms forming surface
c store index position of nearest atom in iatmap
c
	      if(it1(ieps)) then
	        ibnum = ibnum+1
		  inear(1) = iepsmp2(i,j,k,1)
		  inear(2) = iepsmp2(i,j,k,2)
		  inear(3) = iepsmp2(i,j,k,3)
		  inear(4) = iepsmp2(i-1,j,k,1)
		  inear(5) = iepsmp2(i,j-1,k,2)
		  inear(6) = iepsmp2(i,j,k-1,3)
		  iat = 0
		  dmin2 = 1.e6
		  srfx = i
		  srfy = j
		  srfz = k
		  do i1 = 1,6
		    if((inear(i1).ge.1).and.(inear(i1).le.natom))then
		      dist2 = (atcrd(1,inear(i1))-srfx)**2 
     &            + (atcrd(2,inear(i1))-srfy)**2 
     &            + (atcrd(3,inear(i1))-srfz)**2
		      if(dist2.lt.dmin2)then
			  dmin2 = dist2
			  iat = inear(i1)
			end if
		    end if
		  end do
		  iatmap(i,j,k) = iat
	      end if
	    end do
	  end do
	end do
	write(6,*) "number of grid points on molecular surface= ",ibnum
	finish = cputime(start)
	write(6,*) 'time to turn everything in is',finish 
c
c debug
c
c	filnam = 'dump19.map'
c	call dmpeps(filnam)
c
c debug
c
	return
	end

	subroutine setout(atcrd,atrad,natom,iepsmp2,iatmapd,igrid,ngrid,
     &exrad,radprb,scale)
	parameter (nsp=500000)
	parameter (lmtb=15)
	dimension atcrd(3,natom),atrad(natom),xn(3)
	integer iepsmp2(ngrid,ngrid,ngrid,3)
	integer iatmapd(ngrid,ngrid,ngrid) ! kas 8-feb-01
	integer ioff(3,nsp),ismin(3),ismax(3)
	dimension sq1(-lmtb:lmtb),sq2(-lmtb:lmtb),sq3(-lmtb:lmtb)
	dimension rad2a1(-lmtb:lmtb),rad2a2(-lmtb:lmtb),rad2a3(-lmtb:lmtb)
	logical itobig,itest2
c
	itobig = .false. ! kas 1-feb-01 init itobig
c
c find largest possible extended radius
c
	radmax2 = 0.0
	do 691 ix = 1,natom
	  radmax2 = amax1(radmax2,atrad(ix))
691	continue
	temp = amax1(exrad,radprb)
	radmax2 = scale*(radmax2+temp)
	lim = 1 + radmax2
c 
	if(lim.gt.lmtb-1) itobig = .true.
	print *,'lim in setout: ',lim
c 
	if(itobig) goto 7878 
c
c set up list (ioff) of grid indices that lie within possible
c hitting distance of any atom's extended radius- only these need to be check
c later
c
	radtest = (radmax2 + 0.5*sqrt(3.0))**2
	ibox = 0
	do 692 ix = -lim,lim
	  do 693 iy = -lim,lim
	    do 694 iz = -lim,lim
	    dist = ix**2 + iy**2 + iz**2 
	    dist1 = dist + ix + 0.25
	    dist2 = dist + iy + 0.25
	    dist3 = dist + iz + 0.25
	    itest = 0
	    if(dist.lt.radtest) itest = 1
	    if(dist1.lt.radtest) itest = 1
	    if(dist2.lt.radtest) itest = 1
	    if(dist3.lt.radtest) itest = 1 
	    if(itest.eq.1) then 
	    ibox = ibox+1
	    if(ibox.gt.nsp)then
		print *,'ibox in setout: ',ibox
	      print *,'max # of eps pts in check list IOFF exceeded, increase nsp: ',nsp
	      stop
	    end if
	    ioff(1,ibox) = ix
	    ioff(2,ibox) = iy
	    ioff(3,ibox) = iz
	    end if
694	continue
693	continue
692	continue
7878	continue 
	print *,'ibox in setout: ',ibox
c
c loop over all atoms
c 
	num = 0
	do 607 iv = 1, natom
c
c extract coordinates and radius 
c
	  rad =  atrad(iv)
	  xn(1) = atcrd(1,iv)
	  xn(2) = atcrd(2,iv)
	  xn(3) = atcrd(3,iv)
	  if(rad.le.0.) goto 608
c
c scale radius to grid, and calc. various check radii needed later
c
	  rad = rad*scale
	  radp = rad + exrad*scale
	  rad = rad + radprb*scale
	  rad2 = rad*rad
	  radp2 = radp*radp
c
c find upper and lower grid index limits
c that fully enclose sphere of atom
c and ensure they lie within grid
c
c
	  itest2=.false.
        do 9017 k = 1,3
          ismin(k) = (xn(k) - radmax2 - 1.)
	    itest1=ismin(k)
	    ismin(k) = min(ismin(k),igrid)
	    ismin(k) = max(ismin(k),1)
	    if(itest1.ne.ismin(k)) itest2=.true.
          ismax(k) = (xn(k) + radmax2 + 1.)
	    itest1=ismax(k)
	    ismax(k) = min(ismax(k),igrid)
	    ismax(k) = max(ismax(k),1)
	    if(itest1.ne.ismax(k)) itest2=.true.
9017	  continue
c
c
	  if(itest2.or.itobig) then
c
c CASE 1: extended atom radius touches edge of grid
c loop trough all grid indices within limits
c if within Ratom+Rprobe then set iepsmp2 = atom #
c if within Ratom+Rionic then set iatmapd to atom number
c note code uses efficient computation of distances to 3 midlines
c
	    num=num+1
	    rad2a = rad2 - 0.25
          do 9019 iz =  ismin(3),ismax(3)
            do 9020 iy =  ismin(2),ismax(2)
              do 9021 ix =  ismin(1),ismax(1)
                dxyz1 = (ix - xn(1))
                dxyz2 = (iy - xn(2))
                dxyz3 = (iz - xn(3))
                distsq = dxyz1**2 +dxyz2**2 +dxyz3**2 
		    distsq1 = distsq + dxyz1
		    distsq2 = distsq + dxyz2
		    distsq3 = distsq + dxyz3
		    if(distsq1.lt.rad2a) iepsmp2(ix,iy,iz,1)=iv
		    if(distsq2.lt.rad2a) iepsmp2(ix,iy,iz,2)=iv
		    if(distsq3.lt.rad2a) iepsmp2(ix,iy,iz,3)=iv
		    if(distsq.lt.radp2) iatmapd(ix,iy,iz)=iv ! kas 4-may-01
9021		  continue
9020	      continue
9019	    continue
	  else
c
c CASE 2: extended atom radius totaly within grid
c check only grid points stored in IOFF
c precompute x,y,z direction distances to possible grid points (sq1/2/3), 
c and rad2a1/2/3, assuming nearest grid point to atom is at 0,0,0,
c then simple add offsets ix,iy,iz
c if within Ratom+Rprobe then set iepsmp2 = atom #
c if within Ratom+Rionic then set iatmapd = atom #
c
	    rad2a = rad2 - 0.25
	    ixn1=nint(xn(1))
	    iyn1=nint(xn(2))
	    izn1=nint(xn(3))
	    fxn1=ixn1-xn(1)
	    fxn2=iyn1-xn(2)
	    fxn3=izn1-xn(3)
	    do 6020 ix=-lim,lim
	      temp1=ix+fxn1
	      temp2=ix+fxn2
	      temp3=ix+fxn3
	      sq1(ix)=temp1*temp1
	      sq2(ix)=temp2*temp2
	      sq3(ix)=temp3*temp3
	      rad2a1(ix)=rad2a-temp1
	      rad2a2(ix)=rad2a-temp2
	      rad2a3(ix)=rad2a-temp3
6020	    continue 
C$DIR NO_RECURRENCE
	    do 9024 i=1,ibox
	      i1= ioff(1,i)
	      i2= ioff(2,i)
	      i3= ioff(3,i)
	      ix=ixn1+ i1
	      iy=iyn1+ i2
	      iz=izn1+ i3
            distsq = sq1(i1) +sq2(i2) + sq3(i3)
	      if(distsq.lt.rad2a1(i1)) iepsmp2(ix,iy,iz,1)=iv
	      if(distsq.lt.rad2a2(i2)) iepsmp2(ix,iy,iz,2)=iv
	      if(distsq.lt.rad2a3(i3)) iepsmp2(ix,iy,iz,3)=iv
            if(distsq.lt.radp2)   iatmapd(ix,iy,iz)= iv ! kas 8-feb-01
9024	    continue
	  end if
608	  continue ! radius = 0.
607	continue ! end of atom loop
	print *,'atoms in case 1 (num): ',num
c
	return
	end

	subroutine setin(atcrd,atrad,natom,iepsmp2,igrid,ngrid,scale)
	parameter (nsp=500000)
	parameter (lmtb=15)
	dimension atcrd(3,natom),atrad(natom),xn(3)
	integer iepsmp2(ngrid,ngrid,ngrid,3)
	integer ioff(3,nsp),ismin(3),ismax(3)
	dimension sq1(-lmtb:lmtb),sq2(-lmtb:lmtb),sq3(-lmtb:lmtb)
	dimension rad2a1(-lmtb:lmtb),rad2a2(-lmtb:lmtb),rad2a3(-lmtb:lmtb)
	logical itobig,itest2
c-----------------------------------------------------------
	itobig = .false. ! kas 1-feb-01 init itobig
c
c need to regenerate iepsmp2, index list of atom #'s that the midline point is 
c within, so that in reaction field energy calc. in phierg, the point can be moved
c to exact position on atom surface
c
c reset iepsmp2
c
	do i = 1,3
	  do iz = 1,igrid
	    do iy = 1,igrid
	      do ix = 1,igrid
	        iepsmp2(ix,iy,iz,i) = 0
	      end do
	    end do
	  end do
	end do
c
c
c find largest possible extended radius
c
	radmax2 = 0.0
	do 691 ix = 1,natom
	  radmax2 = amax1(radmax2,atrad(ix))
691	continue
	radmax2 = scale*radmax2
	lim = 1+ radmax2
c 
	if(lim.gt.lmtb-1) itobig=.true.
	print *,'lim in setin: ',lim
c 
	if(itobig) goto 7878 
c
c set up list (ioff) of grid indices that lie within possible
c hitting distance of any atom's extended radius- only these need to be check
c later
c
	radtest= (radmax2 + 0.5*sqrt(3.0))**2
	ibox=0
	do 692 ix=-lim,lim
	  do 693 iy=-lim,lim
	    do 694 iz=-lim,lim
	      dist=ix**2 + iy**2 + iz**2 
	      dist1=dist + ix + 0.25
	      dist2=dist + iy + 0.25
	      dist3=dist + iz + 0.25
	      itest=0
	      if(dist.lt.radtest) itest=1
	      if(dist1.lt.radtest) itest=1
	      if(dist2.lt.radtest) itest=1
	      if(dist3.lt.radtest) itest=1 
	      if(itest.eq.1) then 
	        ibox=ibox+1
	        if(ibox.gt.nsp)then
		    print *,'ibox in setin: ',ibox
	          print *,'max # of eps pts in check list IOFF exceeded, increase nsp: ',nsp
	          stop
	        end if
	        ioff(1,ibox)=ix
	        ioff(2,ibox)=iy
	        ioff(3,ibox)=iz
	      end if
694	    continue
693	  continue
692	continue
7878	continue 
	print *,'ibox in setin: ',ibox
c
c loop over all atoms
c 
	num = 0
	do 607 iv=1, natom
c
c extract coordinates and radius 
c
	  rad= atrad(iv)
	  xn(1)=atcrd(1,iv)
	  xn(2)=atcrd(2,iv)
	  xn(3)=atcrd(3,iv)
	  if(rad.le.0.) goto 608
c
c scale radius to grid, and calc. various check radii needed later
c
	  rad = rad*scale
	  rad2 = rad*rad
c
c set dielectric map
c
c find upper and lower grid index limits
c that fully enclose sphere of atom
c ensure they lie within grid
c
c
	  itest2=.false.
        do 9017 k = 1,3
          ismin(k) = (xn(k) - radmax2 - 1.)
	    itest1=ismin(k)
	    ismin(k) = min(ismin(k),igrid)
	    ismin(k) = max(ismin(k),1)
	    if(itest1.ne.ismin(k)) itest2=.true.
          ismax(k) = (xn(k) + radmax2 + 1.)
	    itest1=ismax(k)
	    ismax(k) = min(ismax(k),igrid)
	    ismax(k) = max(ismax(k),1)
	    if(itest1.ne.ismax(k)) itest2=.true.
9017	  continue
c
c
	  if(itest2.or.itobig) then
c
c CASE 1: atom radius touches edge of grid
c loop trough all grid indices within limits
c if within Ratom+Rprobe then set iepsmp2 = atom #
c note code uses efficient computation of distances to 3 midlines
c
	    num=num+1
	    rad2a = rad2 - 0.25
          do 9019 iz =  ismin(3),ismax(3)
            do 9020 iy =  ismin(2),ismax(2)
              do 9021 ix =  ismin(1),ismax(1)
                dxyz1 = (ix - xn(1))
                dxyz2 = (iy - xn(2))
                dxyz3 = (iz - xn(3))
                distsq = dxyz1**2 +dxyz2**2 +dxyz3**2 
		    distsq1 = distsq + dxyz1
		    distsq2 = distsq + dxyz2
		    distsq3 = distsq + dxyz3
		    if(distsq1.lt.rad2a) iepsmp2(ix,iy,iz,1)=iv
		    if(distsq1.lt.rad2a) iepsmp2(ix,iy,iz,1)=iv
		    if(distsq2.lt.rad2a) iepsmp2(ix,iy,iz,2)=iv
9021		  continue
9020	      continue
9019	    continue
	  else
c
c CASE 2: atom radius totaly within grid
c check only grid points stored in IOFF
c precompute x,y,z direction distances to possible grid points (sq1/2/3), 
c and rad2a1/2/3, assuming nearest grid point to atom is at 0,0,0,
c then simple add offsets ix,iy,iz
c if within Ratom+Rprobe then set iepsmp2 = atom #
c
	    rad2a = rad2 - 0.25
	    ixn1=nint(xn(1))
	    iyn1=nint(xn(2))
	    izn1=nint(xn(3))
	    fxn1=ixn1-xn(1)
	    fxn2=iyn1-xn(2)
	    fxn3=izn1-xn(3)
	    do 6020 ix=-lim,lim
	      temp1=ix+fxn1
	      temp2=ix+fxn2
	      temp3=ix+fxn3
	      sq1(ix)=temp1*temp1
	      sq2(ix)=temp2*temp2
	      sq3(ix)=temp3*temp3
	      rad2a1(ix)=rad2a-temp1
	      rad2a2(ix)=rad2a-temp2
	      rad2a3(ix)=rad2a-temp3
6020	    continue 
C$DIR NO_RECURRENCE
	    do 9024 i=1,ibox
	      i1= ioff(1,i)
	      i2= ioff(2,i)
	      i3= ioff(3,i)
	      ix=ixn1+ i1
	      iy=iyn1+ i2
	      iz=izn1+ i3
            distsq = sq1(i1) +sq2(i2) + sq3(i3)
	      if(distsq.lt.rad2a1(i1)) iepsmp2(ix,iy,iz,1)=iv
	      if(distsq.lt.rad2a2(i2)) iepsmp2(ix,iy,iz,2)=iv
	      if(distsq.lt.rad2a3(i3)) iepsmp2(ix,iy,iz,3)=iv
9024	    continue
	  end if
608	  continue ! radius = 0.
607	continue ! end of atom loop
	print *,'atoms in case 1 (num): ',num
c
	return
	end

	subroutine fleps(radprb,ibnum,ibgrd,nsp,atrad,atcrd,natom,igrid,ngrid,
     &scale,iepsmp,iepsmp2)
c
c fills up intersticies between  vanderWaals atom spheres
c that would otherwise be given a dielectric constant of water:
c if probe radius, radprb, not zero, then finds
c every inside grid point that has any outside eps
c points as neighbours, (ie all boundary points)
c and turns all eps points within
c probe radius back to outside dielectric- this generates an eps map
c that approximates the solvent accessible surface
c
	dimension ibgrd(3,nsp),atcrd(3,natom),atrad(natom),mv(6)
	integer*1 iepsmp(ngrid,ngrid,ngrid,3)
	integer iepsmp2(ngrid,ngrid,ngrid,3)
	real*4 dx2(igrid),dy2(igrid),dz2(igrid)
	real*4 dx(igrid),dy(igrid),dz(igrid)
c
c----------------------------------------------------------------
c----------------------------------------------------------------
c
	if(radprb.eq.0.0) return
	write(6,*)' finishing off dielectric map...'
	kas = 1
	if(kas.eq.1)print *,' using kas code in fleps...'
c
c limits of expanded surface
c
	ix1=igrid
	iy1=igrid
	iz1=igrid
	ix65=1
	iy65=1
	iz65=1
	do 543 i=1,ibnum
	  ib1=ibgrd(1,i)
	  ib2=ibgrd(2,i)
	  ib3=ibgrd(3,i)
	  ix1=min(ib1,ix1)
	  iy1=min(ib2,iy1)
	  iz1=min(ib3,iz1)
	  ix65=max(ib1,ix65)
	  iy65=max(ib2,iy65)
	  iz65=max(ib3,iz65)
543	continue
	minp=min(ix1,iy1,iz1)
	maxp=max(ix65,iy65,iz65)
	write(6,*) "expanded surface maximum, minimum (g.u.)=",minp,maxp
c
c scale probe radius for checks
c
	radp = radprb*scale
	radp2 = radp**2
	lim=int(radp+1.5)
	print *,'lim in fleps: ',lim
	rad2a = radp2 - 0.25
c
	ivz=0
	ivz1 = 0
c
c loop over all surface gridpoints
c
	print *,'in fleps ibnum: ',ibnum
	do 50 i=1,ibnum
c
c coordinates
c
	  ib1=ibgrd(1,i)
	  ib2=ibgrd(2,i)
	  ib3=ibgrd(3,i)
	  xn=float(ib1)
	  yn=float(ib2)
	  zn=float(ib3)
c
c atom #'s of neighboring midlines
c
	  mv(1)=iepsmp2(ib1,ib2,ib3,1)
	  mv(2)=iepsmp2(ib1,ib2,ib3,2)
	  mv(3)=iepsmp2(ib1,ib2,ib3,3)
	  mv(4)=iepsmp2(ib1-1,ib2,ib3,1)
	  mv(5)=iepsmp2(ib1,ib2-1,ib3,2)
	  mv(6)=iepsmp2(ib1,ib2,ib3-1,3)
c
c find closest atom
c
	  dism=1.0e6
	  iv=0
	  do j=1,6
	    m=mv(j)
	    if((m.ne.0).and.(m.ne.iv)) then
	      dist=(atcrd(1,m)-xn)**2 + (atcrd(2,m)-yn)**2 
     &      + (atcrd(3,m)-zn)**2
	      if(dist.lt.dism) then
	        iv=m
	        dism=dist
	      end if
	    end if
	  end do
	  if(iv.eq.0) then
c
c count surface points for which a nearest atom isn't found-
c only > 0 if an error somewhere
c
	    ivz=ivz+1
	  else
c
c move surface G'point to exactly on atom's extended surface
c
	    xa=atcrd(1,iv)
	    ya=atcrd(2,iv)
	    za=atcrd(3,iv)
c
	    arad=(xn-xa)**2 + (yn-ya)**2 + (zn-za)**2
	    rad=(atrad(iv)*scale+radp)
	    sfact=rad/(sqrt(arad))
c
	    xn=xa + sfact*(xn-xa)
	    yn=ya + sfact*(yn-ya)
	    zn=za + sfact*(zn-za)
	    ib1=int(xn)
	    ib2=int(yn)
	    ib3=int(zn)
	  end if
c
c set grid index loop limits , and reset every midline < Rprobe from
c moved atom to ouside.
c
	  lim1=max0(ib1-lim,1)
	  lim2=min0(ib1+lim,igrid)
	  lim3=max0(ib2-lim,1)
	  lim4=min0(ib2+lim,igrid)
	  lim5=max0(ib3-lim,1)
	  lim6=min0(ib3+lim,igrid)
c
	  z=float(lim5)-1.0
	  y1=float(lim3)-1.0
	  x1=float(lim1)-1.0
c
	  if(kas.eq.0)then
        do 9019 iz = lim5,lim6 
	    z=z+1.0
	    y=y1
          do 9020 iy = lim3,lim4 
	      y=y+1.0
	      x=x1
            do 9021 ix = lim1,lim2 
	        x=x+1.0
              dxyz1 = (x - xn)
              dxyz2 = (y - yn)
              dxyz3 = (z - zn)
              distsq = dxyz1**2 +dxyz2**2 +dxyz3**2 
		  distsq1 = distsq + dxyz1
		  distsq2 = distsq + dxyz2
		  distsq3 = distsq + dxyz3
		  if(distsq1.lt.rad2a) iepsmp(ix,iy,iz,1)=0
		  if(distsq2.lt.rad2a) iepsmp(ix,iy,iz,2)=0
		  if(distsq3.lt.rad2a) iepsmp(ix,iy,iz,3)=0
9021	      continue
9020	    continue
9019	  continue
c kas code 14-feb-01
	  else
	  do ix = lim1,lim2
	    dx(ix) = rad2a - (ix-xn)
	    dx2(ix) = (ix-xn)**2
	  end do
	  do iy = lim3,lim4
	    dy(iy) = rad2a - (iy-yn)
	    dy2(iy) = (iy-yn)**2
	  end do
	  do iz = lim5,lim6
	    dz(iz) = rad2a - (iz-zn)
	    dz2(iz) = (iz-zn)**2
	  end do
        do iz = lim5,lim6 
          do iy = lim3,lim4 
            do ix = lim1,lim2 
		  dr2 = dx2(ix)+dy2(iy)+dz2(iz)
		  if(dr2.lt.dx(ix)) iepsmp(ix,iy,iz,1)=0
		  if(dr2.lt.dy(iy)) iepsmp(ix,iy,iz,2)=0
		  if(dr2.lt.dz(iz)) iepsmp(ix,iy,iz,3)=0
		end do
	    end do
	  end do
	  end if
c kas code 14-feb-01
50	continue ! end of atom loop
	if(ivz.ne.0) write(6,*) "no. of surface points unassigned (ivz)= ",ivz
	if(ivz1.ne.0) write(6,*) "no. of surface points unassigned (ivz1)= ",ivz1
	return
	end

	subroutine debmemb(zin,zout)
	include 'qdiffpar.h'
c
c set ionic strength in inner region (zin< z < zout)
c and outer region (z > zout) to rionsti, rionsto (debfcti,debfcto)
c respectively
c
	goff = (igrid+1.)/2.
	do iz = 1,igrid
	  zz = (iz - goff)/scale + oldmid(3)
c	  print *,'left: ',iz
	  if(zz.gt.zin)then
	    if(zz.gt.zout)then
c	      print *,'right: ',iz
	      do iy = 1,igrid
	        do ix = 1,igrid
		    if(debmap(ix,iy,iz).ne.0)then
			debmap(ix,iy,iz) = 3 ! kas 8-feb-01
		    end if
	        end do
	      end do
	    else
c	      print *,'middle: ',iz
	      do iy = 1,igrid
	        do ix = 1,igrid
		    if(debmap(ix,iy,iz).ne.0)then
			debmap(ix,iy,iz) = 2 ! kas 8-feb-01
		    end if
	        end do
	      end do
	    end if
	  end if
	end do
	return
	end
