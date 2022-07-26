	subroutine mkmaps(natom,exrad,radprb,debfct,epsin,epsout)
c
c implement anthony nicholls algorithm for making eps, debye map
c as in qdiffxs
c
c------------------------------------------------------------------
	include 'qdiffpar.h'
	parameter (nsp=100000)
	dimension ibgrd(4,nsp)
	dimension iatmap(ngrid,ngrid,ngrid)
	dimension iatmap2(ngrid,ngrid,ngrid)
	dimension xo(3),xn(3)
	logical   it1(0:6),ihs
	dimension iline(65)
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
		debmap(i,j,k)  = debfct
	    end do
	  end do
	end do
c
	call setout(atcrd,atrad,natom,iepsmp2,debmap,igrid,ngrid,
     &exrad,radprb,scale,iatmap,iatmap2)
c
c finish off dielectric map, producing solvent accessible
c volume defined by prob radius
c
	finish=cputime(start)
	write(6,*) 'time to turn everything out is',finish 
c
c check for boundary elements
c
	do i=1,3
	do iz=1,igrid
	do iy=1,igrid
	do ix=1,igrid
	if(iepsmp2(ix,iy,iz,i).ne.0) iepsmp(ix,iy,iz,i)=1
	end do
	end do
	end do
	end do
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
          print *,'max # eps bndry pts exceeded, increase nps: ',nsp
	          stop
	        end if
	        ibgrd(1,ibnum)=i
	        ibgrd(2,ibnum)=j
	        ibgrd(3,ibnum)=k
	        ibgrd(4,ibnum)=iatmap(i,j,k)
	      end if
	    end do
	  end do
	end do
	write(6,*) "number of grid points on expanded surface= ",ibnum
c
c reset iepsmap using these boundary elements
c
c
	ihs = .false.
	if(ihs) then
	  write(6,*) "writing expanded surface data file: hsurf.dat"
	  open(40,file="hsurf.dat")
	  sahen=1.2079/(scale**2)
	  do i=1,ibnum
	    xo(1)=ibgrd(1,i)
	    xo(2)=ibgrd(2,i)
	    xo(3)=ibgrd(3,i)
	    j=ibgrd(4,i)
	    call gtoc(xo,xn)
	    write(40,*) j,xn,sahen
	  end do
	  close(40)
	end if
414	format(i5,3f8.3,f8.3)
c
	call fleps(radprb,ibnum,ibgrd,atrad,atcrd,natom,igrid,ngrid,
     &scale,iepsmp,iepsmp2)
c
c
c have used iepsmp2, now reassign iepsmp2
c
	finish=cputime(start)
	write(6,*) 'redo insides at ',finish 
c
	call setin(atcrd,atrad,natom,iepsmp2,igrid,ngrid,
     &exrad,radprb,scale,iatmap,iatmap2)
c
	finish=cputime(start)
	write(6,*) 'time to turn everything in is',finish 
c 
	  do k=1,igrid
	    do j=1,igrid
	      do i=1,igrid
		  do id = 1,3
		    epsmap(i,j,k,id) = iepsmp(i,j,k,id)*epsin + 
     $		    (1-iepsmp(i,j,k,id))*epsout
		  end do
	      end do
	    end do
	  end do
c	midg = (igrid+1)/2
c	print *,'iepsmap: '
c	print *,' '
c	do j = 1,65
c	 iline(j) = 0
c	end do
c	do i = 1,igrid
c	    do j = 1,igrid
c		iline(j) = iepsmp(i,j,midg,1)
c	    end do
c	    write(6,'(x,65i1)')iline
c	end do
c	print *,' '
	end

	subroutine setout(atcrd,atrad,natom,iepsmp2,debmap,igrid,ngrid,
     &exrad,radprb,scale,iatmap,iatmap2)
	parameter (nsp=100000)
	dimension atcrd(3,natom),atrad(natom),xn(3)
	dimension iepsmp2(ngrid,ngrid,ngrid,3)
	real*4 debmap(ngrid,ngrid,ngrid)
	dimension iatmap(ngrid,ngrid,ngrid),iatmap2(ngrid,ngrid,ngrid)
	integer ioff(3,nsp),ismin(3),ismax(3)
	dimension sq1(-15:15),sq2(-15:15),sq3(-15:15)
	dimension rad2a1(-15:15),rad2a2(-15:15),rad2a3(-15:15)
	logical itobig,itest2
c
	radmax2=0.0
	do 691 ix=1,natom
	radmax2=amax1(radmax2,atrad(ix))
691	continue
	temp=amax1(exrad,radprb)
	radmax2=scale*(radmax2+temp)
	lim= 1+ radmax2
c 
	limmax = 12
	if(lim.gt.limmax) itobig=.true.
c 
	if(itobig) goto 7878 
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
		print *,'ibox in setout: ',ibox
	      print *,'max # eps bndry pts exceeded, increase nps: ',nsp
	      stop
	    end if
	    ioff(1,ibox)=ix
	    ioff(2,ibox)=iy
	    ioff(3,ibox)=iz
	    end if
694	continue
693	continue
692	continue
7878	continue 
c
c 
c set interiors
c
c read data file
c
	do 607 iv=1, natom
c
c restore values
c
	rad= atrad(iv)
	xn(1)=atcrd(1,iv)
	xn(2)=atcrd(2,iv)
	xn(3)=atcrd(3,iv)
	if(rad.le.0.) goto 608
c
c scale radius to grid
c
	  iv10=10*iv+1
	  rad = rad*scale
	  rad5= (rad + 0.5)**2
	  radp = rad + exrad*scale
	  rad = rad + radprb*scale
	  rad4= (rad + 0.5)**2
	  rad2 = rad*rad
	  radp2 = radp*radp
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
9017	continue
c
c
	if(itest2.or.itobig) then
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
		if(distsq.lt.radp2)  debmap(ix,iy,iz) =0.
		if(distsq.lt.rad4)iatmap(ix,iy,iz)=iatmap(ix,iy,iz)+iv10
		if(distsq.lt.rad5)iatmap2(ix,iy,iz)=iatmap2(ix,iy,iz)+iv10
9021		continue
9020	    continue
9019	  continue
		else
	rad2a = rad2 - 0.25
	ixn1=nint(xn(1))
	iyn1=nint(xn(2))
	izn1=nint(xn(3))
	fxn1=ixn1-xn(1)
	fxn2=iyn1-xn(2)
	fxn3=izn1-xn(3)
	rad2ax=rad2a-fxn1
	rad2ay=rad2a-fxn2
	rad2az=rad2a-fxn3
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
6020	continue 
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
        if(distsq.lt.radp2)   debmap(ix,iy,iz)=0.
	if(distsq.lt.rad4)    iatmap(ix,iy,iz)=iatmap(ix,iy,iz)+iv10
	if(distsq.lt.rad5)    iatmap2(ix,iy,iz)=iatmap2(ix,iy,iz)+iv10
9024	continue
		end if
c
608	continue
c
607	continue
c
	return
	end

	subroutine setin(atcrd,atrad,natom,iepsmp2,igrid,ngrid,
     &exrad,radprb,scale,iatmap,iatmap2)
	parameter (nsp=100000)
	dimension atcrd(3,natom),atrad(natom),xn(3)
	dimension iepsmp2(ngrid,ngrid,ngrid,3)
	dimension iatmap(ngrid,ngrid,ngrid),iatmap2(ngrid,ngrid,ngrid)
	integer ioff(3,nsp),ismin(3),ismax(3)
	dimension sq1(-15:15),sq2(-15:15),sq3(-15:15)
	dimension rad2a1(-15:15),rad2a2(-15:15),rad2a3(-15:15)
	logical itobig,itest2
c
	do i=1,3
	do iz=1,igrid
	do iy=1,igrid
	do ix=1,igrid
	iepsmp2(ix,iy,iz,i)=0
	end do
	end do
	end do
	end do
c
	radmax2=0.0
	do 691 ix=1,natom
	radmax2=amax1(radmax2,atrad(ix))
691	continue
	radmax2=scale*radmax2
	lim= 1+ radmax2
c 
	limmax = 12
	if(lim.gt.limmax) itobig=.true.
c 
	if(itobig) goto 7878 
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
		print *,'ibox in setout: ',ibox
	      print *,'max # eps bndry pts exceeded, increase nps: ',nsp
	      stop
	    end if
	    ioff(1,ibox)=ix
	    ioff(2,ibox)=iy
	    ioff(3,ibox)=iz
	    end if
694	continue
693	continue
692	continue
7878	continue 
c
c 
c set interiors
c
c read data file
c
	do 607 iv=1, natom
c
c restore values
c
	rad= atrad(iv)
	xn(1)=atcrd(1,iv)
	xn(2)=atcrd(2,iv)
	xn(3)=atcrd(3,iv)
	if(rad.le.0.) goto 608
c
c scale radius to grid
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
9017	continue
c
c
	if(itest2.or.itobig) then
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
9021		continue
9020	    continue
9019	  continue
		else
	rad2a = rad2 - 0.25
	ixn1=nint(xn(1))
	iyn1=nint(xn(2))
	izn1=nint(xn(3))
	fxn1=ixn1-xn(1)
	fxn2=iyn1-xn(2)
	fxn3=izn1-xn(3)
	rad2ax=rad2a-fxn1
	rad2ay=rad2a-fxn2
	rad2az=rad2a-fxn3
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
6020	continue 
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
9024	continue
		end if
c
608	continue
c
607	continue
c
c
	return
	end

	subroutine fleps(radprb,ibnum,ibgrd,atrad,atcrd,natom,igrid,ngrid,
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
	parameter (nsp=100000)
	dimension ibgrd(4,nsp),atcrd(3,natom),atrad(natom),mv(6)
	dimension iepsmp(ngrid,ngrid,ngrid,3)
	dimension iepsmp2(ngrid,ngrid,ngrid,3)
c
c----------------------------------------------------------------
c----------------------------------------------------------------
c
	if(radprb.eq.0.0) return
	write(6,*)' finishing off dielectric map...'
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
c
	radp = radprb*scale
	radp2 = radp**2
	radp4= (radp-0.5)**2
	rgrid = igrid
	lim=int(radp+1.5)
c
	rad2a = radp2 - 0.25
	ivz=0
c
	do 50 i=1,ibnum
c
	ib1=ibgrd(1,i)
	ib2=ibgrd(2,i)
	ib3=ibgrd(3,i)
	xn=float(ib1)
	yn=float(ib2)
	zn=float(ib3)
c
	mv(1)=iepsmp2(ib1,ib2,ib3,1)
	mv(2)=iepsmp2(ib1,ib2,ib3,2)
	mv(3)=iepsmp2(ib1,ib2,ib3,3)
	mv(4)=iepsmp2(ib1-1,ib2,ib3,1)
	mv(5)=iepsmp2(ib1,ib2-1,ib3,2)
	mv(6)=iepsmp2(ib1,ib2,ib3-1,3)
	dism=1000.
	iat=100000
	do j=1,6
	  m=mv(j)
	  if((m.ne.0).and.(m.ne.iat)) then
	    dist=(atcrd(1,m)-xn)**2 + (atcrd(2,m)-yn)**2 
     &    + (atcrd(3,m)-zn)**2
	    if(dist.lt.dism) then
	      iat=m
	    dism=dist
	    end if
	  end if
	end do
	iv=iat
c
	if(iv.eq.0) then
	  ivz=ivz+1
	else
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
c set loop limits 
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
9021	    continue
9020	  continue
9019	continue
50	continue
	if(ivz.ne.0) write(6,*) "no. of surface points unassigned= ",ivz
	return
	end

	subroutine debmemb(debfct,debfcti,debfcto,zin,zout)
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
		    if(debmap(ix,iy,iz).ne.0.)then
			debmap(ix,iy,iz) = debfcto
		    end if
	        end do
	      end do
	    else
c	      print *,'middle: ',iz
	      do iy = 1,igrid
	        do ix = 1,igrid
		    if(debmap(ix,iy,iz).ne.0.)then
			debmap(ix,iy,iz) = debfcti
		    end if
	        end do
	      end do
	    end if
	  end if
	end do
	return
	end
