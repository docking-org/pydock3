	subroutine mapit(exrad,radprb,scale,natom,atcrd,atrad,ngrid,igrid,iepsmp,iatmap,debmap,iatmapd)
c--------------------------
c implement less scale dependent generation of eps, debye mapping
c--------------------------
	implicit none
c--------------------------
	integer natom,ngrid,igrid
	integer nspt,nsp
	parameter (nsp=9,nspt=nsp**3  )
c--------------------------
c passed arrays
	integer*1 iepsmp(ngrid,ngrid,ngrid,3)
	integer*1 debmap(ngrid,ngrid,ngrid)
	integer iatmap(ngrid,ngrid,ngrid)
	integer iatmapd(ngrid,ngrid,ngrid)
	real*4 atcrd(3,natom),atrad(natom)
	real*4 exrad,radprb,scale
c local arrays
c--------------------------
c various radii all scaled to grid scale
	real*4 satrad(natom) ! atom radius
	real*4 satradph ! atom radius+probe+0.5
	real*4 satrad2p ! atom radius+2*probe
	real*4 satradp(natom) ! atom radius+probe sqq.d
	real*4 sradmq,satradi ! sq. of atom radius-1/4, atom radius+ionic
	real*4 prad,sprad ! probe radius,radius sqd. 
	real*4 irad,sirad,spradmq ! ionic radius, atom+ionic sqd., atom+probe sqd -1/4
c linked lists
	integer ihead(igrid,igrid,igrid) ! entry into linked list of atoms mapped to each grid
	integer ilink(natom) ! linked list of atoms
c various distance stuff
	real*4 dxarr(igrid,3), dxarr2(igrid,3) ! distance arrays for grid point  x,y,z components
	real*4 xat,yat,zat
	real*4 distsq,distsq1,distsq2,distsq3,dmin2
	real roff(3,3) ! offsets to grid line mid points
	logical it1(0:6)
	integer inear(6)
c
	integer iepsmp2(igrid,igrid,igrid,3) ! stores atom numbers  temporarily
	real*4 sphcrd(3,nspt),gsphcrd(3,nspt) ! coords of canonical, translated sphere
c--------------------------
	real*4 start,fin,cputime
	real*4 xran,radmax
	real*4 pi,ang,dang,zed,dzed,crad
c--------------------------
	integer i,j,k,n,i1,j1,k1,n1
	integer i2,j2,k2,n2,inbr
	integer nin,ninter,natout,ndeb,natin,nout,nflip,nstay,ncheck
	integer llim(3),ulim(3),iloc(3),icutl,icutu
	integer iat,jat,nleft,midg,ipt,nlive,nslab
c--------------------------
	character*80 filnam
c--------------------------
	data roff / 0.5, 0.0, 0.0,   0.0, 0.5, 0.0,   0.0, 0.0, 0.5 /
c--------------------------
	pi = 355./113.
	start = cputime(0.0)
c
c clear maps, initialize values
c
	print *,' '
	print *,'in mapit, initializing arrays....'
      do k = 1,igrid
	  dxarr(k,1) = 0.
	  dxarr(k,2) = 0.
	  dxarr(k,3) = 0.
	  dxarr2(k,1) = 0.
	  dxarr2(k,2) = 0.
	  dxarr2(k,3) = 0.
        do j = 1,igrid
          do i = 1,igrid
            iepsmp(i,j,k,1) = 0
            iepsmp(i,j,k,2) = 0
            iepsmp(i,j,k,3) = 0
            iepsmp2(i,j,k,1) = 0
            iepsmp2(i,j,k,2) = 0
            iepsmp2(i,j,k,3) = 0
            debmap(i,j,k)  = 1 
            iatmap(i,j,k)  = 0 
            iatmapd(i,j,k)  = 0 
	      ihead(i,j,k) = 0
          end do
        end do
      end do
	print *,' exrad,radprb,scale ',exrad,radprb,scale ! debug
	print *,'natom, ngrid,igrid: ',natom,ngrid,igrid ! debug
	prad = radprb*scale
	irad = exrad*scale
	print *,'scaled probe ionic radii: ',prad,irad ! debug
c
c loop over all atoms
c
	print *,'generating debye map and solvent accessible volume...'
	natout = 0
	natin = 0
	radmax = -1.e6
	do n = 1,natom
c
c scale radii
c
	  satrad(n) = scale*atrad(n)
	  radmax = max(radmax,satrad(n))
	  satradi = satrad(n)+irad
	  satrad2p = satrad(n)+2.*prad
	  satradph = satrad(n)+prad+0.5
	  satradp(n) = (satrad(n)+prad)**2
c	  print *,'atm,+ionic, +2probe,+probe+0.5, +probe, sqd radii: ',satrad(n),satradi,satrad2p,satradph,satradp(n) ! debug
c
c map to nearest grid point , and set up linked list for efficient check later in reentrant volume
c generation
c
	  ilink(n) = 0
	  do k = 1,3
	    iloc(k) = nint(atcrd(k,n)) ! nearest grid point
	  end do
c	  print *,'n,iloc: ',n,iloc ! debug
	  icutl = - satrad2p
	  icutu = igrid + 2 + satrad2p
c	  print *,'icutl, icutu: ',icutl,icutu ! debug
	  natout = natout + 1
	  if((iloc(1).ge.icutl).and.(iloc(2).ge.icutl).and.(iloc(3).ge.icutl))then
	    if((iloc(1).le.icutu).and.(iloc(2).le.icutu).and.(iloc(3).le.icutu))then
		natin = natin + 1
	      natout = natout - 1
	      do k = 1,3
		  iloc(k) = max(1,iloc(k))
		  iloc(k) = min(igrid,iloc(k))
		end do
c	      print *,'after clipping n,iloc: ',n,iloc ! debug
		ilink(n) = ihead(iloc(1),iloc(2),iloc(3))
		ihead(iloc(1),iloc(2),iloc(3)) = n
	    end if
	  end if
c
c find upper, lower grid indices that enclose atom
c
	  icutl = max(satradi,satradph) + 2
	  do k = 1,3
	   i1 = nint(atcrd(k,n)) - icutl
	   i1 = max(1,i1)
	   llim(k) = min(igrid,i1)
	   i1 = nint(atcrd(k,n)) + icutl
	   i1 = max(1,i1)
	   ulim(k) = min(igrid,i1)
	  end do
c	  print *,'icutl,llim,ulim: ',icutl,llim,ulim ! debug
c
c set up x,y,z component distance arrays, probe and ionic distance cutoffs
c
	  sradmq = satrad(n)**2 - 0.25
	  sirad = satradi**2
	  spradmq = (satrad(n)+prad)**2 - 0.25
c	  print *,'vdw, probe, ionic dist. cutoffs: ',sradmq,spradmq,sirad  ! debug
	  do k = 1,3
	    xat = atcrd(k,n)
	    do i = llim(k),ulim(k)
	      dxarr(i,k) = (i-xat)
	      dxarr2(i,k) = (i-xat)**2
c		print *,'dxarr,dxarr2: ',i,xat,dxarr(i,k),dxarr2(i,k) ! debug
	    end do
	  end do
c
c loop over all grid points 
c if dist < rad+ionic iatmapd = atom #
c if dist < rad iatmap = atom #
c if dist >rad, and < rad+probe, iatmap = -atom #
c
	  do k = llim(3),ulim(3)
	    do j = llim(2),ulim(2)
	      do i = llim(1),ulim(1)
	        distsq = dxarr2(i,1) + dxarr2(j,2) + dxarr2(k,3)
		  distsq1 = distsq + dxarr(i,1)
		  distsq2 = distsq + dxarr(j,2)
		  distsq3 = distsq + dxarr(k,3)
		  if(distsq.lt.sirad)iatmapd(i,j,k) = n
		  if(distsq1.lt.sradmq)then
		    iepsmp2(i,j,k,1) = n
		  else if(distsq1.lt.spradmq)then
		    iepsmp2(i,j,k,1) = -n
		  end if
		  if(distsq2.lt.sradmq)then
		    iepsmp2(i,j,k,2) = n
		  else if(distsq2.lt.spradmq)then
		    iepsmp2(i,j,k,2) = -n
		  end if
		  if(distsq3.lt.sradmq)then
		    iepsmp2(i,j,k,3) = n
		  else if(distsq3.lt.spradmq)then
		    iepsmp2(i,j,k,3) = -n
		  end if
	      end do
	    end do
	  end do
	end do ! end of atom loop
	print *,'max radius (grids): ',radmax ! debug
	print *,'atoms excluded (out of grid), included: ',natout,natin ! debug
c
c debug- count point classes
c
	nin = 0
	nout = 0
	ninter = 0
      do k=1,igrid
        do j=1,igrid
          do i=1,igrid
            do i1=1,3
              if(iepsmp2(i,j,k,i1).gt.0) then
		    nin = nin + 1
	  	  else if(iepsmp2(i,j,k,i1).lt.0)then
		    ninter = ninter + 1
		  else
		    nout = nout + 1
		  end if
            end do
            if(iatmapd(i,j,k).ne.0) ndeb = ndeb + 1
          end do
        end do
      end do
	print *,'nin,ninter,nout,ndeb: ',nin,ninter,nout,ndeb ! debug
c
c debug- count point classes
c
c
c debug check linked list
c
	print *,'linked list i,j,k,jat: ' ! debug
	natin = 0
	do k = 1,igrid
	  do j = 1,igrid
	    do i = 1,igrid
		iat = ihead(i,j,k)
		do while(iat.ne.0)
		  natin = natin + 1
		  print *,i,j,k,iat ! debug
		  iat = ilink(iat)
		end do
	    end do
	  end do
	end do
c debug check linked list
	print *,'natout,natin: ',natout,natin ! debug
	print *,'cpu time (s) used so far in mapit: ',fin
c--------------------------
	if(radprb.le.0.)goto 101 ! skip reentrant volume generation
c--------------------------
c
c generate rentrant volume from inter zone region
c to do this we have to determine fate of ech interzone point:
c Examine possible positions of a probe center that could turn it 'out'
c These possible positions lie within a sphere probe radius centered on point
c examine a 'grid' of points in this sphere, AND are further than probe radius + atom radius 
c from all atoms.  however we only need to check atoms whose center is radius + probe diameter
c from midpoint. here we use linked list of atoms close to each grid point generated earlier
c generate canonical sphere of ns points, radius prad, centerd at 0,0,0
c
	nflip = 0
	nstay = 0
	sprad = prad**2 + 1.e-6
	midg = (nsp+1)/2
	xran = prad/(midg-1)
	print *,'sprad,midg,xran: ',sprad,midg,xran ! debug
	fin = cputime(start)
	print *,'generating probe test location sphere...'
c
c solid sphere of points
c
	goto 111 ! debug
	do i = 1,nspt
	  do k = 1,3
	    sphcrd(k,i) = 0.
	  end do
	end do
	nleft = 1
	do k = 1,nsp
	  zat = xran*(k-midg)
	  do j = 1,nsp
	    yat = xran*(j-midg)
	    do i = 1,nsp
	      xat = xran*(i-midg)
		distsq = zat**2 + yat**2 + xat**2
		if(distsq.le.sprad)then
		  nleft = nleft + 1
		  sphcrd(1,nleft) = xat
		  sphcrd(2,nleft) = yat
		  sphcrd(3,nleft) = zat
		end if
	    end do
	  end do
	end do
111	continue
c
c shell of points
c
      nslab = nsp
      dang = 2*pi/nslab
      nleft = 0
      dzed = 2./nslab
      print *,'nslab,dang,dzed: ',nslab,dang,dzed ! debug
c
      zed = -1.-dzed/2.
      do i = 1,nslab
        zed = zed+dzed
        crad = sqrt(1.-zed**2)
	  print *,'zed, crad: ',zed,crad
        do j=1,nslab
          nleft = nleft+1
          ang = (j-1)*dang
          sphcrd(1,nleft) = prad*crad*sin(ang)
          sphcrd(2,nleft) = prad*crad*cos(ang)
          sphcrd(3,nleft) = prad*zed
        end do
      end do
	print *,'# of points in probe check sphere: ',nleft 
c
c check sphere list
c
	do ipt = 1,nleft
	  print *,'ipt,xyz: ',ipt,(sphcrd(k,ipt),k=1,3) ! debug
	end do
c
c loop over all gridline points in inter zone (iepsmp2 < 0)
c translate sphere of points to gridline point
c eliminate any point closer than atrad+probe rad to any atom
c
	icutl = 2 + 2*prad + radmax
c	print *,'icutl: ',icutl ! debug
	print *,'deciding whether interzone points are in or out....'
	do k = 1,igrid
	  print *,'doing layer: ',k
	  i1 = k - icutl
	  i1 = max(1,i1)
	  llim(3) = min(igrid,i1)
	  i1 = k + icutl
	  i1 = max(1,i1)
	  ulim(3) = min(igrid,i1)
	  do j = 1,igrid
c	    print *,'doing slice: ',j
	    i1 = j - icutl
	    i1 = max(1,i1)
	    llim(2) = min(igrid,i1)
	    i1 = j + icutl
	    i1 = max(1,i1)
	    ulim(2) = min(igrid,i1)
	    do i = 1,igrid
	      do 100 n2 = 1,3
		  if(iepsmp2(i,j,k,n2).lt.0)then
c
c point in inter zone
c
c
c translate sphere of points to center on this grid point
c eliminating any closer than atrad+probe, first checking
c atom that this point 'belongs to'
c
		    iat = abs(iepsmp2(i,j,k,n2))
c		    print *,'iat: ',iat ! debug
		    nlive = 0
		    do ipt = 1,nleft
		      xat = sphcrd(1,ipt) + i + roff(1,n2)
		      yat = sphcrd(2,ipt) + j + roff(2,n2)
		      zat = sphcrd(3,ipt) + k + roff(3,n2)
c			print *,'gsph: ',xat,yat,zat ! debug
			distsq = (xat - atcrd(1,iat))**2 + (yat - atcrd(2,iat))**2 + (zat - atcrd(3,iat))**2 
			if(distsq.gt.satradp(iat))then
			  nlive = nlive + 1
			  gsphcrd(1,nlive) = xat
			  gsphcrd(2,nlive) = yat
			  gsphcrd(3,nlive) = zat
			end if
		    end do
	          i1 = i - icutl
	          i1 = max(1,i1)
	          llim(1) = min(igrid,i1)
	          i1 = i + icutl
	          i1 = max(1,i1)
	          ulim(1) = min(igrid,i1)
c		    print *,'ijk,n2,llim,ulim: ',i,j,k,n2,llim,ulim ! debug
c		    print *,'nlive: ',nlive ! debug
c		    if((nflip.gt.3).and.(nstay.gt.3))goto 300 ! debug- cut short after a few of each for check
c
c loop over all remaining sphere points
c
		    do 200 ipt = 1,nlive
		      xat = gsphcrd(1,ipt)
		      yat = gsphcrd(2,ipt)
		      zat = gsphcrd(3,ipt)
c			ncheck = 0 ! debug
	            do k2 = llim(3),ulim(3)
	              do j2 = llim(2),ulim(2)
	                do i2 = llim(1),ulim(1)
				jat = ihead(i2,j2,k2)
				do while(jat.ne.0)
c				  ncheck = ncheck + 1 ! debug
			  	  distsq = (xat - atcrd(1,jat))**2 + (yat - atcrd(2,jat))**2 + (zat - atcrd(3,jat))**2 
c				  print *,'ipt,next atom, dist: ',ipt,jat,distsq ! debug
c
c point collides, go to next sphere point
c
			  	  if(distsq.lt.satradp(jat))goto 200
		  		  jat = ilink(jat)
				end do
		          end do
		        end do
		      end do
c			print *,'# atoms checked for this point : ',ncheck ! debug
c
c if we reach here then current point survived thru all possible collision
c candidates, so it is a possible probe position, turn the midpoint out
c and go to next grid point
c
			iepsmp2(i,j,k,n2) = 0
			nflip = nflip + 1
			go to 100
200		    continue ! next sphere point
c
c if we reach here, no possible point has survived collision check with atoms
c so there is no possible probe position that could turn this midpoint out0 it must be in
c reentrant volume
c
c	          iepsmp(i,j,k,n2) = iat ! change to positive atom number
		    nstay = nstay + 1
		  end if
100	      continue ! next grid point
	    end do
	  end do
	end do
300	continue
	print *,'number of interzone points confirmed in: ',nstay
	print *,'number of interzone points flipped out: ',nflip
	fin = cputime(start)
	print *,'cpu time (s) used so far in mapit: ',fin

	
101	continue ! here if probe radius zero, skipped reentrant volume
c
c finish up maps
c set iepsmap to inside (1) wherever iespmp2 is set to any atom number,
c i.e. is inside any atom
c
	nin = 0
	nout = 0
	ninter = 0
      do k=1,igrid
        do j=1,igrid
          do i=1,igrid
            do i1=1,3
              if(iepsmp2(i,j,k,i1).gt.0) then
		    nin = nin + 1
		    iepsmp(i,j,k,i1)=1 ! puts out vdw zone
	  	  else if(iepsmp2(i,j,k,i1).lt.0)then
		    ninter = ninter + 1
		    iepsmp(i,j,k,i1)=1 ! puts out inter zone
		  else
		    nout = nout + 1
		  end if
            end do
            if(iatmapd(i,j,k).ne.0) debmap(i,j,k)=0
          end do
        end do
      end do
	print *,'nin,ninter,nout,ndeb: ',nin,ninter,nout,ndeb ! debug
c
	ninter = 0
	do i = 1,5
	  it1(i) = .true.
	end do
	it1(0) = .false.
	it1(6) = .false.
	do k = 2,igrid-1
	  do j = 2,igrid-1
	    do i = 2,igrid-1
            inbr = iepsmp(i,j,k,1) + iepsmp(i,j,k,2) + iepsmp(i,j,k,3) 
     &	     + iepsmp(i-1,j,k,1)+ iepsmp(i,j-1,k,2) + iepsmp(i,j,k-1,3) 
	      if(it1(inbr)) then
	        ninter = ninter+1
		  inear(1) = iepsmp2(i,j,k,1)
		  inear(2) = iepsmp2(i,j,k,2)
		  inear(3) = iepsmp2(i,j,k,3)
		  inear(4) = iepsmp2(i-1,j,k,1)
		  inear(5) = iepsmp2(i,j-1,k,2)
		  inear(6) = iepsmp2(i,j,k-1,3)
		  iat = 0
		  dmin2 = 1.e6
		  do i1 = 1,6
		    if((inear(i1).ge.1).and.(inear(i1).le.natom))then
		      distsq = (atcrd(1,inear(i1))-i)**2 + (atcrd(2,inear(i1))-j)**2 + (atcrd(3,inear(i1))-k)**2
		      if(distsq.lt.dmin2)then
			  dmin2 = distsq
			  iat = inear(i1)
			end if
		    end if
		  end do
		  iatmap(i,j,k) = iat
	      end if
	    end do
	  end do
	end do
	print *,'# of grid points on molecular surface= ',ninter
c
c dump map out for check
c
	filnam = 'mapit.map'
	call dmpeps(filnam)
c
	fin = cputime(start)
	print *,'total mapit cpu time (s): ',fin
c	stop
	end
