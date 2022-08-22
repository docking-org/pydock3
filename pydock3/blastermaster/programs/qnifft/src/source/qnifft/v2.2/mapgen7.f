	subroutine mapgen(exrad,radprb,natom)
c--------------------------
c implement less scale dependent generation of eps, debye mapping
c--------------------------
c	implicit none
	integer ngirdmx
	parameter (ngirdmx=300000)
	include 'qdiffpar.h'
c--------------------------
c	integer natom,ngrid,igrid
c--------------------------
c	real*4 exrad,radprb,scale
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
c various distance stuff
	real*4 dxarr(ngrid,3), dxarr2(ngrid,3) ! distance arrays for grid point  x,y,z components
	real*4 xat,yat,zat,xyz(3)
	real*4 distsq,distsq1,distsq2,distsq3,dmin2,dist
	real roff(3,3) ! offsets to grid line mid points
	logical it1(0:6)
	integer inear(6)
c
	integer iepsmp2(ngrid,ngrid,ngrid,3) ! stores atom numbers  temporarily
	integer ihead(ngrid,ngrid,ngrid) ! flag for interzone point/entry into linked list of girdle points
	integer ilink(ngirdmx) ! linked list of girdle points
	real*4  gxyz(3,ngirdmx) ! linked list of girdle points
	real*4 srat,ratio,cratmx,cratmn,rratmx,rratmn,cdstmx,cdstmn,rdstmx,rdstmn
	real*4 cdstmx1,cdstmn1,slop

c--------------------------
	real*4 start,fin,cputime
	real*4 xran,radmax
	real*4 pi,ang,dang,zed,dzed,crad
c--------------------------
	integer i,j,k,n,i1,j1,k1,n1
	integer i2,j2,k2,n2,inbr,lunit
	integer nin,ninter,ndeb,nout,ncheck,ncont,nreent,ncont1
	integer llim(3),ulim(3),icutl,icutu
	integer iat,jat,nleft,midg,ipt
	integer npoint,nbad,ngood,ngird,nlook
	integer iloc(3)
c--------------------------
	character*80 filnam
	character*10 head
c--------------------------
	data roff / 0.5, 0.0, 0.0,   0.0, 0.5, 0.0,   0.0, 0.0, 0.5 /
c--------------------------
	pi = 355./113.
	start = cputime(0.0)
c
c clear maps, initialize values
c
	print *,' '
	print *,'in mapgen, initializing arrays....'
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
	radmax = -1.e6
	do n = 1,natom
c
c scale radii
c
	  satrad(n) = scale*atrad(n)
	  radmax = max(radmax,satrad(n))
	  satradi = satrad(n)+irad
	  satradph = satrad(n)+prad+0.5
	  satradp(n) = (satrad(n)+prad)
c	  print *,'atom xyz: ',atcrd(1,n),atcrd(2,n),atcrd(3,n),satradp(n) ! debug
c	  print *,'atm,+ionic, +2probe,+probe+0.5, +probe, sqd radii: ',satrad(n),satradi,satrad2p,satradph,satradp(n) ! debug
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
		    if(iepsmp2(i,j,k,1).eq.0)iepsmp2(i,j,k,1) = -n
c		    iepsmp2(i,j,k,1) = -n ! old debug- will overwrite i's vdw with j's inter
		  end if
		  if(distsq2.lt.sradmq)then
		    iepsmp2(i,j,k,2) = n
		  else if(distsq2.lt.spradmq)then
		    if(iepsmp2(i,j,k,2).eq.0)iepsmp2(i,j,k,2) = -n
c		    iepsmp2(i,j,k,2) = -n ! old debug
		  end if
		  if(distsq3.lt.sradmq)then
		    iepsmp2(i,j,k,3) = n
		  else if(distsq3.lt.spradmq)then
		    if(iepsmp2(i,j,k,3).eq.0)iepsmp2(i,j,k,3) = -n
c		    iepsmp2(i,j,k,3) = -n ! old debug
		  end if
	      end do
	    end do
	  end do
	end do ! end of atom loop
	print *,'max radius (grids): ',radmax ! debug
c
c debug- count point classes
c
	nin = 0
	nout = 0
	ninter = 0
	ndeb = 0
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
	fin = cputime(start)
	print *,'cpu time (s) used so far in mapgen: ',fin
c--------------------------
	if(radprb.le.0.)goto 101 ! skip reentrant volume generation
c--------------------------
c
c generate rentrant volume 
c first generate SAS points
c then read them in and turn everything within 1 probe rad out again
c
	lunit=40
c	open(unit=lunit,file='qnifft_sas.usr',status='scratch')
	open(lunit,file='qnifft_sas.usr',status='unknown')
	call sasgen(natom,atcrd,satradp,prad,scale,npoint,lunit,ngird)

	fin = cputime(start)
	print *,'cpu time (s) used so far in mapgen: ',fin
	print *,'reading sas points from qnifft_sas.usr...'
	print *,'expecting # of points: ',npoint
c	open(lunit,file='qnifft_sas.usr',status='old')
	rewind(lunit)
	read(lunit,'(a)')head
	print *,head
c 
c repeat procedure for setting SAS volume in, but use SAS points instead
c of atoms, and set anything within probe rad of a SAS point out
c 
c
	icutl = prad + 2.5
	sradmq = prad**2 - 0.25
	print *,'probe cutoffs: ',icutl,sradmq  ! debug
	ncheck = npoint/10 + 1
	nbad = 0
	ngood = 0
	do n = 1,npoint
	  if(mod(n,ncheck).eq.1)then
	    print *,'point ',n
	  end if
	  read(lunit,223,end=901,err=902)xyz,head
c	  call pntchk(natom,atcrd,satradp,prad,scale,xyz,nbad,ngood) ! debug
c
c find upper, lower grid indices that enclose point + prad
c
	  do k = 1,3
	   i1 = nint(xyz(k)) - icutl
	   i1 = max(1,i1)
	   llim(k) = min(igrid,i1)
	   i1 = nint(xyz(k)) + icutl
	   i1 = max(1,i1)
	   ulim(k) = min(igrid,i1)
	  end do
c	  print *,'icutl,llim,ulim: ',icutl,llim,ulim ! debug
c
c set up x,y,z component distance arrays, probe and ionic distance cutoffs
c
	  do k = 1,3
	    xat = xyz(k)
	    do i = llim(k),ulim(k)
	      dxarr(i,k) = (i-xat)
	      dxarr2(i,k) = (i-xat)**2
c		print *,'dxarr,dxarr2: ',i,xat,dxarr(i,k),dxarr2(i,k) ! debug
	    end do
	  end do
c
c loop over all grid points 
c if dist < rad probe iepsmp = 0
c
	  do k = llim(3),ulim(3)
	    do j = llim(2),ulim(2)
	      do i = llim(1),ulim(1)
	        distsq = dxarr2(i,1) + dxarr2(j,2) + dxarr2(k,3)
		  if(iepsmp2(i,j,k,1).lt.0)then ! interzone midpoint(s) at 1 
		    distsq1 = distsq + dxarr(i,1)
		    if(distsq1.lt.sradmq)then
		      iepsmp2(i,j,k,1) = 0
		    end if
	        end if
		  if(iepsmp2(i,j,k,2).lt.0)then ! interzone midpoint(s) at 1 
		    distsq2 = distsq + dxarr(j,2)
		    if(distsq2.lt.sradmq)then
		      iepsmp2(i,j,k,2) = 0
		    end if
		  end if
		  if(iepsmp2(i,j,k,3).lt.0)then ! interzone midpoint(s) at 1 
		    distsq3 = distsq + dxarr(k,3)
		    if(distsq3.lt.sradmq)then
		      iepsmp2(i,j,k,3) = 0
		    end if
		  end if
	      end do
	    end do
	  end do
	end do ! end of point loop
c	print *,'# of bad, good points: ',nbad,ngood
223   format(3f12.7,a4)
c
c map all girdle points to a linked list
c
	print *,'making linked list of girdle points...'
	rewind(lunit)
	read(lunit,'(a)')head
c	if(ngird.gt.ngirdmx)then
	if(npoint.gt.ngirdmx)then
c	  print *,'WARNING: increase size of girdle point list from, to: ',ngirdmx,ngird
	  print *,'WARNING: increase size of girdle point list from, to: ',ngirdmx,npoint
c	  ngird = ngirdmx
	  npoint = ngirdmx
	end if
c	do n = 1,ngird
	do n = 1,npoint
	  read(lunit,223,end=901,err=902)xyz,head
        ilink(n) = 0
        do k = 1,3
          iloc(k) = nint(xyz(k)) ! nearest grid point
          iloc(k) = max(1,iloc(k))
          iloc(k) = min(igrid,iloc(k))
	    gxyz(k,n) = xyz(k)
        end do
        ilink(n) = ihead(iloc(1),iloc(2),iloc(3))
        ihead(iloc(1),iloc(2),iloc(3)) = n
	end do
	close(lunit)
cc
cc debug check linked list
cc
c      print *,'linked list i,j,k,jat: ' ! debug
c	ngood = 0
c      do k = 1,igrid
c        do j = 1,igrid
c          do i = 1,igrid
c            iat = ihead(i,j,k)
c            do while(iat.ne.0)
c		  ngood = ngood + 1
cc              print *,i,j,k,iat ! debug
c              iat = ilink(iat)
c            end do
c          end do
c        end do
c      end do
cc	print *,'points in girdle,linked list: ',ngird,ngood
c	print *,'points in girdle,linked list: ',npoint,ngood


	fin = cputime(start)
	print *,'cpu time (s) used so far in mapgen: ',fin
	
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
c find boundary points, and assign to an atom if possible, or a girdle point
c
	icutl = prad + 2
	print *,'icutl for girdle search: ',icutl
c	print *,'mapping boundary points to nearest atom/probe center....'
	print *,'mapping boundary points to nearest probe or atom surface....'
	nreen = 0
	nbnd = 0
	ncont = 0
	ncont1 = 0
	nbad = 0
	cratmn = 1.e6
	cratmx = -1.e6
	rratmn = 1.e6
	rratmx = -1.e6
	cdstmn = 1.e6
	cdstmx = -1.e6
	rdstmn = 1.e6
	rdstmx = -1.e6
	cdstmn1 = 1.e6
	cdstmx1 = -1.e6
	slop = 0.5
	print *,'boundary point assignment slop (grids): ',slop
	do i = 1,5
	  it1(i) = .true.
	end do
	it1(0) = .false.
	it1(6) = .false.
c
	do k = 2,igrid-1
	  do j = 2,igrid-1
	    do i = 2,igrid-1
            inbr = iepsmp(i,j,k,1) + iepsmp(i,j,k,2) + iepsmp(i,j,k,3) 
     &	     + iepsmp(i-1,j,k,1)+ iepsmp(i,j-1,k,2) + iepsmp(i,j,k-1,3) 
	      if(it1(inbr)) then
c
c found boundary point, check for nearby atom
c
	        nbnd = nbnd+1
              xyzbnd(1,nbnd) = i
              xyzbnd(2,nbnd) = j
              xyzbnd(3,nbnd) = k
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
c nearest atoms surface
			dist = sqrt(distsq) - satrad(inear(i1))
			distsq = dist**2
c nearest atoms surface
c nearest atom center
		      if(distsq.lt.dmin2)then
			  dmin2 = distsq
			  iat = inear(i1)
			end if
		    end if
		  end do
c search linked list for closest girdle point
c
c		  print *,'finding girdle for orphan point ',i,j,k ! debug
c		  print *,'inear: ',inear
		  ipt = 0
c		  if(ipt.eq.0)then
		  if(iat.eq.0)then
		    dmin2 = 1.e6
	          iloc(1) = i
	          iloc(2) = j
	          iloc(3) = k
	          do k1 = 1,3
                  i1 = iloc(k1) - icutl
                  i1 = max(1,i1)
                  llim(k1) = min(igrid,i1)
                  i1 = iloc(k1) + icutl
                  i1 = max(1,i1)
                  ulim(k1) = min(igrid,i1)
c			llim(k1) = 1 ! debug
c			ulim(k1) = igrid ! debug
	          end do
c		    print *,'searching over ',llim,ulim
		    nlook = 0
                do k2 = llim(3),ulim(3)
                  do j2 = llim(2),ulim(2)
                    do i2 = llim(1),ulim(1)
                      jat = ihead(i2,j2,k2)
                      do while(jat.ne.0)
				nlook = nlook + 1
                        distsq = (i - gxyz(1,jat))**2 + (j - gxyz(2,jat))**2 + (k - gxyz(3,jat))**2
cc nearest probe surface
c				dist = sqrt(distsq) - prad
c				distsq = dist**2
cc			      print *,'looking at point: ',jat,distsq ! debug
cc nearest probe surface
c nearest probe center
                        if(distsq.lt.dmin2)then
		              dmin2 = distsq
				  dist = sqrt(distsq) - prad
		              ipt = jat
		            end if
              		jat = ilink(jat)
	                end do
	              end do
	            end do
	          end do
c		    print *,'found ',i,j,k,' at ',ipt,dmin2
c		    print *,'after searching ',nlook
		  end if
c		  if((ipt.ne.0))then
		  if((ipt.ne.0).and.(abs(dist).lt.slop))then
c closest is girdle point 
c scale grid position to be exactly 1 probe rad from it
                dmin2 = (gxyz(1,ipt)-i)**2 + (gxyz(2,ipt)-j)**2 + (gxyz(3,ipt)-k)**2
                if(dmin2.gt.1.e-6)then
		      nreen = nreen + 1
                  srat = prad/sqrt(dmin2)
			dist = sqrt(dmin2) - prad
c			print *,'srat: ',srat
c			if(srat.lt.0.8)stop
cc
c                 srat = 1. ! debug
cc
                  rratmx = max(rratmx,srat)
                  rratmn = min(rratmn,srat)
			if(dist.gt.rdstmx)then
c			  print *,'new max dist error: ',dist,i,j,k
c			  print *,inear
c			  print *,ipt,gxyz(1,ipt),gxyz(2,ipt),gxyz(3,ipt)
			end if
                  rdstmx = max(rdstmx,dist)
                  rdstmn = min(rdstmn,dist)
                  xyzbnd(1,nbnd) = gxyz(1,ipt) + srat*(i-gxyz(1,ipt)) + 0. ! debug omitt reentrant from rxn
                  xyzbnd(2,nbnd) = gxyz(2,ipt) + srat*(j-gxyz(2,ipt)) + 0. ! debug omitt reentrant from rxn
                  xyzbnd(3,nbnd) = gxyz(3,ipt) + srat*(k-gxyz(3,ipt)) + 0. ! debug omitt reentrant from rxn
                end if
		  else if(iat.ne.0)then
c closest is atom
c scale grid position to be exactly 1 atom rad from it
c
		    iatmap(i,j,k) = iat
                dmin2 = (atcrd(1,iat)-i)**2 + (atcrd(2,iat)-j)**2 + (atcrd(3,iat)-k)**2
                if(dmin2.gt.1.e-6)then
		      ncont = ncont + 1
                  srat = satrad(iat)/sqrt(dmin2)
			dist = sqrt(dmin2) - satrad(iat)
cc
c                 srat = 1. ! debug
cc
		      cratmx = max(cratmx,srat)
		      cratmn = min(cratmn,srat)
                  cdstmx = max(cdstmx,dist)
                  cdstmn = min(cdstmn,dist)
                  xyzbnd(1,nbnd) = atcrd(1,iat) + srat*(i-atcrd(1,iat)) + 0. ! debug omitt reentrant from rxn
                  xyzbnd(2,nbnd) = atcrd(2,iat) + srat*(j-atcrd(2,iat)) + 0. ! debug omitt reentrant from rxn
                  xyzbnd(3,nbnd) = atcrd(3,iat) + srat*(k-atcrd(3,iat)) + 0. ! debug omitt reentrant from rxn
		    end if
		  else
		    nbad = nbad + 1
c		    stop
		    iat = 0
		    dmin2 = 1.e6
cc just look at atom(s) that produced neigboring interzone midpoints (inear < 0)
c		    do i1 = 1,6
c		      if((-inear(i1).ge.1).and.(-inear(i1).le.natom))then
c		        distsq = (atcrd(1,-inear(i1))-i)**2 + (atcrd(2,-inear(i1))-j)**2 + (atcrd(3,-inear(i1))-k)**2
cc nearest atoms surface
c			  dist = sqrt(distsq) - satrad(-inear(i1))
c			  distsq = dist**2
cc nearest atoms surface
cc nearest atom center
c		        if(distsq.lt.dmin2)then
c			    dmin2 = distsq
c			    iat = -inear(i1)
c			  end if
c		      end if
c		    end do
cc end just look at atom(s) that produced neigboring interzone midpoints (inear < 0)
c look at all atoms
		    do i1 = 1,natom
		      distsq = (atcrd(1,i1)-i)**2 + (atcrd(2,i1)-j)**2 + (atcrd(3,i1)-k)**2
c nearest atoms surface
			dist = sqrt(distsq) - satrad(i1)
			distsq = dist**2
c nearest atoms surface
c nearest atom center
		      if(distsq.lt.dmin2)then
			  dmin2 = distsq
			  iat = i1
			end if
		    end do
c look at all atoms
		    if(iat.ne.0)then
c		      iepsmp2(i,j,k,1)   = natom+1 ! debug put oint in vdw zone for display
c		      iepsmp2(i,j,k,2)   = natom+1 ! debug put oint in vdw zone for display
c		      iepsmp2(i,j,k,3)   = natom+1 ! debug put oint in vdw zone for display
c		      iepsmp2(i-1,j,k,1) = natom+1 ! debug put oint in vdw zone for display
c		      iepsmp2(i,j-1,k,2) = natom+1 ! debug put oint in vdw zone for display
c		      iepsmp2(i,j,k-1,3) = natom+1 ! debug put oint in vdw zone for display
cc closest is atom
cc scale grid position to be exactly 1 atom rad from it
cc
c			print *,'ijk: ',i,j,k
c			print *,'iat: ',iat,inear,sqrt(dmin2)
                  dmin2 = (atcrd(1,iat)-i)**2 + (atcrd(2,iat)-j)**2 + (atcrd(3,iat)-k)**2
                  if(dmin2.gt.1.e-6)then
		        ncont1 = ncont1 + 1
                    srat = satrad(iat)/sqrt(dmin2)
			  dist = sqrt(dmin2) - satrad(iat)
cc
c                   srat = 1. ! debug
cc
                    cdstmx1 = max(cdstmx1,dist)
                    cdstmn1 = min(cdstmn1,dist)
                    xyzbnd(1,nbnd) = atcrd(1,iat) + srat*(i-atcrd(1,iat))
                    xyzbnd(2,nbnd) = atcrd(2,iat) + srat*(j-atcrd(2,iat))
                    xyzbnd(3,nbnd) = atcrd(3,iat) + srat*(k-atcrd(3,iat))
		      end if
		    end if
		  end if ! has/not atom neighbor branch
	      end if  ! is boundary point branch
	    end do
	  end do
	end do
	print *,'# of grid points on molecular surface= ',nbnd
	print *,'# of contact molecular surface points assigned to an atom = ',ncont
	print *,'worst contact surface distance error ratios: ',cratmn,cratmx
	print *,'worst contact surface distance errors: ',cdstmn,cdstmx
	print *,'# of reentrant molecular surface points assigned to girdle = ',nreen
	print *,'worst reentrant surface distance error ratios: ',rratmn,rratmx
	print *,'worst reentrant surface distance errors: ',rdstmn,rdstmx
	print *,'# of molecular surface points unassigned = ',nbad
	print *,'# of unassigned molecular surface points assigned to an atom = ',ncont1
	print *,'worst reassigned surface distance errors: ',cdstmn1,cdstmx1
cc debug- put out only problem spots, vdw or reentrant parts
c	nin = 0
c	nout = 0
c	ninter = 0
c      do k=1,igrid
c        do j=1,igrid
c          do i=1,igrid
c            do i1=1,3
c              if(iepsmp2(i,j,k,i1).gt.natom) then
c		    iepsmp(i,j,k,i1)=1 ! puts out problem spots
c		    nin = nin + 1
c	  	  else if(iepsmp2(i,j,k,i1).gt.0)then
c		    nin = nin + 1
c		    iepsmp(i,j,k,i1)=0 ! puts out vdw zone
c	  	  else if(iepsmp2(i,j,k,i1).lt.0)then
c		    ninter = ninter + 1
c		    iepsmp(i,j,k,i1)=1 ! puts out inter zone
c		  else
c		    iepsmp(i,j,k,i1)=0 
c		    nout = nout + 1
c		  end if
c            end do
c          end do
c        end do
c      end do
c	print *,'debug nin,ninter,nout,ndeb: ',nin,ninter,nout,ndeb ! debug
cc
cc dump map out for check
cc
c	filnam = 'mapgen.map'
c	call dmpeps(filnam)
c
	fin = cputime(start)
	print *,'total mapgen cpu time (s): ',fin
	return
c-----------------------------------------------
c	stop
900	print *,'error opening sas file'
	stop
901	print *,'unexpected end of sas file'
	stop
902	print *,'error reading sas file'
	stop
	end
c--------------------------
