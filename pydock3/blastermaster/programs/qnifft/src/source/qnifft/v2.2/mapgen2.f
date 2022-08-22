	subroutine mapgen(exrad,radprb,natom)
c--------------------------
c implement less scale dependent generation of eps, debye mapping
c--------------------------
c	implicit none
	integer ngirdmx
	parameter (ngirdmx=100000)
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
	real*4 distsq,distsq1,distsq2,distsq3,dmin2
	real roff(3,3) ! offsets to grid line mid points
	logical it1(0:6)
	integer inear(6)
c
	integer iepsmp2(ngrid,ngrid,ngrid,3) ! stores atom numbers  temporarily
	integer ihead(ngrid,ngrid,ngrid) ! entry into linked list of girdle points
	integer ilink(ngirdmx) ! linked list of girdle points
	real*4  gxyz(3,ngirdmx) ! linked list of girdle points
	real*4 ratio,ratmx,ratmn

c--------------------------
	real*4 start,fin,cputime
	real*4 xran,radmax
	real*4 pi,ang,dang,zed,dzed,crad
c--------------------------
	integer i,j,k,n,i1,j1,k1,n1
	integer i2,j2,k2,n2,inbr,lunit
	integer nin,ninter,ndeb,nout,ncheck,nassign
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
            iatmapd(i,j,k)  = 0 
            ihead(i,j,k)  = 0 
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
c
c debug- count point classes
c and set iatmap to flag any mid points in interzone that need checking later
c
	nin = 0
	nout = 0
	ninter = 0
	ndeb = 0
      do k=1,igrid
        do j=1,igrid
          do i=1,igrid
		iatmap(i,j,k) = 0
            do i1=1,3
              if(iepsmp2(i,j,k,i1).gt.0) then
		    nin = nin + 1
	  	  else if(iepsmp2(i,j,k,i1).lt.0)then
		    ninter = ninter + 1
		    iatmap(i,j,k) = 1
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
	  call pntchk(natom,atcrd,satradp,prad,scale,xyz,nbad,ngood)
c
c find upper, lower grid indices that enclose atom
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
	  ninter = 0
	  do k = llim(3),ulim(3)
	    do j = llim(2),ulim(2)
	      do i = llim(1),ulim(1)
		  if(iatmap(i,j,k).eq.1)then ! interzone midpoint(s) at 1/2/3 
		    ninter = ninter + 1
	          distsq = dxarr2(i,1) + dxarr2(j,2) + dxarr2(k,3)
		    distsq1 = distsq + dxarr(i,1)
		    distsq2 = distsq + dxarr(j,2)
		    distsq3 = distsq + dxarr(k,3)
		    if(distsq1.lt.sradmq)then
		      iepsmp2(i,j,k,1) = 0
		    end if
		    if(distsq2.lt.sradmq)then
		      iepsmp2(i,j,k,2) = 0
		    end if
		    if(distsq3.lt.sradmq)then
		      iepsmp2(i,j,k,3) = 0
		    end if
		  end if
	      end do
	    end do
	  end do
	end do ! end of point loop
	print *,'# of bad, good points: ',nbad,ngood
223   format(3f12.7,a4)
c
c map all girdle points to a linked list
c
	print *,'making linked list of girdle points...'
	rewind(lunit)
	read(lunit,'(a)')head
	if(ngird.gt.ngirdmx)then
	  print *,'WARNING: increase size of girdle point list from, to: ',ngirdmx,ngird
	  ngird = ngirdmx
	end if
	do n = 1,ngird
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
c
c debug check linked list
c
      print *,'linked list i,j,k,jat: ' ! debug
	ngood = 0
      do k = 1,igrid
        do j = 1,igrid
          do i = 1,igrid
            iat = ihead(i,j,k)
            do while(iat.ne.0)
		  ngood = ngood + 1
c              print *,i,j,k,iat ! debug
              iat = ilink(iat)
            end do
          end do
        end do
      end do
	print *,'points in girdle,linked list: ',ngird,ngood


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
c find boundary points, and assign to an atom if possible
c
	icutl = prad + 2
	print *,'icutl for girdle search: ',icutl
	ngood = 0
	ratmn = 1.e6
	ratmx = -1.e6
c
	ninter = 0
	nassign = 0
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
		  if(iat.ne.0)then
		    nassign = nassign + 1
		  else
c
c search linked list for closest girdle point
c
c		    print *,'finding girdle for orphan point ',i,j,k ! debug
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
	          dmin2 = 1.e6
c		    print *,'searching over ',llim,ulim
		    nlook = 0
                do k2 = llim(3),ulim(3)
                  do j2 = llim(2),ulim(2)
                    do i2 = llim(1),ulim(1)
                      jat = ihead(i2,j2,k2)
                      do while(jat.ne.0)
				nlook = nlook + 1
                        distsq = (i - gxyz(1,jat))**2 + (j - gxyz(2,jat))**2 + (k - gxyz(3,jat))**2
c			      print *,'looking at point: ',jat,sqrt(distsq) ! debug
                        if(distsq.lt.dmin2)then
		              dmin2 = distsq
		              iat = -jat
		            end if
              		jat = ilink(jat)
	                end do
	              end do
	            end do
	          end do
		    ngood = ngood + 1
		    do k1 = 1,3
		      xyz(k1) = gxyz(k1,abs(iat))
		    end do
		    ratio = sqrt(dmin2)/prad
		    ratmx = max(ratmx,ratio)
		    ratmn = min(ratmn,ratio)
c		    print *,'found ',i,j,k,' at ',iat,xyz,ratio
c		    print *,'after searching ',nlook
		  end if
		  iatmap(i,j,k) = iat
	      end if
	    end do
	  end do
	end do
	print *,'# of grid points on molecular surface= ',ninter
	print *,'# of molecular surface points assigned to an atom = ',nassign
	print *,'# of molecular surface points assigned to girdle = ',ngood
	print *,'worst reentrant surface distance error ratios: ',ratmn,ratmx
c
c dump map out for check
c
	filnam = 'mapgen.map'
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
