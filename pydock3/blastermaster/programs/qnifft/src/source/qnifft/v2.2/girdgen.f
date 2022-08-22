	subroutine girdgen(natom,atcrd,atrad,prad,scale,ngird,lunit)
c-----------------------------------------------------------
	implicit none
	integer natom,ncrc,nfct
	parameter (ncrc=200,nfct=190)
c-----------------------------------------------------------
c atom coords, radii
	real*4 atcrd(3,natom),atrad(natom),prad,scale
	real*4 gcen(3),xcen(3)
c circle coordinates
	real*4 crccrd(3,ncrc),crccrd1(3,ncrc)
c flags
	logical incrc(ncrc)
	integer istart(natom),istop(natom)
	integer ipair(2,natom*nfct),itrip(3,natom*nfct)
c-----------------------------------------------------------
	real*4 pi
	real*4 zed,dz,tarea,xmin(3),xmax(3)
	real*4 dang,ang,crad
	real*4 dist2,cdist2,trad2,radmax
	real*4 start,fin
	real*4 secnds
c lot of stuff for circle checking
	real*4 dxij,dyij,dzij,dij,dij2,rij,rij2,pij,qij
	real*4 rdsmij,rddfij,trfact,trxij,tryij,trzij
	real*4 dxy,tiny
	real*4 sintht,costht,sinphi,cosphi
	real*4 rtmt11,rtmt21,rtmt31,rtmt12,rtmt22,rtmt32,rtmt13,rtmt23,rtmt33
	real*4 big
c-----------------------------------------------------------
	real*4 xyz(3,3),rad(3),xyz3(3,2)
c-----------------------------------------------------------
	integer il,i,j,k,n,npt,i1,nhatm,need,npair,ntrip,ntripp
	integer ngird,lunit
	integer j1,k1,ip1,jp1,nleft,n1,itrp,i2,ip2,ndo
c-----------------------------------------------------------
	character point*4,dots*4
	data tiny / 1.e-8 / ! prevent /0
c-----------------------------------------------------------
	pi = 355./113.
	print *,'generating GIRDLE points...'
c	lunit = lunit + 1
c
c initialize 
c sasgen is passed scaled atom rad+probe rad
c
	radmax = 0.
	do i = 1,natom
	  istart(i) = 0
	  istop(i) = 0
	  radmax = max(radmax,atrad(i))
	end do
	print *,'max radius (grids)          : ',radmax
c
c
c SAS point generation
c write to scratch file
c
	dots='DOTS'
c	write(lunit,'(a)')dots
c
c generate canonical circle of points for SAS area
c
	dang = 2*pi/ncrc
	point = '  60' !magenta
c
c cap poles
c
	do j=1,ncrc
	    ang = (j-1)*dang
	    crccrd(1,j) = sin(ang)
	    crccrd(2,j) = cos(ang)
	    crccrd(3,j) = 0.
c	    write(6,*)(crccrd(k,j),k=1,3)
	end do
	print *,'# of circle points, coarsest spacing (grids) : ',ncrc,radmax*dang
	need = 2.*pi*radmax/0.25 ! would like largest spacing to be < 1/4 grid
	if(need.gt.ncrc)then
	  print *,'Recommend increasing ncrc above ',ncrc,' to ',need
	  print *,'if time and memory permit'
	end if
c
c loop over all atoms and find pairs
c
	start = secnds(0.0)
	npair = 0
	do i = 1,natom
	  if(mod(i,1+natom/10).eq.1)print *,'Atom: ',i
	  istart(i) = npair+1
	  if(atrad(i).le.0)goto 700
	  do j = 1,natom
	    if((j.ne.i).and.(atrad(j).gt.0.))then
	      cdist2 = (atrad(i)+atrad(j))**2
	      dist2 = 0.
	      do k = 1,3
	        dist2 = dist2 + (atcrd(k,j)-atcrd(k,i))**2
		  if(dist2.ge.cdist2)goto 203
	      end do
	      if(npair.ge.natom*nfct)then
	        print *,'exceeding max number of atom pairs- increase nfct > ',nfct
		  stop
	      end if
	      npair = npair + 1
	      ipair(1,npair) = i
	      ipair(2,npair) = j
	    end if
203	  continue! with next atom
	  end do
700	  continue ! skip to here if zero radius
	  istop(i) = npair
	end do
c	istart(natom) = npair
	print *,'number of pairs: ',npair
	print *,' '
cc
cc debug - check pair list
cc
c	do n = 1,npair
c	  print *,'npair,i,j: ',n,ipair(1,n),ipair(2,n)
c	end do
c	print *,'atom pair start, stop: '
c	do i = 1,natom
c	  print *,istart(i),istop(i)
c	end do
c	print *,' '
cc
cc debug - check pair list
cc
c	do i = 1,natom
c	  print *,'atom: ',i
c	  do i1 = istart(i),istop(i) ! loop over atom i's pairlist
c	    print *,'i,j: ',ipair(1,i1),ipair(2,i1)
c	  end do
c	end do
	fin = secnds(start)
	print *,'pair time: ',fin
	print *,' '
c
c generate girdle points
c loop over all atom pairs, i,j
c scale, rotate, translate canonical circle of points
c check all pairs of i: k not eq. j, and zap circle points if inside- write survivors
c
c
	ngird = 0
	ndo = 0
	do 100 n = 1,npair
	  i = ipair(1,n)
	  j = ipair(2,n)
	  if(i.ge.j)goto 100
c
c atom separation
c
	  rdsmij = (atrad(i) + atrad(j))**2
	  rddfij = (atrad(i) - atrad(j))**2
	  dxij = atcrd(1,j) - atcrd(1,i)
	  dyij = atcrd(2,j) - atcrd(2,i)
	  dzij = atcrd(3,j) - atcrd(3,i)
	  dij2 = dxij*dxij + dyij*dyij + dzij*dzij ! distance squared
c
c skip if too far apart to intersect, or one engulfs the other
c
	  pij = rdsmij - dij2
	  qij = dij2 - rddfij
	  if((pij.le.0.0).or.(qij.le.0.0))goto 100 
	  ndo = ndo + 1
	  rij = 0.5*sqrt(pij*qij/dij2) ! radius of circle of intersection (coi)
c	  print *,'rij: ',i,j,rij
c
c center of coi
c
	  trfact = 0.5*(atrad(i)**2-atrad(j)**2)/dij2 + 0.5
	  trxij = atcrd(1,i) + trfact*dxij
	  tryij = atcrd(2,i) + trfact*dyij
	  trzij = atcrd(3,i) + trfact*dzij
c	  print *,'translation: ',trxij,tryij,trzij
c
c rotation
c
	  dij = sqrt(dij2)
	  costht = dzij/dij
	  sintht = sqrt(1.0-costht**2)
	  dxy = sqrt(dxij**2 + dyij**2)
	  cosphi = dxij/(dxy + tiny)
	  sinphi = -dyij/(dxy + tiny)
c	  print *,'cos, sin theta: ',costht,sintht
c	  print *,'cos, sin phi: ',cosphi,sinphi
        rtmt11 = sinphi*(1.-costht)*sinphi + costht
        rtmt21 = sinphi*(1.-costht)*cosphi
        rtmt31 = sintht*cosphi
        rtmt12 = sinphi*(1.-costht)*cosphi
        rtmt22 = (1.-costht)*cosphi*cosphi+ costht
        rtmt32 = -sintht*sinphi
        rtmt13 = -sintht*cosphi
        rtmt23 = sintht*sinphi
        rtmt33 = costht
c	  print *,'rotation: '
c	  print *,'---------------------- '
c	  print *,rtmt11,rtmt21,rtmt31
c	  print *,rtmt12,rtmt22,rtmt32
c	  print *,rtmt13,rtmt23,rtmt33
c	  print *,' ----------------------'
c
c transform coi
c
	  do n1 = 1,ncrc
	    incrc(n1) = .true.
	    crccrd1(1,n1) =  trxij + rij*(rtmt11*crccrd(1,n1)+rtmt21*crccrd(2,n1) + rtmt31*crccrd(3,n1))
	    crccrd1(2,n1) =  tryij + rij*(rtmt12*crccrd(1,n1)+rtmt22*crccrd(2,n1) + rtmt32*crccrd(3,n1))
	    crccrd1(3,n1) =  trzij + rij*(rtmt13*crccrd(1,n1)+rtmt23*crccrd(2,n1) + rtmt33*crccrd(3,n1))
	  end do
	  nleft = ncrc
c
c loop over all i's pairs- except j
c
	  do 200 i1 = istart(i),istop(i)
	    ip1 = ipair(2,i1)
	    if(ip1.eq.j) goto 200
c	    print *,'check COI of ',i,j,' against ',ip1
	    trad2 = atrad(ip1)**2
	    do 300 n1 = 1,ncrc
	      if(incrc(n1))then
		  dist2 = 0.
		  do k = 1,3
		    dist2 = dist2 + (crccrd1(k,n1)-atcrd(k,ip1))**2
		    if(dist2.gt.trad2) goto 300 ! not in ipts SAS, survives, go to next point
	        end do
c debug
		  incrc(n1) = .false. ! inside ips SAS, zap point
		  nleft = nleft - 1
c debug
		  if(nleft.eq.0) goto 100 ! nothing left on this COI- skip to next
	      end if
300	    continue ! next sphere point
200	  continue ! next atom
c
c if we reach here, checked all circle points against all atoms that overlap i
c and some points survive-
	  big = 4.
	  do n1 = 1,ncrc
	    if(incrc(n1))then
	     write(lunit,223)(crccrd1(k,n1),k=1,3),point
	     do k = 1,3
	       crccrd1(k,n1) = big*crccrd1(k,n1)
	     end do
	     ngird = ngird + 1
c	     write(lunit,223)'HETATM    1  h    SPH    1    ',(crccrd1(k,n1),k=1,3),'  2.00  1.000'
	    end if
	  end do
c 
100	continue ! next pair COI
	print *,'# of accessible girdle points found: ',ngird
	fin = secnds(start)
	print *,'pairs checked: ',ndo
	print *,'girdle time: ',fin
c
c loop over all pairs and find triples:
c for each atom i, search its pair list for all unique i1,j1 atoms it pairs with
c Then search i1's pair list to see if j1 is found their
c note pair list has indices in order ipair(1,i) < ipair(2,i)
c note i<i1<j1
c
c	lunit = lunit + 1
	print *,'generating triples....'
	ntrip = 0
	ntripp = 0
	do i = 1,natom
	  if(mod(i,1+natom/10).eq.1)print *,'atom: ',i
	  do i1 = istart(i),istop(i) ! loop over atom i's pairlist
	    ip1 = ipair(2,i1)
	    if(ip1.gt.i)then
	      do j1 = i1+1,istop(i)      ! loop over atom i's pairlist for each pair there
		  jp1 = ipair(2,j1)
c	        print *,'my pairs: ',ip1,jp1
	        do k1 = istart(ip1),istop(ip1)
		    if(ipair(2,k1).eq.jp1)then
c	            if(ntrip.ge.natom*nfct)then
c	              print *,'exceeding max number of atom triples- increase nfct > ',nfct
c		        stop
c	            end if
		      ntrip = ntrip + 1
c		      itrip(1,ntrip) = i
c		      itrip(2,ntrip) = ip1
c		      itrip(3,ntrip) = jp1
			do k = 1,3
			  xyz(k,1) = atcrd(k,i)
			  xyz(k,2) = atcrd(k,ip1)
			  xyz(k,3) = atcrd(k,jp1)
			end do
			rad(1) = atrad(i)
			rad(2) = atrad(ip1)
			rad(3) = atrad(jp1)
c			print *,'triple: ',i,ip1,jp1
			call trppnt(xyz,rad,xyz3,itrp)
c			print *,'out of trppnt, itrp: ',itrp
c			print *,'out of trppnt, xyz3: ',xyz3
			if(itrp.eq.2)then
c we've found two triple points- check against i's pair list for collision
c except ip1, jp1
c 
	              do 302 n1 = 1,2
			    incrc(n1) = .true.
	                do 202 i2 = istart(i),istop(i)
	                  ip2 = ipair(2,i2)
	                  if((ip2.eq.ip1).or.(ip2.eq.jp1)) goto 202
c	                  print *,'check TRP of ',i,ip1,jp1,' against ',ip2
	                  trad2 = atrad(ip2)**2
		            dist2 = 0.
		            do k = 1,3
		              dist2 = dist2 + (xyz3(k,n1)-atcrd(k,ip2))**2
		              if(dist2.gt.trad2) goto 202 ! not in ip2s SAS, survives, go to next atom
	                  end do
c if we are here, inside an atom- zap point and skip to next point
				incrc(n1) = .false.
				goto 302
202	                continue ! next atom
302			  continue ! next point
	              do n1 = 1,2
	                if(incrc(n1))then
c	                   write(lunit,223)(crccrd1(k,n1),k=1,3),point !a bug?!
	                   write(lunit,223)(xyz3(k,n1),k=1,3),point
	                   do k = 1,3
	                     xyz3(k,n1) = big*xyz3(k,n1)
	                   end do
	                   ntripp = ntripp + 1
c	                   write(lunit,223)'HETATM    1  P    SPH    2    ',(xyz3(k,n1),k=1,3),'  2.00  1.000'
	                end if
	              end do
			end if
		    end if
	        end do
	      end do
	    end if
	  end do
	end do
cc
cc debug - check triple list
cc
c	do n = 1,ntrip
c	  print *,',ntrip,i,j,k: ',n,itrip(1,n),itrip(2,n),itrip(3,n)
c	end do
cc
	print *,'number of triples: ',ntrip
	print *,'# of accessible triple points found: ',ntripp
	ngird = ngird + ntripp
	print *,'total number of accessible intersection points: ',ngird
	print *,' '
	fin = secnds(start)
	print *,'triple time: ',fin
c

223	format(3f12.7,a4)
c223	format(a30,3f8.3,a13) !debug
c	lunit = lunit - 2
c	stop
	end

	subroutine trppnt(xyz,rad,xyz3,itrp)
c
c generate triple points
c
	implicit none
	real*4 xyz(3,3),rad(3),xyz3(3,2)
	integer itrp
	integer iring(2,3),i,j,k
	data iring / 1,2, 2,3, 3,1 /
c-------------------------------------------------
	real*4 dxij(3),dyij(3),dzij(3),dij(3),dij2(3),rij(3),rij2(3),pij(3),qij(3)
	real*4 rdsmij(3),rddfij(3),trfact(3),trxij(3),tryij(3),trzij(3)
c-------------------------------------------------
	real*4 a0(3),a1(3),a2(3),amag
	real*4 b0(3),b1(3),b2(3),bmag
	real*4 c0(3),c1(3),c2(3),cmag
	real*4 sinta,sintb,dena,denb
	real*4 costa,costb
c-------------------------------------------------
	itrp = 2
c	print *,'in trppnt: '
c	print *,'xyz: ',xyz
c	print *,'rad: ',rad
	dij2(2) = (xyz(1,2)-xyz(1,3))**2 + (xyz(1,2)-xyz(1,3))**2 + (xyz(1,2)-xyz(1,3))**2 ! distance squared
	do k = 1,3,2
	  i = iring(1,k)
	  j = iring(2,k)
c
c atom separation
c
	  rdsmij(k) = (rad(i) + rad(j))**2
	  rddfij(k) = (rad(i) - rad(j))**2
	  dxij(k) = xyz(1,j) - xyz(1,i)
	  dyij(k) = xyz(2,j) - xyz(2,i)
	  dzij(k) = xyz(3,j) - xyz(3,i)
	  dij2(k) = dxij(k)**2 + dyij(k)**2 + dzij(k)**2 ! distance squared
c
c skip if too far apart to intersect, or one engulfs the other
c
	  pij(k) = rdsmij(k) - dij2(k)
	  qij(k) = dij2(k) - rddfij(k)
	  if((pij(k).le.0.0).or.(qij(k).le.0.0))then
	    itrp = 0
	    return
	  end if
	  rij(k) = 0.5*sqrt(pij(k)*qij(k)/dij2(k)) ! radius of circle of intersection (coi)
c
c center of coi
c
	  trfact(k) = 0.5*(rad(i)**2-rad(j)**2)/dij2(k) + 0.5
	  trxij(k) = xyz(1,i) + trfact(k)*dxij(k)
	  tryij(k) = xyz(2,i) + trfact(k)*dyij(k)
	  trzij(k) = xyz(3,i) + trfact(k)*dzij(k)
      end do
c	print *,'in trppnt dij2: ',dij2
c	print *,'in trppnt rij: ',rij
c	print *,'x translation: ',trxij
c	print *,'y translation: ',tryij
c	print *,'z translation: ',trzij
	if(dij2(2).gt.(dij2(1)+dij2(3)))then
c	  print *,'Whoa, we have an angle > 90...'
	  itrp = 0
	  return
	end if
c unit vectors from atom i to j,k
	a1(1) = dxij(1)
	b1(1) = -dxij(3)
	a1(2) = dyij(1)
	b1(2) = -dyij(3)
	a1(3) = dzij(1)
	b1(3) = -dzij(3)
	call norm(a1,amag)
	call norm(b1,bmag)
c normal to ijk atom center plane
	call cross(a1,b1,c1)
	call norm(c1,cmag)
c normal vectors to projection of intersection point in plane
	call cross(c1,a1,a2)
	call cross(b1,c1,b2)
c centers of i's COI's with j,k
	a0(1) = trxij(1) 
	b0(1) = trxij(3) 
	a0(2) = tryij(1) 
	b0(2) = tryij(3) 
	a0(3) = trzij(1) 
	b0(3) = trzij(3) 
c	print *,' a0,a1,a2: ',a0,a1,a2
c	print *,' b0,b1,b2: ',b0,b1,b2
c	print *,'c1: ',c1
c point of intersection
	dena = 0.
	denb = 0.
	do k = 1,3
	  dena = dena + rij(1)*b1(k)*a2(k)
	  denb = denb + rij(3)*b2(k)*a1(k)
	end do
	if((abs(dena).lt.1.e-6).or.(abs(denb).lt.1.e-6))then
c	  print *, 'atoms ijk collinear: ',dena,denb
c	  print *,xyz
c	  print *,rad
	  itrp = 0
	  return
	end if
c	print *, 'dena,denb: ',dena,denb
	costa = 0.
	costb = 0.
	do k = 1,3
	  costa = costa + b1(k)*(b0(k)-a0(k))/dena
	  costb = costb + a1(k)*(a0(k)-b0(k))/denb
	end do
c	print *,'costa: ',costa
c	print *,'costb: ',costb
	if((costa.le.0.).or.(costb.le.0.))then
c	  print *,'triple hole!'
	  itrp = 0
	  return
	end if
	if((costa.ge.1.).or.(costb.ge.1.))then
c	  print *,'triple hole!'
	  itrp = 0
	  return
	end if
	sinta = sqrt(1.-costa**2)
	sintb = sqrt(1.-costb**2)
c	print *,'sinta: ',sinta
c	print *,'sintb: ',sintb
	do k = 1,3
	 xyz3(k,1) = rij(1)*(sinta*c1(k) + costa*a2(k)) + a0(k)
	 xyz3(k,2) = rij(3)*(-sintb*c1(k) + costb*b2(k)) + b0(k)
	end do
c	print *,'xyz3 ',xyz3(1,1),xyz3(2,1),xyz3(3,1)
c	print *,'     ',xyz3(1,2),xyz3(2,2),xyz3(3,2)
c
c debug check
c
	do i = 1,3
	  dij(i) = 0.
	  do k = 1,3
	    dij(i) = dij(i) + (xyz3(k,1)- xyz(k,i))**2
	  end do
	  dij(i) = sqrt(dij(i))-rad(i)
	  if(abs(dij(i)).gt.1.e-4) print *,'*******triple point 1 errors: ',dij(i)
      end do
	do i = 1,3
	  dij(i) = 0.
	  do k = 1,3
	    dij(i) = dij(i) + (xyz3(k,2)- xyz(k,i))**2
	  end do
	  dij(i) = sqrt(dij(i))-rad(i)
	  if(abs(dij(i)).gt.1.e-4) print *,'*******triple point 1 errors: ',dij(i)
      end do

	return
	end
