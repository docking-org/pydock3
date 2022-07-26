	subroutine sasgen(natom,atcrd,atrad,prad,scale,npoint,lunit)
c-----------------------------------------------------------
	implicit none
	integer natom,nsph
	parameter (nsph=1200)
c-----------------------------------------------------------
c atom coords, radii
	real*4 atcrd(3,natom),atrad(natom),prad,scale
	real*4 atarea(natom),gcen(3),xcen(3)
c sphere coordinates
	real*4 sphcrd(3,nsph),sphcrd1(3,nsph)
c flags
	integer insph(nsph),nonh(natom)
c-----------------------------------------------------------
	real*4 pi
	real*4 zed,dz,tarea,xmin(3),xmax(3)
	real*4 dang,ang,crad
	real*4 dist2,cdist2,trad2,radmax
	real*4 start,fin
	real*4 secnds
c-----------------------------------------------------------
	integer il,i,j,k,n,npt,i1,nhatm,nslab,nspt,need
	integer npoint,lunit
	integer ngird
c-----------------------------------------------------------
	character point*4,dots*4
c-----------------------------------------------------------
c-----------------------------------------------------------
c-----------------------------------------------------------
	pi = 355./113.
c
c write to scratch file
	dots='DOTS'
	write(lunit,'(a)')dots
c
c find points at intersection of 3 and 2 atoms
c
	ngird = 0
	call girdgen(natom,atcrd,atrad,prad,scale,ngird,lunit)
	print *,'generating SAS points...'
c
c initialize 
c sasgen is passed scaled atom rad+probe rad
c
	tarea = 0.
	do k = 1,3
	  xmin(k) = 1.e6
	  xmax(k) = -1.e6
	  xcen(k) = 0.
	end do
	radmax = 0.
	do i = 1,natom
	  atarea(i) = 0.
	  nonh(i) = 1
	  radmax = max(radmax,atrad(i))
	  do k = 1,3
	    xmin(k) = min(xmin(k),atcrd(k,i))
	    xmax(k) = max(xmax(k),atcrd(k,i))
	    xcen(k) = xcen(k) + atcrd(k,i)
	  enddo
	end do
c
c center
c
	do k = 1,3
	  gcen(k) = (xmin(k)+xmax(k))/2.
	  xcen(k) = xcen(k)/natom
	end do
	print *,' '
	print *,'Geometric center of molecule: ',gcen
	print *,'average center of molecule  : ',xcen
	print *,'max radius (grids)          : ',radmax
c
c SAS point generation
c
c
c
c generate canonical sphere of points for SAS area
c
	nslab = sqrt(1.*(nsph-2))
	dang = 2*pi/nslab
	dz = 2./nslab
c
	zed = -1.-dz/2.
	point = '  60' !magenta
c
c cap poles
c
	sphcrd(1,1) = 0.
	sphcrd(2,1) = 0.
	sphcrd(3,1) = 1.
	sphcrd(1,2) = 0.
	sphcrd(2,2) = 0.
	sphcrd(3,2) = -1.
	nspt = 2
	do i = 1,nslab
	  zed = zed+dz
	  crad = sqrt(1.-zed**2)
	  do j=1,nslab
	    nspt = nspt+1
	    ang = (j-1)*dang
	    sphcrd(1,nspt) = crad*sin(ang)
	    sphcrd(2,nspt) = crad*cos(ang)
	    sphcrd(3,nspt) = zed
c	    write(lunit,223)(sphcrd(k,n),k=1,3),point
	  end do
	end do
	print *,'# of sphere points, coarsest spacing (grids) : ',nspt,radmax*dz
	need = 2.*radmax/0.25 ! would like largest spacing to be < 1/4 grid
	if(need**2.gt.nspt)then
	  print *,'Recommend increasing nsph above ',nsph,' to ',need**2
	  print *,'if time and memory permit'
	end if
c
c loop over all atoms
c
	start = secnds(0.0)
	npoint = 0
	do i = 1,natom
	  if(mod(i,1+natom/10).eq.1)print *,'Atom: ',i
c	  if(nonh(i).eq.0)goto 700
c
c generate atom i's sphere of points
c
	  do n = 1,nspt
	    insph(n) = 1
	    do k = 1,3
	      sphcrd1(k,n) = sphcrd(k,n)*atrad(i)+atcrd(k,i)
	    end do
	  end do
	  npt = nspt
c 
c AREA/SAS energy
c
c
c loop over all atom pairs
c
	  do j = 1,natom
	    if((j.eq.i).or.(nonh(j).eq.0))goto 203
	    cdist2 = (atrad(i)+atrad(j))**2
	    dist2 = 0.
	    do k = 1,3
	      dist2 = dist2 + (atcrd(k,j)-atcrd(k,i))**2
		if(dist2.gt.cdist2)goto 203
	    end do
c
c loop over remaining sphere points
c
	    trad2 = atrad(j)**2
	    do n = 1,nspt
		if(insph(n).eq.1)then
c
c check if within second atoms sas
c
		  dist2 = 0.
	        do k = 1,3
	          dist2 = dist2 + (sphcrd1(k,n)-atcrd(k,j))**2
		    if(dist2.gt.trad2)goto 204
	        end do
c decrease # of surviving points if a new hit
		  insph(n) = 0
		  npt = npt - 1
		  if(npt.eq.0)goto 205! bail out of pair check
		end if
204		continue ! with next sphere point
	    end do
203	  continue! with next atom
	  end do
205	  continue
	  atarea(i) = npt*4.*pi*atrad(i)**2/nspt
	  tarea = tarea + atarea(i)
c
c write surviving points
c
	  point = ' 60' !magenta
	  do n = 1,nsph
	    if(insph(n).eq.1)then
		npoint = npoint + 1
	      write(lunit,223)(sphcrd1(k,n),k=1,3),point
	    end if
	  end do
700	continue
223	format(3f12.7,a4)
	end do
	fin = secnds(start)
	print *,'number of SAS points: ',npoint
	npoint = npoint + ngird
	print *,'number of SAS+GIRDLE points: ',npoint
	print *,'Area time: ',fin
	tarea = tarea/scale**2
	print *,'total SAS area: ', tarea 

	return
	end
