c---------------------------------------------------
c parameters for qdiff.f and subroutines
c---------------------------------------------------
	parameter (nclist = 3000)
	parameter (nrlist = 1000)
	parameter (nrmax = 1000)
	parameter (ncmax = 3000)
	parameter (natmx = 132000)
	parameter (ngrid = 65)
	parameter (nbndmx = ngrid*ngrid*20)
c---------------------------------------------------
c---------------------------------------------------
c
c	atom names for radii table
c
      dimension atnam(nrmax)		
c
c	residue names for radii table
c
      dimension rnam(nrmax)		
c
c	radius table 
c
      dimension radt(nrmax)		
c
c	links for radius hash table
c
      dimension irlink(nrlist)	
c
c	radii id numbers in hash table
c
      dimension irnumb(nrlist)	
c
c	atom names   for charge table
c
      dimension catnam(ncmax)		
c
c	chain names for charge table
c
      dimension cchn(ncmax)		
c
c	residue names for charge table
c
      dimension crnam(ncmax)		
c
c	residue number for charge table
c
      dimension crnum(ncmax)		
c
c	charge table 
c
      dimension chrgvt(ncmax)		
c
c	links for charge hash table
c
      dimension iclink(nclist)	
c
c	charge entry id numbers in hash table
c
      dimension icnumb(nclist)	
c
c	atom charges for boundary condition
c
	dimension atcrd(3,natmx),atrad(natmx),atcrg(natmx)	
 	real*4 phimap(ngrid,ngrid,ngrid)
 	integer*1 iepsmp(ngrid,ngrid,ngrid,3)
 	integer iatmap(ngrid,ngrid,ngrid)
 	integer iatmapd(ngrid,ngrid,ngrid)
 	integer*1 debmap(ngrid,ngrid,ngrid)
 	real*4 qmap(ngrid,ngrid,ngrid)
	dimension oldmid(3),oldmid1(3)
	dimension qass(4,natmx),qsrf(4,natmx)
	dimension cion(4),zion(4),czion(4),czion2(4) ! kas 8-feb-01
	real*4 epstab(0:1),debtab(0:3),sepstab(0:6) ! kas 8-feb-01
	real*4 xyzbnd(3,nbndmx) ! corrected coords of grid boundary points ! kas 1-jun-01
c---------------------------------------------------
c---------------------------------------------------
      character*1 chn,schn,cchn
      character*3 crnam,rnam,sres,res
      character*4 rnum,snum,crnum
      character*6 atnam,catnam,atm
c---------------------------------------------------
c---------------------------------------------------
	common
     &	/lnkt/  irlink,irnumb,iclink,icnumb,irtot,ictot
     &	/name/  atnam,rnam,catnam,cchn,crnam,crnum
     &	/value/ radt,chrgvt
     &	/maps/  phimap,qmap,iatmap,iatmapd,atcrg,atcrd,atrad
     &	/imaps/  debmap,iepsmp
     &	/scale/ scale,oldmid,scale1,oldmid1,igrid
     &  /charges/ qass,nqass,qsrf,nqsrf	
     & / sion / cion,zion,czion,czion2 ! kas 8-feb-01
     & / tabs / epstab,debtab,sepstab ! kas 8-feb-01
     & / cbnd / xyzbnd,nbnd ! kas 1-jun-01
c---------------------------------------------------

