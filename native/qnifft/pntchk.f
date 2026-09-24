	subroutine pntchk(natom,atcrd,atrad,prad,scale,xyz,nbad,ngood)
c-----------------------------------------------------------
c check putative SAS point xyz to ensure it is outside all atom (rad+probe)
c slow as hell, but sure
c we are passed scaled atom rad+probe rad
c
	implicit none
	integer natom
c-----------------------------------------------------------
c atom coords, radii
	real*4 atcrd(3,natom),atrad(natom),prad,scale,xyz(3)
c-----------------------------------------------------------
	real*4 dist2,cdist2,trad2,slop
c-----------------------------------------------------------
	integer i,j,k,n,nbad,nbada,ngood,ngooda
c-----------------------------------------------------------
	data slop / 0.001 /
c-----------------------------------------------------------
c
c loop over all atoms
c
	nbada = 0
	ngooda = 0
c	print *,'in pntchk xyz, natom: ',xyz,natom
	do i = 1,natom
c
c check if within atoms sas
c
	  trad2 = atrad(i)**2
	  dist2 = 0.
	  do k = 1,3
	    dist2 = dist2 + (xyz(k)-atcrd(k,i))**2
	  end do
	  if((trad2-dist2).gt.slop)then
	    print *,'WARNING point below surface of atom: ',i,dist2,trad2
	    nbada = nbada + 1
	  else if(abs(trad2-dist2).lt.slop)then
c	    print *,'point on surface of atom: ',i,dist2,trad2
	    ngooda = ngooda + 1
	  end if
	end do
	if(nbada.gt.0)nbad = nbad + 1
	if(ngooda.gt.0)ngood = ngood + 1

	return
	end
