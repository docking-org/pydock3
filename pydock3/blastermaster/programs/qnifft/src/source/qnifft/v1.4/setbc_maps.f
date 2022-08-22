	subroutine setbc(ibctyp,qplus,qmin,cqplus,cqmin,epsout,
     &deblen,natom,phifili)
c
c assigns either zero boundary conditions (ibctyp=1)
c quasi-coulombic based on debye dipole qplus at cqplus
c and qmin at cqmin (ibctyp=2)
c focussing (ibctyp=3)(continuance if scales are the same)
c full quasi-coulombic (ibctyp=4)
c or constant external field (ibctyp=5)
c option 2 will be appropriately modified by periodic
c adapted to new phimap format 19 june 1994,kas
c--------------------------------------------------------------
	include 'qdiffpar.h'
c--------------------------------------------------------------
	logical iper(3)
	dimension cqplus(3),cqmin(3),g(3),go(3),c(3)
	character*20 toblbl
	character*10 label
	character*60 title
	character*80 phifili
	character*16 botlbl
	real*4 realv(30)
	integer intv(30)
	character*80 titles(5)
c--------------------------------------------------------------
c prevent divide by zero distance:
	tiny = 1.e-6 
c
c zero option, clear boundary values
c
          do iz=1,igrid
             do iy=1,igrid
                do ix=1,igrid,igrid-1
                   phimap(ix,iy,iz) = 0.0
	    end do
	    end do
	    end do
          do iz=1,igrid
             do iy=1,igrid,igrid-1
                do ix=1,igrid
                   phimap(ix,iy,iz) = 0.0
	    end do
	    end do
	    end do
          do iz=1,igrid,igrid-1
             do iy=1,igrid
                do ix=1,igrid
                   phimap(ix,iy,iz) = 0.0
	    end do
	    end do
	    end do
c
c end of zero option
c
      if(ibctyp.eq.5) then
c
c constant field, 1 kT/e/grid unit, in the x direction
c
          do iz=1,igrid
             do iy=1,igrid
                do ix=1,igrid
                   phimap(ix,iy,iz) = ix
	    end do
	    end do
	    end do
	end if
c
c end external field option
c
      if(ibctyp.eq.2) then
c
c quasi coulombic dipole option
c
	do iz=1,igrid
	  do iy=1,igrid
	    do ix=1,igrid,igrid-1
            dist = (cqplus(1)-ix)**2 + (cqplus(2)-iy)**2
     &	     + (cqplus(3)-iz)**2
 	      dist = sqrt(dist)/scale + tiny
 	      tempp = qplus*exp(-dist/deblen )/(dist*epsout) 
            dist = (cqmin(1)-ix)**2 + (cqmin(2)-iy)**2
     &	     + (cqmin(3)-iz)**2
 	      dist = sqrt(dist)/scale + tiny
 	      tempn = qmin*exp(-dist/deblen )/(dist*epsout) 
            phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempp + tempn
	    end do
	    end do
	    end do
	do iz=1,igrid
	  do iy=1,igrid,igrid-1
	    do ix=1,igrid
            dist = (cqplus(1)-ix)**2 + (cqplus(2)-iy)**2
     &	     + (cqplus(3)-iz)**2
 	      dist = sqrt(dist)/scale + tiny
 	      tempp = qplus*exp(-dist/deblen )/(dist*epsout) 
            dist = (cqmin(1)-ix)**2 + (cqmin(2)-iy)**2
     &	     + (cqmin(3)-iz)**2
 	      dist = sqrt(dist)/scale + tiny
 	      tempn = qmin*exp(-dist/deblen )/(dist*epsout) 
            phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempp + tempn
	    end do
	    end do
	    end do
	do iz=1,igrid,igrid-1
	  do iy=1,igrid
	    do ix=1,igrid
            dist = (cqplus(1)-ix)**2 + (cqplus(2)-iy)**2
     &	     + (cqplus(3)-iz)**2
 	      dist = sqrt(dist)/scale + tiny
 	      tempp = qplus*exp(-dist/deblen )/(dist*epsout) 
            dist = (cqmin(1)-ix)**2 + (cqmin(2)-iy)**2
     &	     + (cqmin(3)-iz)**2
 	      dist = sqrt(dist)/scale + tiny
 	      tempn = qmin*exp(-dist/deblen )/(dist*epsout) 
            phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempp + tempn
	    end do
	    end do
	    end do
      end if 
c
c end of quasi coulombic dipole option
c
      if(ibctyp.eq.4) then
c
c a summation of the potential resulted from each point of charge 
c
	  do ic = 1,natom
	  if(atcrg(ic).ne.0.)then
c	    print *,'atom, q: ',ic,atcrg(ic)
	    do iz=1,igrid
	    do iy=1,igrid
	    do ix=1,igrid,igrid-1
              dist = (atcrd(1,ic)-ix)**2 + (atcrd(2,ic)-iy)**2
     &	     + (atcrd(3,ic)-iz)**2
 	        dist = sqrt(dist)/scale + tiny
	        tempd = atcrg(ic)*exp(-dist/deblen)/(dist*epsout)
              phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempd 
	    end do
	    end do
	    end do
	    do iz=1,igrid
	    do iy=1,igrid,igrid-1
	    do ix=1,igrid
              dist = (atcrd(1,ic)-ix)**2 + (atcrd(2,ic)-iy)**2
     &	     + (atcrd(3,ic)-iz)**2
 	        dist = sqrt(dist)/scale + tiny
	        tempd = atcrg(ic)*exp(-dist/deblen)/(dist*epsout)
              phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempd 
	    end do
	    end do
	    end do
	    do iz=1,igrid,igrid-1
	    do iy=1,igrid
	    do ix=1,igrid
              dist = (atcrd(1,ic)-ix)**2 + (atcrd(2,ic)-iy)**2
     &	     + (atcrd(3,ic)-iz)**2
 	        dist = sqrt(dist)/scale + tiny
	        tempd = atcrg(ic)*exp(-dist/deblen)/(dist*epsout)
              phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempd 
	    end do
	    end do
	    end do
	  end if
	  end do
      end if 
c
c end of the option for the complete summation of potential
c
c
      if(ibctyp.eq.3) then
c 
c focussing option-bc's come from a previous phimap
c
	nunit=18
	open(unit=nunit,status='old',file=phifili,err=900,form='unformatted')
	write(6,*)' '
	write(6,*)'focussing boundary condition read from file'
	write(6,*)phifili
	write(6,*)' '
	call getphi(nunit,ngrid,igrid1,phimap,realv,intv,titles)
	scale1 = realv(1)
	do k = 1,3
	  oldmid1(k) = realv(k+1)
	end do
	close(nunit)
c check to see if this is a continuence
	
	if(igrid1.ne.igrid)then
	  print *,'WARNING: current,focussing grid should be same dimension'
	  print *,igrid,igrid1
	  stop
	end if
	if(scale1.eq.scale) then
	write(6,*) 'scales are the same.' 
	write(6,*) 'therefore assuming this to be a continuence'
	goto 511
	endif
	write(6,*)'original scale (grids/A)      : ',scale1
	write(6,*)'object centre at (A) : ',oldmid1
	write(6,*)' '
c
c check to see that new grid lies within old one that is going to
c provide bc's
c
      iout = 0
	goff = (ngrid+1.)/2.
	do iz=1,igrid,igrid-1
        g(3) = iz
	  do iy=1,igrid,igrid-1
          g(2) = iy
	    do ix=1,igrid,igrid-1
            g(1) = ix
c
c for each new grid corner, calculate old grid coords
c
            call gtoc(g,c)
            do i = 1,3
		  gold = (c(i) - oldmid1(i))*scale1 + goff
              if((gold.le.2.).or.(gold.ge.ngrid-1.))  iout = 1
	      end do
	end do
	end do
	end do
      if(iout.ne.0) then
	  write(6,*)'part of new grid lies outside old grid'
	  write(6,*)'check scaling of both grids'
	  write(6,*)'old scale:'
	  write(6,*)'scale (grids/A)      : ',scale1
	  write(6,*)'object centre at (A) : ',oldmid1
	  write(6,*)'new scale:'
	  write(6,*)'scale (grids/A)      : ',scale
	  write(6,*)'object centre at (A) : ',oldmid
	  stop
      end if
c
c for each boundary point
c convert to real coordinates
c convert to old grid coordinates
c interpolate potential
c note that can use same potential array for boundaries
c since old potentials at boundary are not used for new ones
c
c
c save new grid size, and set temporarily to 65
c
	isgrid = igrid
	igrid = ngrid
	gmid = (isgrid + 1.)/2.
	write(6,*)igrid
      write(6,*)'pulling boundary values out of old potential map...'
	do 9022 iz=2,isgrid-1
        g(3) = iz
	  do 9023 iy=2,isgrid-1
          g(2) = iy
	    do 9024 ix=1,isgrid,isgrid-1
            g(1) = ix
c
c for each new grid side, calculate old grid coords
c
		do 9025 i = 1,3
		  c(i) = (g(i) - gmid)/scale + oldmid(i)
		  go(i) = (c(i) - oldmid1(i))*scale1 + goff

9025		continue
c
c find potential
c
            call phintp(go,phiv)
            phimap(ix,iy,iz) = phiv

9024		    continue
9023		 continue
9022	    continue

	do 9026 iz=2,isgrid-1
        g(3) = iz
	  do 9027 iy=1,isgrid,isgrid-1
          g(2) = iy
	    do 9028 ix=2,isgrid-1
            g(1) = ix
c
c for each new grid side, calculate old grid coords
c
		do 9029 i = 1,3
		  c(i) = (g(i) - gmid)/scale + oldmid(i)
		  go(i) = (c(i) - oldmid1(i))*scale1 + goff
9029		continue
c
c find potential
c
            call phintp(go,phiv)
            phimap(ix,iy,iz) = phiv

9028		    continue
9027		 continue
9026	    continue

	do 9030 iz=1,isgrid,isgrid-1
        g(3) = iz
	  do 9031 iy=2,isgrid-1
          g(2) = iy
	    do 9032 ix=2,isgrid-1
            g(1) = ix
c
c for each new grid side, calculate old grid coords
c
		do 9033 i = 1,3
		  c(i) = (g(i) - gmid)/scale + oldmid(i)
		  go(i) = (c(i) - oldmid1(i))*scale1 + goff
9033		continue
c
c find potential
c
            call phintp(go,phiv)
            phimap(ix,iy,iz) = phiv
c
9032		    continue
9031		 continue
9030	    continue
c restore new grid size
c
	igrid = isgrid
511   end if 
c
c end of focussing option
c
c---------------------------------------------------------------
	midg = (igrid+1)/2
c	write(6,*)' some initial phi values: '
c	write(6,*)' midg,midg,1; midg,midg,igrid '
c	write(6,*)phimap(midg,midg,1), phimap(midg,midg,igrid)
c	write(6,*)' midg,1,midg; midg,igrid,midg '
c	write(6,*)phimap(midg,1,midg),phimap(midg,igrid,midg)
c	write(6,*)' 1,midg,midg; igrid,midg,midg '
c	write(6,*)phimap(1,midg,midg),phimap(igrid,midg,midg)
C---------------------------------------------------------------
	return
900	write(6,*)' no potl map for focussing boundary conditions'
   	stop
	end
