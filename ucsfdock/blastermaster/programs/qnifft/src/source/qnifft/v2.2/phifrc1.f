	subroutine phifrc(ianal,donon,rionst,epsin,epsout,natom,iconc,epkt,rtemp,atfrct,ergt)
c 
c compute rxn field forces
c
c-------------------------------------------------
	include 'qdiffpar.h'
c-------------------------------------------------
c
	integer indxat(natom)
	real*4 atfrcc(3,0:natom) ! coulombic force
	real*4 atfrcr(3,0:natom) ! reaction field force at atoms
	real*4 atfrcb(3,0:natom) ! atoms field force at surface element's atom
	real*4 atfrcd(3,0:natom) ! stress force at surface
	real*4 atfrci(3,0:natom) ! ionic force
	real*4 atfrct(3,natom) ! total solvent force
	real*4 atsrfq(0:natom)   ! induced surface charge mapped to atoms
	real*4 dr(3),srfx(3),de(3),demag
	real*4 df(3),ftrans(3),frot(3),xcen(3),dr1(3)
	real*4 tin(3,3),froti(3)
	real*4 efact,earea,dpres,phi0
	integer npeps(3),ide(3)
	character*8 hour
	logical donon,iconc,ianal
	real*4 calk,boltz
	real*4 dfrcgd(3)
      data cfact  / 6.0e-4  /   !converts moles to ions/cubic angstrom
					 ! ==Na,l/moles,Ang*3
	data boltz / 0.00199 /
c-------------------------------------------------------
	if(.not.ianal)return
	start = cputime(0.0)
c	call time(hour)
	print *,' '
c	write(6,*)'computing forces at: ',hour
	write(6,*)'computing forces...'
c-------------------------------------------------------
	calk = rtemp*boltz
	print *,'kT in kcal/mole: ',calk
	midg = (igrid+1)/2
c
	pi = 355./113.
	if(rionst.le.0.)donon=.false.
c
c collect assigned atomic charge
c
	nqass = 0
	qfix = 0.
	rmidg = (igrid+1)/2.
	do k = 1,3
	  xcen(k) = 0.
	end do
	do i = 1,natom
	  atrad(i) = atrad(i)*scale
	  do k = 1,3
	    xcen(k) = xcen(k) + atcrd(k,i)
	  end do
	  if(atcrg(i).ne.0.)then
	    nqass = nqass + 1
	    do k = 1,3
		qass(k,nqass) = atcrd(k,i)
	    end do
	    indxat(nqass) = i
	    qfix = qfix + atcrg(i)
	    qass(4,nqass) = atcrg(i)
	  end if
	end do
cc debug
c	do i = 1,nqass
c	  print *,'q,i: ',qass(4,i),indxat(i)
c	end do
cc debug
	print *,'total atomic charge: ',qfix
	do k = 1,3
	  xcen(k) = xcen(k)/natom
	end do
	print *,'geometric center: ',xcen
c
c init. forces, charge
c
	do i = 0,natom
	  atsrfq(i) = 0.
	  do k = 1,3
	    atfrcc(k,i) = 0.
	    atfrcd(k,i) = 0.
	    atfrcr(k,i) = 0.
	    atfrcb(k,i) = 0.
	    atfrci(k,i) = 0.
	  end do
	end do
c
c coulombic forces
c
	ergc = 0.
	do i = 1,nqass
	  do j = i+1,nqass
	    do k = 1,3
	      dr(k) = qass(k,i)-qass(k,j)
	    end do
	    r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
	    if(r2.gt.1.e-6)then
	      ergc = ergc + qass(4,i)*qass(4,j)/sqrt(r2)
	      rr3 = 1./(sqrt(r2)*r2)
	      do k = 1,3
	        dfrc = dr(k)*qass(4,i)*qass(4,j)*rr3
	        atfrcc(k,indxat(i)) = atfrcc(k,indxat(i)) + dfrc
	        atfrcc(k,indxat(j)) = atfrcc(k,indxat(j)) - dfrc
	      end do
	    end if
	  end do
	end do
	open(41,file='qnifft.frc')
	write(41,*)'start Coulomb forces (kT/A): '
	do i = 0,natom
	  do k = 1,3
	    atfrcc(k,i) = atfrcc(k,i)*scale**2/epsin
	  end do
	  write(41,'(i6,3f12.6)')i,(atfrcc(k,i),k=1,3)
	end do
c	close(41)
	ergc = ergc/epsin*scale
	write(41,*)'end coulombic energy in uniform Epsin (kT): ',ergc
	print *,'coulombic energy in uniform Epsin (kT): ',ergc
	fin = cputime(start)
	print *,'coulombic forces took ',fin,' sec'
c
c dielectric boundary forces
c
c
c surface charge reaction field energy
c
	ergs = 0.
	qsurf = 0.
	nmove = 0
	nsurf = 0
	nsurfq = 0
	nbound = 0
	efact = (epsin-epsout)/(8*pi)
	do k = 1,3
	  npeps(k) = 0
	end do
	earea = 0.
	dpres = 0.
	print *,'dielectric surface force factor: ',efact
	if(efact.eq.0)goto 301
c	print *,'NOT using iatmap to correct surface points...' ! if debug srat==1 set
	print *,' using iatmap to correct surface points...'
c	print *,' using iatmap to correct surface points for srat 0.8-1.2...'
c	print *,' using gridmapped atom charges for reaction...'
	do i = 2,igrid-1
	  do j = 2,igrid-1
	    do k = 2,igrid-1
		isumeps = iepsmp(i,j,k,1) + iepsmp(i-1,j,k,1)
     &		  + iepsmp(i,j,k,2) + iepsmp(i,j-1,k,2)
     &		  + iepsmp(i,j,k,3) + iepsmp(i,j,k-1,3)
		if((isumeps.ne.0).and.(isumeps.ne.6))nbound = nbound+1
		ide(1) = (iepsmp(i,j,k,1) - iepsmp(i-1,j,k,1))
		ide(2) = (iepsmp(i,j,k,2) - iepsmp(i,j-1,k,2))
		ide(3) = (iepsmp(i,j,k,3) - iepsmp(i,j,k-1,3))
		ibound = abs(ide(1)) + abs(ide(2)) + abs(ide(3))
		demag = sqrt(float(ibound))
	      if(demag.gt.0) then
		  if(qmap(i,j,k).ne.0.)then
		    nsurfq = nsurfq+1
		  else
		    nsurf = nsurf+1
		  end if
		  npeps(ibound) = npeps(ibound) + 1
		  earea = earea + demag
		  srfq = 6.*phimap(i,j,k)-(phimap(i+1,j,k)+phimap(i-1,j,k)
     &	  + phimap(i,j+1,k) + phimap(i,j-1,k)
     &	  + phimap(i,j,k+1) + phimap(i,j,k-1) + qmap(i,j,k)/epsin)
		   qsurf = qsurf + srfq
c
c correction position of surface charge using position of nearest atom
c iatmap contains indices of atoms forming surface
c
		  srfx(1) = i
		  srfx(2) = j
		  srfx(3) = k
		  iat = iatmap(i,j,k)
		  atsrfq(iat) = atsrfq(iat) + srfq
		  if(iat.ne.0)then !found a neighbor atom
		    dmin2 = (atcrd(1,iat)-srfx(1))**2 + (atcrd(2,iat)-srfx(2))**2 + (atcrd(3,iat)-srfx(3))**2
		    if(dmin2.gt.1.e-6)then
			nmove = nmove + 1
		      srat = atrad(iat)/sqrt(dmin2)
cc
c			srat = 1. ! debug
c			if((srat.lt.0.8).or.(srat.gt.1.2))srat = 1. ! debug
cc
			srfx(1) = atcrd(1,iat) + srat*(srfx(1)-atcrd(1,iat))
			srfx(2) = atcrd(2,iat) + srat*(srfx(2)-atcrd(2,iat))
			srfx(3) = atcrd(3,iat) + srat*(srfx(3)-atcrd(3,iat))
		    end if
		  end if
		  fmagnm = 0.
cc magnitude of normal component= Ein(norm).Eout(norm)
		  fmagnm = fmagnm+abs(ide(1))*(phimap(i,j,k)-phimap(i-1,j,k))*(phimap(i+1,j,k)-phimap(i,j,k))
		  fmagnm = fmagnm+abs(ide(2))*(phimap(i,j,k)-phimap(i,j-1,k))*(phimap(i,j+1,k)-phimap(i,j,k))
		  fmagnm = fmagnm+abs(ide(3))*(phimap(i,j,k)-phimap(i,j,k-1))*(phimap(i,j,k+1)-phimap(i,j,k))
c dielectric boundary forces
		  do k1 = 1,3
		    df(k1) = ide(k1)*efact*fmagnm
	          atfrcd(k1,iat) = atfrcd(k1,iat) - df(k1)
		  end do
		  dpres = dpres + efact*fmagnm*demag
c
c loop thru charges and calc reaction forces
c
		  do i1 = 1,nqass
c
c using actual atom charge positions
c
	          do k1 = 1,3
	            dr(k1) = qass(k1,i1)-srfx(k1)
	          end do
	          r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
	          if(r2.gt.1.e-6)then
		      ergs = ergs + qass(4,i1)*srfq/sqrt(r2)
	            rr3 = 1./(sqrt(r2)*r2)
	            do k1 = 1,3
			  dfrc = dr(k1)*qass(4,i1)*srfq*rr3
	              atfrcr(k1,indxat(i1)) = atfrcr(k1,indxat(i1)) + dfrc
	            end do
		    end if
c using actual atom charge positions
cc
cc using grid mapped positions
cc
c		    call gridphi(scale,igrid,srfx,srfq,qass(1,i1),ergsgd,dfrcgd)
c		    ergs = ergs + ergsgd
c	          do k1 = 1,3
c	            atfrcr(k1,indxat(i1)) = atfrcr(k1,indxat(i1)) + dfrcgd(k1)
c		    end do
		  end do
	      end if
          end do   
        end do      
	end do
	earea = earea/scale**2
	dpres = dpres/earea
301	continue ! to here if dielectrics the same
	qsurf = qsurf/(4*pi*scale*epkt)
	qsurfe = (epsin-epsout)*qfix/epsin/epsout/epkt
	print *,'surface points readjusted: ',nmove
	print *,'surface points with, without fixed charge: ',nsurfq,nsurf
	print *,'expected, actual induced surface charge: ',qsurfe,qsurf
	print *,'number of boundary points: ',nbound
	print *,'# eps boundary points,eps area (A**3): ',npeps,earea
	print *,'Mean dielectric pressure(kT/A**3): ',dpres
	ergs = ergs/8./pi
	print *,'dielectric rxn field energy from surface Q (kT): ',ergs
c
c ionic force
c
	earea = 0.
	dpres = 0.
	pout = 0.
	cout = 0.
	ergi = 0.
	do k = 1,3
	  npeps(k) = 0
	end do
	nint = 0
	nout = 0
	if(rionst.eq.0.)goto 401
	start = secnds(0.0)
	dfact = 3.047*sqrt(rtemp*epsout*epkt/80./298.)
	deblen = dfact/sqrt(rionst)
	print *,'debye length (A): ',deblen
	print *,' '
	if(donon)then
	  print *,'computing nonlinear ionic forces'
	else
	  print *,'computing linear ionic forces'
	end if
	print *,' '
c
	dfact = rionst*cfact
	print *,'ionic pressure factor: ',dfact
	do i = 2,igrid-1
	  do j = 2,igrid-1
	    do k = 2,igrid-1
c loop over maps finding ion surface points
c compute inward vector
		if(debmap(i,j,k).eq.1.)then
              ide(1) = (debmap(i+1,j,k) - debmap(i-1,j,k))
		  ide(2) = (debmap(i,j+1,k) - debmap(i,j-1,k))
		  ide(3) = (debmap(i,j,k+1) - debmap(i,j,k-1))
		  ibound = abs(ide(1)) + abs(ide(2)) + abs(ide(3))
		  demag = sqrt(float(ibound))
		  fmag = 0.
		  if(demag.gt.0)then
		    nout = nout + 1
		    iat = 0
		    if(ide(1).ne.0)then
		      iat = max(iat,iatmapd(i+1,j,k))
		      iat = max(iat,iatmapd(i-1,j,k))
		    end if
		    if(ide(2).ne.0)then
		      iat = max(iat,iatmapd(i,j+1,k))
		      iat = max(iat,iatmapd(i,j-1,k))
		    end if
		    if(ide(3).ne.0)then
		      iat = max(iat,iatmapd(i,j,k+1))
		      iat = max(iat,iatmapd(i,j,k-1))
		    end if
		    npeps(ibound) = npeps(ibound)+1
		    earea = earea + demag
		    phi0 = phimap(i,j,k)
		    pout = pout + phi0
		    if(donon)then
			exphi0 = exp(-phi0)
			exphi02 = exphi0*exphi0
			cex = cfact*(cion(1)*(exphi02-1.) 
     &			+ cion(2)*(exphi0-1.) 
     &			+ cion(3)*(-1.+1./exphi02) 
     &			+ cion(4)*(-1.+1./exphi0))
			cnet = cfact*(czion(1)*exphi02 + czion(2)*exphi0 + czion(3)/exphi02 + czion(4)/exphi0)
		    else
		      cex = dfact*phi0**2
		      cnet = -2.0*dfact*phi0
		    end if
cc debug
c			print *,'i,j,k,ide,demag,phi,cex,iat: ',i,j,k,ide,demag,phi0,cex,iat
cc debug
		    dpres = dpres + cex*demag
		    cout = cout + cnet
		    do i1 = 1,3
		      df(i1) = ide(i1)*cex/scale**2
			atfrci(i1,iat) = atfrci(i1,iat) - df(i1)
		    end do
		  else ! outside, but not boundary
		    phi0 = phimap(i,j,k)
		    if(donon)then
			exphi0 = exp(-phi0)
			exphi02 = exphi0*exphi0
			cnet = cfact*(czion(1)*exphi02 + czion(2)*exphi0 + czion(3)/exphi02 + czion(4)/exphi0)
		    else
		      cnet = -2.0*dfact*phi0
		    end if
		    cnet = 0. ! debug
		    cout = cout + cnet
		  end if ! boundary point branch
c
c loop thru charges and calc ionic reaction forces
c
		  if(cnet.ne.0.)then
		    srfx(1) = i
		    srfx(2) = j
		    srfx(3) = k
		    do i1 = 1,nqass
	            do k1 = 1,3
	              dr(k1) = qass(k1,i1)-srfx(k1)
	            end do
	            r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
	            if(r2.gt.1.e-6)then
		        ergi = ergi + qass(4,i1)*cnet/sqrt(r2)
	              rr3 = 1./(sqrt(r2)*r2)
	              do k1 = 1,3
			    dfrc = dr(k1)*qass(4,i1)*cnet*rr3/scale/epsout
	                atfrci(k1,indxat(i1)) = atfrci(k1,indxat(i1)) + dfrc
	              end do
		      end if
		    end do
		  end if
	      else ! inside
		  nint = nint + 1
		end if ! in/out branch
	    end do
	  end do
	end do
	dpres = dpres/earea
	earea = earea/scale**2
	pout = pout/nout
	cout = cout/(scale**3)
	ergi = ergi/scale**2/epsout/2.
401	continue
	print *,'# keps boundary points,eps area: ',npeps,earea
	print *,'Mean ion pressure (kT/A**3): ',dpres
	print *,'Mean outer boundary potl (kT/e): ',pout
	print *,'net outer boundary charge (e): ',cout
	print *,'# inner, outer keps points: ',nint,nout
	print *,'ionic rxn field energy from (surface) rho (kT): ',ergi
	fin = secnds(start)
	print *,'ionic pressure cpu time (s): ',fin
	ergt = ergs + ergi
c
c force write
c
	write(41,*)'start dielectric boundary forces (kT/A): '
	do i = 0,natom
	  write(41,'(i6,3f12.6)')i,(atfrcd(k,i),k=1,3)
	end do
	Write(41,*)'end dielectric rxn field energy from surface Q (kT): ',ergs
	write(41,*)'start Rxn field atomic forces (kT/A): '
	do i = 0,natom
	  do k = 1,3
	    atfrcr(k,i) = atfrcr(k,i)*scale/(4.*pi)
	  end do
	  write(41,'(i6,3f12.6)')i,(atfrcr(k,i),k=1,3)
	end do
	write(41,*)'end '
	write(41,*)'start ionic boundary force (kT/A): '
	do i = 0,natom
	  write(41,'(i6,3f12.6)')i,(atfrci(k,i),k=1,3)
	end do
	write(41,*)'end ionic rxn field energy from (surface) rho (kT): ',ergi
      do k = 1,3
          ftrans(k) = 0.
          frot(k) = 0.
      end do
      do i = 1,natom
        do k = 1,3
	      atfrct(k,i) = atfrcd(k,i) + atfrcr(k,i) + atfrci(k,i)
            dr(k)=atcrd(k,i)-xcen(k)
        end do
        call cross(dr,atfrct(1,i),dr1)
        do k = 1,3
            frot(k) = frot(k) + dr1(k)
        end do
        do k = 1,3
            ftrans(k) = ftrans(k) + atfrct(k,i)
        end do
      end do
      print *,'Net Frot  : ',(frot(k),k=1,3)
      print *,'Net Ftrans: ',(ftrans(k),k=1,3)
      print *,' '
c
c remove net translational, rotational forces
c
c calculate and invert inertia tensor
c
	xx = 0.
	yy = 0.
	zz = 0.
	xy = 0.
	xz = 0.
	yz = 0.
c	print *,'xx,yy,zz,xy,xz,yz: '
	do i = 1,natom
	  do k = 1,3
	    dr(k)=atcrd(k,i)-xcen(k)
	  end do
c	  write(6,'(3f8.3)')(atcrd(k,i),k=1,3)
c	  print *,'dr: ',dr
	  xx = xx + dr(1)*dr(1)
	  yy = yy + dr(2)*dr(2)
	  zz = zz + dr(3)*dr(3)
	  xy = xy + dr(1)*dr(2)
	  xz = xz + dr(1)*dr(3)
	  yz = yz + dr(2)*dr(3)
c	  print *,'xx: ',xx,yy,zz,xy,xz,yz
	end do
	tin(1,1) = yy+zz
	tin(2,2) = xx+zz
	tin(3,3) = xx+yy
	tin(1,2) = -xy
	tin(2,1) = -xy
	tin(1,3) = -xz
	tin(3,1) = -xz
	tin(2,3) = -yz
	tin(3,2) = -yz
c	print *,'xx,yy,zz,xy,xz,yz: '
c	print *,xx,yy,zz,xy,xz,yz
	print *,'Inertia tensor: '
	do j = 1,3
	  write(6,'(''('',3f15.5,'')'')')(tin(k,j),k=1,3)
	end do
c
c invert
c
	call mtdet(tin,det)
	print *,'determinant: ',det
	if(det.gt.1.e-3)then
	  call gaussj(tin,3,3,dr,1,1)
	  print *,'Inverse inertia tensor: '
	  do j = 1,3
	    write(6,'(''('',3f15.5,'')'')')(tin(k,j),k=1,3)
	  end do
	  do k = 1,3
	    froti(k) = frot(1)*tin(k,1)+frot(2)*tin(k,2)+frot(3)*tin(k,3)
	  end do
	  print *,'Inverse torque: '
	  write(6,'(''('',3f15.5,'')'')')(froti(k),k=1,3)
c
c remove net forces
c
	  do i = 1,natom
	    do k = 1,3
		dr(k) = atcrd(k,i) - xcen(k)
	    end do
	    call cross(froti,dr,dr1)
	    do k = 1,3
	      atfrct(k,i) = atfrct(k,i) - ftrans(k)/natom - dr1(k)
	    end do
	  end do
	else
	  do i = 1,natom
	    do k = 1,3
	      atfrct(k,i) = atfrct(k,i) - ftrans(k)/natom
	    end do
	  end do
	end if
c
c check net translational, rotational forces
c
	do k = 1,3
	  ftrans(k) = 0.
	  frot(k) = 0.
	end do
	do i = 1,natom
	  do k = 1,3
	    dr(k)=atcrd(k,i)-xcen(k)
	  end do
	  call cross(dr,atfrct(1,i),dr1)
	  do k = 1,3
	    frot(k) = frot(k) + dr1(k)
	  end do
	  do k = 1,3
		ftrans(k) = ftrans(k) + atfrct(k,i)
	  end do
	end do
	print *,'Residual Frot  : ',(frot(k),k=1,3)
	print *,'Residual Ftrans: ',(ftrans(k),k=1,3)
	print *,' '
c	write(41,*)'start total solvation forces (kT/A): '
	write(41,*)'start total solvation forces (kcal/mole/A): '
	do i = 1,natom
	  do k = 1,3
	    atfrct(k,i) = calk*atfrct(k,i)
	  end do
	  write(41,'(i6,3f12.6)')i,(atfrct(k,i),k=1,3)
	end do
	write(41,*)'end '
	close(41)


	fin = cputime(start)
	print *,'Dielectric Reaction forces took ',fin,' sec'
	start = cputime(0.0)
	return
	end

	subroutine gridphi(scale,igrid,srfx,srfq,qatom,ergsgd,dfrcgd)
c
c map atom charge onto grid just like in qmap assignment,
c then compute interaction with reaction charge
c
	real*4 srfx(3),qatom(4),ergsgd,dfrcgd(3)
	real*4 dr(3),qgrid(4,8)
	integer kb(3),kt(3)
	real*4  dg(3),cg(3),chrgv
c---------------------------------------------------
      iout = 0
	rgrid = igrid
	do k = 1,3
	  if((qatom(k).le.1.).or.(qatom(k).ge.rgrid)) iout = 1
	end do
	if(iout.eq.1)then
c
c just use actual atom charge position
c
	  do k1 = 1,3
	    dr(k1) = qatom(k1)-srfx(k1)
	  end do
	  r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
	  if(r2.gt.1.e-6)then
	    ergsgd = qatom(4)*srfq/sqrt(r2)
	    rr3 = 1./(sqrt(r2)*r2)
	    do k1 = 1,3
	      dfrcgd(k1) = dr(k1)*qatom(4)*srfq*rr3
	    end do
	  endif
	else
	  chrgv = qatom(4)
        do i = 1,3
c
c trucate to nearest grid point
c find upper grid point of box
c find position of charge in box in terms
c of fractional distance along edges
c
          kb(i) = qatom(i)
          kt(i) = kb(i) + 1
          dg(i) = qatom(i) - kb(i)
          cg(i) = 1-dg(i)
        end do
c
c put fractional charges at 8 corners
c
        qgrid(4,1) = cg(1)*cg(2)*cg(3)*chrgv
        qgrid(1,1) = kb(1)
        qgrid(2,1) = kb(2)
        qgrid(3,1) = kb(3)
c
        qgrid(4,2) = dg(1)*cg(2)*cg(3)*chrgv
        qgrid(1,2) = kt(1)
        qgrid(2,2) = kb(2)
        qgrid(3,2) = kb(3)
c
        qgrid(4,3) = cg(1)*dg(2)*cg(3)*chrgv
        qgrid(1,3) = kb(1)
        qgrid(2,3) = kt(2)
        qgrid(3,3) = kb(3)
c
        qgrid(4,4) = dg(1)*dg(2)*cg(3)*chrgv
        qgrid(1,4) = kt(1)
        qgrid(2,4) = kt(2)
        qgrid(3,4) = kb(3)
c
        qgrid(4,5) = cg(1)*cg(2)*dg(3)*chrgv
        qgrid(1,5) = kb(1)
        qgrid(2,5) = kb(2)
        qgrid(3,5) = kt(3)
c
        qgrid(4,6) = dg(1)*cg(2)*dg(3)*chrgv
        qgrid(1,6) = kt(1)
        qgrid(2,6) = kb(2)
        qgrid(3,6) = kt(3)
c
        qgrid(4,7) = cg(1)*dg(2)*dg(3)*chrgv
        qgrid(1,7) = kb(1)
        qgrid(2,7) = kt(2)
        qgrid(3,7) = kt(3)
c
        qgrid(4,8) = dg(1)*dg(2)*dg(3)*chrgv
        qgrid(1,8) = kt(1)
        qgrid(2,8) = kt(2)
        qgrid(3,8) = kt(3)
	  ergsgd = 0.
	  do k1 = 1,3
	    dfrcgd(k1) = 0.
	  end do
	  do i = 1,8
	    do k1 = 1,3
	      dr(k1) = qgrid(k1,i)-srfx(k1)
	    end do
	    r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
	    if(r2.gt.1.e-6)then
	      ergsgd = ergsgd + qgrid(4,i)*srfq/sqrt(r2)
	      rr3 = 1./(sqrt(r2)*r2)
	      do k1 = 1,3
	        dfrcgd(k1) = dfrcgd(k1) + dr(k1)*qgrid(4,i)*srfq*rr3
	      end do
	    endif
	  end do
	end if
				        
	return
	end
