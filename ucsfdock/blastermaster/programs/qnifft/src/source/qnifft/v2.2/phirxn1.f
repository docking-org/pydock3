	subroutine phirxn(ianal,epsin,epsout,natom,epkt,rtemp,frcfil)
c 
c compute rxn potential/field at target site locations
c
c-------------------------------------------------
	include 'qdiffpar.h'
c-------------------------------------------------
c
	real*4 atphix(natom) ! rxn potl
	real*4 atfldx(3,natom) ! rxn field
	real*4 xo(3,natom),chrgv(natom),xyz(3),xg(3)
	character*15 atlab(natom)
	real*4 dr(3),srfx(3)
	real*4 xcen(3),dr1(3)
	real*4 earea
	integer npeps(3),ide(3)
	character*8 hour
	character*80 frcfil
	character*132 line
	logical ianal
c-------------------------------------------------------
	if(.not.ianal)return
	if(epsin.eq.epsout)then
	  print *,'inner, outer dielectrics the same, omitting reaction potl. calc.'
	  return
	endif
	start = cputime(0.0)
c	call time(hour)
	print *,' '
c	write(6,*)'computing rxn field/potentials at: ',hour
c-------------------------------------------------------
	pi = 355./113.
c
c read frcfil to get coords, charges, atom labels of sites
c scale to grid coords
c
	qfix = 0.
	nqass = 0
	do k = 1,3
	  xcen(k) = 0.
	end do
	open(unit=16,file=frcfil)
	open(42,file='reaction.fld')
	do i = 1,12
	  read(16,'(a)')line
c	  print *,line
	  if(i.eq.11)then
	    write(42,'(a)')'         coordinates          charge   reaction potential  field (kT/e/Ang.)'
	  else
	    write(42,'(a)')line
	  end if
	end do
	do i = 1,natom
	  read(16,230)xyz,chrgv(i),phiv,fx,fy,fz,atlab(i)
230     format(8f11.4,1x,a15)
	  qfix = qfix + chrgv(i)
	  call ctog(xyz,xg)
	  do k = 1,3
	    xo(k,i) = xg(k)
	  end do
	  do k = 1,3
	    xcen(k) = xcen(k) + xo(k,i)
	  end do
	  if(chrgv(i).ne.0.)then
	    nqass = nqass + 1
	  end if
	end do
	print *,'total atomic charge in site file: ',qfix
	do k = 1,3
	  xcen(k) = xcen(k)/natom
	end do
	print *,'geometric center: ',xcen
c
c init. field, potl
c
	do i = 1,natom
	  atphix(i) = 0.
	  do k = 1,3
	    atfldx(k,i) = 0.
	  end do
	end do
c
c surface charge reaction field / potential
c
	qsurf = 0.
	nmove = 0
	nsurf = 0
	nsurfq = 0
	nbound = 0
	do k = 1,3
	  npeps(k) = 0
	end do
	earea = 0.
	dpres = 0.
	qfix = 0.
	print *,'using iatmap to correct surface points...'
	do i = 2,igrid-1
	  do j = 2,igrid-1
	    do k = 2,igrid-1
		qfix = qfix + qmap(i,j,k)
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
		  if(iat.ne.0)then !found a neighbor atom
		    dmin2 = (atcrd(1,iat)-srfx(1))**2 + (atcrd(2,iat)-srfx(2))**2 + (atcrd(3,iat)-srfx(3))**2
		    if(dmin2.gt.1.e-6)then
			nmove = nmove + 1
		      srat = atrad(iat)/sqrt(dmin2)
c			print *,'iat, atrad, srat: ',iat,atrad(iat), srat
			srfx(1) = atcrd(1,iat) + srat*(srfx(1)-atcrd(1,iat))
			srfx(2) = atcrd(2,iat) + srat*(srfx(2)-atcrd(2,iat))
			srfx(3) = atcrd(3,iat) + srat*(srfx(3)-atcrd(3,iat))
		    end if
		  end if
c		  print *,'srfxyz after correction: ',srfx
c
c loop thru charges and calc reaction forces/ potls.
c
		  do i1 = 1,natom
	          do k1 = 1,3
	            dr(k1) = xo(k1,i1)-srfx(k1)
	          end do
c		    print *,'xo,dr: ',xo(1,i1),xo(2,i1),xo(3,i1),dr
	          r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
	          if(r2.gt.1.e-6)then
			atphix(i1) = atphix(i1) + srfq/sqrt(r2)
	            rr3 = 1./(sqrt(r2)*r2)
	            do k1 = 1,3
	              atfldx(k1,i1) = atfldx(k1,i1) + dr(k1)*srfq*rr3
	            end do
		    end if
		  end do
	      end if
          end do   
        end do      
	end do
	qfix = qfix/4./pi/scale
	print *,'grid fixed charge of : ',qfix
	qsurf = qsurf/(4*pi*scale*epkt)
	qsurfe = (epsin-epsout)*qfix/epsin/epsout/epkt
	print *,'surface points readjusted: ',nmove
	print *,'surface points with, without fixed charge: ',nsurfq,nsurf
	print *,'expected, actual induced surface charge: ',qsurfe,qsurf
	earea = earea/scale**2
	print *,'number of boundary points: ',nbound
	print *,'# eps boundary points,eps area (A**3): ',npeps,earea
	ergs = 0.
c
c force write
c
	do i = 1,natom
	  atphix(i) = atphix(i)/4./pi
	  ergs = ergs + chrgv(i)*atphix(i)
	  do k = 1,3
	    atfldx(k,i) = atfldx(k,i)*scale/4./pi
	    xg(k) = xo(k,i)
	  end do
	  call gtoc(xg,xyz)
	  write(42,230)xyz,chrgv(i),atphix(i),(atfldx(k,i),k=1,3),atlab(i)
	end do
	write(6,*)' Interaction between reaction potl. and target charges (kT) ',ergs
	write(42,*)' Interaction between reaction potl. and target charges (kT) ',ergs
	write(42,*)' (Excludes Coulombic energy)'
	close(42)
	fin = cputime(start)
	print *,'Reaction field/potl. took ',fin,' sec'
	return
	end
