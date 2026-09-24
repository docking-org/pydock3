	subroutine phierg(ianal,donon,rionst,epsin,epsout,natom,iconc,epkt,rtemp,fasteps)
c analyse phimap energies etc
c handles mixed salt cases- oct 99 VM/KAS
c
c-------------------------------------------------
	include 'qdiffpar.h'
c-------------------------------------------------
c
	dimension cconc(4), tconc(4), tcrg(4)
	dimension trhophi(4), tosm(4)
	integer nepsatm(0:natmx) ! kas debug 13-feb-01
	dimension inear(6)
	character*8 hour
	real*8 ergs,qsurf,sumrr,corrmn
	real*4 corrmx
	logical donon,iconc,ianal,fasteps
      data fpi / 12.566 /
      data fact  / 6.0e-4  /   !converts moles to ions/cubic angstrom
					 ! ==Na,l/moles,Ang*3
c	data conv / 1.78e-3  /	 !converts integration to units of e
c					 !==erg,cm,e,e/kT,Ang,esu,esu
c-------------------------------------------------------
	start = cputime(0.0)
c	call time(hour)
	print *,' '
c	write(6,*)'analysing phimap at: ',hour
c-------------------------------------------------------
	midg = (igrid+1)/2
c
c total of q.phi, includes grid energy
c
	pi = 355./113.
	terg = 0.
	qfix = 0.
	do k = 1,igrid
	  do j = 1,igrid
	    do i = 1,igrid
		qfix = qfix + qmap(i,j,k)
		terg = terg + qmap(i,j,k)*phimap(i,j,k)
	    end do
	  end do
	end do
	terg = terg/8./pi/scale
	print *,'0.5*sum(phi.q) (kT) =  ',terg
	print *,'Interaction of source charges with themselves (kT)',terg
	print *,'(Includes Grid+Coulombic+Solvation energy)'
	print *,'Equals Total energy for Linear PB only'
	qfix = qfix/4./pi/scale
	print *,'For net fixed charge of : ',qfix
	if(.not.ianal)return
	if(rionst.le.0.)donon=.false.
c
c integrate normal field at edge of box to get net charge 
c
	eint = 0.0
	darea = 1./scale**2
	do i = 2,igrid-1
	  do k = 2,igrid-1
	    eint = eint + (phimap(k,1,i) - phimap(k,3,i))
	    eint = eint + (phimap(k,igrid,i) - phimap(k,igrid-2,i))
	    eint = eint + (phimap(1,k,i) - phimap(3,k,i))
	    eint = eint + (phimap(igrid,k,i) - phimap(igrid-2,k,i))
	    eint = eint + (phimap(k,i,1) - phimap(k,i,3))
	    eint = eint + (phimap(k,i,igrid) - phimap(k,i,igrid-2))
	  end do
	end do
	qint = -eint*epsout*darea*scale/fpi/2.
	print *,'net charge from gauss integral: ',qint
	if(rionst.gt.0.)then
c
c initialize energy density integrals
c
	do il = 1,4
          trhophi(il) = 0.0
          tosm(il) = 0.0
          tcrg(il) = 0.0
          tconc(il) = 0.0
        end do
	  tcions = 0.0
	  tcrhophi = 0.0
	  tedotd = 0.0
	  tcosm = 0.0
	  tcpos = 0.
	  tcneg = 0.
        vol = 1./scale**3
c cal z range indicies
	  zlw = (2-midg)/scale+oldmid(3)
	  zup = (igrid-1-midg)/scale+oldmid(3)
	  print *,'integrating over zrange: ',zlw,zup
	  do k = 2,igrid-1
	    do j = 2,igrid-1
	      do i = 2,igrid-1
	        if(debmap(i,j,k).gt.0.)then
	          pot = phimap(i,j,k)
c charge concentrations 
		    if(donon)then 
		      do il = 1, 4
		        cconc(il) = cion(il)*(exp(-zion(il)*pot)-1.)
			  tconc(il) = tconc(il) + cconc(il)
			  tcrg(il) = tcrg(il) + zion(il)*cconc(il)
c entropy contributions
			  trhophi(il) = trhophi(il) + pot*zion(il)*cconc(il)
c osmotic pressure terms
			  tosm(il) = tosm(il) + cconc(il)
		      end do
		    else
		      do il = 1, 4
		        cconc(il) = cion(il)*(-zion(il)*pot)
			  tconc(il) = tconc(il) + cconc(il)
			  tcrg(il) = tcrg(il) + zion(il)*cconc(il)
c entropy contributions
			  trhophi(il) = trhophi(il) + pot*zion(il)*cconc(il)
c osmotic pressure terms
			  tosm(il) = tosm(il) + cconc(il)
		      end do
		    end if
	        end if
c
c electrostatic stress term (E.D)
c
	        iin = (iepsmp(i,j,k,1) + iepsmp(i,j,k,2) 
     &      + iepsmp(i,j,k,3) + iepsmp(i-1,j,k,1)
     &      + iepsmp(i,j-1,k,2) + iepsmp(i,j,k-1,3))
c use eps map
	        if(iin.eq.0) then
		    fldx = (phimap(i+1,j,k) - phimap(i-1,j,k))/2.
		    fldy = (phimap(i,j+1,k) - phimap(i,j-1,k))/2.
		    fldz = (phimap(i,j,k+1) - phimap(i,j,k-1))/2.
		    tedotd = tedotd + (fldx**2 + fldy**2 + fldz**2)
	        end if
            end do   
          end do      
	  end do

	  do il = 1,4
	    tcrg(il) = tcrg(il)*vol*fact
	    tosm(il) = tosm(il)*vol*fact
	    trhophi(il) = trhophi(il)*fact*vol
	    tcions = tcions + tcrg(il)
	    tcosm = tcosm + tosm(il)
	    tcrhophi = tcrhophi + trhophi(il)
	  end do

	  tcedotd = tedotd*vol*scale*scale*epsout/fpi
	  tcpos = (tcrg(1) + tcrg(2))
	  tcneg = (tcrg(3) + tcrg(4))
	  ergion = tcedotd/2. - tcrhophi - tcosm
	  qnrat = tcneg/abs(tcions)
	  qprat = tcpos/abs(tcions)
  
	  calpkt = 1.99*rtemp/1000.
	  print *,'kT = ',calpkt,' kcal at ',rtemp,'K'
	  print *,' '
	  print *,'Q net                       : ',tcions
	  print *,'Q+ excess,ratio             : ',tcpos,qprat
	  print *,'Q- excess,ratio             : ',tcneg,qnrat
	  print *,' '
	  print *,'rho.phi or TdS          (kT): ',tcrhophi
	  print *,'   rho.phi(++)             (kT): ',trhophi(1)
	  print *,'   rho.phi(+)              (kT): ',trhophi(2)
	  print *,'   rho.phi(--)             (kT): ',trhophi(3)
	  print *,'   rho.phi(-)              (kT): ',trhophi(4)
	  print *,'E.D                     (kT): ',tcedotd
	  print *,'dPI                     (kT): ',tcosm
	  print *,'   dPI(++)                 (kT): ',tosm(1)
	  print *,'   dPI(+)                  (kT): ',tosm(2)
	  print *,'   dPI(--)                 (kT): ',tosm(3)
	  print *,'   dPI(-)                  (kT): ',tosm(4)
	  print *,' '
	  print *,'E.D/2 - rho.phi - dPI   (kT): ',ergion
	  print *,' '
	fin = cputime(start)
	print *,'Salt analysis took ',fin,' sec'
	start = cputime(0.0)
	end if
c
c collect assigned charge
c
	nqass = 0
	qfix = 0.
	rmidg = (igrid+1)/2.
	do i = 1,natom
	  if(atcrg(i).ne.0.)then
	    nqass = nqass + 1
	    do k = 1,3
		qass(k,nqass) = atcrd(k,i)
c		qass(k,nqass) = (atcrd(k,i)-rmidg)/scale+oldmid(k)
	    end do
	    qfix = qfix + atcrg(i)
	    qass(4,nqass) = atcrg(i)
c          print *,qass(4,nqass),qass(1,nqass),qass(2,nqass),qass(3,nqass)
	  end if
	end do
	print *,'total atomic charge: ',qfix
c
c coulombic energy
c
	ergc = 0.
	do i = 1,nqass
	  do j = i+1,nqass
	    r2 = (qass(1,i)-qass(1,j))**2 + (qass(2,i)-qass(2,j))**2
     &    + (qass(3,i)-qass(3,j))**2
c	    print *,'r2: ',r2
	    if(r2.gt.1.e-6)ergc = ergc + qass(4,i)*qass(4,j)/sqrt(r2)
	  end do
	end do
c	print *,'epsin: ',epsin
	ergc = ergc/epsin*scale
	print *,'Coulombic energy in uniform Epsin (kT): ',ergc
	fin = cputime(start)
	print *,'Coulombic energy took ',fin,' sec'
	start = cputime(0.0)
c
c surface charge reaction field energy
c
cc debug kas 13-feb-01 check iepsmp2
c	do i = 0,natom
c	  nepsatm(i) = 0
c	end do
cc debug kas 13-feb-01 check iepsmp2
	ergs = 0.
	sumrr = 0.
	corrmn = 0.
	corrmx = 0.
	qsurf = 0.
	sfact = 1./4./pi/scale/epkt
	nmove = 0
	nsurf=0
	nsurfq = 0
	print *,'using list of corrected boundary points ...',nbnd
	ibnd = 0
	do k = 2,igrid-1
	  do j = 2,igrid-1
	    do i = 2,igrid-1
	      iin = (iepsmp(i,j,k,1) + iepsmp(i,j,k,2) 
     &      + iepsmp(i,j,k,3) + iepsmp(i-1,j,k,1)
     &      + iepsmp(i,j-1,k,2) + iepsmp(i,j,k-1,3))
	      if((iin.ne.0).and.(iin.ne.6)) then
		  if(qmap(i,j,k).ne.0.)then
		    nsurfq = nsurfq+1
		  else
		    nsurf = nsurf+1
		  end if
		  srfq = 6.*phimap(i,j,k)-(phimap(i+1,j,k)+phimap(i-1,j,k)
     &	  + phimap(i,j+1,k) + phimap(i,j-1,k)
     &	  + phimap(i,j,k+1) + phimap(i,j,k-1) + qmap(i,j,k)/epsin)
		   qsurf = qsurf + srfq
c
c correction position of surface charge using position of nearest atom
c iepsmp2 contains indices of atoms forming surface
c
		  if(fasteps)then
		    srfx = i
		    srfy = j
		    srfz = k
		    srferr = 0.
                iat = iatmap(i,j,k)
                if(iat.ne.0)then !found a neighbor atom
                  dmin2 = (atcrd(1,iat)-srfx)**2
     &            + (atcrd(2,iat)-srfy)**2
     &            + (atcrd(3,iat)-srfz)**2
                  if(dmin2.gt.1.e-6)then
		        ibnd = ibnd + 1
                    srat = atrad(iat)*scale/sqrt(dmin2)
                    srfx = atcrd(1,iat) + srat*(srfx-atcrd(1,iat))
                    srfy = atcrd(2,iat) + srat*(srfy-atcrd(2,iat))
                    srfz = atcrd(3,iat) + srat*(srfz-atcrd(3,iat))
		        srferr = (srfx-i)**2 + (srfy-j)**2 + (srfz-k)**2
                  end if
                end if
		  else
		    ibnd = ibnd + 1
		    srfx = xyzbnd(1,ibnd)
		    srfy = xyzbnd(2,ibnd)
		    srfz = xyzbnd(3,ibnd)
		    srferr = (srfx-i)**2 + (srfy-j)**2 + (srfz-k)**2
		  end if
		  corrmn = corrmn + srferr
		  corrmx = max(corrmx,srferr)
c		  print *,'srfxyz after correction: ',i,j,k,srfx,srfy,srfz
c
c loop thru charges and calc reaction energy
c
		  do i1 = 1,nqass
	    	    r2 = (qass(1,i1)-srfx)**2 + (qass(2,i1)-srfy)**2 
     &             + (qass(3,i1)-srfz)**2
	          if(r2.gt.1.e-6)ergs = ergs + qass(4,i1)*srfq/sqrt(r2)
	          if(r2.gt.1.e-6)sumrr = sumrr + 1./sqrt(r2)
		  end do
	      end if
          end do   
        end do      
	end do
cc debug
c	print *,'in Phierg, after surface moving iat, nepsatm '
c	ntot = 0
c	do i = 0,natom
c	  ntot = ntot + nepsatm(i)
c	  print *,i,nepsatm(i)
c	end do
c	print *,'ntot: ',ntot
cc debug
	qsurf = qsurf*sfact
	qsurfe = (epsin-epsout)*qfix/epsin/epsout/epkt
	print *,'readjusted surface points: ',ibnd
	print *,'surface points with, without fixed charge: ',nsurfq,nsurf
	print *,'expected, actual induced surface charge: ',qsurfe,qsurf
	ergs = ergs*epkt*sfact*scale/2.
	print *,'Dielectric rxn field energy from surface Q (kT): ',ergs
	sumrr = sumrr*scale/ibnd
	print *,'sum 1/rr q-surface: ',sumrr
	corrmn = sqrt(corrmn/ibnd)
	corrmx = sqrt(corrmx)
	print *,'mean surface position correction (A): ',corrmn
	print *,'max  surface position correction (A): ',corrmx
	fin = cputime(start)
	print *,'Dielectric Reaction energy took ',fin,' sec'
	start = cputime(0.0)
	if(iconc)then
	  print *,'converting potentials to Molar concentrations'
	  do k = 1,igrid
	    do j = 1,igrid
	      do i = 1,igrid
	        if(debmap(i,j,k).gt.0.)then 
		     cnet1 = cion(1)*(exp(-zion(1)*phimap(i,j,k))-1)
		     cnet2 = cion(2)*(exp(-zion(2)*phimap(i,j,k))-1)
		     cnet3 = cion(3)*(exp(-zion(3)*phimap(i,j,k))-1)
		     cnet4 = cion(4)*(exp(-zion(4)*phimap(i,j,k))-1)
		    else
		     cnet1 = 0.
		     cnet2 = 0.
		     cnet3 = 0.
		     cnet4 = 0.
	        end if
		  phimap(i,j,k) = 2*cnet1+cnet2-2*cnet3-cnet4
	      end do
	    end do
	  end do
	end if
	return
	end
