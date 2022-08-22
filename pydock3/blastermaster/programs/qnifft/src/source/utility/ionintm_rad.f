      program ionintm
c kim sharp, 7 jan 1992
c integrates various ion atmosphere quantities
c from potential output of delphi
c combines ionat3 method for E.D integration, and
c ionat4 for rhophi, dPI terms
c includes double layer electrical
c entropic, E.D, osmotic pressure integrations
c and donnan coefficients
c also handles a series of focussing maps
c kas feb 01 modified to handle mixed divalent/monovalents
c kas jan 07 modified to compute radial dist. fcns. of ions
c---------------------------------------------------------
	parameter (nbin = 129 )
      parameter ( isize = 65 )    !grid size
      dimension reps(isize,isize,isize,3)    !fine eps map
      dimension rkeps(isize,isize,isize)    !fine debye map
      dimension phi(isize,isize,isize)	     !potential map
      dimension oldmid(3),oldmid1(3),oldmid2(3) !scale infoo
      dimension xx(3),xlw(3),xlw1(3),xup(3),xup1(3)
	dimension cion(4),zion(4),tosm(4),trhophi(4)
	dimension tcrg(4),cconc(4)
	real*4 radcorr(nbin),dvol
	real*4 radion(7,nbin)
	integer ibin,jbin(nbin)
      character*80 phiname
	character*10 nlbl
      character*20 toplbl
      character*16 botlbl
      character*60 philbl
	logical donon
c--------------------------------------------------------------
      data twopi / 6.2832 /
      data fpi / 12.566 /
	data epsout / 80. /
      data fact  / 6.0e-4  /   !converts moles to ions/cubic angstrom
					 ! ==Na,l/moles,Ang*3
	data conv / 1.78e-3  /	 !converts integration to units of e
					 !==erg,cm,e,e/kT,Ang,esu,esu
	data zion / 2., 1., -2., -1. /
c--------------------------------------------------------------
	print *,'kim sharp, 2 feb 2001'
	print *,'integrates various ion atmosphere quantities'
	print *,'from potential output of delphi'
	print *,'combines ionat3 method for E.D integration, and'
	print *,'ionat4 for rhophi, dPI terms'
	print *,'includes, dS, E.D, osmotic pressure integrations'
	print *,'and donnan coefficients'
	print *,'also handles a series of focussing maps'
	print *,'start with innermost map and work outwards'
	print *,'if molecule is cylindrical, must align it on z axis'
c	print *,'using debye map for E.D integration mask like ionat4'
	print *,'using eps map for E.D integration mask like ionat3'
	print *,'handles mixed salt'
	print *,'c kas jan 07 modified to compute radial dist. fcns. of ions'
	print *,'grid size: ',isize
c--------------------------------------------------------------
	ifoc = 0
	izlck = 0
	do i = 1,nbin
	  jbin(i) = 0
	  radcorr(i) = 0
	  do j = 1,7
	    radion(j,i) = 0.
	  end do
	end do
c----------------------------
101   print *,'C++,C+,C--,C- (M)?'
      read(5,*,err=101)cion
	cbulk = 0.
	do il =1,4
	  cbulk = cbulk + 0.5*cion(il)*zion(il)**2
	end do
	deblen = 3.047/sqrt(cbulk+0.0000001) ! kas 1-feb-01
	print *,'debye length: ',deblen
c-----------------------------
102	print *,'# linear PB ? 1=yes>>'
	read(5,*,err=102)nterm
	donon = .true.
	if(nterm.eq.1)donon = .false.
300	continue
c
c read in compact eps and debye maps, expand and make coarse file
c
	call geteps(gpa2,oldmid2,reps,rkeps)
	print *,' scale, molc. centre coordinates'
	print *,gpa2,oldmid2
c
c read in potential file
c
100	print *,' enter name of potential map file (w/o .phi)'
	read(5,'(a)')phiname
      il = lnblnk(phiname)
      open(unit=40,file=phiname(:il)//'.phi',status='old',
     &form='unformatted',err=100)
	print *,' reading potential file...'
	read(40)toplbl
	print *,toplbl
	read(40)nlbl,philbl
	print *,' name of phi map is:'
	write(6,'(a)')philbl
	read(40)phi
	read(40)botlbl
	read(40)gpa,oldmid
	close(40)
	print *,botlbl
	print *,' scale, molc. centre coordinates from *.phi file'
	print *,gpa,oldmid
c check scales are the same
	isame = 1	
	if(abs(gpa2-gpa).gt.0.0005)isame = 0
	do k = 1,3
	  if(abs(oldmid2(k)-oldmid(k)).gt.0.0005)isame = 0
	end do
	if(isame.eq.0)print *,'WARNING: eps & phi map scales NOT the same!'
	if(ifoc.eq.1)then
c check that new map completely encloses oldmap
	  ienc = 1
	  do k = 1,3
	    xlw(k) = oldmid(k)-0.5*(isize-1.)/gpa
	    xlw1(k) = oldmid1(k)-0.5*(isize-1.)/gpa1
	    xup(k) = oldmid(k)+0.5*(isize-1.)/gpa
	    xup1(k) = oldmid1(k)+0.5*(isize-1.)/gpa1
	    if(xup(k).le.xup1(k))then
		ienc = 0
	    end if
	    if(xlw(k).ge.xlw1(k))then
		ienc = 0
	    end if
	  end do
	  if(ienc.ne.1)then
	    print *,' '
	    print *,'outer map doesnt enclose previous map'
	    print *,'old lows: ',xlw1
	    print *,'new lows: ',xlw
	    print *,'old ups: ',xup1
	    print *,'new ups: ',xup
	    print *,' '
	    stop
	  end if
	else
c save z range of 1st map
	  zlwsv = oldmid(3)-0.5*(isize-2.)/gpa
	  zupsv = oldmid(3)+0.5*(isize-2.)/gpa
	end if
c
c integrate normal field at edge of box to get charge density
c by gauss's law- do x,y edges only also
c
	print *,' integrating normal field flux...'
	eint = 0.0
	exyint = 0.0
	darea = 1./(gpa*gpa)
	do i = 2,isize-1
		do k = 2,isize-1
		  exyint = exyint + (phi(k,1,i) - phi(k,3,i))
		  exyint = exyint + (phi(k,isize,i) - phi(k,isize-2,i))
		  exyint = exyint + (phi(1,k,i) - phi(3,k,i))
		  exyint = exyint + (phi(isize,k,i) - phi(isize-2,k,i))
		end do
	end do
	print *,' integrating normal cylindrical field flux...'
	do i = 2,isize-1
		do k = 2,isize-1
		  eint = eint + (phi(k,1,i) - phi(k,3,i))
		  eint = eint + (phi(k,isize,i) - phi(k,isize-2,i))
		  eint = eint + (phi(1,k,i) - phi(3,k,i))
		  eint = eint + (phi(isize,k,i) - phi(isize-2,k,i))
		  eint = eint + (phi(k,i,1) - phi(k,i,3))
		  eint = eint + (phi(k,i,isize) - phi(k,i,isize-2))
		end do
	end do
	qint = eint*epsout*darea*conv*gpa/fpi/2.
	qxyint = exyint*epsout*darea*conv*gpa/fpi/2.
	print *,' '
	print *,'full, cyl int(E) residual charge: ',qint,qxyint
c
c initialize integrals
c
	do il = 1,4
	  trhophi(il) = 0.0
	  tosm(il) = 0.0
	  tcrg(il) = 0.0
	enddo
      tcions = 0.0
	ttds = 0.0
	tedotd = 0.0
	tedotde = 0.0
	npoint = 0
	nfpoint = 0
      vol = 1./(gpa*gpa*gpa)
	midg = (isize+1)/2
	print *,'grid volume: ',vol
c cal z range indicies
	if((ifoc.eq.0).or.(izlck.eq.0))then
	  klw = 2
	  kup = isize-1
	else
	  izrnge = (zupsv - zlwsv)*gpa/2.
	  klw = midg - izrnge
	  kup = midg + izrnge
	end if
	zlw = (klw-midg)/gpa+oldmid(3)
	zup = (kup-midg)/gpa+oldmid(3)
	print *,'integrating over zrange: ',zlw,zup
	print *,'layers : ', klw,kup
	do k = klw,kup
c	    print *,'layer: ',k ! debug
	    xx(3) = (k-midg)/gpa + oldmid(3)
	    do j = 2,isize-1
	      xx(2) = (j-midg)/gpa + oldmid(2)
		do i = 2,isize-1
	        xx(1) = (i-midg)/gpa + oldmid(1)
c
c find how inside this point is, if an outer map, and
c point is inside previous map, ignore
c
		  iwant = 1
		  if(ifoc.eq.1)then
		    iwant = 0
		    do m = 1,3
			if((xx(m).gt.xup1(m)).or.(xx(m).lt.xlw1(m)))iwant = 1
		    end do
		  end if
		  if(iwant.eq.1)then
		    rad = (i - midg)**2 + (j - midg)**2
		    rad = sqrt(rad)/gpa
		    ibin = int(rad) + 1
c		    if((ibin.le.nbin).and.(k.eq.midg))then
		    if((ibin.le.nbin))then
		      radcorr(ibin) = radcorr(ibin) + vol
c		      radcorr(ibin) = radcorr(ibin) + 1 ! debug
c		      if(rad.le.2.0)then
c		        print *,'i,j,k,rad: ',i,j,k,rad
c		      end if
		    end if
		  end if
c
		  if(iwant.eq.1.and.rkeps(i,j,k).gt.0.)then
		    pot = phi(i,j,k)
		    pot = max(pot,-20.) ! kas 1-feb-01
		    pot = min(pot,20.) ! kas 1-feb-01
		    if(donon)then
			do il = 1,4
			  cconc(il) = cion(il)*(exp(-zion(il)*pot)-1.)
			  tcrg(il) = tcrg(il) + zion(il)*cconc(il)
			end do
		    else
			do il = 1,4
			  cconc(il) = cion(il)*(-zion(il)*pot)
			  tcrg(il) = tcrg(il) + zion(il)*cconc(il)
			end do
		    end if
c
c osmotic pressure term
c enthalpy and entropy contributions
c
		    do il  = 1,4
			trhophi(il) = trhophi(il) + pot*zion(il)*cconc(il)
			tosm(il) = tosm(il) + cconc(il)
		    end do
c
c rad. dist. fcn of ions
c
		    rad = (i - midg)**2 + (j - midg)**2
		    rad = sqrt(rad)/gpa
		    ibin = int(rad) + 1
		    if(ibin.le.nbin)then
		      jbin(ibin) = jbin(ibin) + 1
			do il = 1,4
			  radion(il,ibin) = radion(il,ibin) + (cconc(il) + cion(il))*fact*vol
			end do
		    end if
		  end if
c
c electrostatic stress term (E.D)
c
c use eps map
		  temp = 6. - reps(i,j,k,1) - reps(i,j,k,2) 
     &             - reps(i,j,k,3) - reps(i-1,j,k,1) 
     &             - reps(i,j-1,k,2) - reps(i,j,k-1,3) 
		  ffracv = iwant*temp/6.
		  if(ffracv.eq.1.) then
		    fldx = (phi(i+1,j,k) - phi(i-1,j,k))
		    fldy = (phi(i,j+1,k) - phi(i,j-1,k))
		    fldz = (phi(i,j,k+1) - phi(i,j,k-1))
		    tedotde = tedotde + (fldx**2 + fldy**2 + fldz**2)
		  end if
            end do   
          end do      
	end do
c	print *,'end of integration' ! debug

	do il = 1,4
	  tcrg(il) = tcrg(il)*vol*fact
	  tosm(il) = tosm(il)*vol*fact
	  trhophi(il) = trhophi(il)*fact*vol
	  tcions = tcions + tcrg(il)
	end do
	if(ifoc.eq.0)then
	  gtedotde = tedotde*vol*conv*gpa*gpa*epsout/fpi/4.
	  gtcions = tcions
	  gtrhophi = trhophi(1)+trhophi(2)+trhophi(3)+trhophi(4)
	  gttds = gtrhophi
	  gtosm = tosm(1)+tosm(2)+tosm(3)+tosm(4)
	  gtcpos = tcrg(1)+tcrg(2)
	  gtcneg = tcrg(3)+tcrg(4)
	else
	  gtedotde = gtedotde + tedotde*vol*conv*gpa*gpa*epsout/fpi/4.
	  gtcions = gtcions + tcions
	  gtrhophi = gtrhophi + trhophi(1)+trhophi(2)+trhophi(3)+trhophi(4)
	  gttds = gtrhophi
	  gtosm = gtosm + tosm(1)+tosm(2)+tosm(3)+tosm(4)
	  gtcpos = gtcpos + tcrg(1)+tcrg(2)
	  gtcneg = gtcneg + tcrg(3)+tcrg(4)
	end if
	ergion = gtedotde/2. - gtrhophi - gtosm
	qnrat = gtcneg/abs(gtcions+0.000001) ! kas 2-feb-01 fix / by 0
	qprat = gtcpos/abs(gtcions+0.000001) ! kas 2-feb-01 fix / by 0

	print *,' '
	print *,'Q net                      i: ',gtcions
	print *,'Q+ excess,ratio            i: ',gtcpos,qprat
	print *,'Q- excess,ratio            i: ',gtcneg,qnrat
	print *,'TdS  (kT)                  i: ',gttds
	print *,'rho.phi                    i: ',gtrhophi
	print *,'E.D by epsmap              i: ',gtedotde
c	print *,'E.Dby debye map            i: ',gtedotd
	print *,'dPI                        i: ',gtosm
	print *,'E.D/2 - rho.phi - dPI      i: ',ergion
	print *,' '

c
c get next outer map if required
c
	print *,'continue with next outermost map? 1:yes'
	read(5,'(i1)')iopt
	if(iopt.eq.1)then
	  if(izlck.eq.0)then
	    print *,'clip z integration to original box? 1: yes'
	    read(5,'(i1)')izlck
	  end if
	  ifoc = 1
	  gpa1 = gpa
	  do k = 1,3
	    oldmid1(k) = oldmid(k)
	  end do
	  goto 300
	end if
c	print *,'To get Donnan coefficient divide coion excess'
c	print *,'(Q+ or Q- depending on sign) by total charge in grid'
c
c rad. dist. fcn for ions
c
	open(29,file='ionintm_rad.dat')
	do i = 1,nbin
c
c integrated charge cations, anions, net charge
c
	  radion(5,i) = radion(5,i-1) + (2.*radion(1,i) + radion(2,i))
	  radion(6,i) = radion(6,i-1) + (2.*radion(3,i) + radion(4,i))
	  radion(7,i) = radion(5,i) - radion(6,i)
c
c convert actual concs. to cyl.  symm. equivalent to compare with radion.f analysis of MD
c
	  dvol = (i - 0.5)*twopi*(zup - zlw)
	  radcorr(i) = radcorr(i)/dvol
	  do k = 1,4
	    radion(k,i) = radion(k,i)/dvol/fact
	    if(radcorr(i).gt.0.)then
	      radion(k,i) = radion(k,i)/radcorr(i)
	    end if
	  end do
	  write(29,'(i6,8f10.3,i8,f10.5)')i,dvol,(radion(k,i),k=1,7),jbin(i),radcorr(i)
	end do
	close(29)

      end
c----------------------------------------------------------
	function factor(n)
	if((n.lt.0).or.(n.gt.20)) then
	  print *,' factorial argument out of range'
	else
	  factor = 1.
	  do i = 1,n
	    factor = factor*i
	  end do
	end if
	return
	end

c----------------------------------------------------------
	subroutine geteps(scale,oldmid,reps,rkeps)
c	kim sharp/9 aug 86
c	reads in compressed eps file
c	new format: neps(5,65,65,3) integer*2
c	where first index 1-65 now compressed
c	into 1-5 plus offset into 16 bit words
c	compact format also contains oldmid, the
c	center of the protein in real coordinates
c	thus obviating the scale file.
c	compaction is effected by storing
c	real eps (which take values of 0. and 1.)
c	as bits in a 16 bit word
c	access is via pointers idimx and ioffset
c	thus x arrary indices of reps 0-15 -> word 1
c	16-31 -> word 2 etc
c	
      parameter (isize = 65 ,isize1=isize/16+1 )    !grid size
	
	integer*2 neps,keps,ioffset
	dimension neps(isize1,isize,isize,3),keps(isize1,isize,isize)
      dimension reps(isize,isize,isize,3)    !fine eps map
      dimension rkeps(isize,isize,isize)    !fine debye map
      dimension idimx(isize)		!array of pointers to words
      dimension ioffset(isize)	!array of pointers to bit offsets
	dimension oldmid(3)
	character*80 epsname
c	integer*2 for$ib_set

101   print *, 'name of epsilon file without .eps?'
      read(5,'(a)')epsname
	il=lnblnk(epsname)
      open(unit=10,file=epsname(:il)//'.eps',form='unformatted',
     &status = 'old',err=101)
      print *,'reading compact epsilon file...'
      read(10)kmap, scale, oldmid
      read(10)neps
	if(kmap.eq.1)then
	  print *,'reading debye map too...'
	  read(10)keps
	end if
      close(10)

      print *,' setting up pointers...'
	do ix = 1, isize
	idimx(ix) = ix/16 + 1
	ioffset(ix) = mod(ix,16)
	end do
c-------------------------------------
      print * , ' generating real fine epsilon array...'
	nout = 0
	nin = 0
	do idir = 1, 3
      do iz=1,isize
        do iy=1,isize
          do ix=1,isize
	if(btest(neps(idimx(ix),iy,iz,idir),ioffset(ix)))then 
           reps(ix,iy,iz,idir) = 1.
	     nin = nin + 1
        else
           reps(ix,iy,iz,idir) = 0.
	     nout = nout + 1
        end if
	enddo
	enddo
	enddo
	enddo
	print *,'inside, outside dielectric points: ',nin,nout
      print * , ' generating debye array...'
	nout = 0
	nin = 0
      do iz=1,isize
        do iy=1,isize
          do ix=1,isize
	if(btest(keps(idimx(ix),iy,iz),ioffset(ix)))then 
           rkeps(ix,iy,iz) = 1.
	     nout = nout + 1
        else
           rkeps(ix,iy,iz) = 0.
	     nin = nin + 1
        end if
	enddo
	enddo
	enddo
	print *,'inside, outside debye points: ',nin,nout
      return
      end
