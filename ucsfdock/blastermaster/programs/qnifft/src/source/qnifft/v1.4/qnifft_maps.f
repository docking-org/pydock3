	program qnifft
c------------------------------------------------------------------
c author: kim sharp
c------------------------------------------------------------------
	include 'qdiffpar.h'
c------------------------------------------------------------------
      character*8 hour
      character*9 day
      character*60 toplbl,comment
      character*10 nxtlbl
	character*9 bclab(5)
	character*80 line,filnam
	character*6 head
	character*3 tres
	character*6 tatm
	character*24 crdstr
	character*13 radstr
	logical iper,iconc,isite,iatout,ibios,isph,icen,ianal
	logical iauto,donon,dmmy2
	character*80 radfil,crgfil,pdbfili,pdbfilo,phifili,phifilo
	character*80 epsfil,sitefil,frcfil
	real*4 realv(30)
	integer intv(30)
	character*80 titles(5)

c------------------------------------------------------------------
c------------------------------------------------------------------
	dimension iper(3)
	dimension cmin(3)
	dimension cmax(3)
	dimension cran(3),pmid(3)
	dimension xo(3)
	dimension xn(3)
	dimension ismin(3)
	dimension ismax(3)
	dimension xyzmid(3)
	dimension offset(3)
	dimension cqplus(3),cqmin(3)
	dimension tomid(3,3)   
c	offset of midlines from grid points
c-----------------------------------------------------------------
c------------------------------------------------------------------
c
c conversion to kT/charge at 25 celsius
c
      data epkt,fpi /561.0,12.566/
c
c converts between ionic st (M) and debye lenght
c
	data dfact / 3.047/	
c
c defines offset to midlines along three axes
c
      data tomid /0.5,0.0,0.0, 0.0,0.5,0.0, 0.0,0.0,0.5/
c
	data bclab / 'zero     ','dipolar  ',
     &             'focussing','coulombic',
     &             'ext.field'/
c
c------------------------------------------------------------------
c------------------------------------------------------------------
	write(6,*)'   '
	write(6,*)'------------------QNIFFT----------------------------'
	write(6,*)'*                                                  *'
	write(6,*)'* a program to solve the non-linear P-B equation   *'
	write(6,*)'* 2 dielectric regions, ionic strength, periodic   *'
	write(6,*)'* and focussing boundary conditions. Includes      *'
	write(6,*)'* corrected epsmap & reaction field calculation.   *'
	write(6,*)'* Does nonlinear salt analysis                     *'
	write(6,*)'* uses mike holsts inexact newton method           *'
	write(6,*)'* implemented by jonathan hecht, new input style   *'
	write(6,*)'* and multigridding by kas 6 sept, 1993            *'
	write(6,*)'*variable dimension format for phi,eps maps, 19/6/94*'
	write(6,*)'------------------QNIFFT----------------------------'
	write(6,*)'  '
c------------------------------------------------------------------
c
c print title, description
c
      call date(day)
      call time(hour)
	write(6,*)'  '
      write(6,*)' program run on ',day
      write(6,*)'             at ',hour
	write(6,*)'  '
      print *,'program dimensioned at: ',ngrid
c------------------------------------------------------------------
c
c initialize rad, charge link lists
c
	start1 = cputime(0.0)
	do 9000 i = 1,nrlist
	  irnumb(i) = 0
	  irlink(i) = 0
9000	continue
	do 9001 i = 1,nclist
	  icnumb(i) = 0
	  iclink(i) = 0
9001	continue
c
c read in run parameters
c
	call getpar(perfil,oldmid,scale,border,isize,epsout,epsin,radprb,
     &  exrad,rionst,conv,omega,rtemp,igrid,nitmx,nit0,
     &  nit1,nit2,nit3,iupdate,ichk0,iper,ibctyp,
     &  donon,iconc,ibios,isite,iatout,isph,icen,ianal,
     &  toplbl,radfil,crgfil,pdbfili,pdbfilo,
     &  phifili,phifilo,epsfil,sitefil,frcfil)
      if(igrid.gt.ngrid)then
        print *,'igrid exceeds maximum dimension: ',ngrid
        stop
	endif
c 
c convert ionic strength to debye length, taking account of temperature
c
	dfact = 3.047*sqrt(rtemp*epsout/80./298.)
	if(rionst.ne.0.) then
	  deblen = dfact/sqrt(rionst)
	  print *,'Debye length (A): ',deblen
	else
	  deblen = 1.e6
	  print *,'Debye length infinite'
	end if
c
c scale dielectric so that potentials turn out as kT when 
c charges are in e's and distances are in angstroms
c
	epkt = epkt*298./rtemp
	epsin = epsin/epkT
	epsout = epsout/epkT
	
c
c read radius file and store in hash table
c
	open(unit=11,file=radfil,status='old',err=901)
	print *,'reading radii from file'
	print *,radfil
c
c skip/print comments (start with !) and one column header line 
c
105	read(11,201,end=901)comment
	if(comment(1:1).eq.'!') then
	  write(6,*)comment
	  goto 105
	end if
201   format(a)
	nrdrec = 0
100   continue
	  nrdrec = nrdrec + 1
	  if(nrdrec.gt.nrmax) then
	    write(6,*)' maximum # of radius records exceeded'
	    write(6,*)' increase nrmax'
	    stop 
	  end if
	  read(11,200,err=904,end=300)atm,res,rad
	  call up(atm,6)
	  call elb(atm,6)
	  call up(res,3)
	  call elb(res,3)
	  atnam(nrdrec) = atm
	  rnam(nrdrec)  = res
	  radt(nrdrec)  = rad
	  call rent(atm,res,nrdrec)
	goto 100
300   continue
	close(11)
200   format(A6,A3,F8.3)
	nrdrec = nrdrec - 1
	write(6,*)'# of radius parameter records:',nrdrec
c
c read charge parameter file and store in hash table
c
	open(unit=12,file=crgfil,status='old',err=902)
	print *,'reading charges from file'
	print *,crgfil
c
c skip/print comments (start with !) and one column header line 
c
106	read(12,201,end=901)comment
	if(comment(1:1).eq.'!') then
	  write(6,*)comment
	  goto 106
	end if
	nchrec = 0
101   continue
	  nchrec = nchrec + 1
	  if(nchrec.gt.ncmax) then
	    write(6,*)' maximum # of charge records exceeded'
	    write(6,*)' - increase ncmax'
	    stop 
	  end if
	  read(12,202,err=905,end=301)atm,res,rnum,chn,chrgv
	  call up(atm,6)
	  call elb(atm,6)
	  call up(res,3)
	  call elb(res,3)
	  call up(rnum,4)
	  call elb(rnum,4)
	  call up(chn,1)
	  call elb(chn,1)
	  catnam(nchrec) = atm
	  crnam(nchrec)  = res
	  crnum(nchrec)  = rnum
	  cchn(nchrec)   = chn
	  chrgvt(nchrec) = chrgv
	  call cent(atm,res,rnum,chn,nchrec)
	goto 101
301   continue
	close(12)
202   format(A6,A3,A4,A1,F8.3)
	nchrec = nchrec - 1 
	write(6,*)'# of charge parameter records:',nchrec
c
c atom coordinate file assign radii and charges
c
	write(6,*)'   '
	write(6,*)'assigning radii and charges to atoms...'
	write(6,*)'   '
	natom = 0.
	qnet = 0.
	nqass = 0
	qplus = 0.
	qmin = 0.
	do k = 1,3
	  cqplus(k) = 0.
	  cqmin(k) = 0.
	end do
	if(iatout) then
	  open(unit=19,file=pdbfilo)
	  filnam = ' '
	  write(6,*)'atomic coordinates, charges and radii written to'
	  print *,pdbfilo
	  write(6,*)'   '
	  write(19,*)'HEADER output from qmiff'
	  write(19,*)'HEADER atom radii in columns 55-60'
	  write(19,*)'HEADER atom charges in columns 61-67'
	end if
204	format(a)
205	format(3f8.3)
	open(13,file=pdbfili,err=903,status='old')
	print *,'reading coordinates from file'
	print *,pdbfili
	ncorr = 0
103   continue
        read(13,204,end=303)line
        head = line(1:6)
	  call up(head,6)
c
c skip header lines
c
        if((head.ne.'ATOM  ').and.(head.ne.'HETATM')) then
	    if((iatout).and.(head.eq.'HEADER')) then
		write(19,204)line
	    endif
	    go to 103
	  end if
	  natom = natom + 1
	  if(natom.gt.natmx)then
	    print *,'Max # of atoms exceeded: increase natmx: ',natmx
	    stop
	  end if
	  crdstr = line(31:54)
        read(crdstr,205)atcrd(1,natom),atcrd(2,natom),atcrd(3,natom)
	  atm = line(12:16)
	  res = line(18:20)
	  rnum = line(23:26)
	  chn = line(22:22)
	  call up(atm,6)
	  call elb(atm,6)
c
c kludge to get around those programs that generate atom names of form
c 1hb1 etc, like insight
c
c	  print *,atm(1:1)
	  if((atm(1:1).eq.'1').or.(atm(1:1).eq.'2').or.(atm(1:1).eq.'3'))then
	    ncorr = ncorr + 1
	    tatm = atm(2:6)
	    do i = 1,6
		if(tatm(i:i).eq.' ')then
		  tatm(i:i) = atm(1:1)
c		  print *,'correcting atom name: ',atm,tatm
		  atm = tatm
		  goto 400
		end if
	    end do
400	    continue
	  end if
	  call up(res,3)
	  call elb(res,3)
	  call up(rnum,4)
	  call elb(rnum,4)
	  call up(chn,1)
	  call elb(chn,1)
c
c assign radius, searching for decreasingly specific specification
c ending with generic atom type
c note all atoms must have an assignment
c
	  call rfind(atm,res,ifind,n)
	  if(ifind.eq.0) then
	    tres = '   '
	    call rfind(atm,tres,ifind,n)
	    if(ifind.eq.0) then
		tatm = atm(1:1)//'     '
	      call rfind(tatm,tres,ifind,n)
	      if(ifind.eq.0) then
	        write(6,*)'no radius record found for'
	        write(6,*)line
	        stop
		end if
	    end if
	  end if
	  atrad(natom) = radt(n)
c
c
c assign charge to atom, searching for decreasingly specific specification
c kas,29jn 89- add extra searches with wild card:
C search order:  atom, res, num, chain
c		      x     x    x   x
c			x     x    x   
c			x     x        x
c			x     x
c			x          x   x
c			x          x
c			x              x
c			x
c
c note if no charge record found, is assumed to be 0.0
c
	  call cfind(atm,res,rnum,chn,ifind,n)
	  if(ifind.eq.0) then
	    schn = chn
	    chn = ' '
	    call cfind(atm,res,rnum,chn,ifind,n)
	    if(ifind.eq.0) then
		chn = schn
		snum = rnum
		rnum = '    '
	      call cfind(atm,res,rnum,chn,ifind,n)
	      if(ifind.eq.0) then
	        schn = chn
	        chn = ' '
	        call cfind(atm,res,rnum,chn,ifind,n)
		  if(ifind.eq.0) then
		    chn = schn
		    rnum = snum
		    sres = res
		    res = '   '
	          call cfind(atm,res,rnum,chn,ifind,n)
	          if(ifind.eq.0) then
	            schn = chn
	            chn = ' '
	            call cfind(atm,res,rnum,chn,ifind,n)
	            if(ifind.eq.0) then
		        chn = schn
		        snum = rnum
		        rnum = '    '
	              call cfind(atm,res,rnum,chn,ifind,n)
	              if(ifind.eq.0) then
	                schn = chn
	                chn = ' '
	                call cfind(atm,res,rnum,chn,ifind,n)
			  endif
			endif
		    endif
		  end if
		end if
	    end if
	  end if
c
c net charge etc
c
	  if(ifind.ne.0) then
	    chrgv = chrgvt(n)
	    nqass = nqass + 1
	  else
	    chrgv = 0.
	  end if
	  qnet = qnet + chrgv
	  atcrg(natom) = chrgv
	  if(chrgv.gt.0.0) then
	    qplus = qplus + chrgv
	    do k = 1,3
	      cqplus(k) = cqplus(k) + atcrd(k,natom)*chrgv
	    end do
	  end if
	  if(chrgv.lt.0.0) then
	    qmin = qmin + chrgv
	    do k = 1,3
	      cqmin(k) = cqmin(k) + atcrd(k,natom)*chrgv
	    end do
	  end if
c
c write record to new coordinate file if required, with
c occupancy and temperature factor fields replaced by radius and charge
c
	  if(iatout) then
          write(radstr,206)atrad(natom),atcrg(natom)
	    line(55:67) = radstr
206       format(F6.2,F7.3)
	    if(atrad(natom).ge.0.)write(19,204)line
	  end if

      go to 103
303   continue		
c	end of file
	close(13)
	if(iatout) then
	  close(19)
	end if
	write(6,*)'   '
	write(6,*)'number of atom coordinates read  : ',natom
	if(ncorr.gt.0)then
	  write(6,*)'number of atom names of type 1xxx corrected: ',ncorr
	end if
c
c calculate scale
c
	call sclit(atcrd,atrad,natom,igrid,perfil,scale,oldmid,border,isize,icen)
	midg = (igrid+1)/2
c
c scale atoms to grid coordinates
c
	do i = 1,natom
	  do k = 1,3
	    xo(k) = atcrd(k,i)
	  end do
	  call ctog(xo,xn)
	  do k = 1,3
	    atcrd(k,i) = xn(k)
	  end do
	end do
	if(qplus.ne.0.0) then
	  do k = 1,3
	    cqplus(k) = cqplus(k)/qplus
        end do
	end if
	if(qmin.ne.0.0) then
	  do k = 1,3
	    cqmin(k) = cqmin(k)/qmin
        end do
	end if
	do k = 1,3
	  cqplus(k) = (cqplus(k)-oldmid(k))*scale+midg
	  cqmin(k) = (cqmin(k)-oldmid(k))*scale+midg
      end do
	write(6,*)'total number of charged atoms    : ',nqass
	write(6,*)'net assigned charge              : ',qnet
	write(6,*)'assigned positive charge         : ',qplus
	write(6,*)'centred at (gu) :',cqplus
	write(6,*)'assigned negative charge         : ',qmin
	write(6,*)'centred at (gu) :',cqmin
	write(6,*)'   '
c
c scale 4.pi factor to grid
c
	fpoh = fpi*scale
c
c assign charges
c
	print *,'assigning charges to grid...'
	do i = 1,natom
	  if(atcrg(i).ne.0.)then
	    do k = 1,3
		xo(k) = atcrd(k,i)
	    end do
	    gchrgv = atcrg(i)*fpoh
	    if(isph) then
		call chrgit1(xo,gchrgv)
	    else
	      call chrgit(xo,gchrgv)
	    end if
	  end if
	end do
c
c calculate boundary conditions
c
	write(6,*)'  '
	write(6,*)'setting boundary conditions...'
	write(6,*)'  '
	call setbc(ibctyp,qplus,qmin,cqplus,cqmin,epsout,deblen,natom,phifili)
	debfct = epsout*rionst/(dfact*scale)**2
c
c make eps and debye maps

	call mkmaps(natom,exrad,radprb,debfct,epsin,epsout)
c	call coneps(natom,exrad,radprb,debfct,epsin,epsout)
c
	finish = cputime(start1)
	write(6,*)'  '
	write(6,*)'setup time was (sec) ',finish
	write(6,*)'  '
c
c iterate-
c
	call newt(nit0,nitmx,nit1,nit2,nit3,iper,conv,omega,ichk0,donon)
c write site potentials
c note if concentration option is selected, concentrations
c and conc. gradients are written
c
	nflag = -1
	if(donon)nflag = 1
	if(ianal)call phierg(nflag,rionst,epsin,epsout,natom,iconc,epkt,rtemp)
	epsin = epsin*epkT
	epsout = epsout*epkT
      if(isite) then
        open(unit=15,file=sitefil,status='old',err=903)
	  print *,'reading site potential coordinates from file'
	  print *,sitefil
        open(unit=16,file=frcfil)
	  print *,'writing potentials,fields ,at charge sites to file'
	  print *,frcfil
        write(16,*)'output from QMIFF   '
        write(16,*)'grid size,percent fill:',igrid,perfil
        write(16,*)'inner,outer dielectric:',epsin,epsout
        write(16,*)'ionic strength (M):',rionst
        write(16,*)'ion excl., probe radius:',exrad,radprb
        write(16,*)'level 0, newton iterations:',nit0,nitmx
        write(16,*)'boundary condition:',ibctyp
        write(16,*)'    '
        write(16,*)'title: '
        write(16,*)toplbl
        write(16,233)
        write(16,234)
233     format('         coordinates            charge
     $   potential         field (kT/e/Ang.)')
234     format('     x         y         z       (e)
     &      (kT/e)      Ex       Ey        Ez')
c      123456789 123456789 123456789 123456789 123456789 123456789 123456789
	  etot = 0.
c
c read atom coordinate file
c
	  natom = 0
104     continue
        read(15,204,end=304)line
        head = line(1:6)
	  call up(head,6)
c
c skip header lines
c
        if((head.ne.'ATOM  ').and.(head.ne.'HETATM')) then
	    go to 104
	  end if
	  natom = natom + 1
	  crdstr = line(31:54)
        read(crdstr,205)xo
	  atm = line(12:16)
	  res = line(18:20)
	  rnum = line(23:26)
	  chn = line(22:22)
	  call up(atm,6)
	  call elb(atm,6)
	  call up(res,3)
	  call elb(res,3)
	  call up(rnum,4)
	  call elb(rnum,4)
	  call up(chn,1)
	  call elb(chn,1)
c
c scale atoms to grid space
c
	  call ctog(xo,xn)
c
c assign charge to atom, searching for decreasingly specific specification
c note if no charge record found, is assumed to be 0.0
c
	  call cfind(atm,res,rnum,chn,ifind,n)
	  if(ifind.eq.0) then
	    schn = chn
	    chn = ' '
	    call cfind(atm,res,rnum,chn,ifind,n)
	    if(ifind.eq.0) then
		chn = schn
		snum = rnum
		rnum = '    '
	      call cfind(atm,res,rnum,chn,ifind,n)
	      if(ifind.eq.0) then
	        schn = chn
	        chn = ' '
	        call cfind(atm,res,rnum,chn,ifind,n)
		  if(ifind.eq.0) then
		    chn = schn
		    rnum = snum
		    sres = res
		    res = '   '
	          call cfind(atm,res,rnum,chn,ifind,n)
	          if(ifind.eq.0) then
	            schn = chn
	            chn = ' '
	            call cfind(atm,res,rnum,chn,ifind,n)
	            if(ifind.eq.0) then
		        chn = schn
		        snum = rnum
		        rnum = '    '
	              call cfind(atm,res,rnum,chn,ifind,n)
	              if(ifind.eq.0) then
	                schn = chn
	                chn = ' '
	                call cfind(atm,res,rnum,chn,ifind,n)
			  endif
			endif
		    endif
		  end if
		end if
	    end if
	  end if
c
	  chrgv = 0.0
	  if(ifind.ne.0) then
	    chrgv = chrgvt(n)
	  else
	    chrgv = 0.0
	  end if
c
c calculate potential, field and energy, and output to file
c note field is -grad(phi): prior to jan 94 grad(phi) was output!
c
        call phintp(xn,phiv)
	  xn(1) = xn(1) + 1.
        call phintp(xn,fxu)
	  xn(1) = xn(1) - 2.
        call phintp(xn,fxl)
	  xn(1) = xn(1) + 1.
	  xn(2) = xn(2) + 1.
        call phintp(xn,fyu)
	  xn(2) = xn(2) - 2.
        call phintp(xn,fyl)
	  xn(2) = xn(2) + 1.
	  xn(3) = xn(3) + 1.
        call phintp(xn,fzu)
	  xn(3) = xn(3) - 2.
        call phintp(xn,fzl)
        qphiv = chrgv*phiv
	  etot = etot + qphiv
	  fx = -(fxu-fxl)*scale/2.
	  fy = -(fyu-fyl)*scale/2.
	  fz = -(fzu-fzl)*scale/2.
	  write(16,230)xo,chrgv,phiv,fx,fy,fz
230     format(8G10.3)
        go to 104
304     continue		
c	end of file
	  close(15)
	  write(6,*)'   '
	  write(6,*)'number of atom coordinates read  : ',natom
	  write(6,*)'   '
	  etot = etot/2.
	  write(16,*)'total energy = ',etot,' kt'
	  close(16)
      end if 
c
	if(ibios) then
c
c write phimap in insight format
c
	  open(unit=14,file=phifilo,form='unformatted')
	  write(6,*)'potential map written in INSIGHT format to file'
	  write(6,*)phifilo
	  write(6,*)'  '
	  ivary = 0
	  nbyte = 4
	  intdat = 0
	  xang = 90.
	  yang = 90.
	  zang = 90.
	  intx = igrid - 1
	  inty = igrid - 1
	  intz = igrid - 1
	  xmax = 0.
	  do 9040 k = 1,3
	    temp = abs(oldmid(k))
	    xmax = amax1(xmax,temp)
9040	  continue
	  range = (igrid-1.)/(2.*scale)
	  extent = range + xmax
	  xstart = (oldmid(1)-range)/extent
	  ystart = (oldmid(2)-range)/extent
	  zstart = (oldmid(3)-range)/extent
	  xend = (oldmid(1)+range)/extent
	  yend = (oldmid(2)+range)/extent
	  zend = (oldmid(3)+range)/extent
        write(14)toplbl
	  write(14)ivary,nbyte,intdat,extent,extent,extent,xang,yang,
     &  zang,xstart,xend,ystart,yend,zstart,zend,intx,inty,intz
	  do 9041 k = 1,igrid
	    do 9042 j = 1,igrid
		write(14)(phimap(i,j,k),i=1,igrid)
9042	  continue
9041	  continue
c
	else
c
c write in DELPHI format
c
	  nunit=14
	  realv(1)=scale
	  realv(2)=oldmid(1)
	  realv(3)=oldmid(2)
	  realv(4)=oldmid(3)
	  realv(5)=epsin
	  realv(6)=epsout
	  realv(7)=rionst
	  do k = 1,5
	   titles(k) = ' '
	  end do
	  titles(1)=toplbl
	  if(donon)then
	    intv(1) = 1
	  else
	    intv(1) = 0
	  end if
	  open(unit=nunit,file=phifilo,form='unformatted')
	  call putphi(nunit,ngrid,igrid,phimap,realv,intv,titles)
        close(nunit)
	endif
c
c write dielectric map
c
	do  k = 1,igrid
       do j = 1,igrid
         do i = 1,igrid
            epsmap(i,j,k,1) = iepsmp(i,j,k,1)
            epsmap(i,j,k,2) = iepsmp(i,j,k,2)
            epsmap(i,j,k,3) = iepsmp(i,j,k,3)
            debmap(i,j,k) = debmap(i,j,k)
	    end do
	  end do
	end do
      write(6,*)' writing to compact epsilon file'
	nunit=17
	open(unit=nunit,file=epsfil,form='unformatted')
	call puteps(nunit,ngrid,igrid,epsmap,debmap,realv,intv,titles)
	close(nunit)

	goto 999
901	write(6,*) 'unexpected end or non-existence of radius file'
   	stop
902	write(6,*) 'unexpected end or non-existence of charge file'
   	stop
903	write(6,*) 'unexpected end or non-existence of atom file'
   	stop
904	write(6,*) 'error in reading radius file'
   	stop
905	write(6,*) 'error in reading charge file'
   	stop
999	continue
	finish = cputime(start1)
	write(6,*)'  '
	write(6,*)'total cpu time was (sec) ',finish
	write(6,*)'  '
	write(6,*)'QNIFFT QNIFF QNIF. QNI.. QN... Q....'
	end
