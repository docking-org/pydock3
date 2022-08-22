	program qnifft
c------------------------------------------------------------------
c author: kim sharp
c------------------------------------------------------------------
	include 'qdiffpar.h'
c------------------------------------------------------------------
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
	logical iper,iconc,isite,isitecrg,iatout,ibios,isph,icen,ianal
	logical iauto,donon,dmmy2
	character*80 radfil,crgfil,pdbfili,pdbfilo,phifili,phifilo
	character*80 epsfil,sitefil,frcfil,sitecrgfil

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
c converts between ionic st (M) and debye length
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
c---------------------------------------------------
	data zion / 2.0, 1.0, -2.0, -1.0 /
c
c------------------------------------------------------------------
c------------------------------------------------------------------
	write(6,*)'   '
	write(6,*)'------------------QNIFFT----------------------------'
	write(6,*)'* Version 1.4                                      *'
	write(6,*)'* a program to solve the non-linear P-B equation   *'
	write(6,*)'* 2 dielectric regions, ionic strength, periodic   *'
	write(6,*)'* and focussing boundary conditions.               *'
	write(6,*)'* Does nonlinear salt analysis                     *'
	write(6,*)'* uses mike holsts inexact newton method           *'
	write(6,*)'* implemented by jonathan hecht, new input style   *'
	write(6,*)'* and multigridding by kas 6 sept, 1993            *'
	write(6,*)'* dielectric smoothing added july 95, kas          *'
	write(6,*)'* 2nd, site potl. charge file added. Reads atom    *'
	write(6,*)'* names of form 2H2                                *'
	write(6,*)'* charge anti-aliasing added may 96                *'
	write(6,*)'* Mixed salt NLPB added Oct 99 by Vinod Misra/KAS       *'
	write(6,*)'* Mixed salt NLPB dosnt work with channel ionic strength*'
	write(6,*)'*  differences                                     *'
	write(6,*)'------------------QNIFFT----------------------------'
	write(6,*)'  '
c------------------------------------------------------------------
	print *,'max grid dimension: ',ngrid
	start1 = cputime(0.0)
c
c print title, description
c
	call dtstmp()
c
c read in run parameters
c
	call getpar(perfil,oldmid,scale,border,isize,epsout,epsin,radprb,
     &  exrad,rionst,rionsti,rionsto,zin,zout,conv,omega,rtemp,igrid,nitmx,nit0,
     &  nit1,nit2,nit3,iupdate,ichk0,iper,ibctyp,
     &  donon,iconc,ibios,isite,isitecrg,iatout,isph,icen,ianal,
     &  toplbl,radfil,crgfil,pdbfili,pdbfilo,
     &  phifili,phifilo,epsfil,sitefil,frcfil,sitecrgfil,ismth,smthlen,
     &  cion,bfield)
      print *,'cion: ',cion
      print *,'zion: ',zion
c
c debug:
c
	if(igrid.gt.ngrid)then
	  print *,'ERROR: grid size greater than max of ',ngrid
	  stop
	end if
	midg = (igrid+1)/2
	do k = 1,3
	  offset(k) = 0.
	end do
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
	call getrfil(radfil)
c
c read charge parameter file and store in hash table
c
	call getcfil(crgfil)
c
c atom coordinate file assign radii and charges
c
	dmmy2 = .false.
	call getpdb(pdbfili,pdbfilo,iatout,dmmy2,isitecrg,cqplus,cqmin,qplus,qmin,natom)
c
c calculate scale
c
	call sclit(atcrd,atrad,natom,igrid,perfil,scale,oldmid,border,isize,icen)
c
c scale 4.pi factor to grid
c
	fpoh = fpi*scale
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
	do k = 1,3
	  cqplus(k) = (cqplus(k)-oldmid(k))*scale+midg
	  cqmin(k) = (cqmin(k)-oldmid(k))*scale+midg
      end do
	qnet = qplus + qmin
	write(6,*)'total number of source charged atoms    : ',nqass
	write(6,*)'net assigned source charge              : ',qnet
	write(6,*)'assigned positive charge                : ',qplus
	write(6,*)'centred at (gu) :',cqplus
	write(6,*)'assigned negative charge                : ',qmin
	write(6,*)'centred at (gu) :',cqmin
	write(6,*)'   '
c
c clear qmap before assigning charges
c
      do i=1,ngrid
        do j=1,ngrid
          do k=1,ngrid
            qmap(i,j,k)=0.
          end do
        end do
      end do
c
c assign charges
c
	print *,'assigning charges to grid...'
      start2= cputime(0.0)
      nout=0
      ichrgsa=0
      ichrgit=0
	do i = 1,natom
	  if(atcrg(i).ne.0.)then
	    do k = 1,3
		xo(k) = atcrd(k,i)
	    end do
	    gchrgv = atcrg(i)*fpoh
          xrad=atrad(i)*scale
	    if(isph.and.xrad.ge.0.246) then
            call chrgal(xo,gchrgv,xrad,nout)
            ichrgsa=ichrgsa+1
	    else
	      call chrgit(xo,gchrgv)
            ichrgit=ichrgit+1
	    end if
	  end if
	end do
      finish2= cputime(start2)
c
c print out qmap and report how atoms were charged and cpu time
c
      print *, '# of atoms charged by trilinear interpolation:',ichrgit
      print *, '# of atoms charged by anti-aliasing:',ichrgsa
      print *, 'time to charge grid:',finish2
c
c print alert if atomic radii exceed grid
c
      if (nout.ne.0) then
              print *,' '
              print *, 'WARNING: some atom radii exceed grid'
              print *, 'number of charged atoms out of grid:',nout
              print *, ' '
      end if
c
c calculate boundary conditions
c
	write(6,*)'  '
	write(6,*)'setting boundary conditions...'
	write(6,*)'  '
	call setbc(ibctyp,qplus,qmin,cqplus,cqmin,epsout,deblen,natom,phifili,bfield)
c	debfct = epsout*rionst/(dfact*scale)**2
	if(donon)then
	debfct = fpi/(1636.3*scale**2)
	else
	debfct = 2.*rionst*fpi/(1636.3*scale**2)
	end if
	print *,'debye factor: ',debfct
c
c make eps and debye maps
c
	call mkmaps(natom,exrad,radprb,debfct,epsin,epsout)
	if((rionsti.ne.rionst).or.(rionsto.ne.rionst))then
	  print *,'setting ionic strength in membrane channel and outer solution'
	  print *,'ionic strength in membrane region (zin < z < zout)'
	  print *,'and outer region (z > zout):   ',rionsti,rionsto
	  print *,'boundaries in Z direction (A): ',zin,zout
	  debfcti = 2.*rionsti*fpi/(1636.3*scale**2)
	  debfcto = 2.*rionsto*fpi/(1636.3*scale**2)
	  call debmemb(debfct,debfcti,debfcto,zin,zout)
	else
	  debfcti = debfct
	  debfcto = debfct
	end if
	if(ismth.ne.0)then
        print *,'calling dielectric smoothing routine'
	  call smtheps(ismth,smthlen)
	end if
c
	finish = cputime(start1)
	write(6,*)'  '
	write(6,*)'total setup time was (sec) ',finish
	write(6,*)'  '
c
c iterate-
c
	call newt(nit0,nitmx,nit1,nit2,nit3,iper,conv,omega,ichk0,donon)
c write site potentials
c note if concentration option is selected, concentrations
c and conc. gradients are written
c
	call phierg(ianal,donon,rionst,epsin,epsout,natom,iconc,epkt,rtemp)
c
c beginning site potential part
c
      if(isite) then
c
c if second charge file is assigned then read it, otherwise use
c data from the first charge file
c
	  qtnet = 0.
	  ntqass = 0
	  qtplus = 0.
	  qtmin = 0.
	  if(isitecrg)then
	    print *,' ' 
          print *,'reading second charge file'
          print *,sitecrgfil
	    call getcfil(sitecrgfil)
	  end if
c
c read site-in pdb file
c
	  iatout = .false.
	  print *,' '
	  print *,'reading site potential coordinates from file'
	  print *,sitefil
	  print *,'writing potentials,fields ,at charge sites to file'
	  print *,frcfil
	  open(unit=16,file=frcfil)
	  epsin = epsin*epkT
	  epsout = epsout*epkT
        write(16,*)'output from QNIFFT   '
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
233     format('         coordinates          charge   potential  field (kT/e/Ang.)')
234     format('     x         y       z       (e)     (kT/e)    Ex        Ey       Ez    Atom Res C   #')

c               123456789 123456789 123456789 123456789 123456789 123456789 123456789
	  call getpdb(sitefil,frcfil,iatout,isite,isitecrg,cqplus,cqmin,qplus,qmin,natom)
	  close(16)
      end if 
c
c charge site option end
c
	write(6,*)'  '
	write(6,*)'expanding potential and dielectric maps'
	write(6,*)'to grid of ',ngrid,'...'
	write(6,*)'  '
	call expand(ibios)
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
	  do k = 1,3
	    temp = abs(oldmid(k))
	    xmax = amax1(xmax,temp)
	  end do
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
	  do k = 1,igrid
	    do j = 1,igrid
		write(14)(phimap(i,j,k),i=1,igrid)
	    end do
	  end do
c
	else
c
c write in DELPHI format
c
        nxtlbl = 'potential '
	    open(unit=14,file=phifilo,form='unformatted')
	    filnam = ' '
	    write(6,*)'potential map written to file'
	    write(6,*)phifilo
	    write(6,*)'  '
          write(14)'now starting phimap '
          write(14)nxtlbl,toplbl
          write(14)phimap
          write(14)' end of phimap  '
          write(14)scale,oldmid
          close(14)
	end if
c
c write dielectric map
c
	call wrteps(epsfil,debfct,debfcti,debfcto)

	goto 999
999	continue
	finish = cputime(start1)
	write(6,*)'  '
	write(6,*)'total cpu time was (sec) ',finish
	write(6,*)'  '
	write(6,*)'QNIFFT QNIFF QNIF. QNI.. QN... Q....'
	end
