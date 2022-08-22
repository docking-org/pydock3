	program qnifft22
c------------------------------------------------------------------
c author: kim sharp
c cleaned up feb 01 so compiles on linux/lahey fortran and compaq/f77
c reduce memory requirements feb 01
c cleaned up mkmaps subroutines feb 14, renamed mkepsd
c option to read radii,charges from atm file, apr 01
c------------------------------------------------------------------
	include 'qdiffpar.h'
c------------------------------------------------------------------
c------------------------------------------------------------------
      character*60 toplbl
      character*10 nxtlbl
	logical iper,iconc,isite,isitecrg,iatout,ibios,isph,icen,ianal
	logical iatin
	logical fasteps,force
	logical donon,dmmy2
	character*80 radfil,crgfil,pdbfili,pdbfilo,phifili,phifilo
	character*80 epsfil,sitefil,frcfil,sitecrgfil,parfil
        character*80 pdb2fili
c pdb2fili is the probe points, added by ryan coleman 2012
        logical run_extra_pdb !process the second pdb file passed in additionally

c------------------------------------------------------------------
c------------------------------------------------------------------
	dimension iper(3)
	dimension xo(3)
	dimension xn(3)
	dimension offset(3)
	dimension cqplus(3),cqmin(3)
c-----------------------------------------------------------------
	real*8 atfrct(3,natmx),ergt
c------------------------------------------------------------------
c
c conversion to kT/charge at 25 celsius
c
      data epkt,fpi /561.0,12.566/
c
c converts between ionic st (M) and debye length
c
	data dfact / 3.047/	
c---------------------------------------------------
c	data zion / 2.0, 1.0, -2.0, -1.0 / ! kas 1-feb-01
c
c------------------------------------------------------------------
c------------------------------------------------------------------
	write(6,*)'   '
	write(6,*)'------------------QNIFFT22--------------------------'
	write(6,*)'* Version 2.2                                            *'
	write(6,*)'* a program to solve the non-linear P-B equation         *'
	write(6,*)'* 2 dielectric regions, ionic strength, periodic         *'
	write(6,*)'* and focussing boundary conditions.                     *'
	write(6,*)'* Does nonlinear salt analysis                           *'
	write(6,*)'* uses mike holsts inexact newton method                 *'
	write(6,*)'* implemented by jonathan hecht, new input style         *'
	write(6,*)'* and multigridding by kas 6 sept, 1993                  *'
	write(6,*)'* dielectric smoothing added july 95, kas                *'
	write(6,*)'* 2nd, site potl. charge file added. Reads atom          *'
	write(6,*)'* names of form 2H2                                      *'
	write(6,*)'* charge anti-aliasing added may 96                      *'
	write(6,*)'* Mixed salt NLPB added Oct 99 by Vinod Misra/KAS        *'
	write(6,*)'* Mixed salt NLPB doesnt work with channel ionic diffrnce*'
	write(6,*)'* reduce memory requirements- 7-feb-01                   *'
	write(6,*)'* cleaned up mkmaps 14-feb-01                            *'
	write(6,*)'* read radii, charges from *.atm file 10-apr-01          *'
	write(6,*)'* compute atomic forces               10-apr-01          *'
	write(6,*)'* output reaction field, potl.        27-apr-01          *'
	write(6,*)'* improved surface mapping            18-may-01          *'
	write(6,*)'* include reentrant surface remapping 28-may-01          *'
        write(6,*)'* added parameter for addtl pdb files rgc 7-aug-12       *'
	write(6,*)'------------------QNIFFT22--------------------------'
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
	call getarg(1,parfil)
	call getpar(perfil,oldmid,scale,border,isize,epsout,epsin,radprb,
     &  exrad,rionst,rionsti,rionsto,zin,zout,conv,omega,rtemp,igrid,nitmx,nit0,
     &  nit1,nit2,nit3,iupdate,ichk0,iper,ibctyp,
     &  donon,iconc,ibios,isite,isitecrg,iatout,isph,icen,ianal,
     &  toplbl,radfil,crgfil,pdbfili,pdbfilo,
     &  phifili,phifilo,epsfil,sitefil,frcfil,sitecrgfil,ismth,smthlen,fasteps,force,
     &  cion,bfield,iatin,parfil,pdb2fili)
      print *,'cion: ',cion
c kas 1-feb-01
	zion(1) = 2.0
	zion(2) = 1.0
	zion(3) = -2.0
	zion(4) = -1.0
	do i = 1,4 ! kas 8-feb-01
	  czion(i) = cion(i)*zion(i)
	  czion2(i) = cion(i)*zion(i)**2
	end do ! kas 8-feb-01
c kas 1-feb-01
      print *,'zion: ',zion
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
	epstab(0) = epsout ! kas 8-feb-01
	epstab(1) = epsin ! kas 8-feb-01 
	do i = 0,6
	  sepstab(i) = i*epsin + (6-i)*epsout
	end do
	
c
c read radius file and store in hash table
c
	if(.not.iatin)call getrfil(radfil)
c
c read charge parameter file and store in hash table
c
	if(.not.iatin)call getcfil(crgfil)
c
c atom coordinate file assign radii and charges
c
	dmmy2 = .false.
        if (pdb2fili(1:1) .eq. ' ') then !if empty, no filename, don't add 
          run_extra_pdb = .false. !probe points in the getpdb function
        else !otherwise
          run_extra_pdb = .true. !add additional probe points in 2nd file
        endif
	call getpdb(pdbfili,pdbfilo,iatout,dmmy2,isitecrg,cqplus,cqmin,qplus,
     &      qmin,natom,iatin,pdb2fili,run_extra_pdb)
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
	if(donon)then
	  debfct = fpi/(1636.3*scale**2)
	else
	  debfct = 2.*rionst*fpi/(1636.3*scale**2)
	end if
	debtab(0) = 0. ! kas 8-feb-01
	debtab(1) = debfct ! kas 8-feb-01
c
c make eps and debye maps
c
c	call mapit(exrad,radprb,scale,natom,atcrd,atrad,ngrid,igrid,iepsmp,iatmap,debmap,iatmapd)
	if(fasteps) then
	  print *,'Using fast dielectric mapping-no surface charge correction'
	  call mkepsd(natom,exrad,radprb)
	else
	  print *,'Using accurate dielectric mapping- Please wait!'
	  call mapgen(exrad,radprb,natom)
	end if
	call wrteps(epsfil) ! kas 8-may-01
c	stop ! debug
c	epsfil = 'mapit.eps'
c	call wrteps(epsfil) ! kas 8-may-01
	if((rionsti.ne.rionst).or.(rionsto.ne.rionst))then
	  print *,'setting ionic strength in membrane channel and outer solution'
	  print *,'ionic strength in membrane region (zin < z < zout)'
	  print *,'and outer region (z > zout):   ',rionsti,rionsto
	  print *,'boundaries in Z direction (A): ',zin,zout
	  debfcti = 2.*rionsti*fpi/(1636.3*scale**2)
	  debfcto = 2.*rionsto*fpi/(1636.3*scale**2)
	  debtab(2) = debfcti ! kas 8-feb-01
	  debtab(3) = debfcto ! kas 8-feb-01
	  call debmemb(zin,zout)
	else
	  debfcti = debfct
	  debfcto = debfct
	  debtab(2) = debfct ! kas 8-feb-01
	  debtab(3) = debfct ! kas 8-feb-01
	end if
	print *,'debye factors: ',debtab ! kas 8-feb-01
	print *,'dieletric factors: ',epstab ! kas 8-feb-01
	if(ismth.ne.0)then
	  print *,'SORRY- dielectric smoothing doesnt go with eps mapping!' ! kas 8-feb-01
c	  call smtheps(ismth,smthlen)
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
	call phierg(ianal,donon,rionst,epsin,epsout,natom,iconc,epkt,rtemp,fasteps)
	if(force)then
	call phifrc(ianal,donon,rionst,epsin,epsout,natom,iconc,epkt,rtemp,atfrct,ergt,
     & fasteps)
      end if
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
234     format('     x         y          z          (e)        (kT/e)       Ex        Ey       Ez    Atom Res C   #')

c               123456789 123456789 123456789 123456789 123456789 123456789 123456789
	  call getpdb(sitefil,frcfil,iatout,isite,isitecrg,cqplus,cqmin,qplus,qmin,natom,iatin,sitefil,.false.)
	  close(16)
	  epsin = epsin/epkT
	  epsout = epsout/epkT
c	  call phirxn(ianal,epsin,epsout,natom,epkt,rtemp,frcfil,fasteps)
      end if 
c
c charge site option end
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
	  if(phifilo(1:9).ne."/dev/null")then
          nxtlbl = 'potential '
	    open(unit=14,file=phifilo,form='unformatted')
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
	end if
c
c write dielectric map
c
c	call wrteps(epsfil) ! kas 8-feb-01

	goto 999
999	continue
	finish = cputime(start1)
	write(6,*)'  '
	write(6,*)'total cpu time was (sec) ',finish
	write(6,*)'  '
	write(6,*)'QNIFFT QNIFF QNIF. QNI.. QN... Q....'
	end
