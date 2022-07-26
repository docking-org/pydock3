	subroutine getpar(perfil,oldmid,scale,border,isize,epsout,epsin,radprb,
     &  exrad,rionst,rionsti,rionsto,zin,zout,conv,omega,rtemp,igrid,nitmx,nit0,
     &  nit1,nit2,nit3,iupdate,ichk0,iper,ibctyp,
     &  donon,iconc,ibios,isite,isitecrg,iatout,isph,icen,ianal,
     &  toplbl,radfil,crgfil,pdbfili,pdbfilo,
     &  phifili,phifilo,epsfil,sitefil,frcfil,sitecrgfil,ismth,smthlen,
     &  cion)
c
c flexible input subroutine for delphi
c to add a new parameter: 
c 1: add its keyword to the names data specification in the real, integer
c    logical or string section depending on what type it is.
c    Make sure that the letters before the abbreviation mark *
c    in the name are unique
c 2: give it a variable name of the same type, that will be added to
c    the calling argument list in this subroutine, AND the call statement
c    in the main program that calls this subroutine.  This is the name it
c    will be known by in main program.
c 3: assign the value that was read in to the variable in the section
c    marked ***ASSOCIATE*** below.
c 4: put entries into the qhelp subroutine.
c
c add parameters for smoothing- kas, 26/7/95
c
	parameter (nname=25)
	parameter (klen=21)
c
c names of variables
c type of variable: 1 real,2 integer, 3 logical, 4 ascii
c
	character*21 names(nname,4)
c
c values and defaults
c
	real valrl(nname)
	integer valint(nname)
	logical vallog(nname),ltemp
	character*80 valasc(nname),parfil,deldir
c
	character*80 line,phrase,kword,vword
	character*50 vtemp
	character*21 vtemp1
	character*3  st1
	character*1 dtemp
	logical hasdef(3,nname,4)
c
c variables passed back to qdiff
c
	dimension oldmid(3)
	logical iper(3),donon,iconc,ibios,isite,isitecrg,iatout,isph,icen,ianal
	character*10 sizing,bctyp
	character*60 toplbl
	character*80 label,radfil,crgfil,pdbfili,pdbfilo,phifili,phifilo
	character*80 epsfil,sitefil,frcfil,sitecrgfil
	dimension cion(4)
c
c here is where we define the names of variables as seen in the parameter file
c anything beyond the abbreviation mark * is ignored, including the * itself
c in the match checking, so the letters before * better be unique!
c reals first
c
	data names(1,1)   /'fill*                '/
	data names(2,1)   /'Xcen*ter             '/
	data names(3,1)   /'Ycen*ter             '/
	data names(4,1)   /'Zcen*ter             '/
	data names(5,1)   /'scal*e               '/
	data names(6,1)   /'bord*er_solvent      '/
	data names(7,1)   /'solv*ent_dielectric  '/
	data names(8,1)   /'solu*te_dielectric   '/
	data names(9,1)   /'prob*e_radius        '/
	data names(10,1)  /'ioni*c_radius        '/
	data names(11,1)  /'salt*_concentration  '/
	data names(12,1)  /'conv*ergence         '/
	data names(13,1)  /'rela*xation parameter'/
	data names(14,1)  /'temp*erature         '/
	data names(15,1)  /'leng*th_smooth       '/
	data names(16,1)  /'isalt*_concentration '/
	data names(17,1)  /'osalt*_concentration '/
	data names(18,1)  /'imemb*rane_position  '/
	data names(19,1)  /'omemb*rane_position  '/
	data names(20,1)  /'cat_div*alent        '/
	data names(21,1)  /'cat_mon*ovalent      '/
	data names(22,1)  /'ani_div*alent        '/
	data names(23,1)  /'ani_mon*ovalent      '/
c integers
	data names(1,2)   /'grid*size            '/
	data names(2,2)   /'newt*on_iterations   '/
	data names(3,2)   /'level0*_multi_grid_it'/
	data names(4,2)   /'level1*_multi_grid_it'/
	data names(5,2)   /'level2*_multi_grid_it'/
	data names(6,2)   /'level3*_multi_grid_it'/
	data names(7,2)   /'chec*k_frequency     '/
	data names(8,2)   /'smoo*th_dielectric   '/
c logicals
	data names(1,3)   /'xper*iodic           '/
	data names(2,3)   /'yper*iodic           '/
	data names(3,3)   /'zper*iodic           '/
	data names(4,3)   /'nonl*inear_equation  '/
	data names(5,3)   /'conc*entration output'/
	data names(6,3)   /'insi*ght_format      '/
	data names(7,3)   /'site_pot*entials     '/
	data names(8,3)   /'atom*_file_output    '/
	data names(9,3)   /'sphe*rical_charge_dis'/
	data names(10,3)  /'anal*yse_map         '/
c asci
	data names(1,4)   /'sizi*ng              '/
	data names(2,4)   /'boun*dary_condition  '/
	data names(3,4)   /'titl*e               '/
	data names(4,4)   /'radi*us_file         '/
	data names(5,4)   /'char*ge_file         '/
	data names(6,4)   /'pdb_in*put_file      '/
	data names(7,4)   /'pdb_out*put_file     '/
	data names(8,4)   /'phi_in*put_file      '/
	data names(9,4)   /'phi_out*put_file     '/
	data names(10,4)  /'diel*ectric_map_file '/
	data names(11,4)  /'site_in*put_file     '/
	data names(12,4)  /'site_out*put_file    '/
	data names(13,4)  /'site_ch*arge_file    '/
c--------------------------------------------------------------
	call getarg(1,parfil)
	line = parfil
	call elb(line,80)
	call up(line,80)
	if(line(1:1).eq.' ')call usage()
	if(line(1:4).eq.'HELP')call help()
	if(line(1:3).eq.'ALL')call parlis()
	call parhlp(line)
c
c initalize arrays
c
	do itype = 1,4
	do ind = 1,nname
	  call up(names(ind,itype),klen)
	  call elb(names(ind,itype),klen)
	  vtemp1 = names(ind,itype)
	  i1 = ichar(names(ind,itype)(1:1))
	  if((i1.lt.65).or.(i1.gt.90))then
	    names(ind,itype)(1:1) = '$'
	  end if
	  hasdef(1,ind,itype) = .false.
	  hasdef(2,ind,itype) = .false.
	  hasdef(3,ind,itype) = .false.
	end do
	end do
c process parameter value files, starting with default file
c
	ifile = 1
401	continue
	nkey = 0
	nline = 0
	nperr = 0
	nkerr = 0
	nverr = 0
	nover = 0
	if(ifile.eq.1)then
	  call getenv("DELDIR",deldir)
	  il = lnblnk(deldir)
	  print *,'trying to read default file delphi.def in directory'
	  print *,deldir(1:il)
	  print *,' '
	  print *,'Default parameter file '
	  print *,'----------------------- '
	  print *,' '
	  open (10,file=deldir(1:il)//'/delphi.def',err=400)
	else
	  print *,'trying to read user parameter file ',parfil
	  print *,' '
	  print *,'User parameter file '
	  print *,'----------------------- '
	  print *,' '
	  open (10,file=parfil,err=400)
	end if

c
c read lines until end of file
c
	eof = 0
	do while(eof.eq.0)
	  read(10,'(a)',end=300)line
	  nline = nline + 1
c	  call up(line,80)
	  ifst = 1
c
c parse line til end of line, picking out phrases between delimiters
c which are the space, comma or end of line
c
	  do while(ifst.lt.80)
	    ifst = nnbl(line,80,ifst)
	    if((ifst.le.80).and.((line(ifst:ifst).eq."!").or.
     &    (line(ifst:ifst).eq."#")))ifst=81
	    ilst = ndlm(line,80,ifst)
	    len = ilst-ifst
	    if((ifst.le.80).and.(len.gt.0))then
c
c we have a phrase
c
	      phrase = line(ifst:ilst-1)
c		print *,'phrase: ',phrase(1:len)
c
c find keyword and value in phrase, which are separated by "="
c
	      neq = ifndeq(phrase,len)
c		print *,'neq: ',neq
		if(neq.ge.len)then
		  print *,'cant parse phrase: ',phrase(1:len)
		  nperr = nperr+1
		else
		  len1 = neq-1
		  kword = phrase(1:len1)
		  call up(kword,len1)
		  len2 = len - neq
		  vword = phrase(neq+1:len)
c		  print *,'keyword: ',kword(1:len1)
c		  print *,'value: ',vword(1:len2)
c
c search for keyword in list, upto unique length
c
	        ind = 0
		  itype = 1
		  ierr = 1
		  match = 0
		  do while((itype.le.4).and.(match.eq.0))
		    ind = ind + 1
		    if(ind.gt.nname)then
			ind = 1
			itype = itype+1
		    end if
		    if(itype.le.4)match = ifndmtch(kword,len1,names(ind,itype))
		  end do
		  if(match.eq.0)then
c
c nomatch: give message and skip to next phrase
c
		    print *,'unknown keyword is skipped: ',kword(1:len1)
		    nkerr = nkerr+1
		  else
c
c find match: read value as appropriate variable type
c if no error, set value, if already set, give warning
c
		  if((hasdef(1,ind,itype)).or.(hasdef(2,ind,itype)))then
	print *,'NOTE: overriding previous value for ',names(ind,itype)
			nover = nover+1
		  end if
		  hasdef(ifile,ind,itype) = .true.
		  if(itype.eq.1)then
	          read(vword,*,err=903)valrl(ind)
		  else if(itype.eq.2)then
	          read(vword,*,err=903)valint(ind)
		  else if(itype.eq.3)then
	          read(vword,*,err=903)vallog(ind)
		  else if(itype.eq.4)then
	          read(vword,'(a)',err=903)valasc(ind)
		  end if
		  ierr = 0
c
c read error: give message and skip
c
903		  continue
	        if(ierr.eq.1)then
		    print *,'ERROR reading parameter value: ',vword(1:40)
		    print *,'for key word ',kword(1:klen),' on line ',nline
		    print *,'keyword skipped'
		    nverr = nverr+1
		  else
		     nkey = nkey+1
		  end if
		end if ! valid keyword
		endif !parseable phrase
	    end if  !phrase
	    ifst = ilst+1
	  end do  !end of line
	  goto 301
300	  eof = 1
301	  continue
	end do !end of file
	close(10)
400	continue
	print *,' '
	if(nline.eq.0)then
	  print *,'WARNING: cannot find parameter file, or file is empty'
	else
	  print *,'lines read                          : ',nline
	endif
	print *,'valid keywords found                : ',nkey
	print *,'parsing errors                      : ',nperr
	print *,'invalid parameter keywords          : ',nkerr
	print *,'invalid parameter values            : ',nverr
	print *,'parameter values overidden          : ',nover
	print *,' '
	ifile = ifile+1
	if(ifile.le.2)goto 401
c
c    ***ASSOCIATE*** 
c here is where we associate parameters with variable names used in
c calling program
c
      perfil    = valrl(1)
      oldmid(1) = valrl(2)
      oldmid(2) = valrl(3)
      oldmid(3) = valrl(4)
      scale     = valrl(5)
      border    = valrl(6)
      epsout    = valrl(7)
      epsin     = valrl(8)
      radprb    = valrl(9)
      exrad     = valrl(10)
      rionst    = valrl(11)
      conv      = valrl(12)
      omega     = valrl(13)
	rtemp     = valrl(14)
	smthlen   = valrl(15)
	rionsti   = valrl(16)
	rionsto   = valrl(17)
	zin       = valrl(18)
	zout      = valrl(19)
	do k = 1,4
	  cion(k) = valrl(k+19)
	end do
c
      igrid     = valint(1)
      nitmx     = valint(2)
      nit0      = valint(3)
      nit1      = valint(4)
      nit2      = valint(5)
      nit3      = valint(6)
      ichk0     = valint(7)
	ismth     = valint(8)
c
      iper(1)   = vallog(1)
      iper(2)   = vallog(2)
      iper(3)   = vallog(3)
      donon     = vallog(4)
      iconc     = vallog(5)
      ibios     = vallog(6)
      isite     = vallog(7)
      iatout    = vallog(8)
      isph      = vallog(9)
	ianal     = vallog(10)
c
      sizing    = valasc(1)
      bctyp     = valasc(2)
      toplbl     = valasc(3)
      radfil    = valasc(4)
      crgfil    = valasc(5)
      pdbfili   = valasc(6)
      pdbfilo   = valasc(7)
      phifili   = valasc(8)
      phifilo   = valasc(9)
      epsfil    = valasc(10)
      sitefil   = valasc(11)
      frcfil    = valasc(12)
      sitecrgfil= valasc(13)	
	print *,'values transferred to variables'
c
c make parameters consistent
c
c 
c if sizing=fill, dont need border,scale,centers
c if sizing=border, dont need fill,scale,centers
c if sizing=scale, dont need fill, border- if center missing, default to
c molecule center
c order of precedence is scale, border,fill
c 
	call up(sizing,10)
	icen = .false.
	if(sizing(1:3).eq.'SCA')then
	  isize = 1
	  hasdef(3,1,1) = .true.
	  hasdef(3,6,1) = .true.
	  if((.not.hasdef(1,2,1)).and.(.not.hasdef(2,2,1)))then
	    hasdef(3,2,1) = .true.
	    hasdef(3,3,1) = .true.
	    hasdef(3,4,1) = .true.
	    icen =.true.
	  end if
	  if((.not.hasdef(1,3,1)).and.(.not.hasdef(2,3,1)))then
	    hasdef(3,2,1) = .true.
	    hasdef(3,3,1) = .true.
	    hasdef(3,4,1) = .true.
	    icen =.true.
	  end if
	  if((.not.hasdef(1,4,1)).and.(.not.hasdef(2,4,1)))then
	    hasdef(3,2,1) = .true.
	    hasdef(3,3,1) = .true.
	    hasdef(3,4,1) = .true.
	    icen =.true.
	  end if
	else if(sizing(1:3).eq.'BOR')then
	  isize = 2
	  hasdef(3,1,1) = .true.
	  hasdef(3,2,1) = .true.
	  hasdef(3,3,1) = .true.
	  hasdef(3,4,1) = .true.
	  hasdef(3,5,1) = .true.
	else if(sizing(1:3).eq.'FIL')then
	  isize = 3
	  hasdef(3,2,1) = .true.
	  hasdef(3,3,1) = .true.
	  hasdef(3,4,1) = .true.
	  hasdef(3,5,1) = .true.
	  hasdef(3,6,1) = .true.
	else
	  print *,'illegal scaling option: ',sizing
	  print *,'using explicit SCALE option'
	  sizing='SCALE'
	  isize = 1
	  hasdef(3,1,1) = .true.
	  hasdef(3,6,1) = .true.
	end if
c
	if(ismth.eq.0)then
c don't need smoothing length
	  hasdef(3,15,1)=.true.
	end if
c if not site potential option, dont need 2nd input pdb and output frc files
c if site potential option, and the 2nd charge file is missing, will 
c automatically use the first charge file
	  isitecrg=.false.
	  hasdef(3,13,4)=.true.
	if(.not.isite)then
	  hasdef(3,11,4)=.true.
	  hasdef(3,12,4)=.true.
	else if(hasdef(1,13,4).or.hasdef(2,13,4))then
	  isitecrg=.true.
	  hasdef(3,13,4)=.false.
	end if
c if not focussing, dont need input phi map
	call up(bctyp,10)
	if(bctyp(1:3).eq.'FOC')then
	  ibctyp = 3
	else if(bctyp(1:3).eq.'COU')then
	  ibctyp = 4
	  hasdef(3,8,4)=.true.
	else if(bctyp(1:3).eq.'ZER')then
	  ibctyp = 1
	  hasdef(3,8,4)=.true.
	else if(bctyp(1:3).eq.'FIE')then
	  ibctyp = 5
	  hasdef(3,8,4)=.true.
	else if(bctyp(1:3).eq.'DIP')then
	  ibctyp = 2
	  hasdef(3,8,4)=.true.
	else 
	  print *,'illegal boundary condition option: ',bctyp
	  print *,'using ZERO potential option'
	  ibctyp = 1
	  bctyp = 'ZERO'
	end if
c if not outputting atm map, dont need name
	if(.not.iatout)then
	  hasdef(3,7,4)=.true.
	end if
c if no ionic strength, dont need ion layer thickness, or nonlinear
c conc output false
	rneut = 2.*(cion(1)-cion(3))+(cion(2)-cion(4))
	if(abs(rneut).gt.1.e-6)then
	  print *,'WARNING: not a neutral salt combination: '
	  print *,cion
	  cion(4) = cion(4) + rneut
	  valrl(23) = cion(4) 
	  print *,'setting monovalent anion conc to ',cion(4)
	else
	  print *,'NOTE: Using explicit anion, cation concentrations'
	  print *,'and ignoring salt concentration variable       '
	end if
	rionst = 2.*(cion(1)+cion(3))+0.5*(cion(2)+cion(4))
	valrl(11) = rionst
	if(rionst.le.0.)then
	  rionst = 0.
	  donon=.false.
	  hasdef(3,4,3)=.true.
	  exrad = 0.
	  hasdef(3,10,1)=.true.
	  iconc=.false.
	  hasdef(3,5,3)=.true.
	  hasdef(3,16,1)=.true.
	  hasdef(3,17,1)=.true.
	  hasdef(3,18,1)=.true.
	  hasdef(3,19,1)=.true.
	end if
c
c if no definition for inside and outside ionic strengths
c then set to standard inic strength
c
	if((hasdef(1,16,1).eq..false.).and.(hasdef(2,16,1).eq..false.))then
	  rionsti = rionst
	  zin = 1.e6
	  hasdef(3,16,1)=.true.
	  hasdef(3,17,1)=.true.
	end if
	if((hasdef(1,17,1).eq..false.).and.(hasdef(2,17,1).eq..false.))then
	  rionsto = rionst
	  zout = 1.e6
	  hasdef(3,18,1)=.true.
	  hasdef(3,19,1)=.true.
	end if
c if uniform dielectric, dont need probe radius
	if(epsin.eq.epsout)then
	  radprb = 0.
	  hasdef(3,9,1)=.true.
	end if
c if no level 0 iterations, dont need 1,2,3 level iterations
c if no level 1 iterations, dont need 2,3 level iterations
c if no level 2 iterations, dont need 3 level iterations
	if(nit0.eq.0)then
	  hasdef(3,4,2)=.true.
	  nit1 = 0
	end if
	if(nit1.eq.0)then
	  hasdef(3,5,2)=.true.
	  nit2 = 0
	end if
	if(nit2.eq.0)then
	  hasdef(3,6,2)=.true.
	  nit3 = 0
	end if
c
c print out parameters
c
	nig = 0
	nmiss = 0
	nass = 0
	nblank = 0
	ndflt = 0
	print *,'Summary of parameter assignments'
	print *,'--------------------------------'
	print *,' '
	print *,'Name                  | Def |  Value'
	print *,'----------------------------------------------------------'
	do itype = 1,4
	  if(itype.eq.1)then
	  print *,'                             real parameters'
	  else if(itype.eq.2)then
	  print *,'                          integer parameters'
	  else if(itype.eq.3)then
	  print *,'                          logical parameters'
	  else if(itype.eq.4)then
	  print *,'                            ascii parameters'
	  end if
	  do ind = 1,nname
	    if(names(ind,itype)(1:1).eq.'$')then
	      nblank = nblank+1
	    else
	      if(itype.eq.1)then
	        write(vtemp1,*)valrl(ind)
	      else if(itype.eq.2)then
	        write(vtemp1,*)valint(ind)
	      else if(itype.eq.3)then
	        write(vtemp1,*)vallog(ind)
	      else if(itype.eq.4)then
	        vtemp1 = valasc(ind)
	      end if
	      if(hasdef(3,ind,itype))then
	        nig = nig+1
	        vtemp1 = 'ignored'
	        st1 = ' '
	      else if((.not.hasdef(1,ind,itype)).and.
     &        (.not.hasdef(2,ind,itype)))then
	        nmiss = nmiss + 1
	        vtemp1 = 'missing'
	        st1 = '   '
	      else if((hasdef(1,ind,itype)).and.
     &        (.not.hasdef(2,ind,itype)))then
	        st1 = 'yes '
	        ndflt = ndflt + 1
	      else if(hasdef(2,ind,itype))then
	        st1 = 'no  '
	        nass = nass + 1
	      endif
	      write(6,*)names(ind,itype),' | ',st1,' | ',vtemp1
	    end if
	  end do
	end do
	print *,' '
	print *,'# of assigned parameters            : ',nass
	print *,'# of default parameters             : ',ndflt
	print *,'# of missing parameters             : ',nmiss
	print *,'# of ignored parameters             : ',nig
c	print *,'# of blank parameters               : ',nblank
	print *,' '
	if(nmiss.ne.0)then
	  print *,'stopping run because of missing parameters'
	  stop
	end if
	return
	end

	function nnbl(line,len,ist)
	character*(*) line
	character*1 chr,spc
	data spc / " "/
c
c find next non blank character
c
	nnbl = ist
	chr = line(nnbl:nnbl)
	do while((nnbl.le.len).and.(chr.eq.spc))
	  nnbl = nnbl+1
	  chr = line(nnbl:nnbl)
	end do
	return
	end

	function ndlm(line,len,ist)
	character*(*) line
c
c find next delimiter (space,comma)
c
c
	ndlm = ist
	do while((ndlm.le.len).and.(line(ndlm:ndlm).ne." ")
     &.and.(line(ndlm:ndlm).ne.","))
	  ndlm = ndlm+1
	end do
	return
	end

	function ifndeq(line,len)
	character*(*) line
	ifndeq = 1
	do while((line(ifndeq:ifndeq).ne."=").and.(ifndeq.le.len))
	  ifndeq = ifndeq + 1
	end do
	return
	end
	function ifndmtch(kword,len1,lname)
	character*(*) kword
	character*21 lname
c
c find unique abbreviation marker (*)
c
	iabb = 1
	ifndmtch = 1
	do while((iabb.le.21).and.(lname(iabb:iabb).ne."*"))
	  iabb = iabb + 1
	end do
	i = 1
c
c check for match, bailing out as soon as we fail
c
	do while((i.lt.iabb).and.(i.le.len1))
	  if(kword(i:i).ne.lname(i:i))then
	    ifndmtch = 0
	    return
	  end if
	  i = i+1
	end do
	return
	end

