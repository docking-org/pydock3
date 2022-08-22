	subroutine getpdb(pdbfili,pdbfilo,iatout,isite,isitecrg,cqplus,cqmin,qplus,qmin,natom)
c------------------------------------------------------------------
	include 'qdiffpar.h'
c------------------------------------------------------------------
c------------------------------------------------------------------
	character*80 line
	character*6 head
	character*3 tres
	character*6 tatm
	character*24 crdstr
	character*13 radstr
	logical isite,isitecrg,iatout
	character*80 pdbfili,pdbfilo
c------------------------------------------------------------------
	dimension xo(3)
	dimension xn(3)
	dimension cqplus(3),cqmin(3)
c-----------------------------------------------------------------
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
	  write(6,*)'atomic coordinates, charges and radii written to'
	  print *,pdbfilo
	  write(6,*)'   '
	  write(19,*)'HEADER output from qnifft'
	  write(19,*)'HEADER atom radii in columns 55-60'
	  write(19,*)'HEADER atom charges in columns 61-67'
	end if
	if(isite)then
	    etot = 0.
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
	  atm = line(12:17)
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
	  if(.not.isite)then
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
	  end if
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
	  if(isite)then
	    do k = 1,3
	      xo(k) = atcrd(k,natom)
	    end do
	    call ctog(xo,xn)
c
c calculate potential, field and energy, and output to file
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
	    xn(3) = xn(3) + 1.
          qphiv = chrgv*phiv
	    etot = etot + qphiv
	    fx = -(fxu-fxl)*scale/2.
	    fy = -(fyu-fyl)*scale/2.
	    fz = -(fzu-fzl)*scale/2.
	    write(16,230)xo,chrgv,phiv,fx,fy,fz,line(12:26)
230       format(8f11.4,1x,a15)
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
      write(6,*)'total number of source charged atoms    : ',nqass
      write(6,*)'net assigned source charge              : ',qnet
      write(6,*)'assigned positive charge                : ',qplus
      write(6,*)'centred at (A)  :',cqplus
      write(6,*)'assigned negative charge                : ',qmin
      write(6,*)'centred at (A)  :',cqmin
      write(6,*)'   '
	if(isite)then
	  if(isitecrg) then
	    print *,'Interaction between source and target charges (kT)',etot
	    print *,'(Includes Coulombic+Solvation energy)'
	    write(16,*)'Interaction between source and target charges (kT)',etot
	    write(16,*)'(Includes Coulombic+Solvation energy)'
	  else
	    etot = etot/2.
	    print *,'Interaction of source charges with themselves (kT)',etot
	    print *,'(Includes Grid+Coulombic+Solvation energy)'
	    write(16,*)'Interaction of source charges with themselves (kT)',etot
	    write(16,*)'(Includes Grid+Coulombic+Solvation energy)'
	  end if
	end if
	return

903	write(6,*) 'unexpected end or non-existence of atom file'
   	stop
	end
