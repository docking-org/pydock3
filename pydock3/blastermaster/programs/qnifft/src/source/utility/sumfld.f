	program sumfld
c------------------------------------------------------------
c kim sharp nov 95
c finds sum(q.phi) for each residue in *.fld file
c------------------------------------------------------------
	dimension crg(10000),phi(10000)
	character*60 file1,arg2
	character*80 skip
	character*132 line
	character atnam*6,rsrec*9,rsrec1*9,atrec*15
	integer natom
	logical dosum
c------------------------------------------------------------
c
c greeting
c
	print *,' '
	print *,'finds sum(q.phi) over residues in *.fld file output from qnifft'
	print *,'or average phi for each residue'
	print *,' '
	marg = iargc()
	if(marg.lt.1)then
	  print *,'usage: sumfld field_file {sum|average}'
	  stop
	end if
	call getarg(1,file1)
	dosum = .true.
	if(marg.gt.1)then
	  call getarg(2,arg2)
	  if(arg2(1:3).eq.'ave')dosum = .false.
	  print *,'averaging phi, not summing q.phi over residues...'
	end if
c
200	format(a)
	open(unit=10,file=file1,status='old',err=900)
c
c skip 12 lines of header
c
	do i = 1,12
	  read(10,200,end=900)skip
	  print *,skip
	end do
	nrec = 0
	nres = 0
	qqrs = 0.
	qqtot = 0.
	qphirs = 0
	qphitot = 0.
	natom = 0
	rsrec1 = 'dummy    '
	print *,' '
	if(dosum)then
	  print *,' Summary of net charge, sum of Q.phi (kT) by residue: '
	  print *,'   #  Residue      Qnet     Q.Phi '
	else
	  print *,' Summary of net charge, mean Phi     (kT) by residue: '
	  print *,'   #  Residue      Qnet     <Phi> '
	end if
c
c read next line of file, checking for EOF
c
110	continue
	read(10,201,err=102,end=102)x1,y1,z1,q1,p1,fx,fy,fz,atrec
201	format(8f11.4,1x,a15)
	nrec = nrec + 1
	rsrec = atrec(7:15)
	if(rsrec.ne.rsrec1)then
	  if(nres.gt.0)then
	    if((.not.dosum).and.(natom.gt.0))then
		qphirs = qphirs/natom
	    end if
	    write(6,'(i6,1x,a9,2f9.4)')nres,rsrec1,qqrs,qphirs
	  end if
	  nres = nres + 1
	  rsrec1 = rsrec
	  qqrs = 0.
	  qphirs = 0.
	  natom = 0
	end if
	qqrs = qqrs + q1
	qqtot = qqtot + q1
	natom = natom + 1
	if(dosum)then
	  qphirs = qphirs + q1*p1
	  qphitot = qphitot + q1*p1
	else
	  qphirs = qphirs + p1
	  qphitot = qphitot + p1
	end if
	goto 110
102	continue
	close(10)
	if(nres.gt.0)then
	  if((.not.dosum).and.(natom.gt.0))then
		qphirs = qphirs/natom
	  end if
	  write(6,'(i6,1x,a9,2f9.4)')nres,rsrec1,qqrs,qphirs
	end if
	print *,'# records read       : ',nrec
	if((.not.dosum).and.(nrec.gt.0))then
	  qphitot = qphitot/nres
	  print *,'Net charge, Av.phi   : ',qqtot,qphitot
	else
	print *,'Net charge, Q.phi    : ',qqtot,qphitot
	end if
	stop
900	continue
	print *,'cant find file : ',file1
	stop
	end
