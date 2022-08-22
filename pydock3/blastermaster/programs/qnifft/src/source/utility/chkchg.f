	program chkchg
c check net charge on residues from delphi output
	character line*80,rnum*4,head*6,atm*5,res*3,chn*1
	character rnuml*4,atml*5,resl*3,chnl*1
	character chgstr*7,warning*30
c--------------------------------------------
	chnl = '*'
	resl = '*  '
	atml = '*    '
	rnuml = '*   '
	warning = '   Warning: non-unit charge!'
	tnetq = 0.
	rnetq = 0.
	nrchg = 0
	nptchg = 0.
	print *,'res   #    chn   atoms   net charge  non-integral charge'
103   continue
	natom = 0
      read(5,'(a)',end=303)line
      head = line(1:6)
	call up(head,6)
c
c skip header lines
c
      if((head.ne.'ATOM  ').and.(head.ne.'HETATM')) goto 103
	natom = natom + 1
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
	chgstr = line(61:67)
	read(chgstr,'(f7.3)',err=103)atq
c	print *,line
c	print *,atq
	if((res.ne.resl).or.(chn.ne.chnl).or.(rnum.ne.rnuml))then
c new residue
c	  if((chnl.ne.'*').and.(rnetq.ne.0.))then
	  if(chnl.ne.'*')then
	    inetq = nint(rnetq)
	    dq = (rnetq-inetq)**2
	    if(dq.gt.0.005)then
	    write(6,'(a3,3x,a4,3x,a1,3x,i3,3x,f8.3,a30)')
     &resl,rnuml,chnl,nrchg,rnetq,warning
	    nptchg = nptchg + 1
	    else
	    write(6,'(a3,3x,a4,3x,a1,3x,i3,3x,f8.3)')
     &resl,rnuml,chnl,nrchg,rnetq
	    endif
	  end if
	  tnetq = tnetq + rnetq
	  rnetq = atq
	  nrchg = 1
	  resl = res
	  chnl = chn
	  rnuml = rnum
	else
	  rnetq = rnetq + atq
	  nrchg = nrchg + 1
	end if
c
      go to 103
303   continue		
c	end of file
	inetq = nint(rnetq)
	dq = (rnetq-inetq)**2
	if(dq.gt.0.005)then
	    write(6,'(a3,3x,a4,3x,a1,3x,i3,3x,f8.3,a30)')
     &resl,rnuml,chnl,nrchg,rnetq,warning
	    nptchg = nptchg + 1
	else
	    write(6,'(a3,3x,a4,3x,a1,3x,i3,3x,f8.3)')
     &resl,rnuml,chnl,nrchg,rnetq
	endif
	tnetq = tnetq + rnetq
	print *,'total charge: ',tnetq
	print *,'number of residues with non-integral charge: ',nptchg
	end 
c------------------------------------------------------------
      subroutine up(txt,len)
c
c convert character string to upper case
c
      character*(*) txt
      character*80 save
      character*26 ualpha,lalpha
      data ualpha /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data lalpha /'abcdefghijklmnopqrstuvwxyz'/

      do 9000 i=1,len
        if((txt(i:i).ge.'a').and.(txt(i:i).le.'z')) then
          match = index(lalpha,txt(i:i))
          save(i:i) = ualpha(match:match)
        else
          save(i:i) = txt(i:i)
        end if
c      end do
9000	continue

      txt = save
      return
      end

c----------------------------------------------------------
      subroutine elb(txt,len)
c
c eliminate leading blanks from a character string
c
      character*(*) txt
      character*80 save

      do 9000 i=1,len
        if(txt(i:i).ne.' ') then
          nonb = i
          go to 100
        end if
9000  continue
      return
100   continue
      save = txt(nonb:len)
      txt = save
      return
      end
