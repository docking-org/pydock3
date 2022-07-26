	program qcall
c
c try calling subroutine version of qnifft21
c
	character line*80,parfil*80,head*6
	integer natmx
	parameter (natmx = 32000)
	real*4 atcrd1(3,natmx),atrad1(natmx),atcrg1(natmx)
	real*8 atfrct(3,natmx),ergt
c--------------------------------------------
	call getarg(1,parfil)
	natom = 0
	qnet = 0
103   continue
      read(5,'(a)',end=303)line
      head = line(1:6)
	call up(head,6)
c	print *,'head: ',head ! debug
c
c skip header lines
c
      if((head.ne.'ATOM  ').and.(head.ne.'HETATM')) goto 103
	natom = natom + 1
	if(natom.gt.natmx)then
	  print *,'WARNING max atoms exceeded increase NATMX ',natmx
	  stop
	end if
	read(line(31:54),'(3f8.3)')(atcrd1(k,natom),k=1,3)
	read(line(55:67),'(f6.2,f7.3)')atrad1(natom),atcrg1(natom)
	qnet = qnet + atcrg1(natom)
	goto 103
303	continue
	print *,'atoms read: ',natom
	print *,'net charge: ',qnet
	print *,'calling qnifft.....'
	call qnifft22(natom,atcrd1,atrad1,atcrg1,atfrct,ergt,parfil)
	print *,'total solvation energy: ',ergt
	end 
c------------------------------------------------------------
