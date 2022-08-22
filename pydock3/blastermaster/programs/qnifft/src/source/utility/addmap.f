      program addmap
C	AUTHOR: KIM SHARP
C	REVISED 6 april 90
c     combine 2 potential maps
	parameter (ng = 65)
      dimension  phi1(ng,ng,ng)
      dimension  phi2(ng,ng,ng)
      dimension  phi3(ng,ng,ng)
	dimension oldmid1(3),oldmid2(3)
	data tiny /1.e-6/ !prevent divide by 0.
c------------------------------
100	continue
	print *,' '
	print *,' 1: add 2 maps'
	print *,' 2: map 1 - map 2 '
	print *,' 3: map 1 / map 2 '
	print *,' 4: multiply 2 maps'
	print *,' 5: add a constant to a map '
	print *,' 6: multiply a map by a constant'
	print *,' 9: quit'
	write(6,'(''option>>'',$)')
	read(5,200,err=100,end=900)iopt
200	format(i1)
	if(iopt.eq.9) goto 900
	if((iopt.lt.1).and.(iopt.gt.6)) goto 100
	print  *,'1st map'
	call getmap(phi1,scale1,oldmid1)
	if((iopt.ge.1).and.(iopt.le.4)) then
	  print  *,'2nd map'
	  call getmap(phi2,scale2,oldmid2)
	  call comp(oldmid1,oldmid2,scale1,scale2)
	end if
	if(iopt.eq.1) then
	  do i = 1,ng
	  do j = 1,ng
	  do k = 1,ng
	    phi3(k,j,i) = phi1(k,j,i) + phi2(k,j,i)
	  end do
	  end do
	  end do
	else if(iopt.eq.2) then
	  do i = 1,ng
	  do j = 1,ng
	  do k = 1,ng
	    phi3(k,j,i) = phi1(k,j,i) - phi2(k,j,i)
	  end do
	  end do
	  end do
	else if(iopt.eq.3) then
	  do i = 1,ng
	  do j = 1,ng
	  do k = 1,ng
	    phi3(k,j,i) = (phi1(k,j,i)+tiny)/(phi2(k,j,i)+tiny)
	  end do
	  end do
	  end do
	else if(iopt.eq.4) then
	  do i = 1,ng
	  do j = 1,ng
	  do k = 1,ng
	    phi3(k,j,i) = phi1(k,j,i) * phi2(k,j,i)
	  end do
	  end do
	  end do
	else if(iopt.eq.5) then
	  write(6,'(''additive constant>> '',$)')
	  read(5,*)const
	  do i = 1,ng
	  do j = 1,ng
	  do k = 1,ng
	    phi3(k,j,i) = phi1(k,j,i) + const
	  end do
	  end do
	  end do
	else if(iopt.eq.6) then
	  write(6,'(''multiplicative constant>> '',$)')
	  read(5,*)const
	  do i = 1,ng
	  do j = 1,ng
	  do k = 1,ng
	    phi3(k,j,i) = phi1(k,j,i) * const
	  end do
	  end do
	  end do
	end if
	call putmap(phi3,scale1,oldmid1)
	goto 100
900	continue
	print *,'bye'
	stop
	end
	subroutine comp(oldmid1,oldmid2,scale1,scale2)
	dimension oldmid1(3),oldmid2(3)
	if(scale1.ne.scale2) then
	 print *,' WARNING: map scales are not the same- will use'
	 print *,'first map values'
	else
	  do k = 1,3
	  if(oldmid1(k).ne.oldmid2(k)) then
	   print *,' WARNING: map centers are not the same- will use'
	   print *,'first map values'
	  endif
	  end do
	end if
	return
	end
c--------------------------------------------------------
	subroutine getmap(phi,scale,oldmid)
	parameter (ng = 65)
      dimension  phi(ng,ng,ng)
      character*80 phifil
      character*20 top_label
      character*16 bot_label
      character*80 up_label
	character*60 ntitle
	character*10 head
	dimension oldmid(3)

200	format(a)
100   print *,'full input file name'
      read(5,200)phifil
	il = lnblnk(phifil)
      open(10,file=phifil(:il),status='old',err=100,form='unformatted')
      read(10)top_label
      print *,top_label
      read(10)head,ntitle
      print *,head
      print *,ntitle
      read(10)phi
      read(10)bot_label
      print *,bot_label
      read(10)scale,oldmid
      print *,' scale, oldmid'
      print *,scale, oldmid
      close(10)
	return
      end
c--------------------------------------------------------
	subroutine putmap(phi,scale,oldmid)
	parameter (ng = 65)
      dimension  phi(ng,ng,ng)
      character*80 phifil
      character*20 top_label
      character*16 bot_label
      character*80 up_label
	character*60 ntitle
	character*10 head
	dimension oldmid(3)

      print *,'enter title for map'
      read(5,200)ntitle
	head=' >addmap  '

200	format(a)
      print *,'full output file name'
      read(5,200)phifil
	il = lnblnk(phifil)
      open(10,file=phifil(:il),status='unknown',form='unformatted')
      write(10)' now starting phimap'
      write(10)head,ntitle
      write(10)phi
      write(10) ' end of phi map '
      write(10)scale,oldmid
      close(10)
	return
      end
