        subroutine getcfil(crgfil)
c------------------------------------------------------------------
        include 'qdiffpar.h'
c------------------------------------------------------------------
        character*80 crgfil
	  character*60 comment
c------------------------------------------------------------------
c------------------------------------------------------------------
        do i = 1,nclist
          icnumb(i) = 0
          iclink(i) = 0
        end do
c
c read charge parameter file and store in hash table
c
        open(unit=12,file=crgfil,status='old',err=902)
        print *,'reading charges from file'
        print *,crgfil
c
c skip/print comments (start with !) and one column header line 
c
106     read(12,201,end=905)comment
        if(comment(1:1).eq.'!') then
          write(6,*)comment
          goto 106
        end if
        nchrec = 0
201	format(a)
101   continue
          nchrec = nchrec + 1
          if(nchrec.gt.ncmax) then
            write(6,*)' maximum # of charge records exceeded'
            write(6,*)' - increase ncmax:',ncmax
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
202   format(A6,A3,A4,A1,g8.3)
        nchrec = nchrec - 1 
        write(6,*)'# of charge parameter records:',nchrec
        return
902     write(6,*) 'unexpected end or non-existence of charge file'
        stop
905     write(6,*) 'error in reading charge file'
        stop
        end

c----------------------------------------------------------------
      function ichash(atxt,rtxt,ntxt,ctxt)
c	
c produce hash number from atom and residue name,rsidue number and chain
c name for charge assignment
c	
	include 'qdiffpar.h'

      character*6 atxt
      character*3 rtxt
      character*4 ntxt
      character*1 ctxt
      character*38 string
	integer n,ichash,j
      data string /'* 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      n = 1
      do 9000 i = 1,3
        j = index(string,rtxt(i:i))
c	  print *,'rtxt-> j,n ',j,n ! debug
        n = 5*n + j
c      end do
9000	continue
      do 9001 i = 1,6
        j = index(string,atxt(i:i))
c	  print *,'rtxt-> j,n ',j,n ! debug
c      end do
9001	continue
      do 9002 i = 1,4
        j = index(string,ntxt(i:i))
c	  print *,'rtxt-> j,n ',j,n ! debug
        n = 5*n + j
c      end do
9002	continue
      do 9003 i = 1,1
        j = index(string,ctxt(i:i))
c	  print *,'rtxt-> j,n ',j,n ! debug
        n = 5*n + j
c      end do
9003	continue
c	print *,'ichash n: ',n ! debug
	n = abs(n)
      ichash = mod(n,nclist) + 1
      return
      end

c------------------------------------------------------------------
      subroutine cent(atm,res,rnum,chn,nent)
c
c enter character strings res,atm,rnum,chn into hash table for charges
c by assigning it a number nent
c
	include 'qdiffpar.h'

c
c check to see if there is room
c
c	print *,'ictot,nclist: ',ictot,nclist
	if(ictot.eq.nclist) then
	  write(6,*)'charge list full- increase nclist'
	  stop
	end if

	n = ichash(atm,res,rnum,chn)
	if(icnumb(n).ne.0) then
c
c slot filled
c run down linked list
c
9000	   continue
	   if(iclink(n).eq.0)goto 9001
c         do while(iclink(n).ne.0)
            n = iclink(n)
c         end do            
	   goto 9000
9001	   continue
c
c  search for empty slot
c
         new = 1
9002	   continue
	   if(icnumb(new).eq.0)goto 9003
c         do while(icnumb(new).ne.0)
            new = new + 1
c         end do
    	   goto 9002
9003     continue

c
c found one- add link
c
         iclink(n) = new
         n = new
      end if

c
c slot empty
c fill slot
c
      icnumb(n) = nent
      iclink(n) = 0
	ictot = ictot + 1
      return
      end

c--------------------------------------------------------
      subroutine cfind(atm,res,rnum,chn,ifind,n)
c
c	find entry nres in hash table and check match with res
c
	include 'qdiffpar.h'


      n = ichash(atm,res,rnum,chn)
c	check for match
      ifind = 0
100   continue
c	while no match and not at end of link list
      if(icnumb(n).eq.0) then !no match
        ifind = 0
        return
      end if

      if((res.eq.crnam(icnumb(n))).and.(atm.eq.catnam(icnumb(n)))
     &  .and.(rnum.eq.crnum(icnumb(n))).and.(chn.eq.cchn(icnumb(n)))) 
     &  then
        n = icnumb(n)
        ifind = 1
        return
      else
        if(iclink(n).ne.0) then
          n = iclink(n)
        else
          ifind = 0
          return
        end if
      end if
      go to 100
      end

