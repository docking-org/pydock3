	subroutine getrfil(radfil)
c------------------------------------------------------------------
	include 'qdiffpar.h'
c------------------------------------------------------------------
	character*80 radfil
	character*60 comment
c------------------------------------------------------------------
c
c initialize radius lists
c
	do i = 1,nrlist
	  irnumb(i) = 0
	  irlink(i) = 0
	end do
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

	return
901	write(6,*) 'unexpected end or non-existence of radius file'
   	stop
904	write(6,*) 'error in reading radius file'
   	stop
999	continue
	end

c----------------------------------------------------------------
      function irhash(atxt,rtxt)
c	
c produce hash number from atom and residue name
c for radius assignment
c	
	include 'qdiffpar.h'
c---------------------------------------------------
      character*6 atxt
      character*3 rtxt
      character*38 string
	integer n,irhash,j
      data string /'* 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      n = 1
      do 9000 i = 1,3
        j = index(string,rtxt(i:i))
        n = 5*n + j
c      end do
9000	continue
      do 9001 i = 1,6
        j = index(string,atxt(i:i))
        n = 5*n + j
c      end do
9001	continue
	n = abs(n)
      irhash = mod(n,nrlist) + 1
      return
      end

c------------------------------------------------------------------
      subroutine rent(atm,res,nent)
c
c enter character strings res,atm into hash table for radii
c by assigning it a number nent
c
	include 'qdiffpar.h'
c
c check to see if there is room
c
	if(irtot.eq.nrlist) then
	  write(6,*)' radii list full- increase nrlist'
	  stop 
	end if

	n = irhash(atm,res)
	if(irnumb(n).ne.0) then
c
c slot filled
c run down linked list
c
9000	continue
	if(irlink(n).eq.0)goto 9001
            n = irlink(n)
	goto 9000
9001	continue
c
c  search for empty slot
c
         new = 1
9002	continue
	if(irnumb(new).eq.0)goto 9003
            new = new + 1
	goto 9002
9003	continue
c
c found one- add link
c
         irlink(n) = new
         n = new
      end if

c
c slot empty
c fill slot
c
      irnumb(n) = nent
      irlink(n) = 0
	irtot = irtot + 1
      return
      end

c--------------------------------------------------------
      subroutine rfind(atm,res,ifind,n)
c
c	find entry res in radius hash table and check match 
c
	include 'qdiffpar.h'
c---------------------------------------------------
      n = irhash(atm,res)
c	check for match
      ifind = 0
100   continue
c	while no match and not at end of link list
      if(irnumb(n).eq.0) then !no match
        ifind = 0
        return
      end if
      if((res.eq.rnam(irnumb(n))).and.(atm.eq.atnam(irnumb(n)))) then
        n = irnumb(n)
        ifind = 1
        return
      else
        if(irlink(n).ne.0) then
          n = irlink(n)
        else
          ifind = 0
          return
        end if
      end if
      go to 100
      end
