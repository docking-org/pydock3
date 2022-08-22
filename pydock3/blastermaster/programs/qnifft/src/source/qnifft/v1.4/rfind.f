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

