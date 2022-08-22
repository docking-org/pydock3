      subroutine wrteps(epsfil)
c reformats epsilon array iepsmp to compact array 
c full format: iepsmp(65,65,65,3) integer*4 bit format: neps(5,65,65,3) integer*2
c where first index 1-65 now compressed into 1-5 plus offset into 16 bit words
c oldmid is the c center of the protein in real coordinates
c compaction is effected by storing real eps (which take values of 0. and 1.)
c as bits in a 16 bit word .access is via pointers idimx and ioffset
c thus x arrary indices of reps 0-15 -> word 1 ; 16-31 -> word 2 etc
c 22 may 89, also include debyemap
c	
	include 'qdiffpar.h'
	parameter (isize1 = ngrid/16+1)
c compact epsilon map
	dimension neps(isize1,ngrid,ngrid,3)	
c compact debye map
	dimension keps(isize1,ngrid,ngrid)	
c array of pointers to words
	dimension idimx(ngrid)		
c array of pointers to bit offsets
	dimension ioffset(ngrid)		
	integer*2  neps,ioffset,keps
c 	integer  neps,ioffset
	character*80 epsfil
c---------------------------------------------------------------------

	do 9000 ix = 1, ngrid
	  idimx(ix) = ix/16 + 1
	  ioffset(ix) = mod(ix,16)
9000	continue

      write(6,*)' generating compact epsilon array...'
	nin = 0
	nout = 0
      do 9005 idir = 1, 3
        do 9006 iz=1,ngrid
          do 9007 iy=1,ngrid
            do 9004 ix=1,isize1
              neps(ix,iy,iz,idir) = 0
9004		continue
            do 9008 ix=1,ngrid
              if(iepsmp(ix,iy,iz,idir).eq.1) then
c set bit
		   nin = nin + 1
                neps(idimx(ix),iy,iz,idir) = ibset(
     &          neps(idimx(ix),iy,iz,idir), ioffset(ix))
		  else
		   nout = nout + 1
              end if
9008		continue
9007	    continue
9006	  continue
9005	continue
	print *,'inside, outside dielectric points: ',nin,nout
      write(6,*)' generating compact debye array...'
	nin = 0
	nout = 0
        do 9016 iz=1,ngrid
          do 9017 iy=1,ngrid
            do 9014 ix=1,isize1
              keps(ix,iy,iz) = 0
9014		continue
            do 9018 ix=1,ngrid
              if(debmap(ix,iy,iz).ne.0) then
c set bit
		   nout = nout + 1
                keps(idimx(ix),iy,iz) = ibset(
     &          keps(idimx(ix),iy,iz),ioffset(ix))
		  else
		   nin = nin + 1
              end if
9018		continue
9017	    continue
9016	  continue
	print *,'inside, outside debye points: ',nin,nout

      write(6,*)' writing to compact epsilon file'
	open(unit=17,file=epsfil,form='unformatted')
	write(6,*)' '
	write(6,*)'dielectric map written to file'
	write(6,*)epsfil
	write(6,*)' '
	kmap = 1
	write (17) kmap,scale,oldmid
	write (17) neps
	write (17) keps
	close (17)
      return
	end
