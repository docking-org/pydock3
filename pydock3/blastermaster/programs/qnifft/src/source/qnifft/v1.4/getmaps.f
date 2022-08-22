c
c all routines that read or write unformatted potential and
c epsmaps for delphi/qdiff/qniff etc
c altered to allow for non-65 dimension grids
c
c arrays realv and intv provided to store various parameters
c
c in realv so far: scale, xmid,ymid,ymid,epsin,epsout,rionst
c in intv so far: nonlinear flag
c
	subroutine putphi(nunit,ngrid,igrid,phimap,realv,intv,titles)
	real*4 phimap(ngrid,ngrid,ngrid),realv(30)
	integer intv(30)
	character*80 titles(5)

	print *,'writing phimap in new format...'
	write(nunit)ngrid,igrid
	write(nunit)realv
	write(nunit)intv
	write(nunit)titles
	do i = 1,igrid
	  do j = 1,igrid
	    write(nunit)(phimap(k,j,i),k=1,igrid)
	  end do
	end do
	return
	end 
	
	subroutine getphi(nunit,ngrid,igrid,phimap,realv,intv,titles)
	real*4 phimap(ngrid,ngrid,ngrid),realv(30)
	integer intv(30)
	character*80 titles(5)

	print *,'reading phimap in new format...'
	read(nunit)ngrid1,igrid
	if(igrid.gt.ngrid)then
	  print *,'WARNING: reading map of larger dimension, ',igrid
	  print *,'than array passed to getphi: ',ngrid
	  print *,'increase ngrid'
	  stop
	end if
	read(nunit)realv
	read(nunit)intv
	read(nunit)titles
	do k = 1,5
	  print *,titles(k)
	end do
	do i = 1,igrid
	  do j = 1,igrid
	    read(nunit)(phimap(k,j,i),k=1,igrid)
	  end do
	end do
	return
	end 

	subroutine puteps(nunit,ngrid,igrid,epsmap,debmap,realv,intv,titles)
	real*4 epsmap(ngrid,ngrid,ngrid,3),realv(30)
	real*4 debmap(ngrid,ngrid,ngrid)
	parameter (ngrid1=257,nbyte=ngrid1/16+1)
	integer intv(30)
	integer*2 btrw(nbyte)
	integer*2 idimx(ngrid1),ioffset(ngrid1)
	character*80 titles(5)

	print *,'writing epsmap in new format...'
	if(igrid.gt.ngrid1)then
	  print *,'grid size, ',igrid
	  print *,'too large for packed bit array in puteps: ',ngrid1
	  print *,'increase ngrid1'
	  return
	end if
	write(nunit)ngrid,igrid
	write(nunit)realv
	write(nunit)intv
	write(nunit)titles
	do i = 1,igrid
	  idimx(i) = i/16 + 1
	  ioffset(i) = mod(i,16)
	end do
	nbyte1 = igrid/16+1
	nin = 0
	nout = 0
	print *,'generating compact epsmap...'
	do idir = 1,3
	  do i = 1,igrid
	    do j = 1,igrid
		do k = 1,nbyte1
		  btrw(k) = 0
		end do
		do k = 1,igrid
		  if(epsmap(k,j,i,idir).eq.1.) then
c set bit
		    nin = nin + 1
		    btrw(idimx(k)) = iibset(btrw(idimx(k)), ioffset(k))
		  else
		    nout = nout + 1
              end if
		end do
	      write(nunit)(btrw(k),k=1,nbyte1)
	    end do
	  end do
	end do
	print *,'inside,outside points: ',nin,nout
	print *,'generating compact debmap...'
	nin = 0
	nout = 0
	  do i = 1,igrid
	    do j = 1,igrid
		do k = 1,nbyte1
		  btrw(k) = 0
		end do
		do k = 1,igrid
		  if(debmap(k,j,i).ne.0.) then
c set bit
		    nout = nout + 1
		    btrw(idimx(k)) = iibset(btrw(idimx(k)),ioffset(k))
		  else
		    nin = nin + 1
              end if
		end do
	      write(nunit)(btrw(k),k=1,nbyte1)
	    end do
	  end do
	print *,'inside,outside points: ',nin,nout
	return
	end 

	subroutine geteps(nunit,ngrid,igrid,epsmap,debmap,realv,intv,titles)
	real*4 epsmap(ngrid,ngrid,ngrid,3),realv(30)
	real*4 debmap(ngrid,ngrid,ngrid)
	parameter (ngrid1=257,nbyte=ngrid1/16+1)
	integer intv(30)
	integer*2 btrw(nbyte)
	integer*2 idimx(ngrid1),ioffset(ngrid1)
	character*80 titles(5)

	print *,'reading epsmap in new format...'
	read(nunit)ngrid2,igrid
	if(igrid.gt.ngrid)then
	  print *,'grid size, ',igrid
	  print *,'too large for array passed to geteps: ',ngrid
	  print *,'increase ngrid'
	  stop
	end if
	if(igrid.gt.ngrid1)then
	  print *,'grid size, ',igrid
	  print *,'too large for packed array: ',ngrid1
	  print *,'increase ngrid1'
	  stop
	end if
	read(nunit)realv
	read(nunit)intv
	read(nunit)titles
	do k = 1,5
	  print *,titles(k)
	end do
	do i = 1,igrid
	  idimx(i) = i/16 + 1
	  ioffset(i) = mod(i,16)
	end do
	nbyte1 = igrid/16+1
	nin = 0
	nout = 0
	print *,'regenerating compact epsmap...'
	do idir = 1,3
	  do i = 1,igrid
	    do j = 1,igrid
	      read(nunit)(btrw(k),k=1,nbyte1)
		do k = 1,igrid
		  if(bitest(btrw(idimx(k)),ioffset(k)))then
		    nin = nin + 1
		    epsmap(k,j,i,idir) = 1.
		  else
		    epsmap(k,j,i,idir) = 0.
		    nout = nout + 1
		  end if
		end do
	    end do
	  end do
	end do
	print *,'inside,outside points: ',nin,nout
	print *,'generating  debmap...'
	nin = 0
	nout = 0
	do i = 1,igrid
	  do j = 1,igrid
	    read(nunit)(btrw(k),k=1,nbyte1)
	    do k = 1,igrid
		 if(bitest(btrw(idimx(k)),ioffset(k)))then
		   debmap(k,j,i) = 1.
		   nout = nout + 1
		 else
		   debmap(k,j,i) = 0.
		   nin = nin + 1
		 end if
	    end do
	  end do
	end do
	print *,'inside,outside points: ',nin,nout
	return
	end 
