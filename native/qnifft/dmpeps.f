	subroutine dmpeps(filnam)
c
c dump out eps, deb maps in compact formatted for for easy comparison
c
c------------------------------------------------------------------
	include 'qdiffpar.h'
	character*80 filnam
c------------------------------------------------------------------
	write(6,*)'   '
	write(6,*)'dumping dielectric and debye maps...'
	write(6,*)'   '
	open(file=filnam,unit=41)
	write(41,*)'epsmap ix,iy,ix,idir,new value'
	ilast=-1
	do i=1,3
	do iz=1,igrid
	  do iy=1,igrid
	    do ix=1,igrid
	      if(iepsmp(ix,iy,iz,i).ne.ilast)then
	        ilast = iepsmp(ix,iy,iz,i)
	        write(41,'(5I8)')ix,iy,iz,i,ilast
	      end if
	    end do
	  end do
	end do
	end do
	write(41,*)'debmap ix,iy,ix,new value'
	ilast=-1
	do iz=1,igrid
	  do iy=1,igrid
	    do ix=1,igrid
	      if(debmap(ix,iy,iz).ne.ilast)then
	        ilast = debmap(ix,iy,iz)
	        write(41,'(5I8)')ix,iy,iz,ilast
	      end if
	    end do
	  end do
	end do
	close(41)
	return
	end 
c
