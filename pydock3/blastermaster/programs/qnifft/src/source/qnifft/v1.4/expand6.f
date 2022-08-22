	subroutine expand(ibios)
c
c expand/contract igrid**3 grid to 65**3
c grid to give compatibility with previous
c phimap and epsmap formats using trilinear interpolation
c inext is the next (multigridding) or final grid size to expand to
c-------------------------------------------------------
	include 'qdiffpar.h'
c-------------------------------------------------------
	dimension gc(3)
c-------------------------------------------------------
	logical ibios
c-------------------------------------------------------
	rscale = (igrid-1.)/(65-1.)
	if(igrid.lt.65)then
	  ibeg = 65
	  ifin = 1
	  iint = -1
	else
	  ibeg = 1
	  ifin = 65
	  iint = 1
	end if
	print *,'ibios: ',ibios
	print *,'rscale: ',rscale
c
c do high end first to prevent overwriting
c find small grid values and interpolate into big grid
c only expand phimap if we're not writing in insight format
c
	  if(.not.ibios) then
	    do 9000 iz = ibeg,ifin,iint
	      gc(3) = (iz-1)*rscale + 1.
	      do 9001 iy = ibeg,ifin,iint
	        gc(2) = (iy-1)*rscale + 1.
	        do 9002 ix = ibeg,ifin,iint
	          gc(1) = (ix-1)*rscale + 1.
		    call phintp(gc,phiv)
		    phimap(ix,iy,iz) = phiv
9002		  continue
9001	      continue
9000	    continue
	  end if
c
c dielectric and debye map
c
	do 9003 iz = ibeg,ifin,iint
	  igz = (iz-1)*rscale + 1.
	  do 9004 iy = ibeg,ifin,iint
	    igy = (iy-1)*rscale + 1.
	    do 9005 ix = ibeg,ifin,iint
	      igx = (ix-1)*rscale + 1.
		epsmap(ix,iy,iz,1) = iepsmp(igx,igy,igz,1)
		epsmap(ix,iy,iz,2) = iepsmp(igx,igy,igz,2)
		epsmap(ix,iy,iz,3) = iepsmp(igx,igy,igz,3)
		debmap(ix,iy,iz) = debmap(igx,igy,igz)
9005	    continue
9004	  continue
9003	continue
c
c recalculate scale for expanded grids
c
	scale = scale/rscale
	write(6,*)'new scale is ',scale,' grids/ang'
	return
	end
