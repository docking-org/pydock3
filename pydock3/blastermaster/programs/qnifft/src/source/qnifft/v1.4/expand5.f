	subroutine expand(ibios)
c
c expands igrid**3 grid to ngrid**3
c grid (65**3 default) to give compatibilty with previous
c phimap and epsmap formats using trilinear interpolation
c inext is the next (multigridding) or final grid size to expand to
c-------------------------------------------------------
	include 'qdiffpar.h'
c-------------------------------------------------------
	dimension gc(3)
c-------------------------------------------------------
	logical ibios
c-------------------------------------------------------
	rscale = (igrid-1.)/(ngrid-1.)
c
c do high end first to prevent overwriting
c find small grid values and interpolate into big grid
c only expand phimap if we're not writing in insight format
c
	  if(.not.ibios.and.(igrid.lt.ngrid)) then
	    do 9000 iz = ngrid,1,-1
	      gc(3) = (iz-1)*rscale + 1.
	      do 9001 iy = ngrid,1,-1
	        gc(2) = (iy-1)*rscale + 1.
	        do 9002 ix = ngrid,1,-1
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
	do 9003 iz = ngrid,1,-1
	  igz = (iz-1)*rscale + 1.
	  do 9004 iy = ngrid,1,-1
	    igy = (iy-1)*rscale + 1.
	    do 9005 ix = ngrid,1,-1
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
