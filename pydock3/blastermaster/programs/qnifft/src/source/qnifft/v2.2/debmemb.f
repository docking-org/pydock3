	subroutine debmemb(zin,zout)
	include 'qdiffpar.h'
c
c set ionic strength in inner region (zin< z < zout)
c and outer region (z > zout) to rionsti, rionsto (debfcti,debfcto)
c respectively
c
	goff = (igrid+1.)/2.
	do iz = 1,igrid
	  zz = (iz - goff)/scale + oldmid(3)
c	  print *,'left: ',iz
	  if(zz.gt.zin)then
	    if(zz.gt.zout)then
c	      print *,'right: ',iz
	      do iy = 1,igrid
	        do ix = 1,igrid
		    if(debmap(ix,iy,iz).ne.0)then
			debmap(ix,iy,iz) = 3 ! kas 8-feb-01
		    end if
	        end do
	      end do
	    else
c	      print *,'middle: ',iz
	      do iy = 1,igrid
	        do ix = 1,igrid
		    if(debmap(ix,iy,iz).ne.0)then
			debmap(ix,iy,iz) = 2 ! kas 8-feb-01
		    end if
	        end do
	      end do
	    end if
	  end if
	end do
	return
	end
