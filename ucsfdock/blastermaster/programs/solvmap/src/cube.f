c------------------------------------------------------------------------
c   subroutine cubedata
c     - compute constants for distance cubing
c------------------------------------------------------------------------
	subroutine cubedata(frac, extrm, cblen)

	  include 'cubic.h'

	  real frac, cblen
	  real extrm(2, 3)
	  real off
	  real cmax(3)
          integer i

c         cubes per angstrom
	  cpa = 1./cblen
	  off = 0.1
          do i = 1,3
c           cube minima
            cmin(i) = extrm(1,i) - frac*cblen - off
c           cube maxima
            cmax(i) = extrm(2,i) + frac*cblen + off
c           number of cubes
            cb(i) = (cmax(i) - cmin(i))/cblen
          enddo

          return
	end

c------------------------------------------------------------------------
c   subroutine cube
c     - divide atoms into a cubic grid
c     - pointers to nearby atoms are stored in cbatom
c     - to get to the right section of cbatom, lower and upper bound pointers
c       into cbatom are stored in cblower and cbupper, respectively 
c------------------------------------------------------------------------
	subroutine cube(natm,crd,rad)

	  include 'cubic.h'

          integer natm
	  real crd(3,natm),rad(natm)
c   crd:  atomic x,y,z coordinates
c   rad:  atomic vdw radii
	  integer cblower(0:cb(1),0:cb(2),0:cb(3))
          integer cbupper(0:cb(1),0:cb(2),0:cb(3))
c   cblower: points to the lower bound of this box in cbatom
c   cbupper: points to the upper bound of this box in cbatom
          integer cbatom(1)
c   cbatom: points to all atoms contained in the nearby 27 cubes
	  integer icbn(3,natm)
c   icbn: integer coordinates of each atom inside the cubic grid
          integer icum
c   icum: cummulative counter
          integer i,j,k,ix,iy,iz,jx,jy,jz
          real x, y, z
          integer i_icbn
          integer memalloc

	  i_icbn = memalloc(i_icbn, 4*3*natm)
c         initialize pointers (so default do x=lower,upper does nothing)
	  do k=0,cb(3)
	    do j=0,cb(2)
	      do i=0,cb(1)
		   cblower(i,j,k)=1
		   cbupper(i,j,k)=0
		enddo
	     enddo
	  enddo

c         place atoms into cube boxes
	  do i=1,natm
	    if (rad(i).gt.0.0) then
              x=(crd(1,i)-cmin(1))*cpa
              ix=int(x)
              y=(crd(2,i)-cmin(2))*cpa
              iy=int(y)
              z=(crd(3,i)-cmin(3))*cpa
              iz=int(z)
              if (ix.le.0 .or. iy.le.0 .or. iz.le.0 .or.
     &             ix.ge.cb(1) .or. iy.ge.cb(2) .or. iz.ge.cb(3)) then
c                write(6, *) 'Atom #', i, ' outside cube bounds!'
c                write(6, *) 'x, y, z = ', (crd(k,i), k=1,3)
c                write(6, *) 'ix, iy, iz = ', ix, iy, iz
                icbn(1,i)=-1
              else
c               count this atom in all 27 nearest cubes
c                write(6, *) 'atom # ', i
                do jz=iz-1,iz+1
                  do jy=iy-1,iy+1
                    do jx=ix-1,ix+1
                      cbupper(jx,jy,jz)=cbupper(jx,jy,jz)+1
                    enddo
                  enddo
                enddo
                icbn(1,i)=ix
                icbn(2,i)=iy
                icbn(3,i)=iz
              endif
              call flush(6)
            endif
          enddo

c       calculate cummulative sum into cblower as lower bound for that box
	icum=1
	do iz=0,cb(3)
	   do iy=0,cb(2)
	      do ix=0,cb(1)
c                check if this box has nearby atoms
		 if (cbupper(ix,iy,iz).gt.0) then
		    cblower(ix,iy,iz)=icum
		    icum=icum+cbupper(ix,iy,iz)
		 endif
	      enddo
	   enddo
	enddo

c       fill cbatom with pointers back to the nearby atoms
c       reuse cblower as per-box pointer to next availabe cbatom location

c       do central box first, then 26 more successively distant neighbors
c       0,0,0
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       nn=1, boxes that share a face with the central box
c       -1,0,0
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       1,0,0
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       0,-1,0
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            iy=iy-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       0,1,0
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            iy=iy+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       0,0,-1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            iz=iz-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       0,0,1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            iz=iz+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       nn=2, boxes which share only an edge with the central box
c       -1,0,-1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix-1
            iz=iz-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       1,0,-1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix+1
            iz=iz-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       0,-1,-1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            iy=iy-1
            iz=iz-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       0,1,-1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            iy=iy+1
            iz=iz-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       -1,0,1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix-1
            iz=iz+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       1,0,1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix+1
            iz=iz+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       0,-1,1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            iy=iy-1
            iz=iz+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       0,1,1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            iy=iy+1
            iz=iz+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       -1,-1,0
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix-1
            iy=iy-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       1,-1,0
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix+1
            iy=iy-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       -1,1,0
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix-1
            iy=iy+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       1,1,0
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix+1
            iy=iy+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       nn=3, boxes which share only a corner with the central box
c       -1,-1,-1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix-1
            iy=iy-1
            iz=iz-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       1,-1,-1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix+1
            iy=iy-1
            iz=iz-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       -1,1,-1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix-1
            iy=iy+1
            iz=iz-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       1,1,-1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix+1
            iy=iy+1
            iz=iz-1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       -1,-1,1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix-1
            iy=iy-1
            iz=iz+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       1,-1,1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix+1
            iy=iy-1
            iz=iz+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       -1,1,1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix-1
            iy=iy+1
            iz=iz+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo
c       1,1,1
	do i=1,natm
          ix=icbn(1,i)
          if (ix .gt. 0 .and. rad(i) .gt. 0.0) then
            iy=icbn(2,i)
            iz=icbn(3,i)
            ix=ix+1
            iy=iy+1
            iz=iz+1
            cbatom(cblower(ix,iy,iz))=i
            cblower(ix,iy,iz)=cblower(ix,iy,iz)+1
          endif
	enddo

c       reset cblower array back to lower bound pointers into cbatom
	icum=1
	do iz=0,cb(3)
	   do iy=0,cb(2)
	      do ix=0,cb(1)
		 if(cbupper(ix,iy,iz).gt.0)then
		    cblower(ix,iy,iz)=icum
		    icum=icum+cbupper(ix,iy,iz)
		    cbupper(ix,iy,iz)=icum-1
		 endif
	      enddo
	   enddo
	enddo

	icum=icum-1
	i_icbn= memalloc(i_icbn,0)
      return
      end
