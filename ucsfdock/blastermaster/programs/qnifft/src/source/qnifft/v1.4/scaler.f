c-------------------------------------------------------
      subroutine gtoc(g,c)
c
c converts grid to real coordinates
c
	include 'qdiffpar.h'
c---------------------------------------------------

c
      dimension g(3),c(3)
	goff = (igrid + 1.)/2.
      do 9000 i = 1,3
         c(i) = (g(i) - goff)/scale + oldmid(i)
c      end do
9000	continue
      return
      end

c-------------------------------------------------------
c-------------------------------------------------------
      subroutine ctog(c,g)
c
c converts real to grid coordinates
c
	include 'qdiffpar.h'
c---------------------------------------------------
c
      dimension g(3),c(3)

	goff = (igrid + 1.)/2.
      do 9000 i = 1,3
         g(i) = (c(i) - oldmid(i))*scale + goff
c      end do
9000	continue
      return
      end

	subroutine sclit(atcrd,atrad,natom,igrid,perfil,scale,oldmid,
     &border,isize,icen)
	dimension atcrd(3,natom),atrad(natom)
	dimension oldmid(3),cenmol(3)
	character*15 ascl,axcen,aycen,azcen
	character*79 sclin
	logical icen
c--------------------------------------------------
c
c find extent of molecule, including the radii
c
	cmin1=6000
	cmin2=6000
	cmin3=6000
	cmax1=-6000
	cmax2=-6000
	cmax3=-6000
	do ix=1,natom
	  if(atrad(ix).ge.0.)then
	    cmin1=amin1(cmin1,atcrd(1,ix)-atrad(ix))
	    cmin2=amin1(cmin2,atcrd(2,ix)-atrad(ix))
	    cmin3=amin1(cmin3,atcrd(3,ix)-atrad(ix))
	    cmax1=amax1(cmax1,atcrd(1,ix)+atrad(ix))
	    cmax2=amax1(cmax2,atcrd(2,ix)+atrad(ix))
	    cmax3=amax1(cmax3,atcrd(3,ix)+atrad(ix))
	  end if
	end do
	cran1=cmax1-cmin1
	cran2=cmax2-cmin2
	cran3=cmax3-cmin3 
	bran=amax1(cran1,cran2)
	bran=amax1(bran,cran3)
	cenmol(1)=(cmax1+cmin1)/2
	cenmol(2)=(cmax2+cmin2)/2
	cenmol(3)=(cmax3+cmin3)/2
	midg=(igrid+1)/2
	rmidg=midg
	gran=igrid-1.

	write(6,*)'  '
	if(isize.eq.1)then
c
c scale, oldmid given explicitly
c
	  if(icen)then
	    print *,'scaling by specified scale'
	    print *, 'WARNING: grid center x/y/z incomplete: using molecule center'
	    oldmid(1)=cenmol(1)
	    oldmid(2)=cenmol(2)
	    oldmid(3)=cenmol(3)
	  else
	    print *,'scaling by specified scale and midpoint'
	  endif
	  perfil = scale*100.*bran/gran
	else if(isize.eq.2)then
c
c using minimum border of solvent
c
c	  print *,'using solvent border of (A) : ',border
	  print *,'scaling by solvent border'
	  oldmid(1)=cenmol(1)
	  oldmid(2)=cenmol(2)
	  oldmid(3)=cenmol(3)
	  scale=gran/(bran+2.*border) 
	  perfil = scale*100.*bran/gran
	else if(isize.eq.3)then
c
c using percent fill, put middle of molecule in middle of box
c
c	  print *,'longest dimension of molecule(s) as % of box length: ',perfil
	  print *,'scaling by percent fill'
	  oldmid(1)=cenmol(1)
	  oldmid(2)=cenmol(2)
	  oldmid(3)=cenmol(3)
	  scale=gran*perfil/(100.*bran) 
	end if
c
c calculate borders
c
	bordmn = 1.e6
	brdxup = gran/2./scale - (cmax1 - oldmid(1))
	bordmn = min(brdxup,bordmn)
	brdyup = gran/2./scale - (cmax2 - oldmid(2))
	bordmn = min(brdyup,bordmn)
	brdzup = gran/2./scale - (cmax3 - oldmid(3))
	bordmn = min(brdzup,bordmn)
	brdxlw = gran/2./scale + (cmin1 - oldmid(1))
	bordmn = min(brdxlw,bordmn)
	brdylw = gran/2./scale + (cmin2 - oldmid(2))
	bordmn = min(brdylw,bordmn)
	brdzlw = gran/2./scale + (cmin3 - oldmid(3))
	bordmn = min(brdzlw,bordmn)
c
c write details
c
	write(6,*)'xmin,xmax     (A)   :',cmin1,cmax1
	write(6,*)'ymin,ymax     (A)   :',cmin2,cmax2
	write(6,*)'zmin,zma      (A)   :',cmin3,cmax3
	write(6,*)'x,y,z range   (A)   : ',cran1,cran2,cran3
	write(6,*)'object centre (A)   : ',cenmol
c	write(6,*)'scale   (grids/A)   : ',scale
c	write(6,*)'grid centre   (A)   : ',oldmid
c
c write scale in format that can be redirected to parameter file for input
c
	write(ascl,'(g15.7)')scale
	call elb(ascl,15)
	write(axcen,'(g15.7)')oldmid(1)
	call elb(axcen,15)
	write(aycen,'(g15.7)')oldmid(2)
	call elb(aycen,15)
	write(azcen,'(g15.7)')oldmid(3)
	call elb(azcen,15)
	write(6,*)'scale=',ascl
	write(6,*)'xcen=',axcen,' ycen=',aycen,' zcen=',azcen
	write(6,*)'low,up x borders (A):',brdxlw,brdxup
	write(6,*)'low,up y borders (A):',brdylw,brdyup
	write(6,*)'low,up z borders (A):',brdzlw,brdzup
	print *,'longest dimension of molecule(s) as % of box length: ',perfil
	print *,'minimum solvent border (A) : ',bordmn
	write(6,*)'  '
	return
	end
