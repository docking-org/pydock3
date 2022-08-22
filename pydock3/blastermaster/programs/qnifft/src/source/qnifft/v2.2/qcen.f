	subroutine qcen(natom,atcrd1,atrad1,atcrg1,cqplus,cqmin,qplus,qmin)
c replaces getpdb when qnifft called as a subroutine- transfers atom coords, radii, charge
c------------------------------------------------------------------
	include 'qdiffpar.h'
c------------------------------------------------------------------
	real*4 atcrd1(3,natom),atrad1(natom),atcrg1(natom)
	dimension cqplus(3),cqmin(3)
c-----------------------------------------------------------------
	print *,'transferring xyz,r and q for # atoms ',natom
	qnet = 0.
	nqass = 0
	qplus = 0.
	qmin = 0.
	do k = 1,3
	  cqplus(k) = 0.
	  cqmin(k) = 0.
	end do
	do i = 1,natom
	  print *,'at xyz,r,q: ',atcrd1(1,i),atrad1(i),atcrg1(i)
	  atcrg(i) = atcrg1(i)
	  chrgv = atcrg1(i)
	  atrad(i) = atrad1(i)
	  atcrd(1,i) = atcrd1(1,i)
	  atcrd(2,i) = atcrd1(2,i)
	  atcrd(3,i) = atcrd1(3,i)
	  qnet = qnet + chrgv
	  if(chrgv.gt.0.0) then
	    qplus = qplus + chrgv
	    do k = 1,3
	      cqplus(k) = cqplus(k) + atcrd(k,i)*chrgv
	    end do
	  end if
	  if(chrgv.lt.0.0) then
	    qmin = qmin + chrgv
	    do k = 1,3
	      cqmin(k) = cqmin(k) + atcrd(k,i)*chrgv
	    end do
	  end if
	end do
      if(qplus.ne.0.0) then
        do k = 1,3
          cqplus(k) = cqplus(k)/qplus
        end do
      end if
      if(qmin.ne.0.0) then
        do k = 1,3
          cqmin(k) = cqmin(k)/qmin
        end do
      end if
      write(6,*)'net assigned source charge              : ',qnet
      write(6,*)'assigned positive charge                : ',qplus
      write(6,*)'centred at (A)  :',cqplus
      write(6,*)'assigned negative charge                : ',qmin
      write(6,*)'centred at (A)  :',cqmin
      write(6,*)'   '
	return
	end
