	subroutine smtheps(ismth,smthlen)
c
c   This subroutine is intended to smooth the abruptly changing
c   of dielectric constant value (eps) on the boundary of two
c   dielectric media.  (July 27, 1995) 
c

	include 'qdiffpar.h'
	dimension epstmp(ngrid,ngrid,ngrid,3)
      integer ia, ib, ic, id
	real xtemp, xcoe1, xcoe2
c------------------------------------------------------
c  check eps values before smoothing 
c-------------------------------------------------------         
c      open(12, file='check.value', status='unknown')
c      write(12,*) 'Values before smtheps'
c	write(12,*) '(30,33,33,1)', 561.*epsmap(30,33,33,1)
c	write(12,*) '(30,32,33,2)', 561.*epsmap(30,32,33,2)
c	write(12,*) '(30,33,33,2)', 561.*epsmap(30,33,33,2)
c	write(12,*) '(31,32,33,2)', 561.*epsmap(31,32,33,2)
c	write(12,*) '(31,33,33,2)', 561.*epsmap(31,33,33,2)
c	write(12,*) '(30,33,33,3)', 561.*epsmap(30,33,33,3)
c	write(12,*) '(30,33,32,3)', 561.*epsmap(30,33,32,3)
c	write(12,*) '(31,33,32,3)', 561.*epsmap(31,33,32,3)
c	write(12,*) '(31,33,33,3)', 561.*epsmap(31,33,33,3) 
c      close(12)
c---------------------------------------------------------
c     set # of nearest neighbors to average over
c
	if(ismth.eq.1) then
	  goto 9981
      end if
	if(ismth.eq.2) then
	  goto 1515
      end if
c----------------------------------------------------------------------
c   averaging over 9 points
c
9981	xcoe1=exp(-1./(2.*scale*scale*smthlen*smthlen))
	print *,'neighbor points: 9    coefficient: ', xcoe1
	do i=2, igrid-1
	  do j=2, igrid-1
	    do k=2, igrid-1
c
		  xtemp=(1/epsmap(i,j,k,1)+xcoe1/epsmap(i,j,k,2)
     6	      +xcoe1/epsmap(i,j-1,k,2)+xcoe1/epsmap(i+1,j,k,2)
     6            +xcoe1/epsmap(i+1,j-1,k,2)
     6 	      +xcoe1/epsmap(i,j,k,3)+xcoe1/epsmap(i,j,k-1,3)
     6	      +xcoe1/epsmap(i+1,j,k,3)+xcoe1/epsmap(i+1,j,k-1,3))
		  epstmp(i,j,k,1)=(1.+8.*xcoe1)/xtemp
              xtemp=(1/epsmap(i,j,k,2)+xcoe1/epsmap(i,j,k,1)
     6            +xcoe1/epsmap(i-1,j,k,1)+xcoe1/epsmap(i-1,j+1,k,1)
     6            +xcoe1/epsmap(i,j+1,k,1)
     6            +xcoe1/epsmap(i,j,k,3)+xcoe1/epsmap(i,j,k-1,3)
     6            +xcoe1/epsmap(i,j+1,k,3)+xcoe1/epsmap(i,j+1,k-1,3))
              epstmp(i,j,k,2)=(1.+8.*xcoe1)/xtemp
              xtemp=(1/epsmap(i,j,k,3)+xcoe1/epsmap(i,j,k,2)
     6            +xcoe1/epsmap(i,j,k,1)+xcoe1/epsmap(i,j-1,k,2)
     6            +xcoe1/epsmap(i-1,j,k,1)
     6            +xcoe1/epsmap(i,j,k-1,2)+xcoe1/epsmap(i,j,k-1,1)
     6            +xcoe1/epsmap(i,j-1,k-1,2)+xcoe1/epsmap(i-1,j,k-1,1))
              epstmp(i,j,k,3)=(1.+8.*xcoe1)/xtemp
          enddo
        enddo
      enddo
	goto 1995
c-----------------------------------------------------------------------
c   averaging over 15 points
c
1515	xcoe1=exp(-1./(2.*scale*scale*smthlen*smthlen))
	xcoe2=exp(-1./(scale*scale*smthlen*smthlen))
	print *,'neighbor points: 15    coefficient: ', xcoe1, xcoe2
	do i=2, igrid-1
	  do j=2, igrid-1
	    do k=2, igrid-1
c
		  xtemp=(1/epsmap(i,j,k,1)+xcoe1/epsmap(i,j,k,2)
     6	      +xcoe1/epsmap(i,j-1,k,2)+xcoe1/epsmap(i+1,j,k,2)
     6            +xcoe1/epsmap(i+1,j-1,k,2)
     6 	      +xcoe1/epsmap(i,j,k,3)+xcoe1/epsmap(i,j,k-1,3)
     6	      +xcoe1/epsmap(i+1,j,k,3)+xcoe1/epsmap(i+1,j,k-1,3)
     6            +xcoe2/epsmap(i-1,j,k,1)+xcoe2/epsmap(i,j+1,k,1)
     6            +xcoe2/epsmap(i+1,j,k,1)
     6            +xcoe2/epsmap(i,j-1,k,1)+xcoe2/epsmap(i,j,k+1,1)
     6            +xcoe2/epsmap(i,j,k-1,1))
		  epstmp(i,j,k,1)=(1.+8.*xcoe1+6.*xcoe2)/xtemp
              xtemp=(1/epsmap(i,j,k,2)+xcoe1/epsmap(i,j,k,1)
     6            +xcoe1/epsmap(i-1,j,k,1)+xcoe1/epsmap(i-1,j+1,k,1)
     6            +xcoe1/epsmap(i,j+1,k,1)
     6            +xcoe1/epsmap(i,j,k,3)+xcoe1/epsmap(i,j,k-1,3)
     6            +xcoe1/epsmap(i,j+1,k,3)+xcoe1/epsmap(i,j+1,k-1,3)
     6            +xcoe2/epsmap(i+1,j,k,2)+xcoe2/epsmap(i-1,j,k,2)
     6            +xcoe2/epsmap(i,j+1,k,2)+xcoe2/epsmap(i,j-1,k,2)
     6            +xcoe2/epsmap(i,j,k+1,2)+xcoe2/epsmap(i,j,k-1,2))
              epstmp(i,j,k,2)=(1.+8.*xcoe1+6.*xcoe2)/xtemp
              xtemp=(1/epsmap(i,j,k,3)+xcoe1/epsmap(i,j,k,2)
     6            +xcoe1/epsmap(i,j,k,1)+xcoe1/epsmap(i,j-1,k,2)
     6            +xcoe1/epsmap(i-1,j,k,1)
     6            +xcoe1/epsmap(i,j,k-1,2)+xcoe1/epsmap(i,j,k-1,1)
     6            +xcoe1/epsmap(i,j-1,k-1,2)+xcoe1/epsmap(i-1,j,k-1,1)
     6            +xcoe2/epsmap(i+1,j,k,3)+xcoe2/epsmap(i-1,j,k,3)
     6            +xcoe2/epsmap(i,j+1,k,3)+xcoe2/epsmap(i,j-1,k,3)
     6            +xcoe2/epsmap(i,j,k+1,3)+xcoe2/epsmap(i,j,k-1,3))
              epstmp(i,j,k,3)=(1.+8.*xcoe1+6.*xcoe2)/xtemp
          enddo
        enddo
      enddo
c---------------------------------------------------------------------------	
c   write the averaging results back to the eps array
c
1995	do i=2, igrid-1
	  do j=2, igrid-1
	    do k=2, igrid-1
		do m=1, 3
		  epsmap(i,j,k,m)=epstmp(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
c------------------------------------------------------------
c   check eps values after smooothing
c------------------------------------------------------------
c      open(12, file='check.value', status='unknown')
c      write(12,*) 'Values after smtheps'
c	write(12,*) '(30,33,33,1)', 561.*epsmap(30,33,33,1)
c	write(12,*) '(30,32,33,2)', 561.*epsmap(30,32,33,2)
c	write(12,*) '(30,33,33,2)', 561.*epsmap(30,33,33,2)
c	write(12,*) '(31,32,33,2)', 561.*epsmap(31,32,33,2)
c	write(12,*) '(31,33,33,2)', 561.*epsmap(31,33,33,2)
c	write(12,*) '(30,33,33,3)', 561.*epsmap(30,33,33,3)
c	write(12,*) '(30,33,32,3)', 561.*epsmap(30,33,32,3)
c	write(12,*) '(31,33,32,3)', 561.*epsmap(31,33,32,3)
c	write(12,*) '(31,33,33,3)', 561.*epsmap(31,33,33,3) 
c      close(12)
c-------------------------------------------------------------
      return
	end
