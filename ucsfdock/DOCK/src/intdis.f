********************************************************************************
*--------------------------------intdis----------------------------------------*
* this subroutine calculates the pairwise internal distances of points in a    *
* cartesian coordinate file.    BKS                                            *
* improved only a tiny bit for readability/code style. RGC 09/2011.
* note that many 'temporary' matrices here are returned and used later.
* note that the idis and distmp matrices are identical here, one is modified
*  later by other programs (distmp), the other (idis) is not and is used
* note that ctrtmp seems a bit silly here, but is used later by reord, etc
*  so don't delete it or change it.
* lots of stuff removed. now just does distances between every pair of points
*  RGC 2012
********************************************************************************
      subroutine intdis(numctr, corr, idis, dmax, maxpts)

      implicit none

c XXX: GFortran needs max lengths to come first
c          idis: contains pair-wise internal distances. (return value)
      real, intent(inout) :: dmax
c          dmax: maximum internal distance. (return value)
      integer, intent(in) :: maxpts !maxpts is passed in here
      integer, intent(in) :: numctr 
c        numctr: number of centers. (input)
      real, dimension(3, maxpts), intent(in) :: corr
c          corr: cartesian coordinates. (input)
      real, dimension(maxpts, maxpts), intent(inout) :: idis
c          idis: ligand sph-sph distances

      real sum, sqrtsum
c        sum and sqrtsum are temporary distance labels.
      integer counti, countj, coord !loop variables

c      initialize.
      dmax=0.
      sqrtsum=0.
c clear many matrices of previous values, only to the point we'll use
      do countj = 1, numctr + 1
        do counti = 1, numctr + 1
          idis(counti, countj)   =0.
        enddo
      enddo
c main loop here, compute distances and save them several places
      do counti = 1, numctr-1
        do countj = counti + 1, numctr !funny loop to only do half the matrix
          sum = 0.
          do coord = 1, 3 !over every xyz, compute sum of distances
            sum = sum+(corr(coord, counti)-corr(coord, countj))**2
          enddo
          sqrtsum = sqrt(sum)
          idis(counti, countj) = sqrtsum
          idis(countj, counti) = sqrtsum !copy to both places in matrix
          if (sqrtsum .gt. dmax) then !dmax should be the max value anywhere
            dmax = sqrtsum
          endif
        enddo
      enddo
      return !many values modified in place, see comments above
      end
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
c
