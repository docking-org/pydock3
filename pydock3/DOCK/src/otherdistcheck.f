c checks to see if the spheres have the correct distance between them, within 
c the tolerance. 1-2, 1-3, 1-4, 1-5 etc have been checked already but
c 2-3, 2-4, 2-5, 3-4, 3-5, 4-5 etc need checked.
c still uses disl & disr to look up distances

      logical function otherdistcheck(nodetotal, sphpairs, 
     &          disl, disr, disttol, maxpts, maxctr, maxatmpos)

      implicit none

c XXX: GFortran needs these first
c bunch of maximum parameters
      integer maxpts !maximum number of sphere points
      integer maxctr !maximum number of sphere centers
      integer maxatmpos !maximum atom positions

c first are all the input parameters
      integer nodetotal !match sphere parameters, how many sphere pairs to use
      integer sphpairs(maxctr, 2) !match sphere parameters
      real disl(maxpts, maxpts) !distance table for ligand spheres
      real disr(maxpts, maxpts) !distance table for receptor spheres
      real disttol !distance tolerance
c temporary variables
      integer nodecount, othercount !counts through the pairs
      real ligdist, recdist !temporary distance holders

      otherdistcheck = .true. !default is true
      do nodecount = 2, nodetotal !1 to whatever is already done
        do othercount = nodecount + 1, nodetotal !other count
          ligdist = disl(sphpairs(nodecount, 1), 
     &                   sphpairs(othercount, 1))
          recdist = disr(sphpairs(nodecount, 2), 
     &                   sphpairs(othercount, 2))
          if (abs(ligdist - recdist) .gt. disttol) then !bad
            otherdistcheck = .false. !can quit now, this triple/quadruple 
              !doesn't fulfill all the distance requirements
            return !otherdistcheck as the return value, no need to continue
          endif
        enddo
      enddo

      return !otherdistcheck as the return value
      end
