c call with a list of matching points in iasign (# points is wrtnod)
c maxor, etc are maximum list sizes from max.h
c adds data into comr, coml, and rot after figuring out the translation/rotation
c removed global variables RGC 09/2011
      subroutine set_rotation(wrtnod, iasign, maxor, maxpts, 
     &    maxctr, maxatmpos, spcorr, coords, match)

      use matchtype
      use filenums

      implicit none

c bunch of maximum parameters, these could be pruned way down
      integer maxor !maximum number of orientations we can store
      integer maxpts !maximum number of sphere points
      integer maxctr !maximum number of sphere centers
      integer maxatmpos !maximum atom positions
c first are all the input parameters
      integer wrtnod !match sphere parameters, how many sphere pairs to use
      integer iasign(maxctr, 2) !match sphere parameters
c next are output variables
      real spcorr(3, maxpts) !array of sphere center coordinates for the cluster
      real coords(3, maxatmpos) !transformed coordinates
      type(matcht), intent(inout) :: match

c temporary variables start here
      integer three1, three2 !temporary counters
      integer count !temporary counters
c temporary matrices/vectors to store rotation/translations in
      real cgl(3), cgr(3) !centers of mass
      real rottemp(3, 3) !rotation matrix
      real xor(3, maxctr)     !temporary coordinate storage, receptor spheres
      real xos(3, maxctr)  !temporary coordinate storage, ligand spheres

      match%nmatch = match%nmatch + 1 !we've found another matchset -> orientation, so ++
      if (match%nmatch .gt. maxor) then !too many orientations
        if (match%nmatch .eq. maxor + 1) then !exactly once
          write(OUTDOCK, *)
     &         'too many orientations generated, skipping.'
          write(OUTDOCK, *) 'adjust maxor in max.h', maxor
        endif
      else !only process if room to store it
c       set up for least-squares fit
        call fixpts(wrtnod, spcorr, coords, iasign, xos, xor, maxpts,
     &            maxatmpos, maxctr)
c       least-squares fit this match of atoms to sphere centers
c       every new orientation loops over this point
        call orient(wrtnod, xor, xos, cgr, cgl, rottemp)
c move the data from the temporary variables into the large lists
c nmatch is updated outside of this loop (for now)
        do three1 = 1, 3
          match%coml(three1, match%nmatch) = cgl(three1)
          match%comr(three1, match%nmatch) = cgr(three1)
          do three2 = 1, 3
            match%rot(three1, three2, match%nmatch) = 
     &          rottemp(three1, three2)
          enddo
          do count = 1, wrtnod
            match%xor(three1, count, match%nmatch) = 
     &          xor(three1, count)
            match%xos(three1, count, match%nmatch) =
     &          xos(three1, count)
          enddo
        enddo
      endif
      return
      end
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
