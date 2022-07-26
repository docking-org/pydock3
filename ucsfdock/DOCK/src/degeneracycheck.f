c makes a hash of the lig spheres and rec spheres. checks this against
c previous hashes, if new (not seen before) add to list andreturn true. 
c if seen before, return false.

      logical function degeneracycheck(nodetotal, sphpairs, 
     &          hash, hashcount,
     &          maxpts, maxctr, maxor)

      implicit none

c XXX: GFortran needs these to come first (TS)
c bunch of maximum parameters
      integer maxpts !maximum number of sphere points
      integer maxctr !maximum number of sphere centers
      integer maxor !maximum number of orientations

c first are all the input parameters
      integer nodetotal !match sphere parameters, how many sphere pairs to use
      integer sphpairs(maxctr, 2) !match sphere parameters
      integer hash(maxor, maxctr) !hash of sph numbers
      integer hashcount !how many hashes have been stored

c temporary variables
      integer nodeindex, hashindex, centerindex, tmpindex !temporary counters
      integer tmphash(maxctr) !the temporarily computed hash
      logical match_found, all_matched

      degeneracycheck = .true. !default is true
      do nodeindex = 1, nodetotal !for every node
        tmphash(nodeindex) = sphpairs(nodeindex, 1) + 
     &                       sphpairs(nodeindex, 2) * maxpts !makes uniq number
      enddo
      do hashindex = 1, hashcount !check all previously found matches
        !assume order doesn't matter, check to see if all of tmphash matches
        all_matched = .true.
        do tmpindex = 1, maxctr
          match_found = .false.
          do centerindex = 1, maxctr
            if (tmphash(tmpindex) .eq. 
     &                   hash(hashindex, centerindex)) then
              match_found = .true.
            endif
          enddo
          if (.not. match_found) then
            all_matched = .false.
            exit !break out of the do loop
          endif       
        enddo       
        if (all_matched) then !whole tmphash matched
c          write (6,*) tmphash, "tmp"
c          write (6,*) hash(hashindex, 1), hash(hashindex, 2),
c     &                hash(hashindex, 3), hash(hashindex, 4), "hash"
          degeneracycheck = .false.
          return !quit and return, match found, nothing to do
        endif
      enddo
      if (degeneracycheck) then !no match found, so add the temps 
        hashcount = hashcount + 1
        do centerindex = 1, maxctr !copy them all, even 0s
          hash(hashcount, centerindex) = tmphash(centerindex)
        enddo
      endif

      return !degeneracycheck as the return value, is still true if here
      end
