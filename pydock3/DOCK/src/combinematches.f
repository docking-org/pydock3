c this routine looks at every combination of the sphmatches

      subroutine combinematches(sphmatches, ligsphcoords, recsphcoords,
     &        disttol, disl, disr,
     &        maxpts, maxcol, maxctr, maxor, maxatmpos, 
     &        minnodes, maxnodes, 
     &        match, 
     &        hash, hashcount)

      use matchtype

      implicit none

c XXX: GFORTRAN needs these to come first
c max parameter inputs, see max.h
      integer maxpts, maxcol, maxctr, maxor, maxatmpos
      integer minnodes, maxnodes !controls how many spheres to use for matching
c first all the inputs
      ! XXX: GFortran complains about 2. being a real
      ! BUG: Changing this probably added a bug!!!!! (TS)
      !integer sphmatches(maxpts*maxpts/2., 2) !1 is rec, 2 is lig, each pair is a pair
      integer sphmatches((maxpts*maxpts/2)+1, 2) !1 is rec, 2 is lig, each pair is a pair
        !of spheres that meet the distance & sometimes color requirements
      real ligsphcoords(3, maxpts) !ligand sphere centers
      real recsphcoords(3, maxatmpos) !receptor sphere centers
      real disttol
      real disl(maxpts, maxpts) !distance table for ligand spheres
      real disr(maxpts, maxpts) !distance table for receptor spheres
c outputs!
      type(matcht), intent(inout) :: match
      integer hash(maxor, maxctr) !hash of the spheres for each orientation
      integer hashcount !count of number of things in hash (same as nmatch)
      
c local variables
      integer sphpairs(maxctr, 2)    !put the sphere pairs here lig=1, rec=2
      integer places(maxctr) !where in the list we are, since we walk through
       !the list of sphmatches many times (sphmatches ^ maxnodes) times
       !maxnodes can't be more than maxctr
      integer placeindex !where in the places list we are/loop variable
      integer placeback !backwards index of places, used in loop
      integer nodeindex !count from minnodes to maxnodes
      logical done_next, continue_loop, process_this !used to control complicated loop 

c functions
      logical otherdistcheck !checks positions 2-3, 2-4, 3-4 etc, true if good
      logical degeneracycheck !checks hash of lig, rec spheres

c have to make all possible combinations of the first pair and all other pairs
c all combinations between minnodes and maxnodes should be saved by calling
c set_rotation for each one.
c first version of this algorithm finds all pairs A-B, A-C, A-D, etc
c where A is the firts pair is sphmatches always, and A-B is a valid pair where
c valid means color & distance are satisfied. B-C does not have to be within
c the proper distance in this version of the algorithm.
c minnodes & maxnodes controls how many difference B, C, D, etc. are added
      !add the original pair first
      sphpairs(1, 1) = sphmatches(1, 1)
      sphpairs(1, 2) = sphmatches(1, 2)
      do placeindex = 1, maxnodes !1 is always used, so walk from 2 to maxnodes
        places(placeindex) = placeindex !2=2, 3=3, etc. we want the unique 
        !starting point for each run through of maxnodes
      enddo
c loop runs with 1-2-3-4 first time, then
c                1-2-3-5 then
c                1-2-3-6 (suppose there are only 6) then
c                1-2-4-5 then
c                1-2-4-6 then
c                1-2-5-6 then
c                1-3-4-5 etc. places is the list being described here
      if (sphmatches(places(minnodes), 1) .gt. 0) then !otherwise there aren't
        !enough sphmatches to actually satisfy minnodes, so we can just quit
        continue_loop = .true.
        do while (continue_loop)  !once this is true, it means
          !every possible combination has been done and we can quit
          if ((sphmatches(places(2), 1) .eq. 0) .or. 
     &            (sphmatches(places(2), 2) .eq. 0)) then !we can quit entirely
            continue_loop = .false.
          endif
          if (continue_loop) then !otherwise we're out of the loop
            process_this = .true.
            do placeindex = 2, maxnodes !copy sphmatches into sphpairs
              sphpairs(placeindex, 1) = sphmatches(places(placeindex),1) !ligsph
              sphpairs(placeindex, 2) = sphmatches(places(placeindex),2) !recsph
              if ((sphpairs(placeindex, 1) .eq. 0) .or. 
     &            (sphpairs(placeindex, 2) .eq. 0)) then !we can quit because
                !we have processed every valid set of nodes
                process_this = .false.
c                exit !break out of this do loop
              endif
            enddo
            do nodeindex = minnodes, maxnodes !for every valid number of nodes
              !this actually makes the matching sphere transformation & rotation
              !calculation happen, saving it in coml, comr, and rot
              if (process_this) then !if any 0s, don't do it
                if (otherdistcheck(nodeindex, sphpairs, 
     &               disl, disr, disttol, 
     &               maxpts, maxctr, maxatmpos)) then !okay to do the rest
                  if (degeneracycheck(nodeindex, sphpairs,
     &               hash, hashcount,
     &               maxpts, maxctr, maxor)) then
                    call set_rotation(nodeindex, sphpairs,   
     &                maxor, maxpts, maxctr, maxatmpos, ligsphcoords,  
     &                recsphcoords, match)
                  endif
                endif
              endif
            enddo
            done_next = .false.
            !have to increment the last thing in places, 
            ! or reset many things in placeindex and 
            ! increment a lower-indexed pointer.
            placeback = maxnodes !the first thing we'll try to increment
            do while (.not. done_next)
              places(placeback) = places(placeback) + 1
              if (placeback .eq. 2) then !if we've gotten to placeback=2,
                !then we don't do this, because the loop might finally be over
                done_next = .true.
              else !we want to check and see if we should decrement placeback
                !and start a new loop of the lower number of places
                if (sphmatches(places(placeback), 1) .le. 0) then !gone too far
                  placeback = placeback - 1 !backwards increment here
                  !don't know if we're done yet, so we reenter the do while loop 
                  !and places(placeback) will be incremented
                else
                  done_next = .true. !just incrementing this is fine
                endif
              endif
            enddo
            !after this we may need to move later pointers if we reset earlier 
            do while (placeback .lt. maxnodes)  ! have to advance these,
              !strictly less than since if just the one equal to maxnodes it 
              !has already been incremented
              placeback = placeback + 1 
              !now set the current place to start one higher than the one that 
              !was reset. possibly this happens more than once.
              places(placeback) = places(placeback - 1) + 1
            enddo
          endif
        enddo
      endif

      return
      end
