c constructs adjacency lists of point-point pairs with distances within the
c tolerance (and correctly colored if desired) a la Ewing & Kuntz J Comp Chem
c no histograms/bins/approximations used
c RGC 09/2011
      subroutine sphsfromdists(ligstart, recstart, disl, 
     &         disr, maxlig, maxrec, disttol, 
     &         colormatch, ligsphcolor, recsphcolor, 
     &         colortable, maxpts, maxcol, sphmatches)

      implicit none

c XXX: GFortran needs these first
c max parameter inputs
      integer maxpts, maxcol !max parameters passed in, see max.h
c first all the inputs
      integer ligstart !starting (seed) sphere for ligand side
      integer recstart !starting (seed) sphere for receptor side
      real disl(maxpts, maxpts) !distance table for ligand spheres
      real disr(maxpts, maxpts) !distance table for receptor spheres
      integer maxlig, maxrec !number of entries in the disl & disr matrices
      real disttol !distance tolerance allowed
c color inputs
      logical colormatch !whether to do color matching or not
      integer ligsphcolor(maxpts) !the color of the rigid component spheres
      integer recsphcolor(maxpts) !the color of the receptor matching spheres
      integer colortable(maxcol, maxcol)
c outputs!
      !integer sphmatches(maxpts*maxpts/2., 2) !1 is rec, 2 is lig, each pair is a pair
      integer sphmatches((maxpts*maxpts/2)+1, 2) !1 is rec, 2 is lig, each pair is a pair
        !of spheres that meet the distance & sometimes color requirements

c temporary variables used for counting through the various lists/data structs
      integer matchcount, matchtmp
      integer curligpt, currecpt
      real ligdist, distdiff !difference in the distances between 2 pairs of points

c function used to check color matching
      logical color_check !function that returns true if the colors match

      do matchtmp = 1, (maxpts*maxpts/2)+1 !initialization step, clearing output data
        sphmatches(matchtmp, 1) = 0 !clear the whole array
        sphmatches(matchtmp, 2) = 0 !clear the whole array
      enddo
      matchcount = 0 !how many confirmed matches have been found so far
      matchcount = matchcount + 1 !going to add one
      sphmatches(matchcount, 1) = ligstart !always add this pair
      sphmatches(matchcount, 2) = recstart !always add this pair
      do curligpt = 1, maxlig !for all ligand spheres
        if (curligpt .ne. ligstart) then !except the starting point
          ligdist = disl(ligstart, curligpt) !dist from start to current, lig
          do currecpt = 1, maxrec !for all receptor spheres
            if (currecpt .ne. recstart) then !except the starting point
              distdiff = abs(disr(recstart, currecpt) - ligdist) !|rec - lig|
              !this is the difference in the distances, if below threshold save
              if (distdiff .le. disttol) then
                !still have to check colors maybe
                if (colormatch) then !first check color if we want to
                  if (.not. color_check(ligsphcolor, recsphcolor,
     &     colortable, curligpt, currecpt, maxpts, maxcol)) then
                    cycle !color was not okay, so go to next iteration
                  endif
                endif
                !if we get here it means the colors were okay
                matchcount = matchcount + 1
                sphmatches(matchcount, 1) = curligpt
                sphmatches(matchcount, 2) = currecpt
              endif
            endif
          enddo
        endif
      enddo
      !matchcount is now equal to how many added, sphmatches contains
      !pairs and then 0s once all pairs exhausted

      return
      end
