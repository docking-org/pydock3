c full matching algorithm, no longer uses bins/histograms, makes all possible
c sph-sph matches in terms of only the distance tolerance, a la Ewing & Kuntz.
c like all matching functions, this is called once per pair of initial spheres
c results are added to nmatch, coml, comr, rot
c 09/2011 RGC
c 
      subroutine fullmatch(ligsphcount, ligsphcolor, ligcount, reccount, 
     &   recsphcolor, colortable, colormatch, disttol, disl, disr,
     &   ligsphcoords, recsphcoords, maxlig, maxrec,
     &   maxpts, maxcol, maxor, maxctr, maxatmpos, 
     &   minnodes, maxnodes,
     &   match,
     &   hash, hashcount)

      use matchtype

      implicit none

c XXX: GFortran needs these to come first
c max parameters (see max.h)
      integer maxpts, maxcol, maxor, maxctr, maxatmpos !max parameters
      integer minnodes, maxnodes !control how many to use for matching
c input variables
      integer ligsphcount !number of heavy atoms in rigid/match component or 
         !extended clouds.
      integer ligsphcolor(maxpts) !the color of the rigid component spheres
      integer ligcount !starting (seed) sphere for ligand side
      integer reccount !starting (seed) sphere for receptor side
      integer recsphcolor(maxpts) !the color of the receptor matching spheres
      integer colortable(maxcol, maxcol)
      logical colormatch !whether to do color matching or not
      real disttol !distance tolerance allowed
      real disl(maxpts, maxpts) !distance table for ligand spheres
      real disr(maxpts, maxpts) !distance table for receptor spheres
      real ligsphcoords(3, maxpts) !ligand sphere centers
      real recsphcoords(3, maxatmpos) !receptor sphere centers
c total count of max distances in lig,rec
      integer maxlig, maxrec !number of entries in the disl & disr matrices
c return variables
      type(matcht), intent(inout) :: match !match data
      integer hash(maxor, maxctr) !hash of the spheres for each orientation
      integer hashcount !total number of hashes, should equal nmatch

c temporary variables for use during matching
      integer count !temporary counter variable
      !XXX: Changed the 2. (real) to 2 (integer) to please Gfortran
      ! BUG: This likely introduced some sort of bug
      !integer sphmatches(maxpts*maxpts/2., 2) !maxpts times maxpts is if every pair
      integer sphmatches((maxpts*maxpts/2)+1, 2) !maxpts times maxpts is if every pair
        !hits, this could be improved or made dynamic
        !1 is rec, 2 is lig, each pair is a pair
        !of spheres that meet the distance & sometimes color requirements

      logical color_check !function that returns true if the colors match

c     chemical matching (coloring) check - do this here to not waste time
      if (colormatch) then
        if (.not. color_check(ligsphcolor, recsphcolor, colortable,
     &                 ligcount, reccount, maxpts, maxcol)) then
          return  !color check failed, just quit and return
        endif
      endif
      call sphsfromdists(ligcount, reccount, disl, disr, 
     &           maxlig, maxrec, disttol,
     &           colormatch, ligsphcolor, recsphcolor, 
     &           colortable, maxpts, maxcol, sphmatches)
      call combinematches(sphmatches, ligsphcoords, recsphcoords,
     &           disttol, disl, disr, 
     &           maxpts, maxcol, maxctr, maxor, maxatmpos, 
     &           minnodes, maxnodes, match,
     &           hash, hashcount)

      return
      end
