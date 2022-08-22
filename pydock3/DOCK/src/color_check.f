c checks to see if the input spheres have an acceptable color match or not

      logical function color_check(ligsphcolor, recsphcolor, colortable, 
     &                 ligand_sph, rec_sph, maxpts, maxcol)

      implicit none

      ! XXX: GFortran needs this to come first (TS)
      integer maxpts, maxcol !max.h parameters, max # of points and max # colors
      integer ligsphcolor(maxpts) !the color of the rigid component spheres
      integer recsphcolor(maxpts) !the color of the receptor matching spheres
      integer colortable(maxcol, maxcol)
      integer ligand_sph, rec_sph !the spheres to check for match

      !variables
      integer lcol, rcol !the actual color of the spheres

      color_check = .true.
      lcol = ligsphcolor(ligand_sph)
      rcol = recsphcolor(rec_sph)
c     chemical matching used only if both are colored
      if ((lcol .gt. 0) .and. (rcol .gt. 0)) then
        if (colortable(lcol, rcol) .eq. 0) then
          color_check = .false.
        endif
      endif

      return !color_check as the return value
      end
