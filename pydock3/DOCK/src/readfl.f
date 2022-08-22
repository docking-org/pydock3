c
c keyword subroutines
c written by Michael L. Connolly
c Department of Pharmaceutical Chemistry
c University of California at San Francisco
c Copyright 1993 by the Regents of the University of California
c last revised February 14, 1994
c
c Modified to support force-field score minimization DAG 4/94.
c DAG modifications 7/94:
c      added critical_clusters (ACG)
c      added solvation database support (BKS)
c DAG modifications 8/94:
c      updated keywords to more meaningful things (group decision)
c At some point critical_clusters were removed
c
c read option or parameter file
c
      subroutine readfl(runit, fname, OUTDOCK, phimap0, recdes0,
     &    ligdes0, options0, spheres0, MAXLNS, lines, errflg)

      use errorformats
      use phimaptype
      use optionstype
      use spheres

      implicit none

      integer runit !the file number/read unit/etc
      character*(*) fname
c outdock file specifier, so data can be written out if there are problems
      integer, intent(in) :: OUTDOCK
c these are user-defined types that contain information about grids (usually)
      type(phimap), intent(inout) :: phimap0, recdes0, ligdes0
      type(options), intent(inout) :: options0
      type(spherest), intent(inout) :: spheres0
c where INDOCK is stored...
      integer, intent(in) :: MAXLNS
      character (len=300), dimension(MAXLNS), intent(inout) :: lines
      logical, intent(inout) :: errflg

      integer nlines, fline, lline
      character (len=30) filvrn
      logical fail

      nlines = 0
      fline = 1
      call readln(runit, fname, fail, nlines, MAXLNS, lines)
      if (fail) then
        write(OUTDOCK, HALT2) 'Cannot open parameter file: ', fname
        write(OUTDOCK, HALT0)
        stop
      endif
      lline = nlines
      if (lline .gt. fline) then
        call getvrn(fline, filvrn, MAXLNS, lines)
        if (filvrn .eq. '3.8') then
          fline = fline + 1
c         replace # and following characters by blanks
          call blankc(fline, lline, MAXLNS, lines)
          call readkw(fline, lline, OUTDOCK, phimap0, recdes0,
     &        ligdes0, options0, spheres0, MAXLNS, lines, errflg)
        else
          write(OUTDOCK, HALT2) 
     &      'Unknown version in parameter file: ', fname
          write(OUTDOCK, HALT0)
          stop
        endif
      endif
      close(runit) !actually close input files

      return
      end
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
