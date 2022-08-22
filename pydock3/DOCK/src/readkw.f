c
      subroutine readkw(fline, lline, OUTDOCK, phimap0, recdes0, 
     &    ligdes0, options0, spheres0, MAXLNS, lines, errflg)

      use errorformats
      use phimaptype
      use optionstype
      use spheres

      implicit none

      integer fline, lline
c outdock file specifier, so data can be written out if there are problems
      integer, intent(in) :: OUTDOCK
c these are user-defined types that contain information about grids (usually)
      type(phimap), intent(inout) :: phimap0, recdes0, ligdes0
      type(options), intent(inout) :: options0
      type(spherest), intent(inout) :: spheres0
c where INDOCK is stored...
      integer, intent(in) :: MAXLNS
      character (len=300), dimension(MAXLNS), intent(inout) :: lines
      logical, intent(inout) :: errflg  !problems?

      integer i, n !replace these with meaningful names
      integer lenkey
      integer pos1, pos2, pos3
      character (len=132) args
      character (len=30) keywrd

c     functions:
      integer nflds, nxtfld, nxtsp

      do i = fline, lline
        n = nflds(lines(i))
        if (n .lt. 1) then
          cycle !next keyword
        endif 
c  read keyword
        keywrd = ' '
        pos1 = nxtfld(lines(i), 1)
        pos2 = nxtsp(lines(i), pos1)
        pos3 = nxtfld(lines(i), pos2)
        lenkey = pos2 - pos1
        if (lenkey .gt. 30) then
          errflg = .true.
          write(OUTDOCK, CAUTION1) 'keyword > 30 characters: ', lines(i)
          call doflush(6)
          cycle !next keyword
        endif
        if ((pos3-pos1) .ge. 30) then
          keywrd = lines(i)(pos1:pos1+29)
        else
          keywrd(1:lenkey) = lines(i)(pos1:pos2-1)
        endif
c convert keyword to lower case to simplify recognition
        call tolow(keywrd)

        if (n .le. 1) then
          args = ' '
        else
          args = lines(i)(pos3:132)
        endif
        call dokw(keywrd, args, OUTDOCK, phimap0, recdes0, ligdes0,
     &       options0, spheres0, errflg)
      enddo
      
      return
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
