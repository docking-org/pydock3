c opens indock file?
      subroutine readln(runit, fname, fail, nlines, MAXLNS, lines)

      use errorformats
      use filenums

      implicit none

      integer runit !file handle
      logical fail
      character*(*) fname !what is this??
      integer nlines
      integer, intent(in) :: MAXLNS
      character (len=300), dimension(MAXLNS), intent(inout) :: lines

c     read  file, if present
      open(unit=runit, file=fname, status='old', err=80, 
     &    action='read')
      do while (.true.)
        if (nlines .ge. MAXLNS) then
          write(OUTDOCK, HALT1) 'Too many lines in parameter file'
          write(OUTDOCK, HALT4) 
     &         '  Recompile with higher ',MAXLNS,' in dock.h'
          write(OUTDOCK, HALT0)
          close(runit)
          stop
        endif
        nlines = nlines + 1
        read(runit, '(a132)', end=20) lines(nlines)
      enddo

   20 continue
      nlines = nlines - 1
      close(runit)
      fail = .false.
      return

   80 continue
      fail = .true.
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
