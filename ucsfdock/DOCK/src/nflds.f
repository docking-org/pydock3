c  --counts number of fields in a line           ECM  5/91
c  --modified to support longer lines            DAG  8/94
c
      integer function nflds(line)

      implicit none
      character*(*) line
      integer l,pos1, pos2

c     functions:
      integer nxtfld, nxtsp

      l = len(line)
      nflds = 0
      pos2 = 1
   10 continue
      pos1 = nxtfld(line, pos2)
      if (pos1.eq.l .and. line(l:l).eq.' ') then
        goto 90
      else
        nflds = nflds + 1
      endif
      pos2 = nxtsp(line, pos1)
      if (pos2.eq.l .and. line(l:l).ne.' ') then
        goto 90
      else
        goto 10
      endif
   90 continue
      return
      end
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
c
