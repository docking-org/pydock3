      subroutine getvrn(fline, filvrn, MAXLNS, lines)

      use filenums

      implicit none

      integer fline
      character (len=30) filvrn
      integer, intent(in) :: MAXLNS
      character (len=300), dimension(MAXLNS), intent(inout) :: lines

      integer n
      integer pos1, pos2, pos3, lenwrd
      character (len=30) word1, word2, word3

c  --functions:
      integer nflds, nxtfld, nxtsp

      n = nflds(lines(fline))
      if (n .lt. 3) then
        write(OUTDOCK, *) 'No INDOCK version number supplied'
        filvrn = '3.7'
        return
      endif

c  read three words
      word1 = ' '
      pos1 = nxtfld(lines(fline), 1)
      pos2 = nxtsp(lines(fline), pos1)
      pos3 = nxtfld(lines(fline), pos2)
      if ((pos3-pos1) .ge. 30) then
        word1 = lines(fline)(pos1:pos1+29)
      else
        lenwrd = pos2 - pos1
        word1(1:lenwrd) = lines(fline)(pos1:pos2-1)
      endif
c
      word2 = ' '
      pos1 = nxtfld(lines(fline), pos2)
      pos2 = nxtsp(lines(fline), pos1)
      pos3 = nxtfld(lines(fline), pos2)
      if ((pos3-pos1) .ge. 30) then
        word2 = lines(fline)(pos1:pos1+29)
      else
        lenwrd = pos2 - pos1
        word2(1:lenwrd) = lines(fline)(pos1:pos2-1)
      endif
c
      word3 = ' '
      pos1 = nxtfld(lines(fline), pos2)
      pos2 = nxtsp(lines(fline), pos1)
      pos3 = nxtfld(lines(fline), pos2)
      if ((pos3-pos1) .ge. 30) then
        word3 = lines(fline)(pos1:pos1+29)
      else
        lenwrd = pos2 - pos1
        word3(1:lenwrd) = lines(fline)(pos1:pos2-1)
      endif

      if (word1 .ne. 'DOCK' .and. word1 .ne. 'dock') then
        filvrn = '3.5'
        return
      endif

      filvrn = word2
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
c

