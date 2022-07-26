c jklyu, 20200223
c read one line from rig_frag_index file
c get ligfil, ligname, pos and rig_frag_name
      subroutine readsdiln(runit, OUTDOCK, options0, 
     &        inputstatus, errflg)

      use errorformats
      use optionstype

      implicit none

      integer runit !file handle
      integer, intent(in) :: OUTDOCK
      type(options), intent(inout) :: options0
      logical, intent(inout) :: errflg  !problems?
      integer, intent(inout) ::  inputstatus

      integer i, n !replace these with meaningful names
      integer pos1,pos2,pos3,pos4,pos5,pos6,pos7,pos8
      character (len=300) line

c     functions:
      integer nflds, nxtfld, nxtsp

c     file is already open, read one line
      read(runit, '(a255)', iostat=inputstatus) line
      if (inputstatus .gt. 0) then
        errflg = .true.
        write(OUTDOCK,*) 'errors in reading rig_frag_dock_index'
        return
      else if (inputstatus .eq. -1) then
        write(OUTDOCK,*) 'we reached the end of the rig_frag_dock_index'
        return
      endif

      n = nflds(line)
      if (n .ne. 4) then
        errflg = .true.
        write(OUTDOCK,*) 'sdi does not have 4 column'
        call doflush(6)
        return
      else
        pos1 = nxtfld(line, 1)
        pos2 = nxtsp(line, pos1)
        options0%ligfil = TRIM(line(pos1:pos2-1))
        write(OUTDOCK,*) options0%ligfil
        pos3 = nxtfld(line, pos2)
        pos4 = nxtsp(line, pos3)
        options0%ligname = TRIM(line(pos3:pos4-1))
        write(OUTDOCK,*) options0%ligname
        pos5 = nxtfld(line, pos4)
        pos6 = nxtsp(line, pos5)
c        options0%pos = ichar(line(pos5:pos6-1))
        read(line(pos5:pos6-1), *) options0%pos
        write(OUTDOCK,*) options0%pos
        pos7 = nxtfld(line, pos6)
        pos8 = nxtsp(line, pos7)
        options0%rig_frag_name = TRIM(line(pos7:pos8-1))
        write(OUTDOCK,*) options0%rig_frag_name
      endif
      
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
