      subroutine readvw(fileunit, options0, MAXTYV, sra, srb)

      use errorformats
      use filenums
      use optionstype

      implicit none

      integer, intent(inout) :: fileunit !file handle
      type(options), intent(in) :: options0
      integer, intent(in) :: MAXTYV
      real, dimension(MAXTYV), intent(inout) :: sra, srb

      character (len=80) line
c       line: temporary variable used for reading

      integer nvtyp
c       nvtyp: number of entries in van der Waals parameter file

      open(unit=fileunit, file=options0%vdwfil, status='old', 
     &    action='read')
      nvtyp = 0
 60   read(fileunit, '(a80)', end=70) line
        if (line(1:1) .eq. '!') goto 60
        nvtyp = nvtyp + 1
        if (nvtyp .gt. maxtyv) then
          write(OUTDOCK, HALT1)
     &         'maximum number of vdw types exceeded'
          write(OUTDOCK, HALT1)
     &         '   recompile with larger maxtyv'
          write(OUTDOCK, HALT0)
          stop
        endif
        read(line, '(10x,f8.2,5x,f8.2)') sra(nvtyp), srb(nvtyp)
      goto 60
 70   continue
      close(fileunit)

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

