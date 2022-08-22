c Read color table from ligand database
c
      subroutine read_color_table(spheres0)

      use spheres
      use filenums

      implicit none

      type(spherest), intent(inout) :: spheres0

      character (len=80) line
      integer istat, ligloc
      integer ncol

      ncol = 0
      do while (.true.)
        read(MATCHSPH, '(a80)') line !file already opened... horrible code
        if (len(line) .eq. 0) then
          cycle      ! skip blank lines
        endif
        if (line(1:7) .eq. 'cluster') then
          backspace(1) ! backspace over the cluster header line
          spheres0%nreccl = ncol
          exit !break out of this loop
        endif
        ncol = ncol + 1
        read(line, '(a30)') spheres0%rclnam(ncol)
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
c
