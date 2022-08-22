
      subroutine clread(runit, options0, spheres0, cluster_idx)

      use errorformats
      use optionstype
      use filenums
      use spheres

      implicit none

      integer runit !file handle, already opened
      type(options), intent(in) :: options0
      type(spherest), intent(inout) :: spheres0
      integer cluster_idx !which cluster index to read into spheres0

      integer i, j
      integer iclus(MAXPTS)
c     iclus:  array of spheres in the cluster of interest.
c     idone:  flag for whether the cluster of interest has been read
      integer idone
c     nctemp:  temporary variable for read in of cluster number
      integer nctemp
c     nmax:  temporary variable for number of spheres in cluster
      integer nmax
      character (len=80) line
c     extraLines - whether extraneous lines exist atop cluster file
      logical extralines
      integer critsph !not used anymore, but read in and discarded

      extralines = .false.

      idone = 0
      do while (idone .eq. 0)
        read(runit,'(a80)', end=200) line
        if (line(1:7) .ne. 'cluster') then
          if (.not. extralines) then
            write(OUTDOCK, CAUTION1) 'extraneous lines detected at',
     &        ' top of cluster file'
          endif
          extralines = .true.
          cycle
        endif


        read(line, 51, end=200) nctemp, nmax
 51     format(8x,i5,32x,i5,2x,5i5)
        write(OUTDOCK,'(a,6i4)') 'Receptor spheres: ',nmax
        
c       check for pseudo-cluster (sphere dump) at end of file
        if (nctemp .eq. 0) goto 200
        if (nctemp .eq. cluster_idx) then
          spheres0%nsphr = nmax
          if (spheres0%nsphr .gt. maxpts) then
            write(OUTDOCK, HALT5)
     &           'array maxpts ',maxpts,'  receptor requires ',
     &           spheres0%nsphr
            write(OUTDOCK, HALT0)
            stop
          endif
c         sphere radius and j-atom number not read
c         flagged sphere is interpreted as a differential docking sphere
c         long words, such as 'flagged' are no longer permitted after
c         coordinates, in order to make room for the sphere color
          do i = 1, spheres0%nsphr
            read(runit,57) line
 57         format(a80)
            read(line,58) iclus(i), (spheres0%spcorr(j,i),j=1,3),
     &           critsph, spheres0%scolor(i)
 58         format(i5,3f10.5,13x,i2,i3)
c           convert from sphere cluster file color numbers
c           to merged (ligand+receptor) color numbers
          enddo
          idone = 1
          if (options0%cmatch) then
            call transp(spheres0%nsphr, spheres0%scolor, spheres0)
          endif
        endif
      enddo

c      write(OUTDOCK, 65) clufil
c   65 format(' end of cluster in file: ',a50)
      write(OUTDOCK, 73) cluster_idx, spheres0%nsphr
   73 format(' cluster',i5,' with',i5,' spheres')
      write(OUTDOCK, 74) (iclus(i), i=1, spheres0%nsphr)
   74 format(15i5)
      write(OUTDOCK, *)
      return

c end of file (or sphere dump) reached without success
  200 continue
      write(OUTDOCK, HALT1)
     &     'end of cluster file read, but cluster not found'
      write(OUTDOCK, HALT0)
      stop
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
