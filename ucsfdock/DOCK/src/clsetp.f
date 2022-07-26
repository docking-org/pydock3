c Setup color table
c

      subroutine clsetp(options0, spheres0)

      use errorformats
      use filenums
      use optionstype
      use spheres

      implicit none

      type(options), intent(in) :: options0
      type(spherest), intent(inout) :: spheres0

      integer ligtrn(MAXCOL)
c ligtrn: ligand translation tables (for color numbers)
      integer l, r, m
      integer narray

      integer lkorad, lkpcol !functions

      !this code sets the default colors, which as far as i can tell, are 
      !the only ones ever used. later plan is to modify this to read in from 
      !ligand file if there is anything there (starting with T for ligand type)
      spheres0%nligcl = 7
c     modified by jklyu, 20191216
      spheres0%lclnam(1) = '1' !'positive'
      spheres0%lclnam(2) = '2' !'negative'
      spheres0%lclnam(3) = '3' !'acceptor'
      spheres0%lclnam(4) = '4' !'donor'
      spheres0%lclnam(5) = '5' !'ester_o'
      spheres0%lclnam(6) = '6' !'amide_o'
      spheres0%lclnam(7) = '7' !'neutral'

      narray = 0
      if (spheres0%nligcl .le. 0) then
        write(OUTDOCK, HALT1) 'clsetp:  no ligand atom colors'
        write(OUTDOCK, HALT0)
c       added this error condition back (needs testing) 
        stop
      endif
      if (spheres0%nreccl .le. 0) then
        write(OUTDOCK, HALT1) 'clsetp:  no receptor sphere colors'
        write(OUTDOCK, HALT0)
        stop
      endif
      if (spheres0%nlrmat .le. 0) then
        write(OUTDOCK, HALT1)
     &       'clsetp:  no ligand-receptor color matches specified'
        write(OUTDOCK, HALT0)
        stop
      endif

c add ligand colors to color table
      do l = 1, MAXCOL
        ligtrn(l) = 0
      enddo
      do l = 1, spheres0%nligcl
        ligtrn(l) = lkorad(spheres0%lclnam(l), options0, spheres0)
        if (ligtrn(l) .le. 0) then
          write(OUTDOCK, HALT2) 
     &         'clsetp:  problem merging ligand color: ',
     &         spheres0%lclnam(l)
          write(OUTDOCK, HALT0)
        endif
        write(OUTDOCK, 295) spheres0%lclnam(l), l, ligtrn(l)
c       extra spaces below to match receptor width
  295   format(' ligand   color ',a30,2x,i5,' merged number = ',i5)
      enddo

c     add receptor colors to color table
      do r = 1, MAXCOL
        spheres0%rectrn(r) = 0
      enddo
      do r = 1, spheres0%nreccl
        spheres0%rectrn(r) = lkorad(spheres0%rclnam(r), options0, 
     &      spheres0)
        if (spheres0%rectrn(r) .le. 0) then
          write(OUTDOCK, HALT2) 
     &          'clsetp:  problem merging ligand color: ',
     &                              spheres0%rclnam(r)
          write(OUTDOCK, HALT0)
          stop
        endif
        write(OUTDOCK, 395) spheres0%rclnam(r), r, spheres0%rectrn(r)
  395   format(' receptor color ',a30,2x,i5,' merged number = ',i5)
      enddo
      write(OUTDOCK, *)
      write(OUTDOCK, 460) spheres0%numcol
  460 format(i5,' colors in merged ligand-sphere color table')

c set up color matching table
      do l = 1, MAXCOL
        do r = 1, MAXCOL
          spheres0%lgspcl(l, r) = 0
        enddo
      enddo
      write(OUTDOCK, *) 'this is a run with labeled color matching:'
      do m = 1, spheres0%nlrmat
        l = lkpcol(spheres0%premat(1,m), options0, spheres0)
        if (l .le. 0) then
          write(OUTDOCK, 
     &         HALT2) 'undefined ligand color used in matching: ',
     &         spheres0%premat(1,m)
          write(OUTDOCK, HALT0)
          stop
        endif
        r = lkpcol(spheres0%premat(2,m), options0, spheres0)
        if (r .le. 0) then
          write(OUTDOCK, HALT2) 
     &         'undefined receptor color used in matching: ',
     &         spheres0%premat(2,m)
          write(OUTDOCK, HALT0)
          stop
        endif
c       set matching array element to 1
        spheres0%lgspcl(l, r) = 1
        narray = narray + 1
        write(OUTDOCK, 585) spheres0%colnam(l), spheres0%colnam(r)
  585   format('         match ',a30,2x,a30)
      enddo

      write(OUTDOCK, 740) narray
  740 format(i5,' chemical matches specified')
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
