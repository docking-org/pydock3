c Open (some) files, mostly the non-grid files.
c

      subroutine open_files(options0, spheres0, MAXTYV, sra, srb, fdlig)

      use errorformats
      use optionstype
      use filenums
      use spheres

      implicit none

      type(options), intent(inout) :: options0
      type(spherest), intent(inout) :: spheres0
      integer, intent(in) :: MAXTYV
      real, dimension(MAXTYV), intent(inout) :: sra, srb
      integer (kind=8), intent(inout) :: fdlig !ligand file handle

      integer pos1, pos2, pos3, lenwrd
      integer nflds, nxtfld, nxtsp
      character (len=80) line
      integer istat,istat2,istat3

      call readvw(VDWTYPESFILE, options0, MAXTYV, sra, srb)

      !write(OUTDOCK, *) "readin:", TRIM(options0%ligfil)
      if (.not. options0%sdi) then
        if (options0%ligfil == "-") then
          fdlig = 0
        endif
      endif

      call gzopen(fdlig, 'r', options0%ligfil, istat)
      if (istat .ne. 0) then !problem opening ligand file
        goto 980
      endif

c     sphere cluster file
      open(MATCHSPH, file=options0%clufil, status='old', err=1000, 
     &    action='read')
      read(MATCHSPH, '(a80)') line
      if (line(1:4) .eq. 'DOCK') then
        call read_color_table(spheres0)
        write(OUTDOCK, 622) spheres0%nreccl
  622   format(i5,' colors used in receptor sphere file')
      else
        backspace(1)   ! save cluster header for later
      endif
      return

  980 write(OUTDOCK, HALT1) 'Error opening ligand atom file'
      write(OUTDOCK, HALT0)
      stop

 1000 write(OUTDOCK, HALT1) 'Error opening sphere cluster file'
      write(OUTDOCK, HALT0)
      stop

      end

c     this is for readin in mol2 (in future, solvation file) to calculate single point energy
      subroutine open_files_rescore(options0, MAXTYV, sra, srb,
     &   fdlig,fdsol,fdvdw)

      use errorformats
      use optionstype
      use filenums

      implicit none

      type(options), intent(inout) :: options0
      integer, intent(in) :: MAXTYV
      real, dimension(MAXTYV), intent(inout) :: sra, srb
      integer (kind=8), intent(inout) :: fdlig !ligand file handle mol2
      integer (kind=8), intent(inout) :: fdsol !ligand file handle desolvation
      integer (kind=8), intent(inout) :: fdvdw !ligand file handle vdw

      integer istat,istat2,istat3

      write(OUTDOCK, *) fdlig, fdsol, fdvdw

      call readvw(VDWTYPESFILE, options0, MAXTYV, sra, srb)

      write(OUTDOCK, *) "file name = ", trim(options0%mol2inputfile) 
      call gzopen(fdlig, 'r', trim(options0%mol2inputfile), istat)
      write(OUTDOCK,*) fdlig, istat

      write(OUTDOCK, *) "file name = ", trim(options0%ligsolinputfile) 
      call gzopen(fdsol, 'r', trim(options0%ligsolinputfile), istat2)
      write(OUTDOCK,*) fdsol, istat2

      write(OUTDOCK, *) "file name = ", trim(options0%ligvdwinputfile) 
      call gzopen(fdvdw, 'r', trim(options0%ligvdwinputfile), istat3)
      write(OUTDOCK,*) fdvdw, istat3
      !call gzopen(fdlig, 'r', 'test.mol2.gz', istat)
      
      if (istat .ne. 0) then !problem opening ligand file
        goto 980
      endif
      if (istat2 .ne. 0) then !problem opening ligand file
        goto 980
      endif
      if (istat3 .ne. 0) then !problem opening ligand file
        goto 980
      endif


      return

  980 write(OUTDOCK, HALT1) 'Error opening ligand mol2 file'
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
