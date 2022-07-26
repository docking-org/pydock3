c----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c                             PROGRAM DOCK
c
c       Copyright (C) 1991,2010-12 Regents of the University of California
c                         All Rights Reserved.
c
c       Copyright (C) 2001 Northwestern University Medical School
c                         All Rights Reserved.
c
c     This program implements the docking algorithm of I.D. Kuntz,
c     J Mol Biol 161, 269-288, 1982. 
c
c     Versions 1.0 and 1.1 of the dock program were based extensively on 
c     the work of Robert Sheridan, Renee DesJarlais, and Tack Kuntz.  
c
c     Version 2.0, macrodock, is based on the work of Brian Shoichet
c     and Tack Kuntz, with help from Dale Bodian. 
c
c     Version 3.0, chemdock, is based on the work of Elaine Meng, Brian
c     Shoichet, and Tack Kuntz.
c
c     Version 3.5 was produced from version 3.0 by Mike Connolly and 
c     updated by Daniel A. Gschwend
c
c     Version 3.5.54 evolved from the efforts of David M. Lorber and
c     John J. Irwin and Brian K. Shoichet
c
c     Version 3.6 re-emerged from the ashes due to renewed efforts by
c     Michael M. Mysinger, Ryan G. Coleman and Michael Carchia
c    
c     Version 3.7 Ryan G. Coleman, Michael Carchia & Teague Sterling
c
c     Version 3.7.1 , Trent Balius, Nir London, Ryan G. Coleman, Michael Carchia & Teague Sterling
c
c     Version 3.8. , Jiankun Lyu, Ben Tingle, Reed Stein, Trent E. Balius ...
c-----------------------------------------------------------------------
c

      program dock

      use grid_reader
      use cov_search
      use cluster_search
      !use rig_frag_search
      !se rig_frag_sample
      use score_only_search
      use search
      use search_multi_cluster
      use rescore ! single point calculation
      use errorformats
      use phimaptype
      use solvmaptype
      use vdwmaptype
      use gisttype
      use optionstype
      use gridstype
      use filenums
      use spheres
      use version

      implicit none

      character (len=8) the_time
      character (len=8) the_date
      character (len=80) the_host
      integer index !counter variable
      integer istat !reopen fdlig
c     these variables (argument & fname) are for dynamic INDOCK filename
      character (len=80) argument, fname
      real stime, ftime
      type(flexgrids) :: grids0
      character (len=80) :: phifil
      character (len=80) :: gistfil
      character (len=80) :: gridbmp
      character (len=80) :: gridvdw
      character (len=80) :: solnam
      character (len=80) :: hsolnam
      type(phimap) :: phimap0, phimap_check, phimap_save
      type(phimap) :: recdes0, recdes_check, recdes_save
      type(phimap) :: ligdes0, ligdes_check, ligdes_save !new grids, stop
          !using GB code for ligand desolvation, 
          !use PB code like everything else in DOCK.      
      type(vdwmap) :: vdwmap0, vdwmap_check
      type(solvmap) :: solvmap0, solvmap_check
      type(gistparm) :: gist0, gist_check, rec_des_r0, rec_des_rcheck
      type(options) :: options0
      integer :: which_grid !which grid to process (for reading)
      logical :: sizescheck !if the sizes of multiple grids are the same
      real, dimension(3) :: allminlimits, allmaxlimits ! for outside_grid checking
      real, dimension(3) :: tempminlimits, tempmaxlimits ! for outside_grid checking
      integer, parameter :: MAXTYV = 26 !vdw types, only ever 26
      integer, parameter :: MAXOR = 1000000 !orientations, should be dynamic 
      integer, parameter :: MAXCTR = 4 !how many centers match, 4 is fine
      real, dimension(MAXTYV) :: sra, srb !square roots of vdw parms
      type(spherest) :: spheres0
      integer, parameter :: MAXLNS = 3000 !max lines in INDOCK file
      character (len=300), dimension(MAXLNS) :: lines !where INDOCK is stored
      logical :: errflg !problems while reading INDOCK
      integer (kind=8) :: fdlig, fdsol, fdvdw !ligand file handle
      integer :: sph_idx 
      type(spherest), dimension(10) :: sphere_sets ! array to hold sph clusters
      INTEGER, PARAMETER :: SEEK_SET = 0, SEEK_CUR = 1, SEEK_END = 2
      INTEGER :: fd, offset, ierr

      ierr   = 0
      offset = 5
      fd     = 10

      fdlig = 100
      fdsol = 101
      fdvdw = 102

      call get_time(stime)
! sets defaults for options set in INDOCK
      call setdef(phimap0, options0, spheres0, recdes0, ligdes0) 

      open(unit=OUTDOCK, file='OUTDOCK', action='write')
      !open(unit=OUTDOCK_premin, file='OUTDOCK_premin', action='write')

c     initialize input error flag
      errflg = .false.

c     INDOCK read is directed in readfl to readkw, dokw, and ckparm
      call run_info(the_date, the_time, the_host)
      write(OUTDOCK, '(a,a)') 'DOCK v. 3.8.0, compiled ', COMPDATE
      write(OUTDOCK, '(a,a)') 'DOCK compiled for bitsize ', ARCHSIZE
      write(OUTDOCK, '(a,a)') 'DOCK subversion revision: ', SVNVERSION
      write(OUTDOCK, '(a,a,1x,a,1x,a)') 'CPU, Date, and Time: ',
     &  trim(the_host), trim(the_date), trim(the_time)
      call doflush(OUTDOCK)

      call get_indock(1, argument) !the 1 here means get the first cli argument
      read(argument, '(a)') fname
      index = len_trim(fname)
      if (index .eq. 0) then !use default INDOCK name
        fname = 'INDOCK'
      endif
!     read and check parameters in keyword format INDOCK
      call readfl(INDOCK, fname, OUTDOCK, phimap0, recdes0, ligdes0, 
     &    options0, spheres0, MAXLNS, lines, errflg)
      call ckparm(OUTDOCK, phimap0, recdes0, ligdes0, 
     &    gist0, options0, vdwmap0, rec_des_r0,
     &    MAXOR, MAXCTR, errflg)
      call allocate_flexgrids(grids0, options0%total_receptors, OUTDOCK)
      call setup_grouppart(grids0, options0%group_numbers, 
     &    options0%part_numbers, options0%total_receptors, OUTDOCK)
      call doflush(OUTDOCK) !flush before reading, in case something bad happens
      !copy the parameters read in to a check copy, for checking flex rec
      vdwmap_check = vdwmap0
      phimap_check = phimap0
      recdes_check = recdes0
      ligdes_check = ligdes0
      gist_check = gist0
      rec_des_rcheck = rec_des_r0
      phimap_save = phimap0
      recdes_save = recdes0
      ligdes_save = ligdes0
      do which_grid = 1, options0%total_receptors
        phimap_check = phimap_save !have to do this to allow phigrids to be sometimes implicit 0 (as long as not the first ones)
        recdes_check = recdes_save !have to do this to allow the sometimes implicit0 grids
        call read_grids(grids0%griddata(which_grid),
     &      OUTDOCK, tempminlimits, tempmaxlimits, 
     &      phimap_check, recdes_check, solvmap_check, vdwmap_check,
     &      ligdes_check, gist_check, rec_des_rcheck, 
     &      options0, which_grid)
        if (1 .eq. which_grid) then !copy into 0, which is the one we keep
          phimap0 = phimap_check
          gist0 = gist_check
          rec_des_r0 = rec_des_rcheck
          recdes0 = recdes_check
          ligdes0 = ligdes_check
          solvmap0 = solvmap_check
          vdwmap0 = vdwmap_check
          allminlimits = tempminlimits
          allmaxlimits = tempmaxlimits
          write(OUTDOCK, *) 'min grid limits:', allminlimits
          write(OUTDOCK, *) 'max grid limits:', allmaxlimits
          write(OUTDOCK, *) 'first set of grids read successfully'
          call doflush(OUTDOCK)
        else !we need to check to make sure the grids are all the same size
          ! incase the frist is not defined. 
          gist0 = gist_check
          rec_des_r0 = rec_des_rcheck
          ! we should add gist to this check. 
          sizescheck = gridsize_check(vdwmap0, vdwmap_check, 
     &        phimap0, phimap_check, recdes0, recdes_check, 
     &        solvmap0, solvmap_check, ligdes0, ligdes_check)
          write(OUTDOCK, *) "grid sizes check", which_grid, sizescheck
          if (.not. sizescheck) then
            write(OUTDOCK, HALT1) "grid sizes not equal"
            write(OUTDOCK, HALT0) 
            stop
          endif
        endif
      enddo

!     see if there was an error and act appropriately
      if (errflg) then
        write(OUTDOCK, HALT1) 'Could not parse all INDOCK parameters!'
        write(OUTDOCK, HALT0) 
        stop
      endif
      if (options0%search_flag) then
          call open_files(options0, spheres0, MAXTYV, sra, srb, fdlig)
    
!         set up color table if chemical matching
          if (options0%cmatch) then
            call clsetp(options0, spheres0)
          endif
    
!         read spheres for this sphere cluster
          if (options0%docking_mode .eq. 4) then
              do sph_idx = 1, options0%k_clusters
                  sphere_sets(sph_idx) = spheres0
                  call clread(MATCHSPH, options0, 
     &              sphere_sets(sph_idx), sph_idx)
              enddo
          else
              call clread(MATCHSPH, options0, spheres0, 1)
!             call clread(MATCHSPH, options0, spheres0)
          endif
          close(MATCHSPH) !close the file before proceeding
          call doflush(OUTDOCK)
          !run_search does all ligand reading, scoring and writing

          !check if this is a covalent run
          write(OUTDOCK, *) 'test',options0%docking_mode 
          if (options0%dockovalent .eqv. .true.) then
             write(OUTDOCK, *) 'DOCKovalent'
             call run_cov_search(grids0, phimap0, recdes0, ligdes0, 
     &            gist0, solvmap0, rec_des_r0, vdwmap0, options0, 
     &            allminlimits, allmaxlimits, spheres0, MAXTYV, sra, 
     &            srb, MAXOR, MAXCTR, fdlig)
          else                  ! a regular (non-covalent) run
             if (options0%docking_mode .eq. 0) then
                 write(OUTDOCK, *) 'ClusteringDOCK'
                 call run_cluster_search(grids0, phimap0, recdes0, 
     &                ligdes0, gist0, solvmap0, rec_des_r0, vdwmap0, 
     &                options0, allminlimits, allmaxlimits, spheres0, 
     &                MAXTYV, sra, srb, MAXOR, MAXCTR, fdlig)
c             else if (options0%docking_mode .eq. 1) then
c                 write(OUTDOCK, *) 'rigid fragment sampling'
c                 call run_rig_frag_sample(grids0, phimap0, recdes0, 
c     &                ligdes0, gist0, solvmap0, rec_des_r0, vdwmap0, 
c     &                options0, allminlimits, allmaxlimits, spheres0, 
c     &                MAXTYV, sra, srb, MAXOR, MAXCTR, fdlig)
             else if (options0%docking_mode .eq. 2) then
                 write(OUTDOCK, *) 'read in orientations then score'
                 call run_score_only_search(grids0, phimap0, recdes0,
     &                ligdes0, gist0, solvmap0, rec_des_r0, vdwmap0, 
     &                options0, allminlimits, allmaxlimits, spheres0, 
     &                MAXTYV, sra, srb, MAXOR, MAXCTR, fdlig)
c             else if (options0%docking_mode .eq. 3) then
c                 write(OUTDOCK, *) 'rigid fragment docking'
c                 call run_rig_frag_search(grids0, phimap0, recdes0,
c     &                ligdes0, gist0, solvmap0, rec_des_r0, vdwmap0, 
c     &                options0, allminlimits, allmaxlimits, spheres0, 
c     &                MAXTYV, sra, srb, MAXOR, MAXCTR, fdlig)
             else if (options0%docking_mode .eq. 4) then
                 write(OUTDOCK, *) 'subcluster matching'
                 !run_search does all ligand reading, scoring and writing
                 call run_search_multi_cluster(grids0, phimap0, 
     &            recdes0, ligdes0,
     &            gist0, solvmap0,rec_des_r0, vdwmap0, options0,
     &            allminlimits, allmaxlimits, sphere_sets,
     &            MAXTYV, sra, srb, MAXOR, MAXCTR, fdlig)
c                 do sph_idx = 1, options0%k_clusters
c                   
c                   if (sph_idx > 1) then
c                     ! reopen SDIFILE
c                     open(unit=SDIFILE, file='split_database_index',
c     &                 status='old', action='read')
c                     ! reload next filename
c                     read(SDIFILE,'(a255)') options0%ligfil
c                     ! open next filename
c                     call gzopen(fdlig, 'r', options0%ligfil, 0)
c                     write(OUTDOCK, *) "open the file: ", 
c     &                 TRIM(options0%ligfil)                 
c                   endif
c
c                   call run_search_multi_cluster(grids0, phimap0, 
c     &              recdes0, ligdes0,
c     &              gist0, solvmap0,rec_des_r0, vdwmap0, options0,
c     &              allminlimits, allmaxlimits, sphere_sets(sph_idx),
c     &              MAXTYV, sra, srb, MAXOR, MAXCTR, fdlig, sph_idx)
c                 enddo
             else
                 !run_search does all ligand reading, scoring and writing
                 call run_search(grids0, phimap0, recdes0, ligdes0, 
     &                gist0, solvmap0,rec_des_r0, vdwmap0, options0, 
     &                allminlimits, allmaxlimits, spheres0, 
     &                MAXTYV, sra, srb, MAXOR, MAXCTR,fdlig)
             endif
          endif
      !this is a rescoring run
      else if (options0%rescore_flag) then
         call open_files_rescore(options0, MAXTYV, sra, srb,
     &          fdlig, fdsol, fdvdw  )
         !write(6,*) "I AM HERE (##)"
         call run_rescore(grids0, phimap0, recdes0, ligdes0, gist0,
     &        solvmap0, rec_des_r0, vdwmap0, options0, allminlimits, 
     &        allmaxlimits, MAXTYV, sra, srb, fdlig, fdsol, fdvdw)
c    &        MAXTYV, sra, srb, MAXOR, MAXCTR, fdlig)
         !write(6,*) "I AM HERE (###)"
      else 
         write(6,*) "Error not a valid search_type options. 
     &              (2 = rescore;1= standard)"
         Stop
      endif

      call get_time(ftime) 
      call run_info(the_date, the_time, the_host) !write cleanup stuff out
      write(OUTDOCK, '(a, 1x, a, 1x, a)') 'Date and Time:', 
     &   trim(the_date), 
     &   trim(the_time)
      write(OUTDOCK, '(a19, f14.4, 1x, a8, f11.4)') 
     &   'elapsed time (sec):', ftime - stime, 
     &   ' (hour):', ((ftime - stime)/3600.)
      call doflush(OUTDOCK)
      close(OUTDOCK)
      close(OUTDOCK_premin)
      
      stop 
      end program dock
      
