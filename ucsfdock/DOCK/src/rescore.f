
c-----------------------------------------------------------------------
c     This subroutine rescores a mol2 file of "ligands".
c     writen by Trent E Balius, 2014
c     modification of subroutine search
c     Jiankun modifed to work beter with the minimizer. 
c     TODO:
c     (1) We need to get flexible receptor working with the minimizer.
c        (1.1) to do this we need to restructure the code so that it
c              score molecules on the combinations of the grids not the
c              indevidual grids as it does now.
c        (1.2) we need to work out how we want to do minimization. show
c              the pose be minimized on everycombination or just the
c              best scoring one? 

c-----------------------------------------------------------------------
      module rescore

      implicit none
      contains


c         call run_rescore(grids0, phimap0, recdes0, ligdes0, gist0,
c    &         solvmap0,vdwmap0, options0, allminlimits, allmaxlimits,
c    &         MAXTYV, sra, srb, MAXOR, MAXCTR, fdlig)

      subroutine run_rescore(grids0, phimap0, recdes0, ligdes0, gist0,
     &    solvmap0,rec_des_r0, vdwmap0, options0, allminlimits, 
     &    allmaxlimits,  MAXTYV, sra, srb, fdlig, fdsol, fdvdw)
c    &    MAXTYV, sra, srb, MAXOR, MAXCTR, fdlig)

      use phimaptype
      use gisttype
      use solvmaptype
      use vdwmaptype
      use status
      use optionstype
      use gridstype
      use db2type
      use ligscoretype
      use filenums
      use score_mol
      use atomscoretype
      use matchtype
      use score_conf
      use check_conf_og

c these are the dynamically allocated grids
      type(flexgrids), intent(inout) :: grids0
c these are user-defined types that contain information about grids (usually)
      type(phimap), intent(inout) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(inout) :: gist0, rec_des_r0
      type(solvmap), intent(inout) :: solvmap0
      type(vdwmap), intent(inout) :: vdwmap0
c options
      type(options), intent(in) :: options0
c for outside grid checking
      real, dimension(3) :: allminlimits, allmaxlimits
c vdw parameters
      integer, intent(in) :: MAXTYV
      real, dimension(MAXTYV), intent(in) :: sra, srb !vdw parms
c     integer, intent(in) :: MAXOR
      integer             :: MAXOR
c     integer, intent(in) :: MAXCTR
      integer (kind=8), intent(inout) :: fdlig, fdsol, fdvdw !ligand file handle

c     variables --
      type(db2) :: db2lig !ligand data stored here
      type(ligscoret) :: ligscore !saves many top poses
      type(ligscoret) :: ligscoreeach !saves one top pose for each flexible receptor combination
      type(atomscoret) :: atomscore !recalculate per atom scores
      type(matcht) :: match !stuff about matching
      character (len=80) molout !filename to save
      logical iend
c        iend:  indicator variable for end of ligand file.
c               is equal to zero if end not reached.
      integer molct
c        molct:  the number of ligands searched on this run.
      integer count !temp counter, don't use i,j,n or other one letter nonsense
c best scoring conformation for a given ligand
      integer flush_int
c        # of compounds interval to flush buffer to OUTDOCK
c     real tcorr(3, MAXPTS) !used to hold the rigid component heavy atoms
c     integer tcolor(MAXPTS)
      real telaps, time2, time1 !times
      logical exists
      integer dock_status !uses status.h enums
      integer istat !ligand file status
      integer tempsave, temprank !counter for multiple output
      integer cloudcount, cloudtotal !counter and overall # of times to run
      logical doing_cloud !whether or not the cloud is actually happening
      integer setstart, setend !controls which sets (complete confs) to score
      integer debug1 !temporary variable used in debugging loops
      real curdislim !temporary variable, keeps track of current distance limit
c spheres, one for rec, one lig. the count is how many pairs have been added.
      !integer hash(MAXOR, MAXCTR), hashcount
      integer (kind=8) :: fdmol !output mol2.gz file handle
      integer :: outcount

      real    vdwasum, vdwbsum, eleceel, gistval, solva, solvp, solvd,
     &        rdscore, hescore 
      integer numscored, conf, which_grid
      integer total_iter 



c format lines for output
      character (len=*), parameter :: molexam = 
     &     '(i8,"  molecules examined on this run")'
      character (len=*), parameter :: outdockerr = 
     &     '(i6,1x,a16,i8,1x,i10,1x,f7.2,1x,a)'
      character (len=*), parameter :: outdockline = '(i6,1x,a16,1x,a16,
     &       1x,i8,1x,i10,1x,f7.2,1x,i3,1x,i9,1x,i9,1x,i6,3x,i3,2x,
     &       f7.2,1x,f7.2,1x,f7.2,f7.2,1x,f7.2,1x,f7.2,1x,f7.2,1x,f7.2,
     &       1x,f7.2,f10.2)'

      logical temp_status !temporary status holder
c     --- initialize variables ---
      integer maxpts 

c     these are not needed, only as pace holders
c     integer,dimension(10000,options0%total_receptors) :: gistgridindex
      integer,dimension(100000) :: gistgridindex
      integer,dimension(options0%total_receptors) :: lengistgi ! this is the curent length of list grid point already asigned to the molecule

      if (options0%mol2_minimizer .ne. 0 .and.
     &   (options0%flexible_receptor)) then
          write(OUTDOCK,*) "minimization and flexible receptor "
          write(OUTDOCK,*) "rescoring are not compatable"
          stop
      endif
      if (options0%flexible_receptor .and.
     &    options0%score_each_flex) then
        !call allocate_ligscore(ligscoreeach, grids0%total_combinations,
        call allocate_ligscore(ligscoreeach, options0%total_receptors,
     &      OUTDOCK)
      endif


      maxpts = 500
      match%nmatch = 1 !initialize various things
      dock_status = RIGIDOKAY !see status.f
      outcount = 1
      cloudcount = 0
      cloudtotal = 0 

      !write(*,*) "MAXOR = ", MAXOR 
      MAXOR = 1
      !write(*,*) "MAXOR = ", MAXOR 
      write(*,*) "gistvolume = ", options0%gistvolume 
      write(*,*) "gist_aprox= ", options0%gist_aprox(1)

      call allocate_match(match, 4, OUTDOCK)
      ! we should initialize stuff here. 
      do count = 1,4
         !row 1 
         match%rot(1,1,count) = 1.0 ! diganal
         match%rot(1,2,count) = 0.0
         match%rot(1,3,count) = 0.0
         ! row 2
         match%rot(2,1,count) = 0.0
         match%rot(2,2,count) = 1.0 ! diganal
         match%rot(2,3,count) = 0.0
         ! row 3
         match%rot(3,1,count) = 0.0
         match%rot(3,2,count) = 0.0
         match%rot(3,3,count) = 1.0 ! diganal
         ! translation ? 
         match%coml(1,count) = 0.0
         match%coml(2,count) = 0.0
         match%coml(3,count) = 0.0
         
         match%comr(1,count) = 0.0
         match%comr(2,count) = 0.0
         match%comr(3,count) = 0.0
      enddo

      !call allocate_ligscore(ligscore, options0%nsav, OUTDOCK)
      call allocate_atomscore(atomscore, maxpts, OUTDOCK)
      if (options0%flexible_receptor .and. 
     &    options0%score_each_flex) then
        call doflush(6)
c       call allocate_ligscore(ligscore, grids0%total_combinations, 
c    &      OUTDOCK)
        call allocate_ligscore(ligscore, options0%total_receptors, 
     &      OUTDOCK)
      else
        call allocate_ligscore(ligscore, options0%nsav, OUTDOCK)
      endif
      molout = trim(options0%outfil)//'mol2.gz'
      call gzopen(fdmol, 'w', TRIM(molout), istat) !open output ligand file
      !top of main section of OUTDOCK file
      write(OUTDOCK, *)  options0%rdscale, options0%rec_des_scale
      write(OUTDOCK, '(a,a,a,a)') 
     & '  mol#           id_num     flexiblecode  matched    nscored  ',
     & '  time hac    setnum    matnum   rank cloud    elect +  gist',
     & ' +   vdW + psol +  asol + inter + rec_e + rec_d + r_hyd =',
     & '    Total'
c     write(OUTDOCK,*) OUTDOCK, ' fdmol:', fdmol, ' fdsol:',fdsol,
c    & ' fdvdw:',fdvdw
c     call doflush(OUTDOCK)
      !write(OUTDOCK,*) "I AM HERE (doflush)"
      call doflush(OUTDOCK)
      iend = .false.
      do while (.not. iend) !an actual loop
        total_iter = 0
c       --- read a ligand structure from the dbase ---
        !write(OUTDOCK,*) "I AM HERE (A)"
        call ligread_mol2(db2lig, options0, 
     &        iend, fdlig, fdsol, fdvdw) !read one hierarchy 
c init these here???
        ! initialize rot and trans matrix
        do count = 1,4
           !row 1 
           match%rot(1,1,count) = 1.0 ! diganal
           match%rot(1,2,count) = 0.0
           match%rot(1,3,count) = 0.0
           ! row 2
           match%rot(2,1,count) = 0.0
           match%rot(2,2,count) = 1.0 ! diganal
           match%rot(2,3,count) = 0.0
           ! row 3
           match%rot(3,1,count) = 0.0
           match%rot(3,2,count) = 0.0
           match%rot(3,3,count) = 1.0 ! diganal
           ! translation ? 
           match%coml(1,count) = 0.0
           match%coml(2,count) = 0.0
           match%coml(3,count) = 0.0
           
           match%comr(1,count) = 0.0
           match%comr(2,count) = 0.0
           match%comr(3,count) = 0.0
        enddo
        !write(OUTDOCK,*) "I AM HERE (B)"
        !call doflush(OUTDOCK)
          !stop
          !ligscore%savedcount = 0 !reset the number of poses saved
          ligscore%savedcount = 1 !reset the number of poses saved
          if (options0%flexible_receptor .and. 
     &        options0%score_each_flex) then
            call reset_ligscore(ligscoreeach)
          endif
        !write(OUTDOCK,*) "I AM HERE (C)"
        !call doflush(OUTDOCK)
        setstart = 1
        setend   = db2lig%total_atoms
c         write(OUTDOCK,*) "setstart = ",setstart,
c    &                     "setend   = ",setend
c         call doflush(OUTDOCK)
          
c         call calc_score_mol(dock_status, setstart, setend, 
c    &        db2lig, options0, ligscore, ligscoreeach, grids0,
c    &        vdwmap0, phimap0, recdes0, ligdes0, gist0, solvmap0,
c    &        match, allminlimits, allmaxlimits,
c    &        sra, srb, 
c    &        MAXOR, db2lig%total_confs, db2lig%total_sets, 
c    &        db2lig%total_atoms, db2lig%total_coords,
c    &        MAXTYV)
cc TEB start
c calling the function calc_score_conf durectly here to remove some of
c unneeded lyers  

          conf      = 1
          vdwasum   = 0.0
          vdwbsum   = 0.0
          eleceel   = 0.0
          gistval   = 0.0
          solva     = 0.0
          solvp     = 0.0
          solvd     = 0.0
          hescore   = 0.0
          numscored = 1
c         which_grid = 1
          tempsave = 1
          
          flush_int = options0%input_flush_int !flush after every (this many) ligands

          call get_time(time1)

          do tempsave = 1, options0%total_receptors !go through every grid & 
c         do tempsave = 1, grids0%total_combinations !each receptor combo
c            gistval   = 0.0
             conf      = 1
             vdwasum   = 0.0
             vdwbsum   = 0.0
             eleceel   = 0.0
             gistval   = 0.0
             solva     = 0.0
             solvp     = 0.0
             solvd     = 0.0
             hescore   = 0.0
             numscored = 1

             which_grid = tempsave
             call get_time(time2)
             telaps = time2 - time1
             molct = 1
             !flush_int = 1
             
c            ligscore%pose_reccode(tempsave) = "rescore"
             write(ligscore%pose_reccode(tempsave), 
     &       "(A7,I2.2)") "rescore", tempsave
             ligscore%setsave(tempsave)      = 1
             ligscore%orientsave(tempsave)   = 1
             !write(OUTDOCK,*) "I AM HERE (D2)"
             !call doflush(OUTDOCK)
             temp_status = atom_out_of_grids(conf, db2lig,
     &           allminlimits, allmaxlimits) !do outside grid stuff
             if (temp_status .eqv. .true.) then !outside grid, quit now
               !confstat = OUTSIDEGRIDS !see status.h
                 write(OUTDOCK, *) "out of grid"
                 !gistval = -1000.0 
                 ligscore%pose_es(tempsave) = 10.0
                 ligscore%pose_gi(tempsave) = 10.0
                 ligscore%pose_vs(tempsave) = 10.0
                 ligscore%pose_ps(tempsave) = 10.0
                 ligscore%pose_as(tempsave) = 10.0
                 ligscore%pose_is(tempsave) = 0.0
                 ligscore%pose_rs(tempsave) = 0.0
                 !ligscore%pose_ds(tempsave) = 0.0
                 ligscore%pose_rd(tempsave) = 0.0
                 ligscore%pose_hs(tempsave) = 0.0
                 ligscore%pose_score(tempsave) = 1000.0
             
             else
                lengistgi(tempsave) = 0
                call calc_score_conf(conf, options0, db2lig, grids0,
     &            vdwmap0, phimap0, recdes0, ligdes0, gist0, solvmap0,
     &            rec_des_r0, sra, srb, which_grid,
     &            maxtyv,
     &            vdwasum, vdwbsum, eleceel, gistval, solva, solvp, 
     &            solvd, rdscore, hescore, numscored,
     &            gistgridindex,lengistgi(tempsave),
     &            gistgridindex,lengistgi(tempsave))
             !write(OUTDOCK,*) "I AM HERE (D3)"
             !call doflush(OUTDOCK)
             
               if ((vdwasum - vdwbsum) .gt. options0%bmpvdw) then 
                 write(OUTDOCK, *) "high vdw energy."
                 write(OUTDOCK, *) 
     &           "consider cheaking vdw parms, amsol values"
                 !gistval = -1000.0 
                 ligscore%pose_es(tempsave) = 10.0
                 ligscore%pose_gi(tempsave) = 10.0
                 ligscore%pose_vs(tempsave) = 10.0
                 ligscore%pose_ps(tempsave) = 10.0
                 ligscore%pose_as(tempsave) = 10.0
                 ligscore%pose_is(tempsave) = 0.0
                 ligscore%pose_rs(tempsave) = 0.0
                 !ligscore%pose_ds(tempsave) = 0.0
                 ligscore%pose_rd(tempsave) = 0.0
                 ligscore%pose_hs(tempsave) = 0.0
                 ligscore%pose_score(tempsave) = 1000.0
               else 
                ligscore%pose_es(tempsave) = options0%escale*eleceel
                ligscore%pose_gi(tempsave) = options0%gistscale*gistval
                ligscore%pose_vs(tempsave) = options0%vscale*(vdwasum -
     &                                        vdwbsum)
                ligscore%pose_ps(tempsave) = options0%solvscale * solvp
                ligscore%pose_as(tempsave) = options0%solvscale * solva
c               ligscore%pose_reccode      = "rescore"
c               ligscore%pose_es(tempsave) = eleceel
c               ligscore%pose_gi(tempsave) = gistval
c               ligscore%pose_vs(tempsave) = vdwasum - vdwbsum
c               ligscore%pose_ps(tempsave) = solvp
c               ligscore%pose_as(tempsave) = solva
                ligscore%pose_is(tempsave) = 0.0
c                ligscore%pose_rs(tempsave) = 0.0
                ligscore%pose_rs(tempsave) =
     &               options0%rec_energy(tempsave)
                !ligscore%pose_ds(tempsave) = 0.0
                ligscore%pose_rd(tempsave) = options0%rec_des_scale 
     &                                        * rdscore
                ligscore%pose_hs(tempsave) = 0.0
                ligscore%pose_score(tempsave)= 
     &          ligscore%pose_es(tempsave) +
     &          ligscore%pose_gi(tempsave) + 
     &          ligscore%pose_vs(tempsave) + 
     &          ligscore%pose_ps(tempsave) + 
     &          ligscore%pose_as(tempsave) + 
     &          ligscore%pose_is(tempsave) + 
     &          ligscore%pose_rs(tempsave) +
c    &          ligscore%pose_ds(tempsave) +
     &          ligscore%pose_rd(tempsave) +
     &          ligscore%pose_hs(tempsave)
               endif
             endif
             
             
             !write(OUTDOCK,*) "I AM HERE (D4)"
             !call doflush(OUTDOCK)
             if (options0%mol2_minimizer .ne. 1) then
               write(OUTDOCK, outdockline) molct, db2lig%refcod,
     &               ligscore%pose_reccode(tempsave),
     &               match%nmatch, ligscore%numscored,
c    &               match%nmatch, ligscore%numscored,
     &               telaps, db2lig%total_heavy_atoms,
     &               ligscore%setsave(tempsave),
     &               ligscore%orientsave(tempsave),
     &               tempsave, cloudcount, ligscore%pose_es(tempsave),
     &               ligscore%pose_gi(tempsave),
     &               ligscore%pose_vs(tempsave),
     &               ligscore%pose_ps(tempsave),
     &               ligscore%pose_as(tempsave),
     &               ligscore%pose_is(tempsave),
     &               ligscore%pose_rs(tempsave),
C    &               ligscore%pose_ds(tempsave),
     &               ligscore%pose_rd(tempsave),
     &               ligscore%pose_hs(tempsave),
     &               ligscore%pose_score(tempsave)
             endif
             ! if we are using flexible receptor do not write poses to
             ! the mol2 file.  there is a bug in writing to the mol2
             ! file. also we are not minimizing with rescoring so there
             ! is no point.  
             ! TODO: we need to handle flexible receptor better. 
               call mol2write(options0, db2lig, ligscore, match,
c    &            molct, cloudcount, MAXOR, fdmol, outcount, outcount,
     &            molct, cloudcount, MAXOR, fdmol, tempsave, tempsave,
     &            atomscore, grids0, vdwmap0, phimap0,
     &            recdes0, ligdes0, gist0, solvmap0, rec_des_r0,
     &            sra, srb, MAXTYV, allminlimits, allmaxlimits,
     &            time1, total_iter)
             !write(OUTDOCK,*) "I AM HERE (H)"
          enddo ! loop over all receptors.
             
        !write(OUTDOCK,*) "I AM HERE (H)"
        !write(OUTDOCK,*) total_iter
c       flush buffer to OUTDOCK so user can see progress
        if (mod(molct, flush_int) .eq. 0) then
          call doflush(OUTDOCK)
        endif
        !write(OUTDOCK,*) "I AM HERE (I)"
        !write(OUTDOCK, *), "fdmol,", fdmol 
      enddo
      !write(OUTDOCK,*) "I AM HERE"
      if (iend) then
        write(OUTDOCK, *) "end of file encountered"
        call doflush(OUTDOCK)
      endif
      call gzclose(fdmol, istat) !close ligand output file
      return
      end subroutine run_rescore

      end module rescore

c---------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c
c       Copyright (C) 1991 Regents of the University of California
c                         All Rights Reserved.
c
