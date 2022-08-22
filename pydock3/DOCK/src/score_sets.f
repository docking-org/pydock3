c just score the conformations directly, no hierarchy, no tree necessary
c this is the most complex subroutine, since it deals with clouds/clusters
c and flexible receptor scoring
      module score_sets
 
      implicit none
      contains

      subroutine calc_score_sets(current_match, dock_status, 
     &    setstart, setend, db2lig, options0, ligscore, ligscoreeach, 
     &    grids0, vdwmap0, phimap0, recdes0, ligdes0, gist0, solvmap0,
     &    rec_des_r0, match,
     &    allminlimits, allmaxlimits,
     &    sra, srb, 
     &    maxor, maxconfs, maxsets, maxatm, maxatmpos,
     &    maxtyv)
 
      use status
      use db2type
      use optionstype
      use ligscoretype
      use gridstype
      use phimaptype
      use gisttype
      use solvmaptype
      use transfm_conf
      use vdwmaptype
      use matchtype
      use doall_conf
      use filenums

c constants, set in max.h and passed in
      integer, intent(in) :: maxor !the size of many arrays that store information about
                    !orientations, is checked several times
      integer, intent(in) :: maxconfs !the size of any array having to do with conformations
      integer, intent(in) :: maxsets !how many max sets can be in a ligand file
      integer, intent(in) :: maxatm !how many atoms max can be in a ligand
      integer, intent(in) :: maxatmpos !how many positions all ligand atoms can take
      integer, intent(in) :: maxtyv !how many vdw atom types there are
c variables
      integer, intent(in) :: current_match !which match number this orientation is from
      integer, intent(inout) :: dock_status !dock_status is a enumerated type, see status.h
      integer, intent(in) :: setstart, setend !controls only doing a certain number of sets
      type(db2), intent(inout) :: db2lig !ligand information
      type(options), intent(in) :: options0 !useful options
      type(ligscoret), intent(inout) :: ligscore
      type(ligscoret), intent(inout) :: ligscoreeach
      type(flexgrids), intent(inout) :: grids0 !precomp values
      type(vdwmap), intent(in) :: vdwmap0
      type(phimap), intent(in) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(in) :: gist0, rec_des_r0
      type(solvmap), intent(in) :: solvmap0
      type(matcht), intent(in) :: match
      real, intent(in), dimension(3) :: allminlimits, allmaxlimits !grid limits in xyz
      real, intent(in), dimension(maxtyv) :: sra, srb !square roots of vdw parameters per vdwtype

      integer :: maingrid, flexgrid !which grid to use
      integer :: tempset !keeps track of which set we are scoring
      integer :: loopconf, tempconf, tempconfnum, oriconf !keeps track of which conf we are scoring, oriconf is the origin conformation, which conf preseids the current conf
      integer :: bestset !keeps track of the best set in terms of energy
      integer :: set_status !status.h status of the current set
      integer :: place_to_save !insertion point into list
      !not accessed outside of this function so made here
      logical, dimension(maxsets, options0%total_receptors) :: goodset !for each set overall, used for flexible receptor
      logical :: goodwhole
      logical :: onegood = .false. ! look to see if none are good. 
      real, dimension(maxsets, grids0%total_combinations) :: 
     &     wholevs, wholees, wholegi, wholeas, wholeps, 
     &     wholeds,wholerd, wholehs, wholers, wholescore !score for all the grids (or just 1 if only 1) 
      real, dimension(maxsets, options0%total_receptors) :: 
     &    setvs, setes, setgi, setas, 
     &    setps, setds, setrd, seths, setscore !scores per set (whole position)
      real, dimension(maxconfs, options0%total_receptors) :: 
     &    confvs, confes,confgi,confas,  !one for each scoring function except
     &    confps, confds, confrd, confhs !internal energy and receptor
      !energy do not have a per-conformation score
      integer, dimension(maxconfs, options0%total_receptors) :: 
     &    conf_status !track status of each conf
      character (len=255) :: reccode
      integer, dimension(grids0%group_len) :: temp_combo !list of grids
      integer :: temp_combo_num !temporary counter for combinations
      integer :: tempcoord,tempatom,count
c     ! this will store the list grid point already asigned to the molecule 
c     integer,dimension(options0%total_receptors,10000) :: gistgridindex
      integer,dimension(5000,db2lig%total_confs,
     &     options0%total_receptors) :: 
     &     gistgridindex
c     integer,dimension(100000) :: gistgridindex_cur_rec
      integer,dimension(db2lig%total_confs,
     &    options0%total_receptors) :: lengistgi ! this is the curent length of list grid point already asigned to the molecule
c     just for debugging, jklyu, 2019.03.28
      character (len=*), parameter :: PDBFORM =
     &  '(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,11X,A2,A2)'
      integer :: atomcount,conf,confcount
      integer atomcoord!used to find the right place to put things

      logical logged_tconf, logged_oconf

      maingrid = grids0%group_part(1, 1) !do the invariant or only grid first

c     intialize to zero the lenth of the gist index list.  
      do tempconf = 1, db2lig%total_confs ! loop over all 
        lengistgi(tempconf,maingrid) = 0
        if (options0%flexible_receptor) then !
          do flexgrid = 2, options0%total_receptors !
             lengistgi(tempconf,flexgrid) = 0
          enddo
        endif
      enddo
      !first do the first conf, which is rigid component
      !tempconf = 1 ! db2lig%conf_sorted(1) will equal 1
      tempconf = db2lig%conf_sorted(1) ! may not equal 1
c     write(*,*) "starting tempconf = ", tempconf
      !gistgridindex_cur_rec = gistgridindex(:,maingrid)
c      write(6,*) "I AM in score_sets"
c      write(6,*) "conf_status before is",conf_status(tempconf, maingrid)
      call doflush(6)
      conf_status(tempconf, maingrid) = 
     &    run_conf(tempconf, current_match,
     &    db2lig, options0, grids0, maingrid,
     &    vdwmap0, phimap0, recdes0, ligdes0, gist0, solvmap0,
     &    rec_des_r0, match,
     &    confvs(tempconf, maingrid), 
     &    confes(tempconf, maingrid), 
     &    confgi(tempconf, maingrid), 
     &    confas(tempconf, maingrid), 
     &    confps(tempconf, maingrid), 
     &    confds(tempconf, maingrid),
     &    confrd(tempconf, maingrid),
     &    confhs(tempconf, maingrid),
     &    allminlimits, allmaxlimits,
     &    sra, srb, 
     &    maxor, maxconfs, maxatm, maxatmpos, maxtyv, 
c    &    ligscore%numscored, gistgridindex_cur_rec,
c    &    ligscore%numscored, gistgridindex(:,tempconf,maingrid),
c    &    lengistgi(tempconf,maingrid))
     &    ligscore%numscored,
     &    gistgridindex(:,tempconf,
     &        maingrid),
     &    lengistgi(tempconf,maingrid),
     &    gistgridindex(:,
     &        db2lig%conf_conection(tempconf),
     &        maingrid),
     &    lengistgi(db2lig%conf_conection(tempconf),
     &        maingrid))


c      write(6,*) "conf_status after is", conf_status(tempconf, maingrid)
c      just for debugging turn off the rigidbump check
      if ((conf_status(tempconf, maingrid) .ne. ALLOKAY) .or. 
     &    (confvs(tempconf, maingrid) .gt. 
     &    options0%bmpvdwrigid)) then
        dock_status = BUMPED !the new version of bumping, rely on the vdw score
      else !the rigid component isn't terrible, go ahead and score
        if ((setstart .eq. 1) .and. 
     &      (setend .eq. db2lig%total_sets)) then
          !not doing clouds, just do all confs at once, we'll need them all
          do loopconf = 2, db2lig%total_confs !do all other confs
            !for tempconf, rotate atoms in that conf in place. score them.
            !save scores to confes, ..vs, ..ps, ..as
            !gistgridindex_cur_rec = gistgridindex(:,maingrid)
            tempconf = db2lig%conf_sorted(loopconf)
            oriconf = db2lig%conf_conection(tempconf) 
c           write(*,*) loopconf, "->", tempconf,
c    &                 "->", oriconf
            call doflush(OUTDOCK)
            if (tempconf.eq.0) then
               if (.not. logged_tconf) then ! really annoying when OUTDOCK is cluttered with things like this
                    write(*,*) "total_confs =", db2lig%total_confs
                    write(*,*) "Warning. tempconf = 0"
                    write(*,*) loopconf, "->", tempconf,
     &                 "->", oriconf
                    call doflush(OUTDOCK)
                    logged_tconf = .true.
               endif
               
               tempconf = loopconf
            else if (oriconf.eq.0) then
               if (.not. logged_oconf) then
                    write(*,*) "Warning. oriconf=0"
                    write(*,*) loopconf, "->", tempconf,
     &                 "->", oriconf
                    call doflush(OUTDOCK)
                    logged_oconf = .true.
               endif
               oriconf = tempconf
            endif
            
            conf_status(tempconf, maingrid) = run_conf(
     &          tempconf, current_match,
     &          db2lig, options0, grids0, maingrid,
     &          vdwmap0, phimap0, recdes0, ligdes0, gist0, solvmap0, 
     &          rec_des_r0, match,
     &          confvs(tempconf, maingrid), 
     &          confes(tempconf, maingrid), 
     &          confgi(tempconf, maingrid), 
     &          confas(tempconf, maingrid), 
     &          confps(tempconf, maingrid), 
     &          confds(tempconf, maingrid),
     &          confrd(tempconf, maingrid),
     &          confhs(tempconf, maingrid),
     &          allminlimits, allmaxlimits,
     &          sra, srb, 
     &          maxor, maxconfs, maxatm, maxatmpos, maxtyv, 
c    &          ligscore%numscored,gistgridindex_cur_rec,
c    &          ligscore%numscored,gistgridindex(:,tempconf,maingrid),
c    &          lengistgi(tempconf,maingrid))
     &          ligscore%numscored,
     &          gistgridindex(:,tempconf,
     &              maingrid),
     &          lengistgi(tempconf,maingrid),
     &          gistgridindex(:,
c    &              db2lig%conf_conection(tempconf),
     &              oriconf,
     &              maingrid),
c    &          lengistgi(db2lig%conf_conection(tempconf),
     &          lengistgi(oriconf,
     &              maingrid))

          enddo
          if (options0%flexible_receptor) then !may want to score more
            do loopconf = 1, db2lig%total_confs !do all other confs
              tempconf = db2lig%conf_sorted(loopconf)
c              write(*,*) "tempconf (sorted) = ",tempconf
c     &                  ,"loopconf = ", loopconf      
              if (conf_status(tempconf, maingrid) .eq. ALLOKAY) then
                do flexgrid = 2, options0%total_receptors !all other recs
                  !gistgridindex_cur_rec = gistgridindex(:,flexgrid)
                  conf_status(tempconf, flexgrid) = run_conf(
     &                tempconf, current_match,
     &                db2lig, options0, grids0, flexgrid,
     &                vdwmap0, phimap0, recdes0, ligdes0, gist0,
     &                solvmap0, rec_des_r0, match,
     &                confvs(tempconf, flexgrid), 
     &                confes(tempconf, flexgrid), 
     &                confgi(tempconf, flexgrid), 
     &                confas(tempconf, flexgrid), 
     &                confps(tempconf, flexgrid), 
     &                confds(tempconf, flexgrid),
     &                confrd(tempconf, flexgrid),
     &                confhs(tempconf, flexgrid),
     &                allminlimits, allmaxlimits,
     &                sra, srb, 
     &                maxor, maxconfs, maxatm, maxatmpos, maxtyv, 
c    &                ligscore%numscored,gistgridindex_cur_rec,
     &                ligscore%numscored,gistgridindex(:,tempconf,
     &                    flexgrid),
     &                lengistgi(tempconf,flexgrid),
     &                gistgridindex(:,
     &                    db2lig%conf_conection(tempconf),
     &                    flexgrid),
     &                lengistgi(db2lig%conf_conection(tempconf),
     &                    flexgrid))
                enddo
              endif
            enddo
          endif
          !all conformations have been scored if possible
        else !we only need some of them.. don't do unnecessary work
          !only do for necessary confs, so mark all as undone for now
          do loopconf = 2, db2lig%total_confs !do all other confs
            tempconf = db2lig%conf_sorted(loopconf)
            do flexgrid = 1, options0%total_receptors !for all of them
              conf_status(tempconf, flexgrid) = NOTCHECKED !set to this, check later
            enddo
          enddo
        endif
        do flexgrid = 1, options0%total_receptors !go through every grid & 
                            !set, assemble scores for each combination
          set_status = NOTCHECKED !temporary for this set
          do tempset = setstart, setend !going through every set as well
            !add up scores from each conf to find the best one.
            setvs(tempset, flexgrid) = 0. !set to 0
            setes(tempset, flexgrid) = 0.
            setgi(tempset, flexgrid) = 0.
            setps(tempset, flexgrid) = 0.
            setas(tempset, flexgrid) = 0.
            setds(tempset, flexgrid) = 0.
            setrd(tempset, flexgrid) = 0.
            seths(tempset, flexgrid) = 0.
            !broken/clashed check on next line, if necessary
            if ((db2lig%set_broken(tempset) .eq. 0) .or. 
     &            (options0%check_clash .eqv. .false.)) then 
              if ((db2lig%set_above_strain(tempset) .eq. 0) .or.
     &            (options0%check_strain .eqv. .false.)) then 
                set_status = ALLOKAY
                do tempconfnum = 1, db2lig%set_conf_len(tempset)
                  tempconf = db2lig%set_conf(tempset, tempconfnum)
                  !build on the fly if necessary
                  if (conf_status(tempconf, flexgrid) .eq. 
     &                NOTCHECKED) then !if it hasn't been checked already
                    !gistgridindex_cur_rec = gistgridindex(:,flexgrid)
                    conf_status(tempconf, flexgrid) = run_conf(
     &                  tempconf, current_match,
     &                  db2lig, options0, grids0, flexgrid,
     &                  vdwmap0, phimap0, recdes0, ligdes0, gist0,
     &                  solvmap0, rec_des_r0, match,
     &                  confvs(tempconf, flexgrid), 
     &                  confes(tempconf, flexgrid), 
     &                  confgi(tempconf, flexgrid), 
     &                  confas(tempconf, flexgrid), 
     &                  confps(tempconf, flexgrid), 
     &                  confds(tempconf, flexgrid),
     &                  confrd(tempconf, flexgrid),
     &                  confhs(tempconf, flexgrid),
     &                  allminlimits, allmaxlimits,
     &                  sra, srb, 
     &                  maxor, maxconfs, maxatm, maxatmpos, maxtyv, 
c    &                  ligscore%numscored,gistgridindex_cur_rec,
     &                  ligscore%numscored,gistgridindex(:,tempconf,
     &                      flexgrid),
     &                  lengistgi(tempconf,flexgrid),
     &                  gistgridindex(:,
     &                      db2lig%conf_conection(tempconf),
     &                      flexgrid),
     &                  lengistgi(db2lig%conf_conection(tempconf),
     &                      flexgrid))
                  endif             

                  if (conf_status(tempconf, flexgrid) .ne. ALLOKAY) then
                    set_status = conf_status(tempconf, flexgrid)
                    exit !break out of loop for this set
                  endif
                  setvs(tempset, flexgrid) = setvs(tempset, flexgrid) + 
     &                confvs(tempconf, flexgrid)
                  setes(tempset, flexgrid) = setes(tempset, flexgrid) + 
     &                confes(tempconf, flexgrid)
                  setgi(tempset, flexgrid) = setgi(tempset, flexgrid) + 
     &                confgi(tempconf, flexgrid)
                  setps(tempset, flexgrid) = setps(tempset, flexgrid) + 
     &                confps(tempconf, flexgrid)
                  setas(tempset, flexgrid) = setas(tempset, flexgrid) +
     &                confas(tempconf, flexgrid)
                  setds(tempset, flexgrid) = setds(tempset, flexgrid) + 
     &                confds(tempconf, flexgrid)
                  setrd(tempset, flexgrid) = setrd(tempset, flexgrid) + 
     &                confrd(tempconf, flexgrid)
                  seths(tempset, flexgrid) = seths(tempset, flexgrid) + 
     &                confhs(tempconf, flexgrid)
                enddo
              else !above strain
                set_status = STRAIN
                dock_status = STRAIN
              endif      
            else !broken conf
              set_status = CLASHES !nothing else needs checked, no need to sum
              dock_status = CLASHES ! this is so that if everything clashes then it reports as a clash.
            endif
            goodset(tempset, flexgrid) = .false.
            if (set_status .eq. ALLOKAY) then !each conf scored
              goodset(tempset, flexgrid) = .true.
              !moved scaling to here for grids stuff, internal energy is 
              !done when it is read in
              setvs(tempset, flexgrid) = options0%vscale * 
     &            setvs(tempset, flexgrid)
              setes(tempset, flexgrid) = options0%escale * 
     &            setes(tempset, flexgrid)
              setgi(tempset, flexgrid) = options0%gistscale * 
     &            setgi(tempset, flexgrid)
              setas(tempset, flexgrid) = options0%solvscale * 
     &            setas(tempset, flexgrid)
              setps(tempset, flexgrid) = options0%solvscale * 
     &            setps(tempset, flexgrid)
              setds(tempset, flexgrid) = options0%rdscale * 
     &            setds(tempset, flexgrid)
              setrd(tempset, flexgrid) = options0%rec_des_scale * 
     &            setrd(tempset, flexgrid)
              seths(tempset, flexgrid) = options0%rhscale * 
     &            seths(tempset, flexgrid)
              setscore(tempset, flexgrid) = setvs(tempset, flexgrid) 
     &            + setes(tempset, flexgrid)
     &            + setgi(tempset, flexgrid)
     &            + setps(tempset, flexgrid) + setas(tempset, flexgrid) 
     &            + db2lig%set_energy(tempset) 
     &            + setds(tempset, flexgrid)
     &            + setrd(tempset, flexgrid)
     &            + seths(tempset, flexgrid)
     &            + options0%rec_energy(flexgrid)
            endif
          enddo      
        enddo
        do temp_combo_num = 1, grids0%total_combinations
          temp_combo = grids0%all_combos(temp_combo_num, :)
          reccode = getcode(grids0, temp_combo_num)        
          do tempset = setstart, setend !going through every set as well
            goodwhole = .true.
            wholevs(tempset, temp_combo_num) = 0.
            wholees(tempset, temp_combo_num) = 0.
            wholegi(tempset, temp_combo_num) = 0.
            wholeas(tempset, temp_combo_num) = 0.
            wholeps(tempset, temp_combo_num) = 0.
            wholeds(tempset, temp_combo_num) = 0.
            wholerd(tempset, temp_combo_num) = 0.
            wholehs(tempset, temp_combo_num) = 0.
            wholers(tempset, temp_combo_num) = 0.
            wholescore(tempset, temp_combo_num) = 0.
            do flexgrid = 1, grids0%group_len
              if (goodset(tempset, temp_combo(flexgrid))) then
                wholevs(tempset, temp_combo_num) = 
     &              wholevs(tempset, temp_combo_num) + 
     &              setvs(tempset, temp_combo(flexgrid))
                wholees(tempset, temp_combo_num) = 
     &              wholees(tempset, temp_combo_num) + 
     &              setes(tempset, temp_combo(flexgrid))
               wholegi(tempset, temp_combo_num) =
     &              wholegi(tempset, temp_combo_num) +
     &              setgi(tempset, temp_combo(flexgrid))
                wholeas(tempset, temp_combo_num) = 
     &              wholeas(tempset, temp_combo_num) + 
     &              setas(tempset, temp_combo(flexgrid))
                wholeps(tempset, temp_combo_num) = 
     &              wholeps(tempset, temp_combo_num) + 
     &              setps(tempset, temp_combo(flexgrid))
                wholeds(tempset, temp_combo_num) = 
     &              wholeds(tempset, temp_combo_num) + 
     &              setds(tempset, temp_combo(flexgrid))
                wholerd(tempset, temp_combo_num) = 
     &              wholerd(tempset, temp_combo_num) + 
     &              setrd(tempset, temp_combo(flexgrid))
                wholehs(tempset, temp_combo_num) = 
     &              wholehs(tempset, temp_combo_num) + 
     &              seths(tempset, temp_combo(flexgrid))
                wholers(tempset, temp_combo_num) = 
     &              wholers(tempset, temp_combo_num) + 
     &              options0%rec_energy(temp_combo(flexgrid)) !note difference!!
                wholescore(tempset, temp_combo_num) = 
     &              wholescore(tempset, temp_combo_num) + 
     &              setscore(tempset, temp_combo(flexgrid))
              else
                goodwhole = .false.
                exit !break out, no need to do more
              endif
            enddo
            if (goodwhole) then !don't save if not good
              onegood = .true. ! there is at least one that is good. 

              !all logic about whether or not to save (if score good enough,
              ! etc.) is in the functions called here

              !if dockovalent supply the ligands coords to save
               if (options0%dockovalent .eqv. .true.) then
                  do tempconfnum = 1, db2lig%set_conf_len(tempset)
                     tempconf = db2lig%set_conf(tempset, tempconfnum)
                     do tempcoord = db2lig%conf_coord(1, tempconf), 
     &                    db2lig%conf_coord(2, tempconf)
                        tempatom = db2lig%coord_index(2, tempcoord)
                        do count = 1, 3
                           db2lig%transfm_coords(count, tempatom) = 
     &                          db2lig%coords(count,tempcoord)
                        enddo
                     enddo
                  enddo
               endif
c              write(6,*) current_match, " ", tempset, " ",
c     &            wholescore(tempset, temp_combo_num)," ",
c     &            wholevs(tempset, temp_combo_num)
c              debug, jklyu, 2019.03.28
c              open(unit=OUTPDB, 
c     &        file='test.pdb', action='write')
c             for debugging, transfm the whole set based on the match
c              do confcount = 1, db2lig%set_conf_len(tempset) !sets of complete conformations
c                 conf = db2lig%set_conf(tempset, confcount)
c                 !each conf needs to transform the coordinates to the right place
c                 call run_transfm_conf(conf, current_match, db2lig,
c     &                match, maxor)  !match-picked transformation
c              enddo
c              write(6,"(A18,I6,A10,I4,F6.2,F6.2)") 
c     &            "REMARK match_num: ",
c     &            current_match,
c     &            " set_num: ", tempset,
c     &            wholescore(tempset, temp_combo_num),
c     &            wholevs(tempset, temp_combo_num)
c             the codes below are trying to write out the whole molecule
c              do atomcount = 1, db2lig%total_atoms
c                write(6, PDBFORM) "ATOM  ",
c     &           db2lig%atom_num(atomcount),
c     &           db2lig%atom_name(atomcount), "","LIG", "A", 1,
c     &           "", db2lig%transfm_coords(1, atomcount),
c     &           db2lig%transfm_coords(2, atomcount),
c     &           db2lig%transfm_coords(3, atomcount),
c     &           0.00, 0.00,
c     &           db2lig%atom_name(atomcount), ""
c              enddo
c              write(6,*) "ENDMDL"
c             the codes below are trying to write out the rigid
c             fragment for each molecule
c             frist conf is the rigid fragment
c              do atomcoord = db2lig%conf_coord(1, 1), 
c     &             db2lig%conf_coord(2, 1)
c                atomcount = db2lig%coord_index(2, atomcoord)
c                write(6, PDBFORM) "ATOM  ",
c     &           db2lig%atom_num(atomcount),
c     &           db2lig%atom_name(atomcount), "","LIG", "A", 1,
c     &           "", db2lig%transfm_coords(1, atomcount),
c     &           db2lig%transfm_coords(2, atomcount),
c     &           db2lig%transfm_coords(3, atomcount),
c     &           0.00, 0.00,
c     &           db2lig%atom_name(atomcount), ""
c              enddo
c              write(6,"(A6)") "ENDMDL"

c             write(6,*) tempset, db2lig%set_max_strain(tempset),
c    &                   db2lig%set_total_strain(tempset) 
              call save_set_score(ligscore, current_match, tempset, 
     &            wholevs(tempset, temp_combo_num), 
     &            wholees(tempset, temp_combo_num), 
     &            wholegi(tempset, temp_combo_num), 
     &            wholeas(tempset, temp_combo_num), 
     &            wholeps(tempset, temp_combo_num), 
c     &            db2lig%set_energy(tempset), 
     &            db2lig%set_total_strain(tempset), 
c     &            wholers(tempset, temp_combo_num), 
     &            db2lig%set_max_strain(tempset), 
     &            wholeds(tempset, temp_combo_num), 
     &            wholerd(tempset, temp_combo_num), 
     &            wholehs(tempset, temp_combo_num), 
     &            wholescore(tempset, temp_combo_num), 
     &            reccode, temp_combo_num, options0%nsav,
     &            options0%dockovalent,db2lig%save_cov_coords,
     &            db2lig%total_atoms,db2lig%transfm_coords)
              if (options0%flexible_receptor .and.
     &            options0%score_each_flex) then
                call save_set_score_each_flex(ligscoreeach, 
     &              current_match, tempset, 
     &              wholevs(tempset, temp_combo_num), 
     &              wholees(tempset, temp_combo_num), 
     &              wholegi(tempset, temp_combo_num), 
     &              wholeas(tempset, temp_combo_num), 
     &              wholeps(tempset, temp_combo_num), 
c     &              db2lig%set_energy(tempset), 
     &              db2lig%set_total_strain(tempset), 
c     &              wholers(tempset, temp_combo_num), 
     &              db2lig%set_max_strain(tempset), 
     &              wholeds(tempset, temp_combo_num), 
     &              wholerd(tempset, temp_combo_num), 
     &              wholehs(tempset, temp_combo_num), 
     &              wholescore(tempset, temp_combo_num), 
     &              reccode, temp_combo_num)
              endif
            endif
          enddo
        enddo
        if (onegood) then
           dock_status = ALLOKAY !at least one conformation so this run is good
        endif
      endif
      return
      end subroutine calc_score_sets

      end module score_sets

