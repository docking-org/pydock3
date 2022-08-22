c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
c--------------------------------------------------------------------
c The rigid body minimimizer only modifies the coordinates of the 
c best molecule (scattered throughout xatm).  This routine points to
c the coordinates in xatm to modify, calls transfm and then scores.
c--------------------------------------------------------------------

      real function rigid_min(confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)

      use db2type
      use optionstype
      use gridstype
      use phimaptype
      use gisttype
      use solvmaptype
      use vdwmaptype
      use phiscore
      use gistscore
      use rec_des_score
      use chemscore
      use solvation_score
      use matchtype
      use ligscoretype
      use atomscoretype
      use score_min
      use transfm_conf
      use status
      use check_conf_og

      implicit none

      type(options), intent(in) :: options0 !useful options
      type(db2), intent(inout) :: db2lig !ligand information
      type(matcht), intent(inout) :: match
      integer, intent(inout):: matchnum !for getting rotation matrices
      integer, intent(in) :: MAXOR ! max orientations
      character (len=255), intent(in) :: reccode !flexible receptor code for this pose
      character (len=255) :: temp_combo !flexible receptor code for this pose
      !integer, intent(in) :: combo_num  ! flexible grid combination
      type(flexgrids), intent(inout) :: grids0 !precomp values
      type(vdwmap), intent(in) :: vdwmap0
      type(phimap), intent(in) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(in) :: gist0, rec_des_r0
      type(solvmap), intent(in) :: solvmap0
      integer, intent(in) :: maxtyv !how many vdw atom types there are
!      integer, intent(in) :: which_grid !which grid to use
      integer, intent(in) :: confset !the confset that the caller would like
                                  !scored
      real, dimension(maxtyv), intent(in) :: sra, srb !square roots of
                                                !vdw parameters per vdwtype
c return variables
      real, dimension(db2lig%total_confs) :: vscoreconf, escoreconf,
     & gscoreconf
      real, dimension(db2lig%total_confs) :: ascoreconf, pscoreconf,
     &          rdscoreconf,hscoreconf !scores by reference
c     &          dscoreconf,hscoreconf !scores by reference
      integer :: numscored !keeps track of how many atom
                                           !positions are evaluated
      integer :: confstat, atomcount
      !character (len=255) :: reccode
      !real, intent(inout) :: best_energy
      integer :: confcount, eachconf, flexgrid
      integer :: temp_combo_num, use_comb_num
      integer, dimension(grids0%group_len) :: temp_combo_list !list of grids
      integer :: count, tempatom, tempcoord, gridnum
      integer, intent(inout) :: min_status
      logical :: temp_status
      real :: temp_set_score, temp_rdscore
      real, dimension(3), intent(in) :: allminlimits, allmaxlimits

      real, intent(inout) ::  vscore, escore, gscore,
     &            ascore, pscore,
     &           rscore, rdscore, hscore, setscore
c     &           rscore, dscore, hscore, setscore


c     flexgrid = 1
c     reccode = getcode(grids0, 1)
c     temp_combo_num =1

c     temp_combo = grids0%all_combos(temp_combo_num, :)
c     reccode = getcode(grids0, temp_combo_num)
      !write (6,*) TRIM(reccode)
      
      do temp_combo_num = 1, grids0%total_combinations
          temp_combo = getcode(grids0, temp_combo_num)
          !grids0%all_combos(temp_combo_num, :)
          !write (6,*) TRIM(temp_combo)
          !write (6,*) TRIM(reccode)
          if (TRIM(temp_combo) .eq. TRIM(reccode)) then
              !write (6,*) "same"
              use_comb_num = temp_combo_num
          endif
          !write(6,*) temp_combo_num, TRIM(temp_combo)
      enddo
      !write(6,*) TRIM(reccode)
      !write(6,*) '\n'
      !call doflush(6)
      !stop

      vscore = 0.0
      escore = 0.0
      gscore = 0.0
      ascore = 0.0
      pscore = 0.0
      rscore = 0.0
      rdscore = 0.0
c      dscore = 0.0
      hscore = 0.0
      setscore = 0.0


c    !    temp_combo = grids0%all_combos(temp_combo_num, :)
c    !    reccode = getcode(grids0, temp_combo_num)
c    !    do tempset = setstart, setend !going through every set as well
c    !      goodwhole = .true.
c    !      wholevs(tempset, temp_combo_num) = 0.
c    !      wholees(tempset, temp_combo_num) = 0.
c    !      wholegi(tempset, temp_combo_num) = 0.
c    !      wholeas(tempset, temp_combo_num) = 0.
c    !      wholeps(tempset, temp_combo_num) = 0.
c    !      wholeds(tempset, temp_combo_num) = 0.
c    !      wholerd(tempset, temp_combo_num) = 0.
c    !      wholehs(tempset, temp_combo_num) = 0.
c    !      wholers(tempset, temp_combo_num) = 0.
c    !      wholescore(tempset, temp_combo_num) = 0.
c    !      do flexgrid = 1, grids0%group_len
c    !        if (goodset(tempset, temp_combo(flexgrid))) then

      
      temp_combo_list = grids0%all_combos(use_comb_num, :)
      ! loop over the grids in the receptor conformation, the one that
      ! is the best. 
      do flexgrid = 1, grids0%group_len
         gridnum = temp_combo_list(flexgrid)
c        write(6,*) gridnum
c        call doflush(6)
c     enddo
c     stop

         

         do confcount = 1, db2lig%set_conf_len(confset) !sets of complete conformations
           eachconf = db2lig%set_conf(confset, confcount)
               !each conf needs to transform the coordinates to the right place
                 !reset these temporary variables
           vscoreconf(eachconf) = 0.0
           escoreconf(eachconf) = 0.0
           gscoreconf(eachconf) = 0.0
           ascoreconf(eachconf) = 0.0
           pscoreconf(eachconf) = 0.0
           rdscoreconf(eachconf) = 0.0
c           dscoreconf(eachconf) = 0.0
           hscoreconf(eachconf) = 0.0
          
          
           call run_transfm_conf(eachconf, matchnum, db2lig, match, 
     &          MAXOR)  !match-picked transformation
          
           temp_status = atom_out_of_grids(eachconf, db2lig,
     &       allminlimits, allmaxlimits) !do outside grid stuff
          
           if (temp_status .eqv. .true.) then !outside grid, quit now
               min_status = OUTSIDEGRIDS
               rigid_min = huge(0.0)
               return
           else


            call recalc_score_min(eachconf,
     &           options0, db2lig, grids0, vdwmap0, phimap0,
     &           recdes0, ligdes0, gist0, solvmap0, rec_des_r0, 
     &           sra, srb, gridnum, maxtyv, vscoreconf(eachconf), 
     &           escoreconf(eachconf),
     &           gscoreconf(eachconf),
     &           ascoreconf(eachconf), pscoreconf(eachconf),
     &           rdscoreconf(eachconf),
c     &           dscoreconf(eachconf),
     &           hscoreconf(eachconf),
     &           numscored, confstat)

            if (confstat .ne. ALLOKAY ) then !bumped, quit now
                write(6,*) "bumped in rigid_min"
                min_status = BUMPED
                rigid_min = huge(0.0)
                return
            else
                vscore = vscore + vscoreconf(eachconf)
                escore = escore + escoreconf(eachconf)
                gscore = gscore + gscoreconf(eachconf)
                pscore = pscore + pscoreconf(eachconf)
                ascore = ascore + ascoreconf(eachconf)
                rdscore = rdscore + rdscoreconf(eachconf)
c                dscore = dscore + dscoreconf(eachconf)
                hscore = hscore + hscoreconf(eachconf)
                temp_rdscore = options0%rec_des_scale * rdscore
                temp_set_score = vscore
     &          + escore + pscore
     &          + ascore + temp_rdscore
     &          + hscore
            endif
        endif
         
          
      enddo ! confcount
      
      rscore = rscore + options0%rec_energy(gridnum)

      !write(6,*) gridnum, options0%rec_energy(gridnum), rscore
      !call doflush(6)
      enddo ! flexgrid
      
      !print *, 'total electrostatics before scaling=', escore
      vscore = options0%vscale * vscore
      escore = options0%escale * escore
      gscore = options0%gistscale * gscore
      pscore = options0%solvscale * pscore
      ascore = options0%solvscale * ascore
      rdscore = options0%rec_des_scale * rdscore
c      dscore = options0%rdscale * dscore
      hscore = options0%rhscale * hscore
      !rscore = options0%rec_energy(flexgrid)
      !rscore = options0%rec_energy(use_comb_num)

      !print *, 'total electrostatics after scaling=', escore


      setscore =
     &    vscore
     &  + escore
     &  + gscore
     &  + pscore
     &  + ascore
     &  + rdscore
c     &  + dscore
     &  + hscore
     &  + rscore + db2lig%set_energy(confset)
      !print *, 'setscore=', setscore
c     call doflush(6)

      rigid_min = setscore

      !stop
      return
      end function rigid_min
