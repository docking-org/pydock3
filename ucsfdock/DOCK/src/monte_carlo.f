! monte for saved poses
      subroutine run_montecarlo(options0, db2lig, ligscore,
     &           ligscore_mini, match, match_mini, lignumber, cloudnum,
     &           MAXOR, fdmol, count1, outrank, atomscore, grids0,
     &           vdwmap0, phimap0, recdes0, ligdes0, gist0, solvmap0,
     &           rec_des_r0, sra, srb, maxtyv, allminlimits, 
     &           allmaxlimits, time1, total_iterations)

      use db2type
      use ligscoretype
      use matchtype
      use mol2sav
      use optionstype
      use transfm_conf
      use atomscoretype
      use gridstype
      use phimaptype
      use gisttype
      use solvmaptype
      use vdwmaptype
      use rescore_byatom
      use gisttype
      use status

      implicit none
      
      type(options), intent(in) :: options0
      type(db2), intent(inout) :: db2lig
      type(ligscoret), intent(inout) :: ligscore
      type(ligscoret), intent(inout) :: ligscore_mini
      type(matcht), intent(inout) :: match
      type(matcht), intent(inout) :: match_mini
      integer, intent(in) :: lignumber !ligand number (this run), the number of ligands searched on this run
      integer, intent(in) :: cloudnum !cloud number, or just 1 if clouds not used
      integer, intent(in) :: MAXOR ! max orientations
      integer (kind=8), intent(inout) :: fdmol !ligand output file handle (mol2.gz)
      integer, intent(in) :: count1 ! the count1-th saved pose
      integer, intent(in) :: outrank !output ranking
      type(atomscoret), intent(inout) :: atomscore
      type(flexgrids), intent(inout) :: grids0 !precomp values
      type(vdwmap), intent(in) :: vdwmap0
      type(phimap), intent(in) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(in) :: gist0, rec_des_r0
      type(solvmap), intent(in) :: solvmap0
      ! GForran likes these before they're used
      integer, intent(in) :: maxtyv !how many vdw atom types there are
      real, dimension(maxtyv), intent(in) :: sra, srb !square roots of vdw params
      real :: time1

! temporary variables used here
      character (len=160) molbuf !output buffer
      integer istat, bondcount !file status, temp counters
      integer setcount, confcount, conf, atomcount, thisconf
      integer thisconfcount !more temp counters
      integer matchnum !for getting rotation matrices
      integer arbcount !for walking through arbitrary lines of output
      integer tempcount, tempcoord
      integer :: i,j
      real :: time3, time_elapsed

! for minimization
      integer :: min_status
      real :: best_energy ! for minimization
      integer, intent(inout) :: total_iterations !cumulative number of minimization steps performed
      integer :: temp_combo_num
      character (len=255) :: reccode
      real :: setvs, setes,setgi, setas, setrs, setrd,
     &       setps, setds, seths, setscore !scores (whole position)
      real, intent(in), dimension(3) :: allminlimits, allmaxlimits !grid limits in xyz
      character (len=255) :: min_info

      setcount = ligscore%setsave(count1) ! gets this particular conformation set
      matchnum = ligscore%orientsave(count1)
      ! copy ligscore to ligscore_mini
      ligscore_mini%setsave(count1) = ligscore%setsave(count1)
      ligscore_mini%orientsave(count1) = count1
      ligscore_mini%pose_reccode(count1) = ligscore%pose_reccode(count1)
      total_iterations = 0
      temp_combo_num = 1
c      if ( best_energy .le. options0%min_cut) then
      best_energy=ligscore%pose_score(count1)
c      write(6,*), "best_energy before monte carlo = ",best_energy
      !print *, 'best_energy=', best_energy
      !ligscore_mini%pose_score(count1) = ligscore%pose_score(count1)
      min_status = ALLOKAY
      setvs = 0.0
      setes = 0.0
      setgi = 0.0
      setas = 0.0
      setps = 0.0
      setrs = 0.0
      setrd = 0.0
      setds = 0.0
      seths = 0.0
      setscore = 0.0
      min_info = 'YES'

c     copy the original rotation matrix to match_mini
      do i = 1, 3
           do j = 1, 3
               match_mini%rot(i,j,count1) = match%rot(i,j,matchnum)
           enddo
           match_mini%comr(i,count1) = match%comr(i,matchnum)
           match_mini%coml(i,count1) = match%coml(i,matchnum)
      enddo

c    store rotation matrix before minimization
      do i = 1, 3
           do j = 1, 3
                  match_mini%rotab(i,j,count1) =
     &            match_mini%rot(i,j,count1)
           enddo
           match_mini%comrab(i,count1) =
     &     match_mini%comr(i, count1)
      enddo

      !setcount = ligscore%setsave(count1) ! gets this particular conformation set 
      confcount = db2lig%set_conf_len(setcount)

      call dockmin_mc(best_energy, setcount, options0, db2lig,
     &           match_mini, count1,
     &           MAXOR,
     &           ligscore%pose_reccode(count1),
c    &           reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &            setvs,
     &            setes,
     &            setgi,
     &            setas,
     &            setps,
     &            setrs,
     &            setrd,
c     &            setds,
     &            seths,
     &            setscore, min_status,
     &            allminlimits, allmaxlimits, total_iterations) 

      call get_time(time3)
      time_elapsed = time3 - time1
      if (min_status .eq. OUTSIDEGRIDS) then
        min_info = 'YES'
      endif
        ! insertion sort (by energy) a ligand into list of saved poses
      call save_set_score(ligscore_mini, count1, setcount,
     &        setvs,
     &        setes,
     &        setgi,
     &        setas,
     &        setps,
     &        db2lig%set_energy(setcount),
     &        setrs,
     &        setds,
     &        setrd, ! inconsistent location, fix in future
     &        seths,
     &        setscore,
     &        ligscore%pose_reccode(count1),
     &        ligscore%flexnumsave(count1),
     &        ligscore%savedcount,
     &        options0%dockovalent,db2lig%save_cov_coords,
     &        db2lig%total_atoms,db2lig%transfm_coords)


      return
      end subroutine run_montecarlo
