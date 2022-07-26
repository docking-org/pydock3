c this file holds the things necessary to write out a bunch of specific
c ligand poses. this is in combination with everything in molecule and hierarchy
c that is read in/set for each ligand. here we only need to save
c the coml, comr, rot matrix and set for each pose we want to output. 
c also save the overall energies of each type (v,e,p,a) (plus internal energy i
c and receptor energy r and receptor desolvation d)
c
c also contains 2 routines to save scores into ligscore type classes
      module ligscoretype

      implicit none

      type ligscoret
        logical, dimension(:), allocatable :: good_pose !whether or not this is good (only used in ligscoreeach incarnation)
        integer, dimension(:), allocatable :: setsave !save the set 
        integer, dimension(:), allocatable :: orientsave !save the orientation number
        integer, dimension(:), allocatable :: flexnumsave !save the flexible receptor combination number
        real, dimension(:), allocatable :: pose_score !overall energy
        real, dimension(:), allocatable :: pose_vs !overall van der waals score
        real, dimension(:), allocatable :: pose_es !overall electrostatic score
        real, dimension(:), allocatable :: pose_gi !overall gist energy
        real, dimension(:), allocatable :: pose_ps !overall polar desolvation score
        real, dimension(:), allocatable :: pose_as !overall apolar desolvation score
        real, dimension(:), allocatable :: pose_is !internal energy of each saved pose
        real, dimension(:), allocatable :: pose_rs !receptor energy of each saved pose
        real, dimension(:), allocatable :: pose_ds !receptor desolvation energy of each saved pose
        real, dimension(:), allocatable :: pose_rd !overall rec_d score
        real, dimension(:), allocatable :: pose_hs !hydrophobic effect of receptor 
        character (len=255), dimension(:), allocatable :: pose_reccode !flexible receptor code for this pose
        integer :: savedcount !number saved so far
        integer :: numscored !how many ligand atoms were scored
      end type ligscoret

      contains

!sets variable so that we know the score is bad for each one. don't write bad ones.
      subroutine reset_ligscore(ligsc)

      type(ligscoret), intent(inout) :: ligsc

      integer :: count

      do count = 1, size(ligsc%good_pose, 1)
        ligsc%good_pose(count) = .false.
      enddo
      return
      end subroutine reset_ligscore

!subroutine allocates the things declared here
      subroutine allocate_ligscore(ligsc, posestosave, OUTDOCK)

      type(ligscoret), intent(inout) :: ligsc
      integer, intent(in) :: posestosave
      integer, intent(in) :: OUTDOCK

      integer alloc_stat

      if (allocated(ligsc%good_pose)) then
        deallocate(ligsc%good_pose)
      endif
      if (allocated(ligsc%setsave)) then
        deallocate(ligsc%setsave)
      endif
      if (allocated(ligsc%orientsave)) then
        deallocate(ligsc%orientsave)
      endif
      if (allocated(ligsc%flexnumsave)) then
        deallocate(ligsc%flexnumsave)
      endif
      if (allocated(ligsc%pose_score)) then
        deallocate(ligsc%pose_score)
      endif
      if (allocated(ligsc%pose_vs)) then
        deallocate(ligsc%pose_vs)
      endif
      if (allocated(ligsc%pose_es)) then
        deallocate(ligsc%pose_es)
      endif
      if (allocated(ligsc%pose_gi)) then
        deallocate(ligsc%pose_gi)
      endif
      if (allocated(ligsc%pose_as)) then
        deallocate(ligsc%pose_as)
      endif
      if (allocated(ligsc%pose_ps)) then
        deallocate(ligsc%pose_ps)
      endif
      if (allocated(ligsc%pose_is)) then
        deallocate(ligsc%pose_is)
      endif
      if (allocated(ligsc%pose_rs)) then
        deallocate(ligsc%pose_rs)
      endif
      if (allocated(ligsc%pose_ds)) then
        deallocate(ligsc%pose_ds)
      endif
      if (allocated(ligsc%pose_rd)) then
        deallocate(ligsc%pose_rd)
      endif
      if (allocated(ligsc%pose_hs)) then
        deallocate(ligsc%pose_hs)
      endif
      if (allocated(ligsc%pose_reccode)) then
        deallocate(ligsc%pose_reccode)
      endif
      allocate(ligsc%good_pose(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate good_pose"
        stop
      endif
      allocate(ligsc%setsave(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate setsave"
        stop
      endif
      allocate(ligsc%orientsave(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate orientsave"
        stop
      endif
      allocate(ligsc%flexnumsave(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate flexnumsave"
        stop
      endif
      allocate(ligsc%pose_score(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate pose_score"
        stop
      endif
      allocate(ligsc%pose_vs(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate pose_vs"
        stop
      endif
      allocate(ligsc%pose_es(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate pose_es"
        stop
      endif
      allocate(ligsc%pose_gi(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate pose_gi"
        stop
      endif
      allocate(ligsc%pose_as(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate pose_as"
        stop
      endif
      allocate(ligsc%pose_ps(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate pose_ps"
        stop
      endif
      allocate(ligsc%pose_is(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate pose_is"
        stop
      endif
      allocate(ligsc%pose_rs(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate pose_rs"
        stop
      endif
      allocate(ligsc%pose_ds(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate pose_ds"
        stop
      endif
      allocate(ligsc%pose_rd(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate pose_rd"
        stop
      endif
      allocate(ligsc%pose_hs(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate pose_hs"
        stop
      endif
      allocate(ligsc%pose_reccode(posestosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate pose_reccode"
        stop
      endif
      
      return
      end subroutine allocate_ligscore

c insertion sort (by energy) a ligand into list of saved poses
      subroutine save_set_score(ligsc, currentmatch, setnum, 
     &    thisvs, thises,thisgi, thisas, thisps, 
     &    thisis, thisrs, thisds, thisrd, thishs, 
     &    thisscore, thisreccode, thisrecnum,
     &    nsav,dockovalent,save_cov_coords,total_atoms,
     &    transfm_coords)

      type(ligscoret), intent(inout) :: ligsc
      integer, intent(in) :: currentmatch !current orientation used to generate this score 
      integer, intent(in) :: setnum !current set number used
      real, intent(in) :: thisvs, thises,thisgi, thisas, thisps,
     &     thisis, thisrs, thisds, thisrd, thishs, thisscore !scores per set (whole position)
c how many poses to save max
      character (len=255), intent(in) :: thisreccode !code like 1.2.4
      integer, intent(in) :: thisrecnum !which combination this is (1-totalcombinations)
      integer, intent(in) :: nsav !max number to save, length of many arrays
      logical, intent(in) :: dockovalent !is this a covalent run
      real, dimension(:,:,:), intent(inout) :: save_cov_coords ! array of saved pose coords
      integer, intent(in) :: total_atoms !numer of atoms in this ligand for pose copying
      real, dimension(:,:), intent(in) :: transfm_coords !the current set pose coords

      integer place_to_save !where to insert
      integer tempcount !used to insertion sort
      character (len=255) :: reccodeinsert
      integer count,atcount !use as for loops indices

      !savedcount is the number saved already. nsav is the array dimesions/size
      ! requested by the user
      ! thisscore, thises etc have data
      !figure out where this one goes
      reccodeinsert = adjustl(thisreccode) !adjustl is builtin function, removes whitespace at left
      place_to_save = 0
      
      if (ligsc%savedcount .eq. 0) then 
        place_to_save = 1 !always save if nothing here yet
      else
        if (thisscore .lt. ligsc%pose_score(ligsc%savedcount)) then
          !find the position to insert into
          place_to_save = ligsc%savedcount
          
          !count backwards to find where the energy is better
c         do while ((place_to_save .gt. 1) .and. 
c    &        (thisscore .lt. 
c    &        ligsc%pose_score(place_to_save-1)))
          do while (place_to_save .gt. 1) 
            if (thisscore .ge. 
     &        ligsc%pose_score(place_to_save-1)) then
                 exit ! breaks out of loop
            endif
            place_to_save = place_to_save - 1
          enddo
        !put this next to the last saved one 
        else if ((thisscore .gt.
     &      ligsc%pose_score(ligsc%savedcount)) .and.
     &      (ligsc%savedcount .lt. nsav)) then
          place_to_save = ligsc%savedcount + 1
        endif
      endif
      if (place_to_save .gt. 0) then !0 means don't save
        !move everything back one step
        tempcount = ligsc%savedcount + 1
        if (tempcount .gt. nsav) then
          tempcount = nsav !don't  write over end of array
        endif
        ligsc%savedcount = tempcount !adjust total number saved
        do while (tempcount .gt. place_to_save)
          ligsc%setsave(tempcount) = ligsc%setsave(tempcount - 1)
          ligsc%orientsave(tempcount) = 
     &        ligsc%orientsave(tempcount - 1)
          ligsc%flexnumsave(tempcount) = 
     &        ligsc%flexnumsave(tempcount - 1)
          ligsc%pose_score(tempcount) = 
     &        ligsc%pose_score(tempcount - 1)
          ligsc%pose_vs(tempcount) = ligsc%pose_vs(tempcount - 1)
          ligsc%pose_es(tempcount) = ligsc%pose_es(tempcount - 1)
          ligsc%pose_gi(tempcount) = ligsc%pose_gi(tempcount - 1)
          ligsc%pose_ps(tempcount) = ligsc%pose_ps(tempcount - 1)
          ligsc%pose_as(tempcount) = ligsc%pose_as(tempcount - 1)
          ligsc%pose_is(tempcount) = ligsc%pose_is(tempcount - 1)
          ligsc%pose_rs(tempcount) = ligsc%pose_rs(tempcount - 1)
          ligsc%pose_ds(tempcount) = ligsc%pose_ds(tempcount - 1)
          ligsc%pose_rd(tempcount) = ligsc%pose_rd(tempcount - 1)
          ligsc%pose_hs(tempcount) = ligsc%pose_hs(tempcount - 1)
          ligsc%pose_reccode(tempcount) = 
     &        ligsc%pose_reccode(tempcount - 1)
c          !if dockovalent also keep track of coords array saved for each pose
          if (dockovalent .eqv. .true.) then
             do count=1,3
                do atcount=1,total_atoms
                   save_cov_coords(tempcount,count,atcount) = 
     &             save_cov_coords(tempcount-1,count,atcount)
                enddo
             enddo
          endif
          tempcount = tempcount - 1 !step backwards now, do it again (maybe)
        enddo
        !save the newest ligand data here
        ligsc%setsave(place_to_save) = setnum
        ligsc%orientsave(place_to_save) = currentmatch
        ligsc%flexnumsave(place_to_save) = thisrecnum
        ligsc%pose_score(place_to_save) = thisscore
        ligsc%pose_vs(place_to_save) = thisvs
        ligsc%pose_es(place_to_save) = thises
        ligsc%pose_gi(place_to_save) = thisgi
        ligsc%pose_ps(place_to_save) = thisps
        ligsc%pose_as(place_to_save) = thisas
        ligsc%pose_is(place_to_save) = thisis
        ligsc%pose_rs(place_to_save) = thisrs
        ligsc%pose_ds(place_to_save) = thisds
        ligsc%pose_rd(place_to_save) = thisrd
        ligsc%pose_hs(place_to_save) = thishs
        ligsc%pose_reccode(place_to_save) = reccodeinsert
        !if dockovalent also keep track of coords array saved for each pose
        if (dockovalent .eqv. .true.) then
           do atcount=1,total_atoms
              do count=1,3
                 save_cov_coords(place_to_save,count,atcount) = 
     &                transfm_coords(count,atcount)
              enddo
           enddo
        endif     
      endif
      return
      end subroutine save_set_score

! save the best score for each combination of flexible receptor.
      subroutine save_set_score_each_flex(ligsceach, currentmatch,
     &    setnum, thisvs, thises, thisgi, thisas, thisps, 
     &    thisis, thisrs, thisds,thisrd, thishs,
     &    thisscore, thisreccode, thisrecnum)

      type(ligscoret), intent(inout) :: ligsceach !where all scores are scored
      integer, intent(in) :: currentmatch !current orientation used to generate this score 
      integer, intent(in) :: setnum !current set number used
      real, intent(in) :: thisvs, thises, thisgi, thisas, thisps,
     &      thisis, thisrs, thisds,thisrd, thishs, thisscore !scores per set (whole position)
c how many poses to save max
      character (len=255), intent(in) :: thisreccode !code like 1.2.4
      integer, intent(in) :: thisrecnum !which combination this is (1-totalcombinations)

      character (len=255) :: reccodeinsert !make sure code is correct size

      ! thisscore, thises etc have data
      reccodeinsert = adjustl(thisreccode) !adjustl is builtin function, removes whitespace at left
      !check if this score is better than previous score
      if ((.not. ligsceach%good_pose(thisrecnum)) .or.
     &    (thisscore .lt. ligsceach%pose_score(thisrecnum))) then
        !save the best ligand data for this receptor combination here
        ligsceach%good_pose(thisrecnum) = .true. 
        ligsceach%setsave(thisrecnum) = setnum
        ligsceach%orientsave(thisrecnum) = currentmatch
        ligsceach%flexnumsave(thisrecnum) = thisrecnum
        ligsceach%pose_score(thisrecnum) = thisscore
        ligsceach%pose_vs(thisrecnum) = thisvs
        ligsceach%pose_es(thisrecnum) = thises
        ligsceach%pose_gi(thisrecnum) = thisgi
        ligsceach%pose_ps(thisrecnum) = thisps
        ligsceach%pose_as(thisrecnum) = thisas
        ligsceach%pose_is(thisrecnum) = thisis
        ligsceach%pose_rs(thisrecnum) = thisrs
        ligsceach%pose_ds(thisrecnum) = thisds
        ligsceach%pose_rd(thisrecnum) = thisrd
        ligsceach%pose_hs(thisrecnum) = thishs
        ligsceach%pose_reccode(thisrecnum) = reccodeinsert
      endif
      return
      end subroutine save_set_score_each_flex

      end module ligscoretype

