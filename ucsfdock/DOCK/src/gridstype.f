! all allocatable large grids stored here/
      module gridstype
 
      implicit none
      type grids
        real, allocatable, dimension(:, :, :) :: gistgrid
        logical, allocatable, dimension(:) :: gisttagvoxelarray  ! this is used to tag gist voxels, in the scoring function
        logical, allocatable, dimension(:) :: gistneighborvoxelarray  ! this is used to id gist voxels as neighbor, in the scoring function
        real, allocatable, dimension(:, :, :) :: rec_des_grid
        real, allocatable, dimension(:, :, :) :: hrec_des_grid
        real, allocatable, dimension(:, :, :, :) :: rec_des_precomp
        real, allocatable, dimension(:, :, :, :) :: hrec_des_precomp
        real, allocatable, dimension(:, :, :, :) :: phimap_precomp
        real, allocatable, dimension(:, :, :)    :: abval_precomp
        real, allocatable, dimension(:, :, :, :) :: solgrd_precomp
        real, allocatable, dimension(:, :, :, :) :: hsolgrd_precomp
        real, allocatable, dimension(:, :, :, :) :: recdes_precomp
        real, allocatable, dimension(:, :, :, :) :: ligdes_precomp
      end type grids

      type flexgrids
        type(grids), allocatable, dimension(:) :: griddata
        integer :: group_len !how many groups
        integer :: total_combinations !how many combinations of the groups can be made
        integer, allocatable, dimension(:) :: group_part_len
        integer, allocatable, dimension(:, :) :: group_part !maps from each group
          !to a list of parts
        integer, dimension(:, :), allocatable :: all_combos
      end type flexgrids
    
      contains

      subroutine allocate_flexgrids(grids0, total_rec, OUTDOCK)

      type(flexgrids), intent(inout) :: grids0
      integer, intent(in) :: total_rec !total number of receptor grids
      integer, intent(in) :: OUTDOCK !output file in case of problems

      integer allocstat !temporary status check

      if (allocated(grids0%griddata)) then
        if (size(grids0%griddata) .lt. total_rec) then
          deallocate(grids0%griddata)
        endif
      endif
      if (.not. allocated(grids0%griddata)) then
        allocate(grids0%griddata(total_rec), stat=allocstat)
        if (allocstat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate griddata"
          stop
        endif
      endif
      if (allocated(grids0%group_part_len)) then
        if (size(grids0%group_part_len) .lt. total_rec) then
          deallocate(grids0%group_part_len)
        endif
      endif
      if (.not. allocated(grids0%group_part_len)) then
        allocate(grids0%group_part_len(total_rec), stat=allocstat)
        if (allocstat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate group_part_len"
          stop
        endif
      endif
      if (allocated(grids0%group_part)) then
        if ((size(grids0%group_part, 1) .lt. total_rec) .or.
     &      (size(grids0%group_part, 2) .lt. total_rec)) then
          deallocate(grids0%group_part)
        endif
      endif
      if (.not. allocated(grids0%group_part)) then
        allocate(grids0%group_part(total_rec, total_rec), 
     &      stat=allocstat)
        if (allocstat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate group_part"
          stop
        endif
      endif

      return
      end subroutine allocate_flexgrids

! logic for setting up the rest of the flexgrids type is here.
      subroutine setup_grouppart(grids0, group_numbers, part_numbers,
     &    total_rec, OUTDOCK)

      type(flexgrids), intent(inout) :: grids0
      integer, dimension(:), intent(in) :: group_numbers  
      integer, dimension(:), intent(in) :: part_numbers  
      integer, intent(in) :: total_rec !length of group&part
      integer, intent(in) :: OUTDOCK !output file in case of problems

      integer :: cur_rec, cur_group, cur_combo, cur_part !counters
      integer :: allocstat !for allocating status
      integer :: prev_grouplen !count number of already seen groups
      integer, dimension(:), allocatable :: allprev_group_len

      grids0%group_len = maxval(group_numbers) !how many groups?
      if (allocated(allprev_group_len)) then !always deallocate in case
        deallocate(allprev_group_len)
      endif
      allocate(allprev_group_len(grids0%group_len), 
     &    stat=allocstat)
      if (allocstat .ne. 0) then
        write(OUTDOCK, *)
     &      "ERROR ---> Failed to allocate allprev_group_len"
        stop
      endif
      do cur_rec = 1, total_rec
        grids0%group_part_len(cur_rec) = 0 !initialize all of them jic
      enddo
      do cur_rec = 1, total_rec
        grids0%group_part_len(group_numbers(cur_rec)) = 
     &      grids0%group_part_len(group_numbers(cur_rec)) + 1
        grids0%group_part(group_numbers(cur_rec), 
     &      grids0%group_part_len(group_numbers(cur_rec))) = 
     &      cur_rec
!        write(OUTDOCK,*) group_numbers(cur_rec), 
!     &      grids0%group_part_len(group_numbers(cur_rec)), "=",
!     &      cur_rec
      enddo
      !now setup the total_combinations
      grids0%total_combinations = 1 !always 1 to start
      prev_grouplen = grids0%total_combinations
      do cur_group = 1, grids0%group_len
        grids0%total_combinations = grids0%total_combinations *
     &      grids0%group_part_len(cur_group)
      enddo
      prev_grouplen = grids0%total_combinations
      do cur_group = 1, grids0%group_len
        prev_grouplen = prev_grouplen / grids0%group_part_len(cur_group)
        allprev_group_len(cur_group) = prev_grouplen
      enddo
c      write(OUTDOCK, *) "prev grouplen", allprev_group_len
c      write(OUTDOCK,*) "total combos", grids0%total_combinations
      if (allocated(grids0%all_combos)) then !always deallocate in case
        deallocate(grids0%all_combos)
      endif
      !allocate here so the arrays are the right length
      allocate(grids0%all_combos(grids0%total_combinations, 
     &    grids0%group_len), stat=allocstat)
      if (allocstat .ne. 0) then
        write(OUTDOCK, *)
     &      "ERROR ---> Failed to allocate all_combos"
        stop
      endif
      !actually setup all_combos now
      do cur_group = 1, grids0%group_len
        cur_part = 1
        do cur_combo = 1, grids0%total_combinations
          grids0%all_combos(cur_combo, cur_group) = 
     &        grids0%group_part(cur_group, cur_part)
          if (allprev_group_len(cur_group) .ne. 0) then
            if (mod(cur_combo, allprev_group_len(cur_group)) 
     &          .eq. 0) then
              cur_part = cur_part + 1
              if (cur_part .gt. grids0%group_part_len(cur_group)) then
                cur_part = 1
              endif
            endif
          endif
        enddo
      enddo

      return
      end subroutine setup_grouppart

!this gets a code, or a string representing which grids were used for each combo
      function getcode(grids0, which_combo) result(reccode)

      type(flexgrids), intent(in) :: grids0
      integer, intent(in) :: which_combo

      ! XXX: gfortran chokes on this. Claims it's not a DUMMY variable (TS)
      !character (len=255), intent(out) :: reccode
      character (len=255) :: reccode

      character (len=255) :: tempfull
      integer :: count !temporary counter
      character (len=10) :: tempcode

      reccode = trim(adjustl('')) !start with empty string
      do count = 1, grids0%group_len
        write (tempcode, '(i10)') grids0%all_combos(which_combo, count)
        tempfull = adjustl(trim(reccode))
        if (1 .eq. count) then
          reccode = adjustl(trim(tempcode))
        else
          reccode = tempfull(1:len_trim(tempfull)) // '.' // 
     &        adjustl(trim(tempcode))
        endif
      enddo

      return
      end function getcode

      end module gridstype
