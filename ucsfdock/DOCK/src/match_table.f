! a simple look-up table of matches for different rigid fragments
! matchtable[rig_frag] = match
! jklyu, 20190612
      module match_table

      use matchtype

      implicit none

      type match_table_t
        !integer, parameter :: charlength = 128 ! char length of keys
        !character(len=charlength), dimension(:), allocatable :: keys
        character(len=128), dimension(:), allocatable :: keys
        type(matcht), dimension(:), allocatable :: matches
        integer :: index
c 
      end type match_table_t

      contains

!subroutine allocates matchtable
      subroutine allocate_matchtable(matchtable, maxmatch, OUTDOCK)
      type(match_table_t), intent(inout) :: matchtable
      integer, intent(in) :: maxmatch !how many to allocate maximum
      integer, intent(in) :: OUTDOCK

      integer alloc_stat

      if (allocated(matchtable%keys)) then
        if (size(matchtable%keys, 1) .lt. maxmatch) then
          deallocate(matchtable%keys)
        endif
      endif
      if (.not. allocated(matchtable%keys)) then
        allocate(matchtable%keys(maxmatch), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate matchtable-->keys"
          stop
        endif
      endif

      if (allocated(matchtable%matches)) then
        if (size(matchtable%matches, 1) .lt. maxmatch) then
          deallocate(matchtable%matches)
        endif
      endif
      if (.not. allocated(matchtable%matches)) then
        allocate(matchtable%matches(maxmatch), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate matchtable-->matches"
          stop
        endif

      matchtable%index = 0

      endif
      
      return
      end subroutine allocate_matchtable

c     read in a match from files
      subroutine match_read(matchtable, maxor, OUTDOCK,
     &     rig_frag_name, match_index)

      type(match_table_t), intent(inout) :: matchtable
      character (len=128), intent(in) :: rig_frag_name
      integer, intent(inout) :: match_index
      integer, intent(in) :: maxor !how many to allocate maximum
      integer, intent(in) :: OUTDOCK

      integer alloc_stat
      integer (kind=8) :: fdcomr
      integer (kind=8) :: fdcoml
      integer (kind=8) :: fdrot
      integer istat_comr
      integer istat_coml
      integer istat_rot
      character (len=*), parameter :: COMR = '(3f10.5)'
      character (len=*), parameter :: COML = '(3f10.5)'
      character (len=*), parameter :: ROT  = '(9f8.4)'
      character (len=80) line_comr
      character (len=80) line_coml
      character (len=80) line_rot
      character (len=80) filename_comr
      character (len=80) filename_coml
      character (len=80) filename_rot
      integer curmatch
      logical read_okay
      integer read_stat

      filename_comr = trim(rig_frag_name)//'_comr.gz'
      filename_coml = trim(rig_frag_name)//'_coml.gz'
      filename_rot  = trim(rig_frag_name)//'_rot.gz'

      matchtable%index = matchtable%index + 1
      ! initialize match
      call allocate_match(matchtable%matches(matchtable%index),
     &      maxor, OUTDOCK)
      ! read in from files
c      write(OUTDOCK, *) filename_comr
c      write(OUTDOCK, *) filename_coml
c      write(OUTDOCK, *) filename_rot
      fdcomr = 104
      fdcoml = 105
      fdrot = 106
      call gzopen(fdcomr,'r',filename_comr,istat_comr)
      call gzopen(fdcoml,'r',filename_coml,istat_coml)
      call gzopen(fdrot,'r',filename_rot,istat_rot)
c      write(OUTDOCK, *) istat_comr
c      write(OUTDOCK, *) istat_coml
c      write(OUTDOCK, *) istat_rot
      if (istat_comr .ne. 0) then !problem opening files
         write(OUTDOCK, *) 'Error opening orientation file'
         stop
      endif
      if (istat_coml .ne. 0) then !problem opening files
         write(OUTDOCK, *) 'Error opening orientation file'
         stop
      endif
      if (istat_rot .ne. 0) then !problem opening files
         write(OUTDOCK, *) 'Error opening orientation file'
         stop
      endif
      curmatch = 0
      read_okay = .true.
      do while (read_okay)
        call gzread(fdcomr, line_comr, 80, istat_comr)
        call gzread(fdcoml, line_coml, 80, istat_coml)
        call gzread(fdrot, line_rot, 80, istat_rot)
c        write(OUTDOCK, *) line_rot
        if (istat_comr .le. -1) then 
          !inputend = .true.
          exit !break out of do loop
        endif
        if (istat_coml .le. -1) then !the end of file has been reached
          !inputend = .true.
          exit !break out of do loop
        endif
        if (istat_rot .le. -1) then !the end of file has been reached
          !inputend = .true.
          exit !break out of do loop
        endif
        curmatch = curmatch + 1
        read(line_comr,COMR,iostat=read_stat) 
     &    matchtable%matches(matchtable%index)%comr(1, curmatch),
     &    matchtable%matches(matchtable%index)%comr(2, curmatch),
     &    matchtable%matches(matchtable%index)%comr(3, curmatch)
        read(line_coml,COML,iostat=read_stat)
     &    matchtable%matches(matchtable%index)%coml(1, curmatch),
     &    matchtable%matches(matchtable%index)%coml(2, curmatch),
     &    matchtable%matches(matchtable%index)%coml(3, curmatch)
        read(line_rot,ROT,iostat=read_stat)
     &    matchtable%matches(matchtable%index)%rot(1, 1, curmatch),
     &    matchtable%matches(matchtable%index)%rot(2, 1, curmatch),
     &    matchtable%matches(matchtable%index)%rot(3, 1, curmatch),
     &    matchtable%matches(matchtable%index)%rot(1, 2, curmatch),
     &    matchtable%matches(matchtable%index)%rot(2, 2, curmatch),
     &    matchtable%matches(matchtable%index)%rot(3, 2, curmatch),
     &    matchtable%matches(matchtable%index)%rot(1, 3, curmatch),
     &    matchtable%matches(matchtable%index)%rot(2, 3, curmatch),
     &    matchtable%matches(matchtable%index)%rot(3, 3, curmatch)
c        write(OUTDOCK,ROT)
c     &    matchtable%matches(matchtable%index)%rot(1, 1, curmatch),
c     &    matchtable%matches(matchtable%index)%rot(2, 1, curmatch),
c     &    matchtable%matches(matchtable%index)%rot(3, 1, curmatch),
c     &    matchtable%matches(matchtable%index)%rot(1, 2, curmatch),
c     &    matchtable%matches(matchtable%index)%rot(2, 2, curmatch),
c     &    matchtable%matches(matchtable%index)%rot(3, 2, curmatch),
c     &    matchtable%matches(matchtable%index)%rot(1, 3, curmatch),
c     &    matchtable%matches(matchtable%index)%rot(2, 3, curmatch),
c     &    matchtable%matches(matchtable%index)%rot(3, 3, curmatch)
        if (read_stat .ne. 0) then
          write(OUTDOCK, *)
     &         'The orientation matrices is broken, failed to be read'
          read_okay = .false.
          exit !break out of do loop
        endif 
      enddo
      if (read_okay) then
        match_index = matchtable%index
        matchtable%matches(matchtable%index)%nmatch = curmatch
        matchtable%keys(matchtable%index) = rig_frag_name
      else
        write(OUTDOCK, *)
     &          'failed to read orientation files'
        stop
      endif
      call gzclose(fdcomr, istat_comr)
      call gzclose(fdcoml, istat_coml)
      call gzclose(fdrot, istat_rot)

c  980 write(OUTDOCK, *) 'Error opening orientation file'
c      stop

      end subroutine match_read


c     find a match based on the rigid fragment name
      subroutine match_get(matchtable, maxor, OUTDOCK,
     &     rig_frag_name, match_index)
      type(match_table_t), intent(inout) :: matchtable
      character (len=128), intent(in) :: rig_frag_name
      integer, intent(inout) :: match_index
      integer, intent(in) :: maxor !how many to allocate maximum
      integer, intent(in) :: OUTDOCK

      integer :: local_index
      logical :: found
      found = .false.
c      do local_index = 1, size(matchtable%keys,1)
      do local_index = 1, matchtable%index
        if(TRIM(matchtable%keys(local_index))
     &    == TRIM(rig_frag_name)) then
          match_index = local_index
c          write(OUTDOCK, *) 'found'
          found = .true.
          exit
        endif
      enddo
      if (.not. found) then
        !match_read
c        write(OUTDOCK, *) 'read in new orientation'
        call match_read(matchtable, maxor, OUTDOCK,
     &     rig_frag_name, match_index)
      endif

      return

      end subroutine match_get
      



      end module match_table

