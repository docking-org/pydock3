!contains information about a bunch of matches (matrices to rotate)
      module matchtype

      implicit none

      type matcht
        real :: dmaxl !ligand sphere diameter
        real :: dmaxr !receptor sphere diameter
        integer :: nmatch !how many matches we have currently, how full com,rot
        real, dimension(:, :), allocatable :: coml ! center of ligand mass, is this for translation ?
        real, dimension(:, :), allocatable :: comr ! center of receptor mass??
        real, dimension(:, :, :), allocatable :: rot ! this is for rotation
        real, dimension(:, :, :), allocatable :: rotab !temporary matrix to store rotation matrix
        real, dimension(:, :, :), allocatable :: rotbc !temporary matrix to store new rotation matrix after translation and rotation
        real, dimension(:, :), allocatable :: comrab ! temporary matrix
        real, dimension(:, :, :), allocatable :: xor ! coordinate storage, receptor spheres 
        real, dimension(:, :, :), allocatable :: xos ! coordinate storage, ligand spheres

c 
      end type matcht

      contains

!subroutine allocates the things declared here
      subroutine allocate_match(match, maxor, OUTDOCK)

      type(matcht), intent(inout) :: match
      integer, intent(in) :: maxor !how many to allocate maximum
      integer, intent(in) :: OUTDOCK

      integer alloc_stat




      if (allocated(match%comrab)) then
        if (size(match%comrab, 2) .lt. maxor) then
          deallocate(match%comrab)
        endif
      endif
      if (.not. allocated(match%comrab)) then
        allocate(match%comrab(3, maxor), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate comrab"
          stop
        endif
      endif
      if (allocated(match%rotab)) then
        if (size(match%rotab, 3) .lt. maxor) then
          deallocate(match%rotab)
        endif
      endif
      if (.not. allocated(match%rotab)) then
        allocate(match%rotab(3, 3, maxor), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate rotab"
          stop
        endif
      endif
      if (allocated(match%rotbc)) then
        if (size(match%rotbc, 3) .lt. maxor) then
          deallocate(match%rotbc)
        endif
      endif
      if (.not. allocated(match%rotbc)) then
        allocate(match%rotbc(3, 3, maxor), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate rotbc"
          stop
        endif
      endif
      if (allocated(match%coml)) then
        if (size(match%coml, 2) .lt. maxor) then
          deallocate(match%coml)
        endif
      endif
      if (.not. allocated(match%coml)) then
        allocate(match%coml(3, maxor), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate coml"
          stop
        endif
      endif
      if (allocated(match%comr)) then
        if (size(match%comr, 2) .lt. maxor) then
          deallocate(match%comr)
        endif
      endif
      if (.not. allocated(match%comr)) then
        allocate(match%comr(3, maxor), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate comr"
          stop
        endif
      endif
      if (allocated(match%rot)) then
        if (size(match%rot, 3) .lt. maxor) then
          deallocate(match%rot)
        endif
      endif
      if (.not. allocated(match%rot)) then
        allocate(match%rot(3, 3, maxor), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate rot"
          stop
        endif
      endif
      if (allocated(match%xor)) then
        if (size(match%xor, 3) .lt. maxor) then
          deallocate(match%xor)
        endif
      endif
      if (.not. allocated(match%xor)) then
        ! we need to replace 4 with a variable
        allocate(match%xor(3, 4, maxor), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate xor"
          stop
        endif
      endif
      if (allocated(match%xos)) then
        if (size(match%xos, 3) .lt. maxor) then
          deallocate(match%xos)
        endif
      endif
      if (.not. allocated(match%xos)) then
        ! we need to replace 4 with a variable
        allocate(match%xos(3, 4, maxor), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate xos"
          stop
        endif
      endif
      
      return
      end subroutine allocate_match

c     output match
      subroutine output_match(match, maxor, prefix)
      type(matcht), intent(inout) :: match
      integer, intent(in) :: maxor !how many to allocate maximum
      character (len=80), intent(in) :: prefix
      integer curmatch
      integer istat_comr
      integer istat_coml
      integer istat_rot
      integer istat
      integer (kind=8) :: fdcomr
      integer (kind=8) :: fdcoml
      integer (kind=8) :: fdrot
      character (len=*), parameter :: COMR = '(3f10.5)'
      character (len=*), parameter :: COML = '(3f10.5)'
      character (len=*), parameter :: ROT  = '(9f8.4)'
      character (len=80) buf
      character (len=80) filename_comr
      character (len=80) filename_coml
      character (len=80) filename_rot
      filename_comr = trim(prefix)//'_comr.gz'
      filename_coml = trim(prefix)//'_coml.gz'
      filename_rot  = trim(prefix)//'_rot.gz'
      call gzopen(fdcomr,'w',filename_comr,istat_comr)
      call gzopen(fdcoml,'w',filename_coml,istat_coml)
      call gzopen(fdrot,'w',filename_rot,istat_rot)
      if (match%nmatch .gt. maxor) then !adjust down since not all rotations save
        match%nmatch = maxor
      endif
      do curmatch = 1, match%nmatch !copy rotmatch into rot, etc
        write(buf,COMR) match%comr(1, curmatch),
     &             match%comr(2, curmatch),
     &             match%comr(3, curmatch)
        call gzwrite(fdcomr, buf, istat) 
        write(buf,COML) match%coml(1, curmatch),
     &             match%coml(2, curmatch),
     &             match%coml(3, curmatch)
        call gzwrite(fdcoml, buf, istat) 
        write(buf,ROT) match%rot(1, 1, curmatch),
     &             match%rot(2, 1, curmatch),
     &             match%rot(3, 1, curmatch),
     &             match%rot(1, 2, curmatch),
     &             match%rot(2, 2, curmatch),
     &             match%rot(3, 2, curmatch),
     &             match%rot(1, 3, curmatch),
     &             match%rot(2, 3, curmatch),
     &             match%rot(3, 3, curmatch)
        call gzwrite(fdrot, buf, istat) 
      enddo
      call gzclose(fdcomr, istat_comr)
      call gzclose(fdcoml, istat_coml)
      call gzclose(fdrot, istat_rot)
      end subroutine output_match

      end module matchtype

