!saves by atom scores
      module atomscoretype

      implicit none

      type atomscoret
        real, dimension(:), allocatable :: atom_score !overall energy
        real, dimension(:), allocatable :: atom_vas !attractive van der waals score
        real, dimension(:), allocatable :: atom_vrs !repulsive van der waals score
        real, dimension(:), allocatable :: atom_es !overall electrostatic score
        real, dimension(:), allocatable :: atom_ps !overall polar desolvation score
        real, dimension(:), allocatable :: atom_as !overall apolar desolvation score
        real, dimension(:), allocatable :: atom_ds !receptor desolvation energy of each saved pose
        real, dimension(:), allocatable :: atom_rd !receptor desolvation energy of each saved pose
        real, dimension(:), allocatable :: atom_hs !receptor hydrophobic effect
      end type atomscoret

      contains

!subroutine allocates the things declared here
      subroutine allocate_atomscore(atomsc, atomstosave, OUTDOCK)

      type(atomscoret), intent(inout) :: atomsc
      integer, intent(in) :: atomstosave
      integer, intent(in) :: OUTDOCK !for printing out if memory allocation problems

      integer alloc_stat

      if (allocated(atomsc%atom_score)) then
        deallocate(atomsc%atom_score)
      endif
      if (allocated(atomsc%atom_vas)) then
        deallocate(atomsc%atom_vas)
      endif
      if (allocated(atomsc%atom_vrs)) then
        deallocate(atomsc%atom_vrs)
      endif
      if (allocated(atomsc%atom_es)) then
        deallocate(atomsc%atom_es)
      endif
      if (allocated(atomsc%atom_as)) then
        deallocate(atomsc%atom_as)
      endif
      if (allocated(atomsc%atom_ps)) then
        deallocate(atomsc%atom_ps)
      endif
      if (allocated(atomsc%atom_ds)) then
        deallocate(atomsc%atom_ds)
      endif
      if (allocated(atomsc%atom_rd)) then
        deallocate(atomsc%atom_rd)
      endif
      if (allocated(atomsc%atom_hs)) then
        deallocate(atomsc%atom_hs)
      endif
      allocate(atomsc%atom_score(atomstosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate atom_score"
        stop
      endif
      allocate(atomsc%atom_vas(atomstosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate atom_vas"
        stop
      endif
      allocate(atomsc%atom_vrs(atomstosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate atom_vrs"
        stop
      endif
      allocate(atomsc%atom_es(atomstosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate atom_es"
        stop
      endif
      allocate(atomsc%atom_as(atomstosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate atom_as"
        stop
      endif
      allocate(atomsc%atom_ps(atomstosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate atom_ps"
        stop
      endif
      allocate(atomsc%atom_ds(atomstosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate atom_ds"
        stop
      endif
      allocate(atomsc%atom_rd(atomstosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate atom_rd"
        stop
      endif
      allocate(atomsc%atom_hs(atomstosave), stat=alloc_stat)
      if (alloc_stat .ne. 0) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate atom_hs"
        stop
      endif
      
      return
      end subroutine allocate_atomscore

      end module atomscoretype
