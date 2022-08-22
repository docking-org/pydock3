! everything about the ligand in db2 format read into these data structures
      module db2type
  
      implicit none
      integer, parameter :: MAXADDMATCH = 5 !bad bad bad but necessary for now
      !would need to modify ligand output to write this out which isn't done 
      !yet
      type db2
        character (len=16) :: refcod !reference or zinc code
        !character (len=20) :: refcod !reference or zinc code
        character (len=9) :: protcod !protonation or other additional id
        integer :: total_atoms !number of atoms in ligand
        integer :: total_bonds !number of bonds in ligand
        integer :: total_coords !number of coordinates in the ensemble/ligand
        integer :: total_confs !number of conformations of all groups in ligand
        integer :: total_sets !number of possible combinations of conformations
        real :: total_charge !total charge of ligand
        real :: total_polar_solv !total polar solvation of ligand
        real :: total_apolar_solv !total apolar solvation of ligand
        real :: total_solv !total solvation of ligand
        real :: total_surfarea !total surface area of ligand
        real :: ligcharge !total ligand charge
        integer :: total_heavy_atoms !total number of heavy atoms
        character (len=78) :: molname        ! name of molecule
        character (len=78) :: smiles  !smiles string of ligand
        ! jklyu, added on Aug. 26th, 2019
        character (len=78) :: rig_frag_code ! name of rigid fragment
        character (len=78), dimension(:), allocatable :: arbitrary 
               !any arbitrary information, currently 20 lines maximum
        integer :: mlines !the number of mlines specified by the ligand input db2
        integer, dimension(:), allocatable :: atom_num !atom number (orig numbering)
        character (len=4), dimension(:), allocatable :: atom_name !atom name (C1, C2, O1, etc)
        character (len=5), dimension(:), allocatable :: atom_type !atom type (C.3, N.2, etc) from sybyl
        integer, dimension(:), allocatable :: atom_vdwtype !atom dock vdw type
        integer, dimension(:), allocatable :: atom_color !atom color number
        real, dimension(:), allocatable :: atom_charge !atom partial charge
        real, dimension(:), allocatable :: atom_polsolv !atom polar solvation
        real, dimension(:), allocatable :: atom_apolsolv !atom apolar solvation
        real, dimension(:), allocatable :: atom_totalsolv !atom total solvation
        real, dimension(:), allocatable :: atom_surfarea !atom surface area
        integer, dimension(:), allocatable :: atom_heavy_valence !heavy valence
        !is the number of heavy atoms bonded to each atom
        integer, dimension(:), allocatable :: bond_num !bond number
        integer, dimension(:), allocatable :: bond_start !bond atom (arbitrary, left side)
        integer, dimension(:), allocatable :: bond_end !bond atom (arbitrary, right side)
        character (len=2), dimension(:), allocatable :: bond_type !bond type (1, 2, am, ar)
        !end of bonds & atoms stuff, now hierarchy stuff
        integer, dimension(:, :), allocatable :: coord_index !coord#, atom#, conf# for each coordinate
        real, dimension(:, :), allocatable :: coords!x,y,z for each atom. same order as coord_index.
        real, dimension(:, :), allocatable :: transfm_coords! coordinates transformed into real space
        real, dimension(:, :, :), allocatable :: save_cov_coords! coordinates transformed into real space
        !sets of complete conformations, length
        integer, dimension(:), allocatable :: set_conf_len
        integer, dimension(:), allocatable :: conf_conection ! this will discribe what each conf ( a segment of the molecule) is conected for each set.
        integer, dimension(:), allocatable :: conf_sorted ! gist need the confs to be in a particular order.
        integer, dimension(:, :), allocatable :: set_conf !sets->conf number 
        integer, dimension(:), allocatable :: set_broken !sets->broken
        integer, dimension(:), allocatable :: set_above_strain !sets->above_strain
        real, dimension(:), allocatable :: set_energy !sets-> energy (intern)
        real, dimension(:), allocatable :: set_total_strain !sets-> total_strain
        real, dimension(:), allocatable :: set_max_strain !sets-> max_strain
        integer, dimension(:, :), allocatable :: conf_coord !index into coord 
          !numbers, 1 is start, 2 is end
        !here are the things that read in clusters
        integer :: total_clusters !how many clusters are in file (30 max presently)
        integer, dimension(:), allocatable :: cluster_set_start !points to start of set for each cluster
        integer, dimension(:), allocatable :: cluster_set_end !points to end of set for each cluster
        integer, dimension(:), allocatable :: cluster_match_start !start of add matching spheres
        integer, dimension(:), allocatable :: cluster_match_end !end of add matching spheres
        integer, dimension(:), allocatable :: addmatch_color !color of each matching sphere
        real, dimension(:, :), allocatable :: addmatch_coord !xyz of each matching sphere, x,y,z is first (1,2,3)
      end type db2
     
      contains
  
!this  allocates data for the ligand parts stored here
      subroutine allocate_db2(db2lig, OUTDOCK, total_sets, 
     &    total_coords, total_atoms, total_bonds, total_confs, 
     &    total_mlines, total_clusts, dockovalent, nsav)

      type(db2), intent(inout) :: db2lig
      integer, intent(in) :: OUTDOCK
      integer, intent(in) :: total_sets !for this ligand
      integer, intent(in) :: total_coords !for this ligand
      integer, intent(in) :: total_atoms !for this ligand
      integer, intent(in) :: total_bonds !for this ligand
      integer, intent(in) :: total_confs !for this ligand
      integer, intent(in) :: total_mlines !for this ligand
      integer, intent(in) :: total_clusts !for this ligand
      logical, intent(in) :: dockovalent !to check if this is a covalent run
      integer, intent(in) :: nsav ! how many poses to save for covalent docking

      integer alloc_stat !status check upon allocation

c      write(OUTDOCK, *) "allocating ligand", total_sets, total_coords,
c     &    total_atoms, total_confs
c      call doflush(OUTDOCK)

!strategy is to deallocate if the last ligand didn't have enough things, 
! then if deallocated, allocate anew
      if (allocated(db2lig%set_conf_len)) then
        if (size(db2lig%set_conf_len) .lt. total_sets) then
          deallocate(db2lig%set_conf_len)
        endif
      endif
      if (.not. allocated(db2lig%set_conf_len)) then
        allocate(db2lig%set_conf_len(total_sets), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate set_conf_len"
          stop
        endif
      endif

      if (allocated(db2lig%conf_conection)) then
        if (size(db2lig%conf_conection) .lt. total_confs) then 
          deallocate(db2lig%conf_conection)
        endif
      endif
      if (.not. allocated(db2lig%conf_conection)) then
        allocate(db2lig%conf_conection(total_confs), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate conf_conection"
          stop
        endif
      endif

      if (allocated(db2lig%conf_sorted)) then
        if (size(db2lig%conf_sorted) .lt. total_confs) then
          deallocate(db2lig%conf_sorted)
        endif
      endif
      if (.not. allocated(db2lig%conf_sorted)) then
        allocate(db2lig%conf_sorted(total_confs), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate conf_sorted"
          stop
        endif
      endif


      if (allocated(db2lig%set_conf)) then
        if ((size(db2lig%set_conf, 1) .lt. total_sets) .or.
     &      (size(db2lig%set_conf, 2) .lt. total_atoms)) then
          deallocate(db2lig%set_conf)
        endif
      endif
      if (.not. allocated(db2lig%set_conf)) then
        allocate(db2lig%set_conf(total_sets, total_atoms), 
     &      stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate set_conf"
          stop
        endif
      endif
      if (allocated(db2lig%conf_coord)) then
        if (size(db2lig%conf_coord, 2) .lt. total_confs) then
          deallocate(db2lig%conf_coord)
        endif
      endif
      if (.not. allocated(db2lig%conf_coord)) then
        allocate(db2lig%conf_coord(2, total_confs), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate conf_coord"
          stop
        endif
      endif
      if (allocated(db2lig%arbitrary)) then
        if (size(db2lig%arbitrary) .lt. total_mlines) then
          deallocate(db2lig%arbitrary)
        endif
      endif
      if (.not. allocated(db2lig%arbitrary)) then
        allocate(db2lig%arbitrary(total_mlines), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate arbitrary"
          stop
        endif
      endif
      if (allocated(db2lig%set_broken)) then
        if (size(db2lig%set_broken) .lt. total_sets) then
          deallocate(db2lig%set_broken)
        endif
      endif
      if (.not. allocated(db2lig%set_broken)) then
        allocate(db2lig%set_broken(total_sets), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate set_broken"
          stop
        endif
      endif
      if (allocated(db2lig%set_above_strain)) then
        if (size(db2lig%set_above_strain) .lt. total_sets) then
          deallocate(db2lig%set_above_strain)
        endif
      endif
      if (.not. allocated(db2lig%set_above_strain)) then
        allocate(db2lig%set_above_strain(total_sets), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate set_above_strain"
          stop
        endif
      endif
      if (allocated(db2lig%set_energy)) then
        if (size(db2lig%set_energy) .lt. total_sets) then
          deallocate(db2lig%set_energy)
        endif
      endif
      if (.not. allocated(db2lig%set_energy)) then
        allocate(db2lig%set_energy(total_sets), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate set_energy"
          stop
        endif
      endif
      if (allocated(db2lig%set_total_strain)) then
        if (size(db2lig%set_total_strain) .lt. total_sets) then
          deallocate(db2lig%set_total_strain)
        endif
      endif
      if (.not. allocated(db2lig%set_total_strain)) then
        allocate(db2lig%set_total_strain(total_sets), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate set_total_strain"
          stop
        endif
      endif
      if (allocated(db2lig%set_max_strain)) then
        if (size(db2lig%set_max_strain) .lt. total_sets) then
          deallocate(db2lig%set_max_strain)
        endif
      endif
      if (.not. allocated(db2lig%set_max_strain)) then
        allocate(db2lig%set_max_strain(total_sets), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate set_max_strain"
          stop
        endif
      endif
      if (allocated(db2lig%atom_num)) then
        if (size(db2lig%atom_num) .lt. total_atoms) then
          deallocate(db2lig%atom_num)
        endif
      endif
      if (.not. allocated(db2lig%atom_num)) then
        allocate(db2lig%atom_num(total_atoms), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate atom_num"
          stop
        endif
      endif
      if (allocated(db2lig%atom_name)) then
        if (size(db2lig%atom_name) .lt. total_atoms) then
          deallocate(db2lig%atom_name)
        endif
      endif
      if (.not. allocated(db2lig%atom_name)) then
        allocate(db2lig%atom_name(total_atoms), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate atom_name"
          stop
        endif
      endif
      if (allocated(db2lig%atom_type)) then
        if (size(db2lig%atom_type) .lt. total_atoms) then
          deallocate(db2lig%atom_type)
        endif
      endif
      if (.not. allocated(db2lig%atom_type)) then
        allocate(db2lig%atom_type(total_atoms), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate atom_type"
          stop
        endif
      endif
      if (allocated(db2lig%atom_vdwtype)) then
        if (size(db2lig%atom_vdwtype) .lt. total_atoms) then
          deallocate(db2lig%atom_vdwtype)
        endif
      endif
      if (.not. allocated(db2lig%atom_vdwtype)) then
        allocate(db2lig%atom_vdwtype(total_atoms), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate atom_vdwtype"
          stop
        endif
      endif
      if (allocated(db2lig%atom_color)) then
        if (size(db2lig%atom_color) .lt. total_atoms) then
          deallocate(db2lig%atom_color)
        endif
      endif
      if (.not. allocated(db2lig%atom_color)) then
        allocate(db2lig%atom_color(total_atoms), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate atom_color"
          stop
        endif
      endif
      if (allocated(db2lig%atom_charge)) then
        if (size(db2lig%atom_charge) .lt. total_atoms) then
          deallocate(db2lig%atom_charge)
        endif
      endif
      if (.not. allocated(db2lig%atom_charge)) then
        allocate(db2lig%atom_charge(total_atoms), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate atom_charge"
          stop
        endif
      endif
      if (allocated(db2lig%atom_heavy_valence)) then
        if (size(db2lig%atom_heavy_valence) .lt. total_atoms) then
          deallocate(db2lig%atom_heavy_valence)
        endif
      endif
      if (.not. allocated(db2lig%atom_heavy_valence)) then
        allocate(db2lig%atom_heavy_valence(total_atoms), 
     &      stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to allocate atom_heavy_valence"
          stop
        endif
      endif
      if (allocated(db2lig%atom_polsolv)) then
        if (size(db2lig%atom_polsolv) .lt. total_atoms) then
          deallocate(db2lig%atom_polsolv)
        endif
      endif
      if (.not. allocated(db2lig%atom_polsolv)) then
        allocate(db2lig%atom_polsolv(total_atoms), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate atom_polsolv"
          stop
        endif
      endif
      if (allocated(db2lig%atom_apolsolv)) then
        if (size(db2lig%atom_apolsolv) .lt. total_atoms) then
          deallocate(db2lig%atom_apolsolv)
        endif
      endif
      if (.not. allocated(db2lig%atom_apolsolv)) then
        allocate(db2lig%atom_apolsolv(total_atoms), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate atom_apolsolv"
          stop
        endif
      endif
      if (allocated(db2lig%atom_totalsolv)) then
        if (size(db2lig%atom_totalsolv) .lt. total_atoms) then
          deallocate(db2lig%atom_totalsolv)
        endif
      endif
      if (.not. allocated(db2lig%atom_totalsolv)) then
        allocate(db2lig%atom_totalsolv(total_atoms), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate atom_totalsolv"
          stop
        endif
      endif
      if (allocated(db2lig%atom_surfarea)) then
        if (size(db2lig%atom_surfarea) .lt. total_atoms) then
          deallocate(db2lig%atom_surfarea)
        endif
      endif
      if (.not. allocated(db2lig%atom_surfarea)) then
        allocate(db2lig%atom_surfarea(total_atoms), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate atom_surfarea"
          stop
        endif
      endif
      if (allocated(db2lig%bond_num)) then
        if (size(db2lig%bond_num) .lt. total_bonds) then
          deallocate(db2lig%bond_num)
        endif
      endif
      if (.not. allocated(db2lig%bond_num)) then
        allocate(db2lig%bond_num(total_bonds), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate bond_num"
          stop
        endif
      endif
      if (allocated(db2lig%bond_start)) then
        if (size(db2lig%bond_start) .lt. total_bonds) then
          deallocate(db2lig%bond_start)
        endif
      endif
      if (.not. allocated(db2lig%bond_start)) then
        allocate(db2lig%bond_start(total_bonds), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate bond_start"
          stop
        endif
      endif
      if (allocated(db2lig%bond_end)) then
        if (size(db2lig%bond_end) .lt. total_bonds) then
          deallocate(db2lig%bond_end)
        endif
      endif
      if (.not. allocated(db2lig%bond_end)) then
        allocate(db2lig%bond_end(total_bonds), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate bond_end"
          stop
        endif
      endif
      if (allocated(db2lig%bond_type)) then
        if (size(db2lig%bond_type) .lt. total_bonds) then
          deallocate(db2lig%bond_type)
        endif
      endif
      if (.not. allocated(db2lig%bond_type)) then
        allocate(db2lig%bond_type(total_bonds), stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate bond_type"
          stop
        endif
      endif
      if (allocated(db2lig%coord_index)) then
        if (size(db2lig%coord_index, 2) .lt. total_coords) then
          deallocate(db2lig%coord_index)
        endif
      endif
      if (.not. allocated(db2lig%coord_index)) then
        allocate(db2lig%coord_index(3, total_coords), 
     &      stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate coord_index"
          stop
        endif
      endif
      if (allocated(db2lig%coords)) then
        if (size(db2lig%coords, 2) .lt. total_coords) then
          deallocate(db2lig%coords)
        endif
      endif
      if (.not. allocated(db2lig%coords)) then
        allocate(db2lig%coords(3, total_coords), 
     &      stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate coords"
          stop
        endif
      endif
      if (allocated(db2lig%transfm_coords)) then
        if (size(db2lig%transfm_coords, 2) .lt. total_atoms) then
          deallocate(db2lig%transfm_coords)
        endif
      endif
      if (.not. allocated(db2lig%transfm_coords)) then
        allocate(db2lig%transfm_coords(3, total_atoms),
     &      stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate transfm_coords"
          stop
        endif
      endif

c######################################################      
      if (dockovalent .eqv. .true.) then
         if (allocated(db2lig%save_cov_coords)) then
            if (size(db2lig%save_cov_coords, 2) .lt. total_atoms) then
               deallocate(db2lig%save_cov_coords)
            endif
         endif
         if (.not. allocated(db2lig%save_cov_coords)) then
            allocate(db2lig%save_cov_coords(nsav, 3, total_atoms),
     &           stat=alloc_stat)
            if (alloc_stat .ne. 0) then
               write(OUTDOCK, *)
     &              "ERROR ---> Failed to  allocate save_cov_coords"
               stop
            endif
         endif
      endif
c######################################################

      if (allocated(db2lig%cluster_set_start)) then
        if (size(db2lig%cluster_set_start) .lt. total_clusts) then
          deallocate(db2lig%cluster_set_start)
        endif
      endif
      if (.not. allocated(db2lig%cluster_set_start)) then
        allocate(db2lig%cluster_set_start(total_clusts),
     &      stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate cluster_set_start"
          stop
        endif
      endif
      if (allocated(db2lig%cluster_set_end)) then
        if (size(db2lig%cluster_set_end) .lt. total_clusts) then
          deallocate(db2lig%cluster_set_end)
        endif
      endif
      if (.not. allocated(db2lig%cluster_set_end)) then
        allocate(db2lig%cluster_set_end(total_clusts),
     &      stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate cluster_set_end"
          stop
        endif
      endif
      if (allocated(db2lig%cluster_match_start)) then
        if (size(db2lig%cluster_match_start) .lt. total_clusts) then
          deallocate(db2lig%cluster_match_start)
        endif
      endif
      if (.not. allocated(db2lig%cluster_match_start)) then
        allocate(db2lig%cluster_match_start(total_clusts),
     &      stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate cluster_match_start"
          stop
        endif
      endif
      if (allocated(db2lig%cluster_match_end)) then
        if (size(db2lig%cluster_match_end) .lt. total_clusts) then
          deallocate(db2lig%cluster_match_end)
        endif
      endif
      if (.not. allocated(db2lig%cluster_match_end)) then
        allocate(db2lig%cluster_match_end(total_clusts),
     &      stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate cluster_match_end"
          stop
        endif
      endif
      if (allocated(db2lig%addmatch_color)) then
        if (size(db2lig%addmatch_color) .lt. 
     &      total_clusts * MAXADDMATCH) then
          deallocate(db2lig%addmatch_color)
        endif
      endif
      if (.not. allocated(db2lig%addmatch_color)) then
        allocate(db2lig%addmatch_color(total_clusts * MAXADDMATCH),
     &      stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate addmatch_color"
          stop
        endif
      endif
      if (allocated(db2lig%addmatch_coord)) then
        if (size(db2lig%addmatch_coord, 2) .lt. 
     &      total_clusts * MAXADDMATCH) then
          deallocate(db2lig%addmatch_coord)
        endif
      endif
      if (.not. allocated(db2lig%addmatch_coord)) then
        allocate(db2lig%addmatch_coord(3, total_clusts * MAXADDMATCH),
     &      stat=alloc_stat)
        if (alloc_stat .ne. 0) then
          write(OUTDOCK, *)
     &        "ERROR ---> Failed to  allocate addmatch_coord"
          stop
        endif
      endif

      return
      end subroutine allocate_db2

      end module db2type
