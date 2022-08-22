! all options that are controlled by INDOCK and other parameters
! stored here
      module optionstype
 
      implicit none
      type options
        logical :: search_flag, rescore_flag ! search or single point
        character (len=255) :: mol2inputfile ! mol2 input for rescoring
        character (len=255) :: ligsolinputfile ! solvation input for rescoring
        character (len=255) :: ligvdwinputfile ! vdw input for rescoring
        logical :: hsolflag !whether or not separate hydrogen solvation used
        integer :: ldesolv !ligand desolvation  0=none, 1=full, 2=context-dependent, 3=qnifft desolvation
        real :: bmpvdw, bmpvdwrigid !max vdw allowed for any group
        real :: bmptotal !max total energy allowed for any group
        real :: total_strain !total strain allowed for any conformer
        real :: max_strain !max strain per torsion allowed for any conformer
        integer :: nsav !how many poses to save for each ligand
        integer :: nwrite !how many poses to write for each ligand
        integer :: input_flush_int !how often should we flush to output file.  Ever flush_int lines we will write to output. 
        real :: timeout  !how long to run for any single ligand
        integer :: maxnodes !for matching, what is the max# of nodes (usually 4)
        integer :: minnodes !for matching, what is the min# of nodes (usually 4, sometimes 3)
        integer :: match_method !1 = just use dislim&timeout
                                !2 = use disstep, matchgoal
                                !timeout and dismax
        real :: dislim !maximum distance tolerance allowed between sph pairs
        real :: disstep !how much to add at a time if match_method = 2
        real :: dismax !when to stop if match_method = 2
        integer :: matchgoal !max number of matches to find if match_method = 2
        logical :: check_clash !check clashes if true
        logical :: check_strain !check strain if true
        logical :: cmatch !use chemical matching if true
        real :: escale !number to scale the electrostatics by
        real :: gistscale !number to scale the gist energy by
        real :: gistvolume !user specified gist-voxel volume 
        integer :: gistH ! !if 0 use hydrogens as is, if 1 do not use hydrogens in gist calculations, if 2 replace nonpolar radii with 0.8
        real :: vscale !number to scale the vdw energy by
        real :: iscale !number to scale the internal energy by
        real :: solvscale !number to scale the ligand desolvation terms by
        real :: rescale !number to scale the receptor energy term by
        real :: rdscale !number to scale the receptor desolvation term by
        real :: rec_des_scale !number to scale the receptor desolvation
        real :: rhscale !number to scale the receptor hydrophobic effect term by
        character (len=80) :: clufil !receptor spheres file (match2.sph)
        integer :: natmin, natmax !min & max numbers of atoms in ligands
        logical :: do_premax !if true, experimental premax scoring code will be used
        logical :: do_clusters !if true, clusters/clouds/additional matching spheres will be used
        logical :: astar !if true, experimental rigid first a* search will be used
        character (len=80) :: vdwfil !file containing vdw parameters
        logical :: casen !casen: case sensitive (regarding color names)
        character (len=255) :: outfil !filename of output files (usually test)
        character (len=255) :: ligfil !filename of input ligand file
        character (len=80) :: ligname !store ligname from the 2nd column
        ! of sdi_rig_frag files
        integer :: pos !store the lig position from the 3rd column of
        ! sdi_rig_frag files
        character (len=80) :: rig_frag_name !store the rigid fragment
        ! name from the 4th column of sdi_rig_frag files
        logical :: sdi !split database index, true if being used
        logical :: use_gist ! used in the ligread2 for sorting confs,
                            ! gist expects the confs to be in a certain
                            ! order
        character (len=255) :: use_gist_val ! string the should be yes or no.  
        logical :: per_atom_scores !true if you want extra output that includes per atom scores for each ligand
        logical :: receptor_desolv !true if you want to read in and score using receptor desolvation grids
        logical :: rec_des_r,hrec_des_r ! true if you want to read in REC_DESOLV 
        logical, allocatable, dimension(:) :: rec_d_flag, hrec_d_flag
        logical, allocatable, dimension(:) :: gistflag !true if you want to read in and score using receptor desolvation GIST grids
        !logical :: gist_aprox !true if you want to use an aproxamation of the gist calucluation. per compute using atom volumn. 
        integer, allocatable, dimension(:) :: gist_aprox !0 (no aprox) 1 (uses closest point), 2 uses clostest 8 points. 
        logical :: recdes_valence_penalty !true if you want to correct for valence of atoms
          !each nearby heavy atom is penalized by a certain amount
        real :: save_limit !poses with scores higher than this aren't written
          !use with prospective screens, not with retrospective screens as this
          !can change your results
! covalent docking options here
        logical :: dockovalent
        real :: bondlen
        real :: bondang1
        real :: bondang2
        real :: bondrange
        real :: dihstep
        real :: bondstep
        real :: bondang1range
        real :: bondang1step
        real :: bondang2range
        real :: bondang2step
! docking options here
        integer :: docking_mode

! flexible receptor options here
        logical :: flexible_receptor !whether or not we're using flexrec code
        integer :: total_receptors
        logical :: score_each_flex !true if you want a score for each flexible receptor, false otherwise
        integer :: current_number !which number we are currently reading in
        integer, allocatable, dimension(:) :: rec_numbers !number or ID of each set
        integer, allocatable, dimension(:) :: group_numbers !which group each belongs
        integer, allocatable, dimension(:) :: part_numbers !group_option, i.e. which
            !subpart this is, there should be part number 1,2,etc for each group
        real, allocatable, dimension(:) :: rec_energy !the energy associated with this
            !conformation.
        character (len=255), allocatable, dimension(:) :: solvnames !solvmap names
        character (len=255), allocatable, dimension(:) :: hsolvnames !hsolvmap names
        character (len=255), allocatable, dimension(:) :: phinames !phimap names
        character (len=255), allocatable, dimension(:) :: gistnames !gist names
        character (len=255), allocatable, dimension(:) :: recdesnames !receptor desolvation phimap names
        character (len=255), allocatable, dimension(:) :: recdesrnames 
        character (len=255), allocatable, dimension(:) :: hrecdesnames 
        character (len=255), allocatable, dimension(:) :: vdwnames !vdw grid names
        character (len=255), allocatable, dimension(:) :: bmpnames !bmp map names
            !bump map names is currently only used to read in the vdw grid
            !options like corners, grids per angstrom, etc, since apparently
            !this isn't contained in the vdw grid itself

c       minimization options are here
c       minimize--whether to energy minimize orientation force-field score
c       minimize: 0: none,  1: best mol only
        logical  allocated_seeds
        integer, allocatable :: seeds(:)
        integer :: mc_minimize ! monte carlo minimizer
        integer :: minimize
        integer :: mol2_minimizer
        integer :: output_premin
        integer :: num_mini_1
        integer :: total_iterations ! cumulative number of minimization steps performed
        logical :: check_degeneracy
        integer :: sim_itmax, mc_itmax ! maximum number of simplex/Monte Carlo iterations to perform
        integer :: mc_accpt ! maximum number of Monte Carlo steps that are accepted
        real :: sim_trnstep, mc_trnstep ! maximum stepsize for x, y, and z translation
        real :: sim_rotstep, mc_rotstep ! maximum stepsize for alpha, beta, and gamma auler angles
        real :: mc_temp ! temperature for Monte Carlo
        real :: sim_need_to_restart ! minimum energy change to signal a restart
        real :: sim_cnvrge !  energy convergence tolerance for simplex completion
        real :: min_cut !
        integer :: nvar ! number of variables to minimize
        integer :: mc_iseed
        real :: alpha ! simplex aular angles
        real :: beta ! simplex aular angles
        real :: gamma ! simplex aular angles
        real :: big_energy  ! large energy to initialize finding minimum energy
        integer :: iseed !random seed here
c        logical :: run_pos_solv ! positive polar desolvation will be redistributed
c       logical :: check_clash ! internal distance clashes will be checked and fixed if possible
        integer :: k_clusters ! number of spheres clusters to read in [default 1]  


      end type options

      contains

! this subroutine allocates all the allocatable things about flexible receptors
      subroutine allocate_flexible(options0, OUTDOCK)

      type(options), intent(inout) :: options0 !options for this run
      integer, intent(in) :: OUTDOCK !output file unit, write out any errors

      integer allocate_stat !temporary, used for making sure we have enough memory
      integer idx 
 
      allocate(options0%rec_numbers(options0%total_receptors), 
     &    stat=allocate_stat)
      if (.not. allocated(options0%rec_numbers)) then
        write(OUTDOCK, *) 
     &    "ERROR ---> Failed to dynamically allocate rec_numbers"
        stop
      endif
      allocate(options0%group_numbers(options0%total_receptors), 
     &    stat=allocate_stat)
      if (.not. allocated(options0%group_numbers)) then
        write(OUTDOCK, *) 
     &    "ERROR ---> Failed to dynamically allocate group_numbers"
        stop
      endif
      allocate(options0%part_numbers(options0%total_receptors), 
     &    stat=allocate_stat)
      if (.not. allocated(options0%part_numbers)) then
        write(OUTDOCK, *) 
     &    "ERROR ---> Failed to dynamically allocate part_numbers"
        stop
      endif
      allocate(options0%rec_energy(options0%total_receptors), 
     &    stat=allocate_stat)
      if (.not. allocated(options0%rec_energy)) then
        write(OUTDOCK, *) 
     &    "ERROR ---> Failed to dynamically allocate rec_energy"
        stop
      endif
      allocate(options0%solvnames(options0%total_receptors), 
     &    stat=allocate_stat)
      if (.not. allocated(options0%solvnames)) then
        write(OUTDOCK, *) 
     &    "ERROR ---> Failed to dynamically allocate solvnames"
        stop
      endif
      allocate(options0%hsolvnames(options0%total_receptors), 
     &    stat=allocate_stat)
      options0%hsolvnames = " " ! this is so that gnu works when no hsolv map is specified. (Trent E Balius)
      if (.not. allocated(options0%hsolvnames)) then
        write(OUTDOCK, *) 
     &    "ERROR ---> Failed to dynamically allocate hsolvnames"
        stop
      endif
      allocate(options0%phinames(options0%total_receptors), 
     &    stat=allocate_stat)
      if (.not. allocated(options0%phinames)) then
        write(OUTDOCK, *) 
     &    "ERROR ---> Failed to dynamically allocate phinames"
        stop
      endif
      !endif
c     teb
      allocate(options0%gistnames(options0%total_receptors),
     &    stat=allocate_stat)
      if (.not. allocated(options0%gistnames)) then
        write(OUTDOCK, *)
     &    "ERROR ---> Failed to dynamically allocate gistnames"
        stop
      endif
      allocate(options0%gistflag(options0%total_receptors),
     &    stat=allocate_stat)
      allocate(options0%gist_aprox(options0%total_receptors),
     &    stat=allocate_stat)
      ! intialize gistflag to be all false
      do idx = 1, options0%total_receptors
         options0%gistflag(idx) = .false.
         options0%gist_aprox(idx) = 0
      enddo 
c     teb
      allocate(options0%recdesrnames(options0%total_receptors),
     &    stat=allocate_stat)
      if (.not. allocated(options0%recdesrnames)) then
        write(OUTDOCK, *)
     &    "ERROR --> Failed to dynamically allocate recdesrnames"
        stop
      endif
      allocate(options0%hrecdesnames(options0%total_receptors),
     &    stat=allocate_stat)
      if (.not. allocated(options0%hrecdesnames)) then
        write(OUTDOCK, *)
     &    "ERROR --> Failed to dynamically allocate hrecdesnames"
        stop
      endif
      allocate(options0%rec_d_flag(options0%total_receptors),
     &    stat=allocate_stat)
      allocate(options0%hrec_d_flag(options0%total_receptors),
     &    stat=allocate_stat)
      do idx = 1, options0%total_receptors
        options0%rec_d_flag(idx) = .false.
        options0%hrec_d_flag(idx) = .false.
      enddo
      allocate(options0%recdesnames(options0%total_receptors), 
     &    stat=allocate_stat)
      if (.not. allocated(options0%recdesnames)) then
        write(OUTDOCK, *) 
     &    "ERROR ---> Failed to dynamically allocate recdesnames"
        stop
      endif
      allocate(options0%vdwnames(options0%total_receptors), 
     &    stat=allocate_stat)
      if (.not. allocated(options0%vdwnames)) then
        write(OUTDOCK, *) 
     &    "ERROR ---> Failed to dynamically allocate vdwnames"
        stop
      endif
      allocate(options0%bmpnames(options0%total_receptors), 
     &    stat=allocate_stat)
      if (.not. allocated(options0%bmpnames)) then
        write(OUTDOCK, *) 
     &    "ERROR ---> Failed to dynamically allocate bmpnames"
        stop
      endif

      return
      end subroutine allocate_flexible

      end module optionstype
