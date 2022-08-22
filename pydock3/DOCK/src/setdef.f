c set defaults
c every parameter falls into one of two mutually exclusive groups:
c 1. optional with default
c 2. required without default
c while the INDOCK file can contain either

      subroutine setdef(phimap0, options0, spheres0, recdes0, ligdes0)

      use phimaptype
      use optionstype
      use spheres

      implicit none

      type(phimap), intent(inout) :: phimap0
      type(options), intent(inout) :: options0
      type(spherest), intent(inout) :: spheres0
      type(phimap), intent(inout) :: recdes0
      type(phimap), intent(inout) :: ligdes0

c     letter: used in test of ascii to integer conversion rules.
      character (len=1) letter



c search type selection:
      options0%search_flag = .true.
      options0%rescore_flag = .false.
      options0%mol2inputfile = ''
      options0%ligsolinputfile = ''
      options0%ligvdwinputfile = ''

c more defaults
      
      options0%use_gist = .false.
      options0%use_gist_val = ''
      options0%per_atom_scores = .true.
      options0%receptor_desolv = .false.
      options0%rec_des_r  = .false.
      options0%hrec_des_r  = .false.
      !options0%rec_d_flag(1) = .false.
      !options0%hrec_d_flag(1) = .false.
      !options0%gistflag(1)   = .false.
      !options0%gist_aprox = .false.
      !options0%gist_aprox = 0
      options0%recdes_valence_penalty = .true.
      options0%hsolflag = .false.
      options0%flexible_receptor = .false.     
      !options0%current_number = 0 !how many flexible receptors have been read in
      options0%current_number = 1 !how many flexible receptors have been read in
      options0%dockovalent = .false. !covalent docking off by default
      options0%dihstep     = 20.0
c filenames:
      !options0%gistnames(options0%current_number) = 'undefined'
      !options0%gistnames(options0%current_number)

c defaults added for version 3.7
      options0%save_limit = 9999999.9 !very high score, save all poses
      phimap0%nsize = 193
      recdes0%nsize = 193
      ligdes0%nsize = 193
      options0%escale = 1.0
      options0%gistscale = 1.0
      options0%gistvolume = 0.0
      options0%gistH = 0
      options0%vscale = 1.0
      options0%solvscale = 1.0
      options0%rescale = 1.0 !receptor energy scale
      options0%rdscale = 1.0 !receptor desolvation scale
      options0%rhscale = 1.0 !receptor hydrophobic effect scale
      options0%iscale = 0.0 !internal energy disabled here for now
      options0%nsav = 1
      options0%nwrite = 1
      options0%input_flush_int = 1000
      options0%bmpvdw = 10.0
      options0%bmpvdwrigid = 10.0
      options0%bmptotal = -10.0
      options0%score_each_flex = .false. !save a score for each flexible receptor combination if true
c matching stuff
      options0%match_method = 1 !default matching method 
      options0%dislim = 0.05 !default distance limit (starting if match_method 2)
      options0%disstep = 0.05 !default stepsize, not used unless match_method set to 2
      options0%dismax = 0.5 !default maximum distance tolerance, not used unless mm = 2
      options0%matchgoal = 1000 !how many orientations we want to find, at a mininum
                       !under match_method = 2. if this hasn't been found,
                       !increment dislim by disstep and re-run matching.
                       !if dismax has been reached, also stop
      options0%timeout = 60.0 !in seconds
c chemical matching (coloring) initialization
      options0%cmatch = .false.
      spheres0%nlrmat = 0
      spheres0%nligcl = 0
      spheres0%nreccl = 0
      spheres0%numcol = 0
      spheres0%colrej = 0
      spheres0%coltest= 0
      options0%casen = .true.
      options0%check_clash = .true.
      options0%do_premax = .false.
      options0%astar = .false.

c minimization defaults and degeneracy checking
      options0%allocated_seeds = .false.
      options0%minimize = 0
      options0%mol2_minimizer = 0
      options0%output_premin = 0
      options0%num_mini_1 = 1
      options0%sim_itmax = 500
      options0%sim_trnstep = 0.2
      options0%sim_rotstep = 5.0
      options0%sim_need_to_restart = 1.0
      options0%sim_cnvrge = 0.1
      options0%min_cut = 1.0e15
      options0%nvar = 6 ! 3 Translation + 3 Rotation
      options0%alpha = 1.0
      options0%beta = 0.5
      options0%gamma = 2.0
      options0%big_energy = 9999.999
      options0%iseed = 777

c Monte Carlo minimization defaults
      options0%mc_minimize = 0
      options0%mc_itmax = 500
      options0%mc_accpt = 250
      options0%mc_trnstep = 0.2
      options0%mc_rotstep = 5.0
      options0%mc_temp = 298.15
      options0%mc_iseed = 777

c clustering docking
c docking mode
      !options0%clustering_docking = .false.
      options0%docking_mode = -1

      return
      end
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
c

