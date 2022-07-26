c this is called once per line of INDOCK with keywrd set to the first token
c and args set to the remaining non-whitespace stuff.
      ! XXX: OUTDOCK ambiguous in GFortran (TS)
      !subroutine dokw(keywrd, args, OUTDOCK, phimap0, recdes0, 
      subroutine dokw(keywrd, args, OUTDOCK_h, phimap0, recdes0, 
     &    ligdes0, options0, spheres0, errflg)

      use errorformats
      use phimaptype
      use optionstype
      use spheres
      use filenums

      !XXX: Trying to get GFortran to work
      !implicit none

      character (len=30) keywrd
      character (len=132) args
c outdock file specifier, so data can be written out if there are problems
      integer, intent(in) :: OUTDOCK_h
c these are user-defined types that contain information about grids (usually)
      type(phimap), intent(inout) :: phimap0, recdes0, ligdes0
      type(options), intent(inout) :: options0
      type(spherest), intent(inout) :: spheres0
      logical, intent(inout) :: errflg

      integer pos1, pos2, pos3, lenwrd
      integer nf, jj
      logical done
      character (len=30) :: clnam1
      character (len=30) :: clnam2

c     functions:
      integer nflds, nxtfld, nxtsp
      integer search_type ! 1= standard; 2 = rescore
      integer val ! temp value
      integer inputstatus

c     argument type conversion (integer --> logical)
c     argument validity checking (some here, some in ckparm)
c     (check filenames for alphanumeric?)

      done = .true.
      nf = nflds(args)

      if (keywrd .eq. 'distance_tolerance') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%dislim
      else if (keywrd .eq. 'distance_step') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%disstep
      else if (keywrd .eq. 'distance_maximum') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%dismax
      else if (keywrd .eq. 'nodes_maximum') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%maxnodes
      else if (keywrd .eq. 'nodes_minimum') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%minnodes
      else if (keywrd .eq. 'match_method') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%match_method
      else if (keywrd .eq. 'match_goal') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%matchgoal
      else if (keywrd .eq. 'bump_maximum') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%bmpvdw
      else if (keywrd .eq. 'bump_rigid') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%bmpvdwrigid
      else if (keywrd .eq. 'total_strain') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%total_strain
      else if (keywrd .eq. 'max_strain') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%max_strain
      else if (keywrd .eq. 'mol2_score_maximum') then ! total energy bump. 
        if (nf .lt. 1) goto 940
        read(args, *) options0%bmptotal
      else if (keywrd .eq. 'electrostatic_scale') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%escale
      else if (keywrd .eq. 'rec_des_scale') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%rec_des_scale
      else if (keywrd .eq. 'gist_scale') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%gistscale
      else if (keywrd .eq. 'gist_volume') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%gistvolume
      else if (keywrd .eq. 'k_clusters') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%k_clusters
      else if (keywrd .eq. 'gist_hydrogens') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%gistH = 0
        else if (args(1:1).eq.'0') then
          options0%gistH = 0
          write(6,*) 'gist will use all hydrogens with vdw radii.' 
          write(6,*) 'this is recommended for resocoring, not docking.' 
          write(6,*) 'use gist_hydrogens= 1 or 2 for docking' 
        else if (args(1:1).eq.'1') then
          options0%gistH = 1
          write(6,*) 'gist scoring componant will not use hydrogens.' 
        else if (args(1:1).eq.'2') then
          options0%gistH = 2
          write(6,*) 'gist scoring componant will use radius of 0.8'
          write(6,*) ' for nonpolar hydrogens'
        else
          goto 940 
        endif
      else if (keywrd .eq. 'vdw_scale') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%vscale
      else if (keywrd .eq. 'internal_scale') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%iscale
      else if (keywrd .eq. 'ligand_desolv_scale') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%solvscale
      else if (keywrd .eq. 'save_limit') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%save_limit
      else if (keywrd .eq. 'delphi_nsize') then
        if (nf .lt. 1) goto 940
        read (args, *) phimap0%nsize
      else if (keywrd .eq. 'recdes_nsize') then
        if (nf .lt. 1) goto 940
        read (args, *) recdes0%nsize
      else if (keywrd .eq. 'ligdes_nsize') then
        if (nf .lt. 1) goto 940
        read (args, *) ligdes0%nsize
      else if (keywrd .eq. 'atom_minimum') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%natmin
      else if (keywrd .eq. 'atom_maximum') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%natmax
      else if (keywrd .eq. 'number_save') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%nsav
      else if (keywrd .eq. 'number_write') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%nwrite
      else if (keywrd .eq. 'flush_int') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%input_flush_int
      else if (keywrd .eq. 'timeout') then
        if (nf .lt. 1) goto 940
        read(args,*) options0%timeout
      else if (keywrd .eq. 'chemical_matching') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%cmatch = .true.
        else if (args(1:1) .eq. 'y') then
          options0%cmatch = .true.
        else if (args(1:1) .eq. 'n') then
          options0%cmatch = .false.
        endif
      else if (keywrd .eq. 'score_each_flex') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%score_each_flex = .false.
        else if (args(1:1) .eq. 'y') then
          options0%score_each_flex = .true.
        else if (args(1:1) .eq. 'n') then
          options0%score_each_flex = .false.
        endif
      else if (keywrd .eq. 'case_sensitive') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%casen = .true.
        else if (args(1:1) .eq. 'y') then
          options0%casen = .true.
        else if (args(1:1) .eq. 'n') then
          options0%casen = .false.
        endif
      else if (keywrd .eq. 'match') then
        if (nf .lt. 2) goto 940
        clnam1 = ' '
        pos1=nxtfld(args,1)
        pos2=nxtsp(args,pos1)
        pos3=nxtfld(args,pos2)
        if ((pos3-pos1).ge.30) then
          clnam1=args(pos1:pos1+29)
        else
          lenwrd = pos2 - pos1
          clnam1(1:lenwrd)=args(pos1:pos2-1)
        endif
        call tolow(clnam1)
        clnam2 = ' '
        pos1=nxtfld(args,pos2)
        pos2=nxtsp(args,pos1)
        pos3=nxtfld(args,pos2)
        if ((pos3-pos1).ge.30) then
          clnam2=args(pos1:pos1+29)
        else
          lenwrd = pos2 - pos1
          clnam2(1:lenwrd)=args(pos1:pos2-1)
        endif
        call tolow(clnam2)
        call addmch(clnam1, clnam2, spheres0)
      else
        done = .false.
      endif
      if (done) then
        return !return once we've stopped reading match lines
      endif
      
      if (keywrd .eq. 'ligand_desolvation') then
        if (nf .lt. 1) then
          options0%ldesolv = 0
        else
          call tolow(args)
c         full desolvation
          if (args(1:1) .eq. 'f' ) then
            options0%ldesolv = 1
c         volume-based context-dependent desolvation - MMM 
          else if (args(1:1) .eq. 'v') then
            options0%ldesolv = 2
          else if (args(1:1) .eq. 'q') then !PB/qnifft style ligand desolvation
            options0%ldesolv = 3
c         no desolvation
          else 
            options0%ldesolv = 0
          endif
        endif
      elseif (keywrd .eq. 'search_type') then
        if (nf .lt. 1) goto 940
        read(args, *) search_type
        if (search_type .eq. 1) then
           options0%search_flag = .true.
           options0%rescore_flag = .false.
        else if (search_type .eq. 2) then
           options0%rescore_flag =  .true.
           options0%search_flag = .false.
        else
          write(6,*) 'Warning: not a valid search_type options.' 
          write(6,*) 'standard search (1) will be performed.'
           !Stop
        endif
      ! this is compatable with rescoring only 
      elseif (keywrd .eq. 'mol2file') then
        if (nf .lt. 1) goto 940
        read(args, '(a255)') options0%mol2inputfile 
        write(6,*) 'mol2 input is compatable with search_type = 2'
        write(6,*) 'filename = ', options0%mol2inputfile 
      elseif (keywrd .eq. 'ligsolfile') then
        if (nf .lt. 1) goto 940
        read(args, '(a255)') options0%ligsolinputfile
        write(6,*) 'mol2 input is compatable with search_type = 2'
        write(6,*) 'filename = ', options0%ligsolinputfile
      elseif (keywrd .eq. 'ligvdwfile') then
        if (nf .lt. 1) goto 940
        read(args, '(a255)') options0%ligvdwinputfile
        write(6,*) 'mol2 input is compatable with search_type = 2'
        write(6,*) 'filename = ', options0%ligvdwinputfile
c    & that is rescoring", options0%mol2inputfile
      else if (keywrd .eq. 'receptor_sphere_file') then
        if (nf .lt. 1) goto 940
        read(args,'(a80)') options0%clufil
      else if (keywrd .eq. 'ligand_atom_file') then
        if (nf .lt. 1) goto 940
        read(args,'(a255)') options0%ligfil
        if (options0%ligfil .eq. 'split_database_index') then
          options0%sdi = .true.
          open(unit=SDIFILE, file='split_database_index',
     &         err=980, status='old', action='read')
          read(SDIFILE,'(a255)') options0%ligfil
c     jklyu, 20200223
        else if (options0%ligfil .eq. 'rig_frag_dock_index') then
          options0%sdi = .true.
          open(unit=SDIFILE, file='rig_frag_dock_index',
     &         err=980, status='old', action='read')
          call readsdiln(SDIFILE,OUTDOCK,options0,inputstatus,errflg)
        else
          options0%sdi = .false.
        endif
      else if (keywrd .eq. 'using_gist') then ! this is for sorting
        if (nf .lt. 1) goto 940
        read(args,'(a255)') options0%use_gist_val
        if (options0%use_gist_val .eq. 'yes') then
          options0%use_gist = .true.
        else if (options0%use_gist_val .eq. 'no') then
          options0%use_gist = .false.
        else
          write(6,*) ' using_gist should have a value of yes or no.' 
          write(6,*) ' using_gist is set to false'
          options0%use_gist = .false.
        endif

      else if (keywrd .eq. 'output_file_prefix') then
        if (nf .lt. 1) goto 940
        read(args,'(a255)') options0%outfil
      else if (keywrd .eq. 'vdw_parameter_file') then
        if (nf .lt. 1) goto 940
        read(args, '(a80)') options0%vdwfil
c new keywords
      else if (keywrd .eq. 'check_clashes') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%check_clash = .true.
        else if (args(1:1) .eq. 'n') then
          options0%check_clash = .false.
        else
          options0%check_clash = .true.
        endif
      else if (keywrd .eq. 'check_strain') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%check_strain = .true.
        else if (args(1:1) .eq. 'n') then
          options0%check_strain = .false.
        else
          options0%check_strain = .true.
        endif
      else if (keywrd .eq. 'do_premax') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%do_premax = .false. !off by default
        else if (args(1:1) .eq. 'n') then
          options0%do_premax = .false.
        else
          options0%do_premax = .true.
        endif
      else if (keywrd .eq. 'do_clusters') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%do_clusters = .false. !off by default
        else if (args(1:1) .eq. 'n') then
          options0%do_clusters = .false.
        else
          options0%do_clusters = .true.
        endif
      else if (keywrd .eq. 'astar') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%astar = .false. !off by default
        else if (args(1:1) .eq. 'y') then
          options0%astar = .true.
        else
          options0%astar = .false.
        endif
      else if (keywrd .eq. 'receptor_desolv') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%receptor_desolv = .false. !off by default
        else if (args(1:1) .eq. 'y') then
          options0%receptor_desolv = .true.
        else
          options0%receptor_desolv = .false.
        endif
      else if (keywrd .eq. 'recdes_valence_penalty') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%recdes_valence_penalty = .true. !ON by default
        else if (args(1:1) .eq. 'y') then
          options0%recdes_valence_penalty = .true.
        else
          options0%recdes_valence_penalty = .false.
        endif
      else if (keywrd .eq. 'per_atom_scores') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%per_atom_scores = .true. !off by default
        else if (args(1:1) .eq. 'y') then
          options0%per_atom_scores = .true.
        else
          options0%per_atom_scores = .false.
        endif
!     start of docking_mode options
!     0: ensemble_docking with pre-aligned db2 files
!     1: rigid fragment sampling
!     2: read in orientation matrix then scoring
!     3: read db2 files from certain position then dock
!     4: subcluster docking
!     numbers that is not 0-4: normal docking
      else if (keywrd .eq. 'docking_mode') then
         !call tolow(args)
         if (nf .lt. 1) goto 940
         read(args, *) options0%docking_mode
!     start of covalent docking options
      else if (keywrd .eq. 'dockovalent') then
         call tolow(args)
         if (nf .lt. 1) then
            options0%dockovalent = .false. !off by default
         else if (args(1:1) .eq. 'y') then
            options0%dockovalent = .true.
         else
            options0%dockovalent = .false.
        endif
      else if (keywrd .eq. 'bond_len') then
        if (nf .lt. 1) goto 940
        read (args, *) options0%bondlen
      else if (keywrd .eq. 'bond_ang1') then
        if (nf .lt. 1) goto 940
        read (args, *) options0%bondang1
      else if (keywrd .eq. 'bond_ang2') then
        if (nf .lt. 1) goto 940
        read (args, *) options0%bondang2
      else if (keywrd .eq. 'len_range') then
        if (nf .lt. 1) goto 940
        read (args, *) options0%bondrange
      else if (keywrd .eq. 'dihedral_step') then
        if (nf .lt. 1) goto 940
        read (args, *) options0%dihstep
      else if (keywrd .eq. 'len_step') then
        if (nf .lt. 1) goto 940
        read (args, *) options0%bondstep
      else if (keywrd .eq. 'ang1_range') then
        if (nf .lt. 1) goto 940
        read (args, *) options0%bondang1range
      else if (keywrd .eq. 'ang1_step') then
        if (nf .lt. 1) goto 940
        read (args, *) options0%bondang1step
      else if (keywrd .eq. 'ang2_range') then
        if (nf .lt. 1) goto 940
        read (args, *) options0%bondang2range
      else if (keywrd .eq. 'ang2_step') then
        if (nf .lt. 1) goto 940
        read (args, *) options0%bondang2step
!     start of flexible receptor options, also all grids processed here now
      else if (keywrd .eq. 'flexible_receptor') then
        call tolow(args)
        if (nf .lt. 1) then
          options0%flexible_receptor = .false. !defaults to false
          options0%total_receptors = 1
          !call allocate_flexible(options0, OUTDOCK) ! XXX: Ambiguous
          call allocate_flexible(options0, OUTDOCK_h)
        else if (args(1:1) .eq. 'y') then
          options0%flexible_receptor = .true. !the number of receptors is coming
              !in the next total_receptors keyword
        else if (args(1:1) .eq. 'n') then 
          options0%flexible_receptor = .false.
          options0%total_receptors = 1
          !options0%current_number = 1
          !call allocate_flexible(options0, OUTDOCK) ! XXX: Ambiguous
          call allocate_flexible(options0, OUTDOCK_h)
          options0%rec_energy(1) = 0.0 
        endif
      else if (keywrd .eq. 'total_receptors') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%total_receptors
        !important to now allocate the various flexible receptor parts
        !call allocate_flexible(options0, OUTDOCK) ! XXX: Ambiguous
        call allocate_flexible(options0, OUTDOCK_h)
      else if (keywrd .eq. 'rec_number') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%current_number
        options0%rec_numbers(options0%current_number) = 
     &      options0%current_number
        options0%gistnames(options0%current_number) = "NONE" ! in case gist grid is not needed
        options0%rec_d_flag(options0%current_number) = .false.
        options0%hrec_d_flag(options0%current_number) = .false.
        options0%recdesrnames(options0%current_number) = "NONE" 
        options0%hrecdesnames(options0%current_number) = "NONE"
      else if (keywrd .eq. 'rec_group') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%group_numbers(options0%current_number)
      else if (keywrd .eq. 'rec_group_option') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%part_numbers(options0%current_number)
      else if (keywrd .eq. 'rec_energy') then
        if (nf .lt. 1) goto 940
        read(args, *) options0%rec_energy(options0%current_number)
      else if (keywrd .eq. 'solvmap_file') then
        if (nf .lt. 1) goto 940
        read(args,'(a80)') options0%solvnames(options0%current_number)
      else if (keywrd .eq. 'hydrogen_solvmap_file') then
        if (nf .lt. 1) goto 940
        read(args,'(a80)') options0%hsolvnames(options0%current_number)
      else if (keywrd .eq. 'delphi_file') then
        if (nf .lt. 1) goto 940
        read(args,'(a80)') options0%phinames(options0%current_number)
      else if (keywrd .eq. 'rec_des_r_file') then
        if (nf .lt. 1) goto 940
        options0%rec_d_flag(options0%current_number) = .true.
        read(args,'(a80)')options0%recdesrnames(options0%current_number)
c        options0%rec_des_r  = .true.
        if (options0%recdesrnames(options0%current_number) == "0") then
           options0%rec_d_flag(options0%current_number) = .false.
c            options0%rec_des_r = .false.
        endif
      else if (keywrd .eq. 'hrec_des_r_file') then
        if (nf. lt. 1) goto 940
        options0%hrec_d_flag(options0%current_number) = .true.
        read(args,'(a80)')options0%hrecdesnames(options0%current_number)
c          options0%hrec_des_r  = .true.
        if (options0%hrecdesnames(options0%current_number) == "0") then
            options0%hrec_d_flag(options0%current_number) = .false.
c             options0%hrec_des_r = .false.
        endif
      else if (keywrd .eq. 'gist_file') then
        if (nf .lt. 1) goto 940
        options0%gistflag(options0%current_number) = .true.
        read(args,'(a80)') options0%gistnames(options0%current_number)
        if (options0%gistnames(options0%current_number) == "0") then ! if we want gist to always have a zerro value, this is important for the flexible receptor code. 
           options0%gistflag(options0%current_number) = .false.
        endif
      else if (keywrd .eq. 'gist_aprox') then
        !if (args(1:1) .eq. 'y') then
        !  options0%gist_aprox = .true. !the number of receptors is coming
        !else if (args(1:1) .eq. 'n') then
        !  options0%gist_aprox = .false.
        !endif
        !read(args,'(i1)') val
        read(args,*) val
        options0%gist_aprox(options0%current_number) = val
        !if (val .ge. 0 .and. val .le. 2) then
        if (val .eq. 0 ) then
           write(OUTDOCK_h,*) "no gist aproxamations is used."
        else if (val .eq. 1 ) then   
           write(OUTDOCK_h,*) "gist aproxamation '1' is used
     & (use closest point, for precomputed)."
        else if (val .eq. 2 ) then   
           write(OUTDOCK_h,*) "gist aproxamation '2' is used 
     & (use 8 closest point)."
        else
           write(OUTDOCK_h,*) 'gist_aprox must be 0 (no aprox),
     & 1 (use closest point, for precomputed), 
     & or 2 (use closest 8 points).'
           stop 
        endif
      else if (keywrd .eq. 'recdes_file') then
        if (nf .lt. 1) goto 940
        read(args,'(a80)') options0%recdesnames(options0%current_number)
      else if (keywrd .eq. 'chemgrid_file') then
        if (nf .lt. 1) goto 940
        read(args,'(a80)') options0%vdwnames(options0%current_number)
      else if (keywrd .eq. 'bumpmap_file') then
        if (nf .lt. 1) goto 940
        read(args,'(a80)') options0%bmpnames(options0%current_number)
c     MONTE CARLO MINIMIZATION
      else if (keywrd .eq. 'monte_carlo') then
         call tolow(args)
         if (nf. lt. 1) then
            write(OUTDOCK_h,*) "Monte Carlo is off"
            options0%mc_minimize = 0
         else if (args(1:1) .eq. 'n') then
            options0%mc_minimize = 0
         else
            if (options0%minimize .eq. 1) then
               write(OUTDOCK_h,*) "Monte Carlo cannot 
     &  be used with Simplex"
               stop
            else if (options0%minimize .eq. 0) then
               write(OUTDOCK_h,*) "Monte Carlo is on"
               options0%mc_minimize = 1
            endif
         endif
      else if (keywrd .eq. 'mc_itmax') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%mc_itmax
      else if (keywrd .eq. 'mc_accpt') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%mc_accpt
      else if (keywrd .eq. 'mc_temp') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%mc_temp
      else if (keywrd .eq. 'mc_trnstep') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%mc_trnstep
      else if (keywrd .eq. 'mc_rotstep') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%mc_rotstep
      else if (keywrd .eq. 'mc_iseed') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%mc_iseed
c     SIMPLEX MINIMIZATION
      else if (keywrd .eq. 'minimize') then
        call tolow(args)
        if (nf .lt. 1) then
          write(OUTDOCK_h,*) "minimization is off"
          options0%minimize = 0
        else if (args(1:1) .eq. 'n') then
          options0%minimize = 0
        else
          if (options0%mc_minimize .eq. 1) then
              write(OUTDOCK_h,*) "Simplex cannot be used with
     & Monte Carlo minimization"
              stop
          else 
              write(OUTDOCK_h,*) "minimization is on"
              options0%minimize = 1
          endif
        endif
      else if (keywrd .eq. 'mol2_minimizer') then
        call tolow(args)
        if (nf .lt. 1) then
          write(OUTDOCK_h,*) "mol2_minimizer is off"
          options0%mol2_minimizer = 0
        else if (args(1:1) .eq. 'n') then
          options0%mol2_minimizer = 0
        else
          write(OUTDOCK_h,*) "mol2_minimizer is on"
          options0%mol2_minimizer = 1
        endif
      else if (keywrd .eq. 'output_premin') then
        call tolow(args)
        if (nf .lt. 1) then
          write(OUTDOCK_h,*) "output_premin is off"
          options0%output_premin = 0
        else if (args(1:1) .eq. 'n') then
          options0%output_premin = 0
        else
          write(OUTDOCK_h,*) "output_premin is on"
          options0%output_premin = 1
        endif
      else if (keywrd .eq. 'num_mini_1') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%num_mini_1
      else if (keywrd .eq. 'sim_itmax') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%sim_itmax
      else if (keywrd .eq. 'sim_trnstep') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%sim_trnstep
      else if (keywrd .eq. 'sim_rotstep') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%sim_rotstep
      else if (keywrd .eq. 'sim_need_to_restart') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%sim_need_to_restart
      else if (keywrd .eq. 'sim_cnvrge') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%sim_cnvrge
      else if (keywrd .eq. 'min_cut') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%min_cut
      else if (keywrd .eq. 'big_energy') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%big_energy
      else if (keywrd .eq. 'iseed') then
        call tolow(args)
        if (nf .lt. 1) goto 940
        read (args, *) options0%iseed
c start of deprecated keywords
      else if (keywrd .eq. 'molecules_maximum') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'vdw_minimum') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'vdw_maximum') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'ligand_binsize') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'ligand_overlap') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'receptor_binsize') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'receptor_overlap') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'electrostatic_limit') then
        goto 880 !deprecated
      else if (keywrd .eq. 'chemgrid_file_prefix') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'ratio_minimum') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'normalize_save') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'restart_interval') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'mirror_ligands') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'scoring_option') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'distmap_file') then !bumping gone
        goto 880 ! deprecated !bumping gone
      else if (keywrd .eq. 'ligand_sphere_file') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'mixmatch') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'remove_positive_solvation') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'write_coordinates') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'random_seed') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'focus_type') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'focus_cycles') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'focus_bump') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'critical_clusters') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'cluster_numbers') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'standard_pdb') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'focus_ratio') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'mode') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'ligand_type') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'restart') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'output_hydrogens') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'recombine_fragments') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'interpolate') then
        goto 880 ! deprecated
      else if (keywrd .eq. 'initial_skip') then
        goto 880 !deprecated
      else
        goto 900
      endif

      return

880   continue
      ! XXX: Ambiguous to Gfortran (TS)
      !write(OUTDOCK, WARNING1) 'deprecated keyword ignored: ', keywrd
      write(OUTDOCK_h, WARNING1) 'deprecated keyword ignored: ', keywrd
      return
900   continue
      ! XXX: Ambiguous to Gfortran (TS)
      !write(OUTDOCK, HALT2) 'keyword not recognized: ', keywrd
      write(OUTDOCK_h, HALT2) 'keyword not recognized: ', keywrd
      errflg = .true.
      return
940   continue
      ! XXX: Ambiguous to Gfortran (TS)
      !write(OUTDOCK, HALT2) 'keyword without argument: ', keywrd
      write(OUTDOCK_h, HALT2) 'keyword without argument: ', keywrd
      errflg = .true.
      return
980   continue
      ! XXX: Ambiguous to Gfortran (TS)
      !write(OUTDOCK, HALT2) 'cannot open input file: ', options0%ligfil
      write(OUTDOCK_h, HALT2) 'cannot open input file: ',
     & options0%ligfil
      errflg = .true.
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
