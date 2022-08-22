c----------------------------------------------------------------------
c
c     check for missing parameters
c     check for incorrect parameters
c     check for parameter compatibility

      subroutine ckparm(OUTDOCK, phimap0, recdes0, ligdes0, gist0,
     &    options0, vdwmap0, rec_des_r0,
     &    MAXOR, MAXCTR, errflg)

      use errorformats
      use phimaptype
      use gisttype
      use optionstype
      use vdwmaptype
      use spheres

      implicit none

c outdock file specifier, so data can be written out if there are problems
      integer, intent(in) :: OUTDOCK
c these are user-defined types that contain information about grids (usually)
      type(phimap), intent(inout) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(inout) :: gist0, rec_des_r0
      type(options), intent(inout) :: options0
      type(vdwmap), intent(inout) :: vdwmap0
      integer, intent(in) :: MAXOR, MAXCTR 
      logical, intent(inout) :: errflg !is there a problem?

      character (len=20) solvtype
      integer which_grid !which flexible receptor grid in examination

      if (options0%ldesolv .eq. 1) then
         solvtype = 'full'
      else if (options0%ldesolv .eq. 2) then
         solvtype = 'context-dependent'
      else
         solvtype = 'no'
      endif

      write(OUTDOCK, '(3a)') 'Solvation type: ', trim(solvtype),
     &     ' ligand desolvation'

      if (options0%use_gist) then
          write(OUTDOCK,'(a)') 
     &   'using_gist == yes. We need to Sort the confs in db2, if full',
     &   ' gist score is being used.'
      else
          write(OUTDOCK,'(a)')  
     &   'using_gist == no. Do not sort the confs in db2, if gist is',
     &    ' not used'
      endif

      write(OUTDOCK,'(a)') 
     &     'Ligands will NOT be recombined (NO mix-n-match)'
      if (options0%check_clash .eqv. .true.) then
        write(OUTDOCK,'(a)') 'Internal Distance Clashes checked'
      else
        write(OUTDOCK,'(a)') 'Internal Distance Clashes NOT checked'
      endif
      if (options0%do_premax .eqv. .true.) then
        write(OUTDOCK,'(a)') 
     &     'Experimental premax scoring will be used'
      else
        write(OUTDOCK,'(a)') 
     &     'Experimental premax scoring will NOT be used'
      endif
      if (options0%do_clusters .eqv. .true.) then
        write(OUTDOCK,'(a)') 
     &     'Clusters/Additional Match Spheres scoring will be used'
      else
        write(OUTDOCK,'(a)') 
     &     'Clusters/Additional Match Spheres scoring will not be used'
      endif
      if (options0%astar .eqv. .true.) then
        write(OUTDOCK,'(a)') 
     &     'Experimental best rigid first a* search'
      else
        write(OUTDOCK,'(a)') 'No a* search, normal search used'
      endif

      if (options0%dislim .lt. 0.0) then !distance tolerance cannot be less than 0
        options0%dislim = 0.05
      endif
      if (options0%maxnodes .lt. 3) then
        options0%maxnodes = 3 !maxnodes must be at least 3
      endif
      if (options0%minnodes .lt. 3) then
        options0%minnodes = 3 !minnodes must be at least 3
      endif
      if (options0%minnodes .gt. options0%maxnodes) then
        options0%minnodes = options0%maxnodes !minnodes cannot be greater than maxnodes
      endif
      write(OUTDOCK, 14) options0%dislim, options0%maxnodes,  
     &    options0%minnodes
   14 format(':dislim = ', f7.4, ' maxnodes = ', i4,
     &     ' minnodes = ', i4)
      if (options0%maxnodes .gt. maxctr) then
        write(OUTDOCK, HALT1) 'maxnodes higher than maximum allowed.'
        write(OUTDOCK, HALT3) '  Recompile with maxctr >',
     &      options0%maxnodes
        write(OUTDOCK, HALT0)
        stop
      endif
      write(OUTDOCK, *)
     &    'Matching method number ', options0%match_method, 
     &    ' will be used'
      if (options0%match_method .eq. 2) then
        if (options0%matchgoal .gt. MAXOR) then
          write(OUTDOCK, HALT1) 'matchgoal higher than MAXOR.'
          write(OUTDOCK, HALT3) '  Recompile with MAXOR >', 
     &        options0%matchgoal
          write(OUTDOCK, HALT0)
          stop
        endif
        write(OUTDOCK, 24) options0%disstep, options0%dismax, 
     &      options0%matchgoal
   24   format(':disstep = ', f7.4, ' dismax = ', f7.4,
     &     ' matchgoal = ', i8)
        write(OUTDOCK,*) "timeout = ", options0%timeout
      endif

      if (options0%clufil(1:1) .eq. char(0)) then
        errflg = .true.
        write(OUTDOCK, HALT1) 'Receptor spheres: missing file name'
      else
        write(OUTDOCK,'(a,a)') 'Receptor spheres: ', 
     &      trim(options0%clufil)
      endif

      if (options0%ligfil(1:1) .eq. char(0)) then
        errflg = .true.
        write(OUTDOCK, HALT1) 'Input ligand: missing file name'
      else
        write(OUTDOCK,'(a,a)') 'Input ligand: ', trim(options0%ligfil)
      endif

      if (options0%outfil(1:1) .eq. char(0)) then
        errflg = .true.
        write(OUTDOCK, HALT1) 'missing ligand output file name'
      endif

      write(OUTDOCK, *) 'DelPhi/Qnifft gridsize is ', phimap0%nsize
      if (options0%receptor_desolv) then
        write(OUTDOCK, *) 'Receptor desolvation being used'
        write(OUTDOCK, *) 'Receptor Desolvation/Qnifft gridsize is ', 
     &      recdes0%nsize
        if (options0%recdes_valence_penalty) then
          write(OUTDOCK, *) 'Receptor desolvation heavy atom ', 
     &        'valence penalty being used'
        else
          write(OUTDOCK, *) 'Receptor desolvation heavy atom ', 
     &        'valence penalty NOT being used'
        endif
      else
        write(OUTDOCK, *) 'Receptor desolvation not being used'
      endif

      if (options0%vdwfil(1:1) .eq. char(0)) then
        errflg = .true.
        write(OUTDOCK, HALT1) 'vdW parms:  missing file name'
      else
        write(OUTDOCK,'(a,a)') 'vdW parms: ', trim(options0%vdwfil)
      endif

c      write(OUTDOCK,*) 
c     &    'resolution = ', vdwmap0%grddiv, 
c     &    'number of points = ', vdwmap0%ngrd
      write(OUTDOCK, 19) options0%bmpvdw
      write(OUTDOCK, 20) options0%bmpvdwrigid
      write(OUTDOCK, 21) options0%bmptotal
   19 format(' a maximum of ',f6.1,
     &' in vdw energy score is allowed for any flexible group')
   20 format(' a maximum of ',f6.1,
     &' in vdw energy score is allowed for the whole rigid component')
   21 format(' a maximum of ',f6.1,
     &' in total energy score is allowed for any pose written to mol2.')
      if (options0%natmin .le. 0) then
        options0%natmin = 1 !atom minimum must be 1
      endif
      if (options0%natmax .le. 0) then
        options0%natmax = 80 !atom maximum can't be less than 0
      endif
      if (options0%natmax .gt. maxpts) then
        write(OUTDOCK, WARNING0) 'natmax cannot be greater than maxpts'
        write(OUTDOCK, WARNING2) 
     &     '   resetting natmax to equal maxpts: ', maxpts
        options0%natmax = maxpts
      endif
      write(OUTDOCK, 22) options0%natmin, options0%natmax
   22 format (' natmin = ',i5,'  natmax = ',i5)

      if (options0%nsav .lt. 0) then
        options0%nsav = 1 !must save at least one pose per molecule
      endif
      write(OUTDOCK, 27) options0%nsav
   27 format(' number of poses saved during processing   = ', i5)
      write(OUTDOCK, 28) options0%nwrite
   28 format(' number of poses written    = ', i5)
      write(OUTDOCK, 29) options0%input_flush_int
   29 format(' OUTDOCK will be updated every N numbers of lines, N = '
     & , i5)

      write(OUTDOCK, *)'hydrogens always written out'

      write(OUTDOCK, *) 'trilinear interpolation is always done'

      if (options0%vscale .ne. 1.0) then
        write(OUTDOCK, *)
     &      'the van der Waals E will be scaled by the factor ',
     &      options0%vscale
      endif
      if (options0%escale .ne. 1.0) then
        write(OUTDOCK, *)
     &      'the electrostatic E will be scaled by the factor ',
     &      options0%escale
      endif
      if (options0%rec_des_scale .ne. 1.0) then
        write(OUTDOCK, *)
     &      'the receptor desolvation will be scaled by ',
     &      options0%rec_des_scale
      endif
      if (options0%gistscale .ne. 1.0) then
        write(OUTDOCK, *)
     &      'the gist E will be scaled by the factor ',
     &      options0%gistscale
      endif
      if (options0%gistvolume .eq. 0.0) then ! use dx file information to caculate volume.
        write(OUTDOCK, *)
     &      'gistvolume = 0.0 (defult) --> ',
     &      'Gist-voxel volume will be calulated based on dxfile.'
      else if (options0%gistvolume .gt. 0.0) then ! should always be true.
        write(OUTDOCK, *)
     &      'user specified gist-voxel volume:',
     &      options0%gistvolume
      else 
        errflg = .true.
        write(OUTDOCK, HALT1) 'gist-voxel volume is not resonable.'
      endif

      if (options0%minimize .gt. 0) then
           !write(OUTDOCK, *) "performing minimization . . ." 
           write(OUTDOCK, '(a10,i10)') "min seed = ", options0%iseed 
      endif

c     we should fix this so that it works with min.
      if ((options0%minimize .gt. 0) .and. options0%gistflag(1)) then
          write(OUTDOCK, *) 
     &    'Error: gist does not work with minimization currently.'
          write(OUTDOCK, *)
     &    'either turn off minimization or use blurry gist.'
          errflg = .true.
      endif

c     minimization does not work with covalent docking
      if ((options0%minimize .gt. 0) .and. options0%dockovalent) then
          write(OUTDOCK, *)
     &    'Error: covalent does not work with minimization currently.'
          write(OUTDOCK, *)
     &    'Turn off minimization for covalent docking.'
          errflg = .true.
      endif

      if (options0%mc_minimize .gt. 0) then
         write(OUTDOCK, '(a10,i10)') "min_seed = ", options0%mc_iseed
      endif

c     gist also does not work with covalent docking
      if ((options0%dockovalent) .and. options0%gistflag(1)) then
          write(OUTDOCK, *)
     &    'Error: gist does not work with covalent dock currently.'
          write(OUTDOCK, *)
     &    'turn off gist for covalent docking.'
          errflg = .true.
      endif

      if (options0%solvscale .ne. 1.0) then
        write(OUTDOCK, *)
     &      'the ligand desolvation scaled by the factor ',
     &      options0%iscale
      endif
      if (options0%iscale .ne. 1.0) then
        write(OUTDOCK, *)
     &      'the ligand internal energy scaled by the factor ',
     &      options0%iscale
      endif

      do which_grid = 1, options0%total_receptors
        write(OUTDOCK, *) 'Receptor grids for Group#: ', which_grid
        if (options0%ldesolv .eq. 2) then
          if (options0%solvnames(which_grid)(1:1) .eq. char(0)) then
            errflg = .true.
            write(OUTDOCK, HALT1) 'Desolvation grid: missing file name'
          else
            write(OUTDOCK,'(a,a)') 'Desolvation grid: ', 
     &          trim(options0%solvnames(which_grid))
          endif
          options0%hsolflag = .false. !by default, only use one grid for solvation
          if ((options0%hsolvnames(which_grid)(1:1) .ne. char(0)) .and.
     &        (options0%hsolvnames(which_grid)(1:1) .ne. ' ')) then !sometimes use separate grids for heavy/hydrogen
            options0%hsolflag = .true.
            write(OUTDOCK,'(a,a)') 
     &          'Additional hydrogen desolvation grid: ',
     &          trim(options0%hsolvnames(which_grid))
          endif
        else if (options0%ldesolv .eq. 3) then !use qnifft grids
          write(OUTDOCK, *) 'Ligand Desolvation/Qnifft gridsize is ', 
     &        ligdes0%nsize
          if (options0%solvnames(which_grid)(1:1) .eq. char(0)) then
            write(OUTDOCK, WARNING0) 
     &          'must specify a DelPhi map for ligdes option 3'
            stop
          else
            write(OUTDOCK,'(a,a)') 'Ligand Desolvation grid: ', 
     &          trim(options0%solvnames(which_grid))
          endif
          options0%hsolflag = .false. !qnifft approximation doesn't need 2 grids
        endif
        if (options0%phinames(which_grid)(1:1) .eq. char(0)) then
          write(OUTDOCK, WARNING0) 
     &        'must specify a DelPhi map'
          stop
        else
          write(OUTDOCK,'(a,a)') 'DelPhi grid: ', 
     &        trim(options0%phinames(which_grid))
        endif
c        if (options0%recdesrnames(which_grid)(1:1) .ne. char(0)) then
        if (options0%rec_d_flag(which_grid)) then
c          options0%rec_d_flag(which_grid) = .true.
          write(OUTDOCK,'(a,a)') 'Heavy Receptor Desolvation grid: ',
     &       trim(options0%recdesrnames(which_grid))
        else
c          options0%rec_d_flag(which_grid) = .false. 
          write(OUTDOCK,'(a,a)') 'Rec Desolv is not being used'
        endif
        !if (options0%hrecdesnames(which_grid)(1:1) .ne. char(0)) then
        if (options0%hrec_d_flag(which_grid)) then
          write(OUTDOCK,'(a,a)') 'Hydrogen Receptor Desolvation grid: ',
     &       trim(options0%hrecdesnames(which_grid))
        endif
        if (options0%gistnames(which_grid)(1:1) .eq. '0' .or. 
     & options0%gistnames(which_grid)(1:1) .eq. ' ' .or.
     & options0%gistnames(which_grid) == 'NONE') then
          options0%gistflag(which_grid) = .false. 
          write(OUTDOCK, *) 'Gist grid not used. 0 or " " given.'
        else
          options0%gistflag(which_grid) = .true. 
          write(OUTDOCK,'(a,a)') 'Gist grid: ',
     &        trim(options0%gistnames(which_grid))
        endif
        if (options0%receptor_desolv) then
          write(OUTDOCK, '(a,a)') 'Receptor Desolvation grid: ',
     &        trim(options0%recdesnames(which_grid))
        endif
        if (options0%bmpnames(which_grid)(1:1) .eq. char(0)) then
          errflg = .true.
          write(OUTDOCK, HALT1) 'bumpmap filename: missing'
        else
          write(OUTDOCK,'(a,a)') 'bumpmap filename: ', 
     &        trim(options0%bmpnames(which_grid))
        endif
        if (options0%vdwnames(which_grid)(1:1) .eq. char(0)) then
          errflg = .true.
          write(OUTDOCK, HALT1) 'vdw filename: missing'
        else
          write(OUTDOCK,'(a,a)') 'vdw filename: ', 
     &        trim(options0%vdwnames(which_grid))
        endif
        write(OUTDOCK,'(a,f7.4)') 'receptor energy: ', 
     &      options0%rec_energy(which_grid)
      enddo

      if (options0%score_each_flex) then
        write(OUTDOCK, *) 'The best score for each flexible receptor ',
     &      'combination will be saved'
      endif

      write(OUTDOCK, *) 'Only poses with a score better than ',
     &    options0%save_limit, ' will be written to disk. ',
     &    'Use with caution for retrospective calculations, ',
     &    'Intended use is for prospective calculations. '


      return
      end
c----------------------------------------------------------------------
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
