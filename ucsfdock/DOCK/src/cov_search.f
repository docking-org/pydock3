c     This subroutine based on search.f utilizes the scoring machinary 
c     to score covalently attached poses.
c     The main changes are removal of matching code (covalent poses are 
c     explicitly generated) 
c
c     Modified search.f (RGC 2011) to support covalent docking (NL 2014)
c-----------------------------------------------------------------------
      module cov_search

      implicit none
      contains

      subroutine run_cov_search(grids0, phimap0, recdes0, ligdes0, 
     &    gist0, solvmap0,rec_des_r0, vdwmap0, options0, allminlimits, 
     &    allmaxlimits, spheres0, MAXTYV, sra, srb, MAXOR, MAXCTR, 
     &    fdlig)

      use phimaptype
      use gisttype
      use solvmaptype
      use vdwmaptype
      use status
      use optionstype
      use gridstype
      use db2type
      use ligscoretype
      use matchtype
      use filenums
      use spheres
      use score_mol
      use atomscoretype
      use score_conf
      use check_conf_og

c these are the dynamically allocated grids
      type(flexgrids), intent(inout) :: grids0
c these are user-defined types that contain information about grids (usually)
      type(phimap), intent(inout) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(inout) :: gist0,rec_des_r0
      type(solvmap), intent(inout) :: solvmap0
      type(vdwmap), intent(inout) :: vdwmap0
c options
      type(options), intent(in) :: options0
c for outside grid checking
      real, dimension(3) :: allminlimits, allmaxlimits
c spheres & chemical matching data
      type(spherest), intent(in) :: spheres0
c vdw parameters
      integer, intent(in) :: MAXTYV
      real, dimension(MAXTYV), intent(in) :: sra, srb !vdw parms
      integer, intent(in) :: MAXOR
      integer, intent(in) :: MAXCTR
      integer (kind=8), intent(inout) :: fdlig !ligand file handle
c     variables --
      type(db2) :: db2lig !ligand data stored here
      type(ligscoret) :: ligscore !saves many top poses
      type(ligscoret) :: ligscoreeach !saves one top pose for each flexible receptor combination
      type(atomscoret) :: atomscore !recalculate per atom scores
      type(matcht) :: match !stuff about matching
      character (len=80) molout !filename to save
      logical iend
c        iend:  indicator variable for end of ligand file.
c               is equal to zero if end not reached.
      integer molct
c        molct:  the number of ligands searched on this run.
      integer count !temp counter, don't use i,j,n or other one letter nonsense
c        best scoring conformation for a given ligand
      integer flush_int
c        # of compounds interval to flush buffer to OUTDOCK
      real tcorr(3, MAXPTS) !used to hold the rigid component heavy atoms
      integer tcolor(MAXPTS)
      integer tcount !count of rigid heavy atoms
      integer tcountorig !not counting cloud additional match spheres
      real telaps, time2, time1 !times
      logical exists
      integer dock_status !uses status.h enums
      integer istat !ligand file status
      integer tempsave, temprank !counter for multiple output
      integer cloudcount, cloudtotal !counter and overall # of times to run
      logical doing_cloud !whether or not the cloud is actually happening
      integer setstart, setend !controls which sets (complete confs) to score
      integer centrl !ligand seed (start) sphere
      integer centrr !receptor seed (start) sphere
c variables for match_method choosing, control
      real curdislim !temporary variable, keeps track of current distance limit
      logical more_matches !true if we want to find more matches.
c these are the new degeneracy controls. the 2 hashes are the sets of matching
c spheres, one for rec, one lig. the count is how many pairs have been added.
      integer hash(MAXOR, MAXCTR), hashcount
      integer (kind=8) :: fdmol !output mol2.gz file handle
      integer :: outcount
      logical temp_status

c format lines for output
      character (len=*), parameter :: molexam = 
     &     '(i8,"  molecules examined on this run")'
      character (len=*), parameter :: outdockerr = 
     &     '(i6,1x,a16,i8,1x,i10,1x,f7.2,1x,a)'
      character (len=*), parameter :: outdockline = '(i6,1x,a16,1x,a16,
     &       1x,i8,1x,i10,1x,f7.2,1x,i3,1x,i9,1x,i9,1x,i6,3x,i3,2x,
     &       f7.2,1x,f7.2,1x,f7.2,f7.2,1x,f7.2,1x,f7.2,1x,f7.2,1x,f7.2,
     &       1x,f7.2,f10.2)'

c     dockovalent local variables
      real p1(3),p2(3),p3(3) !receptor spheres - see cov_bond.f
      integer ligi1,ligi2,ligi3,ligi4 !cov atom and neighbour indices
      integer sp !indicates the geometry of cov atom (sp3/sp2)
      real lfi1mag              !magnitude of bins for fi1
      integer lensteps, ang1steps, ang2steps !sampling indices
      integer s,k,v !loop variables
      real a,fi,d          !cov bond parameters see cov_bond.f
      integer lfi      !loop increment variable
      real la,ld       !loop variables to iterate these
      real lb,lfi1         !loop variables to iterate the lig params
      real cov_loc(3)      !tmp vector to hold the covalent location
      integer atnum,resnum !DEBUG indices for pdb printing
      real rigid_part(3, MAXPTS) !tmp storage for rigid part
      real cl(3,MAXOR) ! tmp storage for all lig ensamble coords
      integer rigcount ! number of rigid part atoms including H's
      integer cli,clj,n,p
      integer tempcoord, tempatom 
      integer ligc, db2c, total_iter 
      logical skpr ! ligread2, skip_processing 

c     for scoring
      real asum, bsum, elescore, alscore, plscore, descore, hescore
      real vdwscore, bfact, gistval, rdscore
      !energy do not have a per-conformation score
      integer conf_status
c      integer, dimension(options0%total_receptors) :: 
c     &    conf_status !track status of each conf
      integer :: maingrid

c     Declare local constant Pi
      REAL, PARAMETER :: Pi = 3.1415927
      REAL, PARAMETER :: deg_to_rad = 0.0174532925

c     debug printing formats
 76   format(A,I5,A,A3,A,I6,A,F7.3,A,F7.3,A,F7.3,A,F6.2)
 77   format(A,I5,A,A3,A,I6,A,F7.3,A,F7.3,A,F7.3)
 78   format(A,I4)

      integer ::  gistindexlist(5000)
      integer :: lengistil
      integer :: dih_num, dih_step, dih_final! number of dihideral angles to be sampled.  
c     --- initialize variables ---
      
      skpr = .false.
      !initialize cloudcount to 1 - not implemented
      cloudcount = 1
      call allocate_ligscore(ligscore, options0%nsav, OUTDOCK)
      call allocate_atomscore(atomscore, maxpts, OUTDOCK)
      if (options0%flexible_receptor .and. 
     &    options0%score_each_flex) then
        call allocate_ligscore(ligscoreeach, grids0%total_combinations, 
     &      OUTDOCK)
      endif
      flush_int = 1000 !flush after every (this many) ligands
      molct = 0
      match%nmatch = 0 !initialize various things
      call allocate_match(match, maxor, OUTDOCK)
      iend = .false.

c     open output files
      molout = trim(options0%outfil)//'mol2.gz'
      write(OUTDOCK, '(a,a)') 'output file: ', molout
c     calculation of receptor sphere-sphere distances.
      call intdis(spheres0%nsphr, spheres0%spcorr, spheres0%disr, 
     &    match%dmaxr, MAXPTS)
      write(OUTDOCK, *) 'maximum receptor sphere-sphere distance',
     &    match%dmaxr

c     set the bin magnitude for fi1 angle and other dihideral sampling 
      lfi1mag=options0%dihstep
      dih_step = INT(CEILING(lfi1mag))
      dih_num = INT(CEILING(360.0 / lfi1mag - 1.0))
      dih_final = dih_num * dih_step
      write(OUTDOCK,*) 'dihederal step size: ', 
     & dih_step,  'number of diherals to sample: ', dih_num

      call gzopen(fdmol, 'w', molout, istat) !open output ligand file
      !top of main section of OUTDOCK file
c     write(OUTDOCK, '(a,a,a,a)') 
c    & '  mol#           id_num     flexiblecode  matched    nscored  ',
c    & '  time hac    setnum    matnum   rank cloud    elect +  gist',
c    & ' +   vdW + psol +  asol + inter + rec_e + rec_d + r_hyd =',
c    & '    Total'
        write(OUTDOCK, '(a,a,a,a)') 
     & '  mol#           id_num     flexiblecode  matched    nscored  ',
     & '  time hac    setnum    matnum   rank charge    elect +  gist',
     & ' +   vdW + psol +  asol + tStrain + mStrain + rec_d + r_hyd =',
     & '    Total'
      call doflush(OUTDOCK)

c     set the bin magnitude for fi1 angle
c     lfi1mag=options0%dihstep
c     dih_step = INT(CEILING(lfi1mag))
c     dih_num = INT(CEILING(360.0/ CEILING(lfi1mag)) - 1.0)
c     dih_final = dih_num * dih_step
c     write(OUTDOCK,*) 'dihederal step size: ', 
c    & dih_num,  'number of diherals to sample: ', dih_num

c     set bond length and angels iterators
      lensteps=int((options0%bondrange/options0%bondstep)*2+1)
      ang1steps=int((options0%bondang1range/options0%bondang1step)*2+1)
      ang2steps=int((options0%bondang2range/options0%bondang2step)*2+1)
      call doflush(OUTDOCK)
c     get receptor atoms from matching spheres file
      p1(1)=spheres0%spcorr(1,1) 
      p1(2)=spheres0%spcorr(2,1) 
      p1(3)=spheres0%spcorr(3,1)
      p2(1)=spheres0%spcorr(1,2) 
      p2(2)=spheres0%spcorr(2,2) 
      p2(3)=spheres0%spcorr(3,2)
      p3(1)=spheres0%spcorr(1,3) 
      p3(2)=spheres0%spcorr(2,3) 
      p3(3)=spheres0%spcorr(3,3)

      do while (.not. iend) !an actual loop
c       --- read a ligand structure from the dbase ---
        ligc = 0
        db2c = 0
        call ligread2(db2lig, options0, 
     &      iend, tcorr, tcolor, tcount, fdlig,  ligc, db2c, skpr) !read one hierarchy 

        if (iend) then !if end of file read, no more ligands
          exit !break out of do loop
        endif
        dock_status = NOMATCH
        ligscore%savedcount = 0 !reset the number of poses saved
        ligscore%numscored = 0
        if (options0%flexible_receptor .and. 
     &       options0%score_each_flex) then
           call reset_ligscore(ligscoreeach)
        endif
c       time elapsed init
        telaps = 0.0
        molct = molct + 1

c       process all conformations
        setstart = 1
        setend = db2lig%total_sets !just do them all

        !detect covalent attachment point and its neighbours
        !assumptions:
        !- cov attach point color is 8
        !- cov n1 to n3 are 9,10,11
        !- if 9,10,11 present - assume sp3 geometry
        !- if only 9,10 present - assume sp2 geometry
        ligi1=0; ligi2=0; ligi3=0; ligi4=0;

        !number of atoms in conf 1 - the rigid part
        rigcount = db2lig%conf_coord(2, 1) -
     &                   db2lig%conf_coord(1, 1) + 1

        do count=1,rigcount !rigcount is the number of atoms in rigid component
           if ( db2lig%atom_color(db2lig%coord_index(2,count)) .eq. 8) !assumes covalent atom color is 8
     &        then            
              ligi1=count
           endif
           if ( db2lig%atom_color(db2lig%coord_index(2,count)) .eq. 9) !assumes its neighbour is 9/10/11
     &        then
              ligi2=count
           endif
           if ( db2lig%atom_color(db2lig%coord_index(2,count)) .eq. 10)!assumes its neighbour is 9/10/11
     &        then
              ligi3=count 
           endif
           if ( db2lig%atom_color(db2lig%coord_index(2,count)) .eq. 11)!assumes its neighbour is 9/10/11
     &        then
              ligi4=count 
           endif
        enddo

        !determine geometry of covalent attachmed point 
        !by number of colored neighbours
        if ((ligi1 .eq. 0) .OR. (ligi2 .eq. 0)) then           
          write(OUTDOCK, outdockerr) molct, db2lig%refcod, 0,
     &                   0, 0.0, 'Ligand did not contain covalent 
     &                            colored atoms'
          cycle !skip to next molecule
        else if (ligi3 .eq. 0) then
          sp=1
        else if (ligi4 .eq. 0) then
           sp=2
        else
          sp=3
        endif

        !set covalent atom vdw type to a dummy atom to not get score penalty      
        !zero electrostatics and solvation
        db2lig%atom_vdwtype(db2lig%coord_index(2,ligi1))=25
        db2lig%atom_charge(db2lig%coord_index(2,ligi1))=0.0
        db2lig%atom_polsolv(db2lig%coord_index(2,ligi1))=0.0
        db2lig%atom_apolsolv(db2lig%coord_index(2,ligi1))=0.0

        !store original ligand coords
        do cli=1,3
           do clj=1,db2lig%total_coords
              cl(cli,clj)=db2lig%coords(cli,clj)
           enddo
        enddo

        ! d will contain the current fixed bond length
        d=options0%bondlen-options0%bondrange !set initial cov bond length
        !DEBUG for pdb printing
        atnum=1; resnum=1;
        ! loop 1 - iterate over bond length d 
        do s=1,lensteps 
           la=options0%bondang1-options0%bondang1range !initial cov bond angle
           ! loop 2 - iterate over la (bond angle)
           do k=1,ang1steps     
              ! loop 3 - iterate over lfi 360 deg
              do lfi=0,dih_final, dih_step
              !do lfi=0,340,20 
                 !convert to radians TODO this can be iterated in RAD
                 a=Pi-(la*deg_to_rad)
                 fi=lfi*deg_to_rad
                 ! get the location of the cov atom (in cov_loc)
                 call cov_bond(p1,p2,p3,a,fi,d,cov_loc)
                 atnum=atnum+1
                 !for each bond position sample lig orientations 
                 ! loop 4 - iterate over lb
                 lb=options0%bondang2-options0%bondang2range
                 do v=1,ang2steps
                    ! restore orig rigid coords
                    do tempcoord = db2lig%conf_coord(1, 1), 
     &                   db2lig%conf_coord(2, 1)
                       do count = 1, 3
                           rigid_part(count,tempcoord) = !note that you can use tempcoord here and not 1..rigcount only
     &                      cl(count, tempcoord)   ! because you assume the rigid coordinates are first.
                       enddo
                    enddo

                    !do n=1,18 !iterate over lfi1
                    do n=1,dih_num !iterate over lfi1
                    lfi1=lfi1mag 
                    !only orient the rigid part first
                    call cov_lig(p3,cov_loc,ligi1,ligi2,
     &                 ligi3,ligi4,rigid_part,rigcount,
     &                 lb,lfi1,sp,
     &                 (n.eq.1)) !adjust lb only the first time              
                  
                    !move rigid part into the transfm matrix for scoring
                    do tempcoord = db2lig%conf_coord(1, 1), 
     &                   db2lig%conf_coord(2, 1)
                       tempatom = db2lig%coord_index(2, tempcoord)
                       do count = 1, 3
                          db2lig%transfm_coords(count, tempatom) = 
     &                    rigid_part(count,tempcoord - 
     &                    db2lig%conf_coord(1, 1) + 1)
                       enddo
                    enddo

                    !conf 1 is the rigid part
                    !maxconfs and curr_match also 1 for the time being
                    temp_status = atom_out_of_grids(1, db2lig,
     &                  allminlimits, allmaxlimits) !do outside grid stuff
                    if (temp_status .eqv. .true.) then !outside grid, quit now
                      conf_status = OUTSIDEGRIDS !see status.h
                    else
                       maingrid = grids0%group_part(1, 1)
                       call calc_score_conf(1, options0, db2lig, grids0,
     &                      vdwmap0, phimap0, recdes0, ligdes0, gist0,
     &                      solvmap0,rec_des_r0, sra, srb, maingrid,
     &                      MAXTYV, asum, bsum, 
     &                      elescore, gistval, 
     &                      alscore, plscore, descore,rdscore, hescore,
     &                      ligscore%numscored,gistindexlist,lengistil,
     &                      gistindexlist,lengistil) ! TEB we need to fix this so gist and covalent are compatable
                       vdwscore = asum - bsum

c debug code to print rigid part and color b factor by energy
c                    if (vdwscore .gt. 100) then
c                      bfact = 100.0
c                   else
c                      bfact = vdwscore
c                   endif
c                   write(OUTDOCK,78) "MODEL ",resnum
c                    do p=1,rigcount
c                    write(6,76) 'HETATM',atnum+p,'  ','N ',
c     &                              ' RIG',resnum,'     ',
c     &                       rigid_part(1,p),' ',
c     &                       rigid_part(2,p),' ',
c     &                       rigid_part(3,p),' 1.00',bfact
c                    enddo
c                    resnum = resnum+1
c                    write(OUTDOCK,78) "ENDMDL"
c end debug code

                       !check for rigid part bumps
                       if (vdwscore .gt. options0%bmpvdwrigid) then !this group bumped
                          conf_status = BUMPED
                       else
                          conf_status = ALLOKAY
                       endif
                    endif

                    !if rigid part doesn't bump go ahead and score the rest of the 
                    !molecule
                    if (conf_status .eq. ALLOKAY ) then
                       match%nmatch = 1   
                       !restore orig coords before orienting
                       do cli=1,3
                          do clj=1,db2lig%total_coords
                             db2lig%coords(cli,clj)=cl(cli,clj)
                          enddo
                       enddo
                       !this acutally handles the rotations
                       call cov_lig(p3,cov_loc,ligi1,ligi2,ligi3,ligi4,
     &                           db2lig%coords,
     &                           db2lig%total_coords,lb,n*lfi1mag,
     &                           sp,.true.) 

                       !after orienting calculate the score
                      call calc_score_mol(dock_status, setstart, setend,
     &        db2lig, options0, ligscore, ligscoreeach, grids0,
     &        vdwmap0, phimap0, recdes0, ligdes0, gist0, solvmap0, 
     &        rec_des_r0,match, allminlimits, allmaxlimits,
     &        sra, srb, 
     &        MAXOR, db2lig%total_confs, db2lig%total_sets, 
     &        db2lig%total_atoms, db2lig%total_coords,
     &        MAXTYV) 

                    endif ! if bumped/ok 
                    enddo ! n ; loop5 (dihedral)
                    lb=lb+options0%bondang2step 
                 enddo ! lb ; loop4 (ang2)  
              enddo ! lfi ; loop3 (dihedral)
              la=la+options0%bondang1step !increment la
           enddo ! k ; loop2 (ang1)
           d = d+options0%bondstep !length update
        enddo ! s ; loop1 (lenth)

        !output to OUTDOCK and to output file
        if (dock_status .eq. ALLOKAY) then 
           do tempsave = 1, ligscore%savedcount
             write(OUTDOCK, outdockline) molct, db2lig%refcod, 
     &            ligscore%pose_reccode(tempsave), 
     &            match%nmatch, ligscore%numscored, 
     &            telaps, db2lig%total_heavy_atoms, 
     &            ligscore%setsave(tempsave), 
     &            ligscore%orientsave(tempsave),
     &            tempsave, cloudcount, ligscore%pose_es(tempsave), 
     &            ligscore%pose_gi(tempsave), 
     &            ligscore%pose_vs(tempsave), 
     &            ligscore%pose_ps(tempsave), 
     &            ligscore%pose_as(tempsave), 
     &            ligscore%pose_is(tempsave), 
     &            ligscore%pose_rs(tempsave), 
     &            ligscore%pose_ds(tempsave), 
     &            ligscore%pose_hs(tempsave), 
     &            ligscore%pose_score(tempsave)
           enddo
           do outcount = 1, ligscore%savedcount
              if (ligscore%pose_score(outcount) .lt. 
     &             options0%save_limit .and.
     &             ligscore%pose_score(outcount) .lt.
     &             options0%bmptotal) then
                 total_iter = 0
                 call mol2write(options0, db2lig, ligscore, match,
     &                molct, cloudcount, MAXOR, fdmol, outcount, 
     &                outcount, atomscore, grids0, vdwmap0, phimap0, 
     &                recdes0, ligdes0, gist0, solvmap0,rec_des_r0,
     &                sra, srb, maxtyv, allminlimits, allmaxlimits,
     &                time1, total_iter)
              endif
           enddo
        endif
        enddo !iend loop

        call gzclose(fdmol, istat) !close ligand output file

        return
        end subroutine run_cov_search
      
      end module cov_search
        
c---------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c
c       Copyright (C) 1991 Regents of the University of California
c                         All Rights Reserved.
c
