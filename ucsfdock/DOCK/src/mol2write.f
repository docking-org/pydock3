! writes a mol2 file for one saved pose
      subroutine mol2write(options0, db2lig, ligscore, match, 
     &    lignumber, cloudnum, MAXOR, fdmol, count1, outrank,
     &    atomscore, grids0, vdwmap0, phimap0, recdes0,
     &    ligdes0, gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &    allminlimits, allmaxlimits, time1, total_iterations)


      use db2type
      use ligscoretype
      use matchtype
      use mol2sav
      use optionstype
      use transfm_conf
      use atomscoretype
      use gridstype
      use phimaptype
      use gisttype
      use solvmaptype
      use vdwmaptype
      use rescore_byatom
      use gisttype
      use status

      implicit none

      type(options), intent(in) :: options0
      type(db2), intent(inout) :: db2lig
      type(ligscoret), intent(inout) :: ligscore
      type(matcht), intent(inout) :: match
      integer, intent(in) :: lignumber !ligand number (this run), the number of ligands searched on this run
      integer, intent(in) :: cloudnum !cloud number, or just 1 if clouds not used
      integer, intent(in) :: MAXOR ! max orientations
      integer (kind=8), intent(inout) :: fdmol !ligand output file handle (mol2.gz)
      integer, intent(in) :: count1 ! the count1-th saved pose
      integer, intent(in) :: outrank !output ranking
      type(atomscoret), intent(inout) :: atomscore
      type(flexgrids), intent(inout) :: grids0 !precomp values
      type(vdwmap), intent(in) :: vdwmap0
      type(phimap), intent(in) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(in) :: gist0, rec_des_r0
      type(solvmap), intent(in) :: solvmap0
      ! GForran likes these before they're used
      integer, intent(in) :: maxtyv !how many vdw atom types there are
      real, dimension(maxtyv), intent(in) :: sra, srb !square roots of vdw params
      real :: time1

! temporary variables used here
      character (len=160) molbuf !output buffer
      integer istat, bondcount !file status, temp counters
      integer setcount, confcount, conf, atomcount, thisconf
      integer thisconfcount !more temp counters
      integer matchnum !for getting rotation matrices
      integer arbcount !for walking through arbitrary lines of output
      integer tempcount, tempcoord
      integer :: i,j
      real :: time3, time_elapsed

! for minimization
      integer :: min_status
      real :: best_energy ! for minimization
      integer, intent(inout) :: total_iterations !cumulative number of minimization steps performed
      integer :: temp_combo_num
      character (len=255) :: reccode
      real :: setvs, setes,setgi, setas, setrs,
     &       setps, setrd, seths, setscore !scores (whole position)
      real, intent(in), dimension(3) :: allminlimits, allmaxlimits !grid limits in xyz

      character (len=*), parameter :: outdockline_min = '(i6,1x,a16,1x,
     &       a16,
     &       1x,i8,1x,i10,1x,f7.2,1x,i3,1x,i9,1x,i9,1x,i6,3x,i3,2x,
     &       f7.2,1x,f7.2,1x,f7.2,f7.2,1x,f7.2,1x,f7.2,1x,f7.2,1x,f7.2,
     &       1x,f7.2,f10.2)'

      character (len=255) :: min_info

      !write(6,*),"I AM IN MOL2_WRITE"
      setcount = ligscore%setsave(count1) ! gets this particular conformation set
      matchnum = ligscore%orientsave(count1)
      total_iterations = 0
      temp_combo_num = 1

      best_energy = 1000

      if (options0%mol2_minimizer .ne. 1) then
       if (options0%dockovalent .eqv. .false.) then
         do confcount = 1, db2lig%set_conf_len(setcount) !sets of complete conformations
            conf = db2lig%set_conf(setcount, confcount)
            !each conf needs to transform the coordinates to the right place
            call run_transfm_conf(conf, matchnum, db2lig,
     &           match, MAXOR)  !match-picked transformation
         enddo
       else
         do tempcoord = 1, db2lig%total_atoms
            do tempcount = 1, 3
               db2lig%transfm_coords(tempcount, tempcoord) =
     &              db2lig%save_cov_coords(count1,tempcount,tempcoord)
            enddo
         enddo
       endif

         
      !now atoms in transfm_coords are in correct position for writing out 
      !also good in case they need re-scored for additional output

       if (options0%per_atom_scores) then !have to score each atom
c       call doflush(6)
        call recalc_score_byatom(options0, db2lig, grids0, vdwmap0, 
     &      phimap0, recdes0, ligdes0, gist0, solvmap0, rec_des_r0, 
     &      sra, srb, maxtyv,  atomscore, setcount, 
     &      ligscore%flexnumsave(count1))
       endif

       !call gzwrite(fdmol, " ", istat)
       !write(*, *) 'istat', istat
       !write(*, *) "I AM HERE..", count1
       !write(*, *) "I AM HERE..", count1, ligscore%pose_reccode(count1)
       !call doflush(6)
      !write(molbuf, MOL2AB) '######## MINIMIZATION: no'
       write(molbuf, MOL2B) 'Name', db2lig%refcod
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2B) 'Protonation', db2lig%protcod
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2A) 'SMILES', db2lig%smiles
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2A) 'Long Name', db2lig%molname
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2A) 'FlexRecCode', ligscore%pose_reccode(count1)
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2I) 'Number', lignumber
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2A) 'Ligand Source File', options0%ligfil
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2I) 'Rank', outrank
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2I) 'Setnum', setcount
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2I) 'Matchnum', matchnum
       call gzwrite(fdmol, molbuf, istat)
cc I (Trent Balius) added a colon.
cc I believe this is for debugging purposes. Can it be commented out? 
c      ! output oxr  <-- looks like this is used for printing out receptor
c      ! spheres.
       do i=1,4
         write(molbuf, MOL2X3) 'OXR', match%xor(1,i,matchnum), 
     &          match%xor(2,i,matchnum), match%xor(3,i,matchnum)
         call gzwrite(fdmol, molbuf, istat)
       enddo
c      ! output oxs <-- looks like this is used for printing out ligand
c      ! spheres. 
       do i=1,4
         write(molbuf, MOL2X3) 'OXS', match%xos(1,i,matchnum), 
     &          match%xos(2,i,matchnum), match%xos(3,i,matchnum)
         call gzwrite(fdmol, molbuf, istat)
       enddo

       write(molbuf, MOL2I) 'Cloud', cloudnum
       call gzwrite(fdmol, molbuf, istat)
      !write(molbuf, MOL2F) 'Partition Function', partfunc
      !call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2F) 'Electrostatic', ligscore%pose_es(count1)
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2F) 'Gist', ligscore%pose_gi(count1)
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2F) 'Van der Waals', ligscore%pose_vs(count1)
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2F) 'Ligand Polar Desolv', 
     &    ligscore%pose_ps(count1)
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2F) 'Ligand Apolar Desolv', 
     &    ligscore%pose_as(count1)
       call gzwrite(fdmol, molbuf, istat)
c       write(molbuf, MOL2F) 'Internal Energy', ligscore%pose_is(count1)
       write(molbuf, MOL2F) 'Total Strain', ligscore%pose_is(count1)
       call gzwrite(fdmol, molbuf, istat)
c       write(molbuf, MOL2F) 'Receptor Energy', ligscore%pose_rs(count1)
       write(molbuf, MOL2F) 'Max Strain', ligscore%pose_rs(count1)
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2F) 'Receptor Desolvation', 
     &    ligscore%pose_rd(count1)
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2F) 'Receptor Hydrophobic', 
     &    ligscore%pose_hs(count1)
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2F) 'Total Energy', ligscore%pose_score(count1)
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2F) 'Ligand Charge', db2lig%total_charge
       call gzwrite(fdmol, molbuf, istat)
       do arbcount = 1, db2lig%mlines - 4
        write(molbuf, MOL2A) 'Arbitrary', db2lig%arbitrary(arbcount)
        call gzwrite(fdmol, molbuf, istat)
       enddo
       write(molbuf, MOL2F) 'Ligand Energy', 
     &    db2lig%set_energy(setcount)
       call gzwrite(fdmol, molbuf, istat)
       if (options0%per_atom_scores) then !have to write out scores for each atom
        write(molbuf, '(a,a)') '# num electro vdwattr vdwrepu polsolv ',
     &       'apolsol rdesolv rhydroe totalscore fullpol fullapo' !write header
        call gzwrite(fdmol, molbuf, istat)
      !write(*,*) "I AM HERE writemol2"
        do atomcount = 1, db2lig%total_atoms
          write(molbuf, MOL2ASD) atomcount, 
     &        atomscore%atom_es(atomcount), 
     &        atomscore%atom_vas(atomcount), 
     &        atomscore%atom_vrs(atomcount), 
     &        atomscore%atom_ps(atomcount), 
     &        atomscore%atom_as(atomcount), 
     &        atomscore%atom_rd(atomcount),
     &        atomscore%atom_hs(atomcount), 
     &        atomscore%atom_score(atomcount), 
     &        -db2lig%atom_polsolv(atomcount), !these are - since full score
     &        -db2lig%atom_apolsolv(atomcount) !is solvation, we want desolv
          call gzwrite(fdmol, molbuf, istat)
        enddo
       endif
      !write(*,*) "I AM HERE writemol2"
       write(molbuf, MOL280) " " !empty line
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2MOL) '@<TRIPOS>MOLECULE'
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2C) db2lig%refcod, db2lig%protcod
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2I5) db2lig%total_atoms, db2lig%total_bonds, 
     &    INT(db2lig%total_charge), 0, 0
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL280) " " !empty line
       call gzwrite(fdmol, molbuf, istat) 
       call gzwrite(fdmol, molbuf, istat)
       call gzwrite(fdmol, molbuf, istat)
       write(molbuf, MOL2AB) '@<TRIPOS>ATOM'
       call gzwrite(fdmol, molbuf, istat)
       do atomcount = 1, db2lig%total_atoms
        write(molbuf, mol2atom) db2lig%atom_num(atomcount), 
     &     db2lig%atom_name(atomcount), 
     &     db2lig%transfm_coords(1, atomcount),
     &     db2lig%transfm_coords(2, atomcount), 
     &     db2lig%transfm_coords(3, atomcount), 
     &     db2lig%atom_type(atomcount),
     &     1,'LIG1', db2lig%atom_charge(atomcount)
        call gzwrite(fdmol, molbuf, istat)
       enddo

       write(molbuf, MOL2AB) '@<TRIPOS>BOND'
       call gzwrite(fdmol, molbuf, istat)
       do bondcount = 1, db2lig%total_bonds !easy because this is always the same
        write(molbuf, MOL2BOND) db2lig%bond_num(bondcount), 
     &      db2lig%bond_start(bondcount), db2lig%bond_end(bondcount),
     &      db2lig%bond_type(bondcount)
        call gzwrite(fdmol, molbuf, istat)
       enddo

       !call gzwrite(fdmol, " ", istat)
      
c      call doflush(6)

      

      else if (options0%mol2_minimizer .gt. 0 .and.
     &           best_energy .le. options0%min_cut) then

         best_energy=ligscore%pose_score(count1)
         min_status = ALLOKAY
         setvs = 0.0
         setes = 0.0
         setgi = 0.0
         setas = 0.0
         setps = 0.0
         setrs = 0.0
         setrd = 0.0
         seths = 0.0
         setscore = 0.0
         min_info = 'YES'

         do i = 1, 3
c           store rotation matrix before minimization
              do j = 1, 3
                  match%rotab(i,j,matchnum) = match%rot(i,j,matchnum)
              enddo
              match%comrab(i,matchnum) = match%comr(i, matchnum)
         enddo

         setcount = ligscore%setsave(count1) ! gets this particular conformation set
         confcount = db2lig%set_conf_len(setcount) !sets of complete conformations

          call dockmin_sim(best_energy, setcount,options0, db2lig,
     &           match, matchnum,
     &           MAXOR, ligscore%pose_reccode(count1),
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &            setvs,
     &            setes,
     &            setgi,
     &            setas,
     &            setps,
     &            setrs,
     &            setrd,
     &            seths,
     &            setscore, min_status,
     &            allminlimits, allmaxlimits, total_iterations)

          call get_time(time3)

          time_elapsed = time3 - time1
          !write(*,*) 'total minimization steps : ', total_iterations
          if (min_status .eq. OUTSIDEGRIDS) then
            min_info = 'YES'
          endif
          !print *, 'min_status=', min_info

c         call doflush(6)
          if (options0%dockovalent .eqv. .false.) then
            do confcount = 1, db2lig%set_conf_len(setcount) !sets of complete conformations
                conf = db2lig%set_conf(setcount, confcount)
                !each conf needs to transform the coordinates to the right place
                call run_transfm_conf(conf, matchnum, db2lig,
     &                                match, MAXOR)
            enddo
          else
              do tempcoord = 1, db2lig%total_atoms
                do tempcount = 1, 3
                    db2lig%transfm_coords(tempcount, tempcoord) =
     &                   db2lig%save_cov_coords(
     &                          count1,tempcount,tempcoord)
                enddo
              enddo
          endif

          write(6, outdockline_min) lignumber, db2lig%refcod,
     &            ligscore%pose_reccode(count1),
     &            match%nmatch, ligscore%numscored,
     &            time_elapsed, db2lig%total_heavy_atoms,
     &            setcount,
     &            matchnum,
     &            outrank, cloudnum, setes,
     &            setgi,
     &            setvs,
     &            setps,
     &            setas,
     &            db2lig%set_energy(setcount),
     &            setrs,
     &            setrd,
     &            seths,
     &            setscore

c         call doflush(6)

          if (options0%per_atom_scores) then !have to score each atom
c           call doflush(6)
            call recalc_score_byatom(options0, db2lig, grids0, vdwmap0,
     &          phimap0, recdes0, ligdes0, gist0, solvmap0, rec_des_r0, 
     &          sra, srb, maxtyv,  atomscore, setcount,
     &          ligscore%flexnumsave(count1))
          endif

          write(molbuf, MOL2B) 'MinStatus', min_info
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2B) 'Name', db2lig%refcod
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2B) 'Protonation', db2lig%protcod
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2A) 'SMILES', db2lig%smiles
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2A) 'Long Name', db2lig%molname
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2A) 'FlexRecCode',
     &             ligscore%pose_reccode(count1)
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2I) 'Number', lignumber
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2A) 'Ligand Source File', options0%ligfil
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2I) 'Rank', outrank
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2I) 'Setnum', setcount
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2I) 'Matchnum', matchnum
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2I) 'Cloud', cloudnum
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2F) 'Electrostatic',
     &           setes
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2F) 'Gist',
     &          setgi
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2F) 'Van der Waals',
     &         setvs
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2F) 'Ligand Polar Desolv',
     &          setps
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2F) 'Ligand Apolar Desolv',
     &         setas
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2F) 'Internal Energy',
     &              db2lig%set_energy(setcount)
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2F) 'Receptor Energy',
     &             setrs
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2F) 'Receptor Desolvation',
     &              setrd
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2F) 'Receptor Hydrophobic',
     &         seths
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2F) 'Total Energy',
     &              setscore
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2F) 'Ligand Charge', db2lig%total_charge
          call gzwrite(fdmol, molbuf, istat)
          do arbcount = 1, db2lig%mlines - 4
                write(molbuf, MOL2A) 'Arbitrary',
     &                    db2lig%arbitrary(arbcount)
                call gzwrite(fdmol, molbuf, istat)
          enddo
          write(molbuf, MOL2F) 'Ligand Energy',
     &        db2lig%set_energy(setcount)
          call gzwrite(fdmol, molbuf, istat)
          if (options0%per_atom_scores) then !have to write out scores for each atom
                write(molbuf, '(a,a)')
     &             '# num electro vdwattr vdwrepu polsolv ',
     &             'apolsol rdesolv rhydroe totalscore fullpol fullapo' !write header
                call gzwrite(fdmol, molbuf, istat)
                !write(*,*) "I AM HERE writemol2"
                do atomcount = 1, db2lig%total_atoms
                    write(molbuf, MOL2ASD) atomcount,
     &                  atomscore%atom_es(atomcount),
     &                  atomscore%atom_vas(atomcount),
     &                  atomscore%atom_vrs(atomcount),
     &                  atomscore%atom_ps(atomcount),
     &                  atomscore%atom_as(atomcount),
     &                  atomscore%atom_rd(atomcount),
c     &                  atomscore%atom_ds(atomcount),
     &                  atomscore%atom_hs(atomcount),
     &                  atomscore%atom_score(atomcount),
     &                  -db2lig%atom_polsolv(atomcount), !these are - since full score
     &                  -db2lig%atom_apolsolv(atomcount) !is solvation, we want desolv
                    call gzwrite(fdmol, molbuf, istat)
                enddo
          endif
          !write(*,*) "I AM HERE writemol2"
          write(molbuf, MOL280) " " !empty line
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2MOL) '@<TRIPOS>MOLECULE'
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2C) db2lig%refcod, db2lig%protcod
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2I5) db2lig%total_atoms,
     &           db2lig%total_bonds,
     &           INT(db2lig%total_charge), 0, 0
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL280) " " !empty line
          call gzwrite(fdmol, molbuf, istat)
          call gzwrite(fdmol, molbuf, istat)
          call gzwrite(fdmol, molbuf, istat)
          write(molbuf, MOL2AB) '@<TRIPOS>ATOM'
          call gzwrite(fdmol, molbuf, istat)
          do atomcount = 1, db2lig%total_atoms
                write(molbuf, mol2atom) db2lig%atom_num(atomcount),
     &              db2lig%atom_name(atomcount),
     &              db2lig%transfm_coords(1, atomcount),
     &              db2lig%transfm_coords(2, atomcount),
     &              db2lig%transfm_coords(3, atomcount),
     &              db2lig%atom_type(atomcount),
     &              1,'LIG1', db2lig%atom_charge(atomcount)
                call gzwrite(fdmol, molbuf, istat)
          enddo
          write(molbuf, MOL2AB) '@<TRIPOS>BOND'
          call gzwrite(fdmol, molbuf, istat)
          do bondcount = 1, db2lig%total_bonds !easy because this is always the same
            write(molbuf, MOL2BOND) db2lig%bond_num(bondcount),
     &             db2lig%bond_start(bondcount),
     &              db2lig%bond_end(bondcount),
     &             db2lig%bond_type(bondcount)
            call gzwrite(fdmol, molbuf, istat)
          enddo

         !print *, 'Finished writing MINIMIZED mol2 file'

      endif
      !write(*,*) "I AM HERE writemol2"

      return
      end subroutine mol2write
