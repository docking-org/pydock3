      module score_mol

      implicit none
      contains

      subroutine calc_score_mol(dock_status, setstart, setend,
     &    db2lig, options0, ligscore, ligscoreeach, grids0,
     &    vdwmap0, phimap0, recdes0, ligdes0, gist0,
     &    solvmap0,rec_des_r0, match, allminlimits, allmaxlimits,
     &    sra, srb, 
     &    maxor, maxconfs, maxsets, maxatm, maxatmpos,
     &    maxtyv)

      use status
      use db2type
      use optionstype
      use ligscoretype
      use gridstype
      use phimaptype
      use gisttype
      use solvmaptype
      use vdwmaptype
      use matchtype
      use score_sets

c XXX: GFortran needs these first
c constants, set in max.h and passed in
      integer, intent(in) :: maxor !the size of many arrays that store information about 
                    !orientations, is checked several times 
      integer, intent(in) :: maxconfs !the size of any array having to do with conformations  
      integer, intent(in) :: maxsets !how many max sets can be in a ligand file
      integer, intent(in) :: maxatm !how many atoms max can be in a ligand
      integer, intent(in) :: maxatmpos !how many positions all ligand atoms can take
      integer, intent(in) :: maxtyv !how many vdw atom types there are

      integer dock_status !see status.h for enumerated types
      integer setstart, setend !lets only certain sets (entire confs) be scored
      type(db2), intent(inout) :: db2lig
      type(options), intent(in) :: options0
      type(ligscoret), intent(inout) :: ligscore
      type(ligscoret), intent(inout) :: ligscoreeach
      type(flexgrids), intent(inout) :: grids0 !precomp values
      type(vdwmap), intent(inout) :: vdwmap0
      type(phimap), intent(inout) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(inout) :: gist0, rec_des_r0
      type(solvmap), intent(inout) :: solvmap0
      type(matcht), intent(inout) :: match
      real, dimension(3), intent(in) :: allminlimits !grid limits in xyz
      real, dimension(3), intent(in) :: allmaxlimits !grid limits in xyz
      real, dimension(maxtyv), intent(in) :: sra, srb !square roots of vdw parameters per vdwtype
c temporary variables declared here
      integer temp_status !for passing into score_sets
      integer curmatch !copying rotation matrices around

      if (match%nmatch .gt. maxor) then !adjust down since not all rotations save
        match%nmatch = maxor
      endif
      do curmatch = 1, match%nmatch !copy rotmatch into rot, etc
        temp_status = NOMATCH
        call calc_score_sets(curmatch, temp_status, setstart, setend,
     &      db2lig, options0, ligscore, ligscoreeach, grids0,
     &      vdwmap0, phimap0, recdes0, ligdes0, gist0, solvmap0, 
     &      rec_des_r0, match,
     &      allminlimits, allmaxlimits,
     &      sra, srb, 
     &      maxor, maxconfs, maxsets, maxatm, maxatmpos,
     &      maxtyv)
c     the meaning here is that if this bumped out but a previous orientation
c     found a flexible conformation without bumps, don't suddenly say bumped,
c     everything is still fine and the pointers are still correctly set, scored
        if (dock_status .ne. ALLOKAY) then
          dock_status = temp_status
        endif
      enddo
      return
      end subroutine calc_score_mol

      end module score_mol
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
