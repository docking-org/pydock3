
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

      real function sim_energy(v, best_energy, confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)


      use db2type
      use optionstype
      use gridstype
      use phimaptype
      use gisttype
      use solvmaptype
      use vdwmaptype
      use phiscore
      use gistscore
      use rec_des_score
      use chemscore
      use solvation_score
      use matchtype
      use ligscoretype
      use atomscoretype
      use status

      implicit none

      integer, intent(in) :: maxtyv !how many vdw atom types there are
      integer, intent(in) :: confset !the confset that the caller would like
                                  !scored
      type(matcht), intent(inout) :: match
      integer, intent(in):: matchnum !for getting rotation matrices
      integer, intent(in) :: MAXOR ! max orientations
      character (len=255), intent(in) :: reccode !flexible receptor code for this pose
      type(options), intent(in) :: options0 !useful options
      type(db2), intent(in) :: db2lig !ligand information
      type(flexgrids), intent(inout) :: grids0 !precomp values
      type(vdwmap), intent(in) :: vdwmap0
      type(phimap), intent(in) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(in) :: gist0, rec_des_r0
      type(solvmap), intent(in) :: solvmap0
      real, dimension(maxtyv), intent(in) :: sra, srb !square roots of
                                                !vdw parameters per vdwtype

      real, intent(inout) ::  vscore, escore, gscore,
     &           ascore, pscore,
     &           rscore, rdscore, hscore, setscore
c     &           rscore, dscore, hscore, setscore

      integer, intent(inout) :: min_status
      real, dimension(6) :: v
      real, intent(inout) :: best_energy
      real :: rigid_min
      real, dimension(3), intent(in) :: allminlimits, allmaxlimits


      integer atomstart, atomend !temp/loop variables

c     convert Euler angles to rotation matrix
      call euler2rot(v, match, matchnum)

c     Compute molecular energy during minimization
      best_energy = rigid_min(confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode, 
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
c     &           rscore, dscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)


      sim_energy = best_energy

c     call doflush(6)
      !stop 46
      return
      end
