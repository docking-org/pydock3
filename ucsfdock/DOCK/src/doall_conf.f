c does everything for one specific conformation of a group.
c returns a status from status.h indicating success or what kind of failure
      module doall_conf
 
      implicit none
      contains

      function run_conf(conf, current_match, 
     &              db2lig, options0, grids0, which_grid,
     &              vdwmap0, phimap0, recdes0, ligdes0,gist0,
     &              solvmap0,rec_des_r0, match,
     &              vdwscore, elescore, gistscore, alscore, plscore, 
     &              descore, rdscore, hescore,
     &              allminlimits, allmaxlimits,
     &              sra, srb,
     &              maxor, maxconfs, maxatm, maxatmpos, maxtyv, 
     &              numscored,
     &              gistgridindcur,lengistgicu,
     &              gistgridindcon,lengistgico) 
     &              result(confstat)

      use status
      use db2type
      use optionstype
      use gridstype
      use phimaptype
      use gisttype
      use solvmaptype
      use vdwmaptype
      use matchtype
      use transfm_conf
      use check_conf_og
      use score_conf

c XXX: GFortran needs these to come first (TS)
c sizing constants from max.h
      integer, intent(in) :: maxor !the size of many arrays that store information about
                    !orientations, is checked several times
      integer, intent(in) :: maxconfs !the size of any array having to do with conformations
      integer, intent(in) :: maxatm !how many atoms max can be in a ligand
      integer, intent(in) :: maxatmpos !how many positions all ligand atoms can take
      integer, intent(in) :: maxtyv !how many vdw atom types there are
c input variables    
      integer, intent(in) :: conf !which conf to do
      integer, intent(in) :: current_match !which match set of trans/rots to do
      type(db2), intent(inout) :: db2lig !ligand information
      type(options), intent(in) :: options0 !useful options
      type(flexgrids), intent(inout) :: grids0 !precomp values
      integer, intent(in) :: which_grid !which grid to use
      type(vdwmap), intent(in) :: vdwmap0
      type(phimap), intent(in) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(in) ::gist0, rec_des_r0
      type(solvmap), intent(in) :: solvmap0
      type(matcht), intent(in) :: match
      real, intent(inout) :: vdwscore, elescore, gistscore, alscore, 
     &    plscore, descore,rdscore, hescore !descore is receptor desolvation,
          !hescore is receptor hydrophobic effect
      real, intent(in), dimension(3) :: allminlimits, allmaxlimits !grid limits in xyz
      real, intent(in), dimension(maxtyv) :: sra, srb !square roots of vdw parameters per vdwtype
c return variables 
      integer, intent(inout) :: numscored !keeps track of how many atom positions are evaluated

      ! XXX: GFortran chokes on this (TS)
      !integer, intent(out) :: confstat !return value 
      integer :: confstat !return value

      logical temp_status !temporary status holder
      real asum, bsum !temp holders for vdw energy

      integer count, tempcoord, tempatom

c     ! these arrays are used so that we are not double counting voxals. 
      integer, intent(inout) :: gistgridindcur(5000), ! the current index list to used to determin which voxals are displaced by this conf
     &                          gistgridindcon(5000)  ! the connected conformation, so we know what was displaced last time. 
 
      integer, intent(inout) :: lengistgicu, lengistgico ! this is the curent length of list grid point already asigned to the molecule


 77   format(A,I5,A,A3,A,I6,A,F7.3,A,F7.3,A,F7.3)

      !check if this is a covalent run 
      if (options0%dockovalent .eqv. .false.) then
         !transform this conf into transfm_coords
         call run_transfm_conf(conf, current_match, db2lig,
     &        match, maxor)     !match-picked transformation
      else
         do tempcoord = db2lig%conf_coord(1, conf), 
     &        db2lig%conf_coord(2, conf)
            tempatom = db2lig%coord_index(2, tempcoord)
            do count = 1, 3
               db2lig%transfm_coords(count, tempatom) = 
     &         db2lig%coords(count,tempcoord)
            enddo
         enddo
      endif
      !print *, 'which_grid in doall_conf.f: ', which_grid
      !print *, 'size of grid0: ', size(grids0%griddata)
      ! XXX: GFortran Logicals at (1) must be compared with .eqv. 
      !      instead of .eq.
      !if (temp_status .eqv. .true.) then !outside grid, quit now
c      write(6,*) "I AM in doall_confs."
c      write(6,*) "temp_status is", temp_status
      temp_status = atom_out_of_grids(conf, db2lig,
     &    allminlimits, allmaxlimits) !do outside grid stuff
      if (temp_status .eqv. .true.) then !outside grid, quit now
        confstat = OUTSIDEGRIDS !see status.h
c      write(6,*) "temp_status is", temp_status
      else
c        write(6,*) "running calc_score_conf"
c       write(*,*) "conf =", conf
c       write(*,*) "start,lengistgi=",lengistgi
        call calc_score_conf(conf, options0, db2lig, grids0,
     &      vdwmap0, phimap0, recdes0, ligdes0, gist0, 
     &      solvmap0,rec_des_r0, sra, srb, which_grid,
     &      maxtyv,
     &      asum, bsum, 
     &      elescore, gistscore, alscore, plscore, descore,rdscore,
     &      hescore,numscored,
     &      gistgridindcur,lengistgicu,
     &      gistgridindcon,lengistgico)
c       write(*,*) "stop,lengistgi=",lengistgi
        vdwscore = asum - bsum
        if (vdwscore .gt. options0%bmpvdw) then !this group bumped
          confstat = BUMPED
        else
          confstat = ALLOKAY !see status.h, everything went well
        endif
      endif

      return
      end function run_conf

      end module doall_conf
