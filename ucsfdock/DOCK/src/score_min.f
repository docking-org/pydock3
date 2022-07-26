c this is a scoring function, to give by-atom energies, useful in some cases
      module score_min

      implicit none
      contains

      subroutine recalc_score_min(conf,
     &           options0, db2lig, grids0,vdwmap0, phimap0,
     &           recdes0, ligdes0, gist0, solvmap0, rec_des_r0, 
     &           sra, srb, which_grid, maxtyv, vdwscore, eleceel, 
     &           gistval, solva, solvp, rdval, hescore,
     &           numscored, confstat)

      use status
      use db2type
      use optionstype
      use gridstype
      use phimaptype
      use gisttype
      use solvmaptype
      use vdwmaptype
      use phiscore
      use rec_des_score
      use solvation_score
      use atomscoretype
      use matchtype
      use ligscoretype
      use atomscoretype
      use status
      use chemscore
      
      integer, intent(in) :: conf !the conf that the caller would like scored
      type(options), intent(in) :: options0 !useful options
      integer, intent(in) :: maxtyv !how many vdw atom types there are
      type(db2), intent(in) :: db2lig !ligand information
      type(flexgrids), intent(inout) :: grids0 !precomp values
      type(vdwmap), intent(in) :: vdwmap0
      type(phimap), intent(in) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(in) :: gist0, rec_des_r0
      type(solvmap), intent(in) :: solvmap0
      real, dimension(maxtyv), intent(in) :: sra, srb !square roots of vdw parameters per vdwtype
      integer, intent(inout) :: confstat !return value

      real, intent(inout) :: vdwscore, eleceel, gistval, rdval
      real, intent(inout) :: solva, solvp, hescore !scores by
c      real, intent(inout) :: solva, solvp, solvd, hescore !scores by
                                                          !reference
      integer :: atomstart, atomend, atomnow, atomindex !temp/loop variables
      real :: potent(db2lig%total_atoms) !the electrostatic potential at each atom, 
            !set in phipdb and used in charge but then discarded, so only needed
            !here, doesn't need to be global
      real rec_d_energy(db2lig%total_atoms)
      real recdes_energy(db2lig%total_atoms) !receptor desolv energy, similar to potent
      real ligdes_energy(db2lig%total_atoms) !ligand desolv potential
      real :: vdwasum, vdwbsum
      real :: solvatemp, solvptemp, rdtemp
c      real :: solvatemp, solvptemp, solvdtemp
c      real :: gist !hydrophobic effect
      real :: hscore, hscoretemp !hydrophobic effect
      real :: rec_energy !sum receptor energies
      integer :: numscored !keeps track of how many atom positions are evaluated
      integer, intent(in) :: which_grid ! the grid number 
      integer :: count_grid !which grids to use
      integer :: confcount, atomcount, coordcount !which conf currently being examined


      vdwasum = 0.0
      vdwbsum = 0.0
      eleceel = 0.0
      solva = 0.0
      solvp = 0.0
      rdval = 0.0
c      solvd = 0.0
      !which_grid = 1
      numscored = 0
      atomstart = db2lig%conf_coord(1, conf) !atom index of start
      atomend = db2lig%conf_coord(2, conf) !same but end (2)
      !print *, 'atomstart=', atomstart
      !print *, 'atomend=', atomend
      call calc_chemscore(atomstart, atomend,
     &    db2lig%coord_index, db2lig%transfm_coords, vdwmap0%offset,
     &    vdwmap0%invgrid, vdwmap0%grdpts,
     &    sra, srb, db2lig%atom_vdwtype,
     &    grids0%griddata(which_grid)%abval_precomp, vdwmap0%ngrd,
     &    db2lig%total_atoms, db2lig%total_coords, maxtyv,
     &    vdwasum, vdwbsum, numscored)


c      if ((vdwasum - vdwbsum) .le. options0%bmpvdw) then !otherwise bump this group
c      if ((vdwasum - vdwbsum) .le. 10) then !otherwise bump this group
            if (options0%ldesolv .le. 2) then !none, full or old GB partial
              call calc_solvation_score(atomstart, atomend,
     &            db2lig%atom_vdwtype, db2lig%atom_polsolv,
     &            db2lig%atom_apolsolv, db2lig%coord_index,
     &            db2lig%transfm_coords, solvmap0%sperang,
     &            solvmap0%smax, solvmap0%scadif,
     &            grids0%griddata(which_grid)%solgrd_precomp,
     &            grids0%griddata(which_grid)%hsolgrd_precomp,
     &            options0%hsolflag, options0%ldesolv,
     &            db2lig%total_atoms, db2lig%total_coords, solvp, solva)
c             solvation_score sets solvp and solva
            else !ldesolv = 3 = PB ligand desolvation
              call phipdb(atomstart, atomend, db2lig%coord_index,
     &            db2lig%transfm_coords,
     &            ligdes0%oldmid, ligdes0%phiscale, ligdes0%goff,
     &            grids0%griddata(which_grid)%ligdes_precomp,
     &            db2lig%total_atoms,
     &            db2lig%total_coords, ligdes0%nsize, ligdes_energy) !set to kT
              solva = 0.0
              solvp = charge_squared(atomstart, atomend,
     &            db2lig%coord_index,
     &            db2lig%atom_charge, ligdes_energy, db2lig%total_atoms,
     &            db2lig%total_coords)
            endif
            eleceel = 0.0
            if (allocated(grids0%griddata(which_grid)%phimap_precomp))
     &          then !this check is to not read data when implicit0 grids are used
              call phipdb(atomstart, atomend, db2lig%coord_index,
     &            db2lig%transfm_coords,
     &            phimap0%oldmid, phimap0%phiscale, phimap0%goff,
     &            grids0%griddata(which_grid)%phimap_precomp,
     &            db2lig%total_atoms,
     &            db2lig%total_coords, phimap0%nsize, potent) !sets potent() array for use next
              !print *, 'electrostatics before calc=', eleceel

              eleceel = charge(atomstart, atomend, db2lig%coord_index,
     &            db2lig%atom_charge, potent, db2lig%total_atoms,
     &            db2lig%total_coords)
              !print *, 'electrostatics after calc=', eleceel
              !write(6,*) "calc_score:I AM HERE(4)"
            endif
            rdval = 0.0
            if (allocated(grids0%griddata(which_grid)%rec_des_precomp))
     &          then
              call dx_pb(db2lig%atom_vdwtype, atomstart, atomend,
     &            db2lig%coord_index,
     &            db2lig%transfm_coords, rec_des_r0%orgin,
     &            rec_des_r0%gistspace, rec_des_r0%xnsize,
     &            rec_des_r0%ynsize,rec_des_r0%znsize,
     &            grids0%griddata(which_grid)%rec_des_precomp,
     &            grids0%griddata(which_grid)%hrec_des_precomp,
     &            db2lig%total_atoms, db2lig%total_coords, rec_d_energy)
              rdval = convert_rdsum(atomstart, atomend,
     &            db2lig%coord_index, rec_d_energy,
     &            db2lig%atom_vdwtype, db2lig%total_atoms,
     &            db2lig%total_coords, rec_des_r0%gistspace)
            endif
c            solvd = 0.0 !receptor desolvation is 0. if not requested
c            if (options0%receptor_desolv) then !do receptor desolvation calculation
c              if (allocated(grids0%griddata(which_grid)%recdes_precomp))
c     &            then !the above check is for implicit0 grids
c                call phipdb(atomstart, atomend, db2lig%coord_index,
c     &            db2lig%transfm_coords,
c     &            recdes0%oldmid, recdes0%phiscale, recdes0%goff,
c     &            grids0%griddata(which_grid)%recdes_precomp,
c     &            db2lig%total_atoms,
c     &            db2lig%total_coords, recdes0%nsize, recdes_energy) !set to kT
c                solvd = convert_sum(atomstart, atomend,
c     &             db2lig%coord_index,
c     &            recdes_energy, db2lig%atom_vdwtype,
c     &            db2lig%atom_heavy_valence, db2lig%total_atoms,
c     &            db2lig%total_coords)
c              !write(6,*) "calc_score:I AM HERE(6)"
c              endif
c            endif

            hescore = 0.0 !receptor hydrophobic effecti s 0. if not requested

c      endif
c      endif
      vdwscore=vdwasum-vdwbsum
      confstat = ALLOKAY

c      if (vdwscore .gt. options0%bmpvdw) then
c      if (vdwscore .gt. 10) then
c        confstat = BUMPED
c        !print *, 'confstat=', confstat
c      else 
c        confstat = ALLOKAY
c      end if


      return
      end subroutine recalc_score_min
 
      end module score_min
