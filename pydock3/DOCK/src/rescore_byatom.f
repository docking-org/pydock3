c this is a scoring function, to give by-atom energies, useful in some cases
      module rescore_byatom

      implicit none
      contains

      subroutine recalc_score_byatom(options0, db2lig, grids0,
     &    vdwmap0, phimap0, recdes0, ligdes0, gist0, solvmap0,
     &    rec_des_r0, sra, srb, 
     &    maxtyv, atomscore, setcount, which_combination)

      use db2type
      use optionstype
      use gridstype
      use phimaptype
      use gisttype
      use solvmaptype
      use vdwmaptype
      use phiscore
      use rec_des_score
      use chemscore
      use solvation_score
      use atomscoretype

      ! XXX: GFortran needs this first
      integer, intent(in) :: maxtyv !how many vdw atom types there are
      type(options), intent(in) :: options0 !useful options
      type(db2), intent(in) :: db2lig !ligand information
      type(flexgrids), intent(in) :: grids0 !precomp values
      type(vdwmap), intent(in) :: vdwmap0
      type(phimap), intent(in) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(in) :: gist0, rec_des_r0
      type(solvmap), intent(in) :: solvmap0
      real, dimension(maxtyv), intent(in) :: sra, srb !square roots of vdw parameters per vdwtype
      type(atomscoret), intent(inout) :: atomscore
      integer, intent(in) :: setcount !which set to score  
      integer :: which_combination !which grids to use

      integer atomstart, atomend, atomnow, atomindex !temp/loop variables
      real potent(db2lig%total_atoms) !the electrostatic potential at each atom, 
            !set in phipdb and used in charge but then discarded, so only needed
            !here, doesn't need to be global
      real rec_d_energy(db2lig%total_atoms)
      real recdes_energy(db2lig%total_atoms) !receptor desolv energy, similar to potent
      real ligdes_energy(db2lig%total_atoms) !ligand desolv potential
      real :: vdwasum, vdwbsum, eleceel, asum, bsum, electemp
      real :: solva, solvp, solvd, solvatemp, solvptemp, solvdtemp
      real :: rdval, rdtemp
c      real :: gist !hydrophobic effect
      real :: hscore, hscoretemp !hydrophobic effect
      real :: rec_energy !sum receptor energies
      integer :: numscored !keeps track of how many atom positions are evaluated
      integer :: which_grid, count_grid !which grids to use
      integer :: conf, confcount !which conf currently being examined

      do confcount = 1, db2lig%set_conf_len(setcount)
        conf = db2lig%set_conf(setcount, confcount)
        atomstart = db2lig%conf_coord(1, conf) !atom index of start
        atomend = db2lig%conf_coord(2, conf) !same but end (2)
        do atomnow = atomstart, atomend
          atomindex = db2lig%coord_index(2, atomnow)
          !reset these temporary variables
          vdwasum = 0.0
          vdwbsum = 0.0
          eleceel = 0.0
          !gist = 0.0 ! we will need to add a gist term here in future. 
          solva = 0.0
          solvp = 0.0
c          solvd = 0.0
          rdval = 0.0
          hscore = 0.0
          rec_energy = 0.0
          !for each grid in case we are doing flexible docking
          do count_grid = 1, grids0%group_len
            which_grid = grids0%all_combos(which_combination, 
     &          count_grid)
            rec_energy = rec_energy + options0%rec_energy(which_grid)
            !compute the score for one atom by calling with atomnow, atomnow
            ! instead of atomstart, atomend
            call calc_chemscore(atomnow, atomnow,
     &          db2lig%coord_index, db2lig%transfm_coords, 
     &          vdwmap0%offset, vdwmap0%invgrid, vdwmap0%grdpts,
     &          sra, srb, db2lig%atom_vdwtype,
     &          grids0%griddata(which_grid)%abval_precomp, vdwmap0%ngrd, 
     &          db2lig%total_atoms, db2lig%total_coords, maxtyv,
     &          asum, bsum, numscored)
            vdwasum = vdwasum + asum
            vdwbsum = vdwbsum + bsum
            if (options0%ldesolv .le. 2) then !full, none or GB partial
              call calc_solvation_score(atomnow, atomnow,
     &            db2lig%atom_vdwtype, db2lig%atom_polsolv, 
     &            db2lig%atom_apolsolv, db2lig%coord_index,
     &            db2lig%transfm_coords, solvmap0%sperang, 
     &            solvmap0%smax, solvmap0%scadif,
     &            grids0%griddata(which_grid)%solgrd_precomp, 
     &            grids0%griddata(which_grid)%hsolgrd_precomp, 
     &            options0%hsolflag, options0%ldesolv,
     &            db2lig%total_atoms, db2lig%total_coords, 
     &            solvptemp, solvatemp)
            else !3 = PB
              call phipdb(atomnow, atomnow, db2lig%coord_index, 
     &            db2lig%transfm_coords,
     &            ligdes0%oldmid, ligdes0%phiscale, ligdes0%goff, 
     &            grids0%griddata(which_grid)%ligdes_precomp, 
     &            db2lig%total_atoms,
     &            db2lig%total_coords, ligdes0%nsize, ligdes_energy) 
              solvatemp = 0.0
              solvptemp = charge_squared(atomnow, atomnow, 
     &            db2lig%coord_index, db2lig%atom_charge, 
     &            ligdes_energy, db2lig%total_atoms,
     &            db2lig%total_coords)
            endif
            solva = solva + solvatemp
            solvp = solvp + solvptemp
            electemp = 0.0
            if (allocated(grids0%griddata(which_grid)%phimap_precomp))
     &          then !check for implicit 0 grids
              call phipdb(atomnow, atomnow, db2lig%coord_index, 
     &            db2lig%transfm_coords,
     &            phimap0%oldmid, phimap0%phiscale, phimap0%goff, 
     &            grids0%griddata(which_grid)%phimap_precomp, 
     &            db2lig%total_atoms,
     &            db2lig%total_coords, phimap0%nsize, potent) !sets potent() array for use next
              electemp = charge(atomnow, atomnow, db2lig%coord_index,
     &            db2lig%atom_charge, potent, db2lig%total_atoms, 
     &            db2lig%total_coords)
            endif
            eleceel = eleceel + electemp
            rdtemp = 0.0
            if (allocated(grids0%griddata(which_grid)%rec_des_precomp))
     &          then
              call dx_pb(db2lig%atom_vdwtype, atomnow, atomnow, 
     &            db2lig%coord_index,
     &            db2lig%transfm_coords, rec_des_r0%orgin,
     &            rec_des_r0%gistspace, rec_des_r0%xnsize,
     &            rec_des_r0%ynsize,rec_des_r0%znsize,
     &            grids0%griddata(which_grid)%rec_des_precomp,
     &            grids0%griddata(which_grid)%hrec_des_precomp,
     &            db2lig%total_atoms, db2lig%total_coords, rec_d_energy)
              rdtemp = convert_rdsum(atomnow, atomnow, 
     &            db2lig%coord_index, rec_d_energy,
     &            db2lig%atom_vdwtype, db2lig%total_atoms, 
     &            db2lig%total_coords, rec_des_r0%gistspace)
            endif
            rdval = rdval + rdtemp
            if (options0%receptor_desolv) then !do receptor desolvation calculation
              if (allocated(grids0%griddata(which_grid)%recdes_precomp))
     &            then !implicit 0 checking
                call phipdb(atomnow, atomnow, db2lig%coord_index,
     &              db2lig%transfm_coords,
     &              recdes0%oldmid, recdes0%phiscale, recdes0%goff,
     &              grids0%griddata(which_grid)%recdes_precomp,
     &              db2lig%total_atoms,
     &              db2lig%total_coords, recdes0%nsize, recdes_energy) !set to kT
                solvdtemp = convert_sum(atomnow, atomnow, 
     &              db2lig%coord_index, recdes_energy, 
     &              db2lig%atom_vdwtype, db2lig%atom_heavy_valence,
     &              db2lig%total_atoms, db2lig%total_coords)
              endif
              solvd = solvd + solvdtemp
            !also need to do hydrophobic effect here
            endif
!            write(6,*) "xxx", atomindex, atomnow, which_grid, vdwasum,
!     &          vdwbsum, solva, solvp, eleceel, rec_energy
!            call doflush(6)
          enddo
          !now that each grid is done, sum them all and scale properly
          atomscore%atom_vas(atomindex) = -vdwbsum * options0%vscale
          atomscore%atom_vrs(atomindex) = vdwasum * options0%vscale
          atomscore%atom_es(atomindex) = eleceel * options0%escale
          atomscore%atom_as(atomindex) = solva * options0%solvscale
          atomscore%atom_ps(atomindex) = solvp * options0%solvscale
c          atomscore%atom_ds(atomindex) = solvd * options0%rdscale
          atomscore%atom_rd(atomindex) = rdval * options0%rec_des_scale
          atomscore%atom_hs(atomindex) = hscore * options0%rhscale
          !compute overall sum, does not include receptor energy 
          ! or internal energy since they aren't atomistic
          atomscore%atom_score(atomindex) = 
     & atomscore%atom_vas(atomindex) + atomscore%atom_vrs(atomindex) 
     & + atomscore%atom_es(atomindex) + atomscore%atom_as(atomindex)
     & + atomscore%atom_ps(atomindex) + atomscore%atom_rd(atomindex)
     & + atomscore%atom_hs(atomindex)
c     & + atomscore%atom_ds(atomindex) + atomscore%atom_hs(atomindex)
        enddo
      enddo

      return
      end subroutine recalc_score_byatom
 
      end module rescore_byatom
