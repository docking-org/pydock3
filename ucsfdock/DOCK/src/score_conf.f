c     this is an actual scoring function. input is a group number. output
c     is a bunch of different kinds of scores. it assumes atoms have already
c     been transformed and put into xatm appropriately
      module score_conf

      implicit none
      contains

      subroutine calc_score_conf(conf, options0, db2lig, grids0,
     &    vdwmap0, phimap0, recdes0, ligdes0, gist0, 
     &    solvmap0, rec_des_r0, sra, srb, which_grid,
     &    maxtyv,
     &    vdwasum, vdwbsum, eleceel, gistval, solva, solvp, 
     &    solvd, rdval, hescore, numscored,
     &    gistindexlist, lengistil,
     &    gistindexlcon, lengistilc)

      use db2type
      use optionstype
      use gridstype
      use phimaptype
      use gisttype
      use solvmaptype
      use vdwmaptype
      use rec_des_score
      use phiscore
      use gistscore
      use chemscore
      use solvation_score

c XXX: GFortran needs this to come first
c sizing constants from max.h
      integer, intent(in) :: maxtyv !how many vdw atom types there are
      integer, intent(in) :: conf !the conf that the caller would like scored
      type(options), intent(in) :: options0 !useful options
      type(db2), intent(in) :: db2lig !ligand information
      type(flexgrids), intent(inout) :: grids0 !precomp values
      type(vdwmap), intent(in) :: vdwmap0
      type(phimap), intent(in) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(in) :: gist0, rec_des_r0
      type(solvmap), intent(in) :: solvmap0
      real, dimension(maxtyv), intent(in) :: sra, srb !square roots of vdw parameters per vdwtype
      integer, intent(in) :: which_grid !which grid to use
c return variables
      real, intent(inout) :: vdwasum, vdwbsum, eleceel, gistval, rdval
      real, intent(inout) :: solva, solvp, solvd, hescore !scores by reference
      integer, intent(inout) :: numscored !keeps track of how many atom positions are evaluated

      integer, intent(inout) :: gistindexlist(5000)
      integer, intent(inout) :: lengistil
      integer, intent(inout) :: gistindexlcon(5000)
      integer, intent(inout) :: lengistilc
  
      integer atomstart, atomend !temp/loop variables
      real potent(db2lig%total_atoms) !the electrostatic potential at each atom, 
            !set in phipdb and used in charge but then discarded, so only needed
            !here, doesn't need to be global
      real gistvec(db2lig%total_atoms)
      real rec_d_energy(db2lig%total_atoms)!receptor desolvation energy
      real recdes_energy(db2lig%total_atoms) !the receptor desolvation energy
            !of each atom, in kT (needs converted still)
      real ligdes_energy(db2lig%total_atoms) !ligand desolvation energy, kT


      atomstart = db2lig%conf_coord(1, conf) !atom index of start
      atomend = db2lig%conf_coord(2, conf) !same but end (2)
      call calc_chemscore(atomstart, atomend,
     &    db2lig%coord_index, db2lig%transfm_coords, vdwmap0%offset, 
     &    vdwmap0%invgrid, vdwmap0%grdpts,
     &    sra, srb, db2lig%atom_vdwtype,
     &    grids0%griddata(which_grid)%abval_precomp, vdwmap0%ngrd, 
     &    db2lig%total_atoms, db2lig%total_coords, maxtyv,
     &    vdwasum, vdwbsum, numscored)
      !write(6,*) "calc_score:I AM HERE(1)" 
c     sets vdwasum and vdwbsum (vdw scoring), increments numscored, reads others
      if ((vdwasum - vdwbsum) .le. options0%bmpvdw) then !otherwise bump this group
        if (options0%ldesolv .le. 2) then !none, full or old GB partial
          call calc_solvation_score(atomstart, atomend,
     &        db2lig%atom_vdwtype, db2lig%atom_polsolv, 
     &        db2lig%atom_apolsolv, db2lig%coord_index,
     &        db2lig%transfm_coords, solvmap0%sperang, 
     &        solvmap0%smax, solvmap0%scadif,
     &        grids0%griddata(which_grid)%solgrd_precomp, 
     &        grids0%griddata(which_grid)%hsolgrd_precomp, 
     &        options0%hsolflag, options0%ldesolv,
     &       db2lig%total_atoms, db2lig%total_coords, solvp, solva)
          !write(6,*) "calc_score:I AM HERE(2)" 
c         solvation_score sets solvp and solva
        else !ldesolv = 3 = PB ligand desolvation
          call phipdb(atomstart, atomend, db2lig%coord_index, 
     &        db2lig%transfm_coords,
     &        ligdes0%oldmid, ligdes0%phiscale, ligdes0%goff, 
     &        grids0%griddata(which_grid)%ligdes_precomp, 
     &        db2lig%total_atoms,
     &        db2lig%total_coords, ligdes0%nsize, ligdes_energy) !set to kT
          solva = 0.0
          solvp = charge_squared(atomstart, atomend, 
     &        db2lig%coord_index,
     &        db2lig%atom_charge, ligdes_energy, db2lig%total_atoms, 
     &        db2lig%total_coords)
          !write(6,*) "calc_score:I AM HERE(3)" 
        endif
        eleceel = 0.0
        if (allocated(grids0%griddata(which_grid)%phimap_precomp))
     &      then !this check is to not read data when implicit0 grids are used
          call phipdb(atomstart, atomend, db2lig%coord_index, 
     &        db2lig%transfm_coords,
     &        phimap0%oldmid, phimap0%phiscale, phimap0%goff, 
     &        grids0%griddata(which_grid)%phimap_precomp, 
     &        db2lig%total_atoms,
     &        db2lig%total_coords, phimap0%nsize, potent) !sets potent() array for use next
          eleceel = charge(atomstart, atomend, db2lig%coord_index,
     &        db2lig%atom_charge, potent, db2lig%total_atoms, 
     &        db2lig%total_coords)
          !write(6,*) "calc_score:I AM HERE(4)" 
        endif
        rdval = 0.0
        if (options0%rec_d_flag(which_grid) .and. 
     &      options0%hrec_d_flag(which_grid)) then
c        if (allocated(grids0%griddata(which_grid)%rec_des_precomp))
c     &      then
          call dx_pb(db2lig%atom_vdwtype, atomstart, atomend, 
     &        db2lig%coord_index,
     &        db2lig%transfm_coords, rec_des_r0%orgin, 
     &        rec_des_r0%gistspace, rec_des_r0%xnsize, 
     &        rec_des_r0%ynsize,rec_des_r0%znsize, 
     &        grids0%griddata(which_grid)%rec_des_precomp, 
     &        grids0%griddata(which_grid)%hrec_des_precomp,
     &        db2lig%total_atoms, db2lig%total_coords, rec_d_energy)  
          rdval = convert_rdsum(atomstart, atomend, db2lig%coord_index,
     &        rec_d_energy, db2lig%atom_vdwtype,
     &        db2lig%total_atoms, db2lig%total_coords,
     &        rec_des_r0%gistspace)
        endif
        gistval = 0.0
        if (options0%gistflag(which_grid)) then ! if gist grid is used call the
                                    ! function  other wise set gistval
                                    ! = 0.0
          !if (.not.options0%gist_aprox) then ! this is the most correte but is slow. 
          !write(6,*) "calc_score:gist_aprox = ", options0%gist_aprox 
          if (options0%gist_aprox(which_grid) == 0) then ! this is the most correte but is slow. 
          !write(*,*) "I AM HERE"
          !write(*,*) gist0%orgin, gist0%gistspace, gist0%xnsize
          call gistscorefuc(atomstart, atomend, db2lig%coord_index,
     &        db2lig%transfm_coords,
     &        gist0%orgin, gist0%gistspace,
     &        options0%gistvolume, options0%gistH,
     &        gist0%xnsize,
     &        gist0%ynsize, gist0%znsize,
     &        sra, srb, db2lig%atom_vdwtype,
     &        grids0%griddata(which_grid)%gistgrid,
     &        grids0%griddata(which_grid)%gisttagvoxelarray,
     &        grids0%griddata(which_grid)%gistneighborvoxelarray,
     &        db2lig%total_atoms,
     &        db2lig%total_coords, maxtyv, gistvec, gistval,
     &        gistindexlist, lengistil,
     &        gistindexlcon, lengistilc) ! gistval store the score
          else if(options0%gist_aprox(which_grid) == 1) then ! use the aproxamtion 1point. 
          call gistscorefuc_fast_1(atomstart, atomend, 
     &        db2lig%coord_index,
     &        db2lig%transfm_coords,
     &        gist0%orgin, gist0%gistspace, 
     &        options0%gistvolume,
     &        gist0%xnsize,
     &        gist0%ynsize, gist0%znsize,
     &        sra, srb, db2lig%atom_vdwtype,
     &        grids0%griddata(which_grid)%gistgrid,
     &        db2lig%total_atoms,
     &        db2lig%total_coords, maxtyv, gistvec, gistval) ! gistval store the score

          else if(options0%gist_aprox(which_grid) == 2) then ! use the aproxamtion 8 points.
          call gistscorefuc_fast_8(atomstart, atomend, 
     &        db2lig%coord_index,
     &        db2lig%transfm_coords,
     &        gist0%orgin, gist0%gistspace, 
     &        options0%gistvolume,
     &        gist0%xnsize,
     &        gist0%ynsize, gist0%znsize,
     &        sra, srb, db2lig%atom_vdwtype,
     &        grids0%griddata(which_grid)%gistgrid,
     &        db2lig%total_atoms,
     &        db2lig%total_coords, maxtyv, gistvec, gistval) ! gistval store the score
c          else c we should never get here
c              write(OUTDOCK,*) "Error: "
          endif
          !write(6,*) "calc_score:I AM HERE(5)" 
        endif
c       write(*,*) "gistflag", which_grid,options0%gistflag(which_grid),
c    & options0%gist_aprox(which_grid),gistval
        solvd = 0.0 !receptor desolvation is 0. if not requested
        if (options0%receptor_desolv) then !do receptor desolvation calculation
          if (allocated(grids0%griddata(which_grid)%recdes_precomp)) 
     &        then !the above check is for implicit0 grids
            call phipdb(atomstart, atomend, db2lig%coord_index, 
     &        db2lig%transfm_coords,
     &        recdes0%oldmid, recdes0%phiscale, recdes0%goff, 
     &        grids0%griddata(which_grid)%recdes_precomp, 
     &        db2lig%total_atoms,
     &        db2lig%total_coords, recdes0%nsize, recdes_energy) !set to kT
            solvd = convert_sum(atomstart, atomend, db2lig%coord_index,
     &        recdes_energy, db2lig%atom_vdwtype, 
     &        db2lig%atom_heavy_valence, db2lig%total_atoms, 
     &        db2lig%total_coords)
          !write(6,*) "calc_score:I AM HERE(6)" 
          endif
        endif
        hescore = 0.0 !receptor hydrophobic effecti s 0. if not requested
      endif
      return
      end subroutine calc_score_conf
 
      end module score_conf

