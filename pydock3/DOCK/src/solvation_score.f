      module solvation_score

      implicit none
      contains

c----------------------------------------------------------------------
c  Prorate desolvation using vol. occluded by receptor
c-----------------------------------------------------------------------
      subroutine calc_solvation_score(atomstart, atomend,
     &    atom_vdwtype, atom_polsolv, atom_apolsolv, coord_index, 
     &    transfm_coords, sperang, smax, scadif, 
     &    solgrd_precomp, hsolgrd_precomp, hsolflag, ldesolv,
     &    maxatm, maxatmpos,
     &    solvp, solva)
      
c XXX: GFortran needs this to come first (TS)
c maximum values, from max.h
      integer, intent(in) :: maxatm !how many atoms are in the ligand
      integer, intent(in) :: maxatmpos !max number of positions all atoms can take
c input values
      integer, intent(in) :: atomstart, atomend !which atoms to calculate
      integer, intent(in) :: atom_vdwtype(maxatm) !atom dock vdw type
      real, intent(in) :: atom_polsolv(maxatm) !atom polar solvation
      real, intent(in) :: atom_apolsolv(maxatm) !atom apolar solvation
      integer, intent(in) :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real, intent(in) :: transfm_coords(3, maxatm) !ligand atoms coordinates moved into bs
      integer, intent(in) :: sperang ! grid points/angstrom for desolvation grid.
      real, intent(in) :: smax(3) !maximum real space xyz values
      integer, intent(in) :: scadif(3) !dimensions of the solvation grid
      real, intent(in) :: solgrd_precomp(8, 
     &    0:scadif(1), 0:scadif(2), 0:scadif(3))
      real, intent(in) :: hsolgrd_precomp(8, 
     &    0:scadif(1), 0:scadif(2), 0:scadif(3))
      logical, intent(in) :: hsolflag !whether or not there is a separate hydrogen solvmap
      integer, intent(in) :: ldesolv !0 = none, 1 = full, 2 = partial
c return values
      real, intent(inout) :: solvp, solva !polar and apolar desolvation

      integer atomindex !actual atom index store here
      integer count, lx, ly, lz
      integer nx, ny, nz
      real    slx, sly, slz, xgr, ygr, zgr
      real    a1, a2, a3, a4, a5, a6, a7, a8
      real    volocl ! % desolvation based on sol. volume occluded by receptor.

c     use precomputed values, removed hierarchy knowledge
c     REMOVED OUTSIDE GRID CHECKS

      volocl = 0.
      if (ldesolv .eq. 0) then
        volocl = 0. !no ligand desolvation at all
      else if (ldesolv .eq. 1) then
        volocl = 1. !full ligand desolvation
      endif !otherwise ldesolv = 2, so use partial
      solvp = 0.0
      solva = 0.0
      do count = atomstart, atomend !for every atom
        atomindex = coord_index(2, count) !gets actual atom non-consecutive #
        if (ldesolv .eq. 2) then !only if doing partial ligand desolvation
          slx = smax(1) - transfm_coords(1, atomindex) * sperang
          sly = smax(2) - transfm_coords(2, atomindex) * sperang
          slz = smax(3) - transfm_coords(3, atomindex) * sperang
          nx = int(slx)
          ny = int(sly)
          nz = int(slz)
c         calculate partial cube coordinates of point
          xgr = slx - float(nx)
          ygr = sly - float(ny)
          zgr = slz - float(nz)
c         calculate coefficients of trilinear function
          if (hsolflag .and. ((atom_vdwtype(atomindex) .eq. 6) .or.
     &        (atom_vdwtype(atomindex) .eq. 7))) then !hydrogen atom 
            a8 = hsolgrd_precomp(8, nx, ny, nz)
            a7 = hsolgrd_precomp(7, nx, ny, nz)
            a6 = hsolgrd_precomp(6, nx, ny, nz)
            a5 = hsolgrd_precomp(5, nx, ny, nz)
            a4 = hsolgrd_precomp(4, nx, ny, nz)
            a3 = hsolgrd_precomp(3, nx, ny, nz)
            a2 = hsolgrd_precomp(2, nx, ny, nz)
            a1 = hsolgrd_precomp(1, nx, ny, nz)
          else !heavy atom
            a8 = solgrd_precomp(8, nx, ny, nz)
            a7 = solgrd_precomp(7, nx, ny, nz)
            a6 = solgrd_precomp(6, nx, ny, nz)
            a5 = solgrd_precomp(5, nx, ny, nz)
            a4 = solgrd_precomp(4, nx, ny, nz)
            a3 = solgrd_precomp(3, nx, ny, nz)
            a2 = solgrd_precomp(2, nx, ny, nz)
            a1 = solgrd_precomp(1, nx, ny, nz)
          endif
c         interpolate values read from solvent grid
          volocl = a1*xgr*ygr*zgr + a2*xgr*ygr + a3*xgr*zgr
     &                   + a4*ygr*zgr + a5*xgr + a6*ygr + a7*zgr + a8
          if (volocl .ge. 1.0 ) then !max is 1, means completely desolvated
            write(*,*) "Warning: volocl (desolvation grid) > 1.0",volocl
            volocl = 1.0
          elseif ( volocl .lt. 0.0) then !max is 1, means completely desolvated
            write(*,*) "Warning: volocl (desolvation grid) < 0.0",volocl
            volocl = 0.0
          endif
        endif
        !write(*,*) "volocl (desolvation grid) = ", volocl
        solvp = solvp - (atom_polsolv(atomindex) * volocl)
        solva = solva - (atom_apolsolv(atomindex) * volocl)
      enddo

      return
      end subroutine calc_solvation_score

      end module solvation_score
