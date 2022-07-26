c----------------------------------------------------------------------
c
c  --calculates force field scores.                      ECMeng  11/91
c  --removed eescore, just calculates van der waal scores
c  --removed all hierarchy stuff from scoring
c  --removed global variables. rgc 12/2011.
c------------------------------------------------------------------------
      module chemscore

      ! XXX: GFortran chokes on this when handling abval_precomp types
      !implicit none
      contains

      subroutine calc_chemscore(atomstart, atomend, 
     &    coord_index, transfm_coords, offset, invgrid, grdpts,
     &    sra, srb, atom_vdwtype,
     &    abval_precomp, ngrd, 
     &    maxatm, maxatmpos, maxtyv, 
     &    vdwasum, vdwbsum, numscored)

      use vdwconsts

c XXX: GFortran needs these first (TS)
c max constants from max.h
      integer, intent(in) :: maxatm !how many atoms max can be in a ligand
      integer, intent(in) :: maxatmpos !how many positions all ligand atoms can take
      integer, intent(in) :: maxtyv !how many vdw atom types there are
c input variables
      integer, intent(in) :: atomstart, atomend !atoms to score
      integer, intent(in) :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real, intent(in) :: transfm_coords(3, maxatm) !ligand atoms coordinates moved into bs
      real, intent(in) :: offset(3) !grid xmin, ymin, zmin in angstroms
      real, intent(in) :: invgrid !angstroms per grid
      integer, intent(in) :: grdpts(3) !number of grid points in each dimenion (x,y,z)
      real, intent(in) :: sra(maxtyv), srb(maxtyv) !square roots of vdw parameters per vdwtype
      integer, intent(in) :: atom_vdwtype(maxatm) !the vdw type of each ligand atom
      real, intent(in) :: abval_precomp(2, 8, ngrd) !actual grid values of vdw scores,
          ! precomputed a la the Carchia optimizations
      integer, intent(in) :: ngrd !how many grid points there are, dynamically changes each run
c return values
      real, intent(inout) :: vdwasum, vdwbsum !where to put scores
      integer, intent(inout) :: numscored !keeps track of how many atom positions are evaluated

c temporary variables
      integer atomindex !actual location of atom
      integer count, array_index
      real a1, a2, a3, a4, a5, a6, a7, a8
      integer nx, ny, nz
      real xn(3), xgr, ygr, zgr
c       xn - grid coordinates of ligand atom
c       xgr, ygr, zgr - fractional cube coordinates
      real apt, bpt
c       apt, bpt: interpolated grid values
      real avdw, bvdw
c       avdw: vdw repulsion energy of a ligand atom
c       bvdw: vdw dispersion energy of a ligand atom
c     abval_precomp declared locally so that dimensions of array can be 
c     known.  This is a dynamically allocated array and thus the common 
c     block does not know its size. except it does. so it isn't redeclared.
c      real abval_precomp(2,8,ngrd)

      vdwasum = 0.0
      vdwbsum = 0.0
      do count = atomstart, atomend !for every atom
        atomindex = coord_index(2, count) !gets actual atom non-consecutive #
c       scale atoms to grid space
c       optimization {MC}: use multiplication instead of division
        xn(1) = (transfm_coords(1, atomindex) - offset(1)) * invgrid
        xn(2) = (transfm_coords(2, atomindex) - offset(2)) * invgrid
        xn(3) = (transfm_coords(3, atomindex) - offset(3)) * invgrid

c       find lower left bottom grid point
        nx = int(xn(1))
        ny = int(xn(2))
        nz = int(xn(3))
         
c       calculate fractional cube coordinates of point
        xgr = xn(1) - float(nx)
        ygr = xn(2) - float(ny)
        zgr = xn(3) - float(nz)

c       get index into precomputed interpolation array
        array_index = grdpts(1)*grdpts(2)*nz + 
     &         grdpts(1)*ny + nx + 1

c       calculate vdw a term - speed optimized by {MC}
c       pull coefficients of trilinear function from precomp
        a1 = abval_precomp(AVAL_INDEX, 1, array_index)
        a2 = abval_precomp(AVAL_INDEX, 2, array_index)
        a3 = abval_precomp(AVAL_INDEX, 3, array_index)
        a4 = abval_precomp(AVAL_INDEX, 4, array_index)
        a5 = abval_precomp(AVAL_INDEX, 5, array_index)
        a6 = abval_precomp(AVAL_INDEX, 6, array_index)
        a7 = abval_precomp(AVAL_INDEX, 7, array_index)
        a8 = abval_precomp(AVAL_INDEX, 8, array_index)

c       determine interpolated VDW repulsion factor
        apt = a1*xgr*ygr*zgr + a2*xgr*ygr +
     &         a3*xgr*zgr + a4*ygr*zgr + a5*xgr +
     &         a6*ygr + a7*zgr + a8

c       calculate vdw b term  - speed optimized by {MC}
        a1 = abval_precomp(BVAL_INDEX, 1, array_index)
        a2 = abval_precomp(BVAL_INDEX, 2, array_index)
        a3 = abval_precomp(BVAL_INDEX, 3, array_index)
        a4 = abval_precomp(BVAL_INDEX, 4, array_index)
        a5 = abval_precomp(BVAL_INDEX, 5, array_index)
        a6 = abval_precomp(BVAL_INDEX, 6, array_index)
        a7 = abval_precomp(BVAL_INDEX, 7, array_index)
        a8 = abval_precomp(BVAL_INDEX, 8, array_index)

c       determine interpolated VDW dispersion factor
        bpt = a1*xgr*ygr*zgr + a2*xgr*ygr +
     &         a3*xgr*zgr + a4*ygr*zgr + a5*xgr +
     &         a6*ygr + a7*zgr + a8

c       compute final vdw terms
        avdw = sra(atom_vdwtype(atomindex)) * apt 
        bvdw = srb(atom_vdwtype(atomindex)) * bpt

c       sum vdw terms
        vdwasum = vdwasum + avdw
        vdwbsum = vdwbsum + bvdw
      enddo

      numscored = numscored + 1
      return
      end subroutine calc_chemscore

      end module chemscore
