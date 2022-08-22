!module for phimap scoring, contains phipdb & charge (also desolvation routines)
      module rec_des_score

      implicit none
      contains     

************************************************************************
c   This was originally QNIFFT/PB-based receptor desolvation (didn't
c   work too well, did it? Use 3 dielectrics, future DOCKer. For those
c   about to DOCK, I salute you.)
c   written 10/2016
c   modified to work with precomputed GIST grids
c   Written by Reed M. Stein. 
************************************************************************
      subroutine dx_pb(atom_vdwtype, atomstart, atomend, 
     &    coord_index, transfm_coords, orgin, gistspace,
     &    xnsize, ynsize, znsize, rec_des_precomp, hrec_des_precomp,
     &    maxatm, maxatmpos, score)

c XXX: These need to come first or Gfortran complains (TS)
c max constants, see max.h
      integer, intent(in) :: maxatm !how many atoms max can be in a ligand
      integer, intent(in) :: maxatmpos !how many positions all ligand atoms can take
      integer, intent(in) :: xnsize !grid size for x
      integer, intent(in) :: ynsize !grid size for y
      integer, intent(in) :: znsize !grid size for z
      integer, intent(in) :: atom_vdwtype(maxatm)
c input variables
      integer, intent(in) :: atomstart, atomend !where the atoms and other things are
      integer, intent(in) :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real, intent(in) :: transfm_coords(3, maxatm) !ligand atoms coordinates moved into bs
      real, intent(in) :: orgin(3) !the front,lower left corner of box
      real, intent(in) :: gistspace !angstroms per grid point on dx grid
      real, intent(in) :: rec_des_precomp(8, xnsize, ynsize, znsize) !actual energy values
      real, intent(in) :: hrec_des_precomp(8, xnsize, ynsize, znsize) !actual energy values
c return variable
      real, intent(inout) :: score(maxatm) !score at each atom

      integer atomindex !actual position used
      integer count !atom position counter
      integer nx, ny, nz
      real xn(3)
      real a1, a2, a3, a4, a5, a6, a7, a8
      real xgr, ygr, zgr
      !real outside_grid_val        

      !write(6,*) "calc_score:I AM IN REC_DES_SCORE" 
      !call doflush(6)
      
      !outside_grid_val = -1000000.0
      !write(6,*) "I AM IN REC_DES_SCORE"
      !call doflush(6)

c     loop over all atoms
      do count = atomstart, atomend
        atomindex = coord_index(2, count) !gets actual atom non-consecutive #
c       scale atoms to grid space
        xn(1) = (transfm_coords(1, atomindex) - 
     &       orgin(1)) / gistspace + 1
        xn(2) = (transfm_coords(2, atomindex) -
     &       orgin(2)) / gistspace + 1
        xn(3) = (transfm_coords(3, atomindex) - 
     &       orgin(3)) / gistspace + 1

c       find lower left bottom grid point
        nx = int(xn(1))
        ny = int(xn(2))
        nz = int(xn(3))
c       $mem prefetch   phimap_precomp(1,nx,ny,nz)


c       calculate fractional cube coordinates of point
        xgr = xn(1) - float(nx)
        ygr = xn(2) - float(ny)
        zgr = xn(3) - float(nz)

        !write(6,*) "defining rec_des_precomp"
        !call doflush(6)
        if ((atom_vdwtype(atomindex) .eq. 6) .or.
     &        (atom_vdwtype(atomindex) .eq. 7)) then 
               a8 = hrec_des_precomp(8, nx, ny, nz)
               a7 = hrec_des_precomp(7, nx, ny, nz)
               a6 = hrec_des_precomp(6, nx, ny, nz)
               a5 = hrec_des_precomp(5, nx, ny, nz)
               a4 = hrec_des_precomp(4, nx, ny, nz)
               a3 = hrec_des_precomp(3, nx, ny, nz)
               a2 = hrec_des_precomp(2, nx, ny, nz)
               a1 = hrec_des_precomp(1, nx, ny, nz) 
        else ! heavy atom
               a8 = rec_des_precomp(8, nx, ny, nz)
               a7 = rec_des_precomp(7, nx, ny, nz)
               a6 = rec_des_precomp(6, nx, ny, nz)
               a5 = rec_des_precomp(5, nx, ny, nz)
               a4 = rec_des_precomp(4, nx, ny, nz)
               a3 = rec_des_precomp(3, nx, ny, nz)
               a2 = rec_des_precomp(2, nx, ny, nz)
               a1 = rec_des_precomp(1, nx, ny, nz)
        endif
             
c       calculate potential - speed optimized by {MC}
c       pull coefficients of trilinear function from precomp

       !if (nx .lt. 1.0 .or. nx + 1 .gt. xnsize) then
       !    score(atomindex) = outside_grid_val 
       !    !write(6, *) "x is too big. nx is", nx
       !else if (ny .lt. 1.0 .or. ny + 1 .gt. ynsize) then
       !    score(atomindex) = outside_grid_val 
       !    !write(6, *) "y is too big. ny is", ny
       !else if (nz .lt. 1.0 .or. nz + 1 .gt. znsize) then
       !    score(atomindex) = outside_grid_val 
       !    !write(6, *) "z is too big. nz is", nz
       !else 
c          determine value of potential
           score(atomindex) = a1*xgr*ygr*zgr + a2*xgr*ygr + a3*xgr*zgr
     &        + a4*ygr*zgr + a5*xgr + a6*ygr + a7*zgr + a8
       !endif
      enddo

      return
      end subroutine dx_pb

!function basically converts kT units to kcal/mol, used for receptor desolvation
!uses atom_heavy_valence to reduce the effect by 1/3 for each heavy atom
!connected. atom_heavy_valence should be 0..3 inclusive, above 3 makes no sense
!and neither does below 0. no need to reference options0%recdes_valence_penalty
!since if it isn't set then everything will be 0.
      function convert_rdsum(atomstart, atomend, coord_index,  
     &    kT_energy, atom_vdwtype, maxatm, maxatmpos, gistspace)
     &    result(kcal_energy_sum)

c XXX: These need to come first or Gfortran complains (TS)
c max constants
      integer, intent(in) :: maxatm !how many atoms max can be in a ligand
      integer, intent(in) :: maxatmpos !how many positions all ligand atoms can take

c input variables
      integer, intent(in) :: atomstart, atomend
      integer, intent(in) :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real, intent(in) :: kT_energy(maxatm) !energy in kT at each atom
      real, intent(in) :: gistspace ! dimension of voxel (Ang)
      integer, intent(in) :: atom_vdwtype(maxatm) !atom dock vdw type
c output values
      !real, intent(out) :: kcal_energy_sum ! GFortran chokes on this (TS)
      real :: kcal_energy_sum
      real :: volvox ! volume of voxel ( ang^3) 

      integer atomindex, count !used to find correct atoms
      real this_energy !used to temporarily store energy

      volvox = gistspace * gistspace * gistspace
      kcal_energy_sum = 0.0 !the returned value
      do count = atomstart, atomend !for every atom
        atomindex = coord_index(2, count) !gets actual atom non-consecutive #
        !ignore hydrogen atoms for receptor desolvation, they don't have enough
        !bulk to desolvate, only the heavy atoms do
        this_energy = kT_energy(atomindex) 
!          write (6,*) 'xxxpre', this_energy !debugging
        this_energy = this_energy * (-1.0) * volvox
!          write (6,*) 'xxxpost', this_energy, 
!     &        atom_heavy_valence(atomindex)
        kcal_energy_sum = kcal_energy_sum + this_energy
!        endif
      enddo
!      kcal_energy_sum = kcal_energy_sum * 0.5924 !const to convert to kcal/mole
      return
      end function convert_rdsum

      end module rec_des_score
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
