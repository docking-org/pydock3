      module check_conf_og


      implicit none
      contains

c checks each atom to make sure it is inside the grids. returns .true. 
c immediately if any atom is outside any grid. .false. if all pass
! XXX function check_conf_og(conf, db2lig,
      function atom_out_of_grids(conf, db2lig,
     &    allminlimits, allmaxlimits) result(conf_out)

      use db2type

      integer, intent(in) :: conf
      type(db2), intent(in) :: db2lig !ligand information
      real, intent(in) :: allminlimits(3), allmaxlimits(3) !grid limits in xyz

      ! XXX: Doesn't play nice with GFortran (TS)
      !logical, intent(out) :: conf_out
      logical :: conf_out

      integer tempcoord, tempatom !used to find the right place to put things
      integer count !from 1 to 3 for x to z
      real xyz(3) !constructed on the fly for each atom
    
      do tempcoord = db2lig%conf_coord(1, conf), 
     &    db2lig%conf_coord(2, conf)
        tempatom = db2lig%coord_index(2, tempcoord)
        do count = 1, 3
          xyz(count) = db2lig%transfm_coords(count, tempatom)
        enddo
        if (outside_grids(xyz, allminlimits, allmaxlimits)) then !atom outside
          conf_out = .true. !short circuit out of loop
          return
        endif
      enddo
      conf_out = .false. !all atoms checked and are okay
      return
      end function atom_out_of_grids

c new function to check 1 atom position to see if it outside any grids
c returns .true. if outside any grid and .false. if inside all grids
      function outside_grids(xyz, 
     &    allminlimits, allmaxlimits) result(atom_out)

c input variables
      real, dimension(3), intent(in) :: xyz !xyz coordinates to be checked
      real, dimension(3), intent(in) :: allminlimits !limits in xyz space
      real, dimension(3), intent(in) :: allmaxlimits !limits in xyz space

c output variable
      ! XXX: Doesn't play nice with GFortran
      !logical, intent(out) :: atom_out
      logical :: atom_out

c temporary variables
      integer coord !counts to 3

c     this set of code checks the overall boundary
      do coord = 1,3
        if (xyz(coord) .lt. allminlimits(coord)) then
          atom_out = .true.
          return
        endif
        if (xyz(coord) .gt. allmaxlimits(coord)) then
          atom_out = .true.
          return
        endif
      enddo
c     if we got here, everything is fine
      atom_out = .false.
      return
      end function outside_grids

      end module check_conf_og

