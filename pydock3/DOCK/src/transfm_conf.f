c translates/rotates one conf into position
      module transfm_conf
 
      implicit none
      contains      

      subroutine run_transfm_conf(conf, curmatch, db2lig,
     &    match, maxor)

      use matchtype
      use db2type

      integer, intent(in) :: conf !input conformation
      integer, intent(in) :: curmatch !which match of trans/rots to do
      type(db2), intent(inout) :: db2lig !ligand information
      type(matcht), intent(in) :: match
c max parameters, see max.h
      integer, intent(in) :: maxor !the size of many arrays that store information about
                    !orientations, is checked several times

      integer tempcoord, tempatom !used to find the right place to put things
      integer count !from 1 to 3 for x to z

      !put coords from this conf into transfm_coords
      do tempcoord = db2lig%conf_coord(1, conf), 
     &     db2lig%conf_coord(2, conf)
        tempatom = db2lig%coord_index(2, tempcoord)
        do count = 1, 3
          db2lig%transfm_coords(count, tempatom) = 
     &        match%comr(count, curmatch) + 
     &        match%rot(count, 1, curmatch) * 
     &         (db2lig%coords(1, tempcoord) - match%coml(1, curmatch)) +
     &        match%rot(count, 2, curmatch) * 
     &         (db2lig%coords(2, tempcoord) - match%coml(2, curmatch)) +
     &        match%rot(count, 3, curmatch) * 
     &         (db2lig%coords(3, tempcoord) - match%coml(3, curmatch))
        enddo
      enddo
      !all coords for this conf have been transformed and put in transfm_coords
      return
      end subroutine run_transfm_conf

      end module transfm_conf
