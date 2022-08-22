!module for gist scoring 
      module gistscore

      implicit none
      contains     

************************************************************************
c   This modual is created to caluculate the desolvation of the 
c   Receptor using Grid Inhomogeneous Solvation Theory (GIST)
c   (Nguyen et al. J. Chem. Phys. 137 044101 2012 )
c   This subroutine was writen by Trent Balius
************************************************************************
c   old function was removed. 

************************************************************************
c   This modual is created to caluculate the desolvation of the 
c   Receptor using Grid Inhomogeneous Solvation Theory (GIST)
c   (Nguyen et al. J. Chem. Phys. 137 044101 2012 )
c   This subroutine was writen by Trent Balius with help from Teague Sterling
c   This is to be a smarter implementation.
c   To speed up the code we did the following:
c     Instead of a taged list, we created a boolen array (xnsize*ynsize*znsize)
c     to store weather a grid voxel has been seen. 
c     This means that we do not have to loop over the tage list, instead
c     we just chech using the index. this change made the code ~2/3 faster
c     We do the samething as for the taged list, but for the
c     nieborlist. this list is wipe after each atom. 
c   Additional improvements to try, not implamented: 
c     Wait to sum up all voxels until the end after all atoms have been
c     processed, not for each atom. 
c     Creat a list of voxels for each atom.  
c     loop over each atom and then cheek for duplicates. 
c     if a voxel is shared among multiple atoms then distribut it
C     value equaly amoung all atoms. 
************************************************************************
c     subroutine gistscorefuc_smart(atomstart, atomend, 
      subroutine gistscorefuc(atomstart, atomend, 
     &    coord_index, transfm_coords, orgin, gistspace,
     &    gistvolume,gistH,
     &    xnsize, ynsize, znsize, sra, srb, atom_vdwtype,
     &    gistgrid, tagvoxelarray, neighvoxelarray,
     &    maxatm, maxatmpos, maxtyv, score, totscore,
     &    gistindexlist, lengistil,
     &    gistindexlcon, lengistilc)

c XXX: These need to come first or Gfortran complains (TS)
c max constants, see max.h
      integer, intent(in)    :: maxatm !how many atoms max can be in a ligand
      integer, intent(in)    :: maxatmpos !how many positions all ligand atoms can take
      integer, intent(in)    :: maxtyv !how many vdw atom types there are
      !integer, intent(in)    :: nsize ! grid size 
      integer, intent(in)    :: xnsize ! grid size for x
      integer, intent(in)    :: ynsize ! grid size for y
      integer, intent(in)    :: znsize ! grid size for z
c input variables
      integer, intent(in)    :: atomstart, atomend !where the atoms and other things are
      integer, intent(in)    :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real,    intent(in)    :: transfm_coords(3, maxatm) !ligand atoms coordinates moved into bs
      real,    intent(in)    :: orgin(3) !the orgin is the front,lower left corner of the box
      real,    intent(in)    :: gistspace !angstroms per grid point on gist grid
      real,    intent(in)    :: gistvolume !user defined volume
      integer, intent(in)    :: gistH !if 0 use hydrogens as is, if 1 do not use hydrogens in gist calculations, if 2 replace nonpolar radii with 0.8

      real, intent(in) :: sra(maxtyv), srb(maxtyv) !square roots of vdw parameters per vdwtype
      integer, intent(in) :: atom_vdwtype(maxatm) !the vdw type of each ligand atom

      !real,    intent(in)    :: gistgrid(xnsize, ynsize, znsize) ! values
      real,    intent(in)    :: gistgrid(xnsize, ynsize, znsize) ! values
c return variable
      real,    intent(inout) :: score(maxatm) ! score at each atom
      real,    intent(inout) :: totscore      ! totscore

      integer, PARAMETER     :: neighMax = 1000 ! this the maxium number
                                               ! of voxels per-atom
      !integer                :: tag(3,neighMax*maxatm) ! that this aray is big enough this will remember which grid point have be taged.

      ! Pass this array.  It is intialized when gistgrids are read in. 
      LOGICAL, intent(inout) :: tagvoxelarray(xnsize*ynsize*znsize) !grid index (i,j,k) -> i+(j-1)*xn+(k-1)*xn*yn
                                                                    ! (1,1,1) -> 1 + 0*xn + 0*zn
                                                                    ! (xn,1,1) -> xn + 0 +0
                                                                    ! (1,yn,1)  ->  1+(yn-1)*xn+0 
                                                                    ! (xn,yn,zn) -> xn + (yn-1)*xn + (zn-1)*yn*xn = xn*yn*yn
      LOGICAL, intent(inout) :: neighvoxelarray(xnsize*ynsize*znsize) 
      integer                :: tag(neighMax*maxatm) ! we can remember the index, so that we can easaly wipe tagvoxelarray 
      !integer                :: atom_voxel_list(maxatm,neighMax) !  we will remember the index for each atom, so that we can easaly wipe the gird 

      integer                :: neighlist(3,neighMax) ! neighbor list needs to be big enough to store point in one atom. ever point has 6 neighbors. 
                                                       ! each of those neighbor will have 5 neighbors (not unique).
      integer                :: currentneigh(3,6) ! neighbor list needs to be is big enough store. ever point has 6 neighbors. 
      !integer                :: current_neighlist(3,30*maxatm) ! neighbor list needs to be is big enough store.  
      !integer                :: all_neighlist(3,30*maxatm) ! neighbor list needs to be is big enough store 
      integer                :: tagcount ! counts how meny are taged
      integer                :: neighcount, current, neighlenth ! for neighbor code
      integer                :: ni,nli                              ! neighbor index; neighbor list index
      LOGICAL                :: flagAdd                            !  this flag is for determining wether the neighbor should be added.
      LOGICAL                :: tagflag  
      integer atomindex !actual position used
      integer atomcount !atom position counter
      integer nx, ny, nz, i, idx, voxnum
      real xn(3)
      !real xgr, ygr, zgr, totscore
      real radius, dist, volume

      ! this is for debuging output
      real cx,cy,cz
      character (len=*), parameter :: PDBFORM = 
     &  '(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)'

      integer,     intent(inout) :: gistindexlist(5000)
      integer,     intent(inout) :: lengistil
      integer,     intent(in) :: gistindexlcon(5000)
      integer,     intent(in) :: lengistilc

c      ! here we intailize the array, we moved to this to grid_reader.f 
c      do idx = 1, xnsize*ynsize*znsize
c         if (tagvoxelarray(idx)) then
c            write(*,*) "tagvoxelarray(", idx, ") = ", 
c    &                  tagvoxelarray(idx)
c            call doflush(6)
c         endif
c         tagvoxelarray(idx) = .FALSE.
c      enddo

c     write(6,*) "IN subroutine gistscorefuc . . ."
c     call doflush(6)
c     write(6,*) "orgin = ", orgin(1),",",orgin(2),",",orgin(3)
c     write(6,*) "start=",atomstart, ";stop=",atomend
c     loop over all atoms

c     !! this was to check the grid, for debuging
c     do nx = 1, nsize
c        do ny = 1, nsize
c           do nz = 1, nsize
c              write(6,*) gistgrid(nx,ny,nz)
c           enddo
c        enddo
c     enddo
c     stop

c     write(*,*) "lengistil = ", lengistil
c     write(*,*) "lengistilc = ", lengistilc
c     call doflush(6)
      if (lengistil .gt. 5000) then 
         write(*,*) "Error. . . lengistil is too big."
         stop
      endif

      if (lengistilc .gt. 5000) then
         write(*,*) "lengistilc=", lengistilc
         write(*,*) "Error. . . lengistilc is too big."
         stop
      endif

      tagcount = 1
      lengistil = 0
      tag(1) = 0
      if (lengistilc .gt. 0) then
          do i = 1, lengistilc ! loop over the gistindexlist to remember what we allready look at in the preceding segments of the molcule. 
                               ! Note here we are loop over the connected conf (molcule segement) 
              !neighvoxelarray(gistindexlist(i)) = .true.
              if (gistindexlcon(i) > xnsize*ynsize*znsize) then
                  write(*,*) "lengistilc=", lengistilc
                  write(*,*) "i=", i
                  write(*,*) "gistindexlcon(i)=", gistindexlcon(i)
                  write(*,*) "Warrning... gistindexlcon(i) is too big."
                  call doflush(6)
                  cycle ! go to next grid index value
              else if (gistindexlcon(i) < 1) then
                  write(*,*) "i=", i
                  write(*,*) "gistindexlcon(i)=", gistindexlcon(i)
                  write(*,*) "Warrning... gistindexlcon(i) is small."
                  call doflush(6)
                  cycle ! go to next grid index value
                  !exit ! break out of the loop
                  !stop ! exit program
              endif
              tagvoxelarray(gistindexlcon(i))   = .true.
              tag(tagcount) = gistindexlcon(i)
              tagcount = tagcount + 1
              !if (lengistilc < 200) then ! if the connected molecule is small then include those point in the current 
              if (lengistilc < 1500) then ! if the connected molecule is small then include those point in the current 
                 lengistil = lengistil + 1
                 gistindexlist(lengistil) = gistindexlcon(i)
              else if (i > (lengistilc-1500+1)) then ! if it is biger only take the stuff at the bottum, the bottum should be closer to new frag
                 lengistil = lengistil + 1
                 gistindexlist(lengistil) = gistindexlcon(i)
              endif
          enddo
          !tagcount = lengistilc + 1
      endif
      
      totscore = 0.0

      if (gistvolume == 0.0) then
        volume = gistspace**3 ! this is the volume of a voxel. 
      else
        volume = gistvolume
      endif

c     write(6,*) "volume = ", volume
c     write(6,*) "atomstart=",atomstart,"atomend=",atomend

      do atomcount = atomstart, atomend

        atomindex = coord_index(2, atomcount) !gets actual atom non-consecutive #
        score(atomindex) = 0.0 ! always set atom score to zero frist then replace

        ! if hydrogen then skip it
        if ((gistH.eq.1) .and.((atom_vdwtype(atomindex).eq. 6) 
     &      .or. (atom_vdwtype(atomindex).eq. 7))) then
c           write(*,*) "skip Hydrogen."
            cycle ! go to the next atom
        endif
        ! calculate radii.!
        !     eps = B^2 / A
        !     R = (A/eps)^(1/12) / 2
        ! or  R = (B/(2*eps))^(1/6) / 2
        ! plug in the value
        !     R = (2*A/B)^(1/6) / 2 
        !       = (2^(1/2) * A^(1/2) / B^(1/2))^(1/3) / 2

        if (srb(atom_vdwtype(atomindex)).eq.0.00) then ! this should only be for dummy atom type 25 atomtype 0.   
            radius = 0.1
c       else if ((atom_vdwtype(atomindex).eq. 6)
c    &      .or. (atom_vdwtype(atomindex).eq. 7)) then
        else if ((gistH.eq.2).and.(atom_vdwtype(atomindex).eq. 7)) then ! 7 is nonpolar H
c           radius = (sqrt(2.0)* sra(atom_vdwtype(atomindex)) /
c    &               srb(atom_vdwtype(atomindex)))**(1.0/3.0) / 2
c           write(*,*) "hydrogen radius is ", radius
c           write(*,*) "use 0.8 instead"
            radius = 0.8
        else 
            radius = (sqrt(2.0)* sra(atom_vdwtype(atomindex)) /
     &               srb(atom_vdwtype(atomindex)))**(1.0/3.0) / 2 
        endif

c       write(*,*) "NEW ATOM" 
c       write(6,*) atom_vdwtype(atomindex), 
c    &             sra(atom_vdwtype(atomindex)), 
c    &             srb(atom_vdwtype(atomindex)), radius 

c       tranform ligand coord to the grid 
c       need to add one so the orgin x,y,z coreponds to grid point (1, 1, 1) not (0,0,0)
        xn(1) = (transfm_coords(1, atomindex) - 
     &       orgin(1)) / gistspace + 1
        xn(2) = (transfm_coords(2, atomindex) -
     &       orgin(2)) / gistspace + 1
        xn(3) = (transfm_coords(3, atomindex) - 
     &       orgin(3)) / gistspace + 1

c       write(6,*) "atom:", atomcount," -> ", atomindex

c       write(6,*) transfm_coords(1, atomindex), ",", 
c    &             transfm_coords(2, atomindex), ",",
c    &             transfm_coords(3, atomindex)

c       find closest grid point by rounding the transformed coordnates
        nx = NINT(xn(1))
        ny = NINT(xn(2))
        nz = NINT(xn(3))

        ! check that the atom is not outside the grid
        ! this will iginore the atom for now.
        ! the conformation is removed at later step.  see function gridsize_check in gridreader.f
        if (nx.lt.1.or.nx.gt.xnsize.or.
     &      ny.lt.1.or.ny.gt.ynsize.or.
     &      nz.lt.1.or.nz.gt.znsize)then
                 ! do nothing and
            !exit ! break out of loop 
            cycle ! continue to next atom
        endif
c       write(6,*) "gistgrid:",nx, ny, nz,gistgrid(nx,ny,nz)

c       we want all voxels that overlap with any ligand
c       atom, but do not want to double count.
c       tag insures that a voxel does not repeat for multiple
c       atoms
c
c       see if taged.  
c       Note that the voxel is added to only one atom that the code incountered frist.  This means that the atom decompotion is dependant on the atomic order.  
c       this is not ideal: we would like to distribute the overlaping voxals (that is voxal  to the relivent atoms in future.  

        ! next we need to look at neighboring grid points 
        ! to do this we need the radius which we now have see above. 
        ! we could do this with recursion but lets do it with a while
        ! loop. 
        score(atomindex) = 0.0
        neighcount  = 1 ! how meny elements remain on the list
        current = 1 ! were are we on the list
        neighlenth  = 1 ! the lenght of the list 
        neighlist(1,current) = nx ! the grid point closest to the atom center
        neighlist(2,current) = ny
        neighlist(3,current) = nz
        ! loop over grid points starting from the closted to atom center
        ! and going outward
c       write(6,*) "I AM HERE (1)"
        do while (neighcount > 0) ! while there are still grid points to check wether they are 
c           write(6,*) "I AM HERE (1.1)"
c            get_nieghbors(radius,neighbor_list,neighcount,current,
c     &                    listcount ) ! the function will return a 

            ! 1. get neighbors for curent neighbor
            currentneigh(1,1) = neighlist(1,current)+1 ! forward,x
            currentneigh(2,1) = neighlist(2,current)   ! forward,y
            currentneigh(3,1) = neighlist(3,current)   ! forward,z

            currentneigh(1,2) = neighlist(1,current)-1 ! backward,x
            currentneigh(2,2) = neighlist(2,current)   ! backward,y
            currentneigh(3,2) = neighlist(3,current)   ! backward,z

            currentneigh(1,3) = neighlist(1,current)   ! right,x
            currentneigh(2,3) = neighlist(2,current)+1 ! right,y
            currentneigh(3,3) = neighlist(3,current)   ! right,z

            currentneigh(1,4) = neighlist(1,current)   ! left,x
            currentneigh(2,4) = neighlist(2,current)-1 ! left,y
            currentneigh(3,4) = neighlist(3,current)   ! left,z

            currentneigh(1,5) = neighlist(1,current)   ! up,x
            currentneigh(2,5) = neighlist(2,current)   ! up,y
            currentneigh(3,5) = neighlist(3,current)+1 ! up,z

            currentneigh(1,6) = neighlist(1,current)   ! down,x
            currentneigh(2,6) = neighlist(2,current)   ! down,y
            currentneigh(3,6) = neighlist(3,current)-1 ! down,z

c           write(6,*) "I AM HERE (1.2)"
            do ni =  1,6
c              write(6,*) "I AM HERE (1.3)"
               flagAdd = .TRUE.
               ! 1.a Check that it is not already on the list.
               !     We speed this up in the sameway as 
               !     as the taged list (below) using a bool array. 
c              do nli = 1, neighlenth 
c                  if (currentneigh(1,ni).eq.neighlist(1,nli).AND.
c    &                 currentneigh(2,ni).eq.neighlist(2,nli).AND.
c    &                 currentneigh(3,ni).eq.neighlist(3,nli)) then
c                      flagAdd = .FALSE.
c                      exit ! leave loop early
c                  endif
c              enddo ! nli -- neighbor list index
               voxnum = currentneigh(1,ni) 
     &                  + xnsize*(currentneigh(2,ni)-1)
     &                  + xnsize*ynsize*(currentneigh(3,ni)-1)
               if (voxnum > xnsize*ynsize*znsize) then
                   ! not a grid point index is not on the array
                   ! we should not kept this pose. 
                   !exit ! break out of loop
                   cycle ! continue to next ni
               endif
               if (voxnum.eq.0) then
                   write(*,*) "Warning: voxnum.eq.0"
                   call doflush(6)
                   cycle ! continue to next ni
               endif

               flagAdd = .NOT.(neighvoxelarray(voxnum)) ! if it has been seen before then array val is true and flagAdd should be false

c              if (.NOT.flagAdd) then 
c                  ! if the voxel is alreay seen then continue to the
c                  ! next voxel, break the loop erly. 
c                  !write(*,*) "the voxel is alreay seen"
c                  !call doflush(6)
c                  write(*,*) voxnum, currentneigh(1,ni),
c    &  currentneigh(2,ni), currentneigh(3,ni)
c                  write(*,*) "flagAdd, neighvoxelarray(voxnum)",
c    &  flagAdd, neighvoxelarray(voxnum) 
c                  exit
c              endif
               ! 1.b check that it is not outside the radius of the atom
               !    real dist = sqrt((xn(1) - currentneigh(1,ni))**2
               dist = sqrt((xn(1) - currentneigh(1,ni))**2
     &                   + (xn(2) - currentneigh(2,ni))**2   
     &                   + (xn(3) - currentneigh(3,ni))**2 )  
               ! xn is the center of the atom transformed and
               ! normalized to the grid. we also need to normalize the 
               ! radius by the grid spacing.  
               if (dist > radius/gistspace) then
c                  write(6,*) "xn",xn(1),xn(2),xn(3) 
c                  write(6,*) "currentneigh",currentneigh(1,ni),
c    &                                       currentneigh(2,ni),
c    &                                       currentneigh(3,ni)
c                  write(6,*) "dist", dist,(radius/gistspace)
                   flagAdd = .FALSE.
               endif 
               ! 1.c check that neighboring points are not outside the grid box
               if (currentneigh(1,ni).lt.1.OR.
     &             currentneigh(1,ni).gt.xnsize.OR.
     &             currentneigh(2,ni).lt.1.OR.
     &             currentneigh(2,ni).gt.ynsize.OR.
     &             currentneigh(3,ni).lt.1.OR.
     &             currentneigh(3,ni).gt.znsize)then
                   flagAdd = .False.
                   ! if this is the case we might want to throw
                   ! the pose out.
               endif
               ! 1.d if good then add neighbor below the last element of the list.
               if (flagAdd) then
                   neighlenth = neighlenth + 1
                   if (neighlenth.gt.(neighMax)) then
                       write(6,*) "WARNING: neighlenth = ", neighlenth
                       write(6,*) "radius=",radius, "atomindex=", 
     &                             atomindex 
                       stop
                   endif
                   neighlist(1,neighlenth) = currentneigh(1,ni)
                   neighlist(2,neighlenth) = currentneigh(2,ni)
                   neighlist(3,neighlenth) = currentneigh(3,ni)

                   neighvoxelarray(voxnum) = .TRUE.
               endif
            enddo ! ni -- neighbor index
            ! 2. neighbor count = lenth - curent 
            current = current + 1
            neighcount = neighlenth - current
c           write(6,*) "neighbor info:",neighcount,neighlenth,current
        enddo ! while there are still neighbors that are not outside the of the atom radius 

      ! reset the neighbor bool array
        do nli = 1, neighlenth
c         taging is to insure that voxals do not overlap, no double
c         counting.
c         Note that this code is modified from some of the above code.
          voxnum = neighlist(1,nli) + xnsize*(neighlist(2,nli)-1)
     &                              + xnsize*ynsize*(neighlist(3,nli)-1)
          if (voxnum > xnsize*ynsize*znsize) then
              ! not a grid point index is not on the array
              ! we should not kept this pose. 
              !exit ! break out of loop
              cycle ! continue to next nli
          endif
          neighvoxelarray(voxnum) = .FALSE.
        enddo


      ! now add the grid point values in the atom to the score.
c       write(6,*) "atomcount=", atomcount,"atomindex=",atomindex,
c    &             ": gist point num = ", neighlenth,
c    &             "; radius=",radius  
        ! neighbor list index
        do nli = 1, neighlenth
c         taging is to insure that voxals do not overlap, no double
c         counting.
c         Note that this code is modified from some of the above code.
          voxnum = neighlist(1,nli) + xnsize*(neighlist(2,nli)-1)
     &                              + xnsize*ynsize*(neighlist(3,nli)-1)
          if (voxnum > xnsize*ynsize*znsize) then
              ! not a grid point index is not on the array
              ! we should not kept this pose. 
              !exit ! break out of loop
              cycle ! continue to next nli
          endif
          tagflag = tagvoxelarray(voxnum)
C         removed the do loop over the taged list. 
c         this should speed things up.
          if (.NOT.(tagflag)) then ! if the grid point is taged then do
c                                    not add score
              ! we can speed this up slightly by multipling the volume
              ! at the end. 
              !score = score + (volume of voxel A^3 ) * (voxel energy kcal/mol/A^3)
              score(atomindex) = score(atomindex) + (volume) * 
     &                           gistgrid( neighlist(1,nli), 
     &                                     neighlist(2,nli),
     &                                     neighlist(3,nli))
              !totscore = totscore + score(atomindex) ## this was a bug sum a sum. 
              totscore = totscore + (volume) *
     &                   gistgrid( neighlist(1,nli),
     &                             neighlist(2,nli),
     &                             neighlist(3,nli))

              tag(tagcount) = voxnum
              tagvoxelarray(voxnum) = .TRUE.
              tagcount = tagcount +1
              ! if (tagcount > )
              if (tagcount.gt.(neighMax*maxatm)) then
                  write(6,*) "Error: tagcount = ", tagcount
                  stop
              endif

              ! store index in gistindexlist so that we can remember. 
              lengistil = lengistil+1 ! starts at zerro, so add one frist then set value.
              if (lengistil .gt. 5000) then 
                 write(*,*) "Error. . . lengistil is too big."
                 stop
              endif
              gistindexlist(lengistil) = voxnum

          endif
        enddo ! neighbor list index
c       stop
      enddo ! atom list

      ! We will now wipe the array clean,this is nesarary when we pass
      ! the tagvoxelarray, which is preintialized. 
      do idx = 1, (tagcount -1)
         !write(*,*) "tag(i)=", tag(idx) 
         !call doflush(6)
         !write(*,*) "tagvoxelarray=", tagvoxelarray(tag(idx)) 
         !call doflush(6)
         tagvoxelarray(tag(idx)) = .FALSE.
      enddo

cc    write(6,*) "Tot Gist Score =  ", totscore
cc    write(6,*) "Tot points =  ", tagcount
cc    write(6,*) "Tot atoms =  ",  (atomend-atomstart)
c     write(6,*) "Tot Gist Score =  ", totscore,
c    &           "; Tot points =  ", tagcount,
c    &           "; Tot atoms =  ",  (atomend-atomstart)

c     write(6,*) "OUT subroutine gistscorefuc . . . "

      return
      !end subroutine gistscorefuc_smart
      end subroutine gistscorefuc

************************************************************************
c   This modual is created to caluculate the desolvation of the
c   Receptor using Grid Inhomogeneous Solvation Theory (GIST)
c   (Nguyen et al. J. Chem. Phys. 137 044101 2012 )
c   This is using and approxamation to speed up the caluclation:  
c   only look at the clossest grid point.

c   we can try to use corser grid. 

c
c   ANOTHER IDEA:  
c   Is to use this function on modified grid where
c   each grid poit is a sphere. 
c   sum up voxels in the sphere center on each voxel. 
c   the grids will then no longer be non overlaping.
c   this is a big aproxamation, we may want to scale depending on number
c   of bonds.  that is, if atom has 4 bonds divid by 4.  This will
c   inpart acount for the overlap region, in a mean aproxamation.  
c   This subroutine was writen by Trent Balius
************************************************************************
      subroutine gistscorefuc_fast_1(atomstart, atomend,
     &    coord_index, transfm_coords, orgin, gistspace,
     &    gistvolume,
     &    xnsize, ynsize, znsize, sra, srb, atom_vdwtype,
     &    gistgrid, maxatm, maxatmpos, maxtyv, score, totscore)

      integer, intent(in)    :: maxatm !how many atoms max can be in a ligand
      integer, intent(in)    :: maxatmpos !how many positions all ligand atoms can take
      integer, intent(in)    :: maxtyv !how many vdw atom types there are
      !integer, intent(in)    :: nsize ! grid size
      integer, intent(in)    :: xnsize ! grid size for x
      integer, intent(in)    :: ynsize ! grid size for y
      integer, intent(in)    :: znsize ! grid size for z
c input variables
      integer, intent(in)    :: atomstart, atomend !where the atoms and other things are
      integer, intent(in)    :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real,    intent(in)    :: transfm_coords(3, maxatm) !ligand atoms coordinates moved into bs
      real,    intent(in)    :: orgin(3) !the orgin is the front,lower left corner of the box
      real,    intent(in)    :: gistspace !angstroms per grid point on gist grid
      real,    intent(in)    :: gistvolume !user defined volume

      real, intent(in) :: sra(maxtyv), srb(maxtyv) !square roots of vdw parameters per vdwtype
      integer, intent(in) :: atom_vdwtype(maxatm) !the vdw type of each ligand atom

      !real,    intent(in)    :: gistgrid(xnsize, ynsize, znsize) !
      !values
      real,    intent(in)    :: gistgrid(xnsize, ynsize, znsize) !  values c return variable
      real,    intent(inout) :: score(maxatm) ! score at each atom
      real,    intent(inout) :: totscore      ! totscore

      integer                :: currentneigh(3,6) ! neighbor list needs to be is big enough store. ever point has 6 neighbors.
      integer atomindex !actual position used
      integer atomcount !atom position counter
      integer nx, ny, nz, i
      real xn(3)
      !real xgr, ygr, zgr, totscore
      !real radius, dist !, volume
      real radius, dist , volume

      ! this is for debuging output
      real cx,cy,cz

      !volume = gistspace**3 ! this is the volume of a voxel. 
      ! we might want to change this to be sphere
      ! scale the volume outside of DOCK. (TEB)
      if (gistvolume == 0.0) then
        volume = gistspace**3 ! this is the volume of a voxel. 
      else
        volume = gistvolume
      endif


      totscore = 0.0

      do atomcount = atomstart, atomend

        atomindex = coord_index(2, atomcount) !gets actual atom non-consecutive #
        score(atomindex) = 0.0 ! always set atom score to zero frist then replace

c    Tranform ligand coord to the grid
c    need to add one so the orgin x,y,z coreponds to grid point (1, 1, 1) not (0,0,0)
        xn(1) = (transfm_coords(1, atomindex) -
     &       orgin(1)) / gistspace + 1
        xn(2) = (transfm_coords(2, atomindex) -
     &       orgin(2)) / gistspace + 1
        xn(3) = (transfm_coords(3, atomindex) -
     &       orgin(3)) / gistspace + 1

c       find closest grid point by rounding the transformed coordnates
        nx = NINT(xn(1))
        ny = NINT(xn(2))
        nz = NINT(xn(3))

        ! check that the atom is not outside the grid
        ! this will eginore the atom for now.
        ! the conformation is removed at later step.  see function
        ! gridsize_check in gridreader.f
        if (nx.lt.1.or.nx.gt.xnsize.or.
     &      ny.lt.1.or.ny.gt.ynsize.or.
     &      nz.lt.1.or.nz.gt.znsize)then
                 ! do nothing and
            exit ! leave loop, go to next element
        endif

        score(atomindex) = score(atomindex)+volume*gistgrid(nx,ny,nz)
        totscore = totscore+volume*gistgrid(nx,ny,nz)
        !score(atomindex) = score(atomindex) + gistgrid(nx,ny,nz) !  gist grid should be scaled by volumn outside.  
        !totscore = totscore + gistgrid(nx,ny,nz)
        !consider scaling by valance, this will scale by the number of
        !atoms conected to.  this will inpart correct for double (or
        !multi) counting in a mean way.   

      enddo ! atom list

      end subroutine gistscorefuc_fast_1

************************************************************************
      subroutine gistscorefuc_fast_8(atomstart, atomend,
     &    coord_index, transfm_coords, orgin, gistspace,
     &    gistvolume,
     &    xnsize, ynsize, znsize, sra, srb, atom_vdwtype,
     &    gistgrid, maxatm, maxatmpos, maxtyv, score, totscore)

      integer, intent(in)    :: maxatm !how many atoms max can be in a ligand
      integer, intent(in)    :: maxatmpos !how many positions all ligand atoms can take
      integer, intent(in)    :: maxtyv !how many vdw atom types there are
      !integer, intent(in)    :: nsize ! grid size
      integer, intent(in)    :: xnsize ! grid size for x
      integer, intent(in)    :: ynsize ! grid size for y
      integer, intent(in)    :: znsize ! grid size for z
c input variables
      integer, intent(in)    :: atomstart, atomend !where the atoms and other things are
      integer, intent(in)    :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real,    intent(in)    :: transfm_coords(3, maxatm) !ligand atoms coordinates moved into bs
      real,    intent(in)    :: orgin(3) !the orgin is the front,lower left corner of the box
      real,    intent(in)    :: gistspace !angstroms per grid point on gist grid
      real,    intent(in)    :: gistvolume !angstroms per grid point on gist grid

      real, intent(in) :: sra(maxtyv), srb(maxtyv) !square roots of vdw parameters per vdwtype
      integer, intent(in) :: atom_vdwtype(maxatm) !the vdw type of each ligand atom

      !real,    intent(in)    :: gistgrid(xnsize, ynsize, znsize) !
      !values
      real,    intent(in)    :: gistgrid(xnsize, ynsize, znsize) !  values c return variable
      real,    intent(inout) :: score(maxatm) ! score at each atom
      real,    intent(inout) :: totscore      ! totscore

      integer                :: currentneigh(3,6) ! neighbor list needs to be is big enough store. ever point has 6 neighbors.
      integer atomindex !actual position used
      integer atomcount !atom position counter
      integer nxc, nyc, nzc, nxf, nyf, nzf, i
      real xn(3)
      !real xgr, ygr, zgr, totscore
      !real radius, dist !, volume
      real radius, dist , volume, gs

      ! this is for debuging output
      real cx,cy,cz

      !volume = gistspace**3 ! this is the volume of a voxel. 
      ! we might want to change this to be sphere
      ! scale the volume outside of DOCK. (TEB)
      if (gistvolume == 0.0) then
        volume = gistspace**3 ! this is the volume of a voxel. 
      else
        volume = gistvolume
      endif


      do atomcount = atomstart, atomend

        atomindex = coord_index(2, atomcount) !gets actual atom non-consecutive #
        score(atomindex) = 0.0 ! always set atom score to zero frist then replace

c    Tranform ligand coord to the grid
c    need to add one so the orgin x,y,z coreponds to grid point (1, 1, 1) not (0,0,0)
        xn(1) = (transfm_coords(1, atomindex) -
     &       orgin(1)) / gistspace + 1
        xn(2) = (transfm_coords(2, atomindex) -
     &       orgin(2)) / gistspace + 1
        xn(3) = (transfm_coords(3, atomindex) -
     &       orgin(3)) / gistspace + 1

c       find closest 8 grid points by using ceil and floor the transformed coordnates
        nxc = CEILING(xn(1))
        nyc = CEILING(xn(2))
        nzc = CEILING(xn(3))
        nxf = FLOOR(xn(1))
        nyf = FLOOR(xn(2))
        nzf = FLOOR(xn(3))

        ! check that the atom is not outside the grid
        ! this will iginore the atom for now.
        ! the conformation is removed at later step.  see function
        ! gridsize_check in gridreader.f
        if (nxc.lt.1.or.nxc.gt.xnsize.or.
     &      nyc.lt.1.or.nyc.gt.ynsize.or.
     &      nzc.lt.1.or.nzc.gt.znsize.or.
     &      nxf.lt.1.or.nxf.gt.xnsize.or.
     &      nyf.lt.1.or.nyf.gt.ynsize.or.
     &      nzf.lt.1.or.nzf.gt.znsize)then
                 ! do nothing and
            !exit ! leave loop, go to next element
            cycle ! continue to next atom
        endif

c      (2)----(6)  !  These are the 8 closest grid points 
c     / |    / |   !  to a point (*) 
c   (4)----(8) |   !
c    |  | * |  |   !
c    | (1)----(5)  !
c    |/     |/
c   (3)----(7)

        gs = 0
        gs = gs+volume*gistgrid(nxf,nyf,nzf) ! (1) f,f,f
        gs = gs+volume*gistgrid(nxf,nyf,nzc) ! (2) f,f,c
        gs = gs+volume*gistgrid(nxf,nyc,nzf) ! (3) f,c,f
        gs = gs+volume*gistgrid(nxf,nyc,nzc) ! (4) f,c,c
        gs = gs+volume*gistgrid(nxc,nyf,nzf) ! (5) c,f,f
        gs = gs+volume*gistgrid(nxc,nyf,nzc) ! (6) c,f,c
        gs = gs+volume*gistgrid(nxc,nyc,nzf) ! (7) c,c,f
        gs = gs+volume*gistgrid(nxc,nyc,nzc) ! (8) c,c,c

        score(atomindex) = score(atomindex)+gs
        totscore = totscore+gs
        !score(atomindex) = score(atomindex) + gistgrid(nx,ny,nz) !  gist grid should be scaled by volumn outside.  
        !totscore = totscore + gistgrid(nx,ny,nz)
        !consider scaling by valance, this will scale by the number of
        !atoms conected to.  this will inpart correct for double (or
        !multi) counting in a mean way.   

      enddo ! atom list

      end subroutine gistscorefuc_fast_8


      end module gistscore
c----------------------------------------------------------------------
c
c       Copyright (C) 2014 Trent Balius, Marcus Fisher and Brian K. Shoichet
c                        University of Toronto / UCSF
c                         All Rights Reserved.
c----------------------------------------------------------------------
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
