
c     Jan 7 2014, Trent E Balius modified to read in Water information
c     GIST. 

      module grid_reader !reads all the various grids. also has check routine.

      implicit none
 
c variables used throughout module here, first are the file numbers foreach file
      integer, parameter :: DELPHIGRID = 99
      integer, parameter :: GISTFILE = 8
      integer, parameter :: RDESGRID = 10
      integer, parameter :: HRDESGRID = 11
      integer, parameter :: VDWGRID = 20
      integer, parameter :: SOLGRIDU = 31
      integer, parameter :: HSOLGRIDU = 32
      integer, parameter :: BUMPGRID = 25
c min & max grid limits, temporary and only used here
      real, dimension(3) :: phiminlimits
      real, dimension(3) :: phimaxlimits
      real, dimension(3) :: gistminlimits
      real, dimension(3) :: gistmaxlimits
      real, dimension(3) :: rec_d_rminlimits
      real, dimension(3) :: rec_d_rmaxlimits
      real, dimension(3) :: rdesrminlimits
      real, dimension(3) :: rdesmaxlimits
      real, dimension(3) :: recminlimits
      real, dimension(3) :: recmaxlimits
      real, dimension(3) :: solminlimits
      real, dimension(3) :: solmaxlimits
      real, dimension(3) :: vdwminlimits
      real, dimension(3) :: vdwmaxlimits
c this parameter makes sure the vdw term doesn't get too large
      real, parameter :: toobig = 100000.0

      contains !several subroutines

c this just calls the other 3 functions to read the various grids
      subroutine read_grids(grids0, OUTDOCK, 
     &    allminlimits, allmaxlimits, 
     &    phimap0, recdes0, solvmap0, vdwmap0, ligdes0, gistparm0,
     &    rec_des_r0, options0, which_grid)

      use phimaptype
      use solvmaptype
      use vdwmaptype
      use optionstype
      use gridstype
      use gisttype

      !write(6,*) "I AM READING GRIDS"

c XXX: GFortran needs these to come first
c grid min & max limits
      real, dimension(3), intent(out) :: allminlimits
      real, dimension(3), intent(out) :: allmaxlimits
c these are the grids, where all the data read in from the files will end up.
      type(grids), intent(out) :: grids0
c outdock file specifier, so data can be written out if there are problems
      integer, intent(in) :: OUTDOCK
c types, containing various options
      type(phimap), intent(inout) :: phimap0, recdes0
      type(gistparm), intent(inout) :: gistparm0, rec_des_r0
      type(solvmap), intent(inout) :: solvmap0
      type(vdwmap), intent(inout) :: vdwmap0
      type(phimap), intent(inout) :: ligdes0
      type(options), intent(inout) :: options0
c which grid to read in, since we support multiple flexible receptor grids now
      integer, intent(in) :: which_grid

      call read_vdw(options0, grids0%abval_precomp, OUTDOCK, 
     &    vdwmap0, which_grid)

      call read_delphi(grids0%phimap_precomp, OUTDOCK, phimap0, 
     &   which_grid, options0%phinames(which_grid), 
     &   phiminlimits, phimaxlimits)

      ! if we do not read in the gist grid it remains an unallocated.
      ! just a pointer.
      if (options0%gistflag(which_grid)) then
        call read_GIST_grid(grids0%gistgrid,
     &     grids0%gisttagvoxelarray, 
     &     grids0%gistneighborvoxelarray, 
     &     OUTDOCK, gistparm0, 
     &     which_grid, options0%gistnames(which_grid),
     &     gistminlimits, gistmaxlimits)
      endif

c      if (options0%rec_des_r.and.options0%hrec_des_r) then
      if (options0%rec_d_flag(which_grid) .and.
     &    options0%hrec_d_flag(which_grid)) then
        call read_DX_grid(grids0%rec_des_precomp, 
     &     grids0%rec_des_grid, options0%recdesrnames(which_grid),
     &     grids0%hrec_des_precomp, grids0%hrec_des_grid, 
     &     options0%hrecdesnames(which_grid), OUTDOCK, rec_des_r0, 
     &     which_grid,rec_d_rminlimits, rec_d_rmaxlimits)
      endif

      if (options0%receptor_desolv) then
        call read_delphi(grids0%recdes_precomp, OUTDOCK, recdes0, 
     &     which_grid, 
     &     options0%recdesnames(which_grid),
     &     recminlimits, recmaxlimits)
      endif

      if (options0%ldesolv .eq. 2) then
        call readsol(grids0%solgrd_precomp, grids0%hsolgrd_precomp,
     &       OUTDOCK, solvmap0, options0, which_grid)
      else if (options0%ldesolv .eq. 3) then !PB/qnifft ligand desolvation
        call read_delphi(grids0%ligdes_precomp, 
     &       OUTDOCK, ligdes0, which_grid, 
     &       options0%solvnames(which_grid), 
     &       solminlimits, solmaxlimits)
      endif

      call og_precomputelimits(allminlimits, allmaxlimits, 
     &    options0%ldesolv, options0%receptor_desolv, phimap0%nsize, 
c    &    gistparm0%nsize)
c    &    gistparm0%xnsize) ! replace with min
     &    options0%gistflag(which_grid), 
     &    (options0%rec_d_flag(which_grid)
     &    .and.options0%hrec_d_flag(which_grid))) ! replace with min
 
      return
      end subroutine read_grids


c this subroutine reads the GIST files water for receptor desolvation
c corlation function solute (receptor) and water
c writen by Trent Balius jan 07, 2014

      subroutine read_DX_grid(rec_des_precomp, rec_d_grid, 
     &    rec_des_filename, hrec_des_precomp, hrec_d_grid, 
     &    hrec_des_filename, OUTDOCK, rec_des_r0, which_grid, 
     &    theminlimits, themaxlimits)

      use errorformats
      use gisttype      

      integer, intent(in) :: OUTDOCK
      integer, intent(in) :: which_grid
      character (len=255), intent(in) :: rec_des_filename
      character (len=255), intent(in) :: hrec_des_filename
      real, dimension(3), intent(inout) :: theminlimits, themaxlimits
      type(gistparm), intent(inout) :: rec_des_r0

      real, allocatable, dimension(:, :, :, :), intent(out) ::
     &          rec_des_precomp
      real, allocatable, dimension(:, :, :, :), intent(out) ::
     &          hrec_des_precomp
      real, allocatable, dimension(:, :, :), intent(out) :: rec_d_grid
      real, allocatable, dimension(:, :, :), intent(out) :: hrec_d_grid
      real, allocatable, dimension(:) :: rec_d_gridarray
      real, allocatable, dimension(:) :: hrec_d_gridarray

      integer istat

c     character variables to read form gist_gsw
      character (len=20) :: text20
      character (len=15) :: text15
      character (len=10) :: text10
      character (len=5)  :: text05
      character (len=6)  :: text06
      character (len=60) :: text60
      character (len=*), parameter :: PDBFORM =
     &  '(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)'

      integer nx, ny, nz, coord, tempint, idx
      integer xngrid, yngrid, zngrid, arin ! this is the demitions of the grid 
      integer xnsize, ynsize, znsize ! this is the demitions of the grid 
      integer hxngrid, hyngrid, hzngrid, harin
      real Xorg, Yorg, Zorg ! the orgin of the grid
      real hXorg, hYorg, hZorg ! the orgin of the grid
      real rgrid
      real dx,dy,dz,dummy
      real hdx,hdy,hdz,hdummy
      real a1, a2, a3, a4, a5, a6, a7, a8

        open(RDESGRID, file=trim(rec_des_filename),
     &      status='old', form='formatted', action='read') 
        read(RDESGRID,"(A20,A10,A6,I3,I3,I3)") text20, text10, text06,
     &       xngrid, yngrid, zngrid
        write(OUTDOCK,*) text20,text10,text06,
     &       xngrid, yngrid, zngrid
        read(RDESGRID,*) text10, Xorg, Yorg, Zorg
        read(RDESGRID,*) text10, dx, dummy, dummy
        read(RDESGRID,*) text10, dummy, dy, dummy
        read(RDESGRID,*) text10, dummy, dummy, dz
        read(RDESGRID,'(A60)') text60
        read(RDESGRID,'(A60)') text60

        write(OUTDOCK,*) 'dim =',xngrid,yngrid,zngrid
        write(OUTDOCK,*) 'spacing =',dx,dy,dz
        write(OUTDOCK,*) 'origin ', Xorg,Yorg,Zorg
        
        rec_des_r0%xnsize = xngrid
        rec_des_r0%ynsize = yngrid
        rec_des_r0%znsize = zngrid
        if (dx.ne.dy .or. dx.ne.dz .or. dy.ne.dz) then
          write(OUTDOCK, HALT3) 'Gist step size are inconsitant',
     &    which_grid
          stop
        endif
        rec_des_r0%gistspace = dx
        rec_des_r0%orgin(1) = Xorg
        rec_des_r0%orgin(2) = Yorg
        rec_des_r0%orgin(3) = Zorg
        
        allocate(rec_d_grid(xngrid, yngrid, zngrid),
     &      stat=istat)
        if (istat .ne. 0) then
          write(OUTDOCK, HALT3) 'Failed to allocate rec_d_grid',
     &    which_grid
          stop
        endif

        allocate(rec_d_gridarray(xngrid*yngrid*zngrid),
     &      stat=istat)
        if (istat .ne. 0) then
          write(OUTDOCK, HALT3) 'Failed to allocate rec_d_grid',
     &    which_grid
          stop
        endif

        read(RDESGRID,*) rec_d_gridarray
        close(RDESGRID)

        arin = 1 ! array index
        do nx = 1, xngrid
          do ny = 1, yngrid
            do nz = 1, zngrid
               rec_d_grid(nx,ny,nz) = rec_d_gridarray(arin)
               arin = arin + 1
            enddo
          enddo
        enddo

        do coord = 1,3
          theminlimits(coord) = rec_des_r0%orgin(coord)
        enddo

        ! subtract one because otherwise sometimes it seems to get
        ! unreasonable values.  note that we should sustact one because
        ! of starting at zerro, not one.  
        themaxlimits(1) = (rec_des_r0%orgin(1) +
     &                        (real(rec_des_r0%xnsize-1) *
     &                        rec_des_r0%gistspace))
        themaxlimits(2) = (rec_des_r0%orgin(2) +
     &                        (real(rec_des_r0%ynsize-1) *
     &                        rec_des_r0%gistspace))
        themaxlimits(3) = (rec_des_r0%orgin(3) +
     &                        (real(rec_des_r0%znsize-1) *
     &                        rec_des_r0%gistspace))

        write(OUTDOCK,*) 'maxlimits are',
     &        themaxlimits(1),themaxlimits(2),themaxlimits(3)

        allocate(rec_des_precomp(8,rec_des_r0%xnsize, rec_des_r0%ynsize,
     &      rec_des_r0%znsize), stat=istat)
        if (istat .ne. 0) then
           write(OUTDOCK, HALT3) 'Failed to allocate rdmap', which_grid
           stop
        endif

        allocate(hrec_des_precomp(8,xngrid,yngrid,zngrid), stat=istat)
        if (istat .ne. 0) then
           write(OUTDOCK, HALT3) 'Failed to allocate rdmap', which_grid
           stop
        endif

        open(HRDESGRID, file=trim(hrec_des_filename),
     &      status='old', form='formatted', action='read') 
        read(HRDESGRID,"(A20,A10,A6,I3,I3,I3)") text20, text10, text06,
     &       hxngrid, hyngrid, hzngrid
        read(HRDESGRID,*) text10, hXorg, hYorg, hZorg
        read(HRDESGRID,*) text10, hdx, hdummy, hdummy
        read(HRDESGRID,*) text10, hdummy, hdy, hdummy
        read(HRDESGRID,*) text10, hdummy, hdummy, hdz
        read(HRDESGRID,'(A60)') text60
        read(HRDESGRID,'(A60)') text60
        
        write(OUTDOCK,*) 'hdim =',hxngrid,hyngrid,hzngrid
        write(OUTDOCK,*) 'hspacing =',hdx,hdy,hdz
        write(OUTDOCK,*) 'horigin ', hXorg,hYorg,hZorg        
        call doflush(OUTDOCK)

        if ((hXorg .ne. Xorg) .or. (hYorg .ne. Yorg) .or.
     &      (hZorg .ne. Zorg) .or. (xngrid .ne. hxngrid) .or.
     &      (yngrid .ne. hyngrid) .or. (zngrid .ne. hzngrid) .or.
     &      (dx .ne. hdx) .or. (dy .ne. hdy) .or. (dz .ne. hdz)) then
          write(OUTDOCK, HALT1)
     &        'The hydrogen rec desolv map dimensions or '
          write(OUTDOCK, HALT1)
     &        'or parameters do not match the others'
          stop
        endif
        
        allocate(hrec_d_grid(xngrid, yngrid, zngrid),
     &      stat=istat)
        if (istat .ne. 0) then
          write(OUTDOCK, HALT3) 'Failed to allocate hrec_d_grid',
     &    which_grid
          stop
        endif

        allocate(hrec_d_gridarray(xngrid*yngrid*zngrid),
     &      stat=istat)
        if (istat .ne. 0) then
          write(OUTDOCK, HALT3) 'Failed to allocate hrec_d_grid',
     &    which_grid
          stop
        endif

        read(HRDESGRID,*) hrec_d_gridarray
        close(HRDESGRID)

        arin = 1 ! array index
        do nx = 1, xngrid
          do ny = 1, yngrid
            do nz = 1, zngrid
               hrec_d_grid(nx,ny,nz) = hrec_d_gridarray(arin)
               !write(6,*), hrec_d_grid(nx,ny,nz)
               !call doflush(6)
               arin = arin + 1
            enddo
          enddo
        enddo

        do nx = 1, rec_des_r0%xnsize-1
          do ny = 1, rec_des_r0%ynsize-1
            do nz = 1, rec_des_r0%znsize-1
             a8 = rec_d_grid(nx,ny,nz)
             a7 = rec_d_grid(nx,ny,nz+1) - a8
             a6 = rec_d_grid(nx,ny+1,nz) - a8
             a5 = rec_d_grid(nx+1,ny,nz) - a8
             a4 = rec_d_grid(nx,ny+1,nz+1) - a7 - rec_d_grid(nx,ny+1,nz)
             a3 = rec_d_grid(nx+1,ny,nz+1) - a7 - rec_d_grid(nx+1,ny,nz)
             a2 = rec_d_grid(nx+1,ny+1,nz) - a6 - rec_d_grid(nx+1,ny,nz)
             a1 = rec_d_grid(nx+1,ny+1,nz+1) - a4
     &                - rec_d_grid(nx+1,ny,nz+1)
     &                + rec_d_grid(nx+1,ny,nz)
     &                - rec_d_grid(nx+1,ny+1,nz)
             rec_des_precomp(1, nx, ny, nz) = a1
             rec_des_precomp(2, nx, ny, nz) = a2
             rec_des_precomp(3, nx, ny, nz) = a3
             rec_des_precomp(4, nx, ny, nz) = a4
             rec_des_precomp(5, nx, ny, nz) = a5
             rec_des_precomp(6, nx, ny, nz) = a6
             rec_des_precomp(7, nx, ny, nz) = a7
             rec_des_precomp(8, nx, ny, nz) = a8
            enddo
          enddo
        enddo  
        !write(6,*) 'rec_des has been precomputed'
        !call doflush(6)

        !write(6,*) "now running hydrogen"
        !call doflush(6)
        do nx = 1, xngrid-1
          do ny = 1, yngrid-1
            do nz = 1, zngrid-1
             a8 = hrec_d_grid(nx,ny,nz)
             a7 = hrec_d_grid(nx,ny,nz+1) - a8
             a6 = hrec_d_grid(nx,ny+1,nz) - a8
             a5 = hrec_d_grid(nx+1,ny,nz) - a8
             a4 = hrec_d_grid(nx,ny+1,nz+1)-a7-hrec_d_grid(nx,ny+1,nz)
             a3 = hrec_d_grid(nx+1,ny,nz+1)-a7-hrec_d_grid(nx+1,ny,nz)
             a2 = hrec_d_grid(nx+1,ny+1,nz)-a6-hrec_d_grid(nx+1,ny,nz)
             a1 = hrec_d_grid(nx+1,ny+1,nz+1)-a4
     &                - hrec_d_grid(nx+1,ny,nz+1)
     &                + hrec_d_grid(nx+1,ny,nz)
     &                - hrec_d_grid(nx+1,ny+1,nz)
             hrec_des_precomp(1, nx, ny, nz) = a1
             hrec_des_precomp(2, nx, ny, nz) = a2
             hrec_des_precomp(3, nx, ny, nz) = a3
             hrec_des_precomp(4, nx, ny, nz) = a4
             hrec_des_precomp(5, nx, ny, nz) = a5
             hrec_des_precomp(6, nx, ny, nz) = a6
             hrec_des_precomp(7, nx, ny, nz) = a7
             hrec_des_precomp(8, nx, ny, nz) = a8
             !write(6,*), nx,ny,nz
             !write(6, '(f,f,f,f,f,f,f,f)'),a8,a7,a6,a5,a4,a3,a2,a1
            enddo
          enddo
        enddo
        write(6,*) "hydrogen has been completed"
        call doflush(6)

      return
      end subroutine read_DX_grid

      subroutine read_GIST_grid(gistgrid, tagvoxelarray, 
     &    neighborvoxelarray,OUTDOCK, 
     &    gistparm0, which_grid, gistfilename, 
     &    theminlimits, themaxlimits)

      use errorformats
      use gisttype

      integer, intent(in) :: OUTDOCK
      integer, intent(in) :: which_grid
      character (len=255), intent(in) :: gistfilename
      real, dimension(3), intent(inout) :: theminlimits, themaxlimits
      type(gistparm), intent(inout) :: gistparm0

c     gist  grid
      real, allocatable, dimension(:, :, :), intent(out) :: gistgrid
c     gistgrid stores the x, y, z and value

      real, allocatable, dimension(:) :: gistgridarray
c     gistgridarray stores the value this will be transfered in gistgrid

      logical, allocatable, dimension(:), intent(out) :: tagvoxelarray  ! This is a boolen array.  this will be used in the scoring function
                                                           ! the array is initalized here.   
      logical, allocatable, dimension(:), intent(out) :: 
     & neighborvoxelarray  ! This is a boolen array.  this will be used in the scoring function

      integer istat !temporary status

c     character variables to read form gist_gsw
      character (len=20) :: text20
      character (len=15) :: text15
      character (len=10) :: text10
      character (len=5)  :: text05
      character (len=6)  :: text06
      character (len=60) :: text60
      character (len=*), parameter :: PDBFORM = 
     &  '(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)'

      integer nx, ny, nz, coord, tempint, idx
      real cx, cy, cz ! convert int nx to coor cx
      integer xngrid, yngrid, zngrid, arin ! this is the demitions of the grid 
      real Xorg, Yorg, Zorg ! the orgin of the grid
      real rgrid
      real dx,dy,dz,dummy
c example file formate.  
c object 1 class gridpositions counts 40 40 40
c origin 34.463 28.671 18.222
c delta 0.5 0 0
c delta 0 0.5 0
c delta 0 0 0.5
c object 2 class gridconnections counts 40 40 40
c object 3 class array type float rank 0 items 64000 data follows
c 0 0 0
c 0 0 0
c 0 0 0
c 0 0 0
c 0 0 0
c 0 0 0
c 0 0 0
c 0 0 0
c 0 0 0
c-0.0233089 -0.210245 -0.0425238
c
cfirst, if '0' is the filename, this means we don't actually want to open a file
c or allocate any of the dynamic arrays, we assume the energy of all atoms is 0
c in all places. this is useful for saving memory instead of reading files of 0s
cnote that for many reasons you can't use 0 as the first grid and then use non
czero as the later grids
      if (trim(gistfilename) .ne. '0') then

        open(GISTFILE, file=trim(gistfilename), 
     &      status='old', form='formatted', action='read')
c    &      status='old', form='unformatted', action='read')
        read(GISTFILE,"(A20,A10,A6,I3,I3,I3)") text20, text10, text06,
     &       xngrid, yngrid, zngrid
c       write(OUTDOCK,*) text15,text10,text06,
        write(OUTDOCK,*) text20,text10,text06,
     &       xngrid, yngrid, zngrid
c        stop
        !read(GISTFILE,"(A6,F7.3,F7.3,F7.3)") text10, Xorg, Yorg, Zorg
        !read(GISTFILE,"(A6,F3.1,F2.1,F2.1)") text10, dx, dummy, dummy
        !read(GISTFILE,"(A6,F2.1,F3.1,F2.1)") text10, dummy, dy, dummy
        !read(GISTFILE,"(A6,F2.1,F2.1,F3.1)") text10, dummy, dummy, dz
        read(GISTFILE,*) text10, Xorg, Yorg, Zorg
        read(GISTFILE,*) text10, dx, dummy, dummy
        read(GISTFILE,*) text10, dummy, dy, dummy
        read(GISTFILE,*) text10, dummy, dummy, dz
        read(GISTFILE,'(A60)') text60
        read(GISTFILE,'(A60)') text60
        
        write(OUTDOCK,*) 'dim =',xngrid,yngrid,zngrid
        write(OUTDOCK,*) 'spaceing =',dx,dy,dz
        write(OUTDOCK,*) 'origin ', Xorg,Yorg,Zorg
*        stop

c       if (xngrid.ne.yngrid.or.xngrid.ne.zngrid.or.
c    &   yngrid.ne.zngrid) then
c         write(OUTDOCK, HALT3) 'Gist demitions are not square', 
c    &    which_grid
c         stop
c       endif
        !gistparm0%nsize = xngrid
        gistparm0%xnsize = xngrid
        gistparm0%ynsize = yngrid
        gistparm0%znsize = zngrid
        if (dx.ne.dy .or. dx.ne.dz .or. dy.ne.dz) then
          write(OUTDOCK, HALT3) 'Gist step size are inconsitant', 
     &    which_grid
          stop
        endif
        gistparm0%gistspace = dx
        gistparm0%orgin(1) = Xorg
        gistparm0%orgin(2) = Yorg
        gistparm0%orgin(3) = Zorg

        allocate(gistgrid(xngrid, yngrid, zngrid),
     &      stat=istat)
        if (istat .ne. 0) then
          write(OUTDOCK, HALT3) 'Failed to allocate gistgrid', 
     &    which_grid
          stop
        endif

        allocate(gistgridarray(xngrid*yngrid*zngrid),
     &      stat=istat)
        if (istat .ne. 0) then
          write(OUTDOCK, HALT3) 'Failed to allocate gistgrid', 
     &    which_grid
          stop
        endif

        read(GISTFILE,*) gistgridarray
        close(GISTFILE)

        arin = 1 ! array index
        do nx = 1, xngrid
          do ny = 1, yngrid
            do nz = 1, zngrid
               gistgrid(nx,ny,nz) = gistgridarray(arin)
c              if (( gistgridarray(arin) >  1.0 ).or.
c    &            (  gistgridarray(arin) < -1.0 ))then
               if (( gistgridarray(arin) >  0.2 ).or.
     &            (  gistgridarray(arin) < -0.2 ))then
                  !cx = nx * dx + Xorg
                  !cy = ny * dx + Yorg
                  !cz = nz * dx + Zorg
                  cx = (nx-1) * dx + Xorg ! we should substraced by 1 so that the frist point is at the orgin.
                  cy = (ny-1) * dx + Yorg
                  cz = (nz-1) * dx + Zorg
cc                 write(OUTDOCK,*) nx,", ",ny,", ",nz," :  ", 
cc    &                             cx,", ",cy,", ",cz," :  ", 
cc    &                          gistgrid(nx,ny,nz),
cc    &                          gistgridarray(arin)
cc ATOM      6  HA  LEU A   1      19.042  43.323  53.498  1.00  0.00           H
c                  write(OUTDOCK,PDBFORM)
c     &                            "ATOM  ",1,"O1","", 
c     &                            "HOH","A",1,"", cx,cy,cz,
cc    &                            0.00,0.00,"","O",""
c     &                            0.00,gistgridarray(arin),"","O","" ! put gist in the b-factor column. 
cc    &                            gistgridarray(arin),0.00,"","O","" ! put gist in the occupency column. 
               endif
               arin = arin + 1
            enddo
          enddo
        enddo

        ! calculate the min and max grid points. 
        do coord = 1,3
c         ! the grid limit is the orgin
          theminlimits(coord) = gistparm0%orgin(coord) 
c         ! the max limit is the orgin plus the number of point in 
c         ! each dir time by the spacing between points 40 (grid points) * 1/2 (ang/gridpoint) = 20 Angstroms
c         themaxlimits(coord) = (gistparm0%orgin(coord) + 
c    &                          (real(gistparm0%nsize) * 
c    &                           gistparm0%gistspace))
        enddo


c    ! the max limit is the orgin plus the number of point in 
c    ! each dir time by the spacing between points 40 (grid points) * 1/2 (ang/gridpoint) = 20 Angstroms
c    ! x,y,z size may differ. 
        themaxlimits(1) = (gistparm0%orgin(1) + 
     &                        (real(gistparm0%xnsize) * 
     &                        gistparm0%gistspace))
        themaxlimits(2) = (gistparm0%orgin(2) + 
     &                        (real(gistparm0%ynsize) * 
     &                        gistparm0%gistspace))
        themaxlimits(3) = (gistparm0%orgin(3) + 
     &                        (real(gistparm0%znsize) * 
     &                        gistparm0%gistspace))

        ! Initailize the tagvoxelarray and neighborvoxelarray, these will be used in the scorring
        ! function 
        allocate(tagvoxelarray(xngrid*yngrid*zngrid),
     &        stat=istat)
        allocate(neighborvoxelarray(xngrid*yngrid*zngrid),
     &        stat=istat)
        write(*,*) "Allocate tagvoxelarray :: dim = ",
     &           xngrid*yngrid*zngrid
        do idx = 1, xngrid*yngrid*zngrid
        !do idx = 1, gistparm0%xnsize*gistparm0%ynsize*gistparm0%znsize 
           tagvoxelarray(idx) = .FALSE.
           neighborvoxelarray(idx) = .FALSE.
        enddo
        !deallocate(phigrid)
      else
        !gistparm0%nsize = 0 !only this will actually get checked to see if should use in grid demintion limits
        gistparm0%xnsize = 0 !only this will actually get checked to see if should use in grid demintion limits
        gistparm0%ynsize = 0 !
        gistparm0%znsize = 0 !
        gistparm0%gistspace = 0.0
        gistparm0%orgin(1) = 0.0
        gistparm0%orgin(2) = 0.0
        gistparm0%orgin(3) = 0.0
      endif !the if statement checked to see if the filename was not 0.
!if it was 0, literally do nothing here



      !endif

      return
      end subroutine read_GIST_grid

c  --reads in grids for force field type scoring.          ECM  4/91

c this subroutine reads the delphi/qnifft phimap files for electrostatics
c also reads phimaps for receptor desolvation. uses the same format.
c 
      subroutine read_delphi(phimap_precomp, OUTDOCK, 
     &    phimap0, which_grid, phifilename, 
     &    theminlimits, themaxlimits)

      use errorformats
      use phimaptype

      real, allocatable, dimension(:, :, :, :), intent(out) :: 
     &           phimap_precomp
      integer, intent(in) :: OUTDOCK
      type(phimap), intent(inout) :: phimap0
      integer, intent(in) :: which_grid
      character (len=255), intent(in) :: phifilename
      real, dimension(3), intent(inout) :: theminlimits, themaxlimits

c     temporary delphi grid
      real, allocatable :: phigrid(:, :, :)
      integer istat !temporary status

c     character variables to read in phimap (and thrown away)
      character (len=20) :: toplabel
      character (len=10) :: head
      character (len=60) :: title
      character (len=16) :: botlabel

      integer nx, ny, nz, coord
      real a1, a2, a3, a4, a5, a6, a7, a8
      real rgrid

!first, if '0' is the filename, this means we don't actually want to open a file
! or allocate any of the dynamic arrays, we assume the energy of all atoms is 0
! in all places. this is useful for saving memory instead of reading files of 0s
!note that for many reasons you can't use 0 as the first grid and then use non
!zero as the later grids
      if (trim(phifilename) .ne. '0') then
        allocate(phigrid(phimap0%nsize, phimap0%nsize, phimap0%nsize),
     &      stat=istat)
        if (istat .ne. 0) then
          write(OUTDOCK, HALT3) 'Failed to allocate phimap', which_grid
          stop
        endif
        
        open(DELPHIGRID, file=trim(phifilename), 
     &      status='old', form='unformatted', action='read')
        read(DELPHIGRID) toplabel
        read(DELPHIGRID) head, title 
        read(DELPHIGRID) phigrid  
        read(DELPHIGRID) botlabel
        read(DELPHIGRID) phimap0%phiscale, phimap0%oldmid
        close(DELPHIGRID)

c     dynamically allocate phimap_precomp
        allocate(phimap_precomp(8, phimap0%nsize, phimap0%nsize, 
     &       phimap0%nsize), stat=istat)
        if (istat .ne. 0) then
          write(OUTDOCK, HALT3)
     &        'Failed to allocate phimap_precomp ', which_grid
          stop
        endif

c     Populate precomputed version of phimap table {MC}
c     must do this after phiscale and oldmid have been read in
c     11/9/2009
c     Michael Carchia
c     Restructured subtraction expression to limit possibility of 
c       overflow/underflow and thus excessively large/small final values.
c
c     If problem continues, one can explore in this order:
c       1. change the declaration of a1 to double
c       2. Change a1-a8 to double
c       3. compile this file without --fastsse (probably just --fast instead)
c            If one tries #2, it should be compared against #3 for performance 
c            to see which is faster.

        do nz = 1, phimap0%nsize-1
          do ny = 1, phimap0%nsize-1
            do nx = 1, phimap0%nsize-1 ! faster memory read order?
              a8 = phigrid(nx,ny,nz)
              a7 = phigrid(nx,ny,nz+1) - a8
              a6 = phigrid(nx,ny+1,nz) - a8
              a5 = phigrid(nx+1,ny,nz) - a8
              a4 = phigrid(nx,ny+1,nz+1) - a7 - phigrid(nx,ny+1,nz)
              a3 = phigrid(nx+1,ny,nz+1) - a7 - phigrid(nx+1,ny,nz)
              a2 = phigrid(nx+1,ny+1,nz) - a6 - phigrid(nx+1,ny,nz)
              a1 = phigrid(nx+1,ny+1,nz+1) - a4 - phigrid(nx+1,ny,nz+1) 
     &                + phigrid(nx+1,ny,nz) - phigrid(nx+1,ny+1,nz) 
              phimap_precomp(1, nx, ny, nz) = a1
              phimap_precomp(2, nx, ny, nz) = a2
              phimap_precomp(3, nx, ny, nz) = a3
              phimap_precomp(4, nx, ny, nz) = a4
              phimap_precomp(5, nx, ny, nz) = a5
              phimap_precomp(6, nx, ny, nz) = a6
              phimap_precomp(7, nx, ny, nz) = a7
              phimap_precomp(8, nx, ny, nz) = a8
            enddo
          enddo
        enddo


        rgrid = real(phimap0%nsize) - 1.0 !set for later during outside grids
        phimap0%goff = (phimap0%nsize + 1.)/2.
        do coord = 1,3
          theminlimits(coord) = ((1.-phimap0%goff)/phimap0%phiscale) + 
     &        phimap0%oldmid(coord)
          themaxlimits(coord) = ((rgrid-phimap0%goff)/phimap0%phiscale)
     &        + phimap0%oldmid(coord)
        enddo

        deallocate(phigrid)
      else
        phimap0%nsize = 0 !only this will actually get checked to see if 0
        phimap0%phiscale = 0.0
        phimap0%oldmid(1) = 0.0
        phimap0%oldmid(2) = 0.0
        phimap0%oldmid(3) = 0.0
        phimap0%goff = 0.0
      endif !the if statement checked to see if the filename was not 0.
!if it was 0, literally do nothing here

      return
      end subroutine read_delphi

c  --reads in grids for force field type scoring.          ECM  4/91
c additionally reads in grids for bump-map now (chem.bmp)
      subroutine read_vdw(options0, abval_precomp, OUTDOCK,
     &    vdwmap0, which_grid)

      use vdwconsts     !aval_index, bval_index
      use errorformats
      use optionstype
      use vdwmaptype

      type(options), intent(in) :: options0
      real, allocatable, dimension(:, :, :), intent(out) :: 
     &        abval_precomp     
c outdock file specifier, so data can be written out if there are problems
      integer, intent(in) :: OUTDOCK
c vdw grid parameters
      type(vdwmap), intent(inout) :: vdwmap0
      integer, intent(in) :: which_grid !which grids to read

      real, allocatable :: aval(:), bval(:)
c       aval(), bval() --values stored "at" grid points
c format statements for reading vdw grid
      character (len=17), parameter :: vdw1 = "(A17)"
      character (len=45), parameter :: vdw2 = "(4F8.3, 3I4)"
c other temporary variables
      character (len=17) header
      integer istat
      real a1, a2, a3, a4, a5, a6, a7, a8
      integer counti, countj, countk !aval and grid counters
      integer i1, i2, i3, i4, i5, i6, i7, i8

c this reads in the header of the bump map, used for the chemscore grid sizes
      open(unit=BUMPGRID, file=options0%bmpnames(which_grid), 
     &    status='old')
      read(BUMPGRID, vdw1) header
      read(BUMPGRID, vdw2) vdwmap0%grddiv, 
     &    (vdwmap0%offset(counti), counti = 1, 3), 
     &    (vdwmap0%grdpts(counti), counti = 1, 3)
      vdwmap0%ngrd = vdwmap0%grdpts(1) * vdwmap0%grdpts(2) * 
     &    vdwmap0%grdpts(3)
      close(BUMPGRID)

c     Compute 1/grddiv so it can be used in interp_opt() as multiplication.
c     This is faster than doing the division in interp_opt(), but slightly 
c     less accurate.  Speeds up code ~5%. {MC}
      vdwmap0%invgrid = 1.0 / vdwmap0%grddiv

      allocate(aval(vdwmap0%ngrd), stat=istat)
      if (istat .ne. 0) then
         write(OUTDOCK, HALT3)
     &           'Failed to allocate aval ', which_grid
         stop
      endif
      allocate(bval(vdwmap0%ngrd), stat=istat)
      if (istat .ne. 0) then
         write(OUTDOCK, HALT3)
     &           'Failed to allocate bval ', which_grid
         stop
      endif

c      VDW grid
      call binopen(VDWGRID, options0%vdwnames(which_grid))
      read(VDWGRID) (aval(counti), counti = 1, vdwmap0%ngrd)
      read(VDWGRID) (bval(counti), counti = 1, vdwmap0%ngrd)
      close(VDWGRID)

c     Dynamically allocate abval_precomp array
      allocate(abval_precomp(2, 8, vdwmap0%ngrd), stat=istat)
      if (istat .ne. 0) then
         write(OUTDOCK, HALT3) 
     &        'Failed to allocate abval_precomp ', which_grid
         stop
      endif
c     precompute values for aval, bval  optimized arrays {MC}

c  Optimization to precompute elements of aval, bval arrays
c the maximum values on the grid are rounded so that when interpolated they 
c should never overflow to a negative number. this is controlled by the toobig
c parameter below.
c also computes the grid xyz limits. 

      do countk = 1, vdwmap0%grdpts(3)-1
         do countj = 1, vdwmap0%grdpts(2)-1
           do counti = 1, vdwmap0%grdpts(1)-1 ! faster memory read order.
c            get 1-dimensional indices for the points of interest
             i8 = vdwmap0%grdpts(1)*vdwmap0%grdpts(2)*(countk-1) + 
     &            vdwmap0%grdpts(1)*(countj-1) + counti
             i7 = vdwmap0%grdpts(1)*vdwmap0%grdpts(2)*countk + 
     &            vdwmap0%grdpts(1)*(countj-1) + counti
             i6 = vdwmap0%grdpts(1)*vdwmap0%grdpts(2)*(countk-1) + 
     &            vdwmap0%grdpts(1)*countj + counti
             i5 = vdwmap0%grdpts(1)*vdwmap0%grdpts(2)*(countk-1) + 
     &            vdwmap0%grdpts(1)*(countj-1) + counti+1
             i4 = vdwmap0%grdpts(1)*vdwmap0%grdpts(2)*countk + 
     &            vdwmap0%grdpts(1)*countj + counti
             i3 = vdwmap0%grdpts(1)*vdwmap0%grdpts(2)*countk + 
     &            vdwmap0%grdpts(1)*(countj-1) + counti+1
             i2 = vdwmap0%grdpts(1)*vdwmap0%grdpts(2)*(countk-1) + 
     &            vdwmap0%grdpts(1)*countj + counti+1
             i1 = vdwmap0%grdpts(1)*vdwmap0%grdpts(2)*countk + 
     &            vdwmap0%grdpts(1)*countj + counti+1
c adjust the terms down if they are too big
             if (aval(i8) .gt. toobig) then
               aval(i8) = toobig
             endif               
             if (aval(i7) .gt. toobig) then
               aval(i7) = toobig
             endif               
             if (aval(i6) .gt. toobig) then
               aval(i6) = toobig
             endif               
             if (aval(i5) .gt. toobig) then
               aval(i5) = toobig
             endif               
             if (aval(i4) .gt. toobig) then
               aval(i4) = toobig
             endif               
             if (aval(i3) .gt. toobig) then
               aval(i3) = toobig
             endif               
             if (aval(i2) .gt. toobig) then
               aval(i2) = toobig
             endif               
             if (aval(i1) .gt. toobig) then
               aval(i1) = toobig
             endif               
c            calculate coefficients of trilinear function for aval
             a8 = aval(i8)
             a7 = aval(i7) - a8
             a6 = aval(i6) - a8
             a5 = aval(i5) - a8
             a4 = aval(i4) - a7 - aval(i6) 
             a3 = aval(i3) - a7 - aval(i5) 
             a2 = aval(i2) - a6 - aval(i5) 
             a1 = aval(i1) - a4 - aval(i3) + aval(i5) - aval(i2)
             abval_precomp(AVAL_INDEX, 1, i8) = a1
             abval_precomp(AVAL_INDEX, 2, i8) = a2
             abval_precomp(AVAL_INDEX, 3, i8) = a3
             abval_precomp(AVAL_INDEX, 4, i8) = a4
             abval_precomp(AVAL_INDEX, 5, i8) = a5
             abval_precomp(AVAL_INDEX, 6, i8) = a6
             abval_precomp(AVAL_INDEX, 7, i8) = a7
             abval_precomp(AVAL_INDEX, 8, i8) = a8
c     --calculate coefficients of trilinear function for bval
             if (bval(i8) .gt. toobig) then
               bval(i8) = toobig
             endif               
             if (bval(i7) .gt. toobig) then
               bval(i7) = toobig
             endif               
             if (bval(i6) .gt. toobig) then
               bval(i6) = toobig
             endif               
             if (bval(i5) .gt. toobig) then
               bval(i5) = toobig
             endif               
             if (bval(i4) .gt. toobig) then
               bval(i4) = toobig
             endif               
             if (bval(i3) .gt. toobig) then
               bval(i3) = toobig
             endif               
             if (bval(i2) .gt. toobig) then
               bval(i2) = toobig
             endif               
             if (bval(i1) .gt. toobig) then
               bval(i1) = toobig
             endif               
             a8 = bval(i8)
             a7 = bval(i7) - a8
             a6 = bval(i6) - a8
             a5 = bval(i5) - a8
             a4 = bval(i4) - a7 - bval(i6) 
             a3 = bval(i3) - a7 - bval(i5) 
             a2 = bval(i2) - a6 - bval(i5) 
             a1 = bval(i1) - a4 - bval(i3) + bval(i5) - bval(i2)
             abval_precomp(BVAL_INDEX, 1, i8) = a1
             abval_precomp(BVAL_INDEX, 2, i8) = a2
             abval_precomp(BVAL_INDEX, 3, i8) = a3
             abval_precomp(BVAL_INDEX, 4, i8) = a4
             abval_precomp(BVAL_INDEX, 5, i8) = a5
             abval_precomp(BVAL_INDEX, 6, i8) = a6
             abval_precomp(BVAL_INDEX, 7, i8) = a7
             abval_precomp(BVAL_INDEX, 8, i8) = a8
            enddo
         enddo
      enddo
c compute limits of this grid
      vdwminlimits(1) = vdwmap0%offset(1) ! yes really. offset is poorly named
      vdwminlimits(2) = vdwmap0%offset(2) 
      vdwminlimits(3) = vdwmap0%offset(3)
      vdwmaxlimits(1) = vdwmap0%grdpts(1) * vdwmap0%grddiv + 
     &    vdwmap0%offset(1)
      vdwmaxlimits(2) = vdwmap0%grdpts(2) * vdwmap0%grddiv + 
     &    vdwmap0%offset(2)
      vdwmaxlimits(3) = vdwmap0%grdpts(3) * vdwmap0%grddiv + 
     &    vdwmap0%offset(3)
c deallocate the read in data grids.
      deallocate(aval)
      deallocate(bval)
      return
      end subroutine read_vdw

************************************************************************
*                subroutine readsol                                    * 
c this subroutine reads in the grid of points defined by the volume of *
c a receptor as created by the program solvmap into the array solgrd,  *
c which is used in dock to decide the desolvation penalty of each atom *
c of a given orientation.                                              *
c BKS 2003 and MMM 2008-2009                                           *
c***********************************************************************
      subroutine readsol(solgrd_precomp,
     &    hsolgrd_precomp, OUTDOCK,
     &    solvmap0, options0, which_grid)

      use errorformats
      use solvmaptype
      use optionstype

      real, allocatable, dimension(:, :, :, :), intent(out) :: 
     &           solgrd_precomp
      real, allocatable, dimension(:, :, :, :), intent(out) :: 
     &           hsolgrd_precomp
      integer, intent(in) :: OUTDOCK
      type(solvmap), intent(inout) :: solvmap0
      type(options), intent(inout) :: options0
      integer, intent(in) :: which_grid

c       solgrd: % desolvated pts in space are, based on volume occluded by
c          the receptor
c       hsolgrd: additional % desolvation grid for hydrogen atoms
      real, allocatable, dimension(:, :, :) :: solgrd
      real, allocatable, dimension(:, :, :) :: hsolgrd

      integer      i, j, k, kk, istat !replace these eventually
      character (len=80) line
      integer      hsperang
      real         hsmax(3)
      integer      hscadif(3) !dimensions of hydrogen solvation grid
      integer nx,ny,nz
      real    a1,a2,a3,a4,a5,a6,a7,a8
c format statements for solvmap reading
      character (len=51), parameter :: solhead = '(4(i4,2x),3(f8.3,1x))'
      character (len=79), parameter :: soldata = '(13f6.3)'

c      first try reading integer edge locations for backwards compatibility
      open(SOLGRIDU, file=options0%solvnames(which_grid), status='old')
      read(SOLGRIDU, '(A80)') line
      read(line, solhead) solvmap0%scadif(1), solvmap0%scadif(2), 
     &    solvmap0%scadif(3), solvmap0%sperang, solvmap0%smax(1),
     &    solvmap0%smax(2), solvmap0%smax(3)
      allocate(solgrd(0:solvmap0%scadif(1), 0:solvmap0%scadif(2),
     &    0:solvmap0%scadif(3)), stat=istat)
      if (istat .ne. 0) then
        write(OUTDOCK, HALT3) 'Failed to allocate solgrd ', which_grid
        stop
      endif
      do i = 0, solvmap0%scadif(1)
        do j = 0, solvmap0%scadif(2)
          do k = 0, solvmap0%scadif(3), 13
            if (k + 12 .le. solvmap0%scadif(3)) then !lines have 12 data chunks
              read(SOLGRIDU, soldata) (solgrd(i,j,kk), kk=k, k+12)
            else !sometimes lines don't have all 12 chunks
              read(SOLGRIDU, soldata) (solgrd(i,j,kk), kk=k,
     &             solvmap0%scadif(3))
            endif
          enddo
        enddo
      enddo
      close(SOLGRIDU) !closing files
c     Dynamically allocate solgrd_precomp and hsolgrd_precomp
      allocate(solgrd_precomp(8, 0:solvmap0%scadif(1), 
     &    0:solvmap0%scadif(2), 0:solvmap0%scadif(3)), stat=istat)
      if (istat .ne. 0) then
        write(OUTDOCK, HALT3)
     &         'Failed to allocate solgrd_precomp ', which_grid
        stop
      endif
      if (options0%hsolflag) then
        open(HSOLGRIDU, file=options0%hsolvnames(which_grid), 
     &      status='old')
        read(HSOLGRIDU, solhead) hscadif(1), hscadif(2), hscadif(3),
     &      hsperang, hsmax(1), hsmax(2), hsmax(3)
        if (hscadif(1) .ne. solvmap0%scadif(1) .or. 
     &       hscadif(2) .ne. solvmap0%scadif(2) .or.
     &       hscadif(3) .ne. solvmap0%scadif(3) .or. 
     &       hsperang .ne. solvmap0%sperang .or.
     &       hsmax(1) .ne. solvmap0%smax(1) .or. 
     &       hsmax(2) .ne. solvmap0%smax(2) .or. 
     &       hsmax(3) .ne. solvmap0%smax(3)) then
          write(OUTDOCK, HALT1)
     &          'The hydrogen solvmap grid dimensions or '
          write(OUTDOCK, HALT1)
     &          '   parameters do not match the corresponding'
          write(OUTDOCK, HALT1) '   ones in solvmap.'
          write(OUTDOCK, HALT0) 
          stop           
        endif
        allocate(hsolgrd(0:solvmap0%scadif(1), 0:solvmap0%scadif(2),
     &      0:solvmap0%scadif(3)), stat=istat)
        if (istat .ne. 0) then
          write(OUTDOCK, HALT3) 
     &        'Failed to allocate hsolgrd ', which_grid
          stop
        endif
        do i = 0, solvmap0%scadif(1)
          do j = 0, solvmap0%scadif(2)
            do k = 0, solvmap0%scadif(3), 13
              if (k+12 .le. solvmap0%scadif(3)) then !same as before, 12 chunks per line
                read(HSOLGRIDU, soldata) (hsolgrd(i, j, kk), 
     &              kk = k, k + 12)
              else
                read(HSOLGRIDU, soldata) (hsolgrd(i, j, kk),
     &              kk = k, solvmap0%scadif(3))
              endif
            enddo
          enddo
        enddo
        close(HSOLGRIDU)
        allocate(hsolgrd_precomp(8, 0:solvmap0%scadif(1), 
     &      0:solvmap0%scadif(2), 0:solvmap0%scadif(3)), stat=istat)
        if (istat .ne. 0) then
           write(OUTDOCK, HALT3) 
     &          'Failed to allocate hsolgrd_precomp', which_grid
          stop
        endif
      endif
c     Precalculate the values from the solgrd and hsolgrd to be used by transfm
c     Major speed optimization.{MC}
      if (options0%hsolflag) then
        do nz = 0, solvmap0%scadif(3)-1
          do ny = 0, solvmap0%scadif(2)-1
            do nx = 0, solvmap0%scadif(1)-1
              a8 = hsolgrd(nx,ny,nz)
              a7 = hsolgrd(nx,ny,nz+1) - a8
              a6 = hsolgrd(nx,ny+1,nz) - a8
              a5 = hsolgrd(nx+1,ny,nz) - a8
              a4 = hsolgrd(nx,ny+1,nz+1) - a7 - hsolgrd(nx,ny+1,nz)
              a3 = hsolgrd(nx+1,ny,nz+1) - a7 - hsolgrd(nx+1,ny,nz)
              a2 = hsolgrd(nx+1,ny+1,nz) - a6 - hsolgrd(nx+1,ny,nz)
              a1 = hsolgrd(nx+1,ny+1,nz+1) - a4 - hsolgrd(nx+1,ny,nz+1)
     &             + hsolgrd(nx+1,ny,nz) - hsolgrd(nx+1,ny+1,nz)
              hsolgrd_precomp(1, nx, ny, nz) = a1
              hsolgrd_precomp(2, nx, ny, nz) = a2
              hsolgrd_precomp(3, nx, ny, nz) = a3
              hsolgrd_precomp(4, nx, ny, nz) = a4
              hsolgrd_precomp(5, nx, ny, nz) = a5
              hsolgrd_precomp(6, nx, ny, nz) = a6
              hsolgrd_precomp(7, nx, ny, nz) = a7
              hsolgrd_precomp(8, nx, ny, nz) = a8
            enddo
          enddo
        enddo
      endif
      do nz = 0, solvmap0%scadif(3)-1
        do ny = 0, solvmap0%scadif(2)-1
          do nx = 0, solvmap0%scadif(1)-1
            a8 = solgrd(nx,ny,nz)
            a7 = solgrd(nx,ny,nz+1) - a8
            a6 = solgrd(nx,ny+1,nz) - a8
            a5 = solgrd(nx+1,ny,nz) - a8
            a4 = solgrd(nx,ny+1,nz+1) - a7 - solgrd(nx,ny+1,nz)
            a3 = solgrd(nx+1,ny,nz+1) - a7 - solgrd(nx+1,ny,nz)
            a2 = solgrd(nx+1,ny+1,nz) - a6 - solgrd(nx+1,ny,nz)
            a1 = solgrd(nx+1,ny+1,nz+1) - a4 - solgrd(nx+1,ny,nz+1)
     &           + solgrd(nx+1,ny,nz) - solgrd(nx+1,ny+1,nz)
            solgrd_precomp(1, nx, ny, nz) = a1
            solgrd_precomp(2, nx, ny, nz) = a2
            solgrd_precomp(3, nx, ny, nz) = a3
            solgrd_precomp(4, nx, ny, nz) = a4
            solgrd_precomp(5, nx, ny, nz) = a5
            solgrd_precomp(6, nx, ny, nz) = a6
            solgrd_precomp(7, nx, ny, nz) = a7
            solgrd_precomp(8, nx, ny, nz) = a8
          enddo
        enddo
      enddo
c for outside_grid checking, precomputed limits
      solmaxlimits(1) = (real(solvmap0%smax(1)))/solvmap0%sperang
      solmaxlimits(2) = (real(solvmap0%smax(2)))/solvmap0%sperang
      solmaxlimits(3) = (real(solvmap0%smax(3)))/solvmap0%sperang
      solminlimits(1) = (real(solvmap0%smax(1)) - 
     &    real(solvmap0%scadif(1)))/solvmap0%sperang
      solminlimits(2) = (real(solvmap0%smax(2)) - 
     &    real(solvmap0%scadif(2)))/solvmap0%sperang
      solminlimits(3) = (real(solvmap0%smax(3)) - 
     &    real(solvmap0%scadif(3)))/solvmap0%sperang
      deallocate(solgrd)
      if (options0%hsolflag) then !single line if statements are evil
        deallocate(hsolgrd)
      endif
      return
      end subroutine readsol

c precomputes the intersection of all xyz limits so outside_grids.f
c can just do 1 check instead of several. saves lots of time.
      subroutine og_precomputelimits(allminlimits, allmaxlimits,
c     &    ldesolv, receptor_desolv, nsize, gistxnsize)
     &    ldesolv, receptor_desolv, nsize, gistflag, rec_des_r)

c grid min & max limits
      real, dimension(3), intent(out) :: allminlimits
      real, dimension(3), intent(out) :: allmaxlimits
c various parameters
      integer, intent(in) :: ldesolv !which desolvation method to use
      logical, intent(in) :: receptor_desolv, rec_des_r
      integer, intent(in) :: nsize !phimap nsize, if 0 ignore phi???limits
c      integer, intent(in) :: gistxnsize !gistparm xnsize, if 0 ignore gistlimits, this means the gist grid is not set!
      logical, intent(in) :: gistflag

      integer coord !temporary variable used to walk through coords


      do coord = 1,3
        allminlimits(coord) = vdwminlimits(coord)
        if ((ldesolv .eq. 2) .or. (ldesolv .eq. 3)) then !if solvmap was loaded
          if (solminlimits(coord) .gt. allminlimits(coord)) then
            allminlimits(coord) = solminlimits(coord)
          endif
        endif
        if (receptor_desolv) then !if receptor desolvation
          if (recminlimits(coord) .gt. allminlimits(coord)) then
            allminlimits(coord) = recminlimits(coord)
          endif
        endif
        if (nsize .ne. 0) then
          if (phiminlimits(coord) .gt. allminlimits(coord)) then
            allminlimits(coord) = phiminlimits(coord)
          endif
        endif
        !if (gistxnsize .ne. 0) then
        if (gistflag) then
          if (gistminlimits(coord) .gt. allminlimits(coord)) then
            allminlimits(coord) = gistminlimits(coord)
          endif
        endif
        if (rec_des_r) then
          if (rec_d_rminlimits(coord) .gt. allminlimits(coord)) then
            allminlimits(coord) = rec_d_rminlimits(coord)
          endif
        endif
c       now the max limit
        allmaxlimits(coord) = vdwmaxlimits(coord)
        if ((ldesolv .eq. 2) .or. (ldesolv .eq. 3)) then !if solvmap was loaded
          if (solmaxlimits(coord) .lt. allmaxlimits(coord)) then
            allmaxlimits(coord) = solmaxlimits(coord)
          endif
        endif
        if (receptor_desolv) then !if receptor desolvation used
          if (recmaxlimits(coord) .lt. allmaxlimits(coord)) then
            allmaxlimits(coord) = recmaxlimits(coord)
          endif
        endif
        if (nsize .ne. 0) then
          if (phimaxlimits(coord) .lt. allmaxlimits(coord)) then
            allmaxlimits(coord) = phimaxlimits(coord)
          endif
        endif
        !if (gistxnsize .ne. 0) then
        if (gistflag) then
          if (gistmaxlimits(coord) .lt. allmaxlimits(coord)) then
            allmaxlimits(coord) = gistmaxlimits(coord)
          endif
        endif
        if (rec_des_r) then
          !write(*,*) "I AM HERE"
          if (rec_d_rmaxlimits(coord) .lt. allmaxlimits(coord)) then
            allmaxlimits(coord) = rec_d_rmaxlimits(coord)
          endif
        endif
      enddo

      return
      end subroutine og_precomputelimits
 
! this subroutine checks to make sure the read in grids are the same size
! this must be done for flexible receptors (for now)
! phimap grids are allowed to be size 0 and implicitly 0 everywhere.
!
! **** this might have to be modified to check gist grid.
      function gridsize_check(vdwmap0, vdwmap1, 
     &    phimap0, phimap1, recdes0, recdes1, solvmap0, solvmap1, 
     &    ligdes0, ligdes1)
     &    result(sizescheck)

      use vdwmaptype
      use phimaptype
      use solvmaptype
 
      type(vdwmap), intent(in) :: vdwmap0, vdwmap1
      type(phimap), intent(in) :: phimap0, phimap1, recdes0, recdes1
      type(solvmap), intent(in) :: solvmap0, solvmap1
      type(phimap), intent(in) :: ligdes0, ligdes1
      ! XXX: GFortran complains this is not a DUMMY
      !logical, intent(out) :: sizescheck !return value
      logical :: sizescheck !return value

      integer :: count !variable to run through multidimensional things

      sizescheck = .true. !default, if any check fails change to false
      if (vdwmap0%ngrd .ne. vdwmap1%ngrd) then
        sizescheck = .false.
      endif
      if (vdwmap0%grddiv .ne. vdwmap1%grddiv) then
        sizescheck = .false.
      endif
      if (vdwmap0%invgrid .ne. vdwmap1%invgrid) then
        sizescheck = .false.
      endif
      if (phimap0%phiscale .ne. phimap1%phiscale) then
        !next line checks for 0s in either grid, which makes things okay
        if ((phimap0%nsize .ne. 0) .and. (phimap1%nsize .ne. 0)) then
          sizescheck = .false.
        endif
      endif
      if (phimap0%goff .ne. phimap1%goff) then
        !next line checks for 0s in either grid, which makes things okay
        if ((phimap0%nsize .ne. 0) .and. (phimap1%nsize .ne. 0)) then
          sizescheck = .false.
        endif
      endif
      if (phimap0%nsize .ne. phimap1%nsize) then
        !next line checks for 0s in either grid, which makes things okay
        if ((phimap0%nsize .ne. 0) .and. (phimap1%nsize .ne. 0)) then
          sizescheck = .false.
        endif
      endif
      if (recdes0%phiscale .ne. recdes1%phiscale) then
        !next line checks for 0s in either grid, which makes things okay
        if ((recdes0%nsize .ne. 0) .and. (recdes1%nsize .ne. 0)) then
          sizescheck = .false.
        endif
      endif
      if (recdes0%goff .ne. recdes1%goff) then
        !next line checks for 0s in either grid, which makes things okay
        if ((recdes0%nsize .ne. 0) .and. (recdes1%nsize .ne. 0)) then
          sizescheck = .false.
        endif
      endif
      if (recdes0%nsize .ne. recdes1%nsize) then
        !next line checks for 0s in either grid, which makes things okay
        if ((recdes0%nsize .ne. 0) .and. (recdes1%nsize .ne. 0)) then
          sizescheck = .false.
        endif
      endif
      if (ligdes0%phiscale .ne. ligdes1%phiscale) then
        !next line checks for 0s in either grid, which makes things okay
        if ((ligdes0%nsize .ne. 0) .and. (ligdes1%nsize .ne. 0)) then
          sizescheck = .false.
        endif
      endif
      if (ligdes0%goff .ne. ligdes1%goff) then
        !next line checks for 0s in either grid, which makes things okay
        if ((ligdes0%nsize .ne. 0) .and. (ligdes1%nsize .ne. 0)) then
          sizescheck = .false.
        endif
      endif
      if (ligdes0%nsize .ne. ligdes1%nsize) then
        !next line checks for 0s in either grid, which makes things okay
        if ((ligdes0%nsize .ne. 0) .and. (ligdes1%nsize .ne. 0)) then
          sizescheck = .false.
        endif
      endif
      if (solvmap0%sperang .ne. solvmap1%sperang) then
        sizescheck = .false.
      endif
      do count = 1, 3
        if (vdwmap0%offset(count) .ne. vdwmap1%offset(count)) then
          sizescheck = .false.
        endif
        if (vdwmap0%grdpts(count) .ne. vdwmap1%grdpts(count)) then
          sizescheck = .false.
        endif
        if (phimap0%oldmid(count) .ne. phimap1%oldmid(count)) then
          !next line checks for 0s in either grid, which makes things okay
          if ((phimap0%nsize .ne. 0) .and. (phimap1%nsize .ne. 0)) then
            sizescheck = .false.
          endif
        endif
        if (recdes0%oldmid(count) .ne. recdes1%oldmid(count)) then
          !next line checks for 0s in either grid, which makes things okay
          if ((recdes0%nsize .ne. 0) .and. (recdes1%nsize .ne. 0)) then
            sizescheck = .false.
          endif
        endif
        if (ligdes0%oldmid(count) .ne. ligdes1%oldmid(count)) then
          !next line checks for 0s in either grid, which makes things okay
          if ((ligdes0%nsize .ne. 0) .and. (ligdes1%nsize .ne. 0)) then
            sizescheck = .false.
          endif
        endif
        if (solvmap0%scadif(count) .ne. solvmap1%scadif(count)) then
          sizescheck = .false.
        endif
        if (solvmap0%smax(count) .ne. solvmap1%smax(count)) then
          sizescheck = .false.
        endif
      enddo

      return
      end function gridsize_check

      end module grid_reader
