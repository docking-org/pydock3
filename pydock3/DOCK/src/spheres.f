!module contains data about spheres & color/chemical matching of them
      module spheres
   
      implicit none

      integer, parameter :: MAXPTS = 100 !ugly ugly ugly, max sphere pts
      integer, parameter :: MAXCOL = 20 !more ugliness, max number of colors
      !integer, parameter :: MAXCOL = 12 !more ugliness, max number of colors
      !integer, parameter :: MAXCOLMATCH = 26 !more ugliness, max number of color matches
      integer, parameter :: MAXCOLMATCH = 40 !more ugliness, max number of color matches
      type spherest
        integer :: nsphr !total number of spheres in RECEPTOR, an important distinction
        real, dimension(3, MAXPTS) :: spcorr !array of sphere center coordinates for the cluster
        integer, dimension(MAXPTS) :: scolor !sphere color
        real, dimension(MAXPTS, MAXPTS) :: disl !ligand sph-sph distances
        real, dimension(MAXPTS, MAXPTS) :: disr !receptor sph-sph distances
        integer :: numcol !number of colors
        integer :: nligcl, nreccl !number of ligand, receptor colors
        integer :: nlrmat !number of ligand-receptor matches
        integer :: colrej !number of matches rejected because of color mismatch
        integer :: coltest ! number of matches tested for color mismatch
        integer, dimension(MAXCOL, MAXCOL) :: lgspcl !ligand-sphere color matching table
        integer, dimension(MAXCOL) :: rectrn !receptor translation tables (for color numbers)
        character (len=30), dimension(MAXCOL) :: lclnam !ligand color name
        character (len=30), dimension(MAXCOL) :: rclnam !rec. color name
        character (len=30), dimension(MAXCOL) :: colnam !color name
        character (len=30), dimension(2, MAXCOLMATCH) :: premat !preliminary(input) ligand-sphere matching colors
      end type spherest

      end module spheres



