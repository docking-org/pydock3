c this file contains a bunch of enumerated types (not really, this is fortran)
c the purpose of these is to make the file i/o more understandable
c by assigning each number a name that makes at least some sense
c write(6, *) will get replaced by write(OUTDOCK, *) which at least has an 
c opportunity to be understood
      module filenums

      implicit none

      integer, parameter :: MATCHSPH = 1 ! the matching sphere file for the receptor, has colors
      integer, parameter :: INDOCK  = 51! the INDOCK file, keywords & parameters to use for docking
      integer, parameter :: OUTDOCK = 6! the outdock file
      integer, parameter :: OUTDOCK_premin = 7! the premin outdock file 
      !integer, parameter :: COMR = 9! COMR file
      !integer, parameter :: COML = 10! COML file 
      !integer, parameter :: ROT = 11! ROT file 
      integer, parameter :: OUTPDB = 8! the debugging pdb file 
      integer, parameter :: VDWTYPESFILE = 70 !the file containing vdw types
      integer, parameter :: SDIFILE = 21 !split database index file
      integer, parameter :: RESTARTFILE = 22 ! restart data file

      end module filenums
