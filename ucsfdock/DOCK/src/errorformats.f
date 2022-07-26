! module contains error formats used in multiple places.
      module errorformats

      implicit none

c XXX: Needed to add commas after "ERROR ---> " for GFortran compatibility (TS)

C These are the error format statements used throughout DOCK
C Include this piece of source to make them accessible to a given routine
C DAG 8/94
! RGC 2012 - take out of .h file and put in a module, modern f95

C Cautions - no big deal

      character (len=*), parameter :: CAUTION0 = '("CAUTION ---> ",A)'        
      character (len=*), parameter :: CAUTION1 = '("CAUTION ---> ",2A)'       

C Warnings - of significant interest

      character (len=*), parameter :: WARNING0 = '("WARNING ---> ",A)'        
      character (len=*), parameter :: WARNING1 = '("WARNING ---> ",2A)'       
      character (len=*), parameter :: WARNING2 = 
     &      '("WARNING ---> ",A,I6)'     

C Serious errors - cause program to halt

      character (len=*), parameter :: HALT0 = '("ERROR ---> stops")'
      character (len=*), parameter :: HALT1 = '("ERROR ---> ",A)'       
      character (len=*), parameter :: HALT2 = '("ERROR ---> ",2A)'            
      character (len=*), parameter :: HALT3 = '("ERROR ---> ",A,I6)'          
      character (len=*), parameter :: HALT4 = '("ERROR ---> ",A,I3,A)'        
      character (len=*), parameter :: HALT5 = 
     &      '("ERROR ---> ",A,I3,A,I6)'     

      end module errorformats
