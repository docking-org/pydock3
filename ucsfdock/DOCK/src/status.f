c this file contains a fake enumerated type to track the status of the 
c docking progress throughout the search/match/hierarchy/minimization 
c procedure. 
      module status

      implicit none
      integer, parameter :: MIXMATCHED = -3
      integer, parameter :: CLASHES = -2
      integer, parameter :: STRAIN = -4
      integer, parameter :: NOTCHECKED = -1
      integer, parameter :: NOMATCH = 0
      integer, parameter :: OUTSIDEGRIDS = 1
      integer, parameter :: BUMPED = 2
      integer, parameter :: TIMEDOUT = 3!actually not going to use this
      integer, parameter :: ALLOKAY = 4
      integer, parameter :: RIGIDOKAY = 5
 

      end module status
