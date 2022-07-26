*  
* Mysinger 12/12/2007
*
*  
      SUBROUTINE GZTELL (HANDLE, LOC)
************************************************************************
*
*  Get location in file.
*
*  Arguments:
*     HANDLE     file handle of gzip file
*     LOC        return - location in file (absolute)
************************************************************************
*
      IMPLICIT NONE
      INTEGER*8 HANDLE
      INTEGER LOC
*
      CALL GZTELLER(HANDLE, LOC)
      RETURN
      END
*
