*  
* Mysinger 12/12/2007
*
      SUBROUTINE GZSEEK (HANDLE, OFFSET, ISTAT)
************************************************************************
*
*  Seek to absolute offset within file.
*
*  Arguments:
*     HANDLE     file handle of gzip file
*     OFFSET     absolute offset in file to go to
*     ISTAT      set to absolute offset in file if successful
************************************************************************
*
      IMPLICIT NONE
      INTEGER*8 HANDLE
      INTEGER  ISTAT, OFFSET
*
      CALL GZSEEKER(HANDLE, OFFSET, ISTAT)
      RETURN
      END
*
