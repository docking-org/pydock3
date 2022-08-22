*  
* Mysinger 12/12/2007
*  
      SUBROUTINE GZFLUSH (HANDLE, ISTAT)
************************************************************************
*
*  Flush output to zipped file.
*
*  Arguments:
*     HANDLE     file handle of gzip file
*     ISTAT      set to zero if successful
************************************************************************
*
      IMPLICIT NONE
      INTEGER*8 HANDLE
      INTEGER   ISTAT
*
      CALL GZFLUSHER(HANDLE, ISTAT)
      RETURN
      END
*
