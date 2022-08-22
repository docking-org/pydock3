*  
* $Id: gzclose.F,v 1.1.1.1 2007-12-12 21:55:41 mysinger Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2002/05/03 16:11:54  hvogt
* Imported sources
*
*  
      SUBROUTINE GZCLOSE (HANDLE,istat)
*
*
************************************************************************
*
*  Closes the zipped file
*
*  Argument:
*     HANDLE     file handle of gzip file
*     MODE       string selecting IO mode (r, w)
*     ISTAT      status, =zero if success, #0 otherwise
*
************************************************************************
*
      IMPLICIT NONE
      INTEGER*8 HANDLE
      integer istat

      CALL GZCLOS (HANDLE,ISTAT)
      IF (ISTAT .NE. 0)  THEN
          WRITE (6,'(A)')
     &          'GZCLOSE: Error - Input file could not be closed'
      END IF
*
      END
*
