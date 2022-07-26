*  
* $Id: gzbwrite.F,v 1.1.1.1 2007-12-12 21:55:41 mysinger Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2002/05/03 16:11:54  hvogt
* Imported sources
*
*  
      SUBROUTINE GZBWRITE (HANDLE,BUFFER,NBYTES, ISTAT)
*
*
************************************************************************
*
*  Writes the binary content in BUFFER to a zipped file
*
*  Arguments:
*     HANDLE     file handle of gzip file
*     BUFFER     Array of binary content
*     NBYTES     Number of Bytes to write
*     ISTAT      status, =zero if success, #0 otherwise
*
*
************************************************************************
*
      IMPLICIT NONE
      INTEGER*8 HANDLE

      INTEGER         BUFFER(*)
      INTEGER         NBYTES, ISTAT
*
      CALL GZPUTS (HANDLE,NBYTES,BUFFER,ISTAT)
      IF (ISTAT .LT. 0)  THEN
        WRITE (6,'(A)')
     &          'GZBWRITE: Error while writing to the zipped file'
        STOP
      END IF
*
      END
*
