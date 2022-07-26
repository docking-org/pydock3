*  
* $Id: gzwrite.F,v 1.1.1.1 2007-12-12 21:55:41 mysinger Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2002/05/03 16:11:54  hvogt
* Imported sources
*
*  
      SUBROUTINE GZWRITE (HANDLE,OUTLINE,ISTAT)
*
*
************************************************************************
*
*  Writes a line to a zipped file
*
*  Argument:
*     HANDLE     file handle of gzip file
*     OUTLINE    Character string to write
*     ISTAT      status, =zero if success, #0 otherwise
*
*
************************************************************************
*
      IMPLICIT NONE
      INTEGER*8 HANDLE

      CHARACTER(*)   OUTLINE
      INTEGER         LGN, ISTAT, LNBLNK
*
      LGN = LNBLNK (OUTLINE)
      IF (OUTLINE(LGN:LGN).ne.char(10)) then
        LGN = LGN + 1
        OUTLINE(LGN:LGN) = CHAR(10)
      endif
      CALL GZPUTS (HANDLE,LGN,OUTLINE,ISTAT)
      IF (ISTAT .LE. 0)  THEN
        WRITE (6,'(A)')
     &          'GZWRITE: Error while writing to the zipped file'
        STOP
      END IF
*
      END
*
