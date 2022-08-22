*  
* $Id: gzopen.F,v 1.1.1.1 2007-12-12 21:55:42 mysinger Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2002/05/03 16:11:54  hvogt
* Imported sources
*
*  
      SUBROUTINE GZOPEN (HANDLE, MODE, FILNAM, ISTAT)
*
************************************************************************
*
*  Opens a file for zipped IO
*
*  Arguments:
*     HANDLE     file handle of gzip file
*     MODE       string selecting IO mode (r, w)
*     FILNAM     name of the file (including the null-terminator)
*     ISTAT      status, =zero if success, #0 otherwise
*
*  Input  file: Extension '.gz' added to the file name
*  Output file: Extension '.gz' is checked
*
************************************************************************
*
      IMPLICIT NONE
      INTEGER*8 HANDLE

      CHARACTER(1)   MODE, MYMODE
      CHARACTER(*) FILNAM
      CHARACTER*256 MYFNAM
      INTEGER       FLNLEN, FLNLEN2, ISTAT
*

      ISTAT=0

      MYFNAM = FILNAM
      FLNLEN = INDEX(MYFNAM,' ')
      IF (FLNLEN .EQ. 0) THEN
        FLNLEN = LEN(MYFNAM)
        FLNLEN2 = FLNLEN
      ELSE
        MYFNAM(FLNLEN:FLNLEN) = CHAR(0)
        FLNLEN2 = FLNLEN-1
      END IF
      IF (MODE .EQ. 'r' .OR. MODE .EQ. 'R')  MYMODE='r'
      IF (MODE .EQ. 'w' .OR. MODE .EQ. 'W')  MYMODE='w'
*
      IF (MYMODE .EQ. 'r')  THEN
c        IF (MYFNAM(FLNLEN2-2:FLNLEN2) .EQ. '.gz')  THEN
c rgc changed this so we could read gzipped streams from standard in
c who knows what other effects this may have
          CALL GZIOOP (HANDLE, 'r', MYFNAM, FLNLEN, ISTAT)
          GO  TO  90
c        ELSE
c          WRITE (6,'(A)')
c     &          'GZOPEN: Error - Input file name has no .gz extension'
c          ISTAT = 1
c          GO  TO  90
c        END IF
      ELSE IF (MYMODE .EQ. 'w')  THEN
        CALL GZIOOP (HANDLE,'w',MYFNAM,FLNLEN,ISTAT)
        GO  TO  90
      ELSE
        WRITE (6,'(A)')
     &          'GZOPEN: Error - IO mode is not R or W'
        ISTAT = 1
      END IF
*
   90 RETURN
*
      END
*
