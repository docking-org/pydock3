*  
* $Id: gzbread.F,v 1.1.1.1 2007-12-12 21:55:41 mysinger Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2002/05/03 16:11:54  hvogt
* Imported sources
*
*  
      SUBROUTINE GZBREAD(HANDLE,CHBUF, NBBUF, ISTAT)
************************************************************************
*
*  Reads a line from a zipped ascii file
*
*  Arguments:
*     HANDLE     file handle of gzip file
*     CHBUF      array for storing unzipped info
*     NBBUF      record size, number of bytes of array CHBUF
*     ISTAT      status, =NBREC if ok., =zero if EOF encountered
************************************************************************
*
      IMPLICIT NONE
      INTEGER*8 HANDLE
*
      CHARACTER(*)   CHBUF
      INTEGER         NBBUF, ISTAT
*
      CALL GZGETB (HANDLE,NBBUF,CHBUF,ISTAT)
*
      END
*
