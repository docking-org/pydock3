*  
* $Id: gzread.F,v 1.1.1.1 2007-12-12 21:55:41 mysinger Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2002/05/03 16:11:54  hvogt
* Imported sources
*
*  
      SUBROUTINE GZREAD (HANDLE, CHBUF, NBBUF, ISTAT)
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

      CHARACTER(*)   CHBUF
      character*80 temp
      INTEGER         NBBUF, ISTAT
*
      CALL GZGETS (HANDLE,NBBUF,CHBUF,ISTAT)
      if (istat.gt.0) then
c        chbuf(istat:istat)=char(10)
c        chbuf(istat+1:istat+1)=char(0)
*        istat=istat-1
         temp = chbuf(1:istat-1)
         chbuf = temp
c //char(0)
c        #//char(10)//char(0)
c       chbuf = temp
c       write(*,'(A80)')chbuf
      endif

      END
*
