/*
*  
* $Id: ftn_gzio.c,v 1.1.1.1 2007-12-12 21:55:41 mysinger Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2002/05/03 16:11:54  hvogt
* Imported sources
*
*  
*/
/* ftn_gzio - usage of the zlib compression library by fortran calls
 * Copyright (C) 2001 Harald Vogt
 * derived from gzio
 * Copyright (C) 1995-1998 Jean-loup Gailly.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

extern void exit  OF((int));


/* ===========================================================================
 * fortran interface to open .gz files for read/write
 *
 *    SUBROTINE GZIOOP (FILDES, MODE, FILENAME, LGFNAME, ISTAT)
 *
 *       FILDES     address of file descriptor area
 *       MODE       string selecting IO mode (r, w)
 *       FILENAME   name of the file (including the null-terminator)
 *       LGFNAME    length of the filename string
 *      *ISTAT      status, =zero if success
 */

void gzioop_(gzFile *fildes,char *fmode,char *fname,int *lgfname,int *stat)
{
   char *pttext;
   gzFile file;

   *stat   = -1;

   if(*lgfname==0) {
        perror("GZIOOP: no filename");
        *stat = 107;
        return;
   }
   if(!(pttext=(char*)malloc(*lgfname+8))) {
         fprintf(stderr, "GZIOOP: no malloc possible\n");
         exit(1);
   }
   strncpy(pttext,fname,*lgfname);
   pttext[*lgfname]='\0';

   switch(*fmode) {
     case 'r':
           if (pttext[0] == '-') {
               file = gzdopen(fileno(stdin), "r");
               if (file == NULL) {
                   fprintf(stderr, "GZIOOP: stdin file open error\n");
                   *stat = 108;
                   free(pttext);
                   return;
               }
           } else {
               file = gzopen(pttext, "rb");
               if (file == NULL) {
                   fprintf(stderr, "GZIOOP: input file open error\n");
                   *stat = 108;
                   free(pttext);
                   return;
               }
           }
           break;
     case 'w':
           file = gzopen(pttext, "wb");
           if (file == NULL) {
               fprintf(stderr, "GZIOOP: output file open error\n");
               *stat = 108;
               free (pttext);
               return;
           }
           break;
     default:  fprintf(stderr,"GZIOOP: unknown i/o mode %s\n",fmode);
               *stat = 108;
               free (pttext);
               return;
   }
/* myprint_(file);  */
   *fildes = file;
   *stat   = 0;
   free (pttext);
   return;
}

/*
 *
 *    SUBROTINE GZCLOS (FILDES, MODE, ISTAT)
 *
 *       FILDES     address of file descriptor area
 *       MODE       string selecting IO mode (r, w)
 *      *ISTAT      status, =zero if success
 */

void gzclos_(gzFile *fildes,int *stat)
{
   gzFile file;

   *stat   = -1;
   file    = *fildes;

           gzclose(file);
           if (file == NULL) {
               fprintf(stderr, "GZCLOS: no file defined to close\n");
               exit(1);
           }
   *stat   = 0;
   return;
}

/*
 *
 *    SUBROTINE GZPUTS (FILDES, NBREC, MBUF, ISTAT)
 *
 *       FILDES     address of file descriptor area
 *       NBREC      record size, number of bytes to be written
 *       MBUF       vector to be written
 *      *ISTAT      status, =zero if success
 */

void gzputs_(gzFile *fildes,int *nbrec, char *mbuf, int *stat)
{
   int  nbdn, nbdo;
   gzFile file;

   *stat   = -1;
   file    = *fildes;

   nbdo   = *nbrec;

   if (file == NULL) {
       fprintf(stderr, "GZPUTS: no file defined to write\n");
       exit(1);
   }
   nbdn    = gzwrite(file,  mbuf, nbdo);
   *stat   = nbdn;
   return;
}

/*
 *
 *    SUBROTINE GZPUTB (FILDES, NBREC, MBUF, LEN, ISTAT)
 *
 *       FILDES     address of file descriptor area
 *       NBREC      record size, number of bytes to be written
 *       MBUF       vector to be written
 *      *ISTAT      status, =zero if success
 */

void gzputb_(gzFile *fildes, char *mbuf, int* length, int *stat)
{
   int  nbdn, nbdo;
   gzFile file;

   *stat   = -1;
   file    = *fildes;

   nbdo   = *length;

   if (file == NULL) {
       fprintf(stderr, "GZPUTS: no file defined to write\n");
       exit(1);
   }
   nbdn    = gzwrite(file,  mbuf, nbdo);
   *stat   = nbdn;
   return;
}

/*
 *
 *    SUBROTINE GZGETS (FILDES, NBREC, MBUF, ISTAT)
 *
 *       FILDES     address of file descriptor area
 *       NBREC      record size, number of bytes of array MBUF
 *       MBUF       array for storing unzipped info
 *      *ISTAT      status, =NBREC if ok., =zero if EOF encountered
 */

void gzgets_(gzFile *fildes,int *nbrec, char *mbuf, int *stat)
{
   unsigned len;
//   char *buf;
   gzFile file;

   *stat   = 0;
   file    = *fildes;

   len = *nbrec;
/*   buf = mbuf;
   if (file == NULL) {
       fprintf(stderr, "GZGETS: no file defined to write\n");
       exit(1);
   }
   if (buf == Z_NULL || len <= 0) *stat = 0;
   while (--len > 0 && gzread(file, buf, 1) == 1 && *buf++ != '\n') ;
   *buf = '\0';
   *stat   = buf - mbuf;

*/

   // mjk jji
   if (gzgets(file,mbuf,*nbrec) == Z_NULL) {
     *stat = -1;
   } else {
     *stat=strlen(mbuf);
   }
   // end

   return;
}

/*
 *
 *    SUBROTINE GZGETB (FILDES, NBREC, MBUF, ISTAT)
 *
 *       FILDES     address of file descriptor area
 *       NBREC      record size, number of bytes of array MBUF
 *       MBUF       array for storing unzipped info
 *      *ISTAT      status, =NBREC if ok., =zero if EOF encountered
 */

void gzgetb_(gzFile *fildes,int *nbrec, char *mbuf, int *stat)
{
   int  nbdn, nbdo;
   gzFile file;

   *stat   = -1;
   file    = *fildes;

   nbdo   = *nbrec;

   if (file == NULL) {
       fprintf(stderr, "GZGETB: Error - no input file defined\n");
       exit(1);
   }
   nbdn    = gzread(file,  mbuf, nbdo);
   *stat   = nbdn;
   return;
}

// Michael Mysinger 12/12/2007

void gzteller_(gzFile *file, int *loc)
{
  *loc = (int) gztell(*file);
  return;
}

void gzseeker_(gzFile *file, int *offset, int *istat)
{

  *istat = gzseek(*file, (z_off_t) *offset, SEEK_SET);
  return;
}

void gzflusher_(gzFile *file, int *istat)
{
  *istat = gzflush(*file, Z_SYNC_FLUSH);
  return;
}

// end
