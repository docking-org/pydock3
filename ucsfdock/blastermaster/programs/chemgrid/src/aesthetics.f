        SUBROUTINE aesthetics(output_unit,version,program_name)

C////////////////////////////////////////////////////////////////////////////
C PROGRAMMER:      Daniel A. Gschwend
C REV.DATE:      08-JUN-94
C PURPOSE:      Provide header information at run-time.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IMPLICIT NONE
      INCLUDE            'aesthetics.h'
      INTEGER            output_unit

C      Get the number of arguments on the command-line
      nargs=iargc()

C      Get today's date
        CALL DATE(the_date)

C      Get the current time
        CALL TIME(the_time)

C      Write the program header
      CALL get_hostname()
      WRITE(output_unit,901) program_name(1:INDEX(program_name,' ')-1),
     +                                          version
      WRITE(output_unit,902) the_date,the_time,nm(1:hostname_len)
      WRITE(output_unit,*)

0901      FORMAT( 'Running [35m', A, '[0m, DML v.', F4.2,'...' )
0902      FORMAT( 16X, A9, 2X, A8, 2X, '(', A, ')' )

      RETURN
      END




        SUBROUTINE get_hostname

C////////////////////////////////////////////////////////////////////////////
C PROGRAMMER:   David M. Lorber
C REV.DATE:     04-JUL-99
C PURPOSE:      Obtain name of linux machine at run-time.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IMPLICIT NONE
      INCLUDE            'aesthetics.h'
      integer hostnm

C      Get the name of machine being used
       hostname_len = hostnm(nm)
       hostname_len = findEOS(nm,40)
      RETURN
      END




      INTEGER      FUNCTION findEOS(strng,start)

C////////////////////////////////////////////////////////////////////////////
C PROGRAMMER:      Daniel A. Gschwend
C REV.DATE:      13-JUL-93
C HISTORY:      27-MAY-93 v1.00 
C PURPOSE:      Finds last non-blank character of a string which may 
C            contain spaces. Searches from character position 
C            "start" backwards.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IMPLICIT NONE

      CHARACTER      strng*132
      INTEGER            start

      findEOS=start
      DO WHILE ((findEOS.GT.0).AND.(strng(findEOS:findEOS).EQ.' '))
          findEOS=findEOS-1
      END DO

      RETURN
      END




      FUNCTION isUCLetter(letter)

C////////////////////////////////////////////////////////////////////////////
C PROGRAMMER:      Daniel A. Gschwend
C REV.DATE:      17-JUN-94
C HISTORY:      17-JUN-94 v1.00 
C PURPOSE:      Determines if a character is an uppercase letter.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IMPLICIT NONE
      INCLUDE            'aesthetics.h'

      CHARACTER*1      letter

      isUCLetter=(INDEX(alphabet_uc,letter).NE.0)

      RETURN
      END


      FUNCTION isLCLetter(letter)

C////////////////////////////////////////////////////////////////////////////
C PROGRAMMER:      Daniel A. Gschwend
C REV.DATE:      17-JUN-94
C HISTORY:      17-JUN-94 v1.00 
C PURPOSE:      Determines if a character is a lowercase letter.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IMPLICIT NONE
      INCLUDE            'aesthetics.h'

      CHARACTER*1      letter

      isLCLetter=(INDEX(alphabet_lc,letter).NE.0)

      RETURN
      END

      FUNCTION isDigit(letter)

C////////////////////////////////////////////////////////////////////////////
C PROGRAMMER:      Daniel A. Gschwend
C REV.DATE:      17-JUN-94
C HISTORY:      17-JUN-94 v1.00 
C PURPOSE:      Determines if a character is a digit.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IMPLICIT NONE
      INCLUDE            'aesthetics.h'

      CHARACTER*1      letter

      isDigit=(INDEX(digits,letter).NE.0)

      RETURN
      END


      FUNCTION isHydrogen(atom_name)

C////////////////////////////////////////////////////////////////////////////
C PROGRAMMER:      Daniel A. Gschwend
C REV.DATE:      17-JUN-94
C HISTORY:      17-JUN-94 v1.00 
C PURPOSE:      Determines if a 4 character atom name is a hydrogen.
C            Method:  an atom name containing a capital H NOT preceded
C            by any uppercase letters is considered to be a hydrogen.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IMPLICIT NONE
      INCLUDE            'aesthetics.h'
      INTEGER            i,j

      CHARACTER*4      atom_name

      j=INDEX(atom_name,'H')
      isHydrogen=(j.NE.0)
      IF (isHydrogen) THEN
          i=1
          DO WHILE ((isHydrogen).AND.(i.LT.j))
            isHydrogen=(.NOT.isUCLetter(atom_name(i:i)))
            i=i+1
          END DO
      END IF

      RETURN
      END



      SUBROUTINE getTime(time)

C////////////////////////////////////////////////////////////////////////////
C PROGRAMMER:      Daniel A. Gschwend
C REV.DATE:      17-JUN-94
C HISTORY:      17-JUN-94 v1.00 
C PURPOSE:      Used in computing elapsed CPU time.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      DIMENSION      tarray(2)
      REAL            time,etime,tarray

      time = etime(tarray)

      RETURN
      END



      CHARACTER*10      FUNCTION num2char(num)

C////////////////////////////////////////////////////////////////////////////
C PROGRAMMER:      Daniel A. Gschwend
C REV.DATE:      14-MAY-92
C HISTORY:      14-MAY-92 v1.00 
C PURPOSE:      Converts an integer less than 10^11 to a 10 character string
C            Useful for appending numeric suffixes to filenames, e.g.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IMPLICIT NONE

      INTEGER            length, num, c, dig,d,num2

      num2=num
      DO c=1,10
          d=10-c
          dig=num2/(10**d)
          num2char(c:c)=CHAR( 48+dig )
          num2=num2-dig*10**d
      END DO

      RETURN
      END
