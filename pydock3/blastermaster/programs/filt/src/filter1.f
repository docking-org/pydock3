*******************************************************************************
CCC This program filters the output of the docking programs, using atoms in
CCC the pdb file and atoms from the output (of the docking programs), 
CCC both of which are specified by the
CCC user, as necessary match criteria.   If the atoms specified are within
CCC a certain distance (also user specified) for a particular ligand
CCC orientation in the output file, that orientation is accepted as fulfilling
CCC the filter requirements.  The user is allowed to specify boolean
CCC operators for the pdb atoms (a maximum of two pdb atoms to be linked by 
CCC booleans) and for ligand atoms (he must specify a boolean between each
CCC ligand atom specified).
CCC April, 1987.  Brian Shoichet.
CCC This version adapted for VMS using filtered dUMP as the 'receptor'.
CCC May 14, 1987.
CCC Improved read statement in subroutine NEWFIL allows for end of input file
CCC escape with and 'end='statement.  Arises from too many ligands in command 
CCC file.  June 28, 1987.
CCC July, 2011. Michael Mysinger
CCC Added chain ID handling to GETFIL, INTERNAL, and NEWFILL
CCC This is enough for MakeDOCK to give the chain ID to DMS
*******************************************************************************
	PARAMETER (MAXATM=16,MAXLIG=40000,MAXLIN=300000) 
	CHARACTER*80 RSTRING, LSTRING, CSTRING, REMARK(MAXLIN)
	CHARACTER*80 ANAME, BNAME, CNAME
	CHARACTER*80 INTEXT, MSRUN, RESATM, NANU
	CHARACTER*3 RBOOL(MAXATM), LBOOL(MAXATM)
	CHARACTER*4 RATNUM(MAXATM),LATNUM(MAXATM),DUMMY(MAXLIN)
	CHARACTER*4 ATYPE(MAXLIN),RATYPE(MAXLIN),DUMR(MAXLIN)
	CHARACTER*4 ARES(MAXLIN),RARES(MAXLIN),SUMNAM(MAXLIN/6)
	CHARACTER*4 XATNUM(MAXATM)
	CHARACTER*1 CHAIN(MAXLIN), RCHAIN(MAXLIN),SUMCHAIN(MAXLIN/6)
	INTEGER NRAT(MAXATM),ATNUM(MAXLIN),IRNM(MAXLIN),SUMRES(MAXLIN/6)
	INTEGER FSET,RIRNM(MAXLIN),INCNUM, INUM(MAXLIN), FILNUM, CNTPLS
	REAL CUT(MAXATM),CA(MAXLIN,3),CB(MAXLIN,3)
	REAL RRADS(MAXLIN),LRADS(MAXLIN),SUMWRT(MAXLIG,2,MAXATM/2)
	REAL LOCCUP(MAXLIN),ROCCUP(MAXLIN)
	LOGICAL CUTOFF(MAXLIG,2,MAXATM/2), BCUTOFF(MAXLIG,MAXATM/2)
	LOGICAL FCUTOFF(MAXLIG),LOGALL(MAXATM),NANLOG

  	WRITE(*,5)'FILTER PROGRAM'
 	WRITE(*,10)'THIS PROGRAM FILTERS ATOMIC COORDINATE FILES'
 	WRITE(*,11)'BASED ON DISTANCE AND SEQUENCE CRITERIA'
5	FORMAT(T30,A,/)
10	FORMAT(T12,A)
11	FORMAT(T14,A,/)

CCC   This subroutine handels the user interface including parsing of input
CCC   and conversion to/from capital letters, etc.
	CALL USERIO(ANAME,NANU,INCNUM,MSRUN,RSTRING,MAXATM,RBOOL,NR,NRAT,
     &  LOGALL,RATNUM,NANLOG,BNAME,LSTRING,LBOOL,LATNUM,NL,INTEXT,
     &  CSTRING,CUT,CNAME,LIGBUM,XATNUM,RESATM)
CCC
	OPEN(13,file=ANAME,status='old')
	OPEN(14,file=BNAME,status='old')
	OPEN(15,file=CNAME,status='new')

CCC   Load receptor file (13) into memory.
	CALL GETFIL(DUMR,ATNUM,RATYPE,RARES,RCHAIN,RIRNM,CB,
     &           MAXATM,MAXLIN,IR,13,RRADS,ROCCUP,REMARK)
	WRITE(6,*)'IR IS',IR

CCC  Determine the size of the repeating units in the ligand file.
	CALL REPEAT(ICOUNT,14)
 	WRITE(6,*)'Number of atoms in ligand - ',ICOUNT

CCC  Load ligand file (14) into memory.
	CALL GETFIL(DUMMY,INUM,ATYPE,ARES,CHAIN,IRNM,CA,
     &           MAXATM,MAXLIN,IL,14,LRADS,LOCCUP,REMARK)
	WRITE(6,*)'IL IS',IL

CCC   This IF LOOP (below) is used to decide what sort of 
CCC   filtering operations to preform, based on the user's choices.  
CCC   The options as of this writing are:
CCC   1) Internal filter; 2) Name filter; 3) Normal receptor/ligand filter.
CCC   Within each of these there are further options.  BKS.  Feb 17, 1988.
	n=1
	m=1
	WRITE(6,*)'NR IS-', NR
CCC
	IF(INTEXT.EQ.'Y'.OR.INTEXT.EQ.'y')THEN
		WRITE(6,*)'INTERNAL BEGUN'
		JCUT=1

CCC 		This subroutine is called when one wishes to filter 
CCC 		a file against itself, as when one wants to know 
CCC 		what residues are within a certain distance of a protein 
CCC 		active site, for instance.  If one is doing a 'self-filter'  
CCC 		one does not compare distances to some external file.
		CALL INTERNAL(MAXATM,CUT,MAXLIG,LATNUM,M,IL,JCUT,N,IRNM,
     & 	CUTOFF,DUMMY,ATYPE,CA,CB,IR,MAXLIN,INUM,ARES,CHAIN,ICOUNT,JLIG,
     &	SUMRES,SUMNAM,SUMCHAIN)
CCC
		WRITE(6,*)'INTERNAL COMPLETED'
		NUMLIG=JLIG
		WRITE(6,*)NUMLIG
CCC
	ELSEIF(NANLOG)THEN
		WRITE(6,*)'NAMFIL BEGUN'

CCC 		This subroutine filters on an atom name basis, with the  
CCC 		expectation that the atom name defined by the user will be 
CCC 		repeated periodically, as when the receptor file is really 
CCC		another ligand file.   The distance calculations are done 
CCC		here, and B and FCUTOFFS are then called.
		CALL NAMFIL(NR,RATNUM,MAXATM,NUMLIG,CUT,RBOOL,MAXLIG,LATNUM,
     & 	CUTOFF,IRNM,DUMMY,ATYPE,CA,MAXLIN,INUM,ARES,ICOUNT,LOGALL,
     & 	LIGBUM,ATNUM,CB,INTEXT,SUMRES,RATYPE,IR,IL,FCUTOFF,NL,LBOOL,
     &    BCUTOFF,DUMR,RARES,RIRNM,INCNUM,SUMWRT,RRADS,LRADS,
     &    ROCCUP,LOCCUP,REMARK)
	ELSE
CCC
CCC This subroutine does all of the distance calculations, matching the
CCC appropriate atoms from the receptor and ligands, using their 
CCC corresponding cutoff distances as the filtering criteria.  Booleans 
CCC are not considered except to determine which receptor atoms are 
CCC appropriate (up to two) for each ligand atom.  The filtered sets of
CCC ligands are returned in the array CUTOFF.
		CALL CUTOFFS(NR,NRAT,MAXATM,NUMLIG,CUT,RBOOL,MAXLIG,LATNUM,
     & 	CUTOFF,IRNM,DUMMY,ATYPE,CA,MAXLIN,INUM,ARES,ICOUNT,
     & 	LOGALL,LIGBUM,ATNUM,CB,INTEXT,IL,RATNUM,RATYPE,
     & 	NL,RARES,RIRNM,INCNUM,SUMWRT,IR,REMARK,CNTPLS)
	ENDIF
	WRITE(6,*)'NUMLIG IS-',NUMLIG
	IF (NUMLIG.GT.MAXLIG)THEN
		WRITE(6,*)'TOO MANY LIGANDS, INCREASE PARAMETER MAXLIG'
		WRITE(6,*)'SORRY CHARLIE, PROGRAM BOMBS'
          STOP
	ENDIF
CCC
	IF(.NOT.NANLOG)THEN

CCC This subroutine takes the filtered sets in the array CUTOFF and performs
CCC the boolean operations specified in the receptor command line on them.
CCC Those which meet the boolean conditions are returned in BCUTOFF.
	   CALL BCUTOFFS(RBOOL,NUMLIG,CUTOFF,NL,MAXATM,MAXLIG,BCUTOFF)
CCC
CCC This subroutine takes the filtered sets which have passed the boolean
CCC operations performed on them in BCUTOFFS and subjects them to the
CCC boolean operation specified in the ligand command line.
CCC The ones which pass are put into FCUTOFF, which contains the final
CCC form of the filtered ligands.
CCC
	   CALL FCUTOFFS(NUMLIG,NL,LBOOL,BCUTOFF,MAXATM,MAXLIG,FCUTOFF)
CCC
CCC This subroutine writes those ligands specified by FCUTOFF into CNAME,
CCC a file which contains a list of all those ligands which fulfiled the
CCC various filter conditions.
	   CALL NEWFILL(NUMLIG,FCUTOFF,MAXLIG,BCUTOFF,IRNM,DUMMY,ATYPE, 
     &  CA,MAXLIN,INUM,ARES,CHAIN,ICOUNT,CNAME,INTEXT,SUMRES,SUMNAM,
     &  SUMCHAIN,IL,DUMR,RARES,RIRNM,ATNUM,RATYPE,CB,MSRUN,RRADS,LRADS,
     &  ROCCUP,LOCCUP,REMARK,CNTPLS)

	ENDIF

CCC If filtering by receptor atom names, that is against a 'receptor' which
CCC is a collection of orientations, then do not call F/BCUTOFFS or NEWFIL
CCC since they will have already been called by NAMFIL from inside of 
CCC CUTOFFS.
CCC
	WRITE(6,*)'filtering complete' 
CCC 	STOP
	END
***************************************************************************
*					    USERIO					        *
* 												  *
* This subroutine of the main FILTER1 is responsible for interfacing with *
* the user.  BKS, 03/88									  *
***************************************************************************	
	SUBROUTINE USERIO(ANAME,NANU,INCNUM,MSRUN,RSTRING,MAXATM,RBOOL,NR,
     &  NRAT,LOGALL,RATNUM,NANLOG,BNAME,LSTRING,LBOOL,LATNUM,NL,INTEXT,
     &  CSTRING,CUT,CNAME,LIGBUM,XATNUM,RESATM)
CCC
	CHARACTER*80 RSTRING, LSTRING,CSTRING
	CHARACTER*80 ANAME, BNAME, CNAME
	CHARACTER*80 INTEXT, MSRUN, RESATM, NANU
	CHARACTER*3 RBOOL(MAXATM), LBOOL(MAXATM)
	CHARACTER*4 RATNUM(MAXATM),LATNUM(MAXATM),XATNUM(MAXATM)
	INTEGER NRAT(MAXATM), INCNUM, MAXATM, NL, LIGBUM, NR
	REAL CUT(MAXATM)
	LOGICAL LOGALL, NANLOG

*	WRITE(*,1)'What sort of files are you filtering?'
*	WRITE(*,3)'1. RECEPTOR vs. LIGAND'
*	WRITE(*,3)'2. LIGAND vs. LIGAND'
*	WRITE(*,3)'3. INTERNAL (FILTERING A FILE AGAINST ITSELF)'
*	WRITE(*,3)'4. EXIT'
*	READ(5,19)NANU

*	WRITE(*,1)'Filtering by atom or by residue?'
*	WRITE(*,3)'1. ATOM'
*	WRITE(*,3)'2. RESIDUE'
*	WRITE(*,3)'3. EXIT'

*	WRITE(*,1)'Filtering by specific atoms/res, or by atom/res types?'
* 	WRITE(*,3)'1. ATOM/RES NAMES (eg; CA, LYS)'
*	WRITE(*,3)'2. ATOM/RES NUMBERS (eg; 3415, 317)'
*	WRITE(*,3)'3. EXIT'

1	FORMAT(/,A)
*3 	FORMAT(/,T16,A)
	WRITE(*,1)'name of file to filter by -'
	READ(5,10)ANAME
10	FORMAT(A80)
	WRITE(6,*)ANAME

CCC  NANU allows for filtering against a receptor file which is likely
CCC  to have many repeating atom types that one wishes to filter by, such
CCC  as when the 'receptor' is really just another set of ligands.
	WRITE(6,*)'filter by atom names or numbers in receptor(na/nu)?'
	READ(5,19)NANU

CCC	Converts to upper case.
 	CALL TOUPPER(NANU)

CCC  INCNUM allows the user to only filter for residues/orientations which
CCC  differ in their numbers (ie; residue numbers) by INCNUM.  This is
CCC  useful when filtering a file against itself and only wanting to look at
CCC  i to i+INCNUM relationships.
	WRITE(6,*)'residue/ligand orientation number increment filter?'
	READ(5,21)INCNUM
	WRITE(6,21)INCNUM

CCC  MSRUN is used in NEWFIL to determine what format the output should
CCC  be in.  If MSRUN .EQ. 'y', then right in MS -i format.
	WRITE(6,*)'Filter for MS run?(Y/N)'
	READ(5,25)MSRUN
 	CALL TOUPPER(MSRUN)
	WRITE(6,*)MSRUN

	WRITE(6,*)'Either names or numbers of the atoms to filter by, as '
	WRITE(6,*)'they appear in the pdb, with boolean and/or operators '
	WRITE(6,*)'if such apply between them.  Maximum of one boolean '
	WRITE(6,*)'per two filter atoms.'	
	READ(5,20)RSTRING
 	CALL TOUPPER(RSTRING)

19	FORMAT(A2)
20	FORMAT(A80)
21	FORMAT(I9)

CCC  This subroutine parses the command line into filter atoms and boolean
CCC  operators.
	CALL RPARSE(RSTRING,MAXATM,RBOOL,NR,NRAT,LOGALL,RATNUM,NANU,NANLOG)

	WRITE(6,*)'name of ligand file-'
	READ(5,10)BNAME
	WRITE(6,*)BNAME
CCC
	WRITE(6,*)'names of ligand atoms to match against receptor '
	WRITE(6,*)'filters.  Include boolean and/ors between EACH '
	WRITE(6,*)'atom specified.'
	READ(5,20)LSTRING
 	CALL TOUPPER(LSTRING)

CCC  Parses command line for ligand atoms to filter by
	CALL LPARSE(LSTRING,MAXATM,LBOOL,LATNUM,NL)

	WRITE(6,*)'internal filtering?'
	READ(5,25)INTEXT,RESATM
 	CALL TOUPPER(INTEXT)
 	CALL TOUPPER(RESATM)
	WRITE(6,25)INTEXT,RESATM
25	FORMAT(A1,1X,A1)


	WRITE(6,*)'distance cutoffs for each of the receptor atoms 
     & specified-'
	READ(5,35)CSTRING

CCC  This subroutine parses the distance input string.
	CALL CPARSE(CSTRING,MAXATM,CUT,XATNUM)

*	WRITE(6,35)(CUT(i),i=1,NR)
35	FORMAT(A80)

CCC  The filtering distance are squared to make the distance calculations 
CCC  faster.  
	DO 40 i=1,NR
		CUT(i)=CUT(i)*CUT(i)
		WRITE(6,*)CUT(i)
40	CONTINUE

	WRITE(6,*)'output filename-' 
	READ(5,10)CNAME
	WRITE(6,*)CNAME

	WRITE(6,*)'number of allowed ligand contacts in an all-filter?'
	READ(5,50)LIGBUM
50	FORMAT(I2)
	WRITE(6,*)'LIGBUM IS-',LIGBUM
CCC
	RETURN
	END
************************************************************************
CCC This subroutine converts from lower to upper case all input.
CCC BKS, 1985.
************************************************************************
	SUBROUTINE TOUPPER(CLINE)

	CHARACTER*80 CLINE
	INTEGER      L

CCC 	LEN is an internal function.  Returns the length of a character 
CCC	string.
	L = LEN(CLINE)

CCC	Use ASCII alph-numeric conversion code (ICHAR) to turn 
CCC  into uppercase.
	DO 50 I=1,L
		IF(ICHAR(CLINE(I:I)).GE.97.AND.ICHAR(CLINE(I:I)).LE.122)THEN
			N=ICHAR(CLINE(I:I))-32
			CLINE(I:I)=CHAR(N)	
		ENDIF
50	CONTINUE

	RETURN
	END
*************************************************************************
CCC This subroutine of the program filter1.f parses the command line and
CCC identifies which words in it refer to atom numbers and which refer to
CCC boolean operators.  It also counts how many atoms to filter by.
CCC April, 1987.  Brian Shoichet.
************************************************************************
	SUBROUTINE RPARSE(RSTRING,MAXATM,RBOOL,NR,NRAT,LOGALL,RATNUM,NANU,
     &  NANLOG)


	CHARACTER*80 NANU
	CHARACTER*4 RATNUM(MAXATM)
	CHARACTER*3 RBOOL(MAXATM)
	CHARACTER*80 RSTRING, TEMP  
	INTEGER NRAT(MAXATM)
	LOGICAL LOGALL(MAXATM), NANLOG
 
 
	IBOOL=0
	IF(RSTRING.EQ.' ')GOTO 200

	DO 100 i=1,MAXATM

CCC 		Removes leading blanks from command line
		CALL BLANKS(RSTRING)

CCC		Finds the first word in the command line.
		CALL NBLANK(RSTRING,RATNUM,MAXATM,L,I)

CCC		If filtering a ligand file against a ligand file.
		IF(NANU.EQ.'NA')THEN
			NANLOG=.TRUE.
			IF(RATNUM(I)(1:L-1).EQ.'ALL')THEN
*				WRITE(6,*)'ITS TRUE FOR I=',I
				LOGALL(I)=.TRUE.
			ELSE
				LOGALL(I)=.FALSE.
			ENDIF
		ELSE
			NANLOG=.FALSE.
			IF(RATNUM(I)(1:L-1).EQ.'ALL')THEN
*				WRITE(6,*)'ITS TRUE FOR I=',I
				LOGALL(I)=.TRUE.
			ELSE
				LOGALL(I)=.FALSE.
				READ(RATNUM(I),*)NRAT(i)
  				WRITE(6,*)'RECEPTOR ATOM',I,'IS',RATNUM(i)
			ENDIF
		ENDIF
CCC
CCC If RATNUM(I) = all, then set up the logical LOGALL which in CUTOFFS
CCC determines whether to do an 'all-filter'.  NANU and by extension 
CCC NANLOG are variables, set by the user, which determine whether or not
CCC the 'receptor' file by atom name or number.  One would search by number
CCC if the 'receptor' is just that, a true macromolecule, and by atom name
CCC if it was a colection of small molecule orientations with the same
CCC chemical formula.
CCC
		TEMP=RSTRING(l+1:)
		RSTRING=TEMP 
*		WRITE(6,*)RSTRING

CCC       Sets last BOOL to 'nul'
		IF(RSTRING.EQ.' ')GOTO 190

		CALL BLANKS(RSTRING)

CCC       Parses for boolean command.
		CALL BOOL(RSTRING,RBOOL,MAXATM, I)

100	CONTINUE
CCC
190	RBOOL(i)='NUL'
*	WRITE(6,*)RBOOL(i)
	NR=i
CCC
CCC
	RETURN
200	WRITE(6,*)'No commands: nothing to filter by - program bombs'
	STOP 	
	END
****************************************************************************
CCC This subroutine of the program filter1.f parses the ligand command
CCC line.
CCC April, 1987.  Brian Shoichet.
****************************************************************************
	SUBROUTINE LPARSE(LSTRING,MAXATM,LBOOL,LATNUM,NL)
CCC
CCC
	CHARACTER*4 LATNUM(MAXATM)
	CHARACTER*3 LBOOL(MAXATM)
	CHARACTER*80 LSTRING, TEMP  
 

CCC  Prevents infinite loops.
	IF(LSTRING.EQ.' ')STOP

	DO 100 i=1,MAXATM

CCC 	     Prevents infinite loops.
	     IF(LSTRING.EQ.' ')GOTO 190

CCC 		Removes leading blanks from command line
		CALL BLANKS(LSTRING)

CCC		Finds the first word in the command line.
		CALL NBLANK(LSTRING,LATNUM,MAXATM,L,I)

		WRITE(6,*)'LIGAND ATOM',I,'IS',LATNUM(I)
		TEMP=LSTRING(l+1:)
		LSTRING=TEMP 

CCC       Sets last BOOL to 'nul'
		IF(LSTRING.EQ.' ')GOTO 190

		CALL BLANKS(LSTRING)

CCC		Parses for necessary boolean command.
		CALL BOOL(LSTRING,LBOOL,MAXATM, I)

		IF(LBOOL(I).EQ.'NUL')THEN
			WRITE(6,*)'Must have booleans on ligand command line.'
			STOP
CCC		     INSERT CODE TO LOOP BACK AND CALL FOR COMMAND LINE AGAIN
CCC 			OR EXIT
		ENDIF

100	CONTINUE
CCC
190	LBOOL(i)='NUL'
	NL=i
*	WRITE(6,*)'NL IS',NL
CCC
CCC
	RETURN
200	WRITE(6,*)'No commands: nothing to filter by - program bombs'
	STOP
	END
**********************************************************************
CCC This subroutine of the main FILTER1 is responsible for parsing
CCC the cutoff distances, stored in the array CUT, so as not to force
CCC the user to worry about how he spaces cut distances.  BKS, 07/87.
********************************************************************** 
	SUBROUTINE CPARSE(CSTRING,MAXATM,CUT,XATNUM)
CCC
CCC
	CHARACTER*4  XATNUM(MAXATM)
	CHARACTER*80 CSTRING, TEMP  
	INTEGER 	   I, MAXATM
	REAL 	   CUT(MAXATM)
 
CCC  Prevents infinite loops.
	IF(CSTRING.EQ.' ')GOTO 200

	DO 100 I=1,MAXATM

CCC 		Removes leading blanks from command line
		CALL BLANKS(CSTRING)

CCC		Finds the first word in the command line.
		CALL NBLANK(CSTRING,XATNUM,MAXATM,L,I)

*         WRITE(6,*)'XATNUM IS -',XATNUM(I),I
		READ(XATNUM(I),*)CUT(i)
*		WRITE(6,*)CUT(i)

		TEMP=CSTRING(L+1:)
		CSTRING=TEMP 
*		WRITE(6,*)CSTRING
		IF(CSTRING.EQ.' ')GOTO 190

100	CONTINUE

190	RETURN
200	WRITE(6,*)'No distance to filter with - program bombs'
	STOP
	END
*********************************************************************
*						BLANKS						   *
* 													   *
* This subroutine removes leading blanks from a command line.       *
*********************************************************************
	SUBROUTINE BLANKS(STRING)

	CHARACTER*80 STRING,TEMP

20	IF(STRING(1:1).EQ.' ')THEN
		TEMP=STRING(2:)
		STRING=TEMP
		GOTO 20		
	ENDIF
* 	WRITE(6,*)'STRING IS -',STRING
CCC
	RETURN
	END
********************************************************************
*						NBLANK						  *
*													  *
* This subroutine finds the next word in the command line, where   *
* word means characters delimited by blanks.					  *
********************************************************************
        SUBROUTINE NBLANK(STRING,AATNUM,MAXATM,L,I)

        CHARACTER*80 STRING
        CHARACTER*4  AATNUM(MAXATM)
        INTEGER       MAXATM,L,I

        l=INDEX(STRING,' ')
        AATNUM(I)=STRING(1:l-1)
*       WRITE(6,*)'AATNUM IS-',AATNUM(I)(1:L),I

        RETURN
        END
**************************************************************************
*						BOOL								   *
*														   *
* This subroutine parses for boolean operators on the command line.      *
**************************************************************************
	SUBROUTINE BOOL(STRING,ABOOL,MAXATM,I)

	CHARACTER*80 STRING, TEMP
	CHARACTER*3  ABOOL(MAXATM)
	INTEGER      m, MAXATM, I

	m=INDEX(STRING,' ')
	ABOOL(i)=STRING(1:m-1)

CCC  The next word in the command line is checked for boolean character.  
CCC  If it is not an 'and' an 'or'  or a 'not' it is set to 'nul'.
	IF(ABOOL(i).NE.'AND'.AND.ABOOL(i).NE.'OR '.AND.ABOOL(I)
     & 					.NE.'NOT')THEN
		ABOOL(i)='NUL'
	ELSE
		TEMP=STRING(m+1:)
		STRING=TEMP
	ENDIF
 	WRITE(6,*)'BOOLEAN',I,' IS',ABOOL(i)

	RETURN
	END
**************************************************************************
CCC This subroutine loads the receptor and ligand files into memory.
CCC BKS, Feb 17, 1988.
**************************************************************************
	SUBROUTINE GETFIL(DUM,SATNUM,SATYPE,SARES,CHAIN,SIRNM,SB,
     &           MAXATM,MAXLIN,IS,FILNUM,RADS,OCCUP,REMARK)

	CHARACTER*80 STRING, REMARK(MAXLIN)
	CHARACTER*4 DUM(MAXLIN)
	CHARACTER*4 SARES(MAXLIN),SATYPE(MAXLIN)
	CHARACTER*1 CHAIN(MAXLIN)
	INTEGER SATNUM(MAXLIN),IS,SIRNM(MAXLIN),FILNUM
	REAL SB(MAXLIN,3),RADS(MAXLIN), OCCUP(*)

20	FORMAT(A4,2X,I5,2X,A4,A3,I6,4X,3F8.3,2F6.2)
22	FORMAT(A4,2X,I5,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3,2F6.2)

CCC  The file is read in and stored in various separate arrays.  It is 
CCC  not stored as a simple (A80) type array to save both memory and disk 
CCC  space (when it comes to writing it).
	DO 35 IS=1,MAXLIN
          READ(FILNUM,'(A80)',END=39)STRING
          IF(STRING(1:3).EQ.'ATO'.OR.STRING(1:3).EQ.'HET'.OR.
     &       STRING(1:3).EQ.'TER')THEN
              REMARK(IS)='X'
              READ(STRING,22)DUM(IS),SATNUM(IS),SATYPE(IS),SARES(IS),
     &          CHAIN(IS),SIRNM(IS),(SB(IS,J),J=1,3),OCCUP(IS),RADS(IS)
          ELSE
               IF(FILNUM.EQ.14)REMARK(IS)=STRING
          ENDIF
35	CONTINUE
39	CONTINUE
	
	IF(IS .GE. MAXLIN)THEN
		WRITE(6,*)'File number',FILNUM,'bigger than arrays, increase'
		WRITE(6,*)'parameter MAXLIN at the beginning of the code.'
		STOP
	ENDIF

	IS=IS-1

	RETURN
	END
**************************************************************************
*						REPEAT					 *
*												 *
* This subroutine counts the size of the repeating units in what will    *
* generally be the ligand file.  Repeating units = molecule size in      *
* lines of the input file.     BKS, 03/88.					 *
************************************************************************** 
	SUBROUTINE REPEAT(ACOUNT,FILNUM)

	CHARACTER*3 DUM1,DUM2
	INTEGER FILNUM,ACOUNT

	ACOUNT=0
	REWIND(FILNUM)
40	READ(FILNUM,45,END=100)DUM1
45	FORMAT (A3)
	IF(DUM1.EQ.'TER')THEN
50		READ(FILNUM,45,END=100)DUM2 
*		IF(DUM2.NE.'TER'.AND.DUM2.NE.'REM')THEN
		IF(DUM2.NE.'TER')THEN
			ACOUNT=ACOUNT+1
			GOTO 50
*         ELSEIF(DUM2.EQ.'REM')THEN
*			GOTO 50
		ENDIF
	ELSE
		GOTO 40
	ENDIF
CCC   The purpose of the above is to count how many atom lines there are
CCC   in the pdb file between TER lines.  This information to be used
CCC   in the write loop to the output files as well as in the various 
CCC   CUTOFFS setting subroutines. 

100	CONTINUE
*	WRITE(6,*)'FILNUM IS ',FILNUM
	REWIND(FILNUM)

	RETURN
	END
*************************************************************************
CCC This subroutine of the program filter1 does the initial distance
CCC calculations for all of the filtering atoms in the receptor
CCC command line.  Ligands orientations which pass the specified 
CCC cutoffs (stored in CUT) are stored in the array CUTOFF as logical
CCC .TRUE., .FALSE. if they fail.  
CCC April, 1987.  Brian Shoichet.
************************************************************************
	SUBROUTINE CUTOFFS(NR,NRAT,MAXATM,NUMLIG,CUT,RBOOL,MAXLIG,LATNUM,
     & CUTOFF,IRNM,DUMMY,ATYPE,CA,MAXLIN,INUM,ARES,ICOUNT,LOGALL,
     & LIGBUM,ATNUM,CB,INTEXT,IL,RATNUM,RATYPE,NL,RARES,RIRNM,INCNUM,
     & SUMWRT,IR,REMARK,CNTPLS)

	CHARACTER*3 RBOOL(MAXATM)
	CHARACTER*4 LATNUM(MAXATM),DUMMY(MAXLIN),ATYPE(MAXLIN)
	CHARACTER*4 ARES(MAXLIN),RATYPE(MAXLIN),RARES(MAXLIN)
	CHARACTER*4 RATNUM(MAXATM)
	CHARACTER*80 REMARK(MAXLIN)
	INTEGER NL,ICOUNT,JCUT,j,INCNUM, INUM(MAXLIN),CNTPLS,JUMP
	INTEGER ATNUM(MAXLIN),IR,NRAT(MAXATM),IRNM(MAXLIN),RIRNM(MAXLIN)
	REAL CUT(MAXATM), CA(MAXLIN,3), CB(MAXLIN,3)
	REAL SUMWRT(MAXLIG,2,MAXATM/2)
	LOGICAL CUTOFF(MAXLIG,2,MAXATM/2), LOGALL(MAXATM)

	CNTPLS=1
	m=1
	n=1

	DO 200 j=1,NR
		IF(LOGALL(J))THEN
			JCUT=J

CCC 			This subroutine filters the ligand file 
CCC 			against all atoms in the receptor file.  It is useful 
CCC			for 'bad contact' type calculations.
			CALL ALLFILT(MAXATM,CUT,MAXLIG,LATNUM,M,IL,JCUT,N,
     & 		LIGBUM,CUTOFF,DUMMY,ATYPE,CA,CB,ATNUM,IR,MAXLIN,INUM,
     & 		ARES,ICOUNT,INCNUM,IRNM,RIRNM,NUMLIG,REMARK,CNTPLS)

		ELSE
			NUMLIG=0

CCC            The following two loops set up for use of the 
CCC            ICOUNT iterator in the k loop.  BKS 08/16/87.
			DO 60 LR=1,IR
*				WRITE(6,400)NRAT(J),ATNUM(LR)
				IF(NRAT(J).EQ.ATNUM(LR))GOTO 64
60			CONTINUE

64			CONTINUE
               K=0
               CNTPLS=1
   			DO 65 L=1,ICOUNT
63				K=K+1
                    IF(DUMMY(K).EQ.'TER ')GOTO 63
                    IF(REMARK(K)(1:4).EQ.'REMA')THEN
                         CNTPLS=CNTPLS+1
                         GOTO 63
                    ENDIF
				IF(ATYPE(K).EQ.LATNUM(M))THEN
					GOTO 66
				ENDIF
65			CONTINUE

66 			KLIG=0
               L=K
               JUMP=ICOUNT+CNTPLS

CCC            Here is where ICOUNT comes in handy, allowing the 
CCC            computer to know where to look for the next 
CCC            corrientation's o-ordinates for the atom type of interest.
			DO 100 k=L,IL,JUMP

				KLIG=KLIG+1

CCC 				The distance calculation is performed here.
				SUM = DISTANCE(CA,CB,K,LR,MAXLIN)

*                   WRITE(6,*)ATYPE(K),SUM,KLIG
CCC 				CUTOFF is given logical values for each 
CCC 				ligand orientation.
				CALL SETCUT(SUM,J,CUT,MAXATM,MAXLIG,KLIG,
     &			   			N,M,SUMWRT,CUTOFF)

100			CONTINUE
			NUMLIG=KLIG
		ENDIF

CCC 		Fixes the values of indexes N and M based on whether or
CCC 		not RBOOL(J) is a true boolean or has been set to nul.
		CALL POINTR(N,M,RBOOL,J,MAXATM)

200	CONTINUE
400	FORMAT(A4,2X,A4)
CCC
	RETURN
	END
***************************************************************************
CCC
***************************************************************************
	SUBROUTINE NAMFIL(NR,RATNUM,MAXATM,NUMLIG,CUT,RBOOL,MAXLIG,LATNUM,
     & CUTOFF,IRNM,DUMMY,ATYPE,CA,MAXLIN,INUM,ARES,ICOUNT,LOGALL,
     & LIGBUM,ATNUM,CB,INTEXT,SUMRES,RATYPE,IR,IL,FCUTOFF,NL,LBOOL,
     & BCUTOFF,DUMR,RARES,RIRNM,INCNUM,SUMWRT,RRADS,LRADS,
     & ROCCUP,LOCCUP,REMARK)

	CHARACTER*80 INTEXT, REMARK(MAXLIN)
	CHARACTER*4 DUMR(MAXLIN),RARES(MAXLIN)
	CHARACTER*4 LATNUM(MAXATM),DUMMY(MAXLIN),ATYPE(MAXLIN)
	CHARACTER*4 ARES(MAXLIN),RATNUM(MAXATM),RATYPE(MAXLIN)
	CHARACTER*3 DUM2,RBOOL(MAXATM),LBOOL(MAXATM),DUM1
	INTEGER ATNUM(MAXLIN),IR,RCOUNT,IRNM(MAXLIN),SUMRES(MAXLIN/6) 
	INTEGER RECCNT,IL,NR, M, N, RATRCN, NL,RIRNM(MAXLIN),INCNUM
	INTEGER INUM(MAXLIN), CNTPLS
	REAL    CUT(MAXATM), CA(MAXLIN,3), CB(MAXLIN,3)
	REAL    SUMWRT(MAXLIG,2,MAXATM/2),RRADS(MAXLIN),LRADS(MAXLIN)
	REAL    ROCCUP(*),LOCCUP(*)
	LOGICAL CUTOFF(MAXLIG,2,MAXATM/2),LOGALL(MAXATM),FCUTOFF(MAXLIG) 
	LOGICAL BCUTOFF(MAXLIG,MAXATM/2)

	KLIG=0
	WRITE(6,*)'INCNUM IS-', INCNUM
  	RCOUNT=0
	REWIND(13)

CCC	Determine the size of the repeating units in the 'receptor' file.
	CALL REPEAT(RCOUNT,13)
 	WRITE(6,*)'Number of atoms in receptor - ',RCOUNT

	DO 200 RECCNT=0,IR-1,RCOUNT+1
CCC 		The booolean pointers m and n must be reset 
CCC 		for every receptor orientation since each is treated 
CCC 		has an independent molecule in space having no knowledge 
CCC 		of the rest of the 'receptor' file.
		M=1
		N=1
		DO 150 J=1,NR
			IF(LOGALL(J))THEN
				JCUT=J

CCC 				This subroutine filters the ligand file 
CCC 				against all atoms in the receptor file.  It is 
CCC 				useful for 'bad contact' type calculations.
				CALL ALLFILT(MAXATM,CUT,MAXLIG,LATNUM,M,IL,JCUT,N,
     & 			LIGBUM,CUTOFF,DUMMY,ATYPE,CA,CB,ATNUM,IR,
     & 			MAXLIN,INUM,ARES,ICOUNT,INCNUM,IRNM,NIRNM,NUMLIG
     &              REMARK,CNTPLS)

 			ELSE

CCC 				The following two loops fix the pointers for 
CCC 				the user defined filtering atoms for the receptor 
CCC 				and the ligand, respectively.  BKS 08/16/87.
				DO 60 LR=1,RCOUNT
					IF(RATNUM(J).EQ.RATYPE(LR))GOTO 64
60				CONTINUE

64 				CONTINUE
				RATRCN=LR+RECCNT
				DO 65 L=1,ICOUNT
					IF(ATYPE(L).EQ.LATNUM(M))GOTO 66
65				CONTINUE 

66	 			CONTINUE

	               KLIG=0
				DO 100 k=L,IL,ICOUNT+1
				   	KLIG=KLIG+1

CCC 					The INCNUM IF loop checks for residues or
CCC 					orientations  which are too close by sequence.
				   	IF(RIRNM(RATRCN) .LT. (IRNM(K)-INCNUM))THEN 

CCC 						The distance calculation is performed here.
						SUM = DISTANCE(CA,CB,K,RATRCN,MAXLIN)
*                             WRITE(6,*)'SUM IS-',SUM,RATRCN,K

CCC 						CUTOFF is given logical values for each 
CCC 						ligand orientation.
						CALL SETCUT(SUM,J,CUT,MAXATM,MAXLIG,KLIG,
     &			   			N,M,SUMWRT,CUTOFF)

				   	ELSE
						CUTOFF(KLIG,n,m)=.FALSE.
				   	ENDIF
100				CONTINUE
				NUMLIG=KLIG
			ENDIF
*              WRITE(6,*)'numlig is-',NUMLIG

CCC 			Fixes the values of indexes N and M based on whether or
CCC 			not RBOOL(J) is a true boolean or has been set to nul.
			CALL POINTR(N,M,RBOOL,J,MAXATM)

150		CONTINUE
**		WRITE(6,*)'n and m are -',n,m
*		WRITE(6,*)'RECCNT IS-',RECCNT

		CALL BCUTOFFS(RBOOL,NUMLIG,CUTOFF,NL,MAXATM,MAXLIG,BCUTOFF)

		CALL FCUTOFFS(NUMLIG,NL,LBOOL,BCUTOFF,MAXATM,MAXLIG,FCUTOFF)

		CALL NANFILL(NUMLIG,FCUTOFF,MAXLIG,BCUTOFF,IRNM,DUMMY,ATYPE, 
     &  	CA,MAXLIN,INUM,ARES,ICOUNT,CNAME,INTEXT,SUMRES,IL,DUMR,RARES,
     &  	RIRNM,ATNUM,RATYPE,CB,RECCNT,RCOUNT,SUMWRT,MAXATM,LATNUM,
     &    RRADS,LRADS,ROCCUP,LOCCUP)

200	CONTINUE


	RETURN
	END
***********************************************************************
* 						SETCUT				    *
*											    *
* This subroutine of the main FILTER1 decides whether an atom or group*
* has met the filtering criteria and CUTOFF is set to .TRUE. (if it   *
* has) or .FALSE. (if it hasn't).                                     *
***********************************************************************
	SUBROUTINE SETCUT(SUM,J,CUT,MAXATM,MAXLIG,KLIG,N,M,SUMWRT,CUTOFF)

	INTEGER J, KLIG, N, M, MAXATM, MAXLIG
	REAL SUM, CUT(MAXATM), SUMWRT(MAXLIG,2,MAXATM/2)
	LOGICAL CUTOFF(MAXLIG,2,MAXATM)

CCC 	CUTOFF is given logical values for each ligand orientation.
	IF(SUM.LT.CUT(j))THEN
		CUTOFF(KLIG,n,m)=.TRUE.
*		WRITE(6,*)KLIG,n,m,SUM,CUTOFF(KLIG,n,m)

CCC 		SUMWRT is used in NANFIL to print the distances out w/the 
CCC 		receptor/ligand info.  BKS 02/88. 
		SUMWRT(KLIG,N,M)=SUM
	ELSE
		CUTOFF(KLIG,n,m)=.FALSE.
*		WRITE(6,*)KLIG,n,m,SUM,CUTOFF(KLIG,n,m)
	ENDIF

*	WRITE(6,*)CUTOFF(KLIG,n,m)

	RETURN
	END
*********************************************************************
*						POINTR						   *
*													   *
* This subroutine of the main FILTER1 decides on what the value of  *
* the indexes N and M to the array CUTOFF are based on whether or   *
* not a boolean was called for between two filtering elements of    *
* the receptor command line.  BKS 03/88                             *
*********************************************************************
	SUBROUTINE POINTR(N,M,RBOOL,J,MAXATM)

	CHARACTER*3 RBOOL(MAXATM)
	INTEGER 	  N, M, J, MAXATM

	N=N+1
	IF(RBOOL(j).EQ.'NUL')THEN
		M=M+1
		N=1
	ENDIF
*	WRITE(6,*)RBOOL(j),n,m

	RETURN
	END
*************************************************************************
CCC
*************************************************************************
	SUBROUTINE ALLFILT(MAXATM,CUT,MAXLIG,LATNUM,M,IL,JCUT,N,LIGBUM,
     & CUTOFF,DUMMY,ATYPE,CA,CB,ATNUM,IR,MAXLIN,INUM,ARES,ICOUNT,
     & INCNUM,IRNM,RIRNM,NUMLIG,REMARK,CNTPLS)


	CHARACTER*80 REMARK(MAXLIN)
	CHARACTER*4  LATNUM(MAXATM),DUM,DUMMY(MAXLIN),ATYPE(MAXLIN)
	CHARACTER*4  ARES(MAXLIN)
	CHARACTER*3  DUM2
	INTEGER      ATNUM(MAXLIN),INCNUM,IRNM(MAXLIN),CNTPLS
	INTEGER      RIRNM(MAXLIN),INUM(MAXLIN),JCUT,JUMP
	REAL         CUT(MAXATM), CA(MAXLIN,3), CB(MAXLIN,3)
	LOGICAL      CUTOFF(MAXLIG,2,MAXATM/2)

	IF(LATNUM(M).EQ.'ALL ')THEN
		JLIG=0
		DO 300 J=1,IL,ICOUNT+1
			JLIG=JLIG+1
	          NUMLIG=JLIG
			NBUMP=0
			DO 200 I=1,IR 
*			  IF(IRNM(J) .LT. (RIRNM(I)-INCNUM) .OR. IRNM(J) .GT.
*     &				(RIRNM(I)+INCNUM))THEN
				DO 100 K=J,J+ICOUNT-1

CCC					The distance calculation is performed here.
					SUM = DISTANCE(CA,CB,K,I,MAXLIN)

					IF(SUM.LT.CUT(JCUT))THEN
*                             WRITE(6,*)SUM,K,J,I,(CB(I,II),II=1,3)
						CUTOFF(JLIG,n,m)=.TRUE.
						NBUMP=NBUMP+1
						IF(NBUMP.GT.LIGBUM)THEN
*							WRITE(6,*)CUTOFF(JLIG,n,m), JLIG,n,m
							GOTO 300
						ENDIF

CCC In an 'all-fit' type loop, one only needs LIGNUM .TRUE. in the 100 loop
CCC to pass control back to the 300 loop, which looks at the next ligand
CCC orientation.  ie; one ligand-receptor 'contact' point is enough, one
CCC does not have to explore further in the loop.

					ELSE
						CUTOFF(JLIG,n,m)=.FALSE.
					ENDIF
100				CONTINUE
*			   ENDIF
200			CONTINUE
300		CONTINUE

CCC  The INCNUM loop filters for residues/orientations which are too close
CCC  in the sequence to be included.

	ELSE
		KLIG=0
          K=0
          CNTPLS=1
		DO 350 L=1,ICOUNT
340            K=K+1
               IF(DUMMY(K).EQ.'TER ')GOTO 340
               IF(REMARK(K)(1:4).EQ.'REMA')THEN
                    CNTPLS=CNTPLS+1
                    GOTO 340         
               ENDIF
			IF(ATYPE(K).EQ.LATNUM(M))GOTO 366
350		CONTINUE

366       CONTINUE

		L=K
          JUMP=ICOUNT+CNTPLS
    		DO 500 K=L,IL,JUMP
			KLIG=KLIG+1
	          NUMLIG=KLIG
			NBUMP=0
			DO 400 I=1,IR 

CCC				The distance calculation is performed here.
 				SUM = DISTANCE(CA,CB,K,I,MAXLIN)

*                   WRITE(6,*)SUM,KLIG,ATYPE(K),CA(K,1),CB(I,1)
				IF(SUM.LT.CUT(JCUT))THEN
                         NBUMP=NBUMP+1
					CUTOFF(KLIG,n,m)=.TRUE.
					IF(NBUMP.GT.LIGBUM)GOTO 500
				ELSE
					CUTOFF(KLIG,n,m)=.FALSE.
				ENDIF
400			CONTINUE
500		CONTINUE
	ENDIF
	RETURN
	END
*************************************************************************
*						DISTANCE						       *
*														  *
* The distance calculations between the filtering atom and the groups   *
* to be filtered file is performed here.   BKS 03/88.				  *
*************************************************************************
	FUNCTION DISTANCE(CA,CB,K,LR,MAXLIN)

	INTEGER ISUM, K, LR,MAXLIN
	REAL CA(MAXLIN,3), CB(MAXLIN,3)

	DISTANCE=0
	DO 75 ISUM=1,3
		DISTANCE=DISTANCE+(CA(K,ISUM)-CB(LR,ISUM))**2
75	CONTINUE

	END
*************************************************************************
CCC This subroutine of the program filter1.f is in charge of applying
CCC the boolean information in the ligand command line to the matches 
CCC stored in CUTOFFS.
CCC The list of ligands which survive the boolean operations are
CCC stored in BCUTOFF.
CCC April, 1987.  Brian Shoichet.
*************************************************************************
	SUBROUTINE BCUTOFFS(RBOOL,NUMLIG,CUTOFF,NL,MAXATM,MAXLIG,BCUTOFF)

	INTEGER NUMLIG,NL,MAXLIG,MAXATM
	CHARACTER*3 RBOOL(MAXATM)
	LOGICAL CUTOFF(MAXLIG,2,MAXATM/2), BCUTOFF(MAXLIG,MAXATM/2) 
 
 
CCC  Sets BCUTOFF to be TRUE - is useful for assessing NOT conditions
CCC  for LBOOL when one only has one condition (see FCUTOFF).
	DO 20 I=1,MAXLIG
		DO 10 J=1,MAXATM/2
               BCUTOFF(I,J)=.TRUE.
  10      CONTINUE
  20	CONTINUE

	j=1 
	DO 200 i=1,NL 
*		WRITE(6,*)RBOOL(j)

CCC 		If a boolean operation was called for.
		IF(RBOOL(j).NE.'NUL')THEN
			IF(RBOOL(J-1).EQ.'NOT')THEN
				RBOOL(J-1)='AND'
				DO 50 INOT=1,NUMLIG
					IF(CUTOFF(K,2,INOT))THEN
*						CUTOFF(K,2,INOT)=.TRUE.
 						CUTOFF(K,2,INOT)=.FALSE.
					ELSE
*						CUTOFF(K,2,INOT)=.FALSE.
 						CUTOFF(K,2,INOT)=.TRUE.
					ENDIF
50				CONTINUE
			ENDIF
			DO 100 k=1,NUMLIG

CCC 				Boolean AND performed.
				IF(RBOOL(j-1).EQ.'AND')THEN
					IF(CUTOFF(k,1,i).AND.CUTOFF(k,2,i))THEN
						BCUTOFF(k,i)=.TRUE.
 						WRITE(6,*)BCUTOFF(k,i),k,i
					ELSE
						BCUTOFF(k,i)=.FALSE.
					ENDIF
				ELSE
CCC 				Boolean OR.
					IF(CUTOFF(k,1,i).OR.CUTOFF(k,2,i))THEN
						BCUTOFF(k,i)=.TRUE.
*						WRITE(6,*)BCUTOFF(k,i),k
					ELSE
						BCUTOFF(k,i)=.FALSE.
					ENDIF
				ENDIF
100			CONTINUE
		ELSE
			IF(j.GT.1)THEN
*				WRITE(6,*)RBOOL(j-1)

CCC 				If this list has not already been checked above.
CCC                 ...1990 and i can't figure out what i was smoking...
*				IF(RBOOL(j-1).EQ.'NUL')THEN
					DO 150 k=1,NUMLIG
						BCUTOFF(k,i)=CUTOFF(k,1,i)
*						IF(BCUTOFF(k,i))THEN
**							WRITE(6,*)BCUTOFF(k,i),k
*						ENDIF
150					CONTINUE
*				ENDIF
			ELSE
				DO 175 k=1,NUMLIG
					BCUTOFF(k,i)=CUTOFF(k,1,i)
					IF(BCUTOFF(k,i))THEN
*						WRITE(6,*)BCUTOFF(k,i),k
					ENDIF
175				CONTINUE
			ENDIF
		ENDIF
		j=j+1
200	CONTINUE
	RETURN
	END
**************************************************************************
CCC This subroutine of the program of filter1.f applies the boolean 
CCC operators of the ligand command line to the logical array BCUTOFFS.
CCC Those ligands which satisfy the boolean operations are kept in
CCC the logical array FCUTOFFS.
CCC April, 1987.  Brian Shoichet.
***************************************************************************
	SUBROUTINE FCUTOFFS(NUMLIG,NL,LBOOL,BCUTOFF,MAXATM,MAXLIG,FCUTOFF)

	CHARACTER*3 LBOOL(MAXATM)
	INTEGER NUMLIG
	LOGICAL BCUTOFF(MAXLIG,MAXATM/2), FCUTOFF(MAXLIG),NOTLOG

	IF(NL.GT.1)THEN
	   DO 200 i=1,NL-1
		   IF(LBOOL(I).EQ.'NOT')THEN
*			   WRITE(6,*)'LBOOL(I) IS-',LBOOL(I),I
			   LBOOL(I)='AND'
			   NOTLOG=.TRUE.
			   DO 50 INOT=1,NUMLIG
				   IF(BCUTOFF(INOT,I+1))THEN
					   BCUTOFF(INOT,I+1)=.FALSE.
				   ELSE
					   BCUTOFF(INOT,I+1)=.TRUE.
				   ENDIF
50			   CONTINUE
		   ENDIF

CCC If lbool is not, then those orientations which were less than CUT
CCC should *not* be included, while those that were greater than CUT 
CCC *should* be included, hence the togling of the logicals.  NOTLOG
CCC records that a .not. togling has occured so that LBOOL(i) can be
CCC set back to 'not' from 'and' at the end of each 200 cycle.  
CCC This is only a meaningful change when the 'receptor' is really a
CCC multiple orientations file such that one is doing a ligandA X ligandB
CCC filter, where the size of A and B vary.

		   DO 100 j=1,NUMLIG
			   IF(LBOOL(i).EQ.'AND')THEN
				   IF(BCUTOFF(j,i).AND.BCUTOFF(j,i+1))THEN
					   FCUTOFF(j)=.TRUE.
**					   WRITE(6,*)FCUTOFF(j),j
				   ELSE
				   	   FCUTOFF(j)=.FALSE.
					   BCUTOFF(J,I+1)=.FALSE.
				   ENDIF
			   ELSE
				   IF(BCUTOFF(j,i).OR.BCUTOFF(j,i+1))THEN
					   FCUTOFF(j)=.TRUE.
*					   WRITE(6,*)FCUTOFF(j),j
				   ELSE
					   FCUTOFF(j)=.FALSE.
				   ENDIF
			   ENDIF
100		   CONTINUE
	        IF(NOTLOG)LBOOL(I)='NOT'
200	   CONTINUE
	ELSE
	   DO 300 K=1,NUMLIG
		IF(BCUTOFF(K,1))THEN
			FCUTOFF(K)=.TRUE.
*			WRITE(6,*)'FCUTOFF',FCUTOFF(K),K
		ELSE
			FCUTOFF(K)=.FALSE.
		ENDIF 
300	   CONTINUE
	ENDIF

	RETURN
	END
*************************************************************************
CCC This subroutine of the program filter1.f is responsible for writing
CCC those ligand orientations specified by FCUTOFFS from BNAME into CNAME.
CCC April, 1987.  Brian Shoichet.
*************************************************************************
	SUBROUTINE NEWFILL(NUMLIG,FCUTOFF,MAXLIG,BCUTOFF,IRNM,DUMMY,ATYPE,
     &  CA,MAXLIN,INUM,ARES,CHAIN,ICOUNT,CNAME,INTEXT,SUMRES,SUMNAM,
     &  SUMCHAIN,IL,DUMR,RARES,RIRNM,ATNUM,RATYPE,CB,MSRUN,RRADS,LRADS,
     &  ROCCUP,LOCCUP,REMARK,CNTPLS)

	CHARACTER*80 STRING, REMARK(MAXLIN)
	CHARACTER*80 INTEXT,MSRUN 
	CHARACTER*4 DUMMY(MAXLIN),ATYPE(MAXLIN),ARES(MAXLIN)
	CHARACTER*4 RARES(MAXLIN),RATYPE(MAXLIN),DUMR(MAXLIN)
	CHARACTER*4 SUMNAM(MAXLIN/6)
	CHARACTER*1 CHAIN(MAXLIN), SUMCHAIN(MAXLIN/6)
	CHARACTER*15 CNAME
	CHARACTER*20 DNAME
	INTEGER IRNM(MAXLIN),SUMRES(MAXLIN/6),RIRNM(MAXLIN),INUM(MAXLIN)
	INTEGER ATNUM(MAXLIN), CNTPLS
	REAL CA(MAXLIN,3),RRADS(MAXLIN),LRADS(MAXLIN)
	REAL ROCCUP(*),LOCCUP(*)
	LOGICAL FCUTOFF(MAXLIG)
	DATA NUMWRT /0/

80	FORMAT(A4,2X,I5,2X,A4,A3,I6,4X,3F8.3,2F6.2)
82	FORMAT(A4,2X,I5,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3,2F6.2)
85	FORMAT(A3,I5,1X,'*')
87	FORMAT(A3,I5,1A,1X,'*')
90	FORMAT(A4)

CCC  If not doing an internal filter.  The structure of the data in an
CCC  internal filter is different from a normal filter in that the 'ligands'
CCC  are really just residues of a contiguous structure.  Hence one must 
CCC  have checks for residue termination when writing an internal filter,
CCC  but not when filtering a normal receptor/ligand type construct.
	IF(INTEXT.NE.'Y'.AND.INTEXT.NE.'y')THEN
          DO 200 i=1,NUMLIG
            IF(FCUTOFF(i))THEN
              JCOUNT=I*(ICOUNT+CNTPLS)-(ICOUNT+CNTPLS)+1
              IF(DUMMY(JCOUNT).EQ.'TER ')THEN
                JCOUNT=JCOUNT+1
              ENDIF
*             WRITE(6,*)'JCOUNT IS -',JCOUNT
              
CCC           MSRUN formats output in MS -I format.
              IF(MSRUN .EQ. 'N' .OR. MSRUN .EQ. 'n')THEN
                IF(REMARK(JCOUNT)(1:1).NE.'X')THEN
                  WRITE(15,'(A70)')REMARK(JCOUNT)
                  JCOUNT=JCOUNT+1
                ENDIF
                DO 100 j=JCOUNT,JCOUNT+ICOUNT-1
                  WRITE(15,82)DUMMY(J),INUM(J),ATYPE(J),
     &                 ARES(J),CHAIN(J),IRNM(J),(CA(J,K),K=1,3),
     &                 LOCCUP(J),LRADS(J)
 100            CONTINUE
                WRITE(15,90)DUMMY(J)
              ELSE
                WRITE(15,87)ARES(JCOUNT),INUM(JCOUNT),CHAIN(JCOUNT)
              ENDIF
              NUMWRT=NUMWRT+1
            ENDIF
 200      CONTINUE
          
CCC       If FCUTOFF(i) is .TRUE., then write ligand orientation i.
	ELSE
          DO 275 i=1,NUMLIG
            IF(FCUTOFF(I))THEN
              IF(MSRUN .EQ. 'N' .OR. MSRUN .EQ. 'n')THEN
                DO 250 J=1,IL
CCC               SUMRES checks for proper residue,
CCC               since residue length is variable.
                  IF(IRNM(J).EQ.SUMRES(I))THEN
                    WRITE(15,82)DUMMY(J),INUM(J),
     &                   ATYPE(J),ARES(J),CHAIN(J),IRNM(J),
     &                   (CA(J,K),K=1,3),LOCCUP(J),LRADS(J)
                  ELSEIF(IRNM(J).GT.SUMRES(I))THEN
                    GOTO 275
                  ENDIF
 250            CONTINUE
              ELSE
                WRITE(15,87)SUMNAM(I),SUMRES(I),SUMCHAIN(I)
              ENDIF
              NUMWRT=NUMWRT+1
            ENDIF
 275      CONTINUE
	ENDIF

CCC  If FCUTOFF(i) is .TRUE., then write ligand orientation i.  SUMRES
CCC  checks for proper residue match.

	WRITE(6,*)'TOTAL NUMBER OF ORIENTATIONS WRITEN IS -',NUMWRT
    	RETURN
	END

*************************************************************************
CCC This subroutine of the main FILTER1 filters a receptor file
CCC against itself.  Its purpose in the form it was orriginaly written in
CCC was essentially to find all residues which were within a certain 
CCC radius of another residue specified in a seperate file.  The orriginal
CCC algorithm assumes that only one atom need be within the radius for the
CCC whole receptor to be included.  This is not always what the user will
CCC want and should be modified.
CCC BKS, OCT 87.
*************************************************************************
	SUBROUTINE INTERNAL(MAXATM,CUT,MAXLIG,LATNUM,m,IL,JCUT,n,IRNM,CUTOFF,
     &DUMMY,ATYPE,CA,CB,IR,MAXLIN,INUM,ARES,CHAIN,ICOUNT,JLIG,SUMRES,
     &SUMNAM,SUMCHAIN)

	CHARACTER*4 LATNUM(MAXATM),DUM,DUMMY(MAXLIN),ATYPE(MAXLIN)
	CHARACTER*4 ARES(MAXLIN),SUMNAM(MAXLIN/6)
	CHARACTER*1 RESATM,CHAIN(MAXLIN),SUMCHAIN(MAXLIN/6)
	REAL        CUT(MAXATM), CA(MAXLIN,3), CB(MAXLIN,3) 
	INTEGER     IRNM(MAXLIN),SUMRES(MAXLIN/6),n,m,INUM(MAXLIN)
	LOGICAL     CUTOFF(MAXLIG,2,MAXATM/2)

CCC  JLIG is a residue pointer used for cutoff.  JOLIG stores the previous 
CCC  value for JLIG.
	JLIG=1
	JOLIG=1
	JNLIG=1

	DO 300 I = 1, IL
CCC This IF loop means that it only takes one atom of a residue to be
CCC within CUT for the whole residue to be included.  Once one atom meets
CCC distance criteria, I is advanced until the next residue is being looked
CCC at.  JOLIG is then reset to JLIG... and so on.  JNLIG is a temporary
CCC variable for the CUTOFF pointer JLIG.
          IF (JOLIG .NE. JNLIG) THEN
            IF (IRNM(I) .EQ. IRNM(I-1)) THEN
              IF (RESATM .EQ. 'R') GOTO 300
            ELSE
              JLIG = JLIG + 1
              JOLIG = JNLIG
*             WRITE(6,*) 'JOLIG IS -', JOLIG
            ENDIF
          ENDIF

CCC       Check for non TER line in the PDB.  Won't catch comments.
          IF (DUMMY(I) .NE. ' ') THEN 
            DO 200 J = 1, IR 
CCC           Distance calculation performed here.
              SUM = 0
              DO 75 ISUM = 1, 3
                SUM = SUM + (CA(I,ISUM) - CB(J,ISUM))**2
 75           CONTINUE

*             WRITE (6,*)SUM, CUT(j)
              IF (SUM .LT. CUT(JCUT)) THEN
                CUTOFF(JLIG,n,m) = .TRUE.
					
CCC             SUMRES and SUMNAM are necessary because
CCC             residues can be of variable length.  Used
CCC             in NEWFILL
                SUMRES(JLIG) = IRNM(I)
                SUMNAM(JLIG) = ARES(I)
                SUMCHAIN(JLIG) = CHAIN(I)
                JNLIG = JLIG + 1
                GOTO 300
              ELSEIF ( .NOT. CUTOFF(JLIG,N,M)) THEN
                CUTOFF(JLIG,n,m) = .FALSE.
                SUMRES(JLIG) = IRNM(I)
              ENDIF
 200        CONTINUE
          ENDIF
 300    CONTINUE
CCC
CCC
	RETURN
	END
**************************************************************************
CCC This subroutine is responsible for writing what are essentially ligand
CCC ligand matches to a collated file.  Recall that in NAMFIL what one is
CCC essentially doing is filtering a ligand file against a ligand file.  
CCC The only real trickiness involved in doing this is not writing a ligand
CCC which may match to more than one 'receptor' ligand more than one time.
**************************************************************************
	SUBROUTINE NANFILL(NUMLIG,FCUTOFF,MAXLIG,BCUTOFF,IRNM,DUMMY,ATYPE,
     &  CA,MAXLIN,INUM,ARES,ICOUNT,CNAME,INTEXT,SUMRES,IL,DUMR,RARES,
     &  RIRNM,ATNUM,RATYPE,CB,RECCNT,RCOUNT,SUMWRT,MAXATM,
     &  LATNUM,RRADS,LRADS,ROCCUP,LOCCUP)
CCC
CCC
	CHARACTER*80 INTEXT
	CHARACTER*4 DUMMY(MAXLIN),ATYPE(MAXLIN),ARES(MAXLIN)
	CHARACTER*4 RARES(MAXLIN),RATYPE(MAXLIN),DUMR(MAXLIN),LATNUM(MAXATM)
	CHARACTER*15 CNAME
	CHARACTER*20 DNAME
	INTEGER IRNM(MAXLIN),SUMRES(MAXLIN/6),RIRNM(MAXLIN)
	INTEGER ATNUM(MAXLIN),RECCNT,RCOUNT,RWRITE,INUM(MAXLIN)
	REAL CA(MAXLIN,3),CB(MAXLIN,3),SUMWRT(MAXLIG,2,MAXATM/2)
	REAL RRADS(MAXLIN),LRADS(MAXLIN)
	REAL ROCCUP(*),LOCCUP(*)
	LOGICAL FCUTOFF(MAXLIG),RECWRT

CCC  RECWRT is a logical that keeps track whether the 'receptor' orientation
CCC  has been writen or not.  
	RECWRT=.TRUE.

10	FORMAT('***************************************************')
80	FORMAT(A4,2X,I5,2X,A4,A3,I6,4X,3F8.3,2F6.2)
85	FORMAT(A4,2X,I5,2X,A4,A3,I6,4X,3F8.3,5X,'xxxxx',X,F6.2)
87	FORMAT(A4,2X,I5,2X,A4,A3,I6,4X,3F8.3,2X,F8.3,F6.2)
90	FORMAT('TER')
	DO 200 I=1,NUMLIG

CCC       If FCUTOFF(i) is .TRUE., then write ligand orientation i.
		IF(FCUTOFF(i))THEN
			IF(RECWRT)THEN
				RECWRT=.FALSE.
	               WRITE(15,10)
				DO 95 RWRITE=RECCNT+1,RECCNT+RCOUNT
					WRITE(15,85)DUMR(RWRITE),ATNUM(RWRITE),
     & 				RATYPE(RWRITE),RARES(RWRITE),RIRNM(RWRITE),
     & 				(CB(RWRITE,K),K=1,3),RRADS(RWRITE)
95				CONTINUE
				WRITE(15,90)
			ENDIF

			JCOUNT=((I-1)*(ICOUNT+1))+2
			M=1
			DO 100 j=JCOUNT,JCOUNT+ICOUNT-1
CCC                 print contact distances along w/coords.
				IF(LATNUM(M).EQ.ATYPE(J))THEN
				    WRITE(15,87)DUMMY(J),INUM(J),ATYPE(J),ARES(J),
     & 			    IRNM(J),(CA(J,K),K=1,3),
     &			    SQRT(SUMWRT(I,1,M)),LRADS(J)
                        M=M+1
                    ELSE
                        WRITE(15,87)DUMMY(J),INUM(J),ATYPE(J),ARES(J),
     &                  IRNM(J),(CA(J,K),K=1,3),LOCCUP(J),LRADS(J)        
                    ENDIF
  100          CONTINUE
			WRITE(15,90)
		ENDIF
200	CONTINUE
CCC
	RETURN
	END
