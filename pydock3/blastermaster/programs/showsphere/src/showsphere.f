        PROGRAM showsphere

C////////////////////////////////////////////////////////////////////////////
C PROGRAMMER:	Stuart Oatley, modified by E. Meng, rewritten by D. Gschwend
C REV.DATE:	22-AUG-94
C HISTORY:	19-AUG-94 v2.00 spaghetti code completely rewritten
C				removal of 3.5 format sphere file headers
C		   JUL-94 v1.02 allowed cluster 0, Diana Roe
C		   MAY-94 v1.01 minor bug fixes & modifications, DAG
C		   ???    v1.00 Stuart Oatley
C INPUT:	sphere cluster file
C OUTPUT:	pdb and MS surface file of spheres
C PURPOSE:	depict sphere site for visualization
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        IMPLICIT NONE

	CHARACTER*80	sphfile,pdbfile,srffile,line
	CHARACTER*1	choice
	INTEGER		i,j,lp
	INTEGER		doclus,cluster,sphnum,nsph
	LOGICAL		dosurf,doall,doZero,found
	REAL		xyz(3),radius

	WRITE(6,*) 'Enter name of sphere cluster file:'
	READ(5,'(A80)') sphfile
	WRITE(6,'($,A)') ' Enter cluster number to process (<0 = all): '
	READ(5,*) doclus
	doall=(doclus.LT.0)
	WRITE(6,'($,A)') ' Generate surfaces as well as pdb files (<N>/Y)? '
	READ(5,'(A1)') choice
	dosurf=((choice.EQ.'y').OR.(choice.EQ.'Y'))

	doZero=(doclus.EQ.0)

	IF (doall) THEN
	    WRITE(6,*) 'Enter name for output file prefix:'
	    READ(5,'(A80)') pdbfile
	    lp=INDEX(pdbfile,' ')-1
	    WRITE(6,'($,A)') ' Process clu0 (ALL spheres) (<N>/Y)? '
	    READ(5,'(A1)') choice
	    doZero=((choice.EQ.'Y').OR.(choice.EQ.'y'))
	ELSE
	    WRITE(6,*) 'Enter name for output PDB file name:'
	    READ(5,'(A80)') pdbfile
	    IF (dosurf) THEN
		WRITE(6,*) 'Enter name for output surface file name:'
		READ(5,'(A80)') srffile
	    END IF
	END IF

	found=.FALSE.

	OPEN(UNIT=1,NAME=sphfile,STATUS='OLD')

C	Scan for the next cluster
1050	READ(1,'(A80)',END=1800) line
	DO WHILE (line(1:7).NE.'cluster')
	    READ(1,'(A80)',END=1800) line
	END DO

C	Process the cluster
1100	READ(line,1901,ERR=1810,END=1800) cluster,nsph
	IF ((doall).OR.(cluster.EQ.doclus)) THEN	! We want this cluster

	    IF ((cluster.EQ.0).AND.(.NOT.doZero)) GOTO 1050

	    WRITE(6,*)
	    WRITE(6,1905) cluster,nsph
	    found=.TRUE.

C	    Silly file naming block
	    IF (doall) THEN
		WRITE(pdbfile(lp+1:lp+2),'(I2)') cluster
		IF (cluster.LT.10) pdbfile(lp+1:lp+1)='_'
		srffile=pdbfile
		pdbfile(lp+3:lp+6)='.pdb'
		srffile(lp+3:lp+5)='.ms'
	    END IF

C	    Open files
	    OPEN(UNIT=2,NAME=pdbfile,STATUS='UNKNOWN')
	    IF (dosurf) OPEN(UNIT=3,NAME=srffile,STATUS='UNKNOWN')

C	    Process each sphere
	    DO i=1,nsph
		READ(1,1902,ERR=1830) sphnum,(xyz(j),j=1,3),radius
		WRITE(6,'($,I4,".")') sphnum
		CALL FLUSH(6)

C		Write sphere in pdb file format
		WRITE(2,1903) sphnum,sphnum,(xyz(j),j=1,3)
		WRITE(2,'("TER")') 

C		Write sphere in MS file format (non-QCPE) if asked for
		IF (dosurf) THEN
		    WRITE(3,1904) sphnum,(xyz(j),j=1,3)
		    CALL sphere(sphnum,xyz(1),xyz(2),xyz(3),radius,cluster)
		END IF

	    END DO
	    CLOSE(3)
	    CLOSE(3)

	    WRITE(6,*)
	    IF (.NOT.doall) GOTO 1800	! done

	END IF

	GOTO 1050	! find next cluster in input

	
	    
1800	CONTINUE	
	IF (.NOT.found) GOTO 1820
	STOP

1810	CONTINUE
	WRITE(6,*)
	WRITE(6,*) 'Error!  Sphere cluster header line has format problem.'
	STOP

1820	CONTINUE
	WRITE(6,*)
	WRITE(6,*) 'Error!  Cluster not found:',doclus
	STOP

1830	CONTINUE
	WRITE(6,*)
	WRITE(6,*) 'Error!  Sphere cluster format in error.'
	STOP

1901	FORMAT(8X,I5,32X,I5)					! cluster line
1902	FORMAT(I5,3F10.5,F8.3)					! sphere read
1903	FORMAT('ATOM  ',X,I4,X,' C  ',' SPH',2X,I4,4X,3F8.3)	! PDB atom
1904	FORMAT('SPH ',I4,'   C',3F9.3,X,'A')			! MS atom
1905	FORMAT('Cluster',I3,':',I4,' spheres.')

	END
c------------------------------------------------------------------------

	SUBROUTINE sphere(sphnum,xs,ys,zs,radius,cluster)

	REAL		pi
	PARAMETER	(pi=3.1415926)

	INTEGER		i,j,nsteps
	INTEGER		sphnum,cluster
	REAL		xs,ys,zs,radius
	REAL		ain,ainc,dinc,ang,x,y,z,r,a,steps

	ainc=pi/9.
	dinc=ainc*radius
	DO i=1,9
	    j=i-5
	    ang=REAL(j)*ainc
	    z=zs+radius*SIN(ang)
	    r=radius*COS(ang)
	    ain=2.*pi
	    IF(ABS(r).GT.0.001) ain=dinc/r
	    steps=2.*pi/ain
	    nsteps=steps+0.5
	    ain=2.*pi/nsteps
	    a=0.
	    DO j=1,nsteps
		x=xs+r*COS(a)
		y=ys+r*SIN(a)
		WRITE(3,2901) sphnum,x,y,z,r
		a=a+ain
	    END DO
	END DO

2901	FORMAT('SPH ',I4,'   C',3F9.3,X,'SC0',F7.3)

	RETURN
	END
