
DMSROOT = $(CURDIR)/../..

# LIBDIR is the directory where a 'dms' subdirectory will be created,
# in which the server binary (dmsd) and various support files will be installed
LIBDIR  = $(DMSROOT)/bin

# BINDIR is the directory where the dms binary will be installed
BINDIR  = $(DMSROOT)/bin

# PDBINC is the directory where the 'pdb.h' file from the UCSF PDB-reading
# library got installed.  'pdb.h' is not installed by default by that
# distribution, so you will have to install it by hand into PDBINC
PDBINC  = libpdb

# PDBLIB is the directory where libpdb.a from the UCSF PDB-reading library 
# got installed.
PDBLIB  = libpdb

# MANDIR is the directory where the dms.1 manual page will be installed
MANDIR  = $(LIBDIR)

ifeq ($(shell test -e /bin/cygwin1.dll && echo found),found)
	PROG	= dms.exe
else
	PROG	= dms
endif

# you should not need to modify any of the following
DESTLIB	= ${LIBDIR}/
SERVER	= ${DESTLIB}dmsd
#DEFS	= -DDESTLIB=\"${DESTLIB}\" -DSERVER_PATH=\"$(SERVER)\"  # Removed
DEFS	= -DDESTLIB=\"${DESTLIB}\"
OPT	= -O3 -static
CFLAGS  = ${DEFS} -I${PDBINC} ${OPT} 

OBJS	= compute.o fwritev.o input.o ms.o output.o emalloc.o tokenize.o dms_paths.o

all:	${PROG}
	cd dmsd ; ${MAKE} LIBDIR=${LIBDIR} OPT="${OPT}"

${PROG}: pdb ${OBJS}
	${CC} -L${PDBLIB} -o ${PROG} ${OBJS} -lpdb

pdb:
	cd libpdb ; ${MAKE} OPT="${OPT}"

install: ${PROG} prepinstall maninstall afterinstall
afterinstall: realinstall
realinstall: beforeinstall
	cp ${PROG} ${BINDIR}
	-strip ${BINDIR}/${PROG}

prepinstall:
	if [ -d ${DESTLIB} ]; then true ; else mkdir -pv ${DESTLIB}; fi

maninstall:
	-rm -f ${MANDIR}/dms.1
	cp dms.1 ${MANDIR}

beforeinstall:
	if [ -d ${DESTLIB} ]; then true ; else mkdir -pv ${DESTLIB}; fi
	cd dmsd ; ${MAKE} LIBDIR=${LIBDIR} install

afterinstall:
	cp radii.proto ${DESTLIB}radii

clean:
	/bin/rm -f ${OBJS}
	cd libpdb ; ${MAKE} clean
	cd dmsd ; ${MAKE} clean

spotless:
	/bin/rm -f ${OBJS} ${PROG}
	cd libpdb ; ${MAKE} spotless
	cd dmsd ; ${MAKE} spotless
