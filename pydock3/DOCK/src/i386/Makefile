# Dock linux Makefile
MV = /bin/mv -f
RM = /bin/rm -f
CP = /bin/cp -f

# Set this to the default target binary to be created. 32b or 64b
# It can be overridden with 'make SIZE=xx'
DEFAULT_SIZE = 64
#DEFAULT_COMPILER = pgf
DEFAULT_COMPILER = gnu

#Set default compiler
ifndef COMPILER
  COMPILER = $(DEFAULT_COMPILER)
endif

#set binary size
ifndef SIZE
  SIZE = $(DEFAULT_SIZE)
endif

#To make dock, simply type 'make' on the commandline
all: dock
	$(shell touch $(ALLOBJSSIZEMODDEP))

#include compiler specific variables
include Makefile.gnu
include Makefile.pg
include Makefile.ifort

# Only change stuff below this line if you know what you are doing 
##############################################################################

#To make dock, simply type 'make' on the commandline
all: dock

ifeq ($(SIZE),64)
  ARCH = $(ARCH64)
  OPT += $(OPT64) 
  OPTPROF += $(OPT64)
else 
  ARCH = $(ARCH32)
  OPT += $(OPT32)
  OPTPROF += $(OPT32)
endif


# If user wants to make dock, which shall we make?
# Make the highest speed version!
dock:	dock_nofb

# To generate various profile outputs
profile: pgprofile gprofile

# For use during code development
dev:  profile

COMPDATE = $(shell date '+%Y%m%d')
#SVNVERZ = $(shell svn info Makefile ../*.f | grep Revision | sed 's/Revision: //g' | sort -n | tail -1 | sed 's/Revision: //g')
#modify this line to add a level when you branch - rgc
LIB_DIR = -L../libfgz $(LDLIBS)
DOCK_RUN_DIR = ../dock_feedback_run/
#DOCK_RUN_DIR = dock_feedback_run/


#following line magically gets every object. 
ALLOBJS = $(shell \ls -1 ../*.f | sed 's/\.\.\///g' | sed 's/\.f/\.o /g' | grep -v code4win | grep -v dock.o )

ALLOBJSI = $(patsubst %.o,%_inst$(SIZE).o,$(ALLOBJS))
ALLOBJSSIZE = $(patsubst %.o,%$(SIZE).o,$(ALLOBJS))
ifdef DEBUG
  ALLOBJSSIZEDEBUG = $(patsubst %.o,%debug.o,$(ALLOBJSSIZE))
  ALLOBJSSIZEMODDEP = $(patsubst %.o,%debug.moddep,$(ALLOBJSSIZE))
  DOCKDEBUGOBJ = dock$(SIZE)debug.o
  DOCKDEBUG = dock$(SIZE)d
else
  ALLOBJSSIZEDEBUG = $(patsubst %.o,%nobug.o,$(ALLOBJSSIZE))
  ALLOBJSSIZEMODDEP = $(patsubst %.o,%nobug.moddep,$(ALLOBJSSIZE))
  DOCKDEBUGOBJ = dock$(SIZE)nobug.o
  DOCKDEBUG = dock$(SIZE)
endif
ALLOBJS_FBOPT = $(patsubst %.o,%_fbopt$(SIZE).o,$(ALLOBJS))
ALLOBJS_PROF = $(patsubst %.o,%_prof$(SIZE).o,$(ALLOBJS))
ALLOBJS_GPROF = $(patsubst %.o,%_gprof$(SIZE).o,$(ALLOBJS))
ALLOBJS_VTUNE = $(patsubst %.o,%_vtune$(SIZE).o,$(ALLOBJS))
ALLOBJS_VTUNEI = $(patsubst %.o,%_vtune_inst$(SIZE).o,$(ALLOBJS))
ALLOBJS_VTUNE_FBOPT = $(patsubst %.o,%_vtune_fbopt$(SIZE).o,$(ALLOBJS))

#The setting of FFLAGS commented out below isnt quite right as the profile and 
#vtune versions wont have the correct values, but for the live production version 
#it will be correct.
#fflags not used at the moment
#	sed "s/\((FFLAGS='\).\+')/\1$(FFLAGS)')/" 
version:
	@echo Setting compile date in ../version.f to $(COMPDATE)
	@echo $(SVNVERZ) is the svn revision number
	@sed "s/\(COMPDATE='\)[0-9]\+'/\1$(COMPDATE)'/" ../version.f.template | \
	sed "s/\(ARCHSIZE='\)[0-9]\+'/\1$(SIZE)'/" | \
	sed "s/\(SVNVERSION='\).*[0-9]\+'/\1$(SVNVERZ)'/" > ../version.f 

#Target w/o feedback compilation optimization
dock_nofb: $(ALLOBJSSIZEDEBUG) version $(DOCKDEBUGOBJ) 
	@echo Compiling new date into dock.o
	$(F77) -c $(FFLAGS) ../dock.f -o $(DOCKDEBUGOBJ)
	@echo Linking all objects and creating dock
#	@$(F77) $(FFLAGS) $(DOCKDEBUGOBJ) $(ALLOBJSSIZEDEBUG) $(LIB_DIR) $(FFLIBS) -o $(DOCKDEBUG)
	$(F77) $(FFLAGS) $(LDFLAGS) $(DOCKDEBUGOBJ) $(ALLOBJSSIZEDEBUG) $(LIB_DIR) $(FFLIBS) -o $(DOCKDEBUG)


#Target w/o feedback compilation optimization for Intel VTune
#No Static linking for call graph
dock_vtune: $(ALLOBJS_VTUNE) version dock_vtune$(SIZE).o
	@echo Compiling new date into dock.o
	$(F77) -c $(FFLAGS_VTUNE) ../dock.f -o dock_vtune$(SIZE).o
	@echo Linking all objects and creating dock
	@$(F77) $(FFLAGS_VTUNE) dock_vtune$(SIZE).o $(ALLOBJS_VTUNE) $(LIB_DIR) $(FFLIBS) -o dock_vtune$(SIZE)

dock_vtune_instrumented$(SIZE):  $(ALLOBJS_VTUNEI) version
	@echo Compiling new date into dock.o
	$(F77) -c $(FFLAGS_VTUNE) $(FBGEN) ../dock.f -o dock_vtune_inst$(SIZE).o
	@echo Linking all objects and creating final version of dock 
	$(F77) $(FFLAGS_VTUNE) $(FBGEN) dock_vtune_inst$(SIZE).o $(ALLOBJS_VTUNEI)  $(LIB_DIR) $(FFLIBS)  -o dock_vtune_instrumented$(SIZE)
#	$(F77) $(FFLAGS_VTUNE) $(FBGEN) $(LDFLAGS) dock_vtune_inst$(SIZE).o $(ALLOBJS_VTUNEI)  $(LIB_DIR) $(FFLIBS)  -o dock_vtune_instrumented$(SIZE)

# Generate feedback instrumentation file
vtune_pgfi$(SIZE).out: dock_vtune_instrumented$(SIZE)
	@echo Current number of molecules for this execution:
	@grep molecules_maximum $(DOCK_RUN_DIR)/INDOCK
	@echo
	cd $(DOCK_RUN_DIR) ; ../i386/dock_vtune_instrumented$(SIZE)
	@$(MV) $(DOCK_RUN_DIR)/pgfi.out ./pgfi_vtune$(SIZE).out


dock_vtune_fb: vtune_pgfi$(SIZE).out copy_vtune_fb_output $(ALLOBJS_VTUNE_FBOPT) version dock_vtune_fbopt$(SIZE).o
	@echo Compiling new date into dock_fbout$(SIZE).o
	$(F77) -c $(FFLAGS_VTUNE) $(FBOPT) ../dock.f -o dock_vtune_fbopt$(SIZE).o
	@echo Linking all objects and creating final version of dock 
#	@$(F77) $(FFLAGS_VTUNE) $(FBOPT) $(LDFLAGS) dock_vtune_fbopt$(SIZE).o $(ALLOBJS_VTUNE_FBOPT) $(LIB_DIR) $(FFLIBS) -o dock$(SIZE)
	@$(F77) $(FFLAGS_VTUNE) $(FBOPT)  dock_vtune_fbopt$(SIZE).o $(ALLOBJS_VTUNE_FBOPT) $(LIB_DIR) $(FFLIBS) -o dock_vtune$(SIZE)
	@$(MAKE) fb_cleanup


#Target w/o feedback compilation for profiling for GNU gprof
dock_gprof$(SIZE): $(ALLOBJS_GPROF) version dock_gprof$(SIZE).o
	@echo Compiling new date into dock.o
	$(F77) -c $(FFLAGS_GPROF) ../dock.f -o dock_gprof$(SIZE).o
	@echo Linking all objects and creating dock
#	@$(F77) $(FFLAGS_GPROF)  dock_gprof$(SIZE).o $(ALLOBJS_GPROF) $(LIB_DIR) $(FFLIBS) -o dock_gprof$(SIZE)
	@$(F77) $(FFLAGS_GPROF) $(LDFLAGS) dock_gprof$(SIZE).o $(ALLOBJS_GPROF) $(LIB_DIR) $(FFLIBS) -o dock_gprof$(SIZE)

#gprof
dock_prof$(SIZE): $(ALLOBJS_PROF) version dock_prof$(SIZE).o
	@echo Compiling new date into dock.o
	$(F77) -c $(FFLAGS_PROF) ../dock.f -o dock_prof$(SIZE).o
	@echo Linking all objects and creating dock
	@$(F77) $(FFLAGS_PROF) $(LDFLAGS) dock_prof$(SIZE).o $(ALLOBJS_PROF) $(LIB_DIR) $(FFLIBS) -o dock_prof$(SIZE)
#	@$(F77) $(FFLAGS_PROF) dock_prof$(SIZE).o $(ALLOBJS_PROF) $(LIB_DIR) $(FFLIBS) -o dock_prof$(SIZE)

# Compile a feedback instrumented version
dock_instrumented$(SIZE):  $(ALLOBJSI) version
	@echo Compiling new date into dock.o
	$(F77) -c $(FFLAGS) $(FBGEN) ../dock.f -o dock_inst$(SIZE).o
	@echo Linking all objects and creating final version of dock 
#	$(F77) $(FFLAGS) $(FBGEN) dock_inst$(SIZE).o $(ALLOBJSI)  $(LIB_DIR) $(FFLIBS)  -o dock_instrumented$(SIZE)
	$(F77) $(FFLAGS) $(FBGEN) $(LDFLAGS) dock_inst$(SIZE).o $(ALLOBJSI)  $(LIB_DIR) $(FFLIBS)  -o dock_instrumented$(SIZE)

# Generate feedback instrumentation file
$(FB_FNAME): dock_instrumented$(SIZE)
	@echo
	@echo Running instrumented DOCK to generate feedback data ... Please wait.
	@echo As of 11/12/2009, run was 25 sec on sgehead for molecules_maximum=200
	@echo If runtime is too long, try modifying molecules_maximum in the INDOCK file.
	@echo 
	@echo Current number of molecules for this execution:
	@grep molecules_maximum $(DOCK_RUN_DIR)/INDOCK
	@echo
	cd $(DOCK_RUN_DIR) ; ../i386/dock_instrumented$(SIZE)
	$(MAKE) move_fb_output

# Target for Feedback optimized version
#dock_fb: $(FB_FNAME)  $(ALLOBJS_FBOPT) version dock_fbopt$(SIZE).o
dock_fb: $(FB_FNAME) copy_fb_output $(ALLOBJS_FBOPT) version dock_fbopt$(SIZE).o
	@echo Compiling new date into dock_fbout$(SIZE).o
	$(F77) -c $(FFLAGS) $(FBOPT) ../dock.f -o dock_fbopt$(SIZE).o
	@echo Linking all objects and creating final version of dock 
	@$(F77) $(FFLAGS) $(FBOPT) $(LDFLAGS) dock_fbopt$(SIZE).o $(ALLOBJS_FBOPT) $(LIB_DIR) $(FFLIBS) -o dock$(SIZE)
#	@$(F77) $(FFLAGS) $(FBOPT)  dock_fbopt$(SIZE).o $(ALLOBJS_FBOPT) $(LIB_DIR) $(FFLIBS) -o dock$(SIZE)
	@$(MAKE) fb_cleanup

# Profile 
# Intel VTune
vprofile:  dock_vtune$(SIZE)
#	vtl activity -d 30 -c sampling -c callgraph -master sampling -app ../i386/dock_vtune -moi ../i386/dock_vtune run
#  -master sampling may give segfault from tool  Not sure of the exact cause yet.
	cp dock_vtune$(SIZE)  $(DOCK_RUN_DIR)
	cd $(DOCK_RUN_DIR) ; vtl activity -d 100 -c sampling -of event.cfg  -app ./dock_vtune$(SIZE) -moi ./dock_vtune$(SIZE) run 
	cd $(DOCK_RUN_DIR) ; vtl view -gui &

# Profile 
# Portland Group Profiler
pgprofile:  pgprof.out
# Not sure of a better way to send commands to stdin
	@echo print > pgcommand.txt ; echo quit >> pgcommand.txt; 
	$(PROF) -text -exe dock_prof$(SIZE) < pgcommand.txt > pgprofile.txt
	@$(RM) pgcommand.txt

pgprof.out: dock_prof$(SIZE)
	cd $(DOCK_RUN_DIR) ; ../i386/dock_prof$(SIZE)
	$(MV) $(DOCK_RUN_DIR)/pgprof.out ./

# GNU gprof
# (The generated gprofile.txt has a useful call graph at the bottom, for those looking to learn the dock code)
gprofile:  gmon.out
	$(GPROF) dock_gprof$(SIZE) > gprofile.txt 

gmon.out: dock_gprof$(SIZE)
	cd $(DOCK_RUN_DIR) ; ../i386/dock_gprof$(SIZE)
	$(MV) $(DOCK_RUN_DIR)/gmon.out ./

# symlink the built executable files into the docking/DOCK/bin directory
install:
	find * -type f ! -iname "Make*" -executable -exec ln -s ../src/i386/{} ../../bin/{} \;

clean:
	$(RM) *.moddep *.o *.oo *.ipa *.disasm *~ \
	dock32 dock_instrumented32 dock_gprof32 dock_prof32 dock_vtune32 \
	dock64 dock_instrumented64 dock_gprof64 dock_prof64 dock_vtune64 \
        dock32d dock64d \
	dock_vtune_instrumented32  dock_vtune_instrumented64 \
	pgprof.out pgprofile.txt gprofile.txt gmon.out pgprof.out \
	*.mod
# Not sure why we build after clean
#	$(MAKE) clean_fb_output

# ##################################################
# How to build objects ...
# ##################################################

#auto dependency generation, gets all modules necessary
%$(SIZE)debug.moddep: ../%.f
	@echo -n $*.f $*$(SIZE)debug.moddep ": " > $*$(SIZE)debug.moddep
	@grep "^      use " ../$*.f | awk '{print $$2"$(SIZE)debug.o "$$2".mod "}' | paste -s -d ' ' >> $*$(SIZE)debug.moddep

#Regular .o without instrumentation of any kind.  For 'dock' target
%$(SIZE)debug.o %.mod: ../%.f 
	$(F77) -c $(FFLAGS) ../$*.f -o $*$(SIZE)debug.o
	@objdump -d $*$(SIZE)debug.o > $*$(SIZE)debug.disasm

%$(SIZE)debug.o: ../%.c
	$(CC) -c $(CFLAGS) ../$*.c -o $*$(SIZE)debug.o
	@objdump -d $*$(SIZE)debug.o > $*$(SIZE)debug.disasm

#auto dependency generation, gets all modules necessary
%$(SIZE)nobug.moddep: ../%.f
	@echo -n $*.f $*$(SIZE)nobug.moddep ": " > $*$(SIZE)nobug.moddep
	@grep "^      use " ../$*.f | awk '{print $$2"$(SIZE)nobug.o "$$2".mod "}' | paste -s -d ' ' >> $*$(SIZE)nobug.moddep

%$(SIZE)nobug.o %.mod: ../%.f
	$(F77) -c $(FFLAGS) ../$*.f -o $*$(SIZE)nobug.o
	@objdump -d $*$(SIZE)nobug.o > $*$(SIZE)nobug.disasm

%$(SIZE)nobug.o: ../%.c
	$(CC) -c $(CFLAGS) ../$*.c -o $*$(SIZE)nobug.o
	@objdump -d $*$(SIZE)nobug.o > $*$(SIZE)nobug.disasm

#Feedback instrumented object files, for 'dock_instrumented' target and therefore 'dock_fb' target.
%_inst$(SIZE).o: ../%.f
	$(F77) -c $(FBGEN) $(FFLAGS) ../$*.f -o $*_inst$(SIZE).o

%_inst$(SIZE).o: ../%.c
	$(CC) -c  $(FBGEN) $(CFLAGS) ../$*.c -o $*_inst$(SIZE).o


#Feedback optimized object files, for 'dock_fb' target.
%_fbopt$(SIZE).o: ../%.f
	$(F77) -c  $(FBOPT) $(FFLAGS) ../$*.f -o $*_fbopt$(SIZE).o
	@objdump -d $*_fbopt$(SIZE).o > $*_fbopt$(SIZE).disasm

%_fbopt$(SIZE).o: ../%.c
	$(F77) -c  $(FBOPT) $(CFLAGS) ../$*.c -o $*_fbopt$(SIZE).o
	@objdump -d $*_fbopt$(SIZE).o > $*_fbopt$(SIZE).disasm


# Instrumented object files for profiling 

# Target 'profile', portland group profiler
%_prof$(SIZE).o: ../%.f
	$(F77) -c  $(FFLAGS_PROF) ../$*.f -o $*_prof$(SIZE).o
# Target 'gprofile', GNU gprof profiler
%_gprof$(SIZE).o: ../%.f
	$(F77) -c  $(FFLAGS_GPROF) ../$*.f -o $*_gprof$(SIZE).o
# Target 'vprofile', Intel Vtune
%_vtune$(SIZE).o: ../%.f
	$(F77) -c  $(FFLAGS_VTUNE) ../$*.f -o $*_vtune$(SIZE).o
	@objdump -d $*_vtune$(SIZE).o > $*_vtune$(SIZE).disasm
%_vtune$(SIZE).o: ../%.c
	$(CC) -c  $(CFLAGS_VTUNE) ../$*.c -o $*_vtune$(SIZE).o
	@objdump -d $*_vtune$(SIZE).o > $*_vtune$(SIZE).disasm
%_vtune_inst$(SIZE).o: ../%.f
	$(F77) -c $(FBGEN) $(FFLAGS_VTUNE) ../$*.f -o $*_vtune_inst$(SIZE).o
	@objdump -d $*_vtune_inst$(SIZE).o > $*_vtune_inst$(SIZE).disasm
%_vtune_inst$(SIZE).o: ../%.c
	$(CC) -c  $(FBGEN) $(CFLAGS_VTUNE) ../$*.c -o $*_vtune_inst$(SIZE).o
	@objdump -d $*_vtune_inst$(SIZE).o > $*_vtune_inst$(SIZE).disasm
%_vtune_fbopt$(SIZE).o: ../%.f
	$(F77) -c  $(FBOPT) $(FFLAGS_VTUNE) ../$*.f -o $*_vtune_fbopt$(SIZE).o
	@objdump -d $*_vtune_fbopt$(SIZE).o > $*_vtune_fbopt$(SIZE).disasm
%_vtune_fbopt$(SIZE).o: ../%.c
	$(F77) -c  $(FBOPT) $(CFLAGS_VTUNE) ../$*.c -o $*_vtune_fbopt$(SIZE).o
	@objdump -d $*_vtune_fbopt$(SIZE).o > $*_vtune_fbopt$(SIZE).disasm


############################################################################
# The following targets can be used to determine which is the fastest
# architecture switch to compile with.  Current experiements show no
# real difference other than 32b vs 64b.  This information can be used to
# set ARCH32 and ARCH64, or the binaries can be used directly for specific
# computers.
#
# make compile_all_archs  -  Compiles various dock_<architecture> binaries
# make run_all_archs  - Runs various dock_<architecture> binaries and
#			puts results in dock_feedback_run/arch_results/



# The following 32b architectures wont compile due to lack of support for 
# huge pages.  These are slow older architectures anyway, so skip.
# 	athlon athlonxp k7 p5 p6 piii px-32
TARGET_ARCHS32 = barcelona-32 core2-32 istanbul-32 k8-32 nehalem-32 p6 p7-32 penryn-32 piv shanghai-32 

# The following architecture target doesnt support feedback optimization.
# So skip.
# 	x64
TARGET_ARCHS64 = amd64 amd64e barcelona-64 core2-64 istanbul-64 k8-64 k8-64e nehalem-64 p7-64 penryn-64 px-64 shanghai-64 

compile_all_archs:
	$(foreach A, $(TARGET_ARCHS32), make compile_arch SIZE=32 ARCH='-tp $(A)' AA=$(A);)
	$(foreach A, $(TARGET_ARCHS64), make compile_arch SIZE=64 ARCH='-tp $(A)' AA=$(A);)

run_all_archs:
	mkdir $(DOCK_RUN_DIR)/arch_results
	$(foreach A, $(TARGET_ARCHS32), make run_arch AA=$(A);)
	$(foreach A, $(TARGET_ARCHS64), make run_arch AA=$(A);)
	tail -n 1 $(DOCK_RUN_DIR)/arch_results/OUTDOCK* > summary.txt

clean_all_archs:
	$(foreach A, $(TARGET_ARCHS32), \rm -f dock_$(A);)
	$(foreach A, $(TARGET_ARCHS64), \rm -f dock_$(A);)

# Dont call manually
compile_arch:	clean	
	$(MAKE) dock
	$(MV) dock$(SIZE) dock_$(AA)	

# Dont call manually
run_arch:	
	cd $(DOCK_RUN_DIR) ; ../i386/dock_$(AA)
	$(MV) $(DOCK_RUN_DIR)/OUTDOCK $(DOCK_RUN_DIR)/arch_results/OUTDOCK_$(AA) 

#disable builtin rules or m2c tries to run to make .mod files
.SUFFIXES:
MAKEFLAGS += --no-builtin-rules

-include $(ALLOBJSSIZEMODDEP)



