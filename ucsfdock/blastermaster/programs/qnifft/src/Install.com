#!/bin/csh -f
#
#-----------------------------------------------
echo ' '
echo 'Installing QNIFFT '
echo ' '
setenv DELDIR `pwd`
echo delphi root directory set to $DELDIR
echo 'environment variable DELDIR defined'
echo ' '
#
#-----------------------------------------------
if ( -e $DELDIR/bin ) then
  echo $DELDIR/bin already exists
else
  mkdir $DELDIR/bin
endif
#set path=( $path $DELDIR/bin )
echo ' '
echo path now set to $path
echo ' '
#-----------------------------------------------
echo 'platform: SGI (1) GNU (2) Compaq (3)  Lahey Compiler (4) Other (5)?'
set  platform=$<
if ( $platform == 1 )then
 set compile = 'make -n -f $DELDIR/source/Makefile_sgi -f Makefile'
else if ( $platform == 2 )then
 set compile = 'make -n -f $DELDIR/source/Makefile_gnu -f Makefile'
else if ( $platform == 3 )then
 set compile = 'make -n -f $DELDIR/source/Makefile_compaq -f Makefile'
else if ( $platform == 4 )then
 set compile = 'make -n -f $DELDIR/source/Makefile_lahey -f Makefile'
else
  if ( -e $DELDIR/source/Makefile_other ) then
   set compile = 'make -n -f $DELDIR/source/Makefile_other -f Makefile'
  else
    echo 'You need to create ' $DELDIR/source/Makefile_other
    echo from $DELDIR/source/Makefile_blank to match your compiler/platform specs then compile
    exit
  endif
endif
echo ' '
echo $compile
#-----------------------------------------------
echo ' '
echo compiling utility programs
echo ' '
cd $DELDIR/source/utility
touch *.f
$compile sumfld
$compile ionintm
$compile addmap
$compile chkchg
#make qequil_pdb
echo ' '

#-----------------------------------------------
echo compiling main program
echo ' '
cd $DELDIR/source/qnifft/v2.2
touch *.h
$compile qnifft_small
$compile qnifft_medium
$compile qnifft_large
$compile qnifft_huge
echo ' '

#-----------------------------------------------
echo List of compiled programs:
ls -lt $DELDIR/bin
echo ' '
#-----------------------------------------------

cd $DELDIR
echo ' '
echo 'To run the programs in future place these 2 commands in your log file:'
echo 'setenv DELDIR' $DELDIR
echo 'set path=( $path' $DELDIR'/bin )'
echo ' '
echo ' '
echo 'Finished Installing QNIFFT - run the examples, and go'
echo ' '


exit
