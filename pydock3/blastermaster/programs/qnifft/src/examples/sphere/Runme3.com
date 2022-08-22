#!/bin/csh -f
#
# file: loop.com
#
# examples of loop contructions in shell scripts
# to run qnifft11
#
# here we run a series of qnifft11 runs where we set a parameter
# (the salt concentration)
# to an explicit list of values: this works for any parameter
# type: real, integer, string
#
# we can either define the list of values in this script file:
#
set values="0.1 0.01 0.001"
#
# or we could get them from a file (called here values.lis):
#
#set values=`cat values.lis`
#
# now we loop through the runs
#
foreach curval ($values)
	echo doing run with salt conc. of $curval M
	#
	# copy template parameter file to temporary file:
	#
	cp parm parm1
	#
	# append current new parameter value to temporary parameter file
	#
      echo "salt=$curval" >> parm1
	#
	# run qnifft22, putting output into new logfile
	#
	qnifft22 parm1 > salt$curval.log
	#
	# or putting output into same logfile
	#
	#	qnifft22 parm1 >> saltall.log
end
# clean up
rm parm1
exit(0)
