#
# run example with new delphi input
#
# set environment variable  to point to where default parameter file
# is (which is in this directory for the purposes of this example). 
# normally this would be done
# in your login file, and delphi.def would not here but in some suitable
# directory
#
setenv DELDIR .
#
# run program with existing parameter file
#
qnifft22 parm > sph1.log
#
# alter some parameter and rerun program, putting output into
# a different potential file. 
# one way is to prepare a new parameter file with the editor
# alternatively we can do everthing in a script file:
# we dont want to touch parm so
# we first copy it to a new file, say parm1, and put new line
# at end of parm1. any parameters can be changed this way,
# Thus we can easily run through a series of parameters
# using only one initial parameter file, and a series of script commands
#
# this example shows a focussing run
#
cp parm parm1
echo 'boundary=focus' >> parm1
echo 'grid=65' >> parm1
echo 'fill=20' >> parm1
echo 'phi_in=sph.phi' >> parm1
echo 'phi_out=sph3.phi' >> parm1
qnifft22 parm1 > sph3.log
