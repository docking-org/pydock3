#!/bin/csh -f
alias qnifft11 qnifft22_65
#
# example delphi calculation- shift in pKa of His 40 in rat trypsin
# and its interaction with the rest of the protein
# parm is basic parameter file which gets modified as needed
#
#-------------------------
# reference state energy (isolated His in water
#--------------------------
#
# +his40 alone in water
#
cp parm parm1
echo 'pdb_in=rat_ptb_his40.pdb' >> parm1
echo 'site_in=rat_ptb_his40.pdb' >> parm1
echo 'site_out=tryp0_hisp.fld1' >> parm1
echo 'charge_file=his40.crg' >> parm1
qnifft11 parm1 > tryp0_hisp.log1
#
# check charging assignment done ok
#
chkchg < rat_ptb_his40.atm > chkchg.out
#
# neutral his40 alone in water
#
cp parm parm1
echo 'pdb_in=rat_ptb_his40.pdb' >> parm1
echo 'site_in=rat_ptb_his40.pdb' >> parm1
echo 'site_out=tryp0_hisn.fld1' >> parm1
echo 'charge_file=his40_neut.crg' >> parm1
qnifft11 parm1 > tryp0_hisn.log1
#
# check charging assignment done ok
#
chkchg < rat_ptb_his40.atm >> chkchg.out
#
#
#-------------------------
#  Energy in protein- rest of protein is neutral for
# potential calculation, but charged for interaction calcs
# in this case using 'full' charges i.e. lys, arg +1, glu, asp -1
#--------------------------
#
# +his40
#
cp parm parm1
echo 'pdb_in=rat_ptb.pdb' >> parm1
echo 'site_in=rat_ptb.pdb' >> parm1
echo 'site_out=tryp_hisp.fld1' >> parm1
echo 'charge_file=his40.crg' >> parm1
echo 'site_charge_file=tryp_full.crg' >> parm1
qnifft11 parm1 > tryp_hisp.log1
#
#
# neutral his40 
#
cp parm parm1
echo 'pdb_in=rat_ptb.pdb' >> parm1
echo 'site_in=rat_ptb.pdb' >> parm1
echo 'site_out=tryp_hisn.fld1' >> parm1
echo 'charge_file=his40_neut.crg' >> parm1
echo 'site_charge_file=tryp_full.crg' >> parm1
qnifft11 parm1 > tryp_hisn.log1
#
find:
# pull out pertinent results from log/fld files
grep Interaction *.log1 > trypsin.sum1
# for the interaction between his and rest of tryp charges, program sumfld
# will collect interact energy in kT on a per residue basis 
sumfld tryp_hisp.fld1 > tryp_hisp.int1
sumfld tryp_hisn.fld1 > tryp_hisn.int1
