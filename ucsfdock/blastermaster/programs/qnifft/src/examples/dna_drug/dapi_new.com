#!/bin/csh -f
#
# file: dapi.com
#
# loop contructions in shell scripts to run qnifft on dapi/dna
# drug interaction- note how power of unix script looping is 
# combined with flexibile qnifft input to do all of runs
# ionint.com runs analysis program ionint, to integrate
# ion atmosphere terms for nested focussing potential
# maps- this is required to obtain energies in the non-linear
# pb equation.
#
#
alias qnifft qnifft22_lnx65
alias qnifft qnifft22_193
/bin/rm dapout.log6 dapin.log5
set molc="dap.pdb dna12.pdb dapdna12.pdb"
foreach cmolc ($molc)
	echo doing run with molc $cmolc
	set salt="0.0 0.1"
	foreach csalt ($salt)
		echo doing run with salt conc. of $csalt M
		#
		echo outer run
		#
		cp parm parm2
		echo "boundary=zero" >> parm2
      	echo "salt=$csalt" >> parm2
      	echo "ani_mon=$csalt" >> parm2
      	echo "cat_mon=$csalt" >> parm2
		echo "pdb_in=$cmolc" >> parm2
		echo "phi_out=dapout.phi" >> parm2
		echo "dielectric=dapout.eps" >> parm2
		echo "scale=0.9" >> parm2
		qnifft parm2 >> dapout.log6
		#
		echo inner run
		#
		cp parm parm2
		echo "boundary=focus" >> parm2
      	echo "salt=$csalt" >> parm2
      	echo "ani_mon=$csalt" >> parm2
      	echo "cat_mon=$csalt" >> parm2
		echo "pdb_in=$cmolc" >> parm2
		echo "phi_in=dapout.phi" >> parm2
		echo "phi_out=dapin.phi" >> parm2
		echo "dielectric=dapin.eps" >> parm2
		echo "scale=3.0" >> parm2
		qnifft parm2 >> dapin.log6
#		ionint.com $csalt >> dapin.log6
		ionintm.com $csalt >> dapin.log6
	end
end
#
# check charge assignment was done ok
#
chkchg < temp.atm > chkchg.out
grepit.com dapin.log6 > dapi.sum6
# clean up
#rm parm2
exit(0)
