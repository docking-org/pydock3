#!/bin/csh -f
#
# file: dapi.com
#
alias qnifft11 qnifft12_r8
# loop contructions in shell scripts to run qnifft on dapi/dna
# drug interaction- note how power of unix script looping is 
# combined with flexibile qnifft input to do all of runs
# ionint.com runs analysis program ionint, to integrate
# ion atmosphere terms for nested focussing potential
# maps- this is required to obtain energies in the non-linear
# pb equation.
#
#
/bin/rm dapout.log1 dapin.log1
set molc="dap.pdb dna12.pdb dapdna12.pdb"
foreach cmolc ($molc)
	echo doing run with molc $cmolc
	set salt="0.0 0.01 0.1"
	foreach csalt ($salt)
		echo doing run with salt conc. of $csalt M
		#
		echo outer run
		#
		cp parm parm2
		echo "boundary=zero" >> parm2
      	echo "salt=$csalt" >> parm2
		echo "pdb_in=$cmolc" >> parm2
		echo "phi_out=dapout.phi" >> parm2
		echo "dielectric=dapout.eps" >> parm2
		echo "scale=0.3" >> parm2
		qnifft11 parm2 >> dapout.log1
		#
		echo inner run
		#
		cp parm parm2
		echo "boundary=focus" >> parm2
      	echo "salt=$csalt" >> parm2
		echo "pdb_in=$cmolc" >> parm2
		echo "phi_in=dapout.phi" >> parm2
		echo "phi_out=dapin.phi" >> parm2
		echo "dielectric=dapin.eps" >> parm2
		echo "scale=1.0" >> parm2
		qnifft11 parm2 >> dapin.log1
		ionint.com $csalt >> dapin.log1
	end
end
#
# check charge assignment was done ok
#
chkchg < temp.atm > chkchg.out
grepit.com dapin.log1 > dapi.sum1
# clean up
#rm parm2
exit(0)
