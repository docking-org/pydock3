#!/bin/csh -f
# run test case for sphere in 0 salt, 1M linaer, 1M non-linear
#
setenv RUN _14node12
alias qnifft qnifft22_lnx65
#alias qnifft qnifft22_145c
#
/bin/rm sph$RUN.sum
#goto grepit
/bin/rm sph$RUN.log
#
# 0 Salt
# outer
cp parm parm1
echo 'site_pot=t' >> parm1
echo 'salt=0.' >> parm1
echo 'ani_mon=0.' >> parm1
echo 'cat_mon=0.' >> parm1
echo 'boundary=zero' >> parm1
echo 'fill=20' >> parm1
echo 'phi_out=sphout.phi' >> parm1
echo 'site_in=sph_pot.pdb' >> parm1
echo 'site_out=sph_I0'$RUN'o.fld' >> parm1
qnifft parm1 >> sph$RUN.log
# inner
echo 'boundary=focus' >> parm1
echo 'fill=70' >> parm1
echo 'phi_in=sphout.phi' >> parm1
echo 'phi_out=sphin.phi' >> parm1
echo 'site_out=sph_I0'$RUN'.fld' >> parm1
qnifft parm1 >> sph$RUN.log
#
# 1M Salt, linear
# outer
cp parm parm1
echo 'site_pot=t' >> parm1
echo 'salt=1.' >> parm1
echo 'nonlinear=f' >> parm1
echo 'ani_mon=1.' >> parm1
echo 'cat_mon=1.' >> parm1
echo 'boundary=zero' >> parm1
echo 'fill=20' >> parm1
echo 'site_in=sph_pot.pdb' >> parm1
echo 'site_out=sph_I1lin'$RUN'o.fld' >> parm1
echo 'phi_out=sphout.phi' >> parm1
echo 'diel=sphout.eps' >> parm1
qnifft parm1 >> sph$RUN.log
# inner
echo 'boundary=focus' >> parm1
echo 'fill=70' >> parm1
echo 'phi_in=sphout.phi' >> parm1
echo 'phi_out=sphin.phi' >> parm1
echo 'diel=sphin.eps' >> parm1
echo 'site_out=sph_I1lin'$RUN'.fld' >> parm1
qnifft parm1 >> sph$RUN.log
#
# 1M Salt, non-linear
# outer
cp parm parm1
echo 'site_pot=t' >> parm1
echo 'salt=1.' >> parm1
echo 'nonlinear=t' >> parm1
echo 'ani_mon=1.' >> parm1
echo 'cat_mon=1.' >> parm1
echo 'boundary=zero' >> parm1
echo 'fill=20' >> parm1
echo 'site_in=sph_pot.pdb' >> parm1
echo 'site_out=sph_I1non'$RUN'o.fld' >> parm1
echo 'phi_out=sphout.phi' >> parm1
echo 'diel=sphout.eps' >> parm1
qnifft parm1 >> sph$RUN.log
# inner
echo 'boundary=focus' >> parm1
echo 'fill=70' >> parm1
echo 'phi_in=sphout.phi' >> parm1
echo 'phi_out=sphin.phi' >> parm1
echo 'diel=sphin.eps' >> parm1
echo 'site_out=sph_I1non'$RUN'.fld' >> parm1
qnifft parm1 >> sph$RUN.log
#
grepit:
ionintm.com > ionintm.out
#
grep -a 'grid=' parm1 >> sph$RUN.sum
grep -a rxn sph$RUN.log >> sph$RUN.sum
grep -a Interaction sph$RUN.log >> sph$RUN.sum
grep -a dPI ionintm.out >> sph$RUN.sum
grep -a 'E.D' ionintm.out >> sph$RUN.sum
grep -a 'rho.phi' ionintm.out >> sph$RUN.sum
grep -a cpu sph$RUN.log >> sph$RUN.sum
