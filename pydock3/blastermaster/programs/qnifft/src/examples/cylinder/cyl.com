#!/bin/csh -f
# run test case for cylinder 1e/A, in mixed salt
#
setenv RUN _22
alias qnifft qnifft22_145c
#
/bin/rm cyl$RUN.sum
#goto grepit
/bin/rm cyl$RUN.log
#
# 0.5,0.5,0.5,0.5 Salt
# outer
cp parm parm1
echo 'nonlinear=t' >> parm1
echo 'zper=t' >> parm1
echo 'site_pot=t' >> parm1
echo 'salt=0.' >> parm1
echo 'ani_mon=0.5' >> parm1
echo 'cat_mon=0.5' >> parm1
echo 'ani_div=0.5' >> parm1
echo 'cat_div=0.5' >> parm1
echo 'boundary=zero' >> parm1
echo 'fill=130' >> parm1
echo 'phi_out=cylout.phi' >> parm1
echo 'dielectric=cylout.eps' >> parm1
echo 'site_in=cyl_pot.pdb' >> parm1
echo 'site_out=cyl_I500'$RUN'o.fld' >> parm1
qnifft parm1 >> cyl$RUN.log
# inner
echo 'boundary=focus' >> parm1
echo 'fill=260' >> parm1
echo 'zper=f' >> parm1
echo 'phi_in=cylout.phi' >> parm1
echo 'phi_out=cylin.phi' >> parm1
echo 'dielectric=cylin.eps' >> parm1
echo 'site_out=cyl_I500'$RUN'.fld' >> parm1
qnifft parm1 >> cyl$RUN.log
#
grepit:
ionintm.com > ionintm.out
#
grep -a rxn cyl$RUN.log >> cyl$RUN.sum
grep -a Interaction cyl$RUN.log >> cyl$RUN.sum
grep -a dPI cyl$RUN.log >> cyl$RUN.sum
grep -a 'E.D' cyl$RUN.log >> cyl$RUN.sum
grep -a 'rho.phi' cyl$RUN.log >> cyl$RUN.sum
grep -a dPI ionintm.out >> cyl$RUN.sum
grep -a 'E.D' ionintm.out >> cyl$RUN.sum
grep -a 'rho.phi' ionintm.out >> cyl$RUN.sum
