#!/bin/csh -f

diff qnifft12.f qnifft12_old.f
diff getpar.f getpar_old.f
diff mkmaps.f mkmaps_old.f
