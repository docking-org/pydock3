# CORE-MATH

Correctly rounded `acos`, `atan2`, `sin` and `cos` from the [CORE-MATH project](https://core-math.gitlabpages.inria.fr/)
(`src/binary64/{acos,atan2,sin,cos}` at commit 708e86e), used by dms (`../dms/cr_math.h`).

Correctly rounded results are the same on every platform, whereas those of the platforms' math
libraries can differ in the last bit, which changes some of dms's surface points.

Change: `acos.c` always uses its portable `roundeven_finite`, as `__builtin_roundeven` can
compile to a call to `roundeven()`, which some C libraries (MinGW's) lack.
