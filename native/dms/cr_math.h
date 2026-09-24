/*
 * (pydock3) Use correctly rounded math functions (from CORE-MATH, in native/core-math),
 * so that dms computes the same surface on every platform: the platforms' math libraries
 * can differ in the last bit, which changes some of the surface points.
 */
#include <math.h>

double	cr_acos(double);
double	cr_atan2(double, double);
double	cr_sin(double);
double	cr_cos(double);

#define	acos	cr_acos
#define	atan2	cr_atan2
#define	sin	cr_sin
#define	cos	cr_cos
