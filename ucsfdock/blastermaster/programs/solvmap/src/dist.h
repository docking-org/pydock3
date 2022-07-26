c--------------------------------------------------------------------
c         header file for solmap.  BKS 11/2002 MMM 4/2009
c--------------------------------------------------------------------
        real cutoff
c        parameter (cutoff=10.0)
c  cutoff: the cutoff distance in A for the numerical integration
        real rad_o, rad_c, rad_s, rad_p, rad_n, rad_q
c  rad_o,c,s,p: various atom radii, rad_q is for everything undefined

        common /rads/ rad_o,rad_c,rad_s,rad_p,rad_n,rad_q
