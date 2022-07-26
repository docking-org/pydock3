
c     accessible surface point header

      integer numpts
c     number of evenly spaced surface points to generate on the spheres  
      parameter (numpts=1000)
      real sphpts(3, numpts)

      pointer (i_accstart, accstart)
      pointer (i_accstop, accstop)
      pointer (i_accpts, accpts)

      common /acc/ sphpts,i_accstart,i_accstop,i_accpts
