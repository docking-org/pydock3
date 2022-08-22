c                         PROGRAM PDBTOSPH
c
c  Makes a sphere file for DOCK out of the nonhydrogen atoms in a
c  PDB-format file.
c  ECM May 92.
c
      character*80 pdblin, line
      character*60 pdbfil, sphfil
      real x, y, z
      integer i, count, atmnum
      real radius
      parameter (radius=0.70) 
c
      call getarg(1,pdbfil)
      call getarg (2,sphfil)
      open (unit=1, file=pdbfil, status='old')
      open (unit=2,  status='scratch')
      count=0
   50 read (1, 1000, end=500) pdblin
 1000 format (A80)
      if ((pdblin(1:4) .eq. 'ATOM') .or. 
     +    (pdblin(1:6) .eq. 'HETATM')) then
        if ((pdblin(14:14) .eq. 'H' .or. pdblin(14:14) .eq. 'D')
     +  .or.(pdblin(13:13) .eq. 'H' .or. pdblin(13:13) .eq. 'D'))
     +  go to 50
        count=count + 1
        read (pdblin, 1001) atmnum, x, y, z
 1001   format (6x, I5, 19x, 3F8.3)
        write (2, 1002) atmnum, x, y, z, radius, atmnum
 1002   format (I5, 3F10.5, F8.3, I5)
      endif
      go to 50
  500 continue
      open (unit=3, file=sphfil, status='new')
      write (3, 1003) 'cluster     1   number of spheres in cluster',
     +count
 1003 format (A44, I6)
      rewind (2)
      do 600 i=1, count
        read (2,1000) line
        write (3,1000) line
  600 continue
      close (1)
      close (2)
      close (3)
      close (4)
      end
