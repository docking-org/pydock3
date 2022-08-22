c------------------------------------------------------------------------
c                          PROGRAM SOLVMAP
c
c       Copyright (C) 2010 Regents of the University of California
c                         All Rights Reserved.
c
c  changes by BKS (2003, who can believe it?)
c   - program now calculates an displaced volume (by solute = protein)
c   - based partial desolvation energy.  Based on calculation of GB
c   - dielectric
c  changes by Michael Mysinger (2008)
c   - Change numerical integration boundary condition and normalization 
c         to match GB theory
c   - Extend dielectric bitmap computation outside the solvmap grid by the
c         cutoff size of the numerical integration to obtain full accuracy 
c          at box edges
c   - Pre-compute grid kernel to regain former speed
c   - Read ligand atom size from INSOLV
c  changes by Michael Mysinger (200902)
c   - Dynamic memory allocation
c   - "optimal" selection of expanded dielectric grid size
c  changes by Michael Mysinger (200905)
c   - Interface to surface routines
c------------------------------------------------------------------------

       program solvmap
       include 'pointer.h'
       include 'dist.h'
       include 'acc.h'

       real pdbcrd(3,1), rad(1), solgrid(1), lksolv(1)
c  pdbcrd: pdb coordinates of the atoms (dynamic)
c  rad: atomic radii (dynamic) 
c  solgrid: grid of percent desolvation, based on displaced volume (dynamic)
c  lksolv: lookup table for core solvmap kernel (dynamic)
       character*80 pdbnam, scoren, boxlin
c  pdbnam: name of pdb file containing site to be scored.
c  scoren: name of output file to write the grid.
c  boxlin: name of pdb file of limiting box (optional); allows one to use
c    the total receptor file without forcing the grid to enclose it
       integer perang
c  perang: number of grid points per angstrom
       real invgrid, born, probe, radmax, maxsas
c  invgrid: inverse of perang
c  born: atomic radius of the ligand atom (all atoms assumed to be the same)
c  probe: radius of the solvent probe
c  radmax: maximum receptor atom radius
c  maxsas: maximum solvent accessible surface radius
       integer n(3), anum, lksize
c  n: number of grid points in x,y,z directions.
c  anum: number of atoms 
c  lksize: integration cutoff size (in grid points)
       integer maxn
c  maxn: maximum dimension of solvation and dielectric grids
       integer dedge(2,3)
c  dedge: extra width at each edge of grid outside of solgrid
       real edgesolv(2,3),edgegrid(2,3)
c  edgesolv: minimum (1,*) and maximum (2,*) coordinates of solvation grid
c  edgegrid: minimum (1,*) and maximum (2,*) coordinates of dielectric grid
       logical limit
c  limit: whether or not the user has specified a limiting box
       integer accnum
c  accnum: number of solvent accessible surface points in the receptor
c--------------------------------------------------------------------------

       integer ierr
       integer i, j 
       real diff
       character*80 filename
       character*3 surftype
       integer memalloc
       real tarray(2), etime, stime, ftime

       stime = etime(tarray)

       open(unit=5,file='INSEV',status='old')
       open(unit=6,file='OUTSEV',status='unknown')

c  read input parameters from INSOLV.
       write(6,'(t5,a,$)') 'name of input pdb file:'
       read(5,'(a80)') pdbnam
       write(6,*) pdbnam
       write(6,'(t5,a,$)') 'name of output contact grid file:'
       read(5,'(a80)') scoren
       write(6,*) scoren
   
       write(6,'(t5,a)') 'atomic radii (A):'
       write(6,'(t10,a,$)') '1. for O N C S P Other [1.4 1.3 1.7 2.2'
       write(6,'(a)') '2.2 1.8]?'
       read(5,*) rad_o, rad_n, rad_c, rad_s, rad_p, rad_q
       if (rad_o .le. 0.0) rad_o = 1.4
       if (rad_n .le. 0.0) rad_n = 1.3
       if (rad_c .le. 0.0) rad_c = 1.7
       if (rad_s .le. 0.0) rad_s = 2.2
       if (rad_p .le. 0.0) rad_p = 2.2
       if (rad_q .le. 0.0) rad_q = 1.8
       write(6,*) rad_o, rad_n, rad_c, rad_s, rad_p, rad_q

       write(6,'(t10,a,$)') '2. for water probe (exclusion) [1.4]?'
       read(5,*) probe
       if (probe .lt. 0.0) probe = 1.4
       write(6,*) probe
       write(6,'(t5,a,$)') 'number of grid points/angstrom?(1-10)'
       read(5,*) perang
       if (perang .le. 0) perang = 2
       if (perang .gt. 10) perang = 2
       write(6,*)perang
       limit=.false.
       boxlin(1:10)='          '
       born = 1.4
       cutoff = 10.0
       surftype = 'SEV'
       read(5,'(a80)',end=40) boxlin
       if (boxlin(1:10).eq.'          ') go to 20
       write(6,*) 'pdb file of limiting box:', boxlin
       limit=.true.
 20    continue
       read(5, '(f7.3)', err=40, end=40) born
       if (born .lt. 0.4 .or. born .gt. 10.0) then
         write(6,*) 'WARNING: out of bounds ligand atomic radius ',
     &                 'reset to default of 1.4 A!'
         born = 1.4
       endif
       read(5, '(f7.3)', err=40, end=40) cutoff
       if (cutoff .lt. 4.0 .or. cutoff .gt. 40.0) then
         write(6,*) 'WARNING: out of bounds integration cutoff ',
     &                 'reset to default of 10.0 A!'
         cutoff = 10.0
       endif
       read(5, '(a3)', err=40, end=40) surftype
       if (surftype .ne. 'SEV' .and. surftype .ne. 'MS' .and.
     &     surftype .ne. 'SAS' .and. surftype .ne. 'VDW' ) then
         write(6,*) 'WARNING: unrecognized surface type ',
     &                 'reset to default of SEV!'
         surftype = 'SEV'
       endif
 40    continue
       close(5)

       write(6,*) 'surface type: ', surftype
       write(6,*) 'assumed ligand atomic radius (A):', born

       radmax = max(rad_o,rad_n,rad_c,rad_s,rad_p,rad_q)
       write(6,*) 'maximum receptor atomic radius (A):', radmax
       maxsas = radmax + probe

       open(14,file=pdbnam,status='old')
       open(16,file=scoren,status='unknown')
 
       invgrid = 1.0/real(perang)
       write(6,*) 'grid resolution (A): ', invgrid

       lksize = cutoff*perang
       write(6,*) 'integration cutoff (A):', cutoff

c      store receptor edges in edgegrid temporarily
       call readpdb(14,edgegrid,anum)
       if (limit) then
         ierr=0
         call readbox(15,boxlin,edgesolv,ierr)
         if (ierr.eq.0) then
c          determine dielectric edge boundaries
           do i = 1,3
             diff = perang*(edgesolv(1,i) - (edgegrid(1,i) - maxsas))
             dedge(1, i) = nint(diff)
             if (dedge(1, i) .lt. 0) dedge(1, i) = 0
             if (dedge(1, i) .gt. lksize) dedge(1, i) = lksize
           enddo
           do i = 1,3
             diff = perang*((edgegrid(2,i) + maxsas) - edgesolv(2,i))
             dedge(2, i) = nint(diff)
             if (dedge(2, i) .lt. 0) dedge(2, i) = 0
             if (dedge(2, i) .gt. lksize) dedge(2, i) = lksize
           enddo
         else
           limit=.false.
           write(6,'("WARNING -> ",A)') 'error opening limiting box'
           write(6,'("WARNING -> ",A)') 
     &		'   attempting to use entire input pdb file instead'
         endif
       endif
       if (.not. limit) then
         do i = 1,3
           do j = 1,2
               dedge(j,i) = nint(maxsas*perang)
           enddo
           edgesolv(1,i) = edgegrid(1,i)
           edgesolv(2,i) = edgegrid(2,i)
           edgesolv(3,i) = edgegrid(3,i)
         enddo
       endif

       write(6,*) 'dielectric grid expansion size (points/edge)'
       do i=1,2
         write(6,*) (dedge(i,j), j=1,3)
       enddo
c      store dielectric grid boundaries
       do i = 1,3
         edgegrid(1, i) = edgesolv(1, i) - dedge(1, i)/real(perang)
         edgegrid(2, i) = edgesolv(2, i) + dedge(2, i)/real(perang)
       enddo
       write(6,*)'dielectric grid minimum and maximum x y z (A):'
       do i=1,2
         write(6,*) (edgegrid(i,j), j=1,3)
       enddo

       filename = 'solv_sev.box'
       call mkbox(15,edgesolv,filename)

       filename = 'dielec_sev.box'
c  make pdb-format file of box outlining the solvation grid
       call mkbox(15,edgegrid,filename)

c  calculate number of grid points per dimension
       call gridinfo(edgesolv,perang,n)
       maxn = max(n(1), n(2), n(3)) + 2*lksize

c  generate accessible surface points
       if ( surftype .eq. 'SEV' .or. surftype .eq. 'MS' ) then
         call accsurf(anum,edgegrid,radmax,probe,accnum)
c        debugging output, all accessible dots on the SAS
c         filename = 'acc.pdb'
c         call writeacc(15, anum, accnum, filename)
       endif

       if ( surftype .eq. 'VDW' ) then
c        vdw surface
         write(6,*) 'WARNING: ignoring INSOLV probe radius'
         write(6,*) 'Using only the receptor van der Waals surface!!!'
         call dielec(anum,n,maxn,perang,0.0,edgesolv,dedge)
c         filename = 'vdw.bmp'
c         call writegrid(15,n,maxn,perang,edgegrid,dedge,filename)
       else
c        compute dielectric boundary of receptor
         call dielec(anum,n,maxn,perang,probe,edgesolv,dedge)
c         filename = 'sas.bmp'
c         call writegrid(15,n,maxn,perang,edgegrid,dedge,filename)

         if ( surftype .ne. 'SAS' ) then
c          compute connolly molecular surface of receptor
           call recblot(accnum,n,maxn,perang,probe,edgesolv,dedge)
c           filename = 'ms.bmp'
c           call writegrid(15,n,maxn,perang,edgegrid,dedge,filename)
         endif

       endif


c  precaculate lookup table
       call solvlookup(perang,born,lksize)

       if ( surftype .eq. 'SEV' ) then
c        new solvmap routine
         call walklig(anum,accnum,maxn,lksize,perang,radmax,
     &                  born,probe,n,dedge,edgesolv,edgegrid)
       else
c        old solvmap routine - used for alternate solvmaps
         call score(n,maxn,dedge,lksize)
       endif

c  deallocate old dynamic arrays
       i_pdbcrd = memalloc(i_pdbcrd, 0)
       i_rad = memalloc(i_rad, 0)
       i_grid = memalloc(i_grid, 0)
       i_lksolv = memalloc(i_lksolv, 0)
       i_accstart = memalloc(i_accstart, 0)
       i_accstop = memalloc(i_accstop, 0)
       i_accpts = memalloc(i_accpts, 0)

c  dump grid in format of distmap, with floats for real space grid center
       call bdump(n,maxn,perang,edgesolv)

       filename = 'solv_sev.plt'
       call writeplt(15,n,maxn,perang,edgesolv,filename)

c  deallocate solvation grid
       i_solgrid = memalloc(i_solgrid, 0)

       write(6,*)'successful completion'

       ftime = etime(tarray)
       write(6,'(a19,f14.4)') 'elapsed time (sec): ', ftime-stime

       close(6)
       close(16)
       close(14)
       stop
       end

c------------------------------------------------------------------------
c   subroutine readpdb
c  -rewind and read coordinates (pdbcrd) and radii (radi) into arrays
c  -find the highest and lowest x, y, and z values present
c------------------------------------------------------------------------
       subroutine readpdb(lun, edgegrid, anum)

       include 'pointer.h'
       include 'dist.h'

       integer lun, anum
       real crd(3)
       real edgegrid(2,3)
       integer i, j
       character*80 string
       character*4 atypr
       integer cnum
       real pdbcrd(3,*),rad(*)
       integer memalloc

c      allocate initial arrays
       cnum = 5000
       i_pdbcrd = memalloc(i_pdbcrd, 4*3*cnum)
       i_rad = memalloc(i_rad, 4*cnum)

c      initialize mins to high value, maxs to low value.
       do 10 i = 1,3
         edgegrid(1,i) =  1000.
         edgegrid(2,i) = -1000.
 10    continue

       i = 0
 70    read(lun, '(a80)', end=110) string
       if (string(1:4).ne.'ATOM' .and. string(1:4).ne.'HETA') go to 70
       if (string(13:13).eq.'H' .or. string(13:13).eq.'D' .or.
     &       string (14:14).eq.'H' .or. string(14:14).eq.'D') go to 70
       i = i + 1
 5     format(12x,a4,14x,3f8.3)
       read (string, 5, end=110) atypr, (pdbcrd(j,i), j=1,3)
c       write(6,*) 'coordinates: ', (pdbcrd(j,i), j=1,3)
       do 90 j=1,3 
         if (pdbcrd(j,i) .lt. edgegrid(1,j)) edgegrid(1,j)=pdbcrd(j,i)
         if (pdbcrd(j,i) .gt. edgegrid(2,j)) edgegrid(2,j)=pdbcrd(j,i)
 90    continue
       if (atypr(2:2).eq.'O') then
            rad(i) = rad_o
         elseif (atypr(2:2).eq.'N') then
            rad(i) = rad_n
         elseif(atypr(2:2).eq.'C') then
            rad(i) = rad_c
         elseif(atypr(2:2).eq.'S') then
            rad(i) = rad_s
         elseif(atypr(2:2).eq.'P') then
            rad(i) = rad_p
         else 
            rad(i) = rad_q
       endif
c      grow dynamic atom arrays as needed
       if (i.ge.cnum) then
          cnum = cnum + 5000
          i_pdbcrd = memalloc(i_pdbcrd, 4*3*cnum)
          i_rad = memalloc(i_rad, 4*cnum)
       endif
       go to 70
110    continue
       anum = i
c      shrink dynamic arrays to only use the space needed
       i_pdbcrd = memalloc(i_pdbcrd, 4*3*anum)
       i_rad = memalloc(i_rad, 4*anum)
       write(6,*) 'receptor minimum and maximum x y z (A):'
       do 130 i=1,2
         write(6,*) (edgegrid(i,j), j=1,3)
130    continue
       return
       end

c------------------------------------------------------------------------
c   subroutine gridinfo
c  --calculates number of grid points in each dimension.
c------------------------------------------------------------------------
        subroutine gridinfo(edgesolv,perang,n)

        real edgesolv(2,3)
        integer perang
        integer n(3)

        n(1) = int((edgesolv(2,1) - edgesolv(1,1))*perang)
        n(2) = int((edgesolv(2,2) - edgesolv(1,2))*perang)
        n(3) = int((edgesolv(2,3) - edgesolv(1,3))*perang)
        write(6,*) 'solvation grid x, y, z dimensions:', n(1),n(2),n(3)
        return
        end

c------------------------------------------------------------------------
c   subroutine score
c   - computes percentage desolvation at each solgrid point
c------------------------------------------------------------------------
        subroutine score(n,maxn,dedge,lksize)

        include 'pointer.h'

        integer n(3), dedge(2,3)
        integer maxn,lksize
        integer grid(0:maxn,0:maxn,0:maxn)
c  grid: bump grid (dynamic)
        real solgrid(0:maxn,0:maxn,0:maxn)
c  solgrid: grid of percent desolvation, based on displaced volume (dynamic)
        real lksolv(-lksize:lksize,-lksize:lksize,-lksize:lksize)
c  lksolv: lookup table for core solvmap kernel (dynamic)
        integer i,j,k,l,ii,jj,kk
        integer minb(3),maxb(3),upper
        integer memalloc

c  leave all zero array initialization to compiler and calloc
        i_solgrid = memalloc(i_solgrid, 4*(maxn+1)**3)

 445     format(i4,a,i4,a)
         upper = n(3) + dedge(1,3) + dedge(2,3) + 1
         do 540 k = -dedge(2,3), n(3)+dedge(1,3)
             write(6,445) k+dedge(2,3)+1,' /',upper,' the way there!'
             call flush (6)
             minb(3)=k-lksize
             maxb(3)=k+lksize
c  check solvation bounds, which remain tight at 0 and n(*)
             if (minb(3) .lt. 0 ) minb(3) = 0
             if (maxb(3) .gt. n(3) ) maxb(3) = n(3)
          do 550 j = -dedge(2,2), n(2)+dedge(1,2)
             minb(2)=j-lksize
             maxb(2)=j+lksize
             if (minb(2) .lt. 0 ) minb(2) = 0
             if (maxb(2) .gt. n(2) ) maxb(2) = n(2)
           do 560 i = -dedge(2,1), n(1)+dedge(1,1)
c  dielectric grid is offset +dedge(2,*) so it starts at 0,0,0
            if (grid(i+dedge(2,1),j+dedge(2,2),k+dedge(2,3)).gt.0) then
             minb(1)=i-lksize
             maxb(1)=i+lksize
             if (minb(1) .lt. 0 ) minb(1) = 0
             if (maxb(1) .gt. n(1) ) maxb(1) = n(1)
             do kk = minb(3), maxb(3)
               do jj = minb(2), maxb(2)
                 do ii = minb(1), maxb(1)
                   solgrid(ii,jj,kk) = solgrid(ii,jj,kk) + 
     &                  lksolv(i-ii,j-jj,k-kk)
                 enddo
               enddo
             enddo
            endif
  560      continue
  550     continue
  540    continue

        return
        end

c------------------------------------------------------------------------
c  - precalculates lookup tables over sub-grid for fast calculation
c------------------------------------------------------------------------
        subroutine solvlookup(perang,born,lksize)

        include 'pointer.h'

        integer perang, lksize
        real invgrid, born
        real lksolv(-lksize:lksize,-lksize:lksize,-lksize:lksize)
c  lksolv: lookup table for core solvmap kernel (dynamic)
        integer i, j, k
        real x2,y2,z2
        real invgrid2
        real dV, dist, dlim
        integer memalloc

        i_lksolv = memalloc(i_lksolv, 4*(2*lksize+1)**3)
        
        invgrid = 1.0/real(perang)
        invgrid2 = invgrid**2
        dV = invgrid**3/(4*3.14159265)
        write(6,*) 'dV is', dV
c     normalize solvation grid using born radius of a typical atom
        dV = born*dV
        dlim = born**2

        write (6, *) 'Precalculating solvmap lookup table'
        do k = -lksize, lksize
          z2 = invgrid2*real(k)**2
          do j = -lksize, lksize
            y2 = invgrid2*real(j)**2
            do i = -lksize, lksize
              x2 = invgrid2*real(i)**2
              dist = x2 + y2 + z2
c     only count grid points outside of the vdw of the reference atom
              if (dist .gt. dlim) lksolv(i,j,k) = dV/dist**2
            enddo
          enddo
        enddo

        return
        end

c------------------------------------------------------------------------
c   subroutine bdump
c  -writes out grid in distmap format for compatibility with output by 
c   distmap.f.
c------------------------------------------------------------------------
        subroutine bdump(n,maxn,perang,edgesolv)

        include 'pointer.h'

        integer maxn
        real solgrid(0:maxn,0:maxn,0:maxn)
c  solgrid: grid of percent desolvation, based on displaced volume (dynamic)
        integer nlist,i,j,k,perang,neg, kk
        integer n(3)
        real edgesolv(2,3)
        logical wlist
        data nlist / 1 /

        write(16,740)n(1),n(2),n(3),perang,edgesolv(2,1)*perang,
     &         edgesolv(2,2)*perang,edgesolv(2,3)*perang
 740    format(4(i4,2x),3(f8.3,1x))
         do 550 i = 0,n(1)
           do 350 j = 0,n(2)
             do 150 k = 0,n(3),13
               if (k+12.le.n(3)) then
                  write(16,761) (solgrid(i,j,kk),kk=k,k+12)
               else
                  if (k.le.n(3)-1) then
                    write(16,762) (solgrid(i,j,kk),kk=k,n(3)-1)
                  endif
                  write(16,'(f6.3)') solgrid(i,j,n(3))
               endif
 761              format(13f6.3)    
 762              format(f6.3,$)    
150           continue
350       continue     
550     continue      
        nlist=1
        return
        end

c------------------------------------------------------------------------
      subroutine mkbox(unitno, edgesolv, filename)
c  --makes a pdb format box to show the size and location of the grid.
c------------------------------------------------------------------------
      real boxdim(3), cent(3), edgesolv(2,3)
      integer unitno, i
      character*80 filename
c
      do 10 i=1,3
        cent(i) = (edgesolv(2,i) + edgesolv(1,i))/2.0
        boxdim(i) = edgesolv(2,i) - edgesolv(1,i)
   10 continue
c
      open (unit=unitno, file=filename, status='unknown')
    7 format ('HEADER    CORNERS OF BOX ', 6F8.3)
      write (unitno, 7) (edgesolv(1,i), i=1,3), (edgesolv(2,i), i=1,3)
      write (unitno, 1) 'REMARK    CENTER (X Y Z) ', (cent(i), i=1,3)
    1 format (A25, 3F8.3)
      write (unitno, 2) 'REMARK    DIMENSIONS (X Y Z) ',
     &(boxdim(i), i=1,3)
    2 format (A29, 3F8.3)
      write (unitno, 3) 'ATOM', 1, 'DUA', 'BOX', 1,
     &edgesolv(1,1), edgesolv(1,2), edgesolv(1,3)
      write (unitno, 3) 'ATOM', 2, 'DUB', 'BOX', 1,
     &edgesolv(2,1), edgesolv(1,2), edgesolv(1,3)
      write (unitno, 3) 'ATOM', 3, 'DUC', 'BOX', 1,
     &edgesolv(2,1), edgesolv(1,2), edgesolv(2,3)
      write (unitno, 3) 'ATOM', 4, 'DUD', 'BOX', 1,
     &edgesolv(1,1), edgesolv(1,2), edgesolv(2,3)
      write (unitno, 3) 'ATOM', 5, 'DUE', 'BOX', 1,
     &edgesolv(1,1), edgesolv(2,2), edgesolv(1,3)
      write (unitno, 3) 'ATOM', 6, 'DUF', 'BOX', 1,
     &edgesolv(2,1), edgesolv(2,2), edgesolv(1,3)
      write (unitno, 3) 'ATOM', 7, 'DUG', 'BOX', 1,
     &edgesolv(2,1), edgesolv(2,2), edgesolv(2,3)
      write (unitno, 3) 'ATOM', 8, 'DUH', 'BOX', 1,
     &edgesolv(1,1), edgesolv(2,2), edgesolv(2,3)
    3 format (A4, 6x, I1, 2x, A3, 1x, A3, 5x, I1, 4x, 3F8.3)
      write (unitno, 4) 'CONECT', 1, 2, 4, 5
      write (unitno, 4) 'CONECT', 2, 1, 3, 6
      write (unitno, 4) 'CONECT', 3, 2, 4, 7
      write (unitno, 4) 'CONECT', 4, 1, 3, 8
      write (unitno, 4) 'CONECT', 5, 1, 6, 8
      write (unitno, 4) 'CONECT', 6, 2, 5, 7
      write (unitno, 4) 'CONECT', 7, 3, 6, 8
      write (unitno, 4) 'CONECT', 8, 4, 5, 7
    4 format (A6, 4I5)
      close (unitno)
      return
      end

c------------------------------------------------------------------------
      subroutine readbox(unitno,boxlin,edgesolv,ierr)
c  --reads in a pdb file of a limiting box, as written using the program
c    showbox or the program chemgrid.
c------------------------------------------------------------------------
      real edgesolv(2,3)
      integer unitno, i, j, ierr
      character*80 boxlin,line
c
      ierr=0
      open (unit=unitno, file=boxlin, status='old',err=80)
   10 read (unitno, 1000) line
 1000 format (a80)
      if (line(1:11).eq.'ATOM      1') then
c       read minimum coordinate
        read (line, '(30x,3f8.3)') (edgesolv(1,i), i=1,3)
        go to 10
      else if (line(1:11).eq.'ATOM      7') then
c       read maximum coordinate
        read (line, '(30x,3f8.3)') (edgesolv(2,i), i=1,3)
        go to 90
      else
        go to 10
      endif
   80 continue
      ierr=1
   90 continue
      write(6,*)'box minimum and maximum x y z (A):'
      do 200 i=1,2
        write(6,*) (edgesolv(i,j), j=1,3)
  200 continue
      return
      end
