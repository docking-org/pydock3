c------------------------------------------------------------------------
c
c       Copyright (C) 2010 Regents of the University of California
c                         All Rights Reserved.
c
c                     Solvmap Surface Subroutines
c
c   This file contains solvmap extensions to compute Connolly surfaces
c     and the atomic solvent-excluded volume between protein and ligand.
c     Utility routines are provided to write internal grids to file formats
c     recognized by UCSF Chimera.
c
c
c   Michael Mysinger 200905 Created
c   Michael Mysinger 201005 Add input control on which surface is used
c------------------------------------------------------------------------

c------------------------------------------------------------------------
c   subroutine writeplt
c   - write solvmap as a gopenmol .plt file (viewable in chimera)
c------------------------------------------------------------------------
       subroutine writeplt(lun,n,maxn,perang,edgesolv,filename)
         include 'pointer.h'

         integer lun, n(3), maxn, perang
         real edgesolv(2,3)
         character*80 filename
         real solgrid(0:maxn,0:maxn,0:maxn)
         real scale, oldmid(3)
         integer i, j, k
         integer rank, type

         rank = 3
         type = 10
         open(unit=lun, file=filename, access='sequential',
     &        form='binary', status='unknown')
         write(lun) rank
         write(lun) type
         write(lun) n(3)+1,n(2)+1,n(1)+1
         do i = 3,1,-1
           write(lun) edgesolv(1,i), edgesolv(2,i)
         enddo
         do k = n(3),0,-1
           do j = n(2),0,-1
             do i = n(1),0,-1
               write(lun) solgrid(i,j,k)
             enddo
           enddo
         enddo

         close(lun)
         return
       end

c------------------------------------------------------------------------
c   subroutine writegrid
c   - write dielectric grid as a dock 4,5,6 style bump map (.bmp)
c------------------------------------------------------------------------
       subroutine writegrid(lun,n,maxn,perang,edgegrid,dedge,filename)
         include 'pointer.h'

         integer lun, n(3), maxn, perang
         real edgegrid(2,3)
         integer dedge(2,3)
         character*80 filename
         integer grid(0:maxn,0:maxn,0:maxn)
         integer size, dims(3)
         real invgrid, origin(3)
         integer i,j,k,out

         invgrid = 1.0/real(perang)
         do i = 1,3
           dims(i) = n(i) + dedge(1,i) + dedge(2,i) + 1
c          solvmap grid is hung backwards from edgegrid(2,*)
           origin(i) = edgegrid(2,i) - (dims(i)-1)*invgrid
         enddo
         size = dims(1)*dims(2)*dims(3)

         open(unit=lun, file=filename, access='sequential',
     &        form='binary', status='unknown')

         write(lun) size
         write(lun) invgrid
         write(lun) origin
         write(lun) dims(1), dims(2), dims(3)

c        iterate biggest to smallest index to invert grid back to normal
c        loops nested k,j,i because file is stored in fortran order
         do k = n(3)+dedge(1,3)+dedge(2,3),0,-1
           do j = n(2)+dedge(1,2)+dedge(2,2),0,-1
             do i = n(1)+dedge(1,1)+dedge(2,1),0,-1
               if (grid(i,j,k) .eq. 0) then
                 out = 85
               else if (grid(i,j,k) .gt. 7500) then
                 out = 88
               else if (grid(i,j,k) .gt. 5000) then
                 out = 87
               else if (grid(i,j,k) .gt. 0) then
                 out = 86
               else
                 out = 84 + grid(i,j,k)/5
               endif
               write(lun) char(out)
             enddo
           enddo
         enddo

         close(lun)
         return
       end

c------------------------------------------------------------------------
c   subroutine writeacc
c   - write accessible points array as a pdb file 
c------------------------------------------------------------------------
       subroutine writeacc(lun, anum, accnum, filename)
         include 'acc.h'

         integer lun, anum, accnum
         character*80 filename
         integer accstart(anum), accstop(anum)
         real accpts(3,accnum)
         integer i,j,k

         open(unit=lun, file=filename, status='unknown')

         do i = 1,anum
           do j = accstart(i),accstop(i)
             write(lun,'(a4,i7,a9,i6,4x,3f8.3)') 'ATOM', j, '  C   UNK',
     &             i, (accpts(k,j), k=1,3)
           enddo
         enddo

         close(lun)
         return
       end

c------------------------------------------------------------------------
c   subroutine sphere
c   - generate evenly spaced points on the unit sphere
c   - uses the golden section spiral algorithm
c------------------------------------------------------------------------
       subroutine sphere(numpts, sphpts)
         integer numpts
         real sphpts(3, numpts)
         integer i
         real long, dlong, z, dz, r, pi
         parameter(pi = 3.14159265359)

         dz = 2.0/numpts
         z = 1 - dz/2.0
         dlong = pi*(3.0-sqrt(5.0))
         long = 0.0
         do i = 1,numpts
            r = sqrt(1 - z*z)
            sphpts(1, i) = cos(long)*r
            sphpts(2, i) = sin(long)*r
            sphpts(3, i) = z
            z = z - dz
            long = long + dlong
         enddo
         return
       end

c------------------------------------------------------------------------
c   subroutine accsurf
c   - find the accessible surface points that surround the receptor
c------------------------------------------------------------------------
       subroutine accsurf(anum,edgegrid,radmax,probe,accnum)
         include 'pointer.h'
         include 'cubic.h'
         include 'acc.h'

         integer anum
         real edgegrid(2,3)
         real radmax, probe
         integer accnum

         real apt(3)
         real pdbcrd(3,anum), rad(anum)
         real cbln
         integer i, j, k, n
         integer liml, limu, icb(3), ind
         integer curp, startp, allocnum
         real dx1, dx2, dx3, d2

         integer cblower(1), cbupper(1), cbatom(1)

         integer accstart(1), accstop(1)
         real accpts(3,1)

         integer memalloc
         pointer (i_radp2, radp2)
         real radp2(1)

c        generate evenly spaced unit sphere points
         call sphere(numpts, sphpts)

c        distance divide receptor atoms for fast intersection computations
c        cube length is maximum distance at which two SAS could intersect
         cbln = 2.0*(radmax + probe)
         call cubedata(1.0, edgegrid, cbln)
         i_cblower = memalloc(i_cblower,4*(cb(1)+1)*(cb(2)+1)*(cb(3)+1))
         i_cbupper = memalloc(i_cbupper,4*(cb(1)+1)*(cb(2)+1)*(cb(3)+1))
         i_cbatom = memalloc(i_cbatom,4*27*anum)
         call cube(anum,pdbcrd,rad)

c        allocate arrays for the accessible SAS point data
c        accessible point storage scheme
c         - 2 arrays over atoms contain lower,upper bound 
c             pointers into the accpts array (accstart, accstop)
c         - 1 dynamic array containing accessible points (accpts)
         i_accstart = memalloc(i_accstart,4*anum)
         i_accstop = memalloc(i_accstop,4*anum)
         allocnum = 5000
         i_accpts = memalloc(i_accpts,4*3*allocnum)

         i_radp2 = memalloc(i_radp2,4*anum)
c        initialize start, stop pointers to 1, 0
         do i = 1,anum
           accstart(i) = 1
           accstop(i) = 0
           radp2(i) = (rad(i) + probe)**2
         enddo

         curp = 1
c        generate the accessible sas points
         do 100 i = 1,anum
           startp = curp
           do 110 j = 1,numpts
c            generate sas points
             do k = 1,3
               apt(k) = sphpts(k,j)*(rad(i)+probe) + pdbcrd(k,i)
               icb(k) = (apt(k) - cmin(k))*cpa
               if (icb(k) .le. 0 .or. icb(k) .ge. cb(k)) then
                 goto 110
               endif
             enddo
c            find nearest receptor atoms
             ind = 1+icb(1)+(cb(1)+1)*icb(2)+(cb(1)+1)*(cb(2)+1)*icb(3)
             liml = cblower(ind)
             limu = cbupper(ind)
             do n = liml,limu
               k = cbatom(n)
               if (k .ne. i) then
                 dx1 = pdbcrd(1,k) - apt(1)
                 dx2 = pdbcrd(2,k) - apt(2)
                 dx3 = pdbcrd(3,k) - apt(3)
                 d2 = dx1*dx1 + dx2*dx2 + dx3*dx3
c                if inside any other atom, apt is inaccessible
                 if (d2 .lt. radp2(k)) then
                   goto 110
                 endif
               endif
             enddo
c            otherwise, apt is accessible
             do k = 1,3
               accpts(k, curp) = apt(k)
             enddo
             curp = curp + 1
c            re-allocate accpts as needed
             if (curp + numpts .gt. allocnum) then
               allocnum = allocnum + 5000
               i_accpts= memalloc(i_accpts,4*3*allocnum)
             endif
 110       continue
c          set start and stop pointers
           if (curp .gt. startp) then
             accstart(i) = startp
             accstop(i) = curp - 1
           endif
 100     continue      

         i_radp2 = memalloc(i_radp2,0)
         i_cblower = memalloc(i_cblower,0)
         i_cbupper = memalloc(i_cbupper,0)
         i_cbatom = memalloc(i_cbatom,0)
         accnum = curp-1
         i_accpts = memalloc(i_accpts,4*3*accnum)
         write(6,*) 'finished receptor accessible surface'
         return
       end

c------------------------------------------------------------------------
c   subroutine blot
c   - generic blotting subroutine for making points more inaccessible
c   - if grid <= 0 and within blot then add increment to grid
c   - if grid becomes -1 (inaccessible), set to penalty value
c------------------------------------------------------------------------
       subroutine blot(size,pt,increment,penalty,n,maxn,perang,
     &                 edgesolv,dedge)
         include 'pointer.h'

         real size, pt(3)
         integer increment, penalty
         integer n(3), maxn, perang
         real edgesolv(2,3)
         integer dedge(2,3)

         real invgrid, size2, dist2
         real ptpdb(3)
         integer i,j,k,ip,jp,kp
         integer p,q
         integer minb(3),maxb(3)
         real tdist2, tdist3

c        grid: dielectric grid (dynamic)
         integer grid(0:maxn,0:maxn,0:maxn)

         invgrid = 1.0/real(perang)
         size2 = size*size

c  note that the grids are hung starting at the maximum edge (edgesolv(2,*))
c  higher coordinates in real space = low coordinates in grid space
c  mmm - this is ass backwards, and will mess you up in the following code
c  in particular this switches around the minimum and maximum dedge code

c        find min, max possible grid points
c        check dielectric bounds, from -dedge(2,*) to n(*)+dedge(1,*)
         minb(1)=nint((edgesolv(2,1)-(pt(1)+size))*perang)
         maxb(1)=nint((edgesolv(2,1)-(pt(1)-size))*perang)
         if (minb(1) .lt. -dedge(2,1)) minb(1)=-dedge(2,1)
         if (maxb(1) .gt. n(1)+dedge(1,1)) maxb(1)=n(1)+dedge(1,1)
         minb(2)=nint((edgesolv(2,2)-(pt(2)+size))*perang)
         maxb(2)=nint((edgesolv(2,2)-(pt(2)-size))*perang)
         if (minb(2) .lt. -dedge(2,2)) minb(2)=-dedge(2,2)
         if (maxb(2) .gt. n(2)+dedge(1,2)) maxb(2)=n(2)+dedge(1,2)
         minb(3)=nint((edgesolv(2,3)-(pt(3)+size))*perang)
         maxb(3)=nint((edgesolv(2,3)-(pt(3)-size))*perang)
         if (minb(3) .lt. -dedge(2,3)) minb(3)=-dedge(2,3)
         if (maxb(3) .gt. n(3)+dedge(1,3)) maxb(3)=n(3)+dedge(1,3)

c        loop over all grid points in the cube between minb and maxb.
         do k = minb(3), maxb(3) 
c          dielectric grid is offset +dedge(2,*) so it starts at 0,0,0
           kp = k + dedge(2,3)
           ptpdb(3) = edgesolv(2,3) - k*invgrid
           tdist3 = (pt(3) - ptpdb(3))**2
           do j = minb(2), maxb(2)
             jp = j + dedge(2,2)
             ptpdb(2) = edgesolv(2,2) - j*invgrid
             tdist2 = (pt(2) - ptpdb(2))**2
             do i = minb(1), maxb(1)
               ip = i + dedge(2,1)
c              increment grid if negative and within blot size of pt
               if (grid(ip,jp,kp) .le. 0) then
                 ptpdb(1) = edgesolv(2,1) - i*invgrid
                 dist2 = (pt(1) - ptpdb(1))**2 + tdist2 + tdist3
                 if (dist2 .le. size2) then
                   grid(ip,jp,kp) = grid(ip,jp,kp) + increment
                   if (grid(ip,jp,kp) .eq. -1) then
                     grid(ip,jp,kp) = penalty
                   endif
                 endif
               endif
             enddo 
           enddo
         enddo

         return
       end

c------------------------------------------------------------------------
c   subroutine accblot
c   - generic blotting subroutine for making points more accessible
c   - 0 < grid <= decnum are decremented by decnum and thresh if in blot
c   - grid <= 0 are decremented by thresh only if in blot
c   - if a grid pt become -1 it is adjusted to -2, this saves -1 to be a 
c       sentinel during blot operations
c------------------------------------------------------------------------
       subroutine accblot(size,pt,decrement,threshold,n,maxn,perang,
     &                    edgesolv,dedge)
         include 'pointer.h'

         real size, pt(3)
         integer decrement, threshold
         integer n(3), maxn, perang
         real edgesolv(2,3)
         integer dedge(2,3)

c        grid: dielectric grid (dynamic)
         integer grid(0:maxn,0:maxn,0:maxn)
         real invgrid, size2, dist2
         real ptpdb(3)
         integer i,j,k,ip,jp,kp
         integer p,q
         integer minb(3),maxb(3)
         real tdist2, tdist3

         invgrid = 1.0/real(perang)
         size2 = size*size

c  note that the grids are hung starting at the maximum edge (edgesolv(2,*))
c  higher coordinates in real space = low coordinates in grid space

c        find min, max possible grid points
c        check dielectric bounds, from -dedge(2,*) to n(*)+dedge(1,*)
         minb(1)=nint((edgesolv(2,1)-(pt(1)+size))*perang)
         maxb(1)=nint((edgesolv(2,1)-(pt(1)-size))*perang)
         if (minb(1) .lt. -dedge(2,1)) minb(1)=-dedge(2,1)
         if (maxb(1) .gt. n(1)+dedge(1,1)) maxb(1)=n(1)+dedge(1,1)
         minb(2)=nint((edgesolv(2,2)-(pt(2)+size))*perang)
         maxb(2)=nint((edgesolv(2,2)-(pt(2)-size))*perang)
         if (minb(2) .lt. -dedge(2,2)) minb(2)=-dedge(2,2)
         if (maxb(2) .gt. n(2)+dedge(1,2)) maxb(2)=n(2)+dedge(1,2)
         minb(3)=nint((edgesolv(2,3)-(pt(3)+size))*perang)
         maxb(3)=nint((edgesolv(2,3)-(pt(3)-size))*perang)
         if (minb(3) .lt. -dedge(2,3)) minb(3)=-dedge(2,3)
         if (maxb(3) .gt. n(3)+dedge(1,3)) maxb(3)=n(3)+dedge(1,3)

c        loop over all grid points in the cube between minb and maxb.
         do k = minb(3), maxb(3) 
c          dielectric grid is offset +dedge(2,*) so it starts at 0,0,0
           kp = k + dedge(2,3)
           ptpdb(3) = edgesolv(2,3) - k*invgrid
           tdist3 = (pt(3) - ptpdb(3))**2
           do j = minb(2), maxb(2)
             jp = j + dedge(2,2)
             ptpdb(2) = edgesolv(2,2) - j*invgrid
             tdist2 = (pt(2) - ptpdb(2))**2
             do i = minb(1), maxb(1)
               ip = i + dedge(2,1)
               ptpdb(1) = edgesolv(2,1) - i*invgrid
               dist2 = (pt(1) - ptpdb(1))**2 + tdist2 + tdist3
               if (dist2 .le. size2) then
                 if (grid(ip,jp,kp) .le. threshold) then
                   if (grid(ip,jp,kp) .gt. 0) then
                     grid(ip,jp,kp) = grid(ip,jp,kp) - threshold -
     &                                                 decrement
                   else
                     grid(ip,jp,kp) = grid(ip,jp,kp) - decrement
                   endif
c                  adjust single accessible pt to grid value -2
c                  reserves -1 to indicate no covering accessible points
                   if (grid(ip,jp,kp) .eq. -1) then
                     grid(ip,jp,kp) = -2
                   endif
                 endif
               endif
             enddo
           enddo
         enddo

         return
       end

c------------------------------------------------------------------------
c   subroutine dielec
c   - compute the dielectric bitmap using receptor's SAS as the boundary
c------------------------------------------------------------------------
        subroutine dielec(anum,n,maxn,perang,probe,edgesolv,dedge)

          include 'pointer.h'

          integer anum, n(3), maxn, perang
          real probe, edgesolv(2,3)
          integer dedge(2,3)

          real discut, apt(3)
          integer p, increment, penalty
          integer memalloc

c         receptor coordinates and radii (dynamic)
          real pdbcrd(3,anum),rad(anum)
c         grid: dielectric grid (dynamic)
          integer grid(0:maxn,0:maxn,0:maxn)

c         allocate dynamic array
          i_grid = memalloc(i_grid, 4*(maxn+1)**3)

          increment = 10000
          penalty = 0
c         blot inside receptor atoms SAS
          do p = 1,anum
            discut = rad(p) + probe
            apt(1) = pdbcrd(1,p)
            apt(2) = pdbcrd(2,p)
            apt(3) = pdbcrd(3,p)
            call blot(discut,apt,increment,penalty,n,maxn,perang,
     &                edgesolv,dedge)
          enddo

          write(6,*) 'finished initializing dielectric grid'

          return
        end

c------------------------------------------------------------------------
c   subroutine recblot
c   - create molecular surface by placing water at every accessible point
c   - every grid point within 1.4 A of an accessible point set to 5
c------------------------------------------------------------------------
       subroutine recblot(accnum,n,maxn,perang,probe,edgesolv,dedge)
         include 'pointer.h'
         include 'acc.h'

         integer accnum, n(3), maxn, perang
         real probe, edgesolv(2,3)
         integer dedge(2,3)
         integer grid(0:maxn,0:maxn,0:maxn)
         real invgrid, probe2, dist2
         integer minb(3),maxb(3)
         real ptpdb(3)
         integer p, decrement, threshold
         real accpts(3,accnum)
         real apt(3)

         decrement = 1
         threshold = 10000
c        blot around all accessible points
         do p = 1,accnum
           apt(1) = accpts(1,p)
           apt(2) = accpts(2,p)
           apt(3) = accpts(3,p)
           call accblot(probe,apt,decrement,threshold,n,maxn,perang,
     &                  edgesolv,dedge)
         enddo

         write(6,*) 'finished creating receptor molecular surface grid'

         return
       end

c------------------------------------------------------------------------
c   subroutine walklig 
c   - computes the atomic SEV and then the percentage desolvation for
c       every ligand atom position in the grid
c   - algorithm summary:
c   0) Subroutine recblot is called beforehand to pre-generate 
c        the starting receptor-only coverage for each grid point
c   1) Loop over all grid points, placing ligand atom there
c     2) Find newly inaccessible receptor SAS points
c     3) If none, then immediately use receptor only ms
c     4) Otherwise, find accessible ligand SAS points
c     5) Add additional coverage due to all the accessible ligand SAS 
c          points (by decrementing grid point by 1 for each)
c     6) For newly inaccessible receptor SAS points, uncover 
c          (increment by 1) their grid points. If a grid point becomes
c          completely uncovered (grid = -1), set it to +7000.
c     7) Finally blot inside the ligand vdW surface, incrementing those 
c          points by 4000 if not already positive
c     8) Score the fractional desolvation using that surface
c     9) Reverse each step from 7) to 6) to 5) to regenerate the original
c          receptor-only grid as from 0)
c   10) Move to the next ligand position
c------------------------------------------------------------------------
       subroutine walklig(anum,accnum,maxn,lksize,perang,radmax,
     &                      born,probe,n,dedge,edgesolv,edgegrid)
         include 'pointer.h'
         include 'cubic.h'
         include 'acc.h'

         integer anum, accnum, maxn, lksize, perang
         real radmax, born, probe
         integer n(3), dedge(2,3)
         real edgesolv(2,3), edgegrid(2,3)

         integer i,j,k,ii,jj,kk,ip,jp,kp,m,a,p,q,r
         real ligp, ligp2, invgrid
         real cbln, lig(3)
         integer cblower(1), cbupper(1), cbatom(1)
         integer icb(3), ming(3), maxg(3), minb(3), maxb(3)
         logical reconly
         integer recnum, inaccnum, ind, liml, limu
         real dx1, dx2, dx3, d2

         integer lignum, gpt
         real apt(3)
         integer increment, penalty, decrement, threshold
         integer lower, upper

c        shared dynamic arrays
c        grid: dielectric grid (dynamic)
         integer grid(0:maxn,0:maxn,0:maxn)
c        solgrid: percentage desolvated grid (dynamic)
         real solgrid(0:maxn,0:maxn,0:maxn)
c        lksolv: lookup table for core solvmap kernel (dynamic)
         real lksolv(-lksize:lksize,-lksize:lksize,-lksize:lksize)
c        receptor coordinates and radii (dynamic)
         real pdbcrd(3,anum), rad(anum)
c        accessible points arrays (dynamic)
         integer accstart(anum), accstop(anum)
         real accpts(3,accnum)

c        debugging vars
         character*80 filename
         character*20 ftemp

         integer memalloc
c        local dynamic arrays
         pointer (i_radp2, radp2)
         pointer (i_radcheck, radcheck)
         pointer (i_recptr,recptr)
         pointer (i_inaccptr,inaccptr)
         pointer (i_ligpts,ligpts)
         real radp2(1), radcheck(1)
         integer recptr(1), inaccptr(1)
         real ligpts(3,1)

         ligp = born + probe
         ligp2 = ligp * ligp
         invgrid = 1.0/real(perang)

c        allocate dynamic arrays
         i_solgrid = memalloc(i_solgrid, 4*(maxn+1)**3)
         i_radp2 = memalloc(i_radp2, 4*anum)
         i_radcheck = memalloc(i_radcheck, 4*anum)
         i_recptr = memalloc(i_recptr, 4*anum)
c        XXX - inaccptr could grow by 1000 at a time instead?
         i_inaccptr = memalloc(i_inaccptr, 4*accnum)
         i_ligpts = memalloc(i_ligpts, 4*3*numpts)


c        initialize arrays
         do i = 1,anum
           radp2(i) = (rad(i) + probe)**2
           radcheck(i) = (born + 2.0*probe + rad(i))**2
         enddo

c        distance divide receptor atoms for fast intersection computations
c        cube length is maximum distance for ligand-receptor SAS intersection
         cbln = born + 2.0*probe + radmax
         call cubedata(1.0, edgegrid, cbln)
         i_cblower = memalloc(i_cblower,4*(cb(1)+1)*(cb(2)+1)*(cb(3)+1))
         i_cbupper = memalloc(i_cbupper,4*(cb(1)+1)*(cb(2)+1)*(cb(3)+1))
         i_cbatom = memalloc(i_cbatom,4*27*anum)
         call cube(anum,pdbcrd,rad)

         write (6,*) 'Starting solvmap calculation.'

         upper = (n(3)+1)*(n(2)+1)
c        main loops over all possible grid positions of the ligand atom
         do k = 0, n(3)
          write(6,'2(i4,a)') k+1,' / ',n(3)+1,' the way there!'
          call flush(6)
          lig(3) = edgesolv(2,3) - k*invgrid
          ming(3) = k - lksize
          maxg(3) = k + lksize
          if (ming(3) .lt. -dedge(2,3)) ming(3) = -dedge(2,3)
          if (maxg(3) .gt. n(3)+dedge(1,3)) maxg(3) = n(3)+dedge(1,3)
          icb(3) = (lig(3) - cmin(3))*cpa
          if (icb(3) .le. 0 .or. icb(3) .ge. cb(3)) then
            write (6,*) 'Outside distance cube, aborting!'
            stop
          endif
          do j = 0, n(2)
           lig(2) = edgesolv(2,2) - j*invgrid
           ming(2) = j - lksize
           maxg(2) = j + lksize
           if (ming(2) .lt. -dedge(2,2)) ming(2) = -dedge(2,2)
           if (maxg(2) .gt. n(2)+dedge(1,2)) maxg(2) = n(2)+dedge(1,2)
           icb(2) = (lig(2) - cmin(2))*cpa
           if (icb(2) .le. 0 .or. icb(2) .ge. cb(2)) then
             write (6,*) 'Outside distance cube, aborting!'
             stop
           endif
           do i = 0, n(1)
            lig(1) = edgesolv(2,1) - i*invgrid
            ming(1) = i - lksize
            maxg(1) = i + lksize
            if (ming(1) .lt. -dedge(2,1)) ming(1) = -dedge(2,1)
            if (maxg(1) .gt. n(1)+dedge(1,1)) maxg(1) = n(1)+dedge(1,1)
            icb(1) = (lig(1) - cmin(1))*cpa
            if (icb(1) .le. 0 .or. icb(1) .ge. cb(1)) then
              write (6,*) 'Outside distance cube, aborting!'
              stop
            endif

c           Achieived big speedup by storing the number of accessible
c           SAS points  covering each grid point (stored as -x-1 in 
c           the old dielectric grid itself). An accessible SAS point 
c           covers all grids points within probe radius of it.
c
c           grid values:
c           <=0 are the high dielectric points
c           0:      bulk solvent
c           <=-2:   -x-1 where x is # of accessible SAS points covering  
c
c           >=1 are the low dielectric points
c           10000:  inside receptor-only molecular surface
c           7000:   newly inaccessible point
c           4000-x: inside ligand vdW (coverd by x receptor acc. SAS)

            reconly = .false.
            recnum = 0
            inaccnum = 0
c           use cube to find all nearby receptor atoms 
            ind = 1+icb(1) + (cb(1)+1)*icb(2) +
     &           (cb(1)+1)*(cb(2)+1)*icb(3)
            liml = cblower(ind)
            limu = cbupper(ind)
            do m = liml,limu
             a = cbatom(m)
             dx1 = pdbcrd(1,a) - lig(1)
             dx2 = pdbcrd(2,a) - lig(2)
             dx3 = pdbcrd(3,a) - lig(3)
             d2 = dx1*dx1 + dx2*dx2 + dx3*dx3
c            ligand-receptor SAS intersection check
             if (d2 .lt. radcheck(a)) then
               recnum = recnum + 1
               recptr(recnum) = a
c              loop over accessible points of that receptor atom
               do p = accstart(a),accstop(a)
                dx1 = accpts(1,p) - lig(1)
                dx2 = accpts(2,p) - lig(2)
                dx3 = accpts(3,p) - lig(3) 
                d2 = dx1*dx1 + dx2*dx2 + dx3*dx3
c               ligand-accessible point intersection check
                if (d2 .lt. ligp2) then
c                 a receptor accessible point is now inaccessible
                  inaccnum = inaccnum + 1
                  inaccptr(inaccnum) = p
                endif
               enddo
             endif
            enddo

c           now have a list of intersecting receptor atoms recptr
c             and a list of newly inaccessible points inaccptr

c           fast path when ligand is far away or completely buried
            if (inaccnum .eq. 0) then
              reconly = .true.
              goto 270
            endif

            lignum = 0
c           find ligand accessible points
            do 200 q = 1,numpts
              do p = 1,3
                apt(p) = sphpts(p,q)*(born+probe) + lig(p)
              enddo

              do m = 1,recnum
                r = recptr(m)
                dx1 = pdbcrd(1,r) - apt(1)
                dx2 = pdbcrd(2,r) - apt(2)
                dx3 = pdbcrd(3,r) - apt(3)
                d2 = dx1*dx1 + dx2*dx2 + dx3*dx3
c               if inside any receptor atom, apt is inaccessible
                if (d2 .lt. radp2(r)) then
                  goto 200
                endif
              enddo
c             otherwise, ligand apt is accessible
              lignum = lignum + 1
              do p = 1,3
                ligpts(p, lignum) = apt(p)
              enddo
 200        continue

c           now have a list of the ligand accessible points ligpts

c           mark grid around ligand accessible points
            decrement = 1
            threshold = 0
            do q = 1,lignum
              apt(1) = ligpts(1,q)
              apt(2) = ligpts(2,q)
              apt(3) = ligpts(3,q)
              call accblot(probe,apt,decrement,threshold,n,maxn,perang,
     &                     edgesolv,dedge)
            enddo
c            filename = 'lacc.bmp'
c            call writegrid(15,n,maxn,perang,edgegrid,dedge,filename)

c           increment grid points near newly inaccessible points
            increment = 1
            penalty = 7000
            do m = 1,inaccnum
              q = inaccptr(m)
              apt(1) = accpts(1,q)
              apt(2) = accpts(2,q)
              apt(3) = accpts(3,q)
              call blot(probe,apt,increment,penalty,n,maxn,perang,
     &                  edgesolv,dedge)
            enddo
c            filename = 'inacc.bmp'
c            call writegrid(15,n,maxn,perang,edgegrid,dedge,filename)

c           blot inside the ligand vdW
            increment = 4000
            penalty = 0
            call blot(born,lig,increment,penalty,n,maxn,perang,
     &                edgesolv,dedge)
c            filename = 'lig.bmp'
c            call writegrid(15,n,maxn,perang,edgegrid,dedge,filename)

 270        continue

c           calculate solvmap for this point
            do kk = ming(3),maxg(3)
c             dielectric grid is offset +dedge(2,*) so it starts at 0,0,0
              kp = kk + dedge(2,3)
              do jj = ming(2),maxg(2)
                jp = jj + dedge(2,2)
                do ii = ming(1),maxg(1)
                  ip = ii + dedge(2,1)
                  if (grid(ip,jp,kp) .gt. 0) then
                    solgrid(i,j,k) = solgrid(i,j,k) + 
     &                                 lksolv(ii-i,jj-j,kk-k)
                  endif
                enddo
              enddo
            enddo

            if ( .not. reconly) then
c             reset grid to receptor only status
c             unblot inside ligand vdW
              decrement = 0
              threshold = 4000
              call accblot(born,lig,decrement,threshold,n,maxn,perang,
     &                     edgesolv,dedge)

c            filename = 'relig.bmp'
c            call writegrid(15,n,maxn,perang,edgegrid,dedge,filename)

c             re-decrement grid points near newly inaccessible points 
              decrement = 1
              threshold = 7000
              do m = 1,inaccnum
                q = inaccptr(m)
                apt(1) = accpts(1,q)
                apt(2) = accpts(2,q)
                apt(3) = accpts(3,q)
                call accblot(probe,apt,decrement,threshold,n,maxn,
     &                       perang,edgesolv,dedge)
              enddo
c            filename = 'reinacc.bmp'
c            call writegrid(15,n,maxn,perang,edgegrid,dedge,filename)

c             finally re-increment near ligand accessible points
              increment = 1
              penalty = 0
              do q = 1,lignum
                apt(1) = ligpts(1,q)
                apt(2) = ligpts(2,q)
                apt(3) = ligpts(3,q)
                call blot(probe,apt,increment,penalty,n,maxn,perang,
     &                    edgesolv,dedge)
              enddo

            endif

c            filename = 'reset.bmp'
c            call writegrid(15,n,maxn,perang,edgegrid,dedge,filename)

c           loop to next ligand position
           enddo
          enddo
         enddo

c        deallocate dynamic arrays
         i_cblower = memalloc(i_cblower,0)
         i_cbupper = memalloc(i_cbupper,0)
         i_cbatom = memalloc(i_cbatom,0)
         i_radp2 = memalloc(i_radp2,0)
         i_radcheck = memalloc(i_radcheck,0)
         i_recptr = memalloc(i_recptr, 0)
         i_inaccptr = memalloc(i_inaccptr, 0)
         i_ligpts = memalloc(i_ligpts, 0)

         return
       end
