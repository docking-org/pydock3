c     copyright 1981 UC Regents
c     author ID Kuntz
c     7/17/81 significant mod to condense by atom section
c     mod1: one surface
c     mod2:  spheres lie on proper side of tangent plane
c         only i,j pair with dotproduct
c         of ij vector and normal vector ge.0. are acceptable.
c     11/70 version
c
c     molecules are characterized by spheres.
c     entry point: connolly surface file with normals
c
c     two molecular types are used: receptor(R) and ligand (L)
c     the convex surface of R and the reentrant surface of L are used.
c
c	 edited to handle chain specifiers in ms/dms file 
c                Diana Roe 5/94
c  	read format changed to correspond to format of ms/dms
c				Diana Roe 6/17/94
c	minor mods to deal with HETATM residues indicated by a * in dms output
c				DAG 5/23/96
c	minor mods to fix bug introduced with chain specifier handling
c				DAG 5/28/96
c	bugfix for chain identifiers when residue numbers are > 999
c				Michael M Mysinger (MMM) 1/25/2012!
c
c     implicit none
c
c     parameters --
      integer nsmax
      parameter (nsmax = 300000)
c        nsmax:  maximun number of surface points.
      integer natmax
      parameter (natmax = 30000)
c        natmax:  maximum number of atoms.
      integer normax
      parameter (normax = 250000)
c        normax:  maximum number of surface normals.
      integer ncmax
      parameter(ncmax=200)
c        ncmax:  maximum number of clusters.
      integer nspmax
      parameter(nspmax=90000)
c        nspmax:  maximum number of spheres (bigger than radmin)
c        note: if nspmax > 99999, change format 651 (and in main.f)
c
c     variables --
      integer nat
c        nat:  number of atoms. 
      integer nsur
c        nsur:  number of surface points. 
      integer isph
c        isph:  total number of spheres for receptor.
      integer ncm
c        ncm:  number of clusters.
      integer nn,i,j
c        nn:  do loop variable.
c        i:  do loop variable.
c        j:  do loop variable.
      integer it
c        it: index into itemp (sphere) array.
      integer nspall
c        nspall: number of all spheres big enough to write.
        
      real dotlim
c        dotlim:  limit on dot product between two surface
c                 normals.  product must be greater that dotlim
c                 for a sphere to be generated.
      real radmax
c        radmax:  maximum radius allowed.
      character*1 surftp
c        surftp:  tag for the type of molecule--
c                 R for receptor, and L for ligand.
      character*1 dentag
c        dentag:  tag for lower density of surface points--
c                 X for all.
         character*80 nma
c            nma:  surface file name.
         character*80 nmb
c            nmb:  output file name. 
c
c     array list --
      integer isatm(nsmax)
c        isatm: atom # for surface point.
      integer inorm(nsmax)
c        inorm: snorm type for each surface point.
      integer iatoms(natmax)
c        iatoms:  number of first atom associated with each
c                 sphere saved.
      integer jatoms(natmax)
c        jatoms:  number of second atom associated with each
c                 sphere saved.
      integer itemp(ncmax,natmax)
c        itemp:  temporary store of spheres in cluster.
      integer nsphc(ncmax)
c        nsphc:  number of spheres in each cluster.
      integer nct(natmax)
c        nct:  set to zero when sphere is part of a cluster.
      integer iclus(natmax)
c        iclus:  cluster # as function of atom (line) #
      integer iall(nspmax)
c        iall:  number of first atom associated with each
c                 sphere kept.
      integer jall(nspmax)
c        jall:  number of second atom associated with each
c                 sphere kept.
      real spcors(natmax,3)
c        spcors:  coordinates for sphere centers of each sphere
c                 saved.
      real rads(natmax)
c        rads:  radius of each sphere saved.
      real scor(nsmax,3)
c        scor: surface coordinates.
      real snorm(normax,3)
c        snorm: normal vector NOTE classified by type.
      real allspc(nspmax,3)
c        allspc:  coordinates for sphere centers of each sphere
c                 kept.
      real allrad(nspmax)
c        allrad:  radius of each sphere kept.
      REAL RADMIN
c
c     open unit 5 (INSPH) and unit 6 (OUTSPH)
      open(unit=5,file='INSPH',status='old')
      open(unit=6,file='OUTSPH',status='new')
c     read in sequence
c
   90 continue
      read(5,100) nma
  100 format(a80)
      open(4,file=nma,status='old')
      rewind(4)
      read(5,102) surftp
  102 format(a1)
      if(surftp.ne.'L'.and.surftp.ne.'R') goto 90
c
      read(5,150) dentag
  150 format(a1)
      write(6,152) dentag
  152 format(' density type = ',a1)
c
      write(6,154) nma, surftp
  154 format(' reading  ',a80,'   type   ',a1)
c
      call reader(inorm,scor,isatm,snorm,dentag,normax,nsmax,natmax,
     &            nat,nsur,surftp)
c
c     end of read in
c
c     find spheres
c
      write(6,160) nma
  160 format(' finding spheres for   ',a80)
      read(5,*) dotlim
      if(dotlim.lt.-1..or.dotlim.gt.1.) dotlim = 0.0
      write(6,172) dotlim
  172 format(' dotlim =  ',f8.3)
c
      read(5,*) radmax
      if(radmax.le.0) radmax = 5.0
      write(6,176) radmax
  176 format(' radmax = ',f8.3)
        WRITE(6,*)'Minimum radius of acceptable spheres?'
        READ(5,*)RADMIN
        WRITE(6,*)RADMIN
      IF(RADMIN.GE.RADMAX)THEN
          WRITE(6,*)'RADMAX MUST BE > THAN RADMIN - TRY AGAIN'
          WRITE(6,*)'SORRY CHARLIE, PROGRAM BOMBS'
          STOP
      ENDIF
c
c     'output dump'
      read(5,100) nmb
      write(6,178) nmb
  178 format(' output to  ',a80)
      open(9,file=nmb,status='new')
c
      isph = 0
      call sphere(inorm,isatm,scor,snorm,dotlim,radmax,RADMIN,
     &spcors,iatoms,jatoms,rads,normax,nsmax,natmax,isph,surftp,nat,
     &nsur,nspmax,nspall,iall,jall,allrad,allspc)
c
c
      call clus(rads,spcors,iatoms,jatoms,itemp,nsphc,ncm,natmax,
     &             iclus,nct,isph)
c
      write (9, 500)
  500 format('DOCK 3.5 receptor_spheres')

c     color table is empty - don't write anything
c
c     dump clusters
      do 670 nn = 1,ncm
        write(9,651) nn,nsphc(nn)
  651   format('cluster ',i5,'   number of spheres in cluster ',i5)
c       zero flag and zero color number for all spheres
        do 660 i = 1,nsphc(nn)
          it = itemp(nn,i)
          write(9,652) iatoms(it),(spcors(it,j),j=1,3),rads(it),
     &    jatoms(it),0,0
  652     format(i5,3f10.5,f8.3,i5,i2,i3)
  660   continue
  670 continue
c
c     dump all spheres
      write(9,651) 0,nspall
c     (was 1 originally, but 0 makes more sense to me - MLC)
      do 680 i = 1,nspall
c       zero flag and zero color number for all spheres
        write(9,652) iall(i),(allspc(i,j),j=1,3),allrad(i),jall(i),0,0
  680 continue
      close(9)
c
      end
c-----------------------------------------------------------------------
      subroutine reader(inorm,scor,isatm,snorm,dentag,normax,
     &  nsmax,natmax,nat,nsur,surftp)
c
c     reader accepts connolly surface files and converts them to a
c        stripped file containing only contact or reentrant points
c        as required by file type(surftp), R or L.
c
c     implicit none
c
c     parameters --
      integer nsmax
c        nsmax:  maximun number of surface points.
      integer natmax
c        natmax:  maximum number of atoms.
      integer normax
c        normax:  maximum number of surface normals.
c
c     variables --
      integer nat
c        nat:  number of atoms.
      integer nsur
c        nsur:  number of surface points.
      integer nnorm
c        nnorm: number of normals in snorm.
      integer nres
c        nres:  residue number.
      integer pnres 
c        pnres:  previous residue number
      integer resincr
c        resincr: amount to increment residue to add for chain specifier
      integer i, j, jj, ino
c        i: counter for surface points.
c        j: counter for atoms.
c        jj: do loop variable.
c        ino: do loop variable.
      character*1 chain
c		chain- chain specifier to be read in
      character*1 pchain
c		previous chain specifier read in.
      character*3 resnam
c        resnam: res # for atom.
      character*4 atnam
c        atnam: atom name.
      character*3 tag
c        tag:  indicator in ms file that tells whether point
c              is an atom (A), reentrant surface (SR), contact
c              surface (SC), or saddle surface (SS).
      character*1 surftp
c        surftp:  tag for the type of molecule--
c                 R for receptor, and L for ligand.
      character*1 dentag
c        dentag:  tag for lower density of surface points--
c                 X for all.
      character*80 line
c
c     array list --
      integer isatm(nsmax)
c        isatm: atom # for surface point.
      integer inorm(nsmax)
c        inorm: snorm type for each surface point.
      real scor(nsmax,3)
c        scor: surface coordinates.
      real snorm(normax,3)
c        snorm: normal vector NOTE classified by type.
      real c(3)
c        c: surface point coordinates if tag(1:1) = 'A', or
c           atom coordinates if tag(1:1) = 'S'.
      real an(3)
c        an: coordinates for the surface normal.
c
c     open temp file
c
      open(1,file='temp1.ms',status='new')
      rewind(1)
c
c     temp file for atomic coords
      open(3,file='temp3.atc',status='new')
      rewind(3)
c     index surface list using i; atom list with j
c
      i=0
      j=0
      nnorm=0
      resincr = 0
      pnres = 0
c
c     loop starts here
c
   10 continue
      read(4,'(A)',end=200) line
      read(line,12,err=2000) resnam,nres,atnam,(c(jj),jj=1,3),
     &tag,(an(jj),jj=1,3)
   12 format(a3,i6,a4,f8.3,2f9.3,1x,a3,7x,3f7.3)
	 goto 2001
c
c     section to read chain specifier and add amount to
c     nres to correct for it.
c     adjusted to handle nres > 999 (MMM)
c
2000  read(line,2002,err=2010) resnam,nres,chain,atnam,
     &(c(jj),jj=1,3),tag,(an(jj),jj=1,3)
2002  format(a3,i5,a1,1x,a4,f8.3,2f9.3,1x,a3,7x,3f7.3)
      goto 2040
2010  read(line,2012,err=2014) resnam,nres,chain,atnam,
     &(c(jj),jj=1,3),tag,(an(jj),jj=1,3)
2012  format(a3,i4,a1,1x,a4,f8.3,2f9.3,1x,a3,7x,3f7.3)
      goto 2040
c
c     this section added to deal with HETATM * residues, e.g. "ZN   1A*" (DAG)
c
2014  read(line,2016) resnam,nres,chain,atnam,(c(jj),jj=1,3),
     &tag,(an(jj),jj=1,3)
2016  format(a3,i3,a1,2x,a4,f8.3,2f9.3,1x,a3,7x,3f7.3)


2040  if ( pchain.ne.chain) then
          resincr = pnres
          pchain = chain
      endif
      pnres = nres
      nres=nres+resincr
c
c     end of chain specifier section
c
2001  continue
      if(tag.eq.'A  ') then
         j=j+1
      if(j+1.gt.natmax) then
         write(6,13) natmax
   13    format(' warning exceeded # of atoms ',i5)
         stop
      endif

c
      write(3,15) j,nres,resnam,atnam,(c(jj),jj=1,3)
   15 format(i6,i4,a3,1x,a4,1x,3f8.3)
c
      elseif(dentag.eq.'X'.or.tag(3:3).eq.dentag) then
         i=i+1
c
c     snorm classification 11/70
c
         do 105 ino=1,nnorm
            do 102 jj=1,3
               if(an(jj).ne. snorm(ino,jj)) goto 105
  102       continue
            inorm(i)=ino
            goto 109
c
  105    continue
         nnorm=nnorm +1
         if(nnorm.gt.normax) then
            write(6,106) nnorm,normax
  106       format(' nnorm = ',i5, ' normax = ',i5)
            write(6,107)
  107       format(' normax exceeded, program stops')
            stop
c
         else
            do 108 jj=1,3
               snorm(nnorm,jj)=an(jj)
  108       continue
            inorm(i)=nnorm
         endif
c
  109    continue
c
         isatm(i)=j
         do 110 jj=1,3
            scor(i,jj)=c(jj)
  110    continue
         write(1,12) resnam,nres,atnam,(c(jj),jj=1,3),tag,
     &(an(jj),jj=1,3)
      endif
c
      if(i+1.gt.nsmax) then
         write(6,115) nsmax
  115    format(' warning exceeded # of surface points ',i5)
         stop
      endif
c
      goto 10
c
c     termination
  200 continue
      nat=j
      nsur=i
      write(6,210) nat,nsur
  210 format(' # of atoms =  ',i5,'   # of surf pts =  ',i5)
      close(4)
c
c     reverse signs if surftp is ligand
c
      if(surftp.eq.'L') then
         do 300 i=1,nnorm
         do 300 jj=1,3
            snorm(i,jj)=-snorm(i,jj)
  300    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine sphere(inorm,isatm,c,snorm,dotlim,radmax,RADMIN,
     &spcors,iatoms,jatoms,rads,normax,nsmax,natmax,isph,surftp,nat,
     &nsur,nspmax,nspall,iall,jall,allrad,allspc)
c
c     sphere accepts stripped files from reader
c       it uses an analytical solution to find the largest
c       sphere that lies tangent to a given surface point,i
c       and tangent to a second surface point, j
c       the largest such sphere for each atom is then found
c
c     implicit none
c
c     parameters --
      integer nsmax
c        nsmax:  maximun number of surface points.
      integer natmax
c        natmax:  maximum number of atoms.
      integer normax
c        normax:  maximum number of surface normals.
      integer nspmax
c        nspmax:  maximum number of spheres.
c
c     variables --
      integer nat
c        nat:  number of atoms.
      integer nsur
c        nsur:  number of surface points.
      integer nspall
c        nspall: number of all spheres big enough to write.
      integer nres
c        nres:  residue number.
      integer iatom
c        iatom:  atom number of first surface point.
      integer jatom
c        jatom:  atom number of second surface point.
      integer izero
c        izero:  counter for number of zero surface normal
c                components.
      integer ipair
c        ipair:  store for surface point for smallest sphere
c                so far.
      integer ick
c        ick:  indicator for whethet a sphere should be written
c              ick = 0 don't write, ick = 1 write.
      integer ipt, jpt
c        ipt:  surface point number for atom i.
c        jpt:  surface point number for atom j.
      integer js
c        js:  used in special case calculation.  for special case 2
c             zeros js is the surface normal component that has the
c             value 1 (one).  for the special case 1 zero js is the
c             surface normal component that has the value 0 (zero).
      integer ja, jb
c        ja:  used in special case calculation.  for special case 2
c             zeros ja is one of the 0 (zero) components.  for special
c             case 1 zero ja is one of the non-zero surface normal
c             components.
c        jb:  used in special case calculation.  for special case 2
c             zeros jb is one of the 0 (zero) components.  for special
c             case 1 zero jb is one of the non-zero surface normal
c             components.
      integer iatnew
c        iatnew:  atom number of surface point i.
      integer jres
c        jres:  residue number of residue j.
      integer itest
c        itest:  abs(nres-jres) used to test how far apart in sequence
c                atoms i and j are.
      integer isph
c        isph: counter for spheres.
      integer i, j, jj
c        i: do loop variable.
c        j: do loop variable.
c        jj: do loop variable.
         real s
c        s:  sum of del's.
      real sumi, sumj
c        sumi:  sum of squared coordinates for surface point i.
c        sumj:  sum of squared coordinates for surface point j.
      real dotck
c        dotck:  dot product of surface normal of point i and a
c                vector from point i to point j.  must be less
c                than zero to generate a sphere.
      real dot
c        dot:  dot product of the surface normals of points i and j.
      real dotp
c        dotp:  dot divided by the radius of the sphere from i and j.
c               only saved for the smallest sphere from i to any j.
      real denom
c        denom:  the denominator for the sphere calculation.
      real dotlim
c        dotlim:  limit on dot product between two surface
c                 normals.  product must be greater that dotlim
c                 for a sphere to be generated.
      real radmax
      REAL RADMIN
c        radmax:  maximum radius allowed.
      real slop, bet
c        slop:  intermediate in the sphere calculation for the special
c               cases.  analogous to slope for the usual case.
c        bet:  intermediate in the sphere calculation for the special
c              cases.  analogous to beta for the usual case.
      real rad
c        rad:  raduis of the smallest sphere for a given surface point.
      real csq
c        csq:  used in special case calculation.  csq is the sum of
c              the squares of the non-zero surface normal components.
      real cp
c        cp:  used for special case 2 zeros. this is the new center
c             point coordinate in the dimension with the non-zero
c             surface normal component.
      real cpa, cpb
c        cpa:  used for special case 1 zero. this is the new center
c              point coordinate in one of the dimensions with a
c              non-zero surface normal component.
c        cpb:  used for special case 1 zero. this is the new center
c              point coordinate in one of the dimensions with a
c              non-zero surface normal component.
      real xp, yp, zp
c        xp: temp store of x coordinate of sphere center.
c        yp: temp store of y coordinate of sphere center.
c        zp: temp store of z coordinate of sphere center.
      real dis
c        dis: temp store of radius.
      real rmin, rtmin, rtmina, rtminb, rtminc
c        rtmin:  minimum radius encountered for a surface point.
c        rmin:  minimum radius squared.
c        rtmina:  minimum radius * 3.5.
c        rtminb:  minimum radius * 2.83.
c        rtminc:  minimum radius * 2.0.
      real rtmax
c        rtmax:  maximum sphere radii saved for a particular atom.
      character*3 resnam
c        resnam: res # for atom.
      character*3 rnam
c        rnam: res # for atom used in read of temp1.
      character*4 atnam
c        atnam: atom name.
      character*4 anam
c        anam: atom name used in read of temp1.
      character*3 tag
c        tag:  indicator in ms file that tells whether point
c              is an atom (A), reentrant surface (SR), contact
c              surface (SC), or saddle surface (SS).
      character*1 surftp
c        surftp:  tag for the type of molecule--
c                 R for receptor, and L for ligand.
c
c     array list --
      integer isatm(nsmax)
c        isatm: atom # for surface point.
      integer inorm(nsmax)
c        inorm: snorm type for each surface point.
      integer iatoms(natmax)
c        iatoms:  number of first atom associated with each
c                 sphere saved.
      integer jatoms(natmax)
c        jatoms:  number of second atom associated with each
c                 sphere saved.
      integer iall(nspmax)
c        iall:  number of first atom associated with each
c                 sphere kept.
      integer jall(nspmax)
c        jall:  number of second atom associated with each
c                 sphere kept.
      real spcors(natmax,3)
c        spcors:  coordinates for sphere centers of each sphere
c                 saved.
      real rads(natmax)
c        rads:  radius of each sphere saved.
      real snorm(normax,3)
c        snorm: normal vector NOTE classified by type.
      real c(nsmax,3)
c        c: surface point coordinates (scor)
      real allspc(nspmax,3)
c        allspc:  coordinates for sphere centers of each sphere
c                 kept.
      real allrad(nspmax)
c        allrad:  radius of each sphere kept.
      real an(3)
c        an: coordinates for the surface normal.
      real slope(2), beta(2)
c        slope:  intermediate in the sphere calculation.
c        beta:  intermediate in the sphere calculation.
      real cs(3)
c        cs:  coordinates for the sphere center.
      real del(3)
c        del:  j coordinate minus i coordinate.
      real sn(3)
c        sn:  surface normals for point i.
      real spt(3)
c        spt:  surface point coordinates for point i.
      real sp(3)
c        sp:  sphere center coordinates
      real spcor(3)
c        spcor:  sphere center coordinates for the largest sphere
c                associated with a particular atom.
      real anorm(3)
c        anorm:  surface normal coordinates for point i.
      real acor(3)
c        acor:  atom coordinates.

      logical fpass

c     formulas are the same for ligands and receptors
c     EXCEPT FOR SIGN OF NORMALS (FIXED IN READER)

c     rewind temp file written by reader

      rewind(1)

c     open temp file for writting by sphere
      open(2,file='temp2.sph',status='new')
      rewind(2)

c     initialize number of spheres kept to 0
      nspall = 0
c     outer loop

      do 300 i=1,nsur

c     sums caln each time for 11/70
c     snorms packed into sn for 11/70
c
         sumi=c(i,1)**2+c(i,2)**2+c(i,3)**2
c
         do 103 jj=1,3
            sn(jj)=snorm(inorm(i),jj)
  103    continue
c
c       rmin=radius**2; rtmin=radius
c
         rmin=1000000.
         rtmin=1000.
c     rtminx's for use in cubing section
         rtmina=3.5*rtmin
         rtminb=2.83*rtmin
         rtminc=2.0*rtmin
         iatom=isatm(i)
c
c     check normals for 6 special cases
c
         izero=0
         do 105 jj=1,3
            if(sn(jj).eq.0.) izero=izero+1
  105    continue
         if(izero.eq.1) goto 211
         if(izero.eq.2) goto 201
         if(izero.eq.3) then
           write(6,107) i
  107      format(' SPHGEN: zero normal for surface point  ',i7)
           write (6,108)
  108      format('program stops')
c           if(izero.eq.3) goto 300
           stop
         endif
c
c     otherwise, continue
c
         slope(1)=sn(2)/sn(1)
         slope(2)=sn(3)/sn(1)
         beta(1) =c(i,2)-slope(1)*c(i,1)
         beta(2) =c(i,3)-slope(2)*c(i,1)
c
c     inner loop
c       i,j probably not symmetric
c
         do 150 j=1,nsur
            if(i.eq.j) goto 150
c         skip other surface points on same atom
            if(iatom.eq.isatm(j)) goto 150
c
            sumj=c(j,1)**2+c(j,2)**2+c(j,3)**2
c
c         cubing filter
c
            s=0.
            do 120 jj=1,3
               del(jj)=c(j,jj)-c(i,jj)
               s=s+abs(del(jj))
  120       continue
c
            if(s.gt.rtmina) goto 150
c
c     check for proper side of tangent plane
            dotck=0.
            do 130 jj=1,3
               dotck=dotck+del(jj)*sn(jj)
  130       continue
c
            if(dotck.lt.0.) goto 150
c
            denom=2.*(del(1)+del(2)*slope(1)+del(3)*slope(2))
            if(denom.ne.0.) goto 135
c           write(6,132) i,j
c           format(' convex hull',i5,2x,i5)
            goto 150
c
  135       continue
c
            xp=(sumj-sumi-2.*(del(2)*beta(1)+del(3)*beta(2)))/denom
            yp=slope(1)*xp + beta(1)
            zp=slope(2)*xp + beta(2)
c
            dis=(xp-c(i,1))**2+(yp-c(i,2))**2+(zp-c(i,3))**2
            if(dis.lt.rmin) then
               rmin=dis
               rtmin=sqrt(dis)
               rtmina=3.5*rtmin
               ipair=j
               cs(1)=xp
               cs(2)=yp
               cs(3)=zp
            endif
c
  150    continue
c     end of inner loop
c     skip around special cases
         goto 250
c
c     special cases
c     two surface normals equal zero
c
  201    continue
c
c     find special coord
c
         do 202 jj=1,3
            cs(jj)=c(i,jj)
            if(abs(sn(jj)).gt..99) js=jj
  202    continue
c
         if(js.eq.1) then
            ja=2
            jb=3
         elseif(js.eq.2) then
            ja=1
            jb=3
         else
            ja=1
            jb=2
         endif
c
         csq=c(i,js)**2
c
c     inner loop
         do 209 j=1,nsur
            if(i.eq.j) goto 209
            if(iatom.eq.isatm(j)) goto 209
            do 204 jj=1,3
               del(jj)=c(j,jj)-c(i,jj)
  204       continue
c
            if(abs(del(js)).gt.rtminc) goto 209
c        check for proper side of tangent plane
            dotck=0.
            do 205 jj=1,3
               dotck=dotck+del(jj)*sn(jj)
  205       continue
c
            if(dotck.lt.0.) goto 209
c
c
c     symmetry related points can fail here if del(js)=0.
            if(del(js).eq.0.) goto 209
c
            cp=(c(j,js)**2-csq+del(ja)**2+del(jb)**2)/(2.*del(js))
            dis=(cp-c(i,js))**2
c
            if(dis.lt.rmin) then
               rmin=dis
               rtmin=abs(cp-c(j,js))
               rtminc=2.*rtmin
               ipair=j
               cs(js)=cp
               do 206 jj=1,3
                  if(js.eq.jj) goto 206
                  cs(jj)=c(i,jj)
  206          continue
c
            endif
c
  209    continue
         goto 250
c
c     special case 1 surface noraml zero
c
  211    continue
c
c     find special coord
c
         do 212 jj=1,3
            if(sn(jj).eq.0.) js=jj
  212    continue
         if(js.eq.1) then
            ja=2
            jb=3
         elseif(js.eq.2)then
            ja=1
            jb=3
         elseif(js.eq.3) then
            ja=1
            jb=2
         endif
c
         slop=sn(jb)/sn(ja)
         bet=c(i,jb)-slop*c(i,ja)
c
         csq=c(i,ja)**2+c(i,jb)**2
c
c     inner loop
c
         do 219 j=1,nsur
            if(i.eq.j) goto 219
            if(iatom.eq.isatm(j)) goto 219
            do 214 jj=1,3
               del(jj)=c(j,jj)-c(i,jj)
  214       continue
            s=abs(del(ja))+abs(del(jb))
            if(s.gt.rtminb) goto 219
c
c        check for proper side of tangent plane
            dotck=0.
            do 215 jj=1,3
               dotck=dotck+del(jj)*sn(jj)
  215       continue
c
            if(dotck.lt.0.) goto 219
c
            denom=2.*(del(ja)+slop*del(jb))
c        patch for symmetry related atoms
            if(del(ja).eq.0..and.del(jb).eq.0.)  then
               dis=del(js)**2
               cpa=c(i,ja)
               cpb=c(i,jb)
c
            endif
            if(denom.eq.0.) then
               goto 216
            endif
            cpa=(c(j,ja)**2+c(j,jb)**2-csq+del(js)**2-2.*bet*
     &del(jb))/denom
            cpb=slop*cpa+bet
c
            dis=(cpa-c(i,ja))**2 +(cpb-c(i,jb))**2
c
  216       continue
c
            if(dis.lt.rmin) then
               rmin=dis
               rtmin=sqrt(dis)
               rtminb=2.83*rtmin
               cs(js)=c(i,js)
               cs(ja)=cpa
               cs(jb)=cpb
               ipair=j
            endif
c
  219    continue
c
c     end of inner loop
c     collect results: min radius assoc with i,j, dot prods.
  250    continue
c
c     dotprod is for vectors from i to sphere and sphere to j
c
         dot=0.
         do 251 jj=1,3
            dot=dot+(cs(jj)-c(i,jj))*(c(ipair,jj)-cs(jj))
  251    continue
         dotp=dot/rmin
c
c     read temp1 file to get nres
c
         read(1,260) resnam,nres,atnam,(spt(jj),jj=1,3),tag,
     &(an(jj),jj=1,3)
  260    format(a3,i5,1x,a4,1x,2f8.3,f9.3,1x,a3,7x,3f7.3)
c  260    format(a3,i5,1x,a4,1x,f7.3,2f9.3,1x,a3,7x,3f7.3)
c     store this in inorm
         inorm(i)=nres
c     write this to temp2
         if(dotp.gt.dotlim.AND.RTMIN.LT.999.0) then
            write(2,270) i,nres,rtmin,ipair,dotp,(cs(jj),jj=1,3)
  270       format(2i5,f8.4,i8,4f12.6)
         endif
c
  300 continue
c     end of outer loop
      close(1)
      close(2)
      close(3)
c
c     condense list by atom
c
c     rewind atom coord file
      open(1,file='temp1.ms',status='old')
      open(2,file='temp2.sph',status='old')
      open(3,file='temp3.atc',status='old')
      rewind(1)
      rewind(2)
      rewind(3)
      fpass=.true.
c
c     entry point: read atom file
  350 continue
      read(3,355,end=500) iatom,nres,resnam,atnam,(acor(jj),jj=1,3)
  355 format(i6,i4,a3,1x,a4,1x,3f8.3)
      rtmax=0.
      ick=0
      if(fpass) goto 360
c
c     goto test point
      goto 380
c
c     entry point: read other files
c
  360 continue
      read(1,260,end=490) rnam,nres,anam,(spt(jj),jj=1,3),tag,
     &(an(jj),jj=1,3)
      read(2,270,err=488,end=490)ipt,nres,rad,jpt,dotp,(sp(jj),jj=1,3)
      fpass=.false.
c
c     test point
  380 continue
      iatnew=isatm(ipt)
      if(iatnew.eq.iatom) then
         jres=inorm(jpt)
c all spheres (not just one per atom) bigger than radmin are kept
c in special arrays for writing at end of sphere cluster file (MLC)
         if(rad .ge. radmin) then
            nspall = nspall + 1
            if (nspall .gt. nspmax) then
              write (6,383) nspmax
  383         format('nspmax(',i8,') exceeded, program stops')
              stop
            endif
            iall(nspall) = iatom
            jall(nspall) = isatm(jpt)
            allrad(nspall) = rad
            allspc(nspall,1) = sp(1)
            allspc(nspall,2) = sp(2)
            allspc(nspall,3) = sp(3)
         endif
         itest=abs(nres-jres)
*        if(rad.gt.rtmax.and.(itest.gt.4.or.surftp.eq.'L')) then
         if(rad.gt.rtmax)then
            rtmax=rad
            ick=1
            jatom = isatm(jpt)
c
            do 390 jj=1,3
               spcor(jj)=sp(jj)
               anorm(jj)=an(jj)
  390       continue
c
         endif
      else
c
c     dump old values if ick =1 (some point met tests)
         if(ick.eq.1.and.rtmax.lt.radmax.AND.RTMAX.GE.RADMIN) then
*        if(ick.eq.1.and.rtmax.lt.radmax) then
            isph = isph + 1
            iatoms(isph) = iatom
            jatoms(isph) = jatom
            rads(isph) = rtmax
            spcors(isph,1) = spcor(1)
            spcors(isph,2) = spcor(2)
            spcors(isph,3) = spcor(3)
         endif
c
      endif
c
c     if files are on same atom or surface file is behind read surface file
c
      if(iatnew.le.iatom) goto 360
c
c     if atom file is behind, read atom file new atom
c
      if(iatnew.gt.iatom) goto 350
c
c     end of surface file, force dump if ick=1
  488 continue
      write(6,489)
  489 format(' error reading temp2.sph')
      write(6,270) ipt,nres,rad,jpt,dotp,(sp(jj),jj=1,3)
      backspace(2)
      read(2,270) ipt,nres,rad,jpt,dotp,(sp(jj),jj=1,3)
      write(6,270) ipt,nres,rad,jpt,dotp,(sp(jj),jj=1,3)
  490 continue
      if(ick.eq.1) then
         isph = isph + 1
         iatoms(isph) = iatom
         jatoms(isph) = isatm(jpt)
         rads(isph) = rtmax
         spcors(isph,1) = spcor(1)
         spcors(isph,2) = spcor(2)
         spcors(isph,3) = spcor(3)
      endif
c
c     end of atom file: exit
  500 continue
      close(3,status='delete')
      close(2,status='delete')
      close(1,status='delete')
c
      return
      end
c-----------------------------------------------------------------------
      subroutine clus(rads,spcors,iatoms,jatoms,itemp,nsphc,ncm,
     &                natmax,iclus,nct,isph)
c     finds overlapping spherical pairs (separation <= sum of radii)
c     then uses build up principle to find clusters
c     writes out clusters in order of size (# of spheres)
c
c     implicit none
c
c     parameters --
      integer natmax
c        natmax:  maximun number of spheres allowed.
      integer npmax
      parameter(npmax=200000)
c        npmax:  maximum number of spheres pairs. 
      integer ncmax
      parameter(ncmax=200)
c        ncmax:  maximum number of clusters.
c
c     variables --
      integer isph
c        isph:  total number of spheres for receptor.
      integer npair
c        npair:  number of overlapping sphere pairs.
      integer nclus
c        nclus:  counter for clusters.
      integer ict
c        ict:  counter for number of spheres in cluster.
      integer ncm
c        ncm:  number of clusters.
      integer it
c        it:  counter for number of spheres in each cluster.
      integer nlar
c        nlar:  number of spheres in largest cluster not yet
c               stored in itemp.
      integer np
c        np:  counter for number of overlapping sphere pairs.
      integer nn,n,j,jj,i
c        nn:  do loop variable.
c        n:  do loop variable.
c        j:  do loop variable.
c        jj:  do loop variable.
c        i:  do loop variable.
      real dis
c        dis:  distance between sphere centers to be compared with
c              their combined radii in overlap test.
c
c     array list --
      integer iatoms(natmax)
c        iatoms:  number of first atom associated with sphere.
      integer jatoms(natmax)
c        jatoms:  number of second atom associated with sphere.
      integer iclus(natmax)
c        iclus:  cluster # as function of atom (line) #.
      integer ipair(npmax,2)
c        ipair:  line numbers of pairs of spheres.
      integer nct(natmax)
c        nct:  set to zero when sphere is part of a cluster.
      integer itemp(ncmax,natmax)
c        itemp:  atom numbers per cluster for output.
      integer nsphc(ncmax)
c        nsphc:  number of spheres in each cluster.
      real rads(natmax)
c        rads:  radius of sphere.
      real spcors(natmax,3)
c        spcors:  sphere center coordinates.

c     find spherical pairs
      np=0
      do 110 i=1,isph-1
         do 100 j=i+1,isph
            dis=0.
            do 90 jj=1,3
               dis=dis+(spcors(i,jj)-spcors(j,jj))**2
   90       continue
            dis=sqrt(dis)
c
c     test
            if(dis.le.rads(i)+rads(j)) then
               np=np+1
               if(np.gt.npmax) then
                  write(6,95) np,npmax
   95             format(' WARNING:exceeded max # pairs',i5,
     &'-pairs used',i5,'-maximum pairs allowed')
               endif
               ipair(np,1)=i
               ipair(np,2)=j
            endif
  100    continue
  110 continue
c
      npair=np
c
c     clustering
c        find first cluster
      nclus=1
      i=ipair(1,1)
      j=ipair(1,2)
      iclus(i)=nclus
      iclus(j)=nclus
      nct(nclus)=2
c
c       scan list of pairs for others containing i,j or their partners
c
  199 continue
c
      do 200 np=1,npair
         i=ipair(np,1)
         j=ipair(np,2)
c
         if(iclus(i).eq.0.and.iclus(j).ne.0) iclus(i)=nclus
         if(iclus(i).ne.0.and.iclus(j).eq.0) iclus(j)=nclus
  200 continue
c
c     count members of this cluster
c
      ict=0
      do 300 i=1,isph
         if(iclus(i).eq.nclus) ict=ict+1
  300 continue
c
c     if ict is the same as last round, find the next cluster.
      if(ict.eq.nct(nclus)) goto 399
      nct(nclus)=ict
      goto 199
c     else, the cluster is complete
c
c     find next cluster
  399 continue
c
      do 400 np=1,npair
         i=ipair(np,1)
         j=ipair(np,2)
         if(iclus(i).eq.0) goto 410
  400 continue
c
c     no zero entries: loop completed  goto 500 for dump
      goto 500
c
c     some non assigned spheres left,set up new cluster and continue
c
  410 continue
      nclus=nclus+1
c     check max clusters
      if(nclus.gt.ncmax) write(6,420)
  420 format (' WARNING: exceeds max clusters')
c
      iclus(i)=nclus
      iclus(j)=nclus
c
      nct(nclus)=2
      goto 199
c
c     loop complete
  500 continue
      write(6,510) nclus
  510 format (' clustering is complete  ',i5,'  clusters')
c
      ncm=nclus
c
c     ranking clusters by size
c
      do 700 nn=1,ncm
         nlar=0
c
         do 600 n=1,ncm
            if(nct(n).gt.nlar) then
               nlar=nct(n)
               nclus=n
            endif
  600    continue
c
c     find all members of this cluster
c
         it=0
         do 650 i=1,isph
            if(iclus(i).eq.nclus) then
               it=it+1
               itemp(nn,it)=i
            endif
            nsphc(nn)=it
  650    continue
c
c     set nct(n) to zero so it won't be found again
         nct(nclus)=0
c
  700 continue
      return
      end
