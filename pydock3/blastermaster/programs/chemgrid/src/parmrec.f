c---------------------------------------------------------------------
c
c       Copyright (C) 1991 Regents of the University of California
c                         All Rights Reserved.
c
      subroutine parmrec(recfil, table, vdwfil, unitno)
c
c     --called from CHEMGRID
c          Parmrec reads charges and VDW parameters for receptor
c     atom types from the appropriate files, indexes them via a hash
c     table, and then associates them with the atoms in a given
c     pdb-format receptor file.
c     Much of this code, namely the hashing and lookup routines,
c     has been adapted from the DelPhi code (program qdiffx and
c     subroutines) of Honig et al., version 3.0.
c
c     ECMeng     January 1991
c  4/93 ECM  altered to report residues with nonzero charge to OUTCHEM
c
c  recfil--name of receptor pdb file
c  table--name of the table to be referenced for receptor atom
c    parameters 
c  vdwfil--name of file containing van der Waals parameters
c  unitno--logical unit number to write parameterization information
c    and warnings to
c---------------------------------------------------------------------
      include 'parmrec.h'
c
      character*80 recfil, table, vdwfil
      integer unitno
      integer i, natm, n
c
      character*3 presid
      integer resn, prevn
      real rescrg
c
      nptyp=0
      do 10 i=1,maxtyp
        inum(i)=0
        ilink(i)=0
   10 continue
c
c  --read receptor atom parameter file, index entries via a hash table
c
      open (unit=11, file=table, status='old')
c
  100 read (11, 1000, end=190) line
 1000 format (A80)
      if (line(1:1) .eq. '!') go to 100
      nptyp=nptyp + 1
      if (nptyp .gt. maxtyp) then
        write (6, *)
     &  'maximum number of atom types exceeded'
        write (6, *) 'increase parameter maxtyp'
        stop
      endif
      read (line, 1001) atm(nptyp), res(nptyp), resnum(nptyp),
     &chain(nptyp), crg(nptyp), vdwtyp(nptyp)
 1001 format (A4, 3x, A3, A4, A1, F8.3, 1x, I2)
c
      call enter(atm(nptyp), res(nptyp), resnum(nptyp),
     &chain(nptyp), nptyp)
c
      go to 100
  190 continue
      close (11)
c
c  --read vdw parameter file
c
      open (unit=11, file=vdwfil, status='old')
c
      nvtyp=0
  200 read (11, 1000, end=290) line
      if (line(1:1) .eq. '!') go to 200
      nvtyp=nvtyp + 1
      if (nvtyp .gt. maxtyv) then
        write (6, *) 'maximum number of vdw types exceeded'
        write (6, *) 'increase parameter maxtyv'
        stop
      endif
      read (line, 1002) sra(nvtyp), srb(nvtyp)
 1002 format (10x, F8.2, 5x, F8.2)
      go to 200
  290 continue
      close (11)
c
c  --read receptor pdb file, associate atoms with parameters, write
c    parameters and coordinates out to another file (PDBPARM)
c
      natm=0
      crgtot=0.0
      rescrg=0.0
c
      open (unit=11, file=recfil, status='old')
      open (unit=12, file='PDBPARM', status='new')
      open (unit=13, file='OUTPARM', status='new')
c
   20 read (11, '(A80)', end=990) line
      if (line(1:4) .ne. 'ATOM' .and. line(1:4) .ne. 'HETA') go to 20
      natm=natm + 1
      if (natm .gt. maxatm) then
        write (6, *) 'maximum number of receptor atoms exceeded'
        write (6, *) 'increase parameter maxatm'
        stop
      endif
c
      if (resid.ne.'   ') then
        presid=resid
      else
        presid=sresid
      endif
c
      atom=line(13:16)
      resid=line(18:20)
      chn=line(22:22)
      resno=line(23:26)

      read (resno, *) resn
c
      call find(atom, resid, resno, chn, found, n)
      if (.not. found) then
        schn=chn
        chn=' '
        call find(atom, resid, resno, chn, found, n)
        if (.not. found) then
          chn=schn
          sresno=resno
          resno='    '
          call find(atom, resid, resno, chn, found, n)
          if (.not. found) then
            schn=chn
            chn=' '
            call find(atom, resid, resno, chn, found, n)
            if (.not. found) then
              chn=schn
              resno=sresno
              sresid=resid
              resid='   '
              call find(atom, resid, resno, chn, found, n)
              if (.not. found) then
                schn=chn
                chn=' '
                call find(atom, resid, resno, chn, found, n)
                if (.not. found) then
                  chn=schn
                  sresno=resno
                  resno='    '
                  call find(atom, resid, resno, chn, found, n)
                  if (.not. found) then
                    schn=chn
                    chn=' '
                    call find(atom, resid, resno, chn, found, n)
                    if (.not. found) then
                      write (13, *) 'WARNING--parameters not found for'
                      write (13, *) line(1:27)
                      write (13, '(A18, A21)') 'sqrt(A), sqrt(B), ',
     &                'and charge set to 0.0'
                      write (12, 2000) natm, 0, 0.0, 0.0, 0.0,
     &                line(31:54)
                      go to 20
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif
      endif
      write (12, 2000) natm, vdwtyp(n), sra(vdwtyp(n)), srb(vdwtyp(n)),
     &crg(n), line(31:54)
 2000 format (2I5, 2(1x, F8.2), 1x, F8.3, 1x, A24)
c
      if (natm.eq.1) prevn=resn
      if (resn.ne.prevn) then
        if (abs(rescrg).gt.0.0001) then
          write (13, 2001) ' CHARGED RESIDUE ', presid, prevn, rescrg
 2001 format (A17, A3, I5, F8.3)
        endif
        rescrg = crg(n)
        prevn = resn
      else
        rescrg = rescrg + crg(n)
      endif
c
      crgtot=crgtot + crg(n)
      go to 20
  990 continue
      if (abs(rescrg).gt.0.0001) then
        if (resid.eq.'   ') resid = sresid
        write (13, 2001) ' CHARGED RESIDUE ', resid, resn, rescrg
      endif
      close (11)
      close (12)
      write (13, *) ' '
      write (13, '(A15, F8.3)') 'Total charge = ', crgtot
      close (13)
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine enter(atom, resid, resno, chn, nent)
c
      include 'parmrec.h'
c
c  --enter receptor atom type entries into hash table according to 
c    entry number (sequential number of occurrence within the parameter
c    table)
c
      integer n, new, nent
c
      integer ihash
c
c
c  --get hash number using function ihash
c
      n=ihash(atom, resid, resno, chn)
      if (inum(n) .ne. 0) then
c
c  --slot filled; keep going along linked numbers until zero found
c
  100   continue
        if (ilink(n) .eq. 0) go to 200
        n=ilink(n)
        go to 100
  200   continue
c
c  --find an empty slot and fill it, leaving a trail in ilink()
c
        do 300 new=1,maxtyp
        if (inum(new) .eq. 0) go to 400
  300   continue
  400   continue
        ilink(n)=new
        n=new
      endif
      inum(n)=nent
      ilink(n)=0
      return
      end
c
c--------------------------------------------------------------------
c
      integer function ihash(atxt,rtxt,ntxt,ctxt)
c       
c  --produce a hash number for an atom, using atom name, residue name,
c    residue number, and chain indicator
c       
      include 'parmrec.h'

      character*4 atxt
      character*3 rtxt
      character*4 ntxt
      character*1 ctxt
      character*38 string
      integer n, i, j
      data string /'* 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      n = 1
      do 100 i = 1,3
        j = index(string,rtxt(i:i))
        n = 5*n + j
 100    continue
      do 101 i = 1,4
        j = index(string,atxt(i:i))
        n = 5*n + j
 101    continue
      do 102 i = 1,4
        j = index(string,ntxt(i:i))
        n = 5*n + j
 102    continue
      do 103 i = 1,1
        j = index(string,ctxt(i:i))
        n = 5*n + j
 103    continue
        n = iabs(n)
      ihash = mod(n,maxtyp) + 1
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine find(atom, resid, resno, chn, found, n)
c
c  --use the hash number of a receptor atom to find the appropriate
c    parameters, following links when necessary; check explicitly for
c    a match
c
      include 'parmrec.h'
c
      integer n
      integer ihash
c
      n=ihash(atom, resid, resno, chn)
      found=.false.
  100 continue
      if (inum(n) .eq. 0) then
        found=.false.
        return
      endif
      if ((resid .eq. res(inum(n))) .and. (atom .eq. atm(inum(n)))
     &.and. (resno .eq. resnum(inum(n))) .and. (chn .eq.
     &chain(inum(n)))) then
        n=inum(n)
        found=.true.
        return
      else
        if (ilink(n) .ne. 0) then
          n=ilink(n)
        else
          found=.false.
          return
        endif
      endif
      go to 100
      end
c
c-----------------------------------------------------------------------
