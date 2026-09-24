c----------------------------------------------------------------
c                  distribution version
c     header file for subroutine parmrec          ECMeng   4/91
c----------------------------------------------------------------
      integer maxtyp, nptyp
      parameter (maxtyp=5000)
c  maxtyp--maximum number of entries in 'prot.table' or 'na.table'
c  nptyp--number of entries in 'prot.table' or 'na.table' so far
      integer inum(maxtyp), ilink(maxtyp)
c  inum()--id numbers in hash table
c  ilink()--links for hash table
      character*1 chain(maxtyp), chn, schn
      character*3 res(maxtyp), resid, sresid
      character*4 atm(maxtyp), resnum(maxtyp), atom, resno, sresno
      real crg(maxtyp)
      integer vdwtyp(maxtyp)
c  vdwtyp()--integer vdw type indicators
      integer maxtyv
      parameter (maxtyv=50)
c  maxtyv--maximum number of entries in 'vdw.parms'
      integer nvtyp
c  nvtyp--number of entries in 'vdw.parms' so far
      real sra(maxtyv), srb(maxtyv)
c  sra(), srb()--vdw parameters, sqrt of A and sqrt of B
      integer maxatm
      parameter (maxatm=100000)
c  maxatm--maximum number of receptor atoms
      logical found
      character*80 line
      real crgtot
c
      common
     &/link/ inum, ilink
     &/name/ atm, res, resnum, chain
     &/value/ crg, vdwtyp, sra, srb
c----------------------------------------------------------------
