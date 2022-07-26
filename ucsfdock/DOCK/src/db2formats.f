c this header file just has format lines for use when reading the new db2 
c ligand database files. 


c XXX: GFortran doesn't know that x means just one space. It requires you 
c provide an explicit number of leading spaces (e.g. 1x). Ran %s/,x/,1x/g

      module db2formats

      implicit none

!T ## namexxxx (implicitly assumed to be the standard 7)
      character (len=*), parameter :: DB2NAME = '(2x,i2,1x,a8)' !1000
!M zincname protname #atoms #bonds #xyz #confs #sets #rigid #maxmlines #clusters
      character (len=*), parameter :: DB2M1 =
     &    '(2x,a16,1x,a9,1x,i3,1x,i3,1x,i6,1x,i6,1x,i6,1x,i6,1x,
     &    i6,1x,i6)' !2000
!M charge polar_solv apolar_solv total_solv surface_area
      character (len=*), parameter :: DB2M2 = 
     &    '(2x,f9.4,1x,f10.3,1x,f10.3,1x,f10.3,1x,f9.3)' !2100
!M smiles/longname/arbitrary
      character (len=*), parameter :: DB2M3 = '(2x,a78)' !2200
!A stuff about each atom, 1 per line 
      character (len=*), parameter :: DB2ATOM = 
     &    '(2x,i3,1x,a4,1x,a5,1x,i2,1x,i2,1x,f9.4,1x,f10.3,1x,
     &    f10.3,1x,f10.3,1x,f9.3)' !3000
!B stuff about each bond, 1 per line
      character (len=*), parameter :: DB2BOND = 
     &    '(2x,i3,1x,i3,1x,i3,1x,a2)' !4000
!X coordnumx atomnum confnum x y z 
      character (len=*), parameter :: DB2COORD = 
     &    '(2x,i9,1x,i3,1x,i6,1x,f9.4,1x,f9.4,1x,f9.4)' !5000
!R rigidnum color x y z
      character (len=*), parameter :: DB2RIGID = 
     &    '(2x,i6,1x,i2,1x,f9.4,1x,f9.4,1x,f9.4)' !6000
!C confnum coordstart coordend
      character (len=*), parameter :: DB2CONF = '(2x,i6,1x,i9,1x,i9)' !7000
!S setnum #lines #confs_total broken hydrogens omega_energy
!      character (len=*), parameter :: DB2SET1 = 
!     &    '(2x,i6,1x,i6,1x,i3,1x,i1,1x,i1,1x,f11.3)' !8000
!S setnum #lines #confs_total broken hydrogens totalStrain maxStrain
      character (len=*), parameter :: DB2SET1 = 
     &    '(2x,i6,1x,i6,1x,i3,1x,i1,1x,i1,1x,f11.3)' !8000
!S setnum #lines #confs_total broken hydrogens totalStrain maxStrain
      character (len=*), parameter :: DB2SET1_S = 
     &    '(2x,i6,1x,i6,1x,i3,1x,i1,1x,i1,1x,f11.3,f11.3)' !8000
!S setnum linenum #confs confs [until full column]
      character (len=*), parameter :: DB2SET2 = 
     &    '(2x,i6,1x,i6,1x,i1,1x,i6,1x,i6,1x,i6,1x,i6,
     &    1x,i6,1x,i6,1x,i6,1x,i6)' !8100
!D CLUSID STASET ENDSET ADD(ittional matching spheres count) MST(art) MEN(d)
      character (len=*), parameter :: DB2CLUSTER = 
     &    '(2x,i6,1x,i6,1x,i6,1x,i3,1x,i3,1x,i3)' !9000
!D NUM CO x y z
!reuse DB2RIGID
!E 
!E does not get a format line

      end module db2formats
