c this header file just has format lines for use when reading mol2 formate 
c Trent Balius

c XXX: GFortran doesn't know that x means just one space. It requires you 
c provide an explicit number of leading spaces (e.g. 1x). Ran %s/,x/,1x/g

      module mol2formats

      implicit none


c@<TRIPOS>MOLECULE
c           331902        98
c   15    16     0     0     0
c
c
c
c@<TRIPOS>ATOM
c      1 N1         46.9071    43.2086    26.9254 N.pl3      1          -0.6036
c      2 C1         48.2161    43.4687    27.1887 C.2        1           0.2623
ccc    *    *  *          *          *          **    *       *                *
ccc    7    5  3          11         11         11    5       8                17
ccc                                              1
c@<TRIPOS>BOND
c     1    1    9 1
c     5    2   11 1
c     7    4    9 ar
c     8    4    5 ar
c      *    *    * *
c      7    5    5 2
      !character (len=*), parameter :: MOL2MOLE1 = '(A18,i10)' ! name, number
      character (len=*), parameter :: MOL2MOLE1 = '(1x,a16,1x,a9)' ! name, number
      character (len=*), parameter :: MOL2MOLE2 = '(i5,i6,i6,i6,i6)' ! number of atoms, number of bonds, more info. 
      character (len=*), parameter :: MOL2ATOM = 
c    &    '(i7,a5,3x,3f11.5,1x,a5,a8,f16.5)' ! atom number, atom name, x,y,z coord, atom type, residue name, charge. 
     &    '(i7,a5,3x,3f11.5,1x,a5,a8,6x,f10.5)' ! atom number, atom name, x,y,z coord, atom type, resnum, resname name, charge. 
c    &    '(i7,a5,3x,3f11.5,1x,a5,a8,a8,f8.5)' ! atom number, atom name, x,y,z coord, atom type, residue name, charge. 
      character (len=*), parameter :: MOL2BOND = 
     &    '(i6,2i5,1x,a2)' ! 3 int (bond number,atom orgin, atom target) and 1 string ( bond type)

      end module mol2formats
