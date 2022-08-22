c this file holds the things necessary to write out a bunch of specific
c ligand poses. this is in combination with everything in molecule and hierarchy
c that is read in/set for each ligand. here we only need to save
c the coml, comr, rot matrix and set for each pose we want to output. 
c also save the overall energies of each type (v,e,p,a)
      module mol2sav

      implicit none

c GFortran doesn't know that the X formatting character means a single space,
c it needs an explicit number of spaces (e.g. 1x instead of x)
c ran %s/,x/,1x/g To fix

c these are format lines for mol2 files
! 1020 format('##########',1x,a20,':',1x,a80) !dock6 mol2 file/ chimera viewdock
      character (len=*), parameter :: MOL2A = 
     &    "('##########',1x,a20,':',1x,a80)" !dock6 mol2 file/ chimera viewdock
! 1030 format('##########',1x,a20,':',1x,a20) !dock6 mol2 file/ chimera viewdock
      character (len=*), parameter :: MOL2B = 
     &    "('##########',1x,a20,':',1x,a20)" !dock6 mol2 file/ chimera viewdock
! 1040 format('##########',1x,a20,':',1x,i20) !integer
      character (len=*), parameter :: MOL2I = 
     &    "('##########',1x,a20,':',1x,i20)" !integer
! 1050 format('##########',1x,a20,':',1x,f20.6) !compatibility
      character (len=*), parameter :: MOL2F = 
     &    "('##########',1x,a20,':',1x,f20.6)" !floating point?
!lines for writing out per atom scores
      character (len=*), parameter :: MOL2ASD = !mol2 atom score data
     &    "('#',1x,i3,1x,7(f7.2,1x),f10.2,1x,2(f7.2,1x))" 
! 1100 format(a17) !molecule line
      character (len=*), parameter :: MOL2MOL = "(a17)" !molecule line
! 1150 format(a13) !bond or atom line
      character (len=*), parameter :: MOL2AB = "(a13)" !bond or atom line
! 1200 format(a80) !generic
      character (len=*), parameter :: MOL280 = "(a80)" !generic
! 1250 format(a17,1x,a9) !codes output
      character (len=*), parameter :: MOL2C = "(a17,1x,a9)" !codes output
! 1300 format(i5,1x,i5,1x,i5,1x,i5,1x,i5) !atom, bond, 0,0,0
      character (len=*), parameter :: MOL2I5 = 
     &    "(i5,1x,i5,1x,i5,1x,i5,1x,i5)" !atom, bond, 0,0,0
! 1400 format(1x,i6,1x,a4,4x,f10.4,1x,f10.4,1x,f10.4,
!     &        1x,a5,1x,i6,10x,f7.4) !atom
      character (len=*), parameter :: MOL2ATOM =
     &    "(1x,i6,1x,a4,4x,f10.4,1x,f10.4,1x,f10.4,
     &    1x,a5,1x,i6,2x,a4,2x,f7.4)" !atom  ! modifed to wite out residue names
!    &    1x,a5,1x,i6,10x,f7.4)" !atom
! 1500 format(i6,i5,i5,1x,a2) !bond
      character (len=*), parameter :: MOL2BOND = 
     &    "(i6,i5,i5,1x,a2)" !bond
      character (len=*), parameter :: MOL2X3 = 
     &    "('##########',1x,a20,':',1x,f10.4,1x,f10.4,1x,f10.4)" !3 coordinates

      end module mol2sav
