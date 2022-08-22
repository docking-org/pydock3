c----------------------------------------------------------------------
      subroutine grdout(grdfil, unitno, npts, grddiv, grdpts, offset)
c
c  --called from CHEMGRID
c  --writes out grids; makes a formatted "bump" file and unformatted
c    van der Waals and electrostatics files
c                                                  ECMeng    4/91
c----------------------------------------------------------------------
      include 'chemgrid.h'
c
      character*80 grdfil
      integer i, namend, unitno
c
      namend=80
      do 100 i=2,80
	if (grdfil(i:i) .eq. ' ') then
	  namend=i-1
	  go to 105
        endif
  100 continue
  105 continue
c
    1 format (A17)
    2 format (4F8.3, 3I4)
    3 format (80A1)
      open (unit=unitno, file=grdfil(1:namend)//'.bmp', status='new')
      write (unitno, 1) 'bump map         '
      write (unitno, 2) grddiv, (offset(i), i=1,3), (grdpts(i), i=1,3)
      write (unitno, 3) (bump(i), i=1, npts)
      close (unitno)
      open (unit=unitno, file=grdfil(1:namend)//'.vdw', status='new',
     &form='unformatted')
      write (unitno) (aval(i), i=1, npts)
      write (unitno) (bval(i), i=1, npts)
      close (unitno)
      open (unit=unitno, file=grdfil(1:namend)//'.esp', status='new',
     &form='unformatted')
      write (unitno) (esval(i), i=1, npts)
      close (unitno)
c
      return
      end
c----------------------------------------------------------------------
