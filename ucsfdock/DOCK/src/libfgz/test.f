      program test
      character(1000000) arr
      integer istat,i,lnblnk
      integer*8 handle,handle2

      print*,'testing binary file i/o'
      call gzopen(handle,'r','test_bin.gz',istat)
      print*,'file open'
      call gzbread(handle,arr,1000000,istat)
      print*,'file read'
      call gzclose(handle,istat)
      print*,'file closed'

      call gzopen(handle,'w','testout_bin.gz',istat)
      print*,'file open'
      call gzbwrite(handle,arr,1000000,istat)
      print*,'file written'
      call gzclose(handle,istat)
      print*,'file closed'

      print*,'testing ascii file i/o'
      call gzopen(handle,'r','test_ascii.gz',istat)
      call gzopen(handle2,'w','testout_ascii.gz',istat)
      print*,'opened simulataneously input and output files'

      do
        do i=1,lnblnk(arr)
          arr(i:i)=' '
        enddo
        call gzread(handle,arr,1000000,istat)
        print*,istat,lnblnk(arr),arr(1:istat)
        if (istat.le.0) exit
        call gzwrite(handle2,arr,istat)
      enddo
      call gzclose(handle,istat)
      call gzclose(handle2,istat)

      end
