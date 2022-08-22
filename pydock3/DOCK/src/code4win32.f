c  This is code specific for the PC
c  PC Getarg
      subroutine get_indock(pos, argument)
      use dflib
      integer pos
      character(80) argument
c      call getarg(pos,argument,1)
      argument = 'INDOCK'
      return
      end

c  PC Time
      subroutine get_time(time_check)
      real time_check
      call cpu_time(time_check)
      return
      end

c  PC Flush
      subroutine doflush(unit_num)
      use dfport
      integer unit_num
      call flush(unit_num)
      return
      end

c PC Binary open of SGI binary file
      subroutine binopen(unit_num,filename)
      integer unit_num
      character*60 filename
      open(unit_num, file=filename, status='old', 
     &convert='big_endian',form='unformatted') 
      return
      end

c run information
      subroutine run_info(the_date, the_time, the_host)

      integer hostnm, i
      character*8 the_time
      character*(*) the_date
      character*80 the_host

      call date(the_date)
      call time(the_time)
c WinNT hostname
      
      the_host = 'Win 32-bit'
      return
      end
