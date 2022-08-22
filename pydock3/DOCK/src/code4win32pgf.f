c
c  This is code specific for the SGI
c  SGI Getarg
      subroutine get_indock(pos, argument)
      integer pos
      character*80 argument
      call getarg(pos,argument)
      return
      end

c  SGI Time
c (note:  on some machines, etime is double precision, so the following
c should be declared accordingly...)
      subroutine get_time(time_check)
      dimension tarray(2)
      real time_check, etime, tarray
      time_check=etime(tarray)
      return
      end

c  SGI Flush
      subroutine doflush(unit_num)
      integer unit_num
      call flush(unit_num)
      return
      end

c SGI Binary open
      subroutine binopen(unit_num,filename)
      integer unit_num
      character*80 filename
      open(unit_num, file=filename, status='old', form='unformatted') 
      return
      end

c run information
      subroutine run_info(the_date, the_time, the_host)

      integer hostnm, i
      character*8 the_time
      character*(*) the_date
      character*80 the_host
      character*80 nm

      call date(the_date)
      call time(the_time)
c Linux hostname
      nm = 'x'
      i = 'the hos'
      the_host = nm
      return
      end
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
