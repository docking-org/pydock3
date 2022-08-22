c
c  This is code specific for the SGI
c  SGI Getarg
      subroutine get_indock(pos, argument)
        integer pos
        character (len=80) argument
        call getarg(pos,argument)
        return
      end

c  SGI Time
c (note:  on some machines, etime is double precision, so the following
c should be declared accordingly...)
      subroutine get_time(time_check)
        real time_check, etime, tarray(2)
        time_check = etime(tarray)
        return
      end

c  SGI Flush
      subroutine doflush(unit_num)
        integer , intent(in) :: unit_num
        call flush(unit_num)
        return
      end

c SGI Binary open
      subroutine binopen(funit, filename)
        integer funit
        character (len=255) filename
        open(funit, file=filename, status='old', form='unformatted', 
     &      action='read')
        return
      end

c run information
      subroutine run_info(the_date, the_time, the_host)
        integer hostnm, i
        character (len=8) the_time
        character (len=8) the_date
        character (len=80) the_host
        character (len=80) nm

        ! XXX: Using Fortran95 method of getting date and time
        call date_and_time(DATE=the_date,TIME=the_time)

c       Linux hostname
        nm = 'x'
        i = hostnm(nm)
        the_host = nm
        return
      end
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
