	function cputime(time0)
c
	dimension tarray(2)
c should return elapsed system time
c since time0 as real number of seconds
c
c	cputime = your_system_timer
c eg cray:
c	cputime = rtc()/1.e9 - time0
c vax/convex:
	cputime = secnds(time0)
c	call etime(tarray)
c	cputime = tarray(1) - time0
	end
