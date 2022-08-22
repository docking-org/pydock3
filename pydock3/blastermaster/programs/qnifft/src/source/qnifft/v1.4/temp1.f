	  epsin = epsin*epkT
	  epsout = epsout*epkT
        write(16,*)'output from QNIFFT   '
        write(16,*)'grid size,percent fill:',igrid,perfil
        write(16,*)'inner,outer dielectric:',epsin,epsout
        write(16,*)'ionic strength (M):',rionst
        write(16,*)'ion excl., probe radius:',exrad,radprb
        write(16,*)'level 0, newton iterations:',nit0,nitmx
        write(16,*)'boundary condition:',ibctyp
        write(16,*)'    '
        write(16,*)'title: '
        write(16,*)toplbl
        write(16,233)
        write(16,234)
233     format('         coordinates          charge   potential  field (kT/e/Ang.)')
234     format('     x         y       z       (e)     (kT/e)    Ex        Ey       Ez    Atom Res C   #')

c               123456789 123456789 123456789 123456789 123456789 123456789 123456789
	  etot = 0.
	  do i = 1,natom
	    do k = 1,3
	      xo(k) = atcrd(k,i)
	    end do
	    call ctog(xo,xn)
c
c calculate potential, field and energy, and output to file
c
          call phintp(xn,phiv)
	    xn(1) = xn(1) + 1.
          call phintp(xn,fxu)
	    xn(1) = xn(1) - 2.
          call phintp(xn,fxl)
	    xn(1) = xn(1) + 1.
	    xn(2) = xn(2) + 1.
          call phintp(xn,fyu)
	    xn(2) = xn(2) - 2.
          call phintp(xn,fyl)
	    xn(2) = xn(2) + 1.
	    xn(3) = xn(3) + 1.
          call phintp(xn,fzu)
	    xn(3) = xn(3) - 2.
          call phintp(xn,fzl)
	    xn(3) = xn(3) + 1.
          qphiv = chrgv*phiv
	    etot = etot + qphiv
	    fx = -(fxu-fxl)*scale/2.
	    fy = -(fyu-fyl)*scale/2.
	    fz = -(fzu-fzl)*scale/2.
	    write(16,230)xo,chrgv,phiv,fx,fy,fz,line(12:26)
230       format(8f11.4,x,a15)
	  end do
	  close(15)
	  write(6,*)'number of target atom coordinates read  : ',natom
	  write(6,*)'total number of target charged atoms    : ',ntqass
	  write(6,*)'net assigned target charge              : ',qtnet
	  write(6,*)'assigned positive charge                : ',qtplus
	  write(6,*)'assigned negative charge                : ',qtmin
	  write(6,*)'   '
	  if(isitecrg) then
	    print *,'Interaction between source and target charges (kT)',etot
	    print *,'(Includes Coulombic+Solvation energy)'
	    write(16,*)'Interaction between source and target charges (kT)',etot
	    write(16,*)'(Includes Coulombic+Solvation energy)'
	  else
	    etot = etot/2.
	    print *,'Interaction of source charges with themselves (kT)',etot
	    print *,'(Includes Grid+Coulombic+Solvation energy)'
	    write(16,*)'Interaction of source charges with themselves (kT)',etot
	    write(16,*)'(Includes Grid+Coulombic+Solvation energy)'
	  end if
	  close(16)
