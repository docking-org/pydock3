c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
c
      subroutine sim_simplex(p, y, delta_energy, iter, best_energy,
     &           confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)


c      Multidimensional minimization of the function funk(x) where x is an
c      nvar-dimensional vector, by the downhill simplex method of Nelder 
c      and Mead.  In this case, funk is the energy function.
c      Input is a matrix p whose nvar+1 rows are nvar-dimensional
c      vectors which are the vertices of the starting simplex.  Also input 
c      is the vector y of length nvar+1, whose      components must be 
c      pre-initialized to the values of funk evaluated at the nvar+1 
c      vertices (rows) of p.
c      On output, p and y will have been reset to nvar+1 new points all 
c      within sim_cnvrgE of a minimum function value; iter gives the number
c      of iterations taken, and deltaE gives the energy convergence 
c      reached at completion.
c
c               Taken from Numerical Recipes:  The Art of Scientific
c               Computing, by Press, Flannery, Teukolsky, and Vetterling.
c               1986 by Cambridge University Press, FORTRAN edition
c               pp. 292-293.
c
c      References and further reading:
c            Nelder, J.A., and Mead, R. 1965, Computer Journal, vol. 7,
c                  p. 308.
c            Yarbro, L.A., and Deming, S.N. 1974, Analytica Chim. Acta,
c                  vol. 73, p. 391.
c            Jacoby, S.L.S., Kowalik, J.S., and Pizzo, J.T. 1972, 
c                  Iterative Methods for Nonlinear Optimization Problems
c                  (Englewood Cliffs, N.J.:  Prentice-Hall).
c
      use db2type
      use optionstype
      use gridstype
      use phimaptype
      use gisttype
      use solvmaptype
      use vdwmaptype
      use phiscore
      use gistscore
      use rec_des_score
      use chemscore
      use solvation_score
      use matchtype
      use ligscoretype
      use atomscoretype
      use status

      implicit none

      type(options), intent(in) :: options0 !useful options
      type(db2), intent(inout) :: db2lig !ligand information
      type(matcht), intent(inout) :: match
      integer, intent(in):: matchnum !for getting rotation matrices
      integer, intent(in) :: MAXOR ! max orientations
      character (len=255) :: reccode !flexible receptor code for this pose
      type(flexgrids), intent(in) :: grids0 !precomp values
      type(vdwmap), intent(in) :: vdwmap0
      type(phimap), intent(in) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(in) :: gist0, rec_des_r0
      type(solvmap), intent(in) :: solvmap0
      integer, intent(in) :: maxtyv !how many vdw atom types there are
      integer, intent(in) :: confset !the confset that the caller would like
                                  !scored
      real, dimension(maxtyv), intent(in) :: sra, srb !square roots of
                                                !vdw parameters per vdwtype
c return variables

      integer :: atomstart, atomend !temp/loop variables
      integer, intent(inout) :: min_status

      real, intent(inout) ::  vscore, escore, gscore,
     &           ascore, pscore,
     &           rscore, rdscore, hscore, setscore
c     &           rscore, dscore, hscore, setscore
      integer :: mpts, i, j
      integer :: ilo,ihi,inhi
      logical :: restart

c      nvar: number of variables to minimize
      !integer, parameter :: nvar !=6 ! = 3 translation + 3 rotation
      !real :: sim_cnvrge=0.20
      !integer :: sim_itmax=500
c     variables for experimental jump-out routine
      real, dimension(options0%nvar+1, options0%nvar),
     &      intent(inout):: p ! simplx
      real, dimension(options0%nvar+1), intent(inout) :: y ! simplx_energy
      integer, intent(inout) :: iter
      real, intent(inout) :: delta_energy
      real, dimension(options0%nvar) :: pr
      real, dimension(options0%nvar) :: prr, pbar
      real :: ypr, yprr, best_energy
      real :: sim_energy   ! a function
      real, dimension(3), intent(in) :: allminlimits, allmaxlimits

c      Note that mp is the physical dimension corresponding to the logical
c      dimension mpts; np to nvar.
      mpts = options0%nvar + 1
      restart = (iter .ne. 0)

c      First we must determine which point is the highest (worst), next-
c      highest, and lowest (best).
 2000 ilo = 1
      if (y(1) .gt. y(2)) then
          ihi = 1
          inhi = 2
      else
          ihi = 2
          inhi = 1
      endif

c      loop over points in the simplex
      do i = 1, mpts
          if (y(i) .lt. y(ilo)) then
             ilo = i
          endif
          if (y(i) .gt. y(ihi)) then
            inhi = ihi
            ihi = i
          else if (y(i) .gt. y(inhi)) then
            if (i .ne. ihi) inhi = i
          endif
      enddo

c      compute the fractional range from highest to lowest and return if
c      satisfactory
      delta_energy = abs(y(ihi) - y(ilo))
      if (delta_energy .lt. options0%sim_cnvrge) return
      if (iter .eq. options0%sim_itmax) then
          return
      endif
      if (restart) then
          restart = .false.
      endif


c      begin a new iteration.  compute the vector average of all points
c      except the highest, i.e. the center of the "face" of the simplex
c      across from the high point.  we will subsequently explore the
c      ray from the high point through that center.
      iter = iter + 1
      do j = 1, options0%nvar
          pbar(j) = 0.
      enddo
      do i = 1, mpts
          if (i .ne. ihi) then
              do j = 1, options0%nvar
                  pbar(j) = pbar(j) + p(i,j)
              enddo
          endif
      enddo

c      extrapolate by a factor alpha through the face, i.e. reflect the
c      simplex from the high point.
      do j = 1, options0%nvar
          pbar(j) = pbar(j)/options0%nvar
          pr(j) = (1.+options0%alpha)*pbar(j) - options0%alpha*p(ihi,j)
      enddo

c      evaluate the function at the reflected point.
      ypr = sim_energy(pr, best_energy, confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)


      if (ypr .le. y(ilo)) then

c          gives a result better than the best point, so try an 
c          additional extrapolation by a factor gamma.
          do j = 1, options0%nvar
            prr(j) = options0%gamma*pr(j) + (1.-options0%gamma)*pbar(j)
          enddo

c          check out the function there.
          yprr = sim_energy(prr, best_energy, confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)


          if (yprr .lt. y(ilo)) then

c            the additional extrapolation succeeded, and replaces 
c            the high point.
            do j = 1, options0%nvar
                p(ihi,j) = prr(j)
            enddo
            y(ihi) = yprr

          else

c            the additional extrapolation failed, but we can still use
c            the reflected point
            do j = 1, options0%nvar
                p(ihi,j) = pr(j)
            enddo
            y(ihi) = ypr

          endif

      else if (ypr .ge. y(inhi)) then

c          the refelcted point is worse than the second-highest.  
c          if it's better than the highest, then replace the highest.
          if (ypr .lt. y(ihi)) then
            do j = 1, options0%nvar
                p(ihi,j) = pr(j)
            enddo
            y(ihi) = ypr
          endif

c          but look for an intermediate lower point, in other words,
c          perform a contraction of the simplex along one dimension, then
c          evaluate the function.  
          do j = 1, options0%nvar
            prr(j) = options0%beta*p(ihi,j) + (1.-options0%beta)*pbar(j)
          enddo

          yprr = sim_energy(prr, best_energy, confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)


          if (yprr .lt. y(ihi)) then

c            contraction gives an improvement, so accept it.
            do j = 1, options0%nvar
                p(ihi,j) = prr(j)
            enddo
            y(ihi) = yprr

          else

c            can't seem to get rid of that high point.  better contract
c            around the lowest (best) point.
            do i = 1, mpts
                if (i .ne. ilo) then
                  do j = 1, options0%nvar
                      pr(j) = 0.5 * (p(i,j) + p(ilo,j))
                      p(i,j) = pr(j)
                  enddo

                 y(i) = sim_energy(pr, best_energy, confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)


                endif
            enddo

          endif

      else

c          we arrive here if the original reflection gives a middling point.
c          replace the old high point and continue.
          do j = 1, options0%nvar
            p(ihi,j) = pr(j)
          enddo
          y(ihi) = ypr

      endif

c      test for doneness and next iteration
      goto 2000


      end subroutine sim_simplex
