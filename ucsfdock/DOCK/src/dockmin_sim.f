c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c
      subroutine dockmin_sim(best_energy, confset, options0,
     &           db2lig, match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits, total_iterations)



C////////////////////////////////////////////////////////////////////////////
C PROGRAMMER:      Daniel A. Gschwend
C REV.DATE:      20-APR-94
C HISTORY:      30-JUL-93 v1.00 
C INPUT:      DOCK3 grid information, ligand orientation to be score
C OUTPUT:      Ligand coordinates minimized to local energy minimum on grid
C PURPOSE:      Provides minimization of orientations to be score during
C            run-time operation of DOCK3 so that a) more representative
C            energies may be obtained and b) sampling space may hence be
C            reduced.
C            Rigid-body minimization of interaction energy
C            (6-12 Lennard-Jones + Electrostatics) of two molecules
C            using the simplex minimizer of Nelder & Mead (taken from
C            Numerical Recipes - see subroutine simplex for complete refs).
C            Portions adapted from J.M.Blaney's RGDMIN.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


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


      type(options), intent(inout) :: options0 !useful options
      type(db2), intent(inout) :: db2lig !ligand information
      type(matcht), intent(inout) :: match
      integer, intent(in):: matchnum !for getting rotation matrices
      integer, intent(in) :: MAXOR ! max orientations
      character (len=255), intent(in) :: reccode !flexible receptor code for this pose
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
      integer, intent(inout) :: min_status
      real, intent(inout) ::  vscore, escore, gscore,
     &           ascore, pscore,
     &           rscore, rdscore, hscore, setscore
c     &           rscore, dscore, hscore, setscore

      integer :: atomstart, atomend !temp/loop variables
      real, intent(inout) :: best_energy
      integer, intent(inout) :: total_iterations

      integer :: i, j
      integer :: iter, lowest, irest

c     iter      = number of simplex iterations performed
c     lowest      = index of lowest energy vertex on simplex
c     irest      = number of simplex restarts performed

      real   :: old_energy, new_energy
c     new_energy      = energy at completion of most recent simplex
c     old_energy      = energy at completion of previous simplex

      real, dimension(options0%nvar+1, options0%nvar)  :: simplx ! 3 trans + 3 rot variables
c              at each of nvar+1 vertices  of the simplex
      real, dimension(options0%nvar+1) :: simplx_energy ! energies at
      !          each of the vertices of the simplex
      real, dimension(options0%nvar) :: tmp
      real :: random, ran, min_energy, delta_energy
c      real :: rand, random_num
      real :: random_num, rand_val
      real :: sim_energy !  a function to compute intermolecular interaction
c                    energy, either by continuum or by grid (chemscore)

      real, dimension(3), intent(in) :: allminlimits, allmaxlimits

c     tmp      = temporary array for passing to energy routine
c     sim_trnstep = maximum translation step size
c     sim_rotstep = maximum rotation step size
c     random  = a random number from -1.0 to 1.0
c     ran      = a random number generator function
c     min_energy      = energy of lowest energy vertex
c     bigE      = large energy to initialize finding min_energy
c     deltaE      = energy convergence returned from simplex

      real :: time2, time_vimelapsed
      !integer :: seeds(2)
      !integer, allocatable :: seeds(:) ! moved to options so that it is
      !remembered
      integer :: n, countn 

c timecheck

!      call get_time(time2)
!      time_elapsed = time2 - time1
!      if (time_elapsed .gt. timeout) then
!        time1 = -1E0
!        return
!      endif
c endtimecheck

c     calculate center of geometry
c     to perform a->c transforms the center of mass must be the same
c     in all structures (using the center from the sphere match dl 11/00
c*******************************************************************************
c      minimization section
c*******************************************************************************
c
c     initialize first simplex vertex to starting position:
1297  do j = 1, options0%nvar
          simplx(1,j) = 0.
          tmp(j) = simplx(1,j)

      enddo

      lowest = 1

      simplx_energy(lowest) = sim_energy(tmp, best_energy,confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)

      !print *, 'first_energy=', simplx_energy(lowest)
      iter = 0
      irest = 0
      total_iterations = 0

      !initialize the random generator
      ! get the size of the seed vector and allocate it
      ! note that this vector will be differnt lenth depending on the
      ! compiler. e.g. gfortran, length 8; pgf95, length 34
      if ( .not. options0%allocated_seeds ) then
         call random_seed(size=n)
         allocate(options0%seeds(n))
         ! get a random seed vector form the random_seed function
         call random_seed(get=options0%seeds)
c        ! set the frist element of the random seed vector to that integer
c        options0%seeds(1) = options0%iseed
c        rand_val = 100000 ! I tried different values and this seemed resonable
c        ! each value of the array needs to be different. 
c        do i = 1,n
c            options0%seeds(i) = i * (rand_val - options0%iseed) ! this way of generating values to populate the array behaved well
         ! each value of the array needs to be different. 
         do i = 1,n
             options0%seeds(i) = options0%seeds(i) - options0%iseed ! subtract the seed from each value
         enddo
c        enddo
         options0%allocated_seeds = .true.
      end if
      ! pass this modifed seed vector to the random seed function
      call random_seed(put=options0%seeds)
      ! this seed will now be used to generate random numbers
c     print *, 'Initialize the random generator'

c      *********************************************************
c      *** loop over simplex restarts for full minimization ***
c      *********************************************************

 1300 continue
      old_energy = simplx_energy(lowest) ! get energy of last restart.
                         ! if it is the first time through, old_energy
                         ! becomes the starting position energy

c     initialize first simplex vertex to starting position, either
c       the starting position on the first pass, or the lowest 
c       vertex from previous simplex if a restart:
      do j = 1, options0%nvar
          simplx(1, j) = simplx(lowest, j)
      enddo

c     initialize the rest of the simplex vertices with
c       random displacements, and get energies at each vertex:
      min_energy = old_energy
      do i = 2, options0%nvar+1
        do j = 1, options0%nvar
c         generate random number from -1.0 to 1.0
c         random = 2.0 * (ran(options0%iseed) - 0.5)
          call random_number(random_num)
          random = 2.0 * (random_num - 0.5)
c         debug
c         print *, 'random_num_1=', random_num
c         print *, 'random_num_2=', random
c         multiply it by the allowed step size to generate a vertex
          if (j .le. 3) then
            simplx(i, j) = simplx(1, j) + options0%sim_rotstep * random
          else
            simplx(i, j) = simplx(1, j) + options0%sim_trnstep * random
          endif
          tmp(j) = simplx(i, j)
        enddo

        simplx_energy(i) = sim_energy(tmp, best_energy, confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)
        
         !print *, 'first_simplx_energy=', simplx_energy(i)
        if (simplx_energy(i) .lt. min_energy) then
          min_energy = simplx_energy(i)
          lowest = i
        endif
      enddo
c     perform the minimization

      call sim_simplex(simplx, simplx_energy, delta_energy, iter,
     &            best_energy,
     &           confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)


c     pick lowest energy vertex of simplex:
      min_energy = options0%big_energy
      do i = 1, options0%nvar+1
          if (simplx_energy(i) .lt. min_energy) then
            lowest = i
            min_energy = simplx_energy(i)
          endif
          !print *, 'energy=', simplx_energy(i)
      enddo

c     get final simplex energy
      do j = 1, options0%nvar
        tmp(j) = simplx(lowest, j)
      enddo
      new_energy = sim_energy(tmp, best_energy, confset,
     &            options0, db2lig, match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)

c     check to see if we need to restart

      if ((old_energy-new_energy .ge. options0%sim_need_to_restart)
     & .and. (iter .lt. options0%sim_itmax)) then
          irest = irest + 1
          !print *, '**********RESTART MINIMIZATION**************'
          !call doflush(6)
          goto 1300
      endif
      total_iterations = total_iterations + iter
      !print *, 'irest=', irest
c     call doflush(6)
      !stop 57
c     minimization is complete

      return
      end subroutine dockmin_sim
