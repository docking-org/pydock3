c  Monte Carlo minimization 
c  written by Reed Stein, October 2018
c
c
c
c
c
      subroutine dockmin_mc(old_energy, confset, options0,
     &           db2lig, match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits, mc_total_iterations)



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
                                                !vdw parameters per
                                                !vdwtype

c return variables
      integer, intent(inout) :: min_status
      real, intent(inout) ::  vscore, escore, gscore,
     &           ascore, pscore,
     &           rscore, rdscore, hscore, setscore
c     &           rscore, dscore, hscore, setscore


      integer :: atomstart, atomend !temp/loop variables
      real, intent(inout) :: old_energy
      integer, intent(inout) :: mc_total_iterations

      integer :: i, j
      integer :: mc_iter 
      integer :: acc_iter !number of MC steps that are accepted


c      real   :: old_energy, new_energy, prob
      real :: new_energy, prob, delta_energy, kT
c     new_energy      = energy at completion of most recent simplex
c     old_energy      = energy at completion of previous simplex

      real, dimension(options0%nvar) :: mc_matrix
      real, dimension(options0%nvar) :: tmp
      real :: random
c      real :: rand, random_num
      real :: random_num, rand_step

      real, dimension(3), intent(in) :: allminlimits, allmaxlimits
      real :: rigid_min

      real :: time2, time_vimelapsed
      integer :: n, countn

      ! initialize matrix to starting position
      do i = 1, options0%nvar
         mc_matrix(i) = 0.
      enddo 
      
      mc_total_iterations = 0
      acc_iter = 0
      
      ! initialize the random generator
      ! get the size of the seed vector and allocate it
      ! note that this vector will be different length depending on the
      ! compiler, e.g. gfortran, length 8; pgf95, length 34
      if ( .not. options0%allocated_seeds ) then
         call random_seed(size=n)
         allocate(options0%seeds(n))
         ! get a random seed vector from the random_seed function
         call random_seed(get=options0%seeds)
         do i = 1,n
             options0%seeds(i) = options0%seeds(i) - options0%mc_iseed ! subtract the seed from each value
         enddo
         options0%allocated_seeds = .true.
      end if

      ! pass this modified seed vector to the random seed function
      call random_seed(put=options0%seeds)
      ! this seed will now be used to generate random numbers

c     *****************************************************************
c     *** Metropolis Monte Carlo to explore space for minimization
c     *****************************************************************

 1300 continue
      call random_number(random_num)
      rand_step = CEILING(6.0 * (random_num))  
      do i = 1, options0%nvar
        if (i .eq. rand_step) then 
           call random_number(random_num)
           random = 2.0 * (random_num - 0.5)
           if (i .le. 3) then
              tmp(i) = mc_matrix(i) + options0%mc_rotstep * random
           else
              tmp(i) = mc_matrix(i) + options0%mc_trnstep * random
           endif
        else
           tmp(i) = mc_matrix(i) 
        endif
      enddo

      call euler2rot(tmp, match, matchnum)

      new_energy = rigid_min(confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
c     &           rscore, dscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)

c      write(6,*), "set_score = ",setscore

! if new energy is better than old, accept new matrix

      if (new_energy .le. old_energy) then
         old_energy = new_energy 
         acc_iter = acc_iter + 1
c         write(6,*), "new accepted energy = ",old_energy
         do i = 1, options0%nvar
            mc_matrix(i) = tmp(i)
         enddo
      else if (new_energy .gt. old_energy) then
         delta_energy = new_energy - old_energy
         call random_number(random_num)  
         kT = .001987206504191549 * options0%mc_temp
         prob = EXP(-delta_energy / kT)
         if (prob .ge. random_num) then
            old_energy = new_energy
            do i = 1, options0%nvar
               mc_matrix(i) = tmp(i)
            enddo
         endif 
      endif

      ! check to see if we need to restart
c      write(6,*), old_energy, mc_total_iterations, acc_iter
      if ((mc_total_iterations .lt. options0%mc_itmax) .and. 
     &    (acc_iter .lt. options0%mc_accpt)) then
c         write(6,*), old_energy, mc_total_iterations, acc_iter
         call doflush(6)
         mc_total_iterations = mc_total_iterations + 1
         goto 1300
      else ! need to update setscore so it's consistent

c         write(6,*), "Total MC moves = ",mc_total_iterations
c         write(6,*), "Total Accepted MC moves = ",acc_iter
         call euler2rot(mc_matrix, match, matchnum)
         new_energy = rigid_min(confset,
     &           options0, db2lig,
     &           match, matchnum,
     &           MAXOR, reccode,
     &           grids0, vdwmap0, phimap0, recdes0, ligdes0,
     &           gist0, solvmap0, rec_des_r0, sra, srb, maxtyv,
     &           vscore, escore, gscore, ascore, pscore,
     &           rscore, rdscore, hscore, setscore, min_status,
c     &           rscore, dscore, hscore, setscore, min_status,
     &           allminlimits, allmaxlimits)
      endif

      return
      end subroutine dockmin_mc
