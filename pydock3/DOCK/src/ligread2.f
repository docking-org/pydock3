
c read procedure for new ligand database format db2 rgc
      subroutine ligread2(db2lig, options0, 
     &    inputend, rigid_heavy_coords, 
     &    rigid_colors, rigid_heavy_count, fdlig, ligc, db2c,
     &    skip_processing)

      use db2type
      use optionstype
      use db2formats
      use filenums
      use spheres

      implicit none

      interface
c        logical function inlist(val,list,len) result(flag) 
         function inlist(val,list,len) 
              implicit none

              logical  inlist
              integer, intent(in) :: val,len
              integer, intent(in), allocatable, dimension(:) :: list
         end function inlist
      end interface

      interface
         function tokenize(line,token) 
              implicit none

              integer tokenize
              character (len=80), intent(in) :: line
              character (len=5), intent(in) :: token
         end function tokenize
      end interface

      type(db2), intent(inout) :: db2lig
      type(options), intent(inout) :: options0
      logical, intent(inout) :: inputend !EOF indicator; 1 if EOF reached, returned to call
      real, intent(inout), dimension(3, MAXPTS) :: rigid_heavy_coords !coordinates of heavy atoms in rigid
      integer, intent(inout), dimension(MAXPTS) :: rigid_colors !actually the only used ligand colors
      integer, intent(inout) :: rigid_heavy_count  !how many heavy atoms, used in matching
      integer (kind=8), intent(inout) :: fdlig !ligand file handle
      ! ligand count, incremented whenever an sdi ligand is read, reset to 1 when new sdi file is opened
      integer, intent(inout) :: ligc 
      ! gets incremented when a new db2 file is read
      integer, intent(inout) :: db2c
      logical, intent(in) :: skip_processing ! don't attempt to process or analyze the ligand, just read through it and update counters

      character (len=80) line !used to store one read line at a time
      character (len=78) temparbitrary !used to store arbitrary data temporarily
      character (len=5)  token
      integer inputstatus !used in gzread functions to get status of read
      character control !first character used to control program flow
      integer lines_seen !records how many seen within this segment of the file
      logical read_okay !used to keep track of a good ligand read. read until 
         !we find a good one or the file is empty
      integer lines_count, current, temp_childlen !used in group/conf trees
      integer curset, curidx, seterror, curseterror, curconf ! BT 4/23 used in detecting set/cluster line errors
      integer clustererror, matcherror, cluster_size ! BT 4/23 used for detecting cluster errors
      logical skipcluster ! BT 4/23 also used for detecting cluster errors
      integer children_seen, current_test, line_test !used in group/conf trees 
      integer data_count !used in group/conf trees
      integer temp_start, temp_end !used in reading in confs
      integer mostfreqconf, confsize1, tempconfsize !used in reading in confs
      integer broken_temp, hydro_temp !used in reading in sets
      real energy_temp !used in reading in sets
      real energy_temp_2 !used in reading in sets
      integer tempatom, tempbond, tempotheratom !used during setup for docking
      integer rigid_temp !used while reading in spheres
      integer temp_cluster, temp_match !used while reading in clusters/additional matching spheres
      integer temp_cur_match, temp_match_place !used while reading in clusters/additional matching spheres
      integer maxread, tempread !useful when reading in, implicit loop variable

c     these are all used for determining which confs are connected to
c     which.  This is for gist 
      integer indtc,indts,indsc, indcord, indbond, confnum, atomnum

c     !integer, dimension(5000,2) ::  bondconflist
      integer, allocatable, dimension(:,:) ::  bondconflist
      
      integer, allocatable, dimension(:) :: currentlist
      integer, allocatable, dimension(:) :: newlist
      integer, allocatable, dimension(:) :: seenlist
      integer, allocatable, dimension(:) :: conf_count
      integer, allocatable, dimension(:) :: conf_connect_count ! count the number of connections to other confs 
      integer currentint, newint,seenint, count1, count2
      integer count3
      integer sort_count
      logical flaginlist
      integer allocate_stat  ! for allocating
      integer read_stat ! for line readling
      logical mol_okay ! for the M section
      logical flag_M ! for molecule that follows a broken molecule
c without E, then first M already read in.
      integer broken_line_count, maxbrokenline
      integer countA, countB, countX ! the counter for each section
      integer countC, countS, countR ! the counter for each section

      energy_temp   = 0.0  ! intitize value to zero
      energy_temp_2 = 0.0  

      lines_seen = 0
      line = ' ' !placeholder that forces new line to always be read in
      token = ' '
      flag_M = .false.
      read_okay = .true.
      mol_okay = .true.

      if (skip_processing) then
        do while (line(1:1) .ne. 'E')
          call gzread(fdlig, line, 80, inputstatus)
          if (inputstatus .le. -1) then
            if (.not. options0%sdi) then
              exit
            else
              read(SDIFILE, '(a255)', err=1111, end=1111) 
     &          options0%ligfil !next filename
              call gzopen(fdlig, 'r', options0%ligfil, inputstatus) 
              db2c = db2c + 1
              continue
1111          exit
            endif
          endif
        end do
        ligc = ligc + 1
        return
      endif


      do while (read_okay)
c       write(OUTDOCK,*) "read_okay"
        do while (line(1:1) .ne. 'E') !detects end of molecule
          !fdlig is always the ligand file, already opened for us already
          !gzread is used since files are always gzipped
          !any errors are indicated by inputstatus
c         write(OUTDOCK,*) "while not E"
          if (.not.flag_M) then
            call gzread(fdlig, line, 80, inputstatus) 
c           if (inputstatus .le. 0) then !means the end of file has been reached
            if (inputstatus .le. -1) then !means the end of file has been reached
              inputend = .true.
              exit !break out of do loop
            endif
          else
            !write(OUTDOCK, *) line
c           write(OUTDOCK, *)
c    &               'The first line in M section is already read.'
            flag_M = .false.
          endif
          !if we are here, there is real data in the line buffer, so we should 
          !process it accordingly (first character always indicates line type)
          control = line(1:1) !grab first character only
          if (control .eq. 'M') then !molecule line
c           write(OUTDOCK,*) "line frist char M"
            if (lines_seen .eq. 0) then
              if (tokenize(line,token) .ne. 10) then
                write(OUTDOCK, *)
     &               'The first line in M section is broken', line
c    &               tokenize(line,token)
                read_okay = .false.
                mol_okay  = .false.
                exit !break out of do loop
              endif
              read(line, DB2M1, iostat=read_stat)
     &            db2lig%refcod, db2lig%protcod, 
     &            db2lig%total_atoms, db2lig%total_bonds, 
     &            db2lig%total_coords, db2lig%total_confs, 
     &            db2lig%total_sets, rigid_heavy_count, 
     &            db2lig%mlines, db2lig%total_clusters
              if (read_stat .ne. 0) then
                write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                read_okay = .false.
                exit !break out of do loop
              endif 
              lines_seen = 1
              !want to dynamically allocate storage now
              !write(OUTDOCK,*) "allocate_db2"
              !call doflush(OUTDOCK)
              call allocate_db2(db2lig, OUTDOCK, 
     &            db2lig%total_sets * 3, ! BT 4/23 see main comment
     &            db2lig%total_coords, db2lig%total_atoms, 
     &            db2lig%total_bonds,
     &            db2lig%total_confs, db2lig%mlines, 
     &            db2lig%total_clusters,
     &            options0%dockovalent, options0%nsav)

              ! BT 4/23
              ! We use set_conf_len to determine whether a set exists or not
              ! must initialize this list with zeros beforehand, assuming all sets don't exist until we encounter them
              db2lig%set_conf_len(:) = 0

               if (allocated(currentlist)) then
                   deallocate(currentlist)
               endif
               if (allocated(newlist)) then
                   deallocate(newlist)
               endif
               if (allocated(seenlist)) then
                   deallocate(seenlist)
               endif
               if (allocated(bondconflist)) then
                   deallocate(bondconflist)
               endif
               if (allocated(conf_count)) then
                   deallocate(conf_count)
               endif
               if (allocated(conf_connect_count)) then
                   deallocate(conf_connect_count)
               endif
              allocate(currentlist(10*db2lig%total_confs),
     &                 stat=allocate_stat) 
              allocate(newlist(10*db2lig%total_confs),
     &                 stat=allocate_stat) 
              allocate(seenlist(db2lig%total_confs),
     &                 stat=allocate_stat) 
              allocate(bondconflist(db2lig%total_bonds,2),
     &                 stat=allocate_stat)
              allocate(conf_count(db2lig%total_confs),  ! used to
              ! count the the number occurences of this conf amoung the
              ! sets. 
     &                 stat=allocate_stat)
              allocate(conf_connect_count(db2lig%total_confs), 
     &                 stat=allocate_stat)
              countA = 0
              countB = 0
              countX = 0
              countC = 0
              countS = 0
              countR = 0
            else if (lines_seen .eq. 1) then !the second line
              read(line, DB2M2, iostat=read_stat) db2lig%total_charge, 
     &            db2lig%total_polar_solv, db2lig%total_apolar_solv, 
     &            db2lig%total_solv, db2lig%total_surfarea
              if (read_stat .ne. 0) then
                write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                read_okay = .false.
                exit !break out of do loop
              endif
              lines_seen = 2
            else if (lines_seen .eq. 2) then !the smiles line
              read(line, DB2M3, iostat=read_stat) db2lig%smiles
              if (read_stat .ne. 0) then
                write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                read_okay = .false.
                exit !break out of do loop
              endif
              lines_seen = 3
            else if (lines_seen .eq. 3) then !the long name line
              read(line, DB2M3, iostat=read_stat) db2lig%molname
              if (read_stat .ne. 0) then
                write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                read_okay = .false.
                exit !break out of do loop
              endif
              lines_seen = 4
            else if (lines_seen .ge. 4) then !arbitrary lines
              ! jklyu, modified on Aug. 26th, 2019 
              !read(line, DB2M3, iostat=read_stat) temparbitrary
              read(line, DB2M3, iostat=read_stat) db2lig%rig_frag_code
              temparbitrary = db2lig%rig_frag_code
              !write(OUTDOCK, *) 'rig_frag_code: ', db2lig%rig_frag_code
              if (read_stat .ne. 0) then
                write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                read_okay = .false.
                exit !break out of do loop
              endif
              lines_seen = lines_seen + 1
              !write(OUTDOCK, *) 'ZINC ID: ', db2lig%refcod
              if (lines_seen .gt. db2lig%mlines) then
              write(OUTDOCK, *)
     &               'lines_seen in M section>=mlines', line
                read_okay = .false.
                exit !break out of do loop
              endif
              ! jklyu, modified on Aug. 26th, 2019 
              db2lig%arbitrary(lines_seen - 4) = temparbitrary !put in the list
            endif
            if (lines_seen .eq. db2lig%mlines) then !read all mlines, continue to A section
              lines_seen = 0
            endif
          else if (control .eq. 'A') then !atom line
            lines_seen = lines_seen + 1
            countA = countA + 1
            if (tokenize(line,token) .ne. 10) then
              write(OUTDOCK, *)
     &             'The line in A section is broken', line,
     &             tokenize(line,token)
              read_okay = .false.
              !mol_okay  = .false.
              exit !break out of do loop
            endif
            if (lines_seen .le. db2lig%total_atoms) then
              read(line, DB2ATOM, iostat=read_stat)
     &              db2lig%atom_num(lines_seen), 
     &              db2lig%atom_name(lines_seen), 
     &              db2lig%atom_type(lines_seen),          
     &              db2lig%atom_vdwtype(lines_seen), 
     &              db2lig%atom_color(lines_seen),          
     &              db2lig%atom_charge(lines_seen), 
     &              db2lig%atom_polsolv(lines_seen),          
     &              db2lig%atom_apolsolv(lines_seen),
     &              db2lig%atom_totalsolv(lines_seen), 
     &              db2lig%atom_surfarea(lines_seen)
              if (read_stat .ne. 0) then
                write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                read_okay = .false.
                exit !break out of do loop
              endif
              if (lines_seen .ne. db2lig%atom_num(lines_seen)) then
                write(OUTDOCK, *)
     &              'atom line mismatch file', line, lines_seen
                read_okay = .false.
                exit !break out of do loop
              endif
            else
              write(OUTDOCK, *) 
     &             'poorly formatted line1 in ligand file', line
              read_okay = .false.
              exit !break out of do loop
            endif
            if (lines_seen .eq. db2lig%total_atoms) then
              lines_seen = 0 !reset for bonds
            endif
          else if (control .eq. 'B') then !bond line
            lines_seen = lines_seen + 1
            countB = countB + 1
            if (tokenize(line,token) .ne. 4) then
              write(OUTDOCK, *)
     &             'The line in B section is broken', line,
     &             tokenize(line,token)
              read_okay = .false.
              !mol_okay  = .false.
              exit !break out of do loop
            endif
            if (lines_seen .le. db2lig%total_bonds) then
              read(line, DB2BOND, iostat=read_stat)
     &            db2lig%bond_num(lines_seen), 
     &            db2lig%bond_start(lines_seen), 
     &            db2lig%bond_end(lines_seen),
     &            db2lig%bond_type(lines_seen)
              if (read_stat .ne. 0) then
                write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                read_okay = .false.
                exit !break out of do loop
              endif
              if (lines_seen .ne. db2lig%bond_num(lines_seen)) then
                write(OUTDOCK, *) 
     &               'bond line mismatch file', line, lines_seen
                read_okay = .false.
                exit !break out of do loop
              endif
            else
              write(OUTDOCK, *) 
     &             'poorly formatted line2 in ligand file', line
              read_okay = .false.
              exit !break out of do loop
            endif
            if (lines_seen .eq. db2lig%total_bonds) then
              lines_seen = 0 !reset for coords
            endif
          else if (control .eq. 'X') then !coordinate line
c           write(OUTDOCK, *) line
c           call doflush(OUTDOCK)
            lines_seen = lines_seen + 1
            countX = countX + 1
            if (tokenize(line,token) .ne. 6) then
              write(OUTDOCK, *)
     &             'The line in X section is broken', line,
     &             tokenize(line,token)
              read_okay = .false.
              !mol_okay  = .false.
              exit !break out of do loop
            endif
            if (lines_seen .le. db2lig%total_coords) then
              read(line, DB2COORD, iostat=read_stat)
     &              db2lig%coord_index(1, lines_seen), 
     &              db2lig%coord_index(2, lines_seen), 
     &              db2lig%coord_index(3, lines_seen), 
     &              db2lig%coords(1, lines_seen),
     &              db2lig%coords(2, lines_seen), 
     &              db2lig%coords(3, lines_seen)
              if (read_stat .ne. 0) then
                write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                read_okay = .false.
                exit !break out of do loop
              endif
            else
              write(OUTDOCK, *) 
     &             'poorly formatted line3 in ligand file', line
              read_okay = .false.
              exit !break out of do loop
            endif
            if (lines_seen .ne. db2lig%coord_index(1, lines_seen)) then
              write(OUTDOCK, *) 
     &             'poorly formatted line4 in ligand file', line
              read_okay = .false.
              exit !break out of do loop
            endif
            if (lines_seen .eq. db2lig%total_coords) then
              lines_seen = 0 !reset for rigidmatchingheavy stuff
              current = 0
            endif
          else if (control .eq. 'R') then !rigid matching heavy line
            lines_seen = lines_seen + 1
            countR = countR + 1
            if (tokenize(line,token) .ne. 5) then
              write(OUTDOCK, *)
     &             'The line in R section is broken', line,
     &             tokenize(line,token)
              read_okay = .false.
              !mol_okay  = .false.
              exit !break out of do loop
            endif
            if (lines_seen .le. rigid_heavy_count) then
              read(line, DB2RIGID, iostat=read_stat)
     &              rigid_temp, rigid_colors(lines_seen),
     &              rigid_heavy_coords(1, lines_seen),
     &              rigid_heavy_coords(2, lines_seen),
     &              rigid_heavy_coords(3, lines_seen)
              if (read_stat .ne. 0) then
                write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                read_okay = .false.
                exit !break out of do loop
              endif
            else
              write(OUTDOCK, *) 
     &             'poorly formatted line3r in ligand file', line
              read_okay = .false.
              exit !break out of do loop
            endif
            if (lines_seen .ne. rigid_temp) then
              write(OUTDOCK, *) 
     &             'poorly formatted line4r in ligand file', line
              read_okay = .false.
              exit !break out of do loop
            endif
            if (lines_seen .eq. rigid_heavy_count) then
              lines_seen = 0 !reset for group stuff
              current = 0
            endif
          else if (control .eq. 'C') then !conf line
            !format is much simpler now. just conf# xstart xend
c           write(OUTDOCK, *) line
c           call doflush(OUTDOCK)
            lines_seen = lines_seen + 1
            countC = countC + 1
            if (tokenize(line,token) .ne. 3) then
              write(OUTDOCK, *)
     &             'The line in C section is broken', line,
     &             tokenize(line,token)
              read_okay = .false.
              mol_okay  = .false.
              exit !break out of do loop
            endif
            if (lines_seen .le. db2lig%total_confs) then !read first line
              read(line, DB2CONF, iostat=read_stat)
     &              current, temp_start, temp_end
              if (read_stat .ne. 0) then
                write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                read_okay = .false.
                exit !break out of do loop
              endif
              db2lig%conf_coord(1, current) = temp_start
              db2lig%conf_coord(2, current) = temp_end
            else
              write(OUTDOCK, *) 
     &             'poorly formatted line66 in ligand file', line
              read_okay = .false.
              exit !break out of do loop
            endif
            if (lines_seen .eq. db2lig%total_confs) then
              lines_seen = 0 !reset for S
              current = 0
            endif
          else if (control .eq. 'S') then !set line
            !in this block, lines_seen is reset at the end of each set of confs
            !complicated two part read structure to read trees which can have
            !almost any number of sub-components and are written in 80charlimits
            if ((lines_seen .eq. 0) .and. (current .eq. 0)) then !read first line
              lines_seen = lines_seen + 1
              countS = countS + 1
              !if (tokenize(line,token) .ne. 6) then
              !if (tokenize(line,token) .ne. 7) then
              !  write(OUTDOCK, *)
              !       'The first line in S section is broken', line
              !  read_okay = .false.
              !  exit !break out of do loop
              !endif
              !write(OUTDOCK, *) 'tokenize', tokenize(line,token) 
              !if (tokenize(line,token) .ne. 7) then
              if (tokenize(line,token) .eq. 7) then
                read(line, DB2SET1_S, iostat=read_stat)
     &                current, lines_count, temp_childlen,
     &                broken_temp, hydro_temp, energy_temp,
     &                energy_temp_2
                if (read_stat .ne. 0) then
                  write(OUTDOCK, *)
     &                 'The line is broken, failed to be read', line
                  read_okay = .false.
                  exit !break out of do loop
                endif
c               write(OUTDOCK, *) 'current', current
c               write(OUTDOCK, *) 'tE', energy_temp
c               write(OUTDOCK, *) 'pE', energy_temp_2
                ! BT 4/23
                ! For the most part I don't need to change the set reading code to
                ! adjust for the errors. However, I do need to be able to recognize
                ! duplicates so they are not included in the set count
                ! Confused? ctrl-F "BT 4/23 Main Comment"
                if (db2lig%set_conf_len(current) .ne. 0) then
                  ! this is a duplicated set. No worry, just make a note of it in countS
                  countS = countS - 1
                endif

                ! Do a bounds check here, just in case the numbering anomaly is too extreme
                if (db2lig%total_sets * 3 .lt. current) then
                  write(OUTDOCK, *) 
     &   'Detected an extreme set numbering anomaly, skipping ligand'
                  read_okay = .false.
                  exit
                endif
                db2lig%set_conf_len(current) = temp_childlen
                db2lig%set_broken(current) = broken_temp
                db2lig%set_energy(current) = energy_temp * 
     &              options0%iscale !scale the internal energy here
                db2lig%set_total_strain(current) = energy_temp
                db2lig%set_max_strain(current) = energy_temp_2
                if ((energy_temp .ge. options0%total_strain) .or.
     &            (energy_temp_2 .ge. options0%max_strain)) then
                  db2lig%set_above_strain(current) = 1
                else
                  db2lig%set_above_strain(current) = 0
                endif
                children_seen = 0
              !else if (tokenize(line,token) .ne. 6) then
              else if (tokenize(line,token) .eq. 6) then
                read(line, DB2SET1, iostat=read_stat)
     &                current, lines_count, temp_childlen,
     &                broken_temp, hydro_temp, energy_temp
c               write(*,*) current, lines_count, temp_childlen, 
c    &                     broken_temp, hydro_temp, energy_temp
c               write(OUTDOCK, *) 'current', current
c               write(OUTDOCK, *) 'tE', energy_temp
                if (read_stat .ne. 0) then
                  write(OUTDOCK, *)
     &                 'The line is broken, failed to be read', line
                  read_okay = .false.
                  exit !break out of do loop
                endif

                if (db2lig%set_conf_len(current) .ne. 0) then
                  countS = countS - 1
                endif

                if (db2lig%total_sets * 3 .lt. current) then
                  write(OUTDOCK, *) 
     &   'Detected an extreme set numbering anomaly, skipping ligand'
                  read_okay = .false.
                  exit
                endif

                db2lig%set_conf_len(current) = temp_childlen
                db2lig%set_broken(current) = broken_temp
                db2lig%set_energy(current) = energy_temp * 
     &              options0%iscale !scale the internal energy here
c               write (*,*) db2lig%set_energy(current) 
                db2lig%set_total_strain(current) = energy_temp
                db2lig%set_max_strain(current) = 0.0
                db2lig%set_above_strain(current) = 0.0
                children_seen = 0
              else
                write(OUTDOCK, *)
     &               'The first line in S section is broken', line
                read_okay = .false.
                exit !break out of do loop
              endif
            else
              lines_seen = lines_seen + 1
              !next line figures out how many to read
              maxread = min(8, 
     &            db2lig%set_conf_len(current) - children_seen)
              if (tokenize(line,token) .ne. (3+maxread)) then
                write(OUTDOCK, *)
     &               'The line in S section is broken', line
                read_okay = .false.
                exit !break out of do loop
              endif
              read(line, DB2SET2, iostat=read_stat)
     &              current_test, line_test, data_count, 
     &            (db2lig%set_conf(current, children_seen + tempread), 
     &            tempread = 1, maxread)
              current_test = current_test
              
              if (read_stat .ne. 0) then
                write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                read_okay = .false.
                exit !break out of do loop
              endif
              if ((current .ne. current_test) .or.
     &            (line_test .ne. lines_seen - 1)) then !problem
                write(OUTDOCK, *) 'set line mismatch', line, 
     &              lines_seen, line_test, current, current_test
                read_okay = .false.
                exit !break out of do loop
              endif
              children_seen = children_seen + data_count
            endif
            if (lines_seen .eq. lines_count + 1) then !done with current group
              lines_seen = 0 !reset for next group
              current = 0 !reset for next group
              clustererror = 0
              seterror = 0
              curseterror = 0
              matcherror = 0
              skipcluster = .false.
            endif
          else if (control .eq. 'D') then !set line
            !read the clusters, first a cluster line then a few additional
            !matching sphere lines for each cluster.
            if ((lines_seen .eq. 0)) then !read first line
              lines_seen = lines_seen + 1
              if (tokenize(line,token) .ne. 6) then
                write(OUTDOCK, *)
     &               'The first line in D section is broken', line
                read_okay = .false.
                exit !break out of do loop
              endif

              !current = current - clustererror

c BT 4/23 adjust current index to account for any deleted clusters
              curidx = current + 1 - clustererror 

              read(line, DB2CLUSTER, iostat=read_stat) current_test,
     &              db2lig%cluster_set_start(curidx),
     &              db2lig%cluster_set_end(curidx), temp_cur_match,
     &              db2lig%cluster_match_start(curidx),
     &              db2lig%cluster_match_end(curidx)

c BT 4/23 Main Comment
c New mol2db2 seems to produce db2 files that have bizarre set numbering/duplication errors, i.e:
c M sets = 6
c S 1
c S 2
c S 3
c S 7
c S 5
c S 6
c S 7
c Note that 4 gets skipped over- the amount of set numbering that is skipped depends on how many
c duplicates are generated.
c When the duplicates are discarded, there are in fact as many unique sets as described in the header
c of the db2 file. To compensate for this, we allocate more space than usual in the db2 set arrays
c (i.e set_conf_len, set_conf). This allows us to read in mis-numbered sets (sets with higher # than
c the declared total_sets) without an error. Then, as we create clusters, simply take note of any
c sets whose number has been skipped over, and move the rest of the set arrays over to compensate. We 
c also need to adjust any clusters that reference non-existent sets, and possibly even delete clusters
c as we go along. Seems to work pretty well!
c ctrl-F BT 4/23 to see all sections of code pertinent to this update
              curseterror = 0
              cluster_size = 
     &          (db2lig%cluster_set_end(curidx)
     &          -db2lig%cluster_set_start(curidx))

c Sometimes sets will be normal, but clusters will refer to sets that
c don't exist. Catch this here and compensate
              if ((db2lig%cluster_set_start(curidx).gt.
     &            (db2lig%total_sets*3)).or.
     &            (db2lig%cluster_set_end(curidx).gt.
     &            (db2lig%total_sets*3))) then
                  cluster_size = 0
                  db2lig%cluster_set_start(curidx) = 0
                  db2lig%cluster_set_end(curidx) = 0
              endif
                
              
              do curset = db2lig%cluster_set_start(curidx), 
     &                    db2lig%cluster_set_end(curidx)
                if (db2lig%set_conf_len(curset) .eq. 0) then
                  curseterror = curseterror + 1
                else if ((seterror + curseterror) .gt. 0) then
                  curidx = curset-curseterror-seterror
                  db2lig%set_conf_len(curidx) = 
     &            db2lig%set_conf_len(curset)
                  do curconf = 1, db2lig%set_conf_len(curset)
                    db2lig%set_conf(curidx, curconf) = 
     &              db2lig%set_conf(curset, curconf)
                  end do
c                  db2lig%set_conf(curidx) = 
c     &            db2lig%set_conf(curset)
                  db2lig%set_broken(curidx) = 
     &            db2lig%set_broken(curset)
                  db2lig%set_above_strain(curidx) = 
     &            db2lig%set_above_strain(curset)
                  db2lig%set_total_strain(curidx) = 
     &            db2lig%set_total_strain(curset)
                  db2lig%set_max_strain(curidx) = 
     &            db2lig%set_max_strain(curset)
                  db2lig%set_energy(curidx) = 
     &            db2lig%set_energy(curset)
                endif

              end do
                
              seterror = seterror + curseterror
              curidx = current + 1 - clustererror

              if (curidx .gt. 1) then
                db2lig%cluster_set_start(curidx) =
     &          db2lig%cluster_set_end(curidx-1) + 1
              else
                db2lig%cluster_set_start(curidx) = 1
              endif

              db2lig%cluster_set_end(curidx) =
     &          db2lig%cluster_set_start(curidx) + 
     &          (cluster_size - curseterror)

              if (curseterror .eq. cluster_size) then
                clustererror = clustererror + 1
                matcherror = matcherror + temp_cur_match
                skipcluster = .true.
              else if (matcherror .gt. 0) then
                curidx = current + 1 - clustererror
                db2lig%cluster_match_start(curidx) =
     &          db2lig%cluster_match_start(curidx) - matcherror
                db2lig%cluster_match_end(curidx) =
     &          db2lig%cluster_match_end(curidx) - matcherror
              endif
c End of BT main error correction logic

              if (read_stat .ne. 0) then
                write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                read_okay = .false.
                exit !break out of do loop
              endif

              if ((current + 1 - clustererror) .ne.
     &            (current_test - clustererror))  then !problem
                write(OUTDOCK, *) 'cluster line mismatch', line, 
     &              current, current_test
                read_okay = .false.
                exit !break out of do loop
              endif

              current = current + 1 !increment now, so correct for add match sph
              temp_match_place = db2lig%cluster_match_start(current)
              if (temp_cur_match .eq. 0) then !no matching spheres for a cloud
                lines_seen = 0 !so we read the next line correctly
              endif
            else !lines_seen is 1 or greater, means we're doing the add match spheres
              ! BT 4/23 skip over copying cluster data if cluster is deleted
              if (.not. skipcluster) then
                !reuse the DB2RIGID line for reading matching spheres
                if (tokenize(line,token) .ne. 5) then
                  write(OUTDOCK, *)
     &               'The following line in D section is broken', line
                  read_okay = .false.
                  exit !break out of do loop
                endif
                read(line, DB2RIGID, iostat=read_stat) temp_match, 
     &             db2lig%addmatch_color(temp_match_place),
     &             db2lig%addmatch_coord(1, temp_match_place),
     &             db2lig%addmatch_coord(2, temp_match_place),
     &             db2lig%addmatch_coord(3, temp_match_place)
                if (read_stat .ne. 0) then
                  write(OUTDOCK, *)
     &               'The line is broken, failed to be read', line
                  read_okay = .false.
                  exit !break out of do loop
                endif
                if (temp_match .ne. temp_match_place) then !we have a problem
                  write(OUTDOCK, *) 'addmatchsph line mismatch', line, 
     &              temp_match, temp_match_place
                  read_okay = .false.
                  exit !break out of do loop
                endif
              endif
              temp_match_place = temp_match_place + 1
              if (temp_match_place .gt. 
     &            db2lig%cluster_match_end(current-clustererror)) then !done with reading matching spheres
                lines_seen = 0 !reset this so go back to reading the cluster 
                skipcluster = .false.
                if (current .eq. 
     &            (db2lig%total_clusters-clustererror)) then ! completely done w/ clusters
                  current = 0 !reset for next step 
                  ! BT 4/23 correct total_clusters to reflect actual number
                  ! Unlike total_sets, this number may be reported incorrectly in the header
                  db2lig%total_clusters = 
     &            db2lig%total_clusters - clustererror
                endif
              endif
            endif
          else if (control .eq. 'E') then !end line
            !To solve early E problem
            !Check if we read in enough lines before initializing the
            !mol
            if (countA .ne. db2lig%total_atoms) then
              write(OUTDOCK, *) "countA != total atoms"
              read_okay = .false.
            endif
            if (countB .ne. db2lig%total_bonds) then
              write(OUTDOCK, *) "countB != total bonds"
              read_okay = .false.
            endif
            if (countX .ne. db2lig%total_coords) then
              write(OUTDOCK, *) "countX != total coords"
              read_okay = .false.
            endif
            if (countC .ne. db2lig%total_confs) then
              write(OUTDOCK, *) "countC != total confs"
              read_okay = .false.
            endif
            if (countR .ne. rigid_heavy_count) then
              write(OUTDOCK, *) "countR != total coords"
              read_okay = .false.
            endif
            if (countS .ne. db2lig%total_sets) then
              write(OUTDOCK, *) "countS != total sets"
              read_okay = .false.
            endif
            ligc = ligc + 1
            exit !break out of do loop since we've read in one ligand to dock
          else 
            write(OUTDOCK, *) 
     &           'poorly formatted line5 in ligand file', line, control
            read_okay = .false.
            exit !break out of do loop
          endif
        enddo
        if (.not. inputend) then !not reached end of file yet
          if (.not. read_okay) then !there was a problem while reading
            !advance until end of current broken molecule
            if (.not. mol_okay) then
              write(OUTDOCK, *) "fatal issue with the line length."
            else
              write(OUTDOCK, *) "skipping molecule: ", 
     &            db2lig%refcod, db2lig%protcod
            endif
            broken_line_count = 1
            maxbrokenline = 50000
            call gzread(fdlig, line, 80, inputstatus)
            do while (line(1:1) .ne. 'E' .and. line(1:1) .ne. 'M') !detects end of molecule
              call gzread(fdlig, line, 80, inputstatus)
              !write(OUTDOCK,*) line 
              !write(OUTDOCK,*) "EM check", broken_line_count
c             if (inputstatus .le. 0) then !means the end of file
c               !write(OUTDOCK, *) "istat: ", inputstatus
c               !inputend = .true.
c               !exit !break out of do loop
c             endif
              if (broken_line_count .ge. maxbrokenline) then !no infinite loop
c               write(OUTDOCK, *) "Program terminating . . ."
                write(OUTDOCK, *) "going to next file . . ."
                write(OUTDOCK, *) "exceed the max number of BrokenLines"
     &                            , maxbrokenline
c               write(OUTDOCK, *) "remove broken db2 file from index"
c               write(OUTDOCK, *) "and resubmit"
                inputend = .true.
                read_okay = .false.
                exit ! break out of E loop
                !stop
              endif
              broken_line_count = broken_line_count + 1
            enddo
            if (line(1:1) .eq. 'M') then
                if (tokenize(line,token) .eq. 10) then
                   flag_M = .true.
                endif
                !backspace(fdlig)
            else
                line = ' ' !otherwise loop doesn't restart since it sees E
            endif
            write(OUTDOCK, *) "BrokenLines: ", broken_line_count
            read_okay = .true. !reset and keep reading next ligand
            mol_okay = .true.
            lines_seen = 0 !really, everything has to be reset
            cycle !skips the rest of the do loop, don't process partial ligand
          endif
        endif
        if ((.not. inputend) .and. (read_okay)) then 
          !put all code to prep this ligand for docking here
          !lots of variables need calculated and put in various other places
          !need to count total heavy atoms in ligand--output foreach ligand
          tempatom = 1
          db2lig%total_heavy_atoms = 0
          db2lig%ligcharge = 0.
          do while (tempatom .le. db2lig%total_atoms)
            if ((db2lig%atom_vdwtype(tempatom) .ne. 6) .and. 
     &          (db2lig%atom_vdwtype(tempatom) .ne. 7)) then
              db2lig%total_heavy_atoms = db2lig%total_heavy_atoms + 1
            endif
            db2lig%ligcharge = db2lig%ligcharge + 
     &          db2lig%atom_charge(tempatom)
            tempatom = tempatom + 1
          enddo 
          !now count heavy atom valence, in other words how many heavy atoms
          !is each atom connected to. store in db2lig%atom_heavy_valence
          tempatom = 1
          do while (tempatom .le. db2lig%total_atoms)
            db2lig%atom_heavy_valence(tempatom) = 0. !reset to zero
            tempatom = tempatom + 1
          enddo 
          !now go through bonds. whenever one atom is heavy, increment the other
          !only do this if we want to use the valence later
          !write (6,*) 'xxxrvp', options0%recdes_valence_penalty
          if (options0%recdes_valence_penalty) then
            tempbond = 1
            do while (tempbond .le. db2lig%total_bonds)
              tempatom = db2lig%bond_start(tempbond)
              tempotheratom = db2lig%bond_end(tempbond)
              if ((db2lig%atom_vdwtype(tempatom) .ne. 6) .and. 
     &            (db2lig%atom_vdwtype(tempatom) .ne. 7)) then !tempatom is heavy
                !so tempotheratom gets incremented
                db2lig%atom_heavy_valence(tempotheratom) = 
     &              db2lig%atom_heavy_valence(tempotheratom) + 1
              endif
              !now go the other way
              if ((db2lig%atom_vdwtype(tempotheratom) .ne. 6) .and. 
     &            (db2lig%atom_vdwtype(tempotheratom) .ne. 7)) then !heavy
                !so tempatom gets incremented
                db2lig%atom_heavy_valence(tempatom) = 
     &              db2lig%atom_heavy_valence(tempatom) + 1
              endif
              tempbond = tempbond + 1
            enddo
            !next part caps the heavy atom valence at 3, above this it 
            !doesn't get any receptor desolvation anyway.
            tempatom = 1
            do while (tempatom .le. db2lig%total_atoms)
              if (db2lig%atom_heavy_valence(tempatom) .gt. 5) then
                db2lig%atom_heavy_valence(tempatom) = 5
              endif
              tempatom = tempatom + 1
            enddo 
            !write(OUTDOCK, *) 'hv', db2lig%atom_heavy_valence !debug output
          endif
          ! WHAT I WANT TO DO.  find out which confs (rigid segement) are connected. 
          ! I can do this by looking at the bond information. 
          ! look at which bonds spans two confs. id which segments
          ! contain that atoms in the bonds. 
          ! loop over all of the sets of confs and store which confs
          ! are conected to witch
          ! we have to loop over each set because one set is a full
          ! molcule, there is overlap in atoms numbers for the sets.
          ! for example sets 5 and 6 might have the same atoms. but 5 is
          ! connected to 2 while 6 is connected to 3.  
          do indtc = 1, db2lig%total_confs ! loop over all rigid segements (confs)
             !write(*,*) indtc
             db2lig%conf_conection(indtc) = indtc ! set conf_conection equal to itself
             db2lig%conf_sorted(indtc) = indtc ! set to itself, if options0%use_gist is false
          enddo
c            seenlist(1) = 1
c            seenint = 1 
          ! find the most fequent conf among the sets. 
          ! if there is a tie, pick the one with the most atoms. 

          if (options0%use_gist) then 
          
          do indtc = 1, db2lig%total_confs ! loop all confs
              conf_count(indtc) = 0
          enddo

          do indtc = 1, db2lig%total_confs ! loop all confs
              currentlist(indtc) = 0
          enddo

          do indts = 1, db2lig%total_sets ! loop over all sets of confs
              do indsc = 1, db2lig%set_conf_len(indts)
                  confnum = db2lig%set_conf(indts,indsc)
                  conf_count(confnum) = conf_count(confnum)+1
              enddo
          enddo

          mostfreqconf = 1
          temp_start = db2lig%conf_coord(1,mostfreqconf)
          temp_end = db2lig%conf_coord(2,mostfreqconf)
          confsize1 = temp_end - temp_start
          tempconfsize = 0
          do indtc = 1, db2lig%total_confs ! loop all confs
              if (conf_count(mostfreqconf).lt.conf_count(indtc)) then 
                  mostfreqconf = indtc
                  temp_start = db2lig%conf_coord(1,mostfreqconf)
                  temp_end = db2lig%conf_coord(2,mostfreqconf)
                  confsize1 = temp_end - temp_start
              else if (conf_count(mostfreqconf).eq.conf_count(indtc))
     &        then
                  temp_start = db2lig%conf_coord(1,mostfreqconf)
                  temp_end = db2lig%conf_coord(2,mostfreqconf)
                  tempconfsize = temp_end - temp_start
                  if (tempconfsize.gt.confsize1) then
                       mostfreqconf = indtc
                       confsize1 = tempconfsize
                  endif
              endif   
          enddo
          
c         write(*,*) "mostfreqconf, confsize1 = ", 
c    &      mostfreqconf, confsize1 
          
          
          do indtc = 1, db2lig%total_confs ! loop all confs
              conf_connect_count(indtc) = 0
              db2lig%conf_sorted(indtc) = 0
          enddo


          !sort_count = 2
          sort_count = 1
          count3 = 1
          do indts = 1, db2lig%total_sets ! loop over all sets of confs
c            write(*,*) "set = ", indts
             ! frist I fined what confs are connected with in the set. 
             do indsc = 1, db2lig%set_conf_len(indts) ! loop over each conf in that set. 
c                ! loop over bonds.  when atom is found in bond
                 ! asign the conf.
                 confnum = db2lig%set_conf(indts,indsc)
c                write(*,*) "conf = ", confnum
                 temp_start = db2lig%conf_coord(1,confnum)
                 temp_end = db2lig%conf_coord(2,confnum)
                 do indcord = temp_start, temp_end 
c                  write(*,*) "coord_index", 
c    &                 db2lig%coord_index(:,indcord)
                   atomnum = db2lig%coord_index(2,indcord) ! look up atom number
c                  write(*,*) atomnum
                   ! check if atom is a hydrogen.  If it is, then move
                   ! on to the next atom.  
c                  call doflush(OUTDOCK)
c                  if ((options0%gistH.eq.1).and.
c    &                 ((db2lig%atom_vdwtype(atomnum).eq.7).or.
c    &                 (db2lig%atom_vdwtype(atomnum).eq.6))) then
c                  if  ((db2lig%atom_vdwtype(atomnum).eq.7).or.
c    &                 (db2lig%atom_vdwtype(atomnum).eq.6)) then
c                      write(*,*) indcord, atomnum, "is an H."
c                      cycle ! go to the next coord 
c                  endif
c                  write(*,*) indcord, "->", atomnum
                   ! loop over the bonds
                   do indbond = 1, db2lig%total_bonds
c                     write(*,*) "coord_index_3", 
c    &                   db2lig%coord_index(3,indcord),
c    &                   ";  conf=", indsc,
c    &                   ";  sets=", indts
                      ! if the atom is in the bond then say the
                      ! conf is in the bond
                      if (db2lig%bond_start(indbond) .eq.
     &                    atomnum) then
c                         write(*,*) "in bond (start): conf =", confnum,
c    &     "bond=",indbond,"atom=",atomnum
                         !bondconflist(indbond,1) = indtc
                         bondconflist(indbond,1) = confnum 
c    &                   db2lig%coord_index(3,indcord)
                      endif
                      !else if (db2lig%bond_end(indbond) .eq. 
                      if (db2lig%bond_end(indbond) .eq. 
     &                         atomnum) then
c                         write(*,*) "in bond (end): conf =", confnum,
c    &     "bond=",indbond,"atom=",atomnum
                         !bondconflist(indbond,2) = indtc
                         bondconflist(indbond,2) = confnum 
c    &                   db2lig%coord_index(3,indcord) 
                      endif
                   enddo ! indbond
              
                 enddo ! indcord
             enddo ! indsc -- for each conf in the set.
             ! finished finding what confs are connected with in the set. 

!            write(*,*) "bond info"
             do indbond = 1, db2lig%total_bonds
                if ( bondconflist(indbond,1).ne.bondconflist(indbond,2))
     &             then
!                       write(*,*) "  ", bondconflist(indbond,1),
!     &                      bondconflist(indbond,2)
                       conf_connect_count(bondconflist(indbond,1)) =
     &                 conf_connect_count(bondconflist(indbond,1)) +1   
                       conf_connect_count(bondconflist(indbond,2)) =
     &                 conf_connect_count(bondconflist(indbond,2)) +1   
                endif
             enddo

            
             do indtc = 1, db2lig%total_confs ! loop all confs
                !if (conf_connect_count(indtc).eq.1) then
                !   write(*,*) "Conf with one connection.  ",
     &          !         conf_connect_count(indtc), indtc 
                !endif
                 if (conf_count(indtc).eq.conf_count(mostfreqconf) )then
                 ! check
! if the the conf with just one connection is also as fequent as the
! most frequent. 
                 !    write(*,*) "replace mostfreqconf = ",indtc
                      mostfreqconf = indtc
                 endif
             enddo
             !write(*,*) "mostfreqconf = ", mostfreqconf 
          !enddo
          !   sort_count = 1
             !db2lig%conf_sorted(1) = mostfreqconf ! set conf_conection equal to itself
             


             ! (1) we need to start with overlaping conf (this is the rigid segment
             ! that is oriented).
             ! (2) find everything that is connected  to the conf1
             !     then have db2lig%conf_conection(conf_k) = conf1
             ! (3) then look for every connected to those confs
             !     and have db2lig%conf_conection point to what they are
             !     connected to.  
             !     db2lig%conf_conection(k) = i 
             !     i is always closer to conf1.
             ! (4) to do this we will use a while loop

! I have a list of connections ! find out what is connected to the
! common segment (conf) not always 1 
! then
! find out what is connected to those nodes 
             !currentlist(1) = 1
             !currentlist(1) = 1
             !currentlist(mostfreqconf) = mostfreqconf
             currentlist(1) = mostfreqconf
             currentint = 1 
             !currentint = mostfreqconf

             ! both seenlist and newlist are reused for multiple
             ! ligand database so we need to cleer the list so confusion
             ! does not happen.  
             !do seenint = 2, db2lig%total_confs
             !do seenint = 1, db2lig%total_confs
             !   seenlist(seenint) = 0
             !   newlist(seenint)  = 0
             !enddo
             !!seenlist(1) = 1
             !seenlist(1) = mostfreqconf
             !newlist(1) =  mostfreqconf
             !
             !seenint = 2

!             if (indts.eq.1) then ! if it the frist segement initalize
                                  ! the seenlist
                do seenint = 1, db2lig%total_confs
                   seenlist(seenint) = 0
                enddo
                seenlist(1) = mostfreqconf
                seenint = 2
!             endif ! otherwise keep appending. 

             do count1 = 1, db2lig%total_confs
                newlist(count1)  = 0
             enddo
             newlist(1) =  mostfreqconf

             count1  = 1 
             count2  = 1 
             newint  = 1

             ! if there is only one confermation 
             if (db2lig%total_confs == 1) then
                 count2 = 100002 ! do not enter loop
             endif

             do while((count1.le.currentint).and.(count2<100000))
! ! look for what is connected to the currentlist entry.    
! the vec allow us
! to look up using the that say segment 201 is connected to 1 or 35 is
! connected to 250.    !       note that, alway, the index points to
! segment closer to segment 1.  
                 if (.not.read_okay) then ! if there is a problem with
                                          !the molecule read in then brake out of the loop amediately
                    exit
                 endif
 
                 do indbond = 1, db2lig%total_bonds
                    if (.not.read_okay) then ! problem with molecule, brake out of loop
                       exit
                    endif
                    if ((bondconflist(indbond,2).eq.bondconflist(
     &                indbond,1)).or.((bondconflist(indbond,2).eq.0)
     &                .or.(bondconflist(indbond,1).eq.0))) then
                        cycle ! go to next bond
                    endif
                    if (bondconflist(indbond,1).eq.
     &                  currentlist(count1)) then 
                        ! make sure we have not already seen this conf
                        flaginlist= inlist(bondconflist(indbond,2),
     &                  seenlist,seenint)
                        if (flaginlist) then
c                           write(*,*) "seen ", bondconflist(indbond,2)
                            cycle ! go to next bond
                        endif
                        newint=newint+1
c                       if (newint.gt.db2lig%total_confs)then
c                         write(*,*) "warning:"
c                         write(*,*) newint, ">", db2lig%total_confs
c                       endif
                        !if (newint.gt.(2*db2lig%total_confs))then
                        if (newint.gt.(10*db2lig%total_confs))then
                          write(*,*) bondconflist(indbond,1),
     &                    bondconflist(indbond,2), "bonds with error"
                          write(*,*) "Error. newlist is not big enough"
                          call doflush(OUTDOCK)
                          read_okay = .false.
                          exit
                          !stop
                        endif
                        newlist(newint) = bondconflist(indbond,2)
                        db2lig%conf_conection(bondconflist(
     &                  indbond,2)) = bondconflist(indbond,1)
c                       write (*,*) bondconflist(indbond,1),
c    &          "is connected to", bondconflist(indbond,2)
c                       write (*,*) db2lig%bond_start(indbond),
c    &                              db2lig%bond_end(indbond)
c                       call doflush(OUTDOCK)

                    else if (bondconflist(indbond,2).eq.
     &                  currentlist(count1)) then
                        ! make sure we have not already seen this conf
c                       if (call inlist(bondconflist(indbond,1),
c    &                  seenlist,seencount )) then
                        flaginlist= inlist(bondconflist(indbond,1),
     &                  seenlist,seenint)
                        if (flaginlist) then
c                           write(*,*) "seen ", bondconflist(indbond,1)
                            cycle ! go to next bond
                        endif
                        newint=newint+1
c                       if (newint.gt.db2lig%total_confs)then
c                         write(*,*) "warning: newint > total_confs"
c                         write(*,*) newint, ">", db2lig%total_confs
c                       endif
                        if (newint.gt.(10*db2lig%total_confs))then
                          write(*,*) bondconflist(indbond,1),
     &                    bondconflist(indbond,2), "bonds with error"
                          write(*,*) "Error. newlist is not big enough"
                          call doflush(OUTDOCK)
                          read_okay = .false.
                          !stop
                          exit
                        endif
                        newlist(newint) = bondconflist(indbond,1)
                        db2lig%conf_conection(bondconflist(
     &                      indbond,1)) = bondconflist(indbond,2)
c                       write (*,*) bondconflist(indbond,2),
c    &          "is connected to", bondconflist(indbond,1)
c                       write (*,*) db2lig%bond_start(indbond),
c    &                              db2lig%bond_end(indbond)
c                       call doflush(OUTDOCK)
                    endif
                 enddo ! loop over bonds
                 if (seenint.gt.db2lig%total_confs)then
                     write(*,*) "Error. seenint > db2lig%total_confs"
                     write(*,*) seenint, db2lig%total_confs
                     call doflush(OUTDOCK)
                     read_okay = .false.
                     !stop
                     exit
                 endif
                 if (count1.gt.(10*db2lig%total_confs))then
                     write(*,*) "Error. count1 > 10*db2lig%total_confs "
                     write(*,*) count1, db2lig%total_confs
                     call doflush(OUTDOCK)
                     read_okay = .false.
                     !stop
                     exit
                 endif

                 flaginlist= inlist(currentlist(count1),
     &           seenlist,seenint)
                 if (.not.flaginlist) then ! check that it is not
                                           ! already in the list
                    seenlist(seenint) = currentlist(count1)
                    seenint = seenint+1
                 endif
c                write (*,*) "seen it list", seenlist(seenint) 
                 count1 = count1+1
                
                 if ((count1.eq.(currentint+1)).and.(newint.gt.0)) then
                    currentint  = newint
                    currentlist = newlist
                    newint = 0
                    count1 = 1
                 endif
                 count2=count2+1 
             enddo ! while loop
                 
c            do indtc = 1, db2lig%total_confs ! loop over all rigid segements (confs)
c                if (db2lig%conf_conection(indtc).ne.indtc) then
c                write (*,*) "conf=",indtc, 'connetted_conf=',
c    &                       db2lig%conf_conection(indtc)
c                endif
c                call doflush(OUTDOCK)
c            enddo
c            stop
c            write (*,*) "look at seenlist list" 
             do indtc = 1, db2lig%total_confs ! loop over all rigid segements (confs)
c               write (*,*) "conf=",indtc, 'seenlist=',
c    &                     seenlist(indtc)
                if (seenlist(indtc).eq.0) then
                   cycle
                endif
                flaginlist= inlist(seenlist(indtc),
     &           db2lig%conf_sorted,count3)
                if (flaginlist) then ! if in sort go to next
                   cycle
                endif
                db2lig%conf_sorted(count3) = seenlist(indtc)
                count3 = count3+1
                call doflush(OUTDOCK)
             enddo
          enddo ! indts -- this is each set

c         ! make a list the has the confs in the order to do the conf scoring.
c         ! we will need to score it in a perticular order to insure
c         ! that the gist is scored properly.  in other words score
c         ! conf 1 frist, then everything connected to the conf, then
c         ! everything connected to those confs.
c         do indtc = 1, db2lig%total_confs ! loop over all rigid segements (confs)
c            !db2lig%conf_sorted(indtc) = indtc ! set conf_conection equal to itself
c            db2lig%conf_sorted(indtc) = 0 ! set conf_conection equal to itself
c         enddo

c         !currentlist(1) = 1
c         currentlist(1) = mostfreqconf
c         currentint = 1
c!          seenlist(1) = 1
c         !db2lig%conf_sorted(1) = 1
c         db2lig%conf_sorted(1) = mostfreqconf
c         seenint = 2
c         count1  = 1
c         newint  = 0
c         count2  = 1

c         ! note that because of how conf_conection is constructed 
c         ! we should never revisit a node.  that is, currentlist should
c         ! never have repeting entry.  
c         do while((count1.le.currentint).and.(count2<100000))
c            !write(*,*) "count1=", count1
c            if (.not.read_okay) then
c                exit 
c            endif 
c            do indtc = 1, db2lig%total_confs ! loop over all rigid segements (confs) no need to look at one again
c               !write(*,*) indtc
c               if (currentlist(count1).eq.
c    &              db2lig%conf_conection(indtc)) then
c            !      write(*,*) currentlist(count1), indtc
c                  db2lig%conf_sorted(seenint) = indtc
c                  newint=newint+1
c                  if (newint.gt.(10*db2lig%total_confs))then
c                      write(*,*) "Error. newlist is not big enough"
c                      call doflush(OUTDOCK)
c                      read_okay = .false.
c                      exit
c                      !stop
c                  endif
c                  newlist(newint)=indtc
c                  seenint = seenint+1
c               endif
c            enddo
c            count1=count1+1
c            if ((count1.eq.(currentint+1)).and.(newint.gt.0)) then
c            !if ((count1.eq.(currentint)).and.(newint.gt.0)) then
c                currentint  = newint
c                currentlist = newlist
c                newint = 0
c                count1 = 1
c            !    write(*,*) "currentint",currentint
c            endif
c            count2=count2+1
c         enddo
c         !write out sorted confs:

c for debuging
c         write (*,*) "look at sorted list" 
c         do indtc = 1, db2lig%total_confs ! loop over all rigid segements (confs)
c            write (*,*) "conf=",indtc, 'sorted_conf=',
c    &                  db2lig%conf_sorted(indtc)
c            call doflush(OUTDOCK)
c         enddo
          endif ! use_gist 

          ! we are already at the End of the hierarchy 
          if (.not. read_okay) then !there was a problem while processing confs 
c            write(OUTDOCK, *) "skipping molecule: ",
c     &          db2lig%refcod, db2lig%protcod
            read_okay = .true. !reset and keep reading next ligand
            line = ' ' !otherwise loop doesn't restart since it sees E
            lines_seen = 0 !really, everything has to be reset
            cycle !skips the rest of the do loop, don't process partial ligand
          endif

          
          exit !break out of do loop
        endif
        if (inputend) then
          call gzclose(fdlig, inputstatus) !close ligand file since we hit EOF
          write(OUTDOCK, *) 'close the file: ', TRIM(options0%ligfil)
          if (.not. options0%sdi) then !no split database, multiple ligand files
            exit !break out of do loop since EOF was hit
          else !using split database index, multiple ligand file
c           read_okay = .false.  !terminate main do loop, return to search,
c                                !a new db2 file may be read in. 
            inputend = .false. ! try to read the next file
            inputstatus = -1 ! this in not zero
            db2c = db2c + 1
            ligc = 0
            !options0%ligfil = ""
            ! SDI stands for split_database_file -- this is used to pass
            ! the program multiple db2 file names to be read in.
            read(SDIFILE, '(a255)', err=9999, end=9999) options0%ligfil !next filename
            !write(OUTDOCK, *) "readin:", TRIM(options0%ligfil)
            call gzopen(fdlig, 'r', options0%ligfil, inputstatus) 
            write(OUTDOCK, *) 'open the file: ', TRIM(options0%ligfil)
            !lines_seen = 0
            !write(OUTDOCK, *) "inputstatus:", inputstatus 
9999        continue ! if read function encounters an error it will come
                     ! here, note inputstatus = -1 so the if statement
                     ! will be entered. 
            if (inputstatus .eq. -1) then !
              write(OUTDOCK, *)
     &              ' we reached the end of the 
     &                split_database_index file'
              write(OUTDOCK, *) "close SDIFILE"
              close(SDIFILE) !close the sdi file
              inputend = .true. ! if problems reading the sdi file
              read_okay = .false.  !terminate main do loop, return to search,
            else if (inputstatus .ne. 0) then !open failed
              write(OUTDOCK, *)
     &             'Error in reading in ligand file:', 
     &             TRIM(options0%ligfil)
              !inputend = .true. ! there is a problem reading in db2 file, so terminate
              !write(OUTDOCK, *) "close SDIFILE"
              !close(SDIFILE) !close the sdi file
              write(OUTDOCK, *) "Full program STOP"
              stop
            endif
          endif
        endif
      enddo
      !enddo
      return
      end subroutine ligread2

c     ! this is a function that look to see if some integer is in the list
c     ! of integers: if it is it returns true, if it is not it returns false
c     logical function inlist(val,list,len) result(flag)
      function inlist(val,list,len) result(flag)

      implicit none

      integer, intent(in) :: val,len
      integer, intent(in), allocatable, dimension(:) :: list
      integer             :: i 
      logical             :: flag
      
      flag = .false.
      do i = 1, len
          if (list(i).eq.val) then
            flag = .true.
            exit ! leave list
          end if
      end do  
      end function inlist


      function tokenize(line,token) result(length)

      implicit none

      character (len=80), intent(in) :: line
      character (len=5), intent(in) :: token
      !integer, intent(inout) :: length
      integer               :: length
      logical               :: flag1,flag2
      integer               :: n,i
      !integer              :: pos1,pos2,n,i

      !pos1 = 1
      n = 0
      flag1 = .false.
      flag2 = .false.
      !len = size(line)
      do i=1, len(line)
          !pos2 = INDEX(line(pos1:),token)
          flag1 = flag2
          if (line(i:i).ne.token) then
            flag2 = .false.
          else
            flag2 = .true.
          endif
          !if ((flag1.eq..true.).and.(flag2.eq..false.)) then
          if ((flag1).and.(.not.flag2)) then
            n = n + 1
          endif
      end do
          !if (pos2 .eq. 0) then
            !n = n + 1
            !exit ! leave list
          !end if
          !n = n + 1
          !pos1 = pos2+pos1
      !end do
      length = n
      end function tokenize
