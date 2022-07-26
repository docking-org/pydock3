c----------------------------------------------------------------------
c----------------------------------------------------------------------
c read procedure for ligand from mol2 for rescoring.  Trent E Balius 
c modifed the read procedure for new ligand database format db2 by rgc
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine ligread_mol2(db2lig, options0, 
     &    inputend, fdlig, fdsol, fdvdw)

      use db2type
      use optionstype
      use mol2formats
      use filenums
      use spheres
      use transfm_conf 

      implicit none

      type(db2), intent(inout) :: db2lig
      type(options), intent(inout) :: options0
      logical, intent(inout) :: inputend !EOF indicator; 1 if EOF reached, returned to call
c      integer (kind=8), intent(inout) :: fdlig !ligand file handle
c      integer (kind=8), intent(inout) :: fdsol !ligand file handle
c      integer (kind=8), intent(inout) :: fdvdw !ligand file handle
      integer (kind=8), intent(inout) :: fdlig, fdsol, fdvdw

      character (len=80) line !used to store one read line at a time
      character (len=78) temparbitrary !used to store arbitrary data temporarily
      integer inputstatus !used in gzread functions to get status of read
      character control !first character used to control program flow
      integer lines_seen !records how many seen within this segment of the file
      logical read_okay !used to keep track of a good ligand read. read until 
         !we find a good one or the file is empty
      logical  read_mole, read_atom, read_bond, read_stop ! this is denote were in mol2 file we are
      integer lines_count, current, temp_childlen !used in group/conf trees
      !integer children_seen, current_test, line_test !used in group/conf trees
      integer data_count !used in group/conf trees
      integer temp_start, temp_end !used in reading in confs
      integer tempbond, tempatom, tempotheratom !for valency calculation
      !integer broken_temp, hydro_temp !used in reading in sets
      !real energy_temp !used in reading in sets
      integer molcount,atomcount, bondcount, tempint!,molnum
      character (len=16) molid,molnum
      character (len=20) temp,temp2
c      character (len=80) vdwfile
c      character (len=80) amsolfile
c      vdwfile   = "VDW.txt"
c      amsolfile = "AMSOL.txt"
      
c      write (*,*) "vdwfile = ", vdwfile
c      write (*,*) "amsolfile = ", amsolfile
c      write (OUTDOCK,*) "In function ligread_mol2"
c      call doflush(OUTDOCK)

       lines_seen = 0
       line = '                   ' !placeholder that forces new line to always be read in
       read_okay      = .true.
       read_mole      = .false.
       read_stop      = .false.
       read_atom      = .false.
       read_bond      = .false.
       molcount       = 0
       atomcount      = 1
       bondcount      = 1

c     set some important values fro output in the header. 

       db2lig%refcod         = ""
       db2lig%protcod        = ""
       db2lig%smiles         = ""
       db2lig%molname        = ""
       options0%ligfil       = ""
       db2lig%total_charge   = 0.0
       db2lig%total_surfarea = 100
       db2lig%total_confs    = 1 
       db2lig%total_sets     = 1
       db2lig%total_clusters = 1 
       db2lig%total_heavy_atoms = 0

c      write (OUTDOCK,*) "In function ligread_mol2 (2) "
c      call doflush(OUTDOCK)
        

c     do while (read_okay)
        do while (.NOT.read_stop) !detects end of molecule
          !fdlig is always the ligand file, already opened for us already
          !gzread is used since files are always gzipped
          !any errors are indicated by inputstatus
          call gzread(fdlig, line, 80, inputstatus)
          if (inputstatus .le. 0) then !means the end of file has been reached
            inputend = .true. 
            exit !break out of do loop
          endif

c         write (OUTDOCK,*) line 

          if (line(1:1) .ne. " " .and. line(1:1) .ne. "@" .and.
     & line(1:1) .ne. '#') then 
            !write (OUTDOCK,*) line(1:1) 
c           write (OUTDOCK,*) "cycle:", line 
            call doflush(OUTDOCK)
            cycle
          endif
          
c         if (line(1:1) .ne. " " .and. line(1:1) .ne. "@" ) then 
c             cycle
c         else if (line(1:10) .eq. "##########") then
          if (line(1:10) .eq. "##########") then
              if (read_bond) then ! if ## is incountered after bonds
                                  ! then exit the loop and return the
                                  ! molecule
                 ! break to loop 
                 read_stop = .true. ! if the next mol is incountered then stop reading
                 !exit
              endif 
              cycle
          else if (line(1:18) .eq. "@<TRIPOS>MOLECULE") then
              ! loop will terminate when the next molecule is reached. 
              !if (read_mole ) then 
              molcount  = 0
              read_mole = .true. 
              if (read_bond ) then 
                  !read_stop = .true. ! if this is the next mol is incountered then stop reading
                  !read_mole = .false. 
                  !inputend = .true.  ! this is so that only one mol2 is read in 
                  write (OUTDOCK,*) 'Error: mol2 must have "##########"
     & Between molecules.'
              endif
              read_atom = .false.
              read_bond = .false.
          else if (line(1:13) .eq. "@<TRIPOS>ATOM") then
              atomcount = 1
              read_mole = .false. 
              read_atom = .true.
              read_bond = .false.
              cycle !continue next line in file
          else if (line(1:13) .eq. "@<TRIPOS>BOND") then
              bondcount = 1
              read_mole = .false. 
              read_bond = .true.
              read_atom = .false.
              cycle !continue to next line in file
          endif

          if (read_mole .and. molcount.eq.1) then ! read in atom info
              read(line, MOL2MOLE1) molid, molnum
              db2lig%refcod = molid
              !call REMOVE(molid,db2lig%refcod," ") !molid
              db2lig%protcod = molnum
c             write(OUTDOCK,*) molid, molnum
          else if (read_mole .and. molcount.eq.2) then ! read in number of atoms and bonds
              read(line, MOL2MOLE2) db2lig%total_atoms,
     &                               db2lig%total_bonds,
     &                               tempint,tempint,tempint
c             write(OUTDOCK,*) "INFO FOR ALLOCATION",db2lig%total_atoms,
c    &                   db2lig%total_bonds,
c    &                   tempint,tempint,tempint
c             call doflush(OUTDOCK)
              db2lig%total_coords = db2lig%total_atoms * 3
c             rigid_heavy_count   = db2lig%total_atoms ! this might need to be changed 
              db2lig%mlines       = db2lig%total_atoms +
     &                              db2lig%total_bonds  ! this need to be changed


              !write (*,*) "In function ligread_mol2 (3) "
              !call doflush(OUTDOCK)
              !want to dynamically allocate storage now
              call allocate_db2(db2lig, OUTDOCK, db2lig%total_sets, 
     &            db2lig%total_coords, db2lig%total_atoms, 
     &            db2lig%total_bonds,
     &            db2lig%total_confs, db2lig%mlines, 
     &            db2lig%total_clusters, 
     &            options0%dockovalent, options0%nsav)

              !write (*,*) "In function ligread_mol2 (4) "
              !call doflush(OUTDOCK)
              call read_vdw_file(db2lig,db2lig%total_atoms,fdvdw)
              !write (*,*) "In function ligread_mol2 (5) "
              !call doflush(OUTDOCK)
              call read_amsol_file(db2lig,db2lig%total_atoms,fdsol)
              !write (*,*) "In function ligread_mol2 (6) "
              !call doflush(OUTDOCK)
              tempatom = 1 ! vdwtype index starts at 1 not zero. 
              ! calculate heavy atom count. 
              do while (tempatom .le. db2lig%total_atoms)
                if ((db2lig%atom_vdwtype(tempatom) .ne. 6) .and.
     &              (db2lig%atom_vdwtype(tempatom) .ne. 7)) then !tempatom is heavy
                    db2lig%total_heavy_atoms = 
     &                 db2lig%total_heavy_atoms + 1 ! count heavy atoms
                endif
                tempatom = tempatom + 1
              enddo


              db2lig%set_energy(1)       = 0.0 
              db2lig%set_conf_len(1)     = 1
              db2lig%set_conf(1, 1)      = 1
              db2lig%mlines = 5
              db2lig%arbitrary(1)        = "999"
              ! there is only one conformation 
              ! set as all atoms. 
              db2lig%conf_coord(1, 1)    = 1 !atom index of start 
              db2lig%conf_coord(2, 1)    = db2lig%total_atoms!same but end (2)
          endif

          if (read_stop) then ! read in number of atoms and bonds


              tempatom = 0
              !intialize the valency
              ! call a function 
              ! move this stuff to that function
              do while (tempatom .le. db2lig%total_atoms)
                db2lig%atom_heavy_valence(tempatom) = 0. !reset to zero
                tempatom = tempatom + 1
              enddo
              !now go through bonds. whenever one atom is heavy, increment
              !the other
              !only do this if we want to use the valence later
              !write (6,*) 'xxxrvp', options0%recdes_valence_penalty
              if (options0%recdes_valence_penalty) then
                tempbond = 1
                do while (tempbond .le. db2lig%total_bonds)
                  tempatom = db2lig%bond_start(tempbond)
                  tempotheratom = db2lig%bond_end(tempbond)
                  if ((db2lig%atom_vdwtype(tempatom) .ne. 6) .and.
     &                (db2lig%atom_vdwtype(tempatom) .ne. 7)) then !tempatom is heavy
                    !so tempotheratom gets incremented
                    db2lig%atom_heavy_valence(tempotheratom) =
     &                  db2lig%atom_heavy_valence(tempotheratom) + 1
                  endif
                  !now go the other way
                  if ((db2lig%atom_vdwtype(tempotheratom) .ne. 6) .and.
     &                (db2lig%atom_vdwtype(tempotheratom) .ne. 7)) then !heavy
                    !so tempatom gets incremented
                    db2lig%atom_heavy_valence(tempatom) =
     &                  db2lig%atom_heavy_valence(tempatom) + 1
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
                !write(OUTDOCK, *) 'hv', db2lig%atom_heavy_valence !debug
                !output
              endif
              !read_stop = .false. 
          endif

          if (read_mole) then
              molcount = molcount +1
          endif

          if (read_atom) then ! read in atom information
              if (atomcount .gt. db2lig%total_atoms) then
c                 write(*,*) "Cycle. . . atom"
                 cycle
              endif
c      1 N1         46.9071    43.2086    26.9254 N.pl3      1          -0.6036
              read(line, MOL2ATOM) db2lig%atom_num(atomcount), 
     &              db2lig%atom_name(atomcount),
     &              db2lig%coords(1, atomcount),
     &              db2lig%coords(2, atomcount),
     &              db2lig%coords(3, atomcount),
     &              db2lig%atom_type(atomcount),
     &              temp,
c    &              temp2,
     &              db2lig%atom_charge(atomcount)
              !write(OUTDOCK, * ) temp2
              !call doflush(OUTDOCK)
              db2lig%coord_index(1, atomcount) = atomcount ! coord number
              db2lig%coord_index(2, atomcount) = atomcount ! atom number
              db2lig%coord_index(3, atomcount) = 1 !  conformation number


              db2lig%total_charge  = db2lig%total_charge + 
     &                               db2lig%atom_charge(atomcount)

c for rescoring transfm_coords will equal the orginal coords
              db2lig%transfm_coords(1, atomcount) = 
     %                db2lig%coords(1, atomcount)
              db2lig%transfm_coords(2, atomcount) = 
     %                db2lig%coords(2, atomcount)
              db2lig%transfm_coords(3, atomcount) = 
     %                db2lig%coords(3, atomcount)

c             write(OUTDOCK,*) db2lig%atom_type(atomcount),
c    &              " ",
c    &              db2lig%atom_name(atomcount),
c    &              atomcount,
c    &              db2lig%coords(1,atomcount),
c    &              db2lig%coords(2, atomcount),
c    &              db2lig%coords(3, atomcount),
c    &              db2lig%atom_charge(atomcount)
c             call doflush(OUTDOCK)
              atomcount = atomcount + 1
          endif

          if (read_bond) then
              if (bondcount .gt. db2lig%total_bonds) then
c                write(*,*) "Cycle. . . bonds"
                 call doflush(OUTDOCK)
                 cycle
              endif
              read(line, MOL2BOND) db2lig%bond_num(bondcount),
     &            db2lig%bond_start(bondcount),
     &            db2lig%bond_end(bondcount),
     &            db2lig%bond_type(bondcount)

c             write(OUTDOCK,*) bondcount, db2lig%bond_type(bondcount),
c    &            db2lig%bond_num(bondcount),
c    &            db2lig%bond_start(bondcount),
c    &            db2lig%bond_end(bondcount)

              bondcount = bondcount + 1
          endif

        enddo ! read_stop 
        if (inputend) then
c         write(*,*) "closing files."
          call gzclose(fdlig, inputstatus) !close ligand file since we hit EOF
          call gzclose(fdsol, inputstatus) !close ligand sol file since we hit EOF
          call gzclose(fdvdw, inputstatus) !close ligand vdw file since we hit EOF
        endif
      end
c end subroten ligreadmol2

c function to remove a char
      subroutine REMOVE(ostr,nstr,rchar)
          character (len=20), intent(in) :: ostr
          character (len=20), intent(out) :: nstr
          character, intent(in)  :: rchar ! replace char
          integer i
          integer ls1
          integer ls2

          ls1 = len_trim(ostr)
          ls2 = 0
          do i = 1,ls1
             if(ostr(i:i).ne.rchar) then
                ls2 = ls2 + 1
                nstr(ls2:ls2) = ostr(i:i)
             endif
          enddo
          return 
      end


c     subroutine cal_valancy(db2lig,atom_num,fdvdw)
c     ! this need to be added
c     ! this is if we want to add in code to calculate vdw type on the
c     ! fly and not read them in from a file.  
c     end ! cal_valancy

c     subroutine read_vdw_file(file_name,atom_vdwtype,atom_num)
      subroutine read_vdw_file(db2lig,atom_num,fdvdw)
          use db2type
          use filenums
c         character (len=80) , intent(in) :: file_name
          integer, intent(in) :: atom_num
          type(db2), intent(inout) :: db2lig
c         real, dimension(atom_num), intent(inout) :: atom_vdwtype
          character (len=80):: tempstr
          integer tempint, stat, i
          integer inputstatus !used in gzread functions to get status of read
          integer (kind=8), intent(inout) :: fdvdw !ligand file handle
          character (len=80) :: line = ''

c         integer, parameter :: VDWFILE = 83
c         write(*,*) file_name
c         write(OUTDOCK, *) "In read_vdw_file"
c         call doflush(OUTDOCK)

c         read in the name of the molecule
          call gzread(fdvdw, line, 80, inputstatus)
          if (line(1:1) .ne. '#') then
              write(OUTDOCK, *) line
              write(OUTDOCK, *) "Error: looks like vdw parameters
     & have too meny lines, may not have name line. expects ## name."
              stop
          end if
c         write(OUTDOCK,*) line
c         call doflush(OUTDOCK)

c         open(VDWFILE, file=file_name, status='old', 
c    &        action='read')

          do i = 1,atom_num ! note index starts at 1 
             !(A6,F7.3,F7.3,F7.3)
             !read(VDWFILE,"(I2 A2 A5 I2)",IOSTAT=stat)  tempint,
c            read(fdvdw,*,IOSTAT=stat)  tempint,
c    &                      tempstr,tempstr,db2lig%atom_vdwtype(i) 
             call gzread(fdvdw, line , 80, inputstatus)
c            write(OUTDOCK, *) line

             if (line(1:1) .eq. '#') then
                 write(OUTDOCK, *) line
                 write(OUTDOCK, *) "Error: looks like vdw parameters
     & have too few lines"
                 stop
             end if
                 
             read(line,*)  tempint,
     &                      tempstr,tempstr,db2lig%atom_vdwtype(i) 
             
             if (stat .lt. 0) then
                 WRITE(OUTDOCK,*) "ERROR: rescore vdw type file is not
     & right may be too long"
                 exit
             endif
c            WRITE(*,*) tempint, tempstr, db2lig%atom_vdwtype(i)
          enddo
c         close(VDWFILE)

      end ! end function read_vdw_file

c     subroutine read_amsol_file(file_name,atom_polsolv,atom_apolsolv,
c    &                        atom_totalsolv, atom_surfarea, atom_num)
      subroutine read_amsol_file(db2lig,atom_num,fdsol)
          use db2type
          use filenums
c         character (len=80) , intent(in) :: file_name
          integer, intent(in) :: atom_num
          type(db2), intent(inout) :: db2lig
c         real, dimension(atom_num), intent(inout) :: atom_polsolv
c         real, dimension(atom_num), intent(inout) :: atom_apolsolv
c         real, dimension(atom_num), intent(inout) :: atom_totalsolv
c         real, dimension(atom_num), intent(inout) :: atom_surfarea
          integer tempint, stat, i
          character (len=8) :: tempstr
          real :: tempreal
          integer, parameter :: AMSOLFILE = 84
          integer inputstatus !used in gzread functions to get status of read
          integer (kind=8), intent(inout) :: fdsol !ligand file handle

          character (len=80) :: line = ''
          character (len=*), parameter :: AMSOLFORMAT =
     &    '(i2,1x,a3,f10.3,f11.3,f11.3,f10.3)'

c         write(OUTDOCK,*) "In read_amsol_file"
c         call doflush(OUTDOCK)
c         write(*,*) file_name

c         open(AMSOLFILE, file=file_name, status='old',
c    &        action='read')


c 1  N1     0.380     -1.210     -0.830     6.230
c 2  C1     0.110      3.300      3.410    12.970
c 3  N2    -3.510     -2.200     -5.710    11.270
c 4  C2    -0.130      0.790      0.660     7.230
c 5  C3    -0.530      0.470     -0.060    10.050
c 6  C4    -0.200      0.470      0.270     9.960
c 7  C5     0.000      0.470      0.470    10.030
c 8  C6    -0.070      0.470      0.400    10.070
c 9  C7     0.030      0.780      0.810     7.310
c10  H1    -2.230     -0.590     -2.820     8.960
c11  H2    -1.610      0.280     -1.330     8.060
c12  H3     0.590      0.280      0.870     8.060
c13  H4    -0.220      0.280      0.060     8.060
c14  H5    -0.600      0.280     -0.320     8.060
c15  H6    -0.430      0.280     -0.150     8.060

c         read in the name of the molecule
          call gzread(fdsol, line, 80, inputstatus)
          if (line(1:1) .ne. '#') then
              write(OUTDOCK, *) line
              write(OUTDOCK, *) "Error: looks like amsol parameters
     & have too meny lines, may not have name line. expects ## name."
              stop
          end if

          if (inputstatus .le. 0) then !means the end of file has been reached
               write(OUTDOCK, *) 'Error: at end of the file.'
               stop
          endif
c         write(OUTDOCK,*) line,'\n'

          do i = 1,atom_num
             call gzread(fdsol, line, 80, inputstatus)
             if (line(1:1) .eq. '#') then
                 write(OUTDOCK, *) line
                 write(OUTDOCK, *) "Error: looks like amsol parameters
     & have too few lines"
                 stop
             end if
c            read(line,AMSOLFORMAT) tempint,
             read(line,*) tempint,
     &                      tempstr, db2lig%atom_polsolv(i),
c    &                      atom_apolsolv(i),atom_totalsolv(i),tempreal
     &                      db2lig%atom_apolsolv(i),
     &                      db2lig%atom_totalsolv(i),
     &                      db2lig%atom_surfarea(i)
             if (inputstatus .le. 0) then
                 WRITE(*,*) "ERROR: rescore vdw type file is not right
     & may be  too long"
                 exit
             endif
c            WRITE(OUTDOCK,*) tempint, tempstr, db2lig%atom_polsolv(i),
c    &                  db2lig%atom_apolsolv(i),
c    &                  db2lig%atom_totalsolv(i),
c    &                  db2lig%atom_surfarea(i)  
c            call doflush(OUTDOCK)
          enddo

c         close(AMSOLFILE)

      end

cc this function is to idenify the dock atom type from the sybyl type
cc this function is not working yet. 
c      subroutine get_vdw_type( numofa, numofb, atom_num, atom_type, 
c     &                         bond_num, bond_start, bond_end, vdw_list)
c
c
c        integer, intent(in) :: numofa
c        integer, intent(in) :: numofb
c        integer, dimension(numofa), intent(in) :: atom_num
c        character (len=5), dimension(numofa), intent(in) :: atom_type
c        integer, dimension(numofb), intent(in) :: bond_num !bond number
c        integer, dimension(numofb), intent(in) :: bond_start !bond atom (arbitrary, left side)
c        integer, dimension(numofb), intent(in) :: bond_end !bond atom (arbitrary, right side)
c
c        integer, dimension(numofa), intent(in) :: vdw_list
c        integer, dimension(numofa) :: at_bo_count  ! atom bond count : cardinality of the atom
c        integer, dimension(numofa,5) :: atom_atlist ! atom atom list : the atoms connected to atom i, a max of 5 atoms
c        integer i, j
c
c        write(*,*) "In get_vdw_type "
c        write(*,*) numofa, numofb
c
c        do i = 1, numofa
c           at_bo_count(i) = 0
cc          write(*,*) i, at_bo_count(i)
c        end do
c
cc       do i = 1, numofa
cc          write(*,*) i, atom_type(i), atom_num(i)
cc       end do
c
c 
cc       write(*,*) "In get_vdw_type "
c
c        ! loop over the bonds 
c        ! find the which atoms are bound to atom i 
c        ! count how meny atoms are bound to atom i
c        do i = 1, numofb
cc          write(*,*) i 
cc          write(*,*) "BOND INFO", bond_num(i), bond_start(i), 
cc    &                bond_end(i)
c           if (at_bo_count(bond_start(i)).gt.5 .or.
c     &         at_bo_count(bond_end(i)).gt.5) then
c               write(*,*) "ERROR"
c               stop
c           endif
c           if (bond_start(i).lt.(numofa+1)) then
c              at_bo_count(bond_start(i)) =at_bo_count(bond_start(i))+ 1
c              atom_atlist(bond_start(i),at_bo_count(bond_start(i))) = 
c     &            bond_end(i)
c           else 
c              write(*,*) "ERROR"
c              stop
c           endif
c           if (bond_end(i).lt.(numofa+1)) then
c              at_bo_count(bond_end(i)) = at_bo_count(bond_end(i)) + 1
c              atom_atlist(bond_end(i),at_bo_count(bond_end(i))) =  
c     &            bond_start(i)
c           else
c              write(*,*) "ERROR"
c              stop
c           endif
c        end do
c
c        do i = 1, numofa
cc          write(*,*) "atom ", i, atom_type(i), "has", at_bo_count(i),
cc    &                "bonds"
c           if (atom_type(i).eq."H") then
c               if (at_bo_count(i).ne.1) then
c                  write(*,*) "ERROR: H should only have on bond"
c               end if
c               !if (atom_type(atom_atlist(i,1)).eq."C.ar")
c           end if
c           do j = 1, at_bo_count(i) 
c              write(*,*) "   ", atom_atlist(i,j),
c     &          atom_type(atom_atlist(i,j))
c           end do
c           
c        end do
c      end 

c  1        888.79        24.81    sp2 and sp C
c  2       1586.37        35.05    CH3 (united atom)
c  3       1128.12        27.96    CH2 (united atom)
c  4        769.72        21.49    CH (united atom)
c  5        533.20        16.16    sp3 C
c  6          0.37         0.31    H on polar atom
c  7         85.37         4.13    H on C
c  8        735.31        24.25    sp2 and sp N
c  9        725.70        20.26    quaternary sp3 N
c 10        888.79        24.81    sp3 N
c 11        480.19        20.72    sp2 O
c 12        500.18        19.68    sp3 O
c 13       2454.77        46.86    P
c 14       1831.79        40.48    S
c 15        251.02        11.92    F
c 16       2194.13        46.37    Cl
c 17       3885.92        66.31    Br
c 18       6817.37        92.86    I
c 19        235.24        10.15    Na+ (unhydrated), K+
c 20         51.92         5.73    Mg++, Li+, Al+++, M++ (except Ca++)
c 21        339.55        14.65    Ca++
c 22       1971.47        37.63    Cl- (unhydrated)
c 23        762.07        24.38    Lennard-Jones water particle
c 24       3885.92        66.31    Si (same numbers as Br)
c 25          0.00         0.00    Du/LP (same numbers as H on polar atom)
c 26         53.67         7.33    Zn2+ (Austin)


c     ! is_polar()
c     subroutine is_polar(string) 

c        H     !hydrogen
c        H.spc !hydrogen in Single Point Charge (SPC) water model
c        H.t3p !hydrogen in Transferable intermolecular Potential (TIP3P) water model

c        O.3   !oxygen sp3
c        O.2   !oxygen sp2
c        O.co2 !oxygen in carboxylate and phosphate groups
c        O.spc !oxygen in Single Point Charge (SPC) water model
c        O.t3p !oxygen in Transferable Intermolecular Potential (TIP3P) water model

c        N.1   !nitrogen sp
c        N.2   !nitrogen sp2
c        N.3   !nitrogen sp3
c        N.4   !nitrogen sp3 positively charged

c        N.am  !nitrogen amide
c        N.ar  !nitrogen aromatic
c        N.pl3 !nitrogen trigonal planar

c        S.2 !sulfur sp2
c        S.3 !sulfur sp3
c        S.O !sulfoxide sulfur
c        S.O2 !sulfone sulfur
c        
c        C.3   !carbon sp3
c        C.2   !carbon sp2
c        C.1   !carbon sp
c        C.ar  !carbon aromatic
c        C.cat !carbocation (C+) used only in a guadinium group

c        LP    !lone pair
c        Du    !dummy atom
c        Du.C  !dummy carbon
c        Any   !any atom

c        Hal   !halogen
c        F     !fluorine
c        Cl    !chlorine
c        I     !iodine
c        Br    !bromine

c        Het   !heteroatom = N, O, S, P
c        Hev   !heavy atom (non hydrogen)

c        Li    !lithium
c        Se    !selenium
c        Mo    !molybdenum
c        Na    !sodium
c        Mg    !magnesium
c        Al    !aluminum
c        Si    !silicon
c        K     !potassium
c        Ca    !calcium
c        Cr.th !chromium (tetrahedral)
c        Cr.oh !chromium (octahedral)
c        Co.oh ! cobalt (octahedral)
c        Mn    !manganese
c        Fe    !iron 
c        P.3   !phosphorous sp3
c        Cu    !copper
c        Zn    !zinc
c        Sn    !tin 
c    



