#!/usr/bin/perl -w
# makespheres1.cli.pl
#actually able to specify inputs/outputs in this version
# makespheres1.cli.pl sph/match.sph sph/sph grids/rec.crg sph/match1.sph 40
# 40 is the minimum number of spheres in output. 
#  set to lower for smaller sites. higher for bigger sites. seriously.
# This is a script to make spheres for docking in an automated non-biased method
# Make sure the INDOCK file uses the output file sph/match1.sph for delphi spheres
# Input files are sph/match.sph, sph/sph from sphgen, and grids/rec.pdb
# Usage: makespheres2.pl       (must have subdirectories sph and grids)
# NB: watch out for "off-by-one bug" (OBOB) for line numbers; most arrays start with 0
##########################################################################################
# Pseudocode:
#   Read in sph/match.sph file (derived from xtal-lig.mol2) and determine center of coordinates
#   If user indicated that ligand is crystallographic, then flag spheres for output
#   Read in file sph/sph cluster 0
#   Remove spheres farther than 10 angstroms from center of coords of ligand (match.sph)
#   Read in grids/rec.crg (receptor atom coordinates file, including H's for polar interactions)
#   Remove spheres farther than 7 angstroms from any receptor atom
#   Remove spheres closer than 1.5 angstroms from any receptor atom (including crystallographic spheres)
#   Flag spheres according to proximity and directionality to polar receptor atoms
#      by determining receptor atom type, calculate theta to basetom, and find best valid theta closest to mean
#      if multiple polar receptor atoms are near a single sphere
#      valid angle is defined by the reference: http://prowl.rockefeller.edu/aainfo/hbonds.html
#      lower and upper angle cutoff values are at 95% level (1.96 times standard deviation from the mean)
#   Remove spheres that have contacts with other spheres too close (less than 1.5 A) -- nonbiased method, preference to polar
#   Remove spheres that are <0.8 angstroms apart in 3D space (keep 1 of the 2 spheres), preference to polar
#   Prune spheres to delete any disconnected islands (connectivity within 2.5 A)
#    and if too few spheres (<40), repeat continuity with value successively increased by 0.5 until continuity >4.0
#   If too many spheres (>120 for delphi calculation): with 25% preference to polar from averageweight
#     if crystallographic then remove spheres farthest from:
#        any ligand atom at 1-weight
#        center of ligand at 1/8-weight
#        center of spheres at 1/8-weight
#     if not crystallographic then remove spheres farthest from:
#        center of receptor at 1-weight
#        center of ligand at 1-weight
#        center of spheres at 1-weight
#   Print out results to file sph/match1.sph
##########################################################################################
# Flagging system: each coordinate line has (x,y,z) in array positions [1], [2], [3]
# The sphere flags in array position [6] as follows:
#   0 is a plain sphere that can be deleted
#   1 is a crystallographic sphere (should not be deleted ever, until the very end if >120 spheres)
#   2 is a sphere that is continuous with the center of sphere distribution and any crystallographic spheres
#   3 is a crystallographic sphere that happens to be the center sphere as well
# The sphere flags in array position [7] as follows:
#   -1 near neither polar nor nonpolar receptor atoms (initial value)
#   -2 near nonpolar receptor atom and not near a polar receptor atom
#   <=180 is theta from sphere to polar receptor atom and base atom
# The sphere value in array position [8] is the weighted distance to sort by
#   in eliminating spheres when >120 exist at the end
# The sphere value in array position [9] is the difference from the mean for the angle
#   which is unsed in calculating which sphere is better than another (closer to its mean)
##########################################################################################
# Current coloring system:
#    positive                       (1)   e.g., sphere near 2 Asp/Glu
#    negative                       (2)   e.g., sphere near 2 Arg/Lys
#    acceptor                       (3)   e.g., sphere near backbone HN, Asn/Gln HN, Ser/Thr/Tyr HO
#    donor                          (4)   e.g., sphere near backbone O, Asn/Gln O
#    ester_o                        (5)
#    amide_o                        (6)
#    neutral                        (7)   e.g., sphere near nonpolar C only (no polar atoms)
#    not_neutral                    (8)   e.g., sphere near both positive and negative, pos_don and neg_acc
#    positive_or_donor              (9)   e.g., sphere near Asp/Glu O 
#    negative_or_acceptor           (10)  e.g., sphere near Arg/Lys HN
#    neutral_or_acceptor_or_donor   (11)  e.g., sphere near S-related atoms (Cys, Cyx, Met S and HG)
#    donacc                         (12)  e.g., sphere near both donor & acceptor receptor atoms
##########################################################################################
# TWEAKS: arbitrary parameters can be modified for custom sphere generation
#   $M is margin to keep spheres from center of ligand coordinates (DEFAULT is 10 angstroms)
#   $R is margin to keep spheres from receptor atoms (DEFAULT is 7 angstroms)
#   $Rclose is margin to keep spheres away from receptor atoms (DEFAULT is 1.2 angstroms)
#   $polardist is the polar distance cutoff between sphere and heavy receptor atom that is polar (DEFAULT is 3.3 angstroms)
#   $Hpolardist is the polar distance cutoff between sphere and H receptor atom that is polar (DEFAULT is 2.5 angstroms)
#   $nonpolardist is the nonpolar distance cutoff between sphere and receptor atom that is nonpolar (DEFAULT is 4.5 angstroms)
#   $gridsize is space between spheres to keep (DEFAULT is 1.5 angstroms)
#   $tooclose is distance for second pass to eliminate spheres too close to each other (DEFAULT is 0.8 angstroms)
#   $continuity is distance between spheres to be considered continuous (DEFAULT is 3.0 angstroms)
#   $contincrement is the increment value to increase continuity if fewer than $minspheres found (DEFAULT is 0.5 angstroms)
#   $maxcontinuity is the maximum permitted value continuity can be when incrementally increased (DEFAULT is 4.5 angstroms)
#   $numrepeats is how many iterations to search for continuity of points (DEFAULT is 10 times) this is more than enough
#   $minspheres is arbitrary value of 40 used to keep continuity finding enough spheres, especially open/flat/groove sites
#   $maxspheres is arbitrary value of 120 used near end of script to remove too many spheres for delphi 
#   $polarbenefit is percent benefit for polar sphere from averageweight at end of sphere selection (DEFAULT is -0.25)
##########################################################################################
# Written by Austin N. Kirschner on February 3, 2003 (2-3-2003), and updated in March 2003 and July-August 2003

# modified by Trent Balius on 2014/02/10 to use reduce atom names.  

# CONSTANTS defined here
$M = 12;              # Margin value to keep spheres close to center of ligand
$R = 7;               # Margin value to keep spheres close to receptor atoms
$Rclose = 1.2;        # Margin value to keep spheres away from receptor atoms
$polardist = 3.3;     # Polar distance cuttoff between sphere and heavy receptor atom that is polar
$Hpolardist = 2.5;    # Polar distance cuttoff between sphere and H receptor atom that is polar
$nonpolardist = 4.5;  # Nonpolar distance cuttoff between sphere and carbon receptor atom
$gridsize = 1.5;      # Typical distance between sphere points, i.e., only one point kept in each grid unit area
$tooclose = 0.8;      # Second pass elimination of spheres that are too close, after nonbiased grid method
$continuity = 3.0;    # How far apart points can be and still considered continuous for delphi 
$contincrement = 0.5; # Increment to increase continuity if fewer than $minspheres found based on continiuty
$maxcontinuity = 4.5; # Maxmimum value continuity is permitted to be when it needs to be increased
$numrepeats = 10;     # How many times point-continuity should be checked to grow spheres to keep
$minspheres = $ARGV[4];     # Minimum number of spheres to have for output , niu
$maxspheres = 120;     # Maximum number of spheres to have for delphi input, ARBITRARY based on DOCK, niu
$polarbenefit =-0.25; # Percent benefit from averageweight for final selection of spheres

##################### READ LIGAND ATOMS ##################################################
# Read in match.sph coordinates
#$LIG = "sph/match.sph";
$LIG = $ARGV[0]; #now able to specify on command line
open (LIG, "<$LIG") or die "My error, cannot open LIG: $!\n";
while (<LIG>) {                     # read in LIG lines
        @ligarray = split;
        push @ligitems, [ @ligarray ];
}
close (LIG);
for ($x=1; $x<=$#ligitems; $x++) {   # start of coordinates line
      push @ligxyzcoords, [ @{ $ligitems[$x] } [ 1..3 ] ];
}
print "There are ".($#ligxyzcoords+1)." ligand heavy atoms\n";
if (($#ligxyzcoords+1)==0) { print "ERROR: There are 0 spheres\n"; exit; }
# Determine center of mass of match.sph
$ligavgX = $ligxyzcoords[0][0];
$ligavgY = $ligxyzcoords[0][1];
$ligavgZ = $ligxyzcoords[0][2];
for ($line = 1; $line <= $#ligxyzcoords; $line++) {
   $ligavgX = $ligavgX + $ligxyzcoords[$line][0];
   $ligavgY = $ligavgY + $ligxyzcoords[$line][1];
   $ligavgZ = $ligavgZ + $ligxyzcoords[$line][2];
}
$ligavgX = $ligavgX/($#ligxyzcoords+1);
$ligavgY = $ligavgY/($#ligxyzcoords+1);
$ligavgZ = $ligavgZ/($#ligxyzcoords+1);
print "Ligand center of coords (x y z): ".$ligavgX." ".$ligavgY." ".$ligavgZ."\n";

##################### SET LIGAND ATOMS AS SPHERES IF CRYSTALLOGRAPHIC ####################
$sphxyzcoords[0] = 0;    # Initialize sphxyzcoords array
# Check if flag .useligsph is "on" for crystallographic ligand, and if so keep spheres for final output
$uselig = "off";   # default value, (initialized for later algorithms to check)
$USELIGSPH = ".useligsph";
if (! -e $USELIGSPH) {
   print "useligsph flag does not exist.\n";
 }else{
#   $uselig = `cat .useligsph`;
   $uselig = "off";  # niu, for delphi input, ligand information is not necessary
   chomp $uselig;
   if ($uselig ne "on") {
      print "useligsph not on\n";
    }else{
      print "Using crystallographic ligand as spheres for output\n";
      pop @sphxyzcoords;
      for ($i = 0; $i<=$#ligxyzcoords; $i++) {
         $sphxyzcoords[$i][0] = "crystal";
         $sphxyzcoords[$i][1] = $ligxyzcoords[$i][0];
         $sphxyzcoords[$i][2] = $ligxyzcoords[$i][1];
         $sphxyzcoords[$i][3] = $ligxyzcoords[$i][2];
         $sphxyzcoords[$i][4] = 0;
         $sphxyzcoords[$i][5] = 0;
         $sphxyzcoords[$i][6] = 1;   # special FLAG for crystallographic ligand spheres
      }
   }
}

##################### READ SPHERES FROM SPHGEN ###########################################
#$SPH = "sph/sph";
$SPH = $ARGV[1]; #specify via command line now
open (SPH, "<$SPH") or die "My error, cannot open SPH: $!\n";
while (<SPH>) {                     # read in SPH lines
   @spharray = split;
   push @sphitems, [ @spharray ];
}
close (SPH);
# Put cluster 0 into array elements
$line = 0;
$numspherelines = 0;  # number of lines to expect containing sphere coords after "cluster"
@temparray = [0];
if ($#sphxyzcoords != 0) { $initialsphxyzsize=$#sphxyzcoords+1; } else { $initialsphxyzsize=0; }  # if have crystallographic spheres
if (@sphxyzcoords eq @temparray) { pop @sphxyzcoords; }   # make sure sphxyzcoords array is completely empty
while ( $line <= $#sphitems ) { 
   if ( $sphitems[$line][1] eq "0" ) {  # starting line of cluster 0 followed by coords
      $numspherelines = $sphitems[$line][7];    # Item 7 on "cluster" text line is num of spheres
      print "Starting number of spheres from ".$SPH." cluster 0: ".$numspherelines."\n";
      $line++;   # Move to first line with coordinates after header cluster 0 line
      while ( $line <= $#sphitems ){     # Cluster 0 goes to the end of the file
         push @sphxyzcoords, [ @{ $sphitems[$line] } ];
         $line++;
      }
   }
   $line++;
}
for ($i=$initialsphxyzsize; $i<=$#sphxyzcoords; $i++) { # set spheres read from sph/sph to non-crystallographic flag
   $sphxyzcoords[$i][6] = 0;
}
print "FINISHED READING IN FILE ".$SPH." CLUSTER 0 \n";
# for $array_ref ( @sphxyzcoords ) {
#    print "\t [ @$array_ref ], \n";
# }

##################### ELIMINATE SPHERES TOO FAR FROM CENTER OF LIGAND ####################
# Trim cluster 0 spheres to boundary of points within $M angstroms of center of ligand coords
$bad = 0;            # boolean for whether or not a coordinate line should be kept
$line = 0;
while ($line <= $#sphxyzcoords) {
   if ( ($sphxyzcoords[$line][1] < $ligavgX-$M ) || ( $sphxyzcoords[$line][1] > $ligavgX+$M) ) {$bad = 1;}
   if ( ($sphxyzcoords[$line][2] < $ligavgY-$M ) || ( $sphxyzcoords[$line][2] > $ligavgY+$M) ) {$bad = 1;}
   if ( ($sphxyzcoords[$line][3] < $ligavgZ-$M ) || ( $sphxyzcoords[$line][3] > $ligavgZ+$M) ) {$bad = 1;}
   if ( ($bad == 1) && ($sphxyzcoords[$line][6] != 1)) {  # not splicing out a crystallographic ligand
      splice(@sphxyzcoords, $line, 1);
      $line=$line-1;
   }
   $line++;
   $bad = 0;
}
print "Number of sphere points after cutting off spheres too far from center of ligand (".$M." angstroms): ".($#sphxyzcoords+1)."\n";

##################### READ RECEPTOR ATOMS ################################################
# Read in receptor rec.crg coordinates
#$REC = "grids/rec.crg";
$REC = $ARGV[2]; #specify via commandline
open (REC, "<$REC") or die "My error, cannot open REC.CRG: $!\n";
while (<REC>) {                     # read in REC lines
        @recarray = split;
        push @recitems, [ @recarray ];
}
close (REC);
for ($x=0; $x<=$#recitems; $x++) {
   if (($recitems[$x][0] eq "ATOM") || ($recitems[$x][0] eq "HETATM")) {
       if ($recitems[$x][5] =~ m/\./) {
           push @recxyzcoords, [ @{ $recitems[$x] } [ 5..7, 2..4 ] ];   # coordinates [0..2] then atomname,residue,num [3..5]
           push @{$rec{$recitems[$x][4]}}, $recxyzcoords[$#recxyzcoords];  # hash of rec indexed by residue number
       } else {
           push @recxyzcoords, [ @{ $recitems[$x] } [ 6..8, 2..3, 5 ] ];   # coordinates [0..2] then atomname,residue,num [3..5]
           push @{$rec{$recitems[$x][5]}}, $recxyzcoords[$#recxyzcoords];  # hash of rec indexed by residue number
       }
      if (defined ($recentries{$recxyzcoords[$#recxyzcoords][5]})) {
         $recentries{$recxyzcoords[$#recxyzcoords][5]}++;
       }else{
         $recentries{$recxyzcoords[$#recxyzcoords][5]} = 0;
      }
   }
}

##################### ELIMINATE SPHERES TOO FAR AND TOO CLOSE TO ANY RECEPTOR ATOM #######
# Remove spheres that have contacts with any receptor atom greater than $R cutoff
# Remove spheres that have contacts with any receptor atom closer than $Rclose cutoff
for ($i = 0; $i <= $#sphxyzcoords; $i++) { #These nested loops determines smallest distance of each sphere to any receptor atom
   $dist=sqrt((($recxyzcoords[0][0]-$sphxyzcoords[$i][1])**2)+(($recxyzcoords[0][1]-$sphxyzcoords[$i][2])**2)+(($recxyzcoords[0][2]-$sphxyzcoords[$i][3])**2));
   for ($j = 1; $j <= $#recxyzcoords; $j++) {
       $distB=sqrt((($recxyzcoords[$j][0]-$sphxyzcoords[$i][1])**2)+(($recxyzcoords[$j][1]-$sphxyzcoords[$i][2])**2)+(($recxyzcoords[$j][2]-$sphxyzcoords[$i][3])**2));
      if ($distB < $dist) {
         $dist = $distB;
      }
   }
   if (($dist > $R) || ($dist < $Rclose) ){ # may splice out a crystallographic sphere
      if ($sphxyzcoords[$i][6]==1) {print "WARNING: removing a crystallographic sphere that is distance $dist angstroms from closest receptor atom\n";}
      splice(@sphxyzcoords, $i, 1);
      $i=$i-1;
   }
}
print "Number of spheres after removing too far (>$R angstroms) and too close to receptor (<$Rclose angstroms): ".($#sphxyzcoords+1)."\n";

##################### FLAG SPHERES NEAR POLAR RECEPTOR ATOMS #############################
# Flag spheres in position [7] with 0 for near nonpolar and not polar, 1 if near polar, default -1 if neither
for ($i=0; $i<=$#sphxyzcoords; $i++) {   # initialize entries hash table
   $polarentries{$i} = -1;
   $nonpolarentries{$i} = -1;
   $sphxyzcoords[$i][7]=-1;   # initialize flag as neither polar nor nonpolar
   $sphxyzcoords[$i][9]=1000; # initialize distance from the mean angle
}
for ($i=0; $i<=$#sphxyzcoords; $i++) {
   for ($j=0; $j<=$#recxyzcoords; $j++) {
      $dist=sqrt((($recxyzcoords[$j][0]-$sphxyzcoords[$i][1])**2)+(($recxyzcoords[$j][1]-$sphxyzcoords[$i][2])**2)+(($recxyzcoords[$j][2]-$sphxyzcoords[$i][3])**2));
      if (($dist<=$polardist)&&((substr($recxyzcoords[$j][3],0,1)eq"N")||(substr($recxyzcoords[$j][3],0,1)eq"O")||(substr($recxyzcoords[$j][3],0,1)eq"S"))){
         push @{$polar{$i}}, $recxyzcoords[$j];
         $polarentries{$i}++;
      }
      if (($dist<=$Hpolardist)&&(substr($recxyzcoords[$j][3],0,1)eq"H")) {
         push @{$polar{$i}}, $recxyzcoords[$j];
         $polarentries{$i}++;
      }
      if (($dist<=$nonpolardist)&&(substr($recxyzcoords[$j][3],0,1)eq"C")){
         push @{$nonpolar{$i}}, $recxyzcoords[$j];
         $nonpolarentries{$i}++;
      }
   }
}
$numpolar=0;
$numnonpolar=0;
foreach ( keys %nonpolar) {
   $nonp = $_;
   $found = 0;
   foreach ( keys %polar) {
      if ($nonp eq $_) { $found = 1; }
   }
   if ($found == 0) { $sphxyzcoords[$nonp][7] = -2; $numnonpolar++;}  # sphere $nonp is not in polar hash table, and is in nonpolar table
}
# Define math required for directionality calculation
my $pi = 3.14159265359;
sub rad_to_deg { ($_[0]/$pi) * 180 }
sub acos { atan2( sqrt(1 - $_[0] * $_[0]), $_[0] ) }

# Subroutine costheta returns cos(theta) based on residue number and atom name for polar and base
sub theta($) {
   ($basename) = @_;
   @polatom = ($polar{$_}[$i][0], $polar{$_}[$i][1], $polar{$_}[$i][2]);
   for ($j=0; $j<=$recentries{($polar{$sphnum}[$i][5])}; $j++) {
      if ($rec{($polar{$sphnum}[$i][5])}[$j][3] eq "$basename") {
         @baseatom = ($rec{($polar{$sphnum}[$i][5])}[$j][0], $rec{($polar{$sphnum}[$i][5])}[$j][1], $rec{($polar{$sphnum}[$i][5])}[$j][2]);
      }
   }
   @base2pol = ($baseatom[0]-$polatom[0], $baseatom[1]-$polatom[1], $baseatom[2]-$polatom[2]);
   @pol2sph = ($sphxyzcoords[$sphnum][1]-$polatom[0], $sphxyzcoords[$sphnum][2]-$polatom[1], $sphxyzcoords[$sphnum][3]-$polatom[2]);
   $lengthV1=sqrt((($base2pol[0])**2)+(($base2pol[1])**2)+(($base2pol[2])**2));
   $lengthV2=sqrt((($pol2sph[0])**2)+(($pol2sph[1])**2)+(($pol2sph[2])**2));
   @unitV1 = ($base2pol[0]/$lengthV1, $base2pol[1]/$lengthV1, $base2pol[2]/$lengthV1);
   @unitV2 = ($pol2sph[0]/$lengthV2, $pol2sph[1]/$lengthV2, $pol2sph[2]/$lengthV2);
   $cosangle = (($unitV1[0]*$unitV2[0])+($unitV1[1]*$unitV2[1])+($unitV1[2]*$unitV2[2]));
   $angle=rad_to_deg(acos($cosangle));
   return $angle;
}

# Assign cos(theta) value to sphere based on position of polar atom and base receptor atom
foreach (keys %polar) {
   $i=0;
   while ($i<=$polarentries{$_}) {
      $sphnum=$_;
      #if ($polar{$_}[$i][3] eq "HN") {                                                                                                                                         
      if ($polar{$_}[$i][3] eq "H") { # REDUCE                                                                                                                                  
         $newangle=theta("N");
         if((abs($newangle-155.21)<$sphxyzcoords[$_][9])&&($newangle>=129.2792)&&($newangle<=181.1408)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-155.21);} # Mean 155.21 Std. Dev. 13.23
      }
      elsif ($polar{$_}[$i][3] eq "O") {
         $newangle=theta("C");
         if((abs($newangle-140.70)<$sphxyzcoords[$_][9])&&($newangle>=104.4596)&&($newangle<=176.9404)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-140.70);} # Mean 140.70 Std. Dev. 18.49
      }
      #elsif (($polar{$_}[$i][3] eq "HNE")&&($polar{$_}[$i][4]eq"ARG")) {                                                                                                       
      elsif (($polar{$_}[$i][3] eq "HE")&&($polar{$_}[$i][4]eq"ARG")) { # REDUCE                                                                                                
         $newangle=theta("NE");
         if((abs($newangle-150.13)<$sphxyzcoords[$_][9])&&($newangle>=121.4356)&&($newangle<=178.8244)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-150.13);} # Mean 150.13 Std. Dev. 14.64
      }
      #elsif ((($polar{$_}[$i][3] eq "HN11")||($polar{$_}[$i][3] eq "HN12"))&&($polar{$_}[$i][4]eq"ARG")) {
      elsif ((($polar{$_}[$i][3] eq "HH11")||($polar{$_}[$i][3] eq "HH12"))&&($polar{$_}[$i][4]eq"ARG")) { # REDUCE
         $newangle=theta("NH1");
         if((abs($newangle-147.29)<$sphxyzcoords[$_][9])&&($newangle>=114.8520)&&($newangle<=179.7280)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-147.29);} # Mean 147.29 Std. Dev. 16.55
      }
      #elsif ((($polar{$_}[$i][3] eq "HN21")||($polar{$_}[$i][3] eq "HN22"))&&($polar{$_}[$i][4]eq"ARG")) {
      elsif ((($polar{$_}[$i][3] eq "HH21")||($polar{$_}[$i][3] eq "HH22"))&&($polar{$_}[$i][4]eq"ARG")) { # REDUCE
         $newangle=theta("NH2");
         if((abs($newangle-146.21)<$sphxyzcoords[$_][9])&&($newangle>=114.4384)&&($newangle<=177.9816)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-146.21);} # Mean 146.21 Std. Dev. 16.21
      }
      #elsif ((($polar{$_}[$i][3] eq "HND1")||($polar{$_}[$i][3] eq "HND2"))&&($polar{$_}[$i][4]eq"ASN")) {
      elsif ((($polar{$_}[$i][3] eq "HD21")||($polar{$_}[$i][3] eq "HD22"))&&($polar{$_}[$i][4]eq"ASN")) { # REDUCE
         $newangle=theta("ND2");
         if((abs($newangle-151.25)<$sphxyzcoords[$_][9])&&($newangle>=114.5784)&&($newangle<=187.9216)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-151.25);} # Mean 151.25 Std. Dev. 18.71
      }
      elsif (($polar{$_}[$i][3] eq "OD1")&&($polar{$_}[$i][4]eq"ASN")) {
         $newangle=theta("CG");
         if((abs($newangle-131.66)<$sphxyzcoords[$_][9])&&($newangle>=97.3796)&&($newangle<=165.9404)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-131.66);} # Mean 131.66 Std. Dev. 17.49
      }
      elsif (($polar{$_}[$i][3] eq "OD1")&&($polar{$_}[$i][4]eq"ASP")) {
         $newangle=theta("CG");
         if((abs($newangle-124.66)<$sphxyzcoords[$_][9])&&($newangle>=90.0464)&&($newangle<=159.2736)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-124.66);} # Mean 124.66 Std. Dev. 17.66
      }
      elsif (($polar{$_}[$i][3] eq "OD2")&&($polar{$_}[$i][4]eq"ASP")) {
         $newangle=theta("CG");
         if((abs($newangle-121.82)<$sphxyzcoords[$_][9])&&($newangle>=92.9884)&&($newangle<=150.6516)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-121.82);} # Mean 121.82 Std. Dev. 14.71
      }
      elsif (($polar{$_}[$i][3] eq "HG")&&($polar{$_}[$i][4]eq"CYS")) {
         $newangle=theta("SG");
         if((abs($newangle-154.82)<$sphxyzcoords[$_][9])&&($newangle>=128.9284)&&($newangle<=180.7116)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-154.82);} # Mean 154.82 Std. Dev. 13.21
      }
      elsif (($polar{$_}[$i][3] eq "SG")&&($polar{$_}[$i][4]eq"CYS")) {
         $newangle=theta("CB");
         if((abs($newangle-112.75)<$sphxyzcoords[$_][9])&&($newangle>=87.1720)&&($newangle<=138.3280)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-112.75);} # Mean 112.75 Std. Dev. 13.05
      }
      elsif (($polar{$_}[$i][3] eq "SG")&&($polar{$_}[$i][4]eq"CYX")) {
         $newangle=theta("CB");
         if((abs($newangle-147.05)<$sphxyzcoords[$_][9])&&($newangle>=127.5480)&&($newangle<=166.5520)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-147.05);} # Mean 147.05 Std. Dev. 9.95
      }
      elsif (($polar{$_}[$i][3] eq "OE1")&&($polar{$_}[$i][4]eq"GLU")) {
         $newangle=theta("CD");
         if((abs($newangle-123.35)<$sphxyzcoords[$_][9])&&($newangle>=86.4824)&&($newangle<=160.2176)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-123.35);} # Mean 123.35 Std. Dev. 18.81
      }
      elsif (($polar{$_}[$i][3] eq "OE2")&&($polar{$_}[$i][4]eq"GLU")) {
         $newangle=theta("CD");
         if((abs($newangle-123.83)<$sphxyzcoords[$_][9])&&($newangle>=90.4512)&&($newangle<=157.2088)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-123.83);} # Mean 123.83 Std. Dev. 17.03
      }
      #elsif ((($polar{$_}[$i][3] eq "HNE1")||($polar{$_}[$i][3] eq "HNE2"))&&($polar{$_}[$i][4]eq"GLN")) {
      elsif ((($polar{$_}[$i][3] eq "HE21")||($polar{$_}[$i][3] eq "HE22"))&&($polar{$_}[$i][4]eq"GLN")) { # REDUCE
         $newangle=theta("NE2");
         if((abs($newangle-147.62)<$sphxyzcoords[$_][9])&&($newangle>=108.714)&&($newangle<=186.5260)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-147.62);} # Mean 147.62 Std. Dev. 19.85
      }
      elsif (($polar{$_}[$i][3] eq "OE1")&&($polar{$_}[$i][4]eq"GLN")) {
         $newangle=theta("CD");
         if((abs($newangle-129.76)<$sphxyzcoords[$_][9])&&($newangle>=96.1068)&&($newangle<=163.4132)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-129.76);} # Mean 129.76 Std. Dev. 17.17
      }
      #elsif (($polar{$_}[$i][3] eq "HND")&&(($polar{$_}[$i][4]eq"HIS")||($polar{$_}[$i][4]eq"HID")||($polar{$_}[$i][4]eq"HIP"))) {
      # Note that here HIS is a tarted (more polar) residue
      elsif (($polar{$_}[$i][3] eq "HD1")&&(($polar{$_}[$i][4]eq"HIS")||($polar{$_}[$i][4]eq"HID")||($polar{$_}[$i][4]eq"HIP"))) {
         $newangle=theta("ND1");
         if((abs($newangle-149.52)<$sphxyzcoords[$_][9])&&($newangle>=116.6508)&&($newangle<=182.3892)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-149.52);} # Mean 149.52 Std. Dev. 16.77
      }
      #elsif (($polar{$_}[$i][3] eq "HNE")&&(($polar{$_}[$i][4]eq"HIE")||($polar{$_}[$i][4]eq"HIP"))) {#
      elsif (($polar{$_}[$i][3] eq "HE2")&&(($polar{$_}[$i][4]eq"HIE")||($polar{$_}[$i][4]eq"HIP"))) {# REDUCE name
         $newangle=theta("NE2");
         if((abs($newangle-152.56)<$sphxyzcoords[$_][9])&&($newangle>=117.6132)&&($newangle<=187.5068)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-152.56);} # Mean 152.56 Std. Dev. 17.83
      }
      elsif (($polar{$_}[$i][3] eq "ND1")&&($polar{$_}[$i][4]eq"HIS")) {
         $newangle=theta("CG");  # ring has symmetry so CG is acceptable
         if((abs($newangle-130.34)<$sphxyzcoords[$_][9])&&($newangle>=110.8968)&&($newangle<=149.7832)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-130.34);} # Mean 130.34 Std. Dev. 9.92
      }
      elsif (($polar{$_}[$i][3] eq "NE2")&&($polar{$_}[$i][4]eq"HIS")) {
         $newangle=theta("CD2");  # ring has symmetry so CD2 is acceptable
         if((abs($newangle-131.83)<$sphxyzcoords[$_][9])&&($newangle>=111.6224)&&($newangle<=152.0376)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-131.83);} # Mean 131.83 Std. Dev. 10.31
      }
      #elsif((($polar{$_}[$i][3]eq"HNZ1")||($polar{$_}[$i][3]eq"HNZ2")||($polar{$_}[$i][3]eq"HNZ3"))&&($polar{$_}[$i][4]eq"LYS")){
      elsif((($polar{$_}[$i][3]eq"HZ1")||($polar{$_}[$i][3]eq"HZ2")||($polar{$_}[$i][3]eq"HZ3"))&&($polar{$_}[$i][4]eq"LYS")){ # REDUCE
         $newangle=theta("NZ");
         if((abs($newangle-139.81)<$sphxyzcoords[$_][9])&&($newangle>=102.1388)&&($newangle<=177.4812)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-139.81);} # Mean 139.81 Std. Dev. 19.22
      }
      elsif (($polar{$_}[$i][3] eq "SD")&&($polar{$_}[$i][4]eq"MET")) {
         $newangle=theta("CG");  # CG is probably pointing for a better angle than CE
         if((abs($newangle-139.81)<$sphxyzcoords[$_][9])&&($newangle>=92.3324)&&($newangle<=169.0076)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-139.81);} # Mean 139.81 Std. Dev. 19.22
      }
      #elsif (($polar{$_}[$i][3] eq "HOG")&&($polar{$_}[$i][4]eq"SER")) {
      elsif (($polar{$_}[$i][3] eq "HG")&&($polar{$_}[$i][4]eq"SER")) { # REDUCE
         $newangle=theta("OG");
         if((abs($newangle-159.76)<$sphxyzcoords[$_][9])&&($newangle>=130.2032)&&($newangle<=189.3168)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-159.76);} # Mean 159.76 Std. Dev. 15.08
      }
      elsif (($polar{$_}[$i][3] eq "OG")&&($polar{$_}[$i][4]eq"SER")) {
         $newangle=theta("CB");
         if((abs($newangle-124.55)<$sphxyzcoords[$_][9])&&($newangle>=93.2292)&&($newangle<=155.8708)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-124.55);} # Mean 124.55 Std. Dev. 15.98
      }
      #elsif (($polar{$_}[$i][3] eq "HOG")&&($polar{$_}[$i][4]eq"THR")) {
      elsif (($polar{$_}[$i][3] eq "HG")&&($polar{$_}[$i][4]eq"THR")) { # REDUCE
         $newangle=theta("OG1");
         if((abs($newangle-163.00)<$sphxyzcoords[$_][9])&&($newangle>=135.8148)&&($newangle<=190.1852)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-163.00);} # Mean 163.00 Std. Dev. 13.87
      }
      elsif (($polar{$_}[$i][3] eq "OG1")&&($polar{$_}[$i][4]eq"THR")) {
         $newangle=theta("CB");
         if((abs($newangle-125.41)<$sphxyzcoords[$_][9])&&($newangle>=97.0488)&&($newangle<=153.7712)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-125.41);} # Mean 125.41 Std. Dev. 14.47
      }
      #elsif (($polar{$_}[$i][3] eq "HNE")&&($polar{$_}[$i][4]eq"TRP")) {
      elsif (($polar{$_}[$i][3] eq "HE1")&&($polar{$_}[$i][4]eq"TRP")) { # REDUCE
         $newangle=theta("NE1");
         if((abs($newangle-153.50)<$sphxyzcoords[$_][9])&&($newangle>=124.4332)&&($newangle<=182.5668)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-153.50);} # Mean 153.50 Std. Dev. 14.83
      }
      #elsif (($polar{$_}[$i][3] eq "HOH")&&($polar{$_}[$i][4]eq"TYR")) {
      elsif (($polar{$_}[$i][3] eq "HH")&&($polar{$_}[$i][4]eq"TYR")) { # REDUCE
         $newangle=theta("OH");
         if((abs($newangle-147.32)<$sphxyzcoords[$_][9])&&($newangle>=112.3928)&&($newangle<=182.2472)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-147.32);} # Mean 147.32 Std. Dev. 17.82
      }
      elsif (($polar{$_}[$i][3] eq "OH")&&($polar{$_}[$i][4]eq"TYR")) {
         $newangle=theta("CZ");
         if((abs($newangle-119.55)<$sphxyzcoords[$_][9])&&($newangle>=94.8540)&&($newangle<=144.2460)) {$sphxyzcoords[$_][7]=$newangle; $sphxyzcoords[$_][9]=abs($newangle-119.55);} # Mean 119.55 Std. Dev. 12.60
      }
      $i++;
   }
   if (($sphxyzcoords[$_][7]!=-1)&&($sphxyzcoords[$_][7]!=-2)) { 
      $numpolar++;
   }
}
print "Number of spheres polar: $numpolar, nonpolar: $numnonpolar, out of current total: ".($#sphxyzcoords+1)."\n";

##################### ELIMINATE SPHERES TOO CLOSE TO EACH OTHER ##########################
# Trim spheres so that none are closer than $gridsize angstroms from each other
# Pseudocode: (1) Establish grid system by finding extents of x y z sphere coordinates and using $gridsize 3D boxes
# Pseudocode: (1) Establish grid system by finding extents of x y z sphere coordinates and using $gridsize 3D boxes
# (2) nested loop through all grid points in x y z, checking for all spheres within current $gridsize box
# (3) keep the sphere with the lowest original line number and splice out others (except crystallographic spheres)
$maxX1 = $sphxyzcoords[0][1];
$maxX2 = $sphxyzcoords[0][1];
$maxY1 = $sphxyzcoords[0][2];
$maxY2 = $sphxyzcoords[0][2];
$maxZ1 = $sphxyzcoords[0][3];
$maxZ2 = $sphxyzcoords[0][3];
for ($i = 1; $i <= $#sphxyzcoords; $i++) {
   if ( $sphxyzcoords[$i][1] < $maxX1 ) { $maxX1 = $sphxyzcoords[$i][1]; } 
   if ( $sphxyzcoords[$i][1] > $maxX2 ) { $maxX2 = $sphxyzcoords[$i][1]; }
   if ( $sphxyzcoords[$i][2] < $maxY1 ) { $maxY1 = $sphxyzcoords[$i][2]; }
   if ( $sphxyzcoords[$i][2] > $maxY2 ) { $maxY2 = $sphxyzcoords[$i][2]; }
   if ( $sphxyzcoords[$i][3] < $maxZ1 ) { $maxZ1 = $sphxyzcoords[$i][3]; }
   if ( $sphxyzcoords[$i][3] > $maxZ2 ) { $maxZ2 = $sphxyzcoords[$i][3]; }
}
print "Extents of grid to search are X1 X2 Y1 Y2 Z1 Z2: ".$maxX1." ".$maxX2." ".$maxY1." ".$maxY2." ".$maxZ1." ".$maxZ2."\n";
for ($x = $maxX1; $x <= $maxX2; $x=$x+$gridsize) {        # grid x coordinates
   for ($y = $maxY1; $y <= $maxY2; $y=$y+$gridsize) {     # grid y coordinates
      for ($z = $maxZ1; $z <= $maxZ2; $z=$z+$gridsize) {  # grid z coordinates
         @gridbox = [0];
         $index=0;
         for ($i=0; $i<=$#sphxyzcoords; $i++) {
            if(($sphxyzcoords[$i][1]>=$x-($gridsize/2))&&($sphxyzcoords[$i][1]<$x+($gridsize/2))&&($sphxyzcoords[$i][2]>=$y-($gridsize/2))&&($sphxyzcoords[$i][2]<$y+($gridsize/2))&&($sphxyzcoords[$i][3]>=$z-($gridsize/2))&&($sphxyzcoords[$i][3]<$z+($gridsize/2))){  # Found sphere in current grid box
               $gridbox[$index][0] = $i;     # Keep track of line number from original sphxyzcoords for splicing out below
               $gridbox[$index][1] = sqrt((($sphxyzcoords[$i][1]-$x)**2)+(($sphxyzcoords[$i][2]-$y)**2)+(($sphxyzcoords[$i][3]-$z)**2)); #dist2center
               $gridbox[$index][2] = $sphxyzcoords[$i][6];   # Associate status of sphere crystallographic flag
               $gridbox[$index][3] = $sphxyzcoords[$i][7];   # Associate status of sphere color flag
               $gridbox[$index][4] = $sphxyzcoords[$i][9];   # Associate value from the mean angle
               $index++;
            } 
         }  # Finished finding all spheres within current grid box; put them all in @gridbox array
         if ($#gridbox >= 1) {  # Found more than 1 sphere in current grid box
            @sortedgridbox = sort{$a->[1] <=> $b->[1]} @gridbox; # sorted spheres in gridbox by ascending distance to center of grid box
            $crystal=0;
            $polar=0;
            for ($i=0; $i<=$#sortedgridbox; $i++) { # count number of crystallographic and polar spheres
               if ($sortedgridbox[$i][2]==1) {    # Found crystallographic sphere
                  $crystal++;    # counter for crystallographic spheres in current gridbox
               }
               if (($sortedgridbox[$i][3]!=-1)&&($sortedgridbox[$i][3]!=-2)) { # Found polar sphere
                  $polar++;      # counter for polar spheres in current gridbox
               }
            }
            if ($crystal>0) {  # found crystal in box, so splice out all non-crystallographic spheres
               for ($line=0; $line<=$#sortedgridbox; $line++) {
                  if ($sortedgridbox[$line][2]!=1) {
                     splice (@sphxyzcoords,$sortedgridbox[$line][0],1);
                     for ($j=($line+1); $j<=$#sortedgridbox; $j++) {  # renumber lines since spliced one out
                        if ($sortedgridbox[$j][0] > $sortedgridbox[$line][0]) {$sortedgridbox[$j][0]--;} #decrease line num if spliced is before 
                     }
                  }
               }
            }elsif ($polar>0) { # no crystal but found polar, so splice out all non-polar spheres and keep best polar sphere
               @polarsortedgridbox = sort{$a->[4] <=> $b->[4]} @sortedgridbox; # sort spheres in gridbox by ascending difference from the mean angle theta
               for ($line=1; $line<=$#polarsortedgridbox; $line++) { # splice out all but the sphere closest to mean theta
                  splice (@sphxyzcoords,$polarsortedgridbox[$line][0],1);
                  for ($j=($line+1); $j<=$#polarsortedgridbox; $j++) {  # renumber lines since spliced one out
                     if ($polarsortedgridbox[$j][0]>$polarsortedgridbox[$line][0]){$polarsortedgridbox[$j][0]--;} #decrease line num if splice is before
                  }
               }
            }else{  # splice out all spheres keeping only 1 closest to center
               for ($line=1; $line<=$#sortedgridbox; $line++) {
                  splice (@sphxyzcoords,$sortedgridbox[$line][0],1);
                  for ($j=($line+1); $j<=$#sortedgridbox; $j++) {  # renumber lines since spliced one out
                     if ($sortedgridbox[$j][0]>$sortedgridbox[$line][0]) {$sortedgridbox[$j][0]--;} #decrease line num if splice is before
                  }
               }
            }
         }
      }
   }
}
print "Number of spheres after removing spheres too close to each other (approximately <".$gridsize." angstroms): ".($#sphxyzcoords+1)."\n";

##################### ELIMINATE SPHERES INDIVIDUALLY TOO CLOSE TO EACH OTHER #############
# Second pass through spheres that remain to do sweeps around each sphere in 3D
# If find that 2 spheres are within <0.5 A from each other then delete one of them
for ($i=0; $i<=$#sphxyzcoords; $i++) {
   for ($j=0; $j<=$#sphxyzcoords; $j++) {
      if(($j!=$i)&&(sqrt((($sphxyzcoords[$i][1]-$sphxyzcoords[$j][1])**2)+(($sphxyzcoords[$i][2]-$sphxyzcoords[$j][2])**2)+(($sphxyzcoords[$i][3]-$sphxyzcoords[$j][3])**2))<$tooclose)) {
         if (($sphxyzcoords[$j][6]!=1)&&(($sphxyzcoords[$j][7]==-1)||($sphxyzcoords[$j][7]==-2))) {  # eliminate j if not crystallographic and not polar
            splice (@sphxyzcoords,$j,1);
            $j--;
            if ($i>$j) { $i--; }    # adjusts for splicing changes index of records after splice point
          }elsif (($sphxyzcoords[$i][6]!=1)&&(($sphxyzcoords[$i][7]==-1)||($sphxyzcoords[$j][7]==-2))) {  # eliminate i if not crystallographic and not polar
            splice (@sphxyzcoords,$i,1);
            $i--;
            if ($j>$i) { $j--; }    # adjusts for splicing changes index of records after splice point
          }elsif (($sphxyzcoords[$j][6]!=1)&&($sphxyzcoords[$j][9]>=$sphxyzcoords[$i][9])) {  # eliminate j if not crystallographic and worse polar farther from mean angle theta
            splice (@sphxyzcoords,$j,1);
            $j--;
            if ($i>$j) { $i--; }    # adjusts for splicing changes index of records after splice point
          }elsif (($sphxyzcoords[$i][6]!=1)&&($sphxyzcoords[$i][9]>$sphxyzcoords[$j][9])) {  # eliminate i if not crystallographic and worse polar farther from mean angle theta
            splice (@sphxyzcoords,$i,1);
            $i--;
            if ($j>$i) { $j--; }    # adjusts for splicing changes index of records after splice point
          }elsif ($sphxyzcoords[$j][6]!=1) {  # eliminate j if not crystallographic
            splice (@sphxyzcoords,$j,1);
            $j--;
            if ($i>$j) { $i--; }    # adjusts for splicing changes index of records after splice point
          }elsif ($sphxyzcoords[$i][6]!=1) {  # eliminate i if not crystallographic
            splice (@sphxyzcoords,$i,1);
            $i--;
            if ($j>$i) { $j--; }    # adjusts for splicing changes index of records after splice point
          }else{                    # both i and j are crystallographic, but <0.8 from each other
            print "WARNING: spheres $i and $j are both crystallographic and are <$tooclose to each other (not removing either of them)\n";
         }
      }
   }
}
print "Number of spheres after second pass removing spheres absolutely <$tooclose from each other: ".($#sphxyzcoords+1)."\n";

##################### ELIMINATE SPHERES THAT ARE NOT CONTINUOUS WITH CENTER OF LIGAND ####
if ($uselig ne "on") { # If not crystallographic ligand, then determine point closest to center of spheres 
   print "No crystallographic spheres present, so determining sphere nearest center of ligand (x y z): ".$ligavgX." ".$ligavgY." ".$ligavgZ."\n";
   $closesttocenter = 0;
   $disttocenter=sqrt((($sphxyzcoords[0][1]-$ligavgX)**2)+(($sphxyzcoords[0][2]-$ligavgY)**2)+(($sphxyzcoords[0][3]-$ligavgZ)**2));
   for ($i = 1; $i <= $#sphxyzcoords; $i++) {
      $newdist=sqrt((($sphxyzcoords[$i][1]-$ligavgX)**2)+(($sphxyzcoords[$i][2]-$ligavgY)**2)+(($sphxyzcoords[$i][3]-$ligavgZ)**2));
      if ($newdist < $disttocenter) {
         $closesttocenter = $i;
         $disttocenter = $newdist;
      }
   }
   if ($sphxyzcoords[$closesttocenter][6]!=1) {
      $sphxyzcoords[$closesttocenter][6] = 2; # center point kept in output, along with crystallographic spheres already flagged
   } else { $sphxyzcoords[$closesttocenter][6] = 3; } # point is both crystallographic and center of spheres
   print "Center point: ".$closesttocenter." xyz: ".$sphxyzcoords[$closesttocenter][1]." ".$sphxyzcoords[$closesttocenter][2]." ".$sphxyzcoords[$closesttocenter][3]."\n";
}else{ print "Crystallographic spheres are present, so doing continuity based on crystallographic spheres\n"; }
# Do successive cycles expanding points to keep based on continuity from continuous points (center point) and from any crystallographic sphere
do {  # outermost do loop if not enough spheres are continuous, then increase continuity parameter
   $numspheres=0;
   for ($i = 0; $i <= $#sphxyzcoords; $i++) {
      if ($sphxyzcoords[$i][6]!=0) {   # counting how many spheres are flagged as continuous
         $numspheres++;
      }
   }
   $repeats=1;
   while (($numspheres<=$maxspheres) && ($repeats <= $numrepeats)) {  # $numrepeats = 10 is more than enough to find all continuous points
      for ($k = 0; $k <= $#sphxyzcoords; $k++) {
         $i = 0;
         while (($sphxyzcoords[$k][6]==0) && ($i <= $#sphxyzcoords)) {
            if(($i!=$k)&&(sqrt((($sphxyzcoords[$k][1]-$sphxyzcoords[$i][1])**2)+(($sphxyzcoords[$k][2]-$sphxyzcoords[$i][2])**2)+(($sphxyzcoords[$k][3]-$sphxyzcoords[$i][3])**2))<=$continuity)&&($sphxyzcoords[$i][6]>=1)){  # flag [6] >= 1 means it is near crystallographic or continuous (center) point or both
               $sphxyzcoords[$k][6]=2;  # sphere $k flagged continuous since it is near a continuous sphere (flag 2) or crystallographic sphere (flag 1)
            }
            $i++;
         }
      }
      $repeats++;
   }
   $numspheres=0;
   for ($i = 0; $i <= $#sphxyzcoords; $i++) {
      if ($sphxyzcoords[$i][6]!=0) {   # counting how many spheres are flagged as continuous
         $numspheres++;
      }
   }
   if ( $numspheres < $minspheres ) { 
      $continuity = $continuity + $contincrement;
      print "Only ".$numspheres." spheres are continuous, so increasing continuity by ".$contincrement." to be ".$continuity."\n";
   }
} until (($numspheres >= $minspheres) || ($continuity > $maxcontinuity));
# splice out those not flagged as being continuous
$i = 0;
$currentlength=$#sphxyzcoords;
while ($i <= $currentlength) {
   if ($sphxyzcoords[$i][6]==0) {
      splice(@sphxyzcoords,$i,1);
      $i--;
      $currentlength--;
   }
   $i++;
}
print "After continuity checking, number of spheres is: ".($#sphxyzcoords+1)."\n";

##################### ELIMINATE SPHERES IF TOO MANY REMAIN  ##############################
# If more than permitted number of spheres for DOCK, eliminate spheres based on those farthest away
if ($#sphxyzcoords+1 > $maxspheres) {       # Arbitrary checkpoint for delphi to have no more than 120 spheres
   # Determine center of mass of match.sph
   $sphavgX = $sphxyzcoords[0][1];
   $sphavgY = $sphxyzcoords[0][2];
   $sphavgZ = $sphxyzcoords[0][3];
   for ($line = 1; $line <= $#sphxyzcoords; $line++) {
      $sphavgX = $sphavgX + $sphxyzcoords[$line][1];
      $sphavgY = $sphavgY + $sphxyzcoords[$line][2];
      $sphavgZ = $sphavgZ + $sphxyzcoords[$line][3];
   }
   $sphavgX = $sphavgX/($#sphxyzcoords+1);
   $sphavgY = $sphavgY/($#sphxyzcoords+1);
   $sphavgZ = $sphavgZ/($#sphxyzcoords+1);
   print "Spheres center of coords (x y z): ".$sphavgX." ".$sphavgY." ".$sphavgZ."\n";
   # Calculate "weighted distance" to eliminate spheres
   if ($uselig eq "on") {  # Eliminate spheres depending on having crystallographic ligand present
      print "Crystallographic ligand, so keeping spheres closest to ANY ligand atom, center of ligand (1/8), and center of spheres (1/8)\n";
      for ($i = 0; $i <= $#sphxyzcoords; $i++) { #These nested loops determines distance of each sphere to closest ligand atom
          $dist=sqrt((($ligxyzcoords[0][0]-$sphxyzcoords[$i][1])**2)+(($ligxyzcoords[0][1]-$sphxyzcoords[$i][2])**2)+(($ligxyzcoords[0][2]-$sphxyzcoords[$i][3])**2));
          for ($j = 1; $j <= $#ligxyzcoords; $j++) {
             $distB=sqrt((($ligxyzcoords[$j][0]-$sphxyzcoords[$i][1])**2)+(($ligxyzcoords[$j][1]-$sphxyzcoords[$i][2])**2)+(($ligxyzcoords[$j][2]-$sphxyzcoords[$i][3])**2));
             if ($distB < $dist) {
                $dist = $distB;
             }
          } 
          $sphxyzcoords[$i][8]=$dist;   # Weighted to 1 importance
      }
      # Uncomment lines below to add weight to ANY receptor atoms
#      for ($i = 0; $i <= $#sphxyzcoords; $i++) { #These nested loops determines smallest distance of each sphere to any receptor atom 
#          $dist=sqrt((($recxyzcoords[0][0]-$sphxyzcoords[$i][1])**2)+(($recxyzcoords[0][1]-$sphxyzcoords[$i][2])**2)+(($recxyzcoords[0][2]-$sphxyzcoords[$i][3])**2));
#          for ($j = 1; $j <= $#recxyzcoords; $j++) {
#             $distB=sqrt((($recxyzcoords[$j][0]-$sphxyzcoords[$i][1])**2)+(($recxyzcoords[$j][1]-$sphxyzcoords[$i][2])**2)+(($recxyzcoords[$j][2]-$sphxyzcoords[$i][3])**2));
#             if ($distB < $dist) {
#                $dist = $distB;
#             }
#          }
#          $sphxyzcoords[$i][8]=$sphxyzcoords[$i][8]+$dist;   # Weighted at 1-whole importance
#      }
      for ($i = 0; $i <= $#sphxyzcoords; $i++) { # determine distance of each sphere to center of ligand atoms
          $dist=sqrt((($ligavgX-$sphxyzcoords[$i][1])**2)+(($ligavgY-$sphxyzcoords[$i][2])**2)+(($ligavgZ-$sphxyzcoords[$i][3])**2));
          $sphxyzcoords[$i][8]=$sphxyzcoords[$i][8]+($dist/8); # Weighted at 1/8 of importance
      }
      for ($i = 0; $i <= $#sphxyzcoords; $i++) { # determine distance of each sphere to center of spheres
          $dist=sqrt((($sphavgX-$sphxyzcoords[$i][1])**2)+(($sphavgY-$sphxyzcoords[$i][2])**2)+(($sphavgZ-$sphxyzcoords[$i][3])**2));
          $sphxyzcoords[$i][8]=$sphxyzcoords[$i][8]+($dist/8); # Weighted at 1/8 of importance
      }
       ########################################################################################################
       # ALTERNATE DISTRIBUTIONS OF SPHERES ARE GENERATED BY MAKING THE FOLLOWING MODIFICATIONS:
       # setting item 7 can determine how points are weighted when cut off being >120 spheres
       # it can be set to the minimum distance to ANY receptor atom (use recxyzcoords in $distB= above)
       # it can be set to the minimum distance to ANY ligand atom (use ligxyzcoords in $distB= above) [current method]
       # it can be set to the distance to center of ligand (use ligavgX, ligavgY, ligavgZ in $distB= above)
       # it can be set to the distance to center of receptor (calculate and use recavgX, recavgY, recavgZ in $distB= above)
       # or it can be a summation of any of these (+), multiplication (*), division (/), squared (**2), or weighted average (constant * each term), etc.
       ########################################################################################################
   }else{   # No crystallographic spheres are present
      print "Not crystallographic ligand, so keeping spheres closest to center of receptor, center of ligand, and center of spheres\n";
      # Uncomment lines below to add weight to ANY receptor atoms
#      for ($i = 0; $i <= $#sphxyzcoords; $i++) { #These nested loops determines distance of each sphere to ANY receptor atom
#          $dist=sqrt((($recxyzcoords[0][0]-$sphxyzcoords[$i][1])**2)+(($recxyzcoords[0][1]-$sphxyzcoords[$i][2])**2)+(($recxyzcoords[0][2]-$sphxyzcoords[$i][3])**2));
#          for ($j = 1; $j <= $#recxyzcoords; $j++) {
#             $distB=sqrt((($recxyzcoords[$j][0]-$sphxyzcoords[$i][1])**2)+(($recxyzcoords[$j][1]-$sphxyzcoords[$i][2])**2)+(($recxyzcoords[$j][2]-$sphxyzcoords[$i][3])**2));
#             if ($distB < $dist) {
#                $dist = $distB;
#             }
#          }
#          $sphxyzcoords[$i][8]=$dist;   # Weighted at 1-whole importance
#      }
      $recavgX = $recxyzcoords[0][0];
      $recavgY = $recxyzcoords[0][1];
      $recavgZ = $recxyzcoords[0][2];
      for ($line = 1; $line <= $#recxyzcoords; $line++) {
         $recavgX = $recavgX + $recxyzcoords[$line][0];
         $recavgY = $recavgY + $recxyzcoords[$line][1];
         $recavgZ = $recavgZ + $recxyzcoords[$line][2];
      }
      $recavgX = $recavgX/($#recxyzcoords+1);
      $recavgY = $recavgY/($#recxyzcoords+1);
      $recavgZ = $recavgZ/($#recxyzcoords+1);
      print "Receptor center of coords (x y z): ".$recavgX." ".$recavgY." ".$recavgZ."\n";
      for ($i = 0; $i <= $#sphxyzcoords; $i++) { # determine distance of each sphere to center of receptor atoms
          $dist=sqrt((($recavgX-$sphxyzcoords[$i][1])**2)+(($recavgY-$sphxyzcoords[$i][2])**2)+(($recavgZ-$sphxyzcoords[$i][3])**2));       
          $sphxyzcoords[$i][8]=$dist; # Weighted at 1 importance                             
      }
      for ($i = 0; $i <= $#sphxyzcoords; $i++) { # determine distance of each sphere to center of ligand atoms
          $dist=sqrt((($ligavgX-$sphxyzcoords[$i][1])**2)+(($ligavgY-$sphxyzcoords[$i][2])**2)+(($ligavgZ-$sphxyzcoords[$i][3])**2));
          $sphxyzcoords[$i][8]=$sphxyzcoords[$i][8]+$dist; # Weighted at 1 importance
      }
      for ($i = 0; $i <= $#sphxyzcoords; $i++) { # determine distance of each sphere to center of spheres
          $dist=sqrt((($sphavgX-$sphxyzcoords[$i][1])**2)+(($sphavgY-$sphxyzcoords[$i][2])**2)+(($sphavgZ-$sphxyzcoords[$i][3])**2));
          $sphxyzcoords[$i][8]=$sphxyzcoords[$i][8]+$dist; # Weighted at 1 importance
      }
   }
   # Calculate average of the "weighted distances" calculated above to give priority to polar spheres
   $avgweight=0;
   for ($i=0; $i<=$#sphxyzcoords; $i++) {
      $avgweight=$avgweight+$sphxyzcoords[$i][8];
   }
   $avgweight = $avgweight/($#sphxyzcoords+1);
   print "Average weighted-distance for final elimination is: $avgweight\n";
   for ($i=0; $i<=$#sphxyzcoords; $i++) {  # Decrease weighted distance for polar spheres
      if (($sphxyzcoords[$i][7]!=-1)&&($sphxyzcoords[$i][7]!=-2)) {
         # print "Adjusting weight for polar sphere ".($i+1)." from $sphxyzcoords[$i][8]";
         $sphxyzcoords[$i][8]=$sphxyzcoords[$i][8]+($avgweight*$polarbenefit);
         # print " to $sphxyzcoords[$i][8]\n";
      }
   }
   # Begin eliminating spheres based on "weighted distances" calculated above
   @tempsortedarray = sort{$b->[8] <=> $a->[8]} @sphxyzcoords; # Sort by descending weighted distance
   $index = 0;
   while (($#tempsortedarray+1 > $maxspheres) && ($index <= $#tempsortedarray)) { # remove spheres flagged continuous
      if ($tempsortedarray[$index][6]==2) { # flag is only continuous, not crystallographic
         splice(@tempsortedarray, $index, 1);
         $index--;
      }
      $index++;
   }
   while ( $#tempsortedarray+1 > $maxspheres) {  # remove ANY spheres until less than MAXSPHERES
      print "WARNING: removing a crystallographic sphere since ".($#tempsortedarray+1).">$maxspheres max spheres are present\n";
      splice(@tempsortedarray, 0, 1);
   }
   @sphxyzcoords = @tempsortedarray;          # Reinitialize @sphxyzcoords
}

##################### OUTPUT SPHERES TO FILE sph/match1.sph ##############################
$numspheres = $#sphxyzcoords+1;
print "Final number of output spheres is: ".$numspheres."\n";
$numcryst=0;
$numpolar=0;
$numnonpolar=0;
$numneither=0;
for ($i=0; $i<=$#sphxyzcoords; $i++) {
   if (($sphxyzcoords[$i][6] == 1) || ($sphxyzcoords[$i][6]==3)) { $numcryst++; }
   if (($sphxyzcoords[$i][7] != -1) && ($sphxyzcoords[$i][7]!=-2)) { $numpolar++; }
   if ($sphxyzcoords[$i][7] == -2) { $numnonpolar++; }
   if ($sphxyzcoords[$i][7] == -1) { $numneither++; }
}
print "Number of polar spheres: ".$numpolar."\n";
print "Number of nonpolar spheres: ".$numnonpolar."\n";
print "Number of spheres that are near neither polar nor nonpolar receptor atoms: ".$numneither."\n";
print "Final number of spheres that are from crystallographic ligand: ".$numcryst."\n";
# output match1.sph file containing dock coloring header and then sphere coordinates
#$MATCHSPH = "sph/match1.sph";
$MATCHSPH = $ARGV[3]; #specify via commandline
open (MATCHSPH, ">$MATCHSPH");
print (MATCHSPH "DOCK 5.2 ligand_atoms\n");
print (MATCHSPH "positive                       (1)\n");
print (MATCHSPH "negative                       (2)\n");
print (MATCHSPH "acceptor                       (3)\n");
print (MATCHSPH "donor                          (4)\n");
print (MATCHSPH "ester_o                        (5)\n");
print (MATCHSPH "amide_o                        (6)\n");
print (MATCHSPH "neutral                        (7)\n");
print (MATCHSPH "not_neutral                    (8)\n");
print (MATCHSPH "positive_or_donor              (9)\n");
print (MATCHSPH "negative_or_acceptor           (10)\n");
print (MATCHSPH "neutral_or_acceptor_or_donor   (11)\n");
print (MATCHSPH "donacc                         (12)\n");
printf (MATCHSPH "cluster     1   number of spheres in cluster %5.f\n",$numspheres);
$outline=1;
for ($line = 0; $line <= $#sphxyzcoords; $line++ ) {
    printf (MATCHSPH "%5.f%10.5f%10.5f%10.5f%8.3f%5.f 0  0          \n",(9000+$outline),($sphxyzcoords[$line][1]),($sphxyzcoords[$line][2]),($sphxyzcoords[$line][3]),($sphxyzcoords[$line][4]),($sphxyzcoords[$line][5]));
    $outline++;
}
close (MATCHSPH);
