#!/usr/bin/perl -w
# makebox.pl
# This is a script to make box
# file created from input spheres file sph/match2.sph
# Usage: in project directory with subdirectories sph and grids, "makebox.pl projectdirectory"
# Pseudocode:
#  Read in Sphere file (sph/match.sph)
#  Trim sphere file into array of only 3D coordinates
#  Start loop that checks box size less than maximum Chemgrid can handle (~1.7 million points)
#    initialize variables, arrays, margin, etc.
#    read all 3D coordinates and create custom box tailored to margin around them
#    if box is too large and eliminate farthest point from the spheres
#  If box is smaller than 50000 points, then increase box size based on receptor atoms
#  Print out results to file /grids/box
# Written by Austin N. Kirschner 12/2002, modified 7/9/2003 ANK;

#$projectdir = $ARGV[0];
$SPH = $ARGV[0]; #actually let user specify things (sphere file)
$REC = $ARGV[1]; #" (receptor pdb file)
$BOX = $ARGV[2]; #" (output box file)
$M   = $ARGV[3]; # Margin value should be at least 6.2 angstroms 

# Read in sphere file (sph/match2.sph)
#$SPH = "$projectdir/sph/match2.sph";
open (SPH, "<$SPH") or die "My error, cannot open SPH: $!\n";
while (<SPH>) {                     # read in SPH lines
        @spharray = split;
        push @sphitems, [ @spharray ];
}
$start = 0;   # flag for starting reading coordinates
for ($i=0; $i<=$#sphitems; $i++) {
   if ($start == 1) {
      push @xyzcoords, [ @{ $sphitems[$i] } [ 1..3 ] ];  # only 3D coordinates
   }
   if ($sphitems[$i][0] eq "cluster") {$start = 1;} # coordinates start next 
}

# start creating box by defining arbitrary parameters
#$M = 10;             # Margin value should be at least 6.2 angstroms 
$PPA = 3;            # Points Per Angstrom
$MAXBOX = 1700000;   # Maximum number of points in box
$MINBOX = 50000;   # Minimum number of points in box #down to 50000 from 150000, rgc 2012
$MAXSIDE = 50.000;   # Maximum length of a box side (determined by 3*side for distmap, limited by 150 for dock)
$CVP = $MAXBOX + 1;  # Current Value of Points is total volume of box in absolute points
$I = 1;              # Increment for adjusting size of box if too small

while ( $CVP > $MAXBOX ) {  # creates initial box and handles if box is too large (input spheres too far apart)
   #Create box using diagonal corners
   @D1 = ([ 0, 0, 0]);  # Three-dimensional Diagonal corner of box
   @D2 = ([ 0, 0, 0]);  # Three-dimensional Diagonal corner of box

   #Initialize box corners
   $initX1 = $xyzcoords[0][0] - $M;
   $initY1 = $xyzcoords[0][1] - $M;
   $initZ1 = $xyzcoords[0][2] - $M;
   @D1 = ( $initX1, $initY1, $initZ1 );
   $initX2 = $xyzcoords[0][0] + $M;
   $initY2 = $xyzcoords[0][1] + $M;
   $initZ2 = $xyzcoords[0][2] + $M;
   @D2 = ( $initX2, $initY2, $initZ2 );

   #Read all xyz coords and resize box accordingly
   for ($line = 1; $line <= $#xyzcoords; $line++) {
      if ( ($D1[0]+$M<=$xyzcoords[$line][0]) && ($xyzcoords[$line][0]<=$D2[0]-$M) ) {;}
      elsif ($xyzcoords[$line][0]<$D1[0]+$M) { print "X- "; $D1[0]=$xyzcoords[$line][0] - $M; }
      elsif ($D2[0]-$M<$xyzcoords[$line][0]) { print "X+ "; $D2[0]=$xyzcoords[$line][0] + $M; }

      if ( ($D1[1]+$M<=$xyzcoords[$line][1]) && ($xyzcoords[$line][1]<=$D2[1]-$M) ) {;}
      elsif ($xyzcoords[$line][1]<$D1[1]+$M) { print "Y- "; $D1[1]=$xyzcoords[$line][1] - $M; }
      elsif ($D2[1]-$M<$xyzcoords[$line][1]) { print "Y+ "; $D2[1]=$xyzcoords[$line][1] + $M; }

      if ( ($D1[2]+$M<=$xyzcoords[$line][2]) && ($xyzcoords[$line][2]<=$D2[2]-$M) ) {;}
      elsif ($xyzcoords[$line][2]<$D1[2]+$M) { print "Z- "; $D1[2]=$xyzcoords[$line][2] - $M; }
      elsif ($D2[2]-$M<$xyzcoords[$line][2]) { print "Z+ "; $D2[2]=$xyzcoords[$line][2] + $M; }
   }

   # Calculate center of box coordinates
   $centX = ($D2[0]+$D1[0])/2;
   $centY = ($D2[1]+$D1[1])/2;
   $centZ = ($D2[2]+$D1[2])/2;

   $CVP = ( abs($D2[0]-$D1[0]) * abs($D2[1]-$D1[1]) * abs($D2[2]-$D1[2]) * $PPA * $PPA * $PPA );
   # print "Box has: ".$CVP." 3D points in it\n"; 

   if ( $CVP > $MAXBOX ) {
      print "Volume points of box is too large >".$MAXBOX.", so replacing farthest point with current center of box\n";
      # calculate distances of each point in three-dimenstional space to center of box
      @xyzdistances = ([ 0 ]);
      $largestdistline = 0;
      for ($line = 0; $line <= $#xyzcoords; $line++) {
         $xyzdistances[$line] = sqrt((($centX-$xyzcoords[$line][0])**2)+(($centY-$xyzcoords[$line][1])**2)+(($centZ-$xyzcoords[$line][2])**2));
         # find the line that has the point farthest from the center of the box
         if ( $xyzdistances[$line] > $xyzdistances[$largestdistline] ) {
            $largestdistline = $line;
         }
      }
      $xyzcoords[$largestdistline][0] = $centX;  # replacing farthest points with center of box points
      $xyzcoords[$largestdistline][1] = $centY;  # this lessens the effect of points too far away
      $xyzcoords[$largestdistline][2] = $centZ;  # but too many too far may make skewed box
   }
}

if ($CVP < $MINBOX ) {   # handles if box is too small, increase x y and z by 1.0 until large enough
   print "\nAfter creating box around spheres, box is too small (".$CVP." < ".$MINBOX."), so increasing its size\n";
   # Calculate center of receptor coords                                     
   #$REC = "$projectdir/rec.pdb";    # Using rec.crg includes the spheres appended to the end of the file
   open (REC, "<$REC") or die "My error, cannot open REC.PDB: $!\n";
   while (<REC>) {                     # read in REC lines
      @recarray = split;
      push @recitems, [ @recarray ];
   }
   for ($x=0; $x<=$#recitems; $x++) {
      if (($recitems[$x][0] eq "ATOM") || ($recitems[$x][0] eq "HETATM")) {
        #push @recxyzcoords, [ @{ $recitems[$x] } [ 5..7 ] ];
        if ($recitems[$x][5] =~ m/\./) {
            push @recxyzcoords, [ @{ $recitems[$x] } [ 5..7 ] ];
        } else {
            push @recxyzcoords, [ @{ $recitems[$x] } [ 6..8 ] ];
        }
      }
   }
   $avgX = $recxyzcoords[0][0];                                                         
   $avgY = $recxyzcoords[0][1];                                                         
   $avgZ = $recxyzcoords[0][2];                                                         
   for ($line = 1; $line <= $#recxyzcoords; $line++) {                                  
      $avgX = $avgX + $recxyzcoords[$line][0];                                          
      $avgY = $avgY + $recxyzcoords[$line][1];                                          
      $avgZ = $avgZ + $recxyzcoords[$line][2];                                          
   }                                                                                 
   $avgX = $avgX/$#recxyzcoords+1;                                                        
   $avgY = $avgY/$#recxyzcoords+1;                                                        
   $avgZ = $avgZ/$#recxyzcoords+1;                                                        
   print "Center of receptor atoms is (x y z): ".$avgX." ".$avgY." ".$avgZ."\n";
   $count = 1;  # just in case it cannot increase box on all 3 sides, it should escape the while loop below
   while (($CVP < $MINBOX ) && ($count <= 200)) {
      $count++;
      $centX = abs($D2[0]+$D1[0])/2;
      $centY = abs($D2[1]+$D1[1])/2;
      $centZ = abs($D2[2]+$D1[2])/2;
      print "Center of current box is (x y z): ".$centX." ".$centY." ".$centZ."\n";
      print "Diagonals of current box are (D1: x y z, D2: x y z): ".$D1[0]." ".$D1[1]." ".$D1[2]." ".$D2[0]." ".$D2[1]." ".$D2[2]."\n";
      @centD1Xwall = ( $D1[0], ($D1[1]+$D2[1])/2, ($D1[2]+$D2[2])/2);
      @centD2Xwall = ( $D2[0], ($D1[1]+$D2[1])/2, ($D1[2]+$D2[2])/2);
      @centD1Ywall = ( ($D1[0]+$D2[0])/2, $D1[1], ($D1[2]+$D2[2])/2);
      @centD2Ywall = ( ($D1[0]+$D2[0])/2, $D2[1], ($D1[2]+$D2[2])/2);
      @centD1Zwall = ( ($D1[0]+$D2[0])/2, ($D1[1]+$D2[1])/2, $D1[2]);
      @centD2Zwall = ( ($D1[0]+$D2[0])/2, ($D1[1]+$D2[1])/2, $D2[2]);
      # print "Center of box wall D1Xwall: ".$centD1Xwall[0].",".$centD1Xwall[1].",".$centD1Xwall[2]."\n";
      # print "Center of box wall D2Xwall: ".$centD2Xwall[0].",".$centD2Xwall[1].",".$centD2Xwall[2]."\n";
      # print "Center of box wall D1Ywall: ".$centD1Ywall[0].",".$centD1Ywall[1].",".$centD1Ywall[2]."\n";
      # print "Center of box wall D2Ywall: ".$centD2Ywall[0].",".$centD2Ywall[1].",".$centD2Ywall[2]."\n";
      # print "Center of box wall D1Zwall: ".$centD1Zwall[0].",".$centD1Zwall[1].",".$centD1Zwall[2]."\n";
      # print "Center of box wall D2Zwall: ".$centD2Zwall[0].",".$centD2Zwall[1].",".$centD2Zwall[2]."\n";
      $D1Xflag = 0;
      $D2Xflag = 0;
      $D1Yflag = 0;
      $D2Yflag = 0;
      $D1Zflag = 0;
      $D2Zflag = 0;
      $i = 0;
      while ( ($i<=$#recxyzcoords) && (($D1Xflag==0)||($D2Xflag==0)||($D1Yflag==0)||($D2Yflag==0)||($D1Zflag==0)||($D2Zflag==0)) ) {
         if(($D1Xflag==0)&&(sqrt((($recxyzcoords[$i][0]-$centD1Xwall[0])**2)+(($recxyzcoords[$i][1]-$centD1Xwall[1])**2)+(($recxyzcoords[$i][2]-$centD1Xwall[2])**2))<$M/2)){$D1Xflag=1;}
         if(($D2Xflag==0)&&(sqrt((($recxyzcoords[$i][0]-$centD2Xwall[0])**2)+(($recxyzcoords[$i][1]-$centD2Xwall[1])**2)+(($recxyzcoords[$i][2]-$centD2Xwall[2])**2))<$M/2)){$D2Xflag=1;}
         if(($D1Yflag==0)&&(sqrt((($recxyzcoords[$i][0]-$centD1Ywall[0])**2)+(($recxyzcoords[$i][1]-$centD1Ywall[1])**2)+(($recxyzcoords[$i][2]-$centD1Ywall[2])**2))<$M/2)){$D1Yflag=1;}
         if(($D2Yflag==0)&&(sqrt((($recxyzcoords[$i][0]-$centD2Ywall[0])**2)+(($recxyzcoords[$i][1]-$centD2Ywall[1])**2)+(($recxyzcoords[$i][2]-$centD2Ywall[2])**2))<$M/2)){$D2Yflag=1;}
         if(($D1Zflag==0)&&(sqrt((($recxyzcoords[$i][0]-$centD1Zwall[0])**2)+(($recxyzcoords[$i][1]-$centD1Zwall[1])**2)+(($recxyzcoords[$i][2]-$centD1Zwall[2])**2))<$M/2)){$D1Zflag=1;}
         if(($D2Zflag==0)&&(sqrt((($recxyzcoords[$i][0]-$centD2Zwall[0])**2)+(($recxyzcoords[$i][1]-$centD2Zwall[1])**2)+(($recxyzcoords[$i][2]-$centD2Zwall[2])**2))<$M/2)){$D2Zflag=1;}
         $i++;
      }
      $i=0;
      $closestD1X = $M;
      $closestD2X = $M;
      $closestD1Y = $M;
      $closestD2Y = $M;
      $closestD1Z = $M;
      $closestD2Z = $M;
      while(($i<=$#recxyzcoords)&&(($D1Xflag==0)||($D2Xflag==0)||($D1Yflag==0)||($D2Yflag==0)||($D1Zflag==0)||($D2Zflag==0))){
         if(sqrt((($recxyzcoords[$i][0]-$centD1Xwall[0])**2)+(($recxyzcoords[$i][1]-$centD1Xwall[1])**2)+(($recxyzcoords[$i][2]-$centD1Xwall[2])**2))<$closestD1X){$closestD1X=sqrt((($recxyzcoords[$i][0]-$centD1Xwall[0])**2)+(($recxyzcoords[$i][1]-$centD1Xwall[1])**2)+(($recxyzcoords[$i][2]-$centD1Xwall[2])**2))}
         if(sqrt((($recxyzcoords[$i][0]-$centD2Xwall[0])**2)+(($recxyzcoords[$i][1]-$centD2Xwall[1])**2)+(($recxyzcoords[$i][2]-$centD2Xwall[2])**2))<$closestD2X){$closestD2X=sqrt((($recxyzcoords[$i][0]-$centD2Xwall[0])**2)+(($recxyzcoords[$i][1]-$centD2Xwall[1])**2)+(($recxyzcoords[$i][2]-$centD2Xwall[2])**2))}
         if(sqrt((($recxyzcoords[$i][0]-$centD1Ywall[0])**2)+(($recxyzcoords[$i][1]-$centD1Ywall[1])**2)+(($recxyzcoords[$i][2]-$centD1Ywall[2])**2))<$closestD1Y){$closestD1Y=sqrt((($recxyzcoords[$i][0]-$centD1Ywall[0])**2)+(($recxyzcoords[$i][1]-$centD1Ywall[1])**2)+(($recxyzcoords[$i][2]-$centD1Ywall[2])**2))}
         if(sqrt((($recxyzcoords[$i][0]-$centD2Ywall[0])**2)+(($recxyzcoords[$i][1]-$centD2Ywall[1])**2)+(($recxyzcoords[$i][2]-$centD2Ywall[2])**2))<$closestD2Y){$closestD2Y=sqrt((($recxyzcoords[$i][0]-$centD2Ywall[0])**2)+(($recxyzcoords[$i][1]-$centD2Ywall[1])**2)+(($recxyzcoords[$i][2]-$centD2Ywall[2])**2))}
         if(sqrt((($recxyzcoords[$i][0]-$centD1Zwall[0])**2)+(($recxyzcoords[$i][1]-$centD1Zwall[1])**2)+(($recxyzcoords[$i][2]-$centD1Zwall[2])**2))<$closestD1Z){$closestD1Z=sqrt((($recxyzcoords[$i][0]-$centD1Zwall[0])**2)+(($recxyzcoords[$i][1]-$centD1Zwall[1])**2)+(($recxyzcoords[$i][2]-$centD1Zwall[2])**2))}
         if(sqrt((($recxyzcoords[$i][0]-$centD2Zwall[0])**2)+(($recxyzcoords[$i][1]-$centD2Zwall[1])**2)+(($recxyzcoords[$i][2]-$centD2Zwall[2])**2))<$closestD2Z){$closestD2Z=sqrt((($recxyzcoords[$i][0]-$centD2Zwall[0])**2)+(($recxyzcoords[$i][1]-$centD2Zwall[1])**2)+(($recxyzcoords[$i][2]-$centD2Zwall[2])**2))}
         $i++;
      }
      if ($D1Xflag==1) { $D1[0]=$D1[0]-$I; print "D1X- ";}
      if ($D2Xflag==1) { $D2[0]=$D2[0]+$I; print "D2X+ ";}
      if ($D1Yflag==1) { $D1[1]=$D1[1]-$I; print "D1Y- ";}
      if ($D2Yflag==1) { $D2[1]=$D2[1]+$I; print "D2Y+ ";}
      if ($D1Zflag==1) { $D1[2]=$D1[2]-$I; print "D1Z- ";}
      if ($D2Zflag==1) { $D2[2]=$D2[2]+$I; print "D2Z+ ";}
      if (($D1Xflag==0)&&($closestD1X==$M)) { $D1[0]=$D1[0]+($I); print "D1X+ ";}
      if (($D2Xflag==0)&&($closestD2X==$M)) { $D2[0]=$D2[0]-($I); print "D2X- ";}
      if (($D1Yflag==0)&&($closestD1Y==$M)) { $D1[1]=$D1[1]+($I); print "D1Y+ ";}
      if (($D2Yflag==0)&&($closestD2Y==$M)) { $D2[1]=$D2[1]-($I); print "D2Y- ";}
      if (($D1Zflag==0)&&($closestD1Z==$M)) { $D1[2]=$D1[2]+($I); print "D1Z+ ";}
      if (($D2Zflag==0)&&($closestD2Z==$M)) { $D2[2]=$D2[2]-($I); print "D2Z- ";}
      print "\n";

      # Check if any sides are >50 A (becomes side*3 for distmap, and maximum set by current dock compilation)
      while (abs($D2[0]-$D1[0]) > $MAXSIDE) {
         $D2[0]=$D2[0]-0.00025;
         $D1[0]=$D1[0]+0.00025;
      }
      while (abs($D2[1]-$D1[1]) > $MAXSIDE) {
         $D2[1]=$D2[1]-0.00025;
         $D1[1]=$D1[1]+0.00025;
      }
      while (abs($D2[2]-$D1[2]) > $MAXSIDE) {
         $D2[2]=$D2[2]-0.00025;
         $D1[2]=$D1[2]+0.00025;
      }

      # Check if adjustments to all three sides (x y z) are enough to make box larger
      $CVP = ( abs($D2[0]-$D1[0]) * abs($D2[1]-$D1[1]) * abs($D2[2]-$D1[2]) * $PPA * $PPA * $PPA );
      print "**Total grid points in box is now: ".$CVP."\n";
      print "**BOX SIZE is now: ".abs($D2[0]-$D1[0])." ".abs($D2[1]-$D1[1])." ".abs($D2[2]-$D1[2])."\n";
   }
   print "Center of receptor atoms is (x y z): ".$avgX." ".$avgY." ".$avgZ."\n";
   $centX = ($D2[0]+$D1[0])/2;
   $centY = ($D2[1]+$D1[1])/2;
   $centZ = ($D2[2]+$D1[2])/2;
   print "Center of current box is (x y z): ".$centX." ".$centY." ".$centZ."\n";
   print "Diagonals of current box is (D1: x y z, D2: x y z): ".$D1[0]." ".$D1[1]." ".$D1[2]." ".$D2[0]." ".$D2[1]." ".$D2[2]."\n";
}
$D1[0] = sprintf ("%8.3f", $D1[0]);
$D1[1] = sprintf ("%8.3f", $D1[1]);
$D1[2] = sprintf ("%8.3f", $D1[2]);
$D2[0] = sprintf ("%8.3f", $D2[0]);
$D2[1] = sprintf ("%8.3f", $D2[1]);
$D2[2] = sprintf ("%8.3f", $D2[2]);
print "Limit diagonals to 8.3f :".$D1[0]." ".$D1[1]." ".$D1[2]." ".$D2[0]." ".$D2[1]." ".$D2[2]."\n";

$centX = ($D2[0]+$D1[0])/2;
$centY = ($D2[1]+$D1[1])/2;
$centZ = ($D2[2]+$D1[2])/2;
close (SPH);  # done reading and managing coordinates as input

# output box file containing center of box coordinates, lengths of box, diagonal coordinates
#$BOX = "$projectdir/grids/box";
open (BOX, ">$BOX");
print (BOX "HEADER    CORNERS OF BOX ");
printf (BOX "%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f",$D1[0],$D1[1],$D1[2],$D2[0],$D2[1],$D2[2]);
print (BOX "\nREMARK    CENTER (X Y Z) ");
printf (BOX "%8.3f%8.3f%8.3f",$centX,$centY,$centZ);
print (BOX "\nREMARK    DIMENSIONS (X Y Z) ");
printf (BOX "%8.3f%8.3f%8.3f",abs($D2[0]-$D1[0]),abs($D2[1]-$D1[1]),abs($D2[2]-$D1[2]));
print (BOX "\nATOM      1  DUA BOX     1    ");
printf (BOX "%8.3f%8.3f%8.3f",$D1[0],$D1[1],$D1[2]);
print (BOX "\nATOM      2  DUB BOX     1    ");
printf (BOX "%8.3f%8.3f%8.3f",$D2[0],$D1[1],$D1[2]);
print (BOX "\nATOM      3  DUC BOX     1    ");
printf (BOX "%8.3f%8.3f%8.3f",$D2[0],$D1[1],$D2[2]);
print (BOX "\nATOM      4  DUD BOX     1    ");
printf (BOX "%8.3f%8.3f%8.3f",$D1[0],$D1[1],$D2[2]);
print (BOX "\nATOM      5  DUE BOX     1    ");
printf (BOX "%8.3f%8.3f%8.3f",$D1[0],$D2[1],$D1[2]);
print (BOX "\nATOM      6  DUF BOX     1    ");
printf (BOX "%8.3f%8.3f%8.3f",$D2[0],$D2[1],$D1[2]);
print (BOX "\nATOM      7  DUG BOX     1    ");
printf (BOX "%8.3f%8.3f%8.3f",$D2[0],$D2[1],$D2[2]);
print (BOX "\nATOM      8  DUH BOX     1    ");
printf (BOX "%8.3f%8.3f%8.3f",$D1[0],$D2[1],$D2[2]);
print (BOX "\nCONECT    1    2    4    5\n");
print (BOX "CONECT    2    1    3    6\n");
print (BOX "CONECT    3    2    4    7\n");
print (BOX "CONECT    4    1    3    8\n");
print (BOX "CONECT    5    1    6    8\n");
print (BOX "CONECT    6    2    5    7\n");
print (BOX "CONECT    7    3    6    8\n");
print (BOX "CONECT    8    4    5    7\n");
close (BOX);
# print "Completed makebox.pl script.\n";
