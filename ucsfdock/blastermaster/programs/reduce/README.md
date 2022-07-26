For more on reduce, see:

http://kinemage.biochem.duke.edu/software/readme.reduce.html

Static copy of page copied here:


README.reduce	Home Page: Richardsons' Laboratory	J. Michael Word - 6/2001


Contents:

Introduction
New Features
Running Reduce
Command Line Help
Fixing an Orientation
Cliques
Program Information
Other Sources of Help


Introduction
This describes the pre-compiled executable for Reduce version 2.14, developed by J. Michael Word. Reduce is a (UNIX/Linux/Windows) program for adding hydrogens to a Protein DataBank (PDB) molecular structure file. It was developed at Duke University in the lab of David and Jane Richardson. Reduce is described in:
Word, et. al. (1999) Asparagine and Glutamine: Using Hydrogen Atom Contacts in the Choice of Side-chain Amide Orientation, J. Mol. Biol. 285, 1733-1747.

Both proteins and nucleic acids can have hydrogens added. HET groups can also be processed as long as the atom connectivity is provided. A slightly modified version of the connectivity table provided by the PDB is included (as a compressed text file). Reduce will attempt to look for this information, in the standard location, "/usr/local/reduce_het_dict.txt". If you place this file somewhere other than /usr/local/ or want to refer to a different connection file for a particular project, use the -DB command line flag to point to the new connection file, or set an environmental variable prior to running Reduce:

setenv REDUCE_HET_DICT connectionfilename


New Features
5/31/01 SegmentID to Chain Mapping

It is now possible to process a PDBfile where multiple chains are identified by their segment ID and the chain ID is missing (certain xtallographic model building and refinement packages do this). A new flag -SEGIDmap "seg,c..." allows you to define a set of segid->chainid mappings. Reduce will then add the chainids and process the PDBfile. In the string of seg,c,seg,c mapping info, use "_" to represent a blank.



Running Reduce
In most circumstances, the recommended command when using Reduce to add hydrogens to a PDB file and standardize the bond lengths of existing hydrogens is:

reduce -build coordfile.pdb > coordfileH.pdb
which includes the optimization of adjustable groups (OH, SH, NH3+, Met-CH3, and Asn, Gln and His sidechain orientation). When speed is important, the -build option can be dropped; hydrogens will still be added, but not His side-chain NH hydrogens, and side-chains will not be flipped. For even greater speed, but even less accuracy, adding -nooh and -noadj will skip the OH and SH hydrogens and eliminate optimization altogether. Input is from the specified PDB format coordinate file and the new, updated PDB coordinates are written to "standard output", here redirected to a file with the '>' symbol.

Disulfides, covalent modifications, and connection of the ribose-phosphate nucleic acid backbone, are recognized and any hydrogens eliminated by bonding are skipped. When an amino acid main-chain nitrogen is not connected to the preceeding residue or some other group, Reduce treats it as the N-terminus and constructs an NH3+ only if the residue number is less than or equal to an ajustable limit (1, by default). Otherwise, it considers the residue to be the observable beginning of an actually-connected fragment and does not protonate the nitrogen. Reduce does not protonate carboxylates (including the C-terminus) bacause it does not specifically consider pH, instead modeling a neutral environment.

Hydrogens are positioned with respect to the covalently bonded neighbors and these are identified by name. Non-standard atom names are the primary cause of missing or misplaced hydrogens. If Reduce tries to process a file which contains hydrogens with non-standard names, the existing hydrogens may not be recognized and may interfere with the generation of new hydrogens. The solution may be to remove existing hydrogens before further processing.

Hydrogens can be removed from a pdb format file with Reduce.

reduce -trim 1abcH > 1abc
This can be used, for example, to update the orientation of Asn/Gln/His side chains where the H atoms are not wanted; first build the hydrogens and then trim them back out. Trimming can occasionally be fooled if a hydrogen has been given a non-standard name. The most common example of this comes from left-justified atom names: gamma hydrogens masquerade as mercury atoms! In this case, manual editing may be required.



Internal program help text
The following brief description of the command line flags is displayed with the -h flag:

$ reduce -h

reduce: version 2.15 10/4/01, Copyright 1997-2001, J. Michael Word
arguments: [-flags] filename or -

Adds hydrogens to a PDB format file and writes to standard output.
(note: By default, HIS sidechain NH protons are not added. See -BUILD)

Flags:
-Trim             remove (rather than add) hydrogens

-NOOH             remove hydrogens on OH and SH groups
-OH               add hydrogens on OH and SH groups (default)

-HIS              create NH hydrogens on HIS rings
-FLIPs            allow complete ASN, GLN and HIS sidechains to flip
                        (usually used with -HIS)
-NOHETh           do not attempt to add NH proton on Het groups
-ROTNH3           allow lysine NH3 to rotate (default)
-NOROTNH3         do not allow lysine NH3 to rotate
-ROTEXist         allow existing rotatable groups (OH, SH, Met-CH3) to rotate
-ROTEXOH          allow existing OH & SH groups to rotate
-ALLMEthyls       allow all methyl groups to rotate
-ONLYA            only adjust 'A' conformations (default)
-ALLALT           process adjustments for all conformations
-NOROTMET         do not rotate methionine methyl groups
-NOADJust         do not process any rot or flip adjustments

-BUILD            add H, including His sc NH, then rotate and flip groups
                  (except for pre-existing methionine methyl hydrogens)
                  (same as: -OH -ROTEXOH -HIS -FLIP)

-Keep             keep bond lengths as found
-NBonds#          remove dots if cause within n bonds (default=3)
-Model#           which model to process (default=1)
-Nterm#           max number of nterm residue (default=1)
-DENSity#.#       dot density (in dots/A^2) for VDW calculations (default=16)
-RADius#.#        probe radius (in A) for VDW calculations (default=0)
-OCCcuttoff#.#    occupancy cutoff for adjustments (default=0.01)
-H2OBcuttoff#.#   B-factor  cutoff for water atoms (default=40)
-H2OOCCcuttoff#.# occupancy cutoff for water atoms (default=0.66)
-PENalty#.#       fraction of std. bias towards original orientation (default=1)
-HBREGcuttoff#.#  over this gap regular HBonds bump (default=0.6)
-HBCHargedcut#.#  over this gap charged HBonds bump (default=0.4)
-BADBumpcut#.#    at this gap a bump is 'bad' (default=0.4)
-SEGIDmap "seg,c..."  assign chainID based on segment identifier field
-Xplor            use Xplor conventions for naming polar hydrogens
-NOCon            drop conect records
-LIMIT#           max num iter. for exhaustive search (default=100000)
-NOTICKs          do not display the set orientation ticker during processing
-SHOWSCore        display scores for each orientation considered during processing
-FIX "filename"   if given, file specifies orientations for adjustable groups
-DB "filename"    file to search for het info
                        (default="/usr/local/reduce_het_dict.txt")
note: can also redirect with unix environment variable: REDUCE_HET_DICT

-Quiet            do not write extra info to the console
-REFerence        display citation reference
-Help             more extensive description of command line arguments


Fixing an Orientation
At times it is useful to control the flip state or rotation angle of an adjustable group when adding hydrogens, either because the correct orientation has already been established, allowing the optimization time to be reduced, or because a non-optimal orientation is sought.

One of the command line flags (-fix myfile.txt) takes a file containing information about which conformation to set for one or more adjustable groups. The colon delimited format is similar to the orientation data that Reduce prints in the header file:

action:residueID:comment
(one line for each group to be fixed) and because spacing matters in the residue identifier string, the easiest way to produce this file is to copy and edit USER MOD records from Reduce output. The action can be one of three kinds, depending on residue type: O to leave in the original orientation, F to flip the orientation, and R# to rotate a dihedral to an angle of #deg. Using either O or F with His sidechains allows the protonation state to vary; to specify a particular orientation and protonation state use F# where # is the number of the state (1, 2 or 3 for the original orientation with H (1) only on NE2, (2) only on ND1, or (3) doubly protonated; 4-6 for the corresponding three flipped states).



Cliques
The current version of Reduce uses brute-force enumeration to optimize the conformations of ajustable groups. If a 'clique' of adjustible groups is too large (> ≈7) this sort of search technique is inadequate--the enumeration will be abandoned and these groups will be left in their original conformations. The cuttoff point is based on the total number of permutations, which the user can control with the -limit# option. Although we are considering more powerful search techniques for these situations, some work-around strategies have been developed.

First check to see if distinct chainIDs are provided for each chain. Reduce does not support files which specify chain information only in the segID field and can get confused.

Examination of the clique may reveal that the orientations of one or more groups are obvious; for instance, they may interact with obligate H-bond donors or acceptors. By fixing the orientation of these groups (as described above), the total number of permutations is reduced. This is especially effective if it breaks the clique into smaller sum-cliques or singletons.

An alternative way to break up cliques is to rotate all the methionine CH3s and lysine/N-terminus NH3+s in an initial pass, then keep them fixed in a second pass.

reduce -nooh inputfile | reduce -build -norotmet -norotnh3 - > outputfileH
The single dash towards the end of the command line tells Reduce to read data piped ('|') from the first pass rather than from a file. A few NH3+ H-bonds may have inferior geometry with this two pass approach but the result is otherwise comparable to using -build alone and can be combined with the previous approach, if neccessary. With this technique, unusual cliques requiring many hours to process have been converted into several smaller problems which were all solved in a matter of minutes.



Program Information
The latest version of Reduce is available at:    http://kinemage.biochem.duke.edu

Version 2 of Reduce is a C++ program and written on a Silicon Graphics workstation. For low level processing of PDB file records it uses Pdb++, our modified version of a class developed as part of Midas by researchers at the University of California, San Francisco.

If you use Reduce, we would appreciate any comments you send our way.
e-mail: 
J. Michael Word or 
David C. Richardson
URL: http://kinemage.biochem.duke.edu
Biochemistry Department
Duke University
Durham, NC USA 27710

