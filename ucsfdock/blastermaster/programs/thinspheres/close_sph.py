
#################################################################################################################
## Writen by Trent Balius in the Shoichet Lab, UCSF in 2015
## modifed in 2016
## modified in 4/2019 by Reed Stein to add trim_spheres function - this is to reduce redundancy of spheres
# because the sphere layer needs to be uniform, but QNIFFT crashes with large atom numbers
# Increase MS density to increase number of spheres for small thin sphere radii 
#################################################################################################################
import sys
import sph_lib as sph
#import mol2 
import pdb_lib as pdb
#import numpy as np

def trim_sph(sph_list, sph_rad):
     
    for i in range(len(sph_list)-1):
        #if sph_list[i][1] == True:
        if (sph_list[i][1]):
            for j in range(i+1, len(sph_list)):
                #if sph_list[j][1] == True:
                if (sph_list[j][1]):
                    dist = (sph_list[i][0].X - sph_list[j][0].X)**2 + (sph_list[i][0].Y - sph_list[j][0].Y)**2 + (sph_list[i][0].Z - sph_list[j][0].Z)**2
                    if dist <= (sph_rad**2.0) / 2.0:
                       sph_list[j][1] = False

    final_sph_list = []
    for sph in sph_list:
        #if sph[1] == True:
        if (sph[1]):
            final_sph_list.append(sph[0])

    return(final_sph_list)

def distance_sph_pdb(sphs,pdbatoms,dt): 

    sphlist = []
    for sph in sphs:
        for atom in pdbatoms:
            d2 = (atom.X - sph.X)**2 + (atom.Y - sph.Y)**2 + (atom.Z - sph.Z)**2
            if d2 < float(dt)**2.0: 
                sphlist.append([sph, True])
                break
           
    return sphlist

def main():
    if len(sys.argv) != 6: # if no input
       print ("ERROR: there need to be 5 inputs: sph inputfilename, pdb inputfilename, outputfilename, distance, sphere_radius")
       return

    fileinputsph = sys.argv[1]
    fileinputpdb = sys.argv[2]
    fileoutput   = sys.argv[3]
    distance     = float(sys.argv[4])
    sph_rad      = float(sys.argv[5])

    print(('input_sph =' + fileinputsph))
    print(('input_pdb =' + fileinputpdb))
    print(('output =' + fileoutput))
    print(('distance = %6.3f'%distance))
    print(('Sphere radius = %6.3f'%sph_rad))

    list = sph.read_sph(fileinputsph,"A","A")
    pdblist = pdb.read_pdb(fileinputpdb)    
    list2 = distance_sph_pdb(list, pdblist,distance)
    final_list = trim_sph(list2, sph_rad)
    sph.write_sph(fileoutput,final_list)

main()
