
#################################################################################################################
## Writen by Trent Balius in the Shoichet Lab, UCSF in 2015
#################################################################################################################
import sys

class sphere:
    def __init__(self,index,X,Y,Z,radius,atomnum,critical_cluster,sphere_color):
        self.index            = int(index)
        self.X                = float(X)
        self.Y                = float(Y)
        self.Z                = float(Z)
        self.radius           = float(radius)
        self.atomnum          = int(atomnum)
        self.critical_cluster = int(critical_cluster)
        self.sphere_color     = int(sphere_color)
def byIndex(x, y):
    return cmp(x.index, y.index)

def Are_equal(x,y):
    if (x.X == y.X and x.Y == y.Y and x.Z == y.Z): 
        return True
    return False 

#The number of the atom with which surface point j (second point used to generate the sphere) is associated.
#The critical cluster to which this sphere belongs.
#The sphere color. The color is simply an index into the color table that was specified in the header. Therefore, 1 corresponds to the first color in the header, 2 for the second, etc. 0 corresponds to unlabeled. 

def in_list(val, list):
    for ele in list:
        if Are_equal(val,ele):
           return True
    return False

def remove_val(val,index, list):
    # remove val from list
    list_ele = [] ## list of elements the same as value. 
    for i in range(index+1,len(list)): ## We assume that only element below in the list can be equal.  
                                       ## this is valid because we start at the begin.
         if Are_equal(val,list[i]):    ## we chech does the frist element have dup, 
            list_ele.append(i);        ## does the second, 3rd, and so on. 
                                       ## so, no earlyer element will be a duplicate of the val. 
    count = 0
    for i in list_ele:
        list.pop(i-count)
        count=count+1 # index will change everytime somthing is removed.
    
# this function will remove duplicates from the list
# the duplicates have the same X,Y,Z coordenates.
def remove_dullicates(list):
    index = 0 
    for ele in list:
        remove_val(ele,index, list)
        index = index+1
        #print len(list)
    

#FORMAT: (I5, 3F10.5, F8.3, I5, I2, I3) 
# ccluster is the chosen cluster

def read_sph(filename,ccluster,color):
    sphere_list = []
    insph = open(filename,'r')
    flag_cluster = False # this flag determins if the sphere is writen to list
    
    for line in insph:
       #if line[0:4] == 'DOCK' or line[0:4]=='clus':
       if line[0:4] == 'DOCK':
          continue
       elif line[0:4]=='clus':
          cluster = int(line[7:16] )
          #print cluster
          if ccluster == 'A': 
             flag_cluster = True
          elif int(ccluster) == cluster:
             print(("cluster", int(ccluster), cluster))
             flag_cluster = True
          else:
             flag_cluster = False
       elif not (line[0:5].replace(' ','').isdigit()):
            #print line[0:5], line[0:5].replace(' ','').isdigit()
            print (line)
       else:
          index   = int(line[0:5])
          x       = float(line[5:15])
          y       = float(line[15:25])
          z       = float(line[25:35])
          r       = float(line[35:43])
          if (r == 0.0): 
              r = 0.5
              print ("radius of 0.0 detected.  changed to 0.5.")
          
          atomnum = int(line[43:48])
          #print '\'' + line[48:50] + '\'' 
          if (line[48:50] != '  '):
              clust = int(line[48:50])
          else: 
              clust = 0
          if (line[50:53] != '   '):
             col = int(line[50:53])
          else:
             col = 0

          if ( color == 'A' or int(color) == col):
             flag_color = True
          else:
             flag_color = False
          
          tmp_sphere = sphere(index,x,y,z,r,atomnum,clust,col)

          if (flag_cluster and flag_color): 
          # only put sphere on list if it is in a cluster of instrested
          # and if the color is the same
             sphere_list.append(tmp_sphere)

    #sphere_list.append(tmp_sphere)
    #sphere_list.sort(byIndex)
    #sphere_list.sort(key = lambda a,b: cmp(a.index, b.index))
    sphere_list.sort(key = lambda a: a.index)
    #sphere_list.sort(key = byIndex)
    ## remove duplicates:
    remove_dullicates(sphere_list)

    if len(sphere_list) == 0:
       print ("there is a problem")

    return sphere_list

def write_sph(filename,spheres):
    outsph = open(filename,'w')
    outsph.write("DOCK spheres generated from read_write_sph.py\n")
    outsph.write("cluster     1   number of spheres in cluster %3d\n" % len(spheres))
    for i in range(len(spheres)):
        outsph.write("%5d%10.5f%10.5f%10.5f%8.3f%5d%2d%3d\n" %
                     (i+1,round(spheres[i].X,3),
                          round(spheres[i].Y,3),
                          round(spheres[i].Z,3),
                          #spheres[i].radius,i+1,0,0) )
                          #spheres[i].radius,spheres[i].atomnum,0,0) )
                          spheres[i].radius,spheres[i].atomnum,0,spheres[i].sphere_color) )

