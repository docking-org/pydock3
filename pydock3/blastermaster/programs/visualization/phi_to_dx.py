import struct
import array
import sys
import string
import os
import math
import copy

# Written by Reed Stein, Shoichet lab, 2016
# took parts from blastermaster and Trent Balius' GIST tools. 
# Updated for python3 by Brendan Hall


class phi(object):

    def __init__(self, phiFileName=False, is64=False, gridSizes=(193,)):
        '''reads the phi file from disk'''
        if gridSizes is None:
            gridSizes = (193,)
        if gridSizes[0] is None:
            gridSizes = (193,)
        self.oldmid = [0., 0., 0.]
        self.__minsmaxs = None
        self.__boundaries = None
        if phiFileName:  # otherwise just creating an empty phi map for writing
            for gridSize in gridSizes:
                try:
                    phiFile = open(phiFileName, 'rb')  # b is for binary, r is for read
                    tempArray = array.array('f')
                    junk = struct.unpack('4s', phiFile.read(4))
                    check = struct.unpack('4s', phiFile.read(4))[0].decode('utf-8')
                    if check == "now ":  # this changed, but this is now correct
                        #print "32bit phimap"
                        pass
                    else:
                        #print "64bit phimap"
                        is64 = True
                    if not is64:
                        temptop = struct.unpack('16s', phiFile.read(16))[0].decode('utf-8')
                        self.toplabel = check + temptop
                        #print(self.toplabel)
                    else:
                        temptop = struct.unpack('20s', phiFile.read(20))[0].decode('utf-8')
                        self.toplabel = temptop
                        #print(self.toplabel)
                    #print "toplabel:", self.toplabel
                    junk = struct.unpack('8s', phiFile.read(8))
                    if is64:
                        junk = struct.unpack('8s', phiFile.read(8))
                    self.head = struct.unpack('10s', phiFile.read(10))[0].decode('utf-8')
                    # print "head:", self.head
                    self.title = struct.unpack('60s', phiFile.read(60))[0].decode('utf-8')
                    #print "title:", self.title
                    junk = struct.unpack('8s', phiFile.read(8))
                    if is64:
                        junk = struct.unpack('8s', phiFile.read(8))
                    #next line raises error if grid too big
                    #GxGxG -> packed into an array xyz order samplePhi = array.array('f')
                    tempArray.fromfile(phiFile, gridSize**3)
                    tempArray.byteswap()
                    #for count in xrange(gridSize**3):
                    #  bats = phiFile.read(4) #raw characters
                    #  blah = struct.unpack('>f', bats)[0] #always big-endian
                    #  tempArray.append(blah)
                    junk = struct.unpack('8s', phiFile.read(8))
                    if is64:
                        junk = struct.unpack('8s', phiFile.read(8))
                    self.gridDimension = gridSize
                    self.phiArray = tempArray
                    break  # read successfully, just go on and read the last bits
                except EOFError:
                    phiFile.close()
            self.botlabel = struct.unpack('16s', phiFile.read(16))[0].decode('utf-8')
            #print "botlabel:", self.botlabel
            junk = struct.unpack('8s', phiFile.read(8))
            if is64:
                junk = struct.unpack('8s', phiFile.read(8))
            #>ffff on next line forces big-endian reading
            (self.scale, self.oldmid[0], self.oldmid[1], self.oldmid[2],) = \
                struct.unpack('>ffff', phiFile.read(16))
      
            # print("scale, oldmid:", self.scale, self.oldmid)
            # print("the scale is", (1/self.scale)*gridSize)
            #box_corner_x = self.oldmid[0]-(((1/self.scale)*gridSize)/2)
            #box_corner_y = self.oldmid[1]-(((1/self.scale)*gridSize)/2)
            #box_corner_z = self.oldmid[2]-(((1/self.scale)*gridSize)/2)
            #box_corner = [box_corner_x, box_corner_y, box_corner_z]
            #print(box_corner)
            junk = struct.unpack('4s', phiFile.read(4))
            phiFile.close()


    def get_stuff(self):
        return self.scale, self.oldmid, self.gridDimension 
  
    #def get_oldmid(self)

    def get_phi_values(self):
        phi_list = []
        for phi in self.phiArray:
            phi_list.append(float(phi))
        return phi_list

def write_out_dx_file(file,xn,yn,zn,dx,dy,dz,origin,values):
    fileh = open(file,'w')
    #object 1 class gridpositions counts 40 40 40
    #origin 35.31 27.576 18.265
    #delta 0.5 0 0
    #delta 0 0.5 0
    #delta 0 0 0.5
    #object 2 class gridconnections counts 40 40 40
    #object 3 class array type float rank 0 items 64000 data follows

    fileh.write('object 1 class gridpositions counts %d %d %d\n' % (xn,yn,zn))
    #fileh.write('origin %6.2f %6.2f %6.2f\n' % (origin[0],origin[1],origin[2]))
    fileh.write('origin %7.4f %7.4f %7.4f\n' % (origin[0],origin[1],origin[2]))
    fileh.write('delta %6.3f 0 0\n' % dx)
    fileh.write('delta 0 %6.3f 0\n' % dy)
    fileh.write('delta 0 0 %6.3f\n' % dz)
    fileh.write('object 2 class gridconnections counts %d %d %d\n' % (xn,yn,zn))
    fileh.write('object 3 class array type float rank 0 items %d data follows\n' % len(values))

    count = 1
    for value in values:
        if (value == 0.0):
            fileh.write('%d' % 0)
        else:
            fileh.write('%f' % value)
            # print newline after 3rd number.
        if (count == 3):
            fileh.write('\n')
            count = 0
            # print space after number but not at the end of the line.
        else:
            fileh.write(' ')
        count = count + 1

    # if the last line has less than 3 numbers then print the a newline.
    if (count < 3):
        fileh.write('\n')
    fileh.close()



def construct_box(dx_file_name, origin, spacing, gridSize, phi_list):

    matrix = []
    for x in range(gridSize):
        y_list = []
        for y in range(gridSize):
            z_list = []
            for z in range(gridSize):
                z_list.append(0)
            y_list.append(z_list)
        matrix.append(y_list)

    index = 0
    for i in range(gridSize):
        for j in range(gridSize):
            for k in range(gridSize):
                matrix[i][j][k] = phi_list[index]
                index += 1
    #print(matrix)


    write_out_list = []
  
    for k in range(gridSize):
        for j in range(gridSize):
            for i in range(gridSize):
                write_out_list.append(matrix[i][j][k])	


    dx = spacing[0]
    dy = spacing[1]
    dz = spacing[2]

    write_out_dx_file(dx_file_name,gridSize,gridSize,gridSize,dx,dy,dz,origin,write_out_list)


def calculate_box_parameters(scale, oldmid, gridSize):

    grid_spacing = float(1/scale)
    side_len = float(grid_spacing*gridSize)

    x = float(oldmid[0]-(side_len/2))
    y = float(oldmid[1]-(side_len/2))
    z = float(oldmid[2]-(side_len/2))

    origin = [x, y, z]

    dx = grid_spacing
    dy = grid_spacing
    dz = grid_spacing
    spacing = [dx, dy, dz]

    return origin, spacing


def determine_grid_size(input_phi_file):

  file_size = os.path.getsize(input_phi_file)
  grid_bytes = file_size - 162  # 162 is number of fixed bytes in a grid file
  grid_points = grid_bytes / 4  # 4 bytes per float
  grid_size = grid_points ** (1.0/3.0)  # Cube root of grid points is size
  return int(grid_size+1)


def create_trim_electrostatics_dx(trim_electrostatics_phi, dx_file_name):
    phiSize = determine_grid_size(trim_electrostatics_phi)
    # print(phiSize)
    phiData = phi(trim_electrostatics_phi, gridSizes=(phiSize,))
    scale, oldmid, gridSize = phiData.get_stuff()
    phi_list = phiData.get_phi_values()

    origin, spacing = calculate_box_parameters(scale, oldmid, gridSize)

    construct_box(dx_file_name, origin, spacing, gridSize, phi_list)
