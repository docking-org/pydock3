import struct
import sys
import array

# written by Reed Stein
# Converted to python3 by Brendan Hall


def write_out_dx_file(file,xn,yn,zn,dx,dy,dz,origin,values):
    # print("I AM HERE in write_out_dx_file")
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
    fileh.write('delta %6.6f 0 0\n' % dx)
    fileh.write('delta 0 %6.6f 0\n' % dy)
    fileh.write('delta 0 0 %6.6f\n' % dz)
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


#def read_bump(bump_file):
#
#        bump_open = open(bump_file,'r')
#        bump_read = bump_open.readlines()
#        bump_open.close()
#
#        x_coords = []
#        y_coords = []
#        z_coords = []
#
#        for line in bump_read[0:2]:
#                #print(line)
#                line = line.strip().split()
#                spacing = '0.200'
#                if line[0] == spacing:
#                        print(line)
#                        box_corner_x = float(line[1])
#                        box_corner_y = float(line[2])
#                        box_corner_z = float(line[3])
#                        x_dim = int(line[4])
#                        y_dim = int(line[5])
#                        z_dim = int(line[6])
#
#	return(x_dim, y_dim, z_dim, [box_corner_x, box_corner_y, box_corner_z])

def construct_box(dx_file_name, origin, spacing, xn, yn, zn, values):
        matrix = []
        for x in range(xn):
            y_list = []
            for y in range(yn):
                z_list = []
                for z in range(zn):
                    z_list.append(0)
                y_list.append(z_list)
            matrix.append(y_list)

        #index = 0
        #for k in range(zn):
        #        for j in range(yn):
        #                for i in range(xn):
        #                        matrix[i][j][k] = values[index]
        #                        index += 1
        index = 0
        for i in range(xn):
            for j in range(yn):
                for k in range(zn):
                    matrix[i][j][k] = values[index]
                    index += 1

        #write_out_list = []
        #for i in range(xn):
        #	for j in range(yn):
        #		for k in range(zn):
        #			write_out_list.append(matrix[i][j][k])

        write_out_list = []
        for i in range(xn-1,-1,-1):
            for j in range(yn-1,-1,-1):
                for k in range(zn-1,-1,-1):
                    write_out_list.append(matrix[i][j][k])


        #write_out_list = []
        #for k in range(zn):
        #        for j in range(yn):
        #                for i in range(xn):
        #			write_out_list.append(matrix[i][j][k])

        dx = spacing
        dy = spacing
        dz = spacing
        write_out_dx_file(dx_file_name,xn,yn,zn,dx,dy,dz,origin,write_out_list)
        #write_out_dx_file(dx_file_name,zn,yn,xn,dx,dy,dz,[origin[2],origin[1],origin[0]],write_out_list)
        #write_out_dx_file(dx_file_name,xn,yn,zn,dx,dy,dz,origin,values)

def read_bump(bump_file):

        bump_open = open(bump_file,'r')
        bump_read = bump_open.readlines()
        bump_open.close()

        x_coords = []
        y_coords = []
        z_coords = []


        line = bump_read[0]
        splitline = line.split()

        x_dim = int(splitline[0])
        y_dim = int(splitline[1])
        z_dim = int(splitline[2])
        spaceperang        = float(splitline[3])
        space = 1.0/spaceperang
        box_upper_corner_x_space = float(splitline[4]) # these values are not in angstroms but grid spaceing. 
        box_upper_corner_y_space = float(splitline[5])
        box_upper_corner_z_space = float(splitline[6])

        upper_corner_x = (box_upper_corner_x_space ) * space
        upper_corner_y = (box_upper_corner_y_space ) * space
        upper_corner_z = (box_upper_corner_z_space ) * space

        botom_corner_x = (box_upper_corner_x_space - x_dim) * space
        botom_corner_y = (box_upper_corner_y_space - y_dim) * space
        botom_corner_z = (box_upper_corner_z_space - z_dim) * space

        # print(x_dim, y_dim, z_dim, space, box_upper_corner_x_space, box_upper_corner_y_space, box_upper_corner_z_space)
        # print(botom_corner_x, botom_corner_y, botom_corner_z)
        # print(upper_corner_x, upper_corner_y, upper_corner_z)
        #exit()

        values = []
        for line in bump_read[1:]:
            splitline = line.split()
            for val in splitline:
                values.append(float(val))

        # print(len(values), (x_dim + 1) * (y_dim+1)* (z_dim+1))

        #values = []	
        #count = 0 
        #for z in z_coords:
        #        for y in y_coords:
        #                for x in x_coords:
        #                        values.append(bump_string[count])
        #			#if bump_string[count] == "F":
        #			#	values.append(float(1))
        #			#else:
        #			#	values.append(float(0))
        #			count += 1
#
        return(values, [x_dim+1, y_dim+1, z_dim+1], [botom_corner_x, botom_corner_y, botom_corner_z ], space)
        

def create_ligdesolv_dx(lig_desolv_heavy, lig_desolv_dx):
        values, dims, origin, spacing = read_bump(lig_desolv_heavy)
        construct_box(lig_desolv_dx, origin, spacing, dims[0], dims[1], dims[2], values)

