import struct
import sys
import array

# written by Reed Stien


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


# def read_bump(bump_file):

#     bump_open = open(bump_file,'r')
#     bump_read = bump_open.readlines()
#     bump_open.close()

#     x_coords = []
#     y_coords = []
#     z_coords = []

#     for line in bump_read[0:2]:
#             #print(line)
#             line = line.strip().split()
#             spacing = '0.200'
#             if line[0] == spacing:
#                     print(line)
#                     box_corner_x = float(line[1])
#                     box_corner_y = float(line[2])
#                     box_corner_z = float(line[3])
#                     x_dim = int(line[4])
#                     y_dim = int(line[5])
#                     z_dim = int(line[6])

# 	return(x_dim, y_dim, z_dim, [box_corner_x, box_corner_y, box_corner_z])

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

    index = 0
    for k in range(zn):
        for j in range(yn):
            for i in range(xn):
                matrix[i][j][k] = values[index]
                index += 1

    write_out_list = []
    for i in range(xn):
        for j in range(yn):
            for k in range(zn):
                write_out_list.append(matrix[i][j][k])

    dx = spacing
    dy = spacing
    dz = spacing
    write_out_dx_file(dx_file_name,xn,yn,zn,dx,dy,dz,origin,write_out_list)

def read_bump(bump_file):

    bump_open = open(bump_file,'r')
    bump_read = bump_open.readlines()
    bump_open.close()

    x_coords = []
    y_coords = []
    z_coords = []

    for line in bump_read[0:2]:
        line = line.strip().split()
        spacing = '0.200' 
        if line[0] == spacing:
            box_corner_x = float(line[1])
            box_corner_y = float(line[2])
            box_corner_z = float(line[3])
            x_dim = int(line[4])
            y_dim = int(line[5])
            z_dim = int(line[6])

    x_coords = [box_corner_x]
    y_coords = [box_corner_y]
    z_coords = [box_corner_z]

    bump_string = ''	

    for line in bump_read[2:]:
        line = line.strip()
        bump_string+=line

    spacing = float(0.200) 
    for i in range(x_dim-1):
        new_val = box_corner_x+spacing
        x_coords.append(new_val)
        box_corner_x = new_val

    for i in range(y_dim-1):
        new_val = box_corner_y+spacing
        y_coords.append(new_val)
        box_corner_y = new_val

    for i in range(z_dim-1):
        new_val = box_corner_z+spacing
        z_coords.append(new_val)
        box_corner_z = new_val


    values = []	
    count = 0 
    for z in z_coords:
        for y in y_coords:
            for x in x_coords:
                if bump_string[count] == "F":
                    values.append(float(1))
                else:
                    values.append(float(0))
                count += 1

    return(values, [x_dim, y_dim, z_dim], [x_coords[0], y_coords[0], z_coords[0]])
    

def create_vdw_dx(vdw_file, bmp_file, vdw_repuslive_dx, vdw_attractive_dx, vdw_dx):
    values, dims, origin = read_bump(bmp_file)
    x_dim, y_dim, z_dim = dims

    spacing = 0.200
    vdwFile = open(vdw_file, 'rb')  # b is for binary, r is for read
    tempArray = array.array('f')
    tempArray.fromfile(vdwFile, x_dim * y_dim * z_dim)
    tempArray.byteswap()

    construct_box(vdw_repuslive_dx, origin, spacing, x_dim, y_dim, z_dim, tempArray)
    tempArray = array.array('f')
    tempArray.fromfile(vdwFile, x_dim * y_dim * z_dim)
    tempArray.byteswap()
    construct_box(vdw_attractive_dx, origin, spacing, x_dim, y_dim, z_dim, tempArray)

    construct_box(vdw_dx, origin, spacing, x_dim, y_dim, z_dim, values)
