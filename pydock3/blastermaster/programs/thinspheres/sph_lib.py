from pydock3.util import get_hexdigest_of_persistent_md5_hash_of_tuple


class Sphere(object):
    def __init__(self, index, X, Y, Z, radius, atomnum, critical_cluster, sphere_color):
        self.index = int(index)
        self.X = float(X)
        self.Y = float(Y)
        self.Z = float(Z)
        self.radius = float(radius)
        self.atomnum = int(atomnum)
        self.critical_cluster = int(critical_cluster)
        self.sphere_color = int(sphere_color)


def cmp(a, b):
    return bool(a > b) - bool(a < b)


def by_index(x, y):
    return cmp(x.index, y.index)


def are_equal(a, b):
    return a.X == b.X and a.Y == b.Y and a.Z == b.Z


# The number of the atom with which surface point j (second point used to generate the sphere) is associated.
# The critical cluster to which this sphere belongs.
# The sphere color. The color is simply an index into the color table that was specified in the header. Therefore, 1 corresponds to the first color in the header, 2 for the second, etc. 0 corresponds to unlabeled.


def in_list(val, list):
    for ele in list:
        if are_equal(val, ele):
            return True
    return False


def remove_val(val, index, list):
    # remove val from list
    list_ele = []  ## list of elements the same as value.
    for i in range(
        index + 1, len(list)
    ):  ## We assume that only element below in the list can be equal.
        ## this is valid because we start at the begin.
        if are_equal(val, list[i]):  ## we chech does the frist element have dup,
            list_ele.append(i)  ## does the second, 3rd, and so on.
            ## so, no earlyer element will be a duplicate of the val.
    count = 0
    for i in list_ele:
        list.pop(i - count)
        count = count + 1  # index will change everytime somthing is removed.


# this function will remove duplicates from the list
# the duplicates have the same X,Y,Z coordinates.
def remove_duplicates(l):
    index = 0
    for ele in l:
        remove_val(ele, index, l)
        index = index + 1


# FORMAT: (I5, 3F10.5, F8.3, I5, I2, I3)
def read_sph(filename, chosen_cluster, color):
    sphere_list = []
    insph = open(filename, "r")
    flag_cluster = False  # this flag determines if the sphere is writen to list

    for line in insph:
        if line[0:4] == "DOCK":
            continue
        elif line[0:4] == "clus":
            cluster = int(line[7:16])
            if chosen_cluster == "A":
                flag_cluster = True
            elif int(chosen_cluster) == cluster:
                flag_cluster = True
            else:
                flag_cluster = False
        elif not (line[0:5].replace(" ", "").isdigit()):
            pass
        else:
            index = int(line[0:5])
            x = float(line[5:15])
            y = float(line[15:25])
            z = float(line[25:35])
            r = float(line[35:43])
            if r == 0.0:
                r = 0.5

            atomnum = int(line[43:48])
            if line[48:50] != "  ":
                clust = int(line[48:50])
            else:
                clust = 0
            if line[50:53] != "   ":
                col = int(line[50:53])
            else:
                col = 0

            if color == "A" or int(color) == col:
                flag_color = True
            else:
                flag_color = False

            tmp_sphere = Sphere(index, x, y, z, r, atomnum, clust, col)

            if flag_cluster and flag_color:
                # only put sphere on list if it is in a cluster of interested
                # and if the color is the same
                sphere_list.append(tmp_sphere)

    # sphere_list.append(tmp_sphere)
    # sphere_list.sort(byIndex)
    # sphere_list.sort(key = lambda a,b: cmp(a.index, b.index))
    sphere_list.sort(key=lambda a: a.index)
    # sphere_list.sort(key = byIndex)
    # remove duplicates:
    remove_duplicates(sphere_list)

    if len(sphere_list) == 0:
        raise Exception("sphere list is empty")

    return sphere_list


def write_sph(filename, spheres):
    outsph = open(filename, "w")
    outsph.write("DOCK spheres generated from read_write_sph.py\n")
    outsph.write("cluster     1   number of spheres in cluster %3d\n" % len(spheres))
    for i in range(len(spheres)):
        outsph.write(
            "%5d%10.5f%10.5f%10.5f%8.3f%5d%2d%3d\n"
            % (
                i + 1,
                round(spheres[i].X, 3),
                round(spheres[i].Y, 3),
                round(spheres[i].Z, 3),
                spheres[i].radius,
                spheres[i].atomnum,
                0,
                spheres[i].sphere_color,
            )
        )
