import struct
import array
import os
import math
import copy
import itertools
import logging

import numpy as np


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


# format follows
#       character*20 toplabel
#       character*10 head,character*60 title
#       real*4 phi(65,65,65) #or now 193,193,193
#       character*16 botlabel
#       real*4 scale, oldmid(3)


NUM_BYTES_PER_GRID_FILE = 162
NUM_BYTES_PER_FLOAT = 4


def grid_size_from_file_size(file_size):

    grid_bytes = file_size - NUM_BYTES_PER_GRID_FILE
    grid_points = grid_bytes / NUM_BYTES_PER_FLOAT
    grid_size = int(np.ceil(np.cbrt(grid_points)))
    return grid_size


class Phi(object):
    def __init__(self, phi_file=None, is_64=False, grid_sizes=None):
        """reads the phi file from disk"""
        if grid_sizes is None:
            grid_sizes = (None,)
        self.oldmid = [0.0, 0.0, 0.0]
        self.__minsmaxs = None
        if phi_file is not None:  # otherwise just creating an empty phi map for writing
            for grid_size in grid_sizes:
                if grid_size is None:
                    grid_size = grid_size_from_file_size(os.stat(phi_file.path).st_size)
                    logger.debug(f"Determined size to be {grid_size}")
                try:
                    f = open(phi_file.path, "rb")
                    temp_array = array.array("f")
                    _ = struct.unpack("4s", f.read(4))
                    (check,) = struct.unpack("4s", f.read(4))
                    if str(check) == "b'now '":  # this changed, but this is now correct
                        logger.debug("32-bit phimap")
                        pass
                    else:
                        logger.debug("64-bit phimap")
                        is_64 = True
                    if not is_64:
                        (temptop,) = struct.unpack("16s", f.read(16))
                        self.toplabel = check + temptop
                    else:
                        (temptop,) = struct.unpack("20s", f.read(20))
                        self.toplabel = temptop
                    logger.debug(f"toplabel: {self.toplabel}")
                    _ = struct.unpack("8s", f.read(8))
                    if is_64:
                        _ = struct.unpack("8s", f.read(8))
                    (self.head,) = struct.unpack("10s", f.read(10))
                    logger.debug(f"head: {self.head}")
                    (self.title,) = struct.unpack("60s", f.read(60))
                    logger.debug(f"title: {self.title}")
                    _ = struct.unpack("8s", f.read(8))
                    if is_64:
                        _ = struct.unpack("8s", f.read(8))
                    # next line raises error if grid too big
                    # GxGxG -> packed into an array xyz order
                    temp_array.fromfile(f, grid_size**3)
                    temp_array.byteswap()

                    self.grid_dimension = grid_size
                    self.phi_array = temp_array
                    break  # read successfully, just go on and read the last bits

                except EOFError:
                    f.close()  # TODO: should a return come after this?

            _ = struct.unpack("8s", f.read(8))
            (self.botlabel,) = struct.unpack("16s", f.read(16))
            logger.debug(f"botlabel: {self.botlabel}")
            _ = struct.unpack("8s", f.read(8))
            if is_64:
                _ = struct.unpack("8s", f.read(8))
            # >ffff on next line forces big-endian reading
            (
                self.scale,
                self.oldmid[0],
                self.oldmid[1],
                self.oldmid[2],
            ) = struct.unpack(">ffff", f.read(16))
            logger.debug(f"\tscale: {self.scale}\n\toldmid: {self.oldmid}")
            _ = struct.unpack("4s", f.read(4))
            f.close()

    def write(self, phi_file=None):
        """write data to member data structure manually,
        then call this to write to file
        the pad lines reproduce the binary padding of an original
        fortran formatted phi file"""
        if phi_file is not None:  # do nothing if no filename given
            out_array = copy.deepcopy(self.phi_array)
            out_array.byteswap()  # switch endianness back, only for writing
            with open(
                phi_file.path, "wb"
            ) as f:  # TODO: b may be unnecessary, have to check
                f.write(struct.pack("4b", 0, 0, 0, 20))  # pad
                f.write(struct.pack("20s", self.toplabel))
                f.write(struct.pack("8b", 0, 0, 0, 20, 0, 0, 0, 70))  # pad
                f.write(struct.pack("10s", self.head))
                f.write(struct.pack("60s", self.title))
                f.write(struct.pack("4b", 0, 0, 0, 70))  # pad, always same
                f.write(struct.pack(">l", len(out_array) * 4))  # diff. pad sometimes
                out_array.tofile(f)  # array
                f.write(struct.pack(">l", len(out_array) * 4))  # diff. pad sometimes
                f.write(struct.pack("4b", 0, 0, 0, 16))  # pad, always same
                f.write(struct.pack("16s", self.botlabel))
                f.write(struct.pack("8b", 0, 0, 0, 16, 0, 0, 0, 16))  # pad
                f.write(
                    struct.pack(
                        ">ffff",
                        self.scale,
                        self.oldmid[0],
                        self.oldmid[1],
                        self.oldmid[2],
                    )
                )
                # > on previous line forces big-endian writing
                f.write(struct.pack("4b", 0, 0, 0, 16))  # pad

    def trim_phi(self, new_mid_indices, new_size):
        """for a new center index and a desired cubic grid size, trim the current
        phimap and return the new trimmed phimap"""
        plus_minus = (new_size - 1) / 2  # how many to add or subtract from the center
        new_phi = Phi()
        new_phi.oldmid = self.get_xyz_list(new_mid_indices)  # only change of these data
        new_phi.toplabel = self.toplabel
        new_phi.head = self.head
        new_phi.title = self.title
        new_phi.botlabel = self.botlabel
        new_phi.scale = self.scale
        # the phi_array does change
        new_phi.phi_array = array.array("f")
        new_phi.phi_array.fromlist(
            [
                self.get_value(old_index_x, old_index_y, old_index_z)
                if np.all(
                    [
                        0 <= old_index < self.grid_dimension
                        for old_index in [old_index_x, old_index_y, old_index_z]
                    ]
                )
                else 0.0
                for old_index_z, old_index_y, old_index_x in itertools.product(
                    *[
                        range(
                            int(new_mid_index - plus_minus),
                            int(new_mid_index + plus_minus + 1),
                        )
                        for new_mid_index in reversed(new_mid_indices)
                    ]
                )
            ]
        )

        return new_phi

    def get_mins_maxs(self):
        """finds the positions of the extreme grid corners"""
        if self.__minsmaxs is None:
            mins = [
                center - ((self.grid_dimension - 1.0) / (2.0 * self.scale))
                for center in self.oldmid
            ]
            maxs = [
                center + ((self.grid_dimension - 1.0) / (2.0 * self.scale))
                for center in self.oldmid
            ]
            self.__minsmaxs = mins, maxs
        return self.__minsmaxs

    def get_xyz_list(self, xyz):
        """changes list to x,y,z calls get_xyz"""
        return self.get_xyz(*xyz)

    def get_xyz(self, x_ind, y_ind, z_ind):
        """returns the xyz coordinate of the center of the box"""
        mins, maxs = self.get_mins_maxs()
        gap = 1.0 / self.scale
        return mins[0] + (x_ind * gap), mins[1] + (y_ind * gap), mins[2] + (z_ind * gap)

    def get_value(self, x_ind, y_ind, z_ind):
        """for a given set of indices, return the value in the array"""
        index = int(
            (z_ind * (self.grid_dimension**2.0))
            + (y_ind * self.grid_dimension)
            + x_ind
        )
        return self.phi_array[index]

    def subtract(self, other):
        """subtract other from self, destructively write over self"""
        self.modify(other, -1)

    def add(self, other):
        """add other to self, destructively write over self."""
        self.modify(other, 1)

    def modify(self, other, change):
        """modify other to self, destructively write over self. allows +-/etc
        presume without checking that grids are compatible (same mid etc)"""
        self.phi_array = [
            this_phi_value + (other_phi_value * change)
            for this_phi_value, other_phi_value in zip(self.phi_array, other.phi_array)
        ]

    def get_indices(self, pt):
        """helper function to find the box a point is in"""
        mins, maxs = self.get_mins_maxs()
        grid_size = 1.0 / self.scale
        return tuple(
            int(math.floor((pt_i - minimum) / grid_size))
            for pt_i, minimum in zip(pt, mins)
        )

    def trim_to_box_center_and_size(self, corners, center):
        """given a box, find the new center and size of a valid phimap based on
        this current phimap"""
        # find the midpoint and corners
        center_indices = self.get_indices(center)
        onecorner = self.get_indices(corners[:3])
        twocorner = [coord + 1 for coord in self.get_indices(corners[3:6])]
        # phimap grid can only be cubic
        biggest_dimension = 0
        if twocorner[1] - onecorner[1] > twocorner[0] - onecorner[0]:
            biggest_dimension = 1
        if (
            twocorner[2] - onecorner[2]
            > twocorner[biggest_dimension] - onecorner[biggest_dimension]
        ):
            biggest_dimension = 2
        new_size = twocorner[biggest_dimension] - onecorner[biggest_dimension]
        if 0 == (new_size % 2):  # if size is even, that's not allowed, so,
            new_size += 1  # make it odd
        return center_indices, new_size

    def trim_to_box(self, corners, center):
        """given a box (see box.py) trim so that the box is enclosed but not more.
        returns the new trimmed phimap"""
        center_indices, new_size = self.trim_to_box_center_and_size(corners, center)
        return self.trim_phi(center_indices, new_size), center_indices, new_size


def add(input_grid, add_this_grid, output_grid, phi_size=None):
    """does the addition, all read code in phi.py"""
    phi_data = Phi(input_grid, grid_sizes=(phi_size,))
    phi_data_2 = Phi(add_this_grid, grid_sizes=(phi_size,))
    phi_data.add(phi_data_2)
    phi_data.write(output_grid)


def subtract(input_grid, subtract_this_grid, output_grid, phi_size=None):
    """does the subtraction, all read code in phi.py"""
    phi_data = Phi(input_grid, grid_sizes=(phi_size,))
    phi_data_2 = Phi(subtract_this_grid, grid_sizes=(phi_size,))
    phi_data.subtract(phi_data_2)
    phi_data.write(output_grid)


def read_box_file(box_file):
    """reads a 'box' file from docking, used for constructing other grids.
    since things outside this box can't be scored, we don't need to save that data
    for reading into DOCK. ridiculously allows any odd-numbered gridsize."""
    corners = []
    center = []
    dimensions = []
    with open(box_file.path, "r") as f:
        for line in f.readlines():
            if line.find("CORNERS") > 0:
                corners = [float(item) for item in line.split()[4:]]
            elif line.find("CENTER") > 0:
                center = [float(item) for item in line.split()[5:]]
            elif line.find("DIMENSIONS") > 0:
                dimensions = [float(item) for item in line.split()[5:]]
    return corners, center, dimensions


def trim(input_phi_file, box_file, output_phi_file):
    """reads a phimap, trims it to be just outside the box, writes the new map
    out. returns the new cubic size."""
    phi_data = Phi(input_phi_file)
    corners, center, dimensions = read_box_file(box_file)
    trimmed_phi_data, new_center, new_size = phi_data.trim_to_box(corners, center)
    trimmed_phi_data.write(output_phi_file)
    return new_size, new_center
