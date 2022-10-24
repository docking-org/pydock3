#!/bin/env python
"""Generate delphi sphere pool, near a given distance from surface, using rec.ms.

Michael Mysinger 201202 Created
Jiankun Lyu, 201701, added parms to control sphere closeness to surface and radius size.
"""

import os
import sys
import logging
import os.path as op
from optparse import OptionParser

DISTANCE = 1.8
# DISTANCE = 1.4
SIZE = 1.9
# SIZE = 1.5


class ScriptError(Exception):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)


# FORTRAN FORMAT: (I5, 3F10.5, F8.3, I5, I2, I3)
def format_sphere_line(atom_num, sphere, size=SIZE):
    """Format a line of a DOCK .sph file"""
    return "%5d%10.5f%10.5f%10.5f%8.3f%5d%2d%3d" % (
        atom_num,
        sphere[0],
        sphere[1],
        sphere[2],
        size,
        atom_num,
        0,
        0,
    )


def thin_spheres(in_f, out_f, distance=DISTANCE, size=SIZE):
    """Generate delphi sphere pool near a given distance using rec.ms."""
    print(("Distance from surface specified: %.2f" % distance))
    spheres = []
    for line in in_f:
        if line[40] == "S":
            splits = line.split()
            # point = [float(x) for x in splits[3:6]]
            # normal = [float(x) for x in splits[8:11]]
            point = [float(line[13:21]), float(line[21:30]), float(line[30:39])]
            normal = [
                float(line[43:50]),
                float(line[50:57]),
                float(line[57:64]),
                float(line[64:71]),
            ]
            # print (line)
            # print (point)
            # print (normal)
            sphere = [p + distance * n for p, n in zip(point, normal)]
            # Handle odd chain id placement
            try:
                atom_num = int(splits[1])
            except (ValueError):
                atom_num = int(splits[1][:-1])
            spheres.append((atom_num, sphere))
    out_f.write("cluster     0   number of spheres in cluster %5d\n" % len(spheres))
    for atom_num, sphere in spheres:
        out_f.write(format_sphere_line(atom_num, sphere, size=size) + "\n")


def handleio(infile=None, outfile=None, distance=DISTANCE, size=SIZE):
    """I/O handling for the script."""
    if infile is None:
        in_f = sys.stdin
    else:
        in_f = open(infile, "r")
    if outfile is None:
        out_f = sys.stdout
    else:
        out_f = open(outfile, "w")
    try:
        try:
            thin_spheres(in_f, out_f, distance=distance, size=size)
        except (ScriptError, message):
            logging.error(message)
            return message.value
    finally:
        in_f.close()
        out_f.close()
    return 0


def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    description = (
        "Generate delphi sphere pool, near a given distance from surface, using rec.ms."
    )
    usage = "%prog [options]"
    version = "%prog: version 201202 - created by Michael Mysinger"
    parser = OptionParser(usage=usage, description=description, version=version)
    parser.set_defaults(infile=None, outfile=None, distance=DISTANCE, size=SIZE)
    parser.add_option("-i", "--infile", help="input file (default: stdin)")
    parser.add_option("-o", "--outfile", help="output file (default: stdout)")
    parser.add_option(
        "-d",
        "--distance",
        type="float",
        help="sphere distance from receptor surface (default: %default)",
    )
    # added by Jiankun Lyu
    parser.add_option(
        "-s", "--size", type="float", help="sphere radius (default: %default)"
    )

    options, args = parser.parse_args(args=argv[1:])
    if len(args):
        parser.error(
            "program takes no positional arguments.\n"
            + "  Use --help for more information."
        )
    # return handleio(infile=options.infile, outfile=options.outfile,
    #                distance=options.distance)
    return handleio(
        infile=options.infile,
        outfile=options.outfile,
        distance=options.distance,
        size=options.size,
    )


if __name__ == "__main__":
    sys.exit(main(sys.argv))
