"""Python port of makebox.smallokay.pl (Austin N. Kirschner 2002-2003; rgc 2012).

Makes the grid box around the ligand spheres, padded by `margin` angstroms on each side.
If the box holds too many grid points, the sphere farthest from the box center is moved
to the center until it fits; if it holds too few, each side is grown (or shrunk) by one
angstrom per cycle depending on whether receptor atoms are near that side.

Ported line by line from the Perl so that the box file is identical.
"""
import math

from pydock3.blastermaster.programs.perl_semantics import field, num, read_tokens

POINTS_PER_ANGSTROM = 3
MAX_BOX_POINTS = 1700000
MIN_BOX_POINTS = 50000
MAX_SIDE = 50.0  # angstroms
INCREMENT = 1  # angstroms per side per growing cycle
MAX_GROWING_CYCLES = 200


def _dist(p, q):
    return math.sqrt((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2)


def _num_points(d1, d2):
    return (
        abs(d2[0] - d1[0]) * abs(d2[1] - d1[1]) * abs(d2[2] - d1[2])
        * POINTS_PER_ANGSTROM * POINTS_PER_ANGSTROM * POINTS_PER_ANGSTROM
    )


def _read_sphere_coords(spheres_path):
    """Coordinates on every line after the first "cluster" line."""
    coords = []
    started = False
    for tokens in read_tokens(spheres_path):
        if started:
            coords.append([num(field(tokens, k)) for k in (1, 2, 3)])
        if field(tokens, 0) == "cluster":
            started = True
    return coords


def _read_receptor_coords(receptor_path):
    coords = []
    for tokens in read_tokens(receptor_path):
        if field(tokens, 0) in ("ATOM", "HETATM"):
            # coordinates are tokens 5-7, or 6-8 if there is a chain ID
            columns = (5, 6, 7) if "." in field(tokens, 5) else (6, 7, 8)
            coords.append([num(field(tokens, k)) for k in columns])
    return coords


def _box_around(coords, margin):
    """Diagonal corners of the box enclosing all `coords` plus `margin`."""
    d1 = [coords[0][k] - margin for k in range(3)]
    d2 = [coords[0][k] + margin for k in range(3)]
    for p in coords[1:]:
        for k in range(3):
            if d1[k] + margin <= p[k] <= d2[k] - margin:
                pass
            elif p[k] < d1[k] + margin:
                d1[k] = p[k] - margin
            elif d2[k] - margin < p[k]:
                d2[k] = p[k] + margin
    return d1, d2


def _grow_box(d1, d2, receptor_coords, margin):
    """Move each side of a too-small box out if receptor atoms are near it, in if none are close."""
    num_points = _num_points(d1, d2)
    cycle = 1
    while num_points < MIN_BOX_POINTS and cycle <= MAX_GROWING_CYCLES:
        cycle += 1
        # center of each wall, in the order: -X, +X, -Y, +Y, -Z, +Z
        walls = []
        for k in range(3):
            for corner in (d1, d2):
                center = [(d1[m] + d2[m]) / 2 for m in range(3)]
                center[k] = corner[k]
                walls.append(center)

        near = [False] * 6  # a receptor atom is within margin/2 of the wall center
        i = 0
        while i < len(receptor_coords) and not all(near):
            for w, wall in enumerate(walls):
                if not near[w] and _dist(receptor_coords[i], wall) < margin / 2:
                    near[w] = True
            i += 1

        closest = [margin] * 6  # distance of the closest receptor atom, if closer than margin
        if not all(near):
            for atom in receptor_coords:
                for w, wall in enumerate(walls):
                    d = _dist(atom, wall)
                    if d < closest[w]:
                        closest[w] = d

        for w in range(6):
            k, sign = w // 2, (-1 if w % 2 == 0 else 1)
            corner = d1 if w % 2 == 0 else d2
            if near[w]:
                corner[k] += sign * INCREMENT
        for w in range(6):
            k, sign = w // 2, (-1 if w % 2 == 0 else 1)
            corner = d1 if w % 2 == 0 else d2
            if not near[w] and closest[w] == margin:
                corner[k] -= sign * INCREMENT

        for k in range(3):
            while abs(d2[k] - d1[k]) > MAX_SIDE:
                d2[k] = d2[k] - 0.00025
                d1[k] = d1[k] + 0.00025

        num_points = _num_points(d1, d2)


def make_box(spheres_path, receptor_path, box_path, margin):
    margin = float(margin)
    coords = _read_sphere_coords(spheres_path)

    num_points = MAX_BOX_POINTS + 1
    while num_points > MAX_BOX_POINTS:
        d1, d2 = _box_around(coords, margin)
        center = [(d2[k] + d1[k]) / 2 for k in range(3)]
        num_points = _num_points(d1, d2)
        if num_points > MAX_BOX_POINTS:  # replace the point farthest from the center by the center
            distances = [_dist(center, p) for p in coords]
            farthest = 0
            for i in range(len(coords)):
                if distances[i] > distances[farthest]:
                    farthest = i
            coords[farthest] = list(center)

    if num_points < MIN_BOX_POINTS:
        _grow_box(d1, d2, _read_receptor_coords(receptor_path), margin)

    # corners are rounded to 3 decimals before the center and dimensions are computed
    d1 = [float("%8.3f" % v) for v in d1]
    d2 = [float("%8.3f" % v) for v in d2]
    center = [(d2[k] + d1[k]) / 2 for k in range(3)]
    dimensions = [abs(d2[k] - d1[k]) for k in range(3)]

    def xyz(*values):
        return "".join("%8.3f" % v for v in values)

    (x1, y1, z1), (x2, y2, z2) = d1, d2
    lines = [
        "HEADER    CORNERS OF BOX " + xyz(x1, y1, z1, x2, y2, z2),
        "REMARK    CENTER (X Y Z) " + xyz(*center),
        "REMARK    DIMENSIONS (X Y Z) " + xyz(*dimensions),
        "ATOM      1  DUA BOX     1    " + xyz(x1, y1, z1),
        "ATOM      2  DUB BOX     1    " + xyz(x2, y1, z1),
        "ATOM      3  DUC BOX     1    " + xyz(x2, y1, z2),
        "ATOM      4  DUD BOX     1    " + xyz(x1, y1, z2),
        "ATOM      5  DUE BOX     1    " + xyz(x1, y2, z1),
        "ATOM      6  DUF BOX     1    " + xyz(x2, y2, z1),
        "ATOM      7  DUG BOX     1    " + xyz(x2, y2, z2),
        "ATOM      8  DUH BOX     1    " + xyz(x1, y2, z2),
        "CONECT    1    2    4    5",
        "CONECT    2    1    3    6",
        "CONECT    3    2    4    7",
        "CONECT    4    1    3    8",
        "CONECT    5    1    6    8",
        "CONECT    6    2    5    7",
        "CONECT    7    3    6    8",
        "CONECT    8    4    5    7",
    ]
    with open(box_path, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
