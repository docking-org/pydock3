"""Python port of makespheres1.cli.pl and makespheres3.cli.pl.

Written in Perl by Austin N. Kirschner (2003), modified by Trent Balius (2014) to use
reduce atom names. The two scripts share one algorithm and differ only in their settings:
makespheres3 makes the matching spheres (it keeps the crystallographic ligand atoms as
spheres) and makespheres1 makes the low-dielectric spheres for electrostatics (it doesn't).

Algorithm, on the spheres of sphgen's cluster 0:
  1. drop spheres farther than `margin` from the ligand center, farther than 7 A from the
     receptor, or closer than 1.2 A to it
  2. flag spheres near polar receptor atoms with the hydrogen-bond angle they make, or as
     nonpolar if only near carbons
  3. thin them on a grid of `gridsize` (preferring crystallographic, then polar spheres),
     then remove spheres closer than `tooclose` to another one
  4. keep spheres continuous with the crystallographic spheres (or the sphere closest to
     the ligand center), relaxing continuity until there are at least `min_spheres`
  5. if more than `max_spheres` remain, drop those farthest from the ligand/spheres

Ported line by line from the Perl so that the output is identical, including its quirks
(noted in comments). The one deliberate difference: polar spheres are visited in sorted
order, as makespheres1 did, rather than makespheres3's random Perl hash order.
"""
import math

from pydock3.blastermaster.programs.perl_semantics import field, num, read_tokens

RECEPTOR_MAX_DIST = 7.0
RECEPTOR_MIN_DIST = 1.2
POLAR_DIST = 3.3  # to N, O, S receptor atoms
H_POLAR_DIST = 2.5  # to polar H receptor atoms
NONPOLAR_DIST = 4.5  # to C receptor atoms
CONTINUITY = 3.0
CONTINUITY_INCREMENT = 0.5
MAX_CONTINUITY = 4.5
NUM_CONTINUITY_REPEATS = 10
POLAR_BENEFIT = -0.25  # fraction of the average weighted distance given to polar spheres
PI = 3.14159265359

# Sphere flags (index 6)
PLAIN, CRYSTAL, CONTINUOUS, CRYSTAL_AND_CENTER = 0, 1, 2, 3
# Polarity (index 7): an H-bond angle, or
NEITHER, NONPOLAR = -1, -2

# Hydrogen-bond geometry for polar receptor atoms, checked in this order:
# (atom names, residue names or None for any, base atom, mean angle, min angle, max angle)
POLAR_ATOM_RULES = [
    (("H",), None, "N", 155.21, 129.2792, 181.1408),
    (("O",), None, "C", 140.70, 104.4596, 176.9404),
    (("HE",), ("ARG",), "NE", 150.13, 121.4356, 178.8244),
    (("HH11", "HH12"), ("ARG",), "NH1", 147.29, 114.8520, 179.7280),
    (("HH21", "HH22"), ("ARG",), "NH2", 146.21, 114.4384, 177.9816),
    (("HD21", "HD22"), ("ASN",), "ND2", 151.25, 114.5784, 187.9216),
    (("OD1",), ("ASN",), "CG", 131.66, 97.3796, 165.9404),
    (("OD1",), ("ASP",), "CG", 124.66, 90.0464, 159.2736),
    (("OD2",), ("ASP",), "CG", 121.82, 92.9884, 150.6516),
    (("HG",), ("CYS",), "SG", 154.82, 128.9284, 180.7116),
    (("SG",), ("CYS",), "CB", 112.75, 87.1720, 138.3280),
    (("SG",), ("CYX",), "CB", 147.05, 127.5480, 166.5520),
    (("OE1",), ("GLU",), "CD", 123.35, 86.4824, 160.2176),
    (("OE2",), ("GLU",), "CD", 123.83, 90.4512, 157.2088),
    (("HE21", "HE22"), ("GLN",), "NE2", 147.62, 108.714, 186.5260),
    (("OE1",), ("GLN",), "CD", 129.76, 96.1068, 163.4132),
    (("HD1",), ("HIS", "HID", "HIP"), "ND1", 149.52, 116.6508, 182.3892),
    (("HE2",), ("HIE", "HIP"), "NE2", 152.56, 117.6132, 187.5068),
    (("ND1",), ("HIS",), "CG", 130.34, 110.8968, 149.7832),
    (("NE2",), ("HIS",), "CD2", 131.83, 111.6224, 152.0376),
    (("HZ1", "HZ2", "HZ3"), ("LYS",), "NZ", 139.81, 102.1388, 177.4812),
    (("SD",), ("MET",), "CG", 139.81, 92.3324, 169.0076),
    (("HG",), ("SER",), "OG", 159.76, 130.2032, 189.3168),
    (("OG",), ("SER",), "CB", 124.55, 93.2292, 155.8708),
    (("HG",), ("THR",), "OG1", 163.00, 135.8148, 190.1852),
    (("OG1",), ("THR",), "CB", 125.41, 97.0488, 153.7712),
    (("HE1",), ("TRP",), "NE1", 153.50, 124.4332, 182.5668),
    (("HH",), ("TYR",), "OH", 147.32, 112.3928, 182.2472),
    (("OH",), ("TYR",), "CZ", 119.55, 94.8540, 144.2460),
]

SPH_HEADER = """DOCK 5.2 ligand_atoms
positive                       (1)
negative                       (2)
acceptor                       (3)
donor                          (4)
ester_o                        (5)
amide_o                        (6)
neutral                        (7)
not_neutral                    (8)
positive_or_donor              (9)
negative_or_acceptor           (10)
neutral_or_acceptor_or_donor   (11)
donacc                         (12)
"""


def make_matching_spheres(gridsize, tooclose, max_spheres, ligand_spheres_path, all_spheres_path, receptor_path, out_path):
    """makespheres3.cli.pl: spheres for DOCK's matching, including the crystallographic ligand atoms."""
    _make_spheres(
        ligand_spheres_path, all_spheres_path, receptor_path, out_path,
        margin=10.0, gridsize=float(gridsize), tooclose=float(tooclose),
        min_spheres=20, max_spheres=float(max_spheres), use_ligand=True,
    )


def make_low_dielectric_spheres(ligand_spheres_path, all_spheres_path, receptor_path, out_path, min_spheres):
    """makespheres1.cli.pl: spheres filling the binding site for the low-dielectric electrostatics region."""
    _make_spheres(
        ligand_spheres_path, all_spheres_path, receptor_path, out_path,
        margin=12.0, gridsize=1.5, tooclose=0.8,
        min_spheres=float(min_spheres), max_spheres=120, use_ligand=False,
    )


def _xyz(sphere):
    return sphere[1], sphere[2], sphere[3]


def _dist(p, q):
    return math.sqrt((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2)


def _mean(points):
    return [sum(p[k] for p in points) / len(points) for k in range(3)]


def _read_ligand_coords(path):
    """Coordinates of every line after the first (the cluster line)."""
    return [[num(field(tokens, k)) for k in (1, 2, 3)] for tokens in read_tokens(path)[1:]]


def _read_cluster_0(path):
    """Spheres of sphgen's cluster 0, which runs to the end of the file.

    Each sphere is a list: [number, x, y, z, radius, atom number, flag, polarity, weighted
    distance, distance of the H-bond angle from its mean].
    """
    items = read_tokens(path)
    spheres = []
    line = 0
    while line < len(items):
        if field(items[line], 1) == "0":
            line += 1
            while line < len(items):
                tokens = items[line]
                spheres.append([field(tokens, 0)] + [num(field(tokens, k)) for k in range(1, 6)] + [0, 0, 0, 0])
                line += 1
        line += 1
    return spheres


def _read_receptor_atoms(path):
    """Receptor atoms as [x, y, z, atom name, residue name, residue number]."""
    atoms = []
    for tokens in read_tokens(path):
        if field(tokens, 0) in ("ATOM", "HETATM"):
            if "." in field(tokens, 5):  # no chain ID
                x, y, z, name, residue, number = (field(tokens, k) for k in (5, 6, 7, 2, 3, 4))
            else:
                x, y, z, name, residue, number = (field(tokens, k) for k in (6, 7, 8, 2, 3, 5))
            atoms.append([num(x), num(y), num(z), name, residue, number])
    return atoms


def _make_spheres(ligand_spheres_path, all_spheres_path, receptor_path, out_path, margin, gridsize, tooclose, min_spheres, max_spheres, use_ligand):
    ligand_coords = _read_ligand_coords(ligand_spheres_path)
    if not ligand_coords:
        raise ValueError(f"There are 0 spheres in {ligand_spheres_path}")
    ligand_center = _mean(ligand_coords)

    spheres = []
    if use_ligand:
        spheres = [["crystal", x, y, z, 0, 0, CRYSTAL, 0, 0, 0] for x, y, z in ligand_coords]
    num_crystal = len(spheres)
    if len(spheres) == 1:  # quirk of the original: a lone crystallographic sphere is dropped
        spheres.pop()
        num_crystal = 0
    spheres += _read_cluster_0(all_spheres_path)
    for sphere in spheres[num_crystal:]:
        sphere[6] = PLAIN

    # 1. drop spheres too far from the ligand center (crystallographic ones are kept)
    line = 0
    while line < len(spheres):
        too_far = any(
            spheres[line][k + 1] < ligand_center[k] - margin or spheres[line][k + 1] > ligand_center[k] + margin
            for k in range(3)
        )
        if too_far and spheres[line][6] != CRYSTAL:
            del spheres[line]
            line -= 1
        line += 1

    # ... and too far from or too close to the receptor (crystallographic ones included)
    receptor_atoms = _read_receptor_atoms(receptor_path)
    receptor_atoms_by_residue_number = {}
    for atom in receptor_atoms:
        receptor_atoms_by_residue_number.setdefault(atom[5], []).append(atom)
    i = 0
    while i < len(spheres):
        dist = min(_dist(atom, _xyz(spheres[i])) for atom in receptor_atoms)
        if dist > RECEPTOR_MAX_DIST or dist < RECEPTOR_MIN_DIST:
            del spheres[i]
            i -= 1
        i += 1

    # 2. flag spheres near polar / nonpolar receptor atoms
    polar_atoms = {}  # sphere index -> nearby polar receptor atoms
    nonpolar_atoms = {}
    for sphere in spheres:
        sphere[7] = NEITHER
        sphere[9] = 1000
    for i, sphere in enumerate(spheres):
        for atom in receptor_atoms:
            dist = _dist(atom, _xyz(sphere))
            element = atom[3][:1]
            if dist <= POLAR_DIST and element in ("N", "O", "S"):
                polar_atoms.setdefault(i, []).append(atom)
            if dist <= H_POLAR_DIST and element == "H":
                polar_atoms.setdefault(i, []).append(atom)
            if dist <= NONPOLAR_DIST and element == "C":
                nonpolar_atoms.setdefault(i, []).append(atom)
    for i in nonpolar_atoms:
        if i not in polar_atoms:
            spheres[i][7] = NONPOLAR

    # The base atom is looked up by residue number; if it is missing, the original reuses
    # the base atom of the previous lookup (hence this state, and the fixed visiting order).
    base_atom = [0.0, 0.0, 0.0]

    def hydrogen_bond_angle(sphere, polar_atom, base_name):
        nonlocal base_atom
        for atom in receptor_atoms_by_residue_number.get(polar_atom[5], []):
            if atom[3] == base_name:
                base_atom = atom[:3]
        base_to_polar = [base_atom[k] - polar_atom[k] for k in range(3)]
        polar_to_sphere = [sphere[k + 1] - polar_atom[k] for k in range(3)]
        length1 = math.sqrt(base_to_polar[0] ** 2 + base_to_polar[1] ** 2 + base_to_polar[2] ** 2)
        length2 = math.sqrt(polar_to_sphere[0] ** 2 + polar_to_sphere[1] ** 2 + polar_to_sphere[2] ** 2)
        unit1 = [base_to_polar[k] / length1 for k in range(3)]
        unit2 = [polar_to_sphere[k] / length2 for k in range(3)]
        cos_angle = unit1[0] * unit2[0] + unit1[1] * unit2[1] + unit1[2] * unit2[2]
        acos = math.atan2(math.sqrt(1 - cos_angle * cos_angle), cos_angle)
        return (acos / PI) * 180

    for i in sorted(polar_atoms, key=str):  # Perl `sort keys` sorts the indices as strings
        sphere = spheres[i]
        for polar_atom in polar_atoms[i]:
            name, residue = polar_atom[3], polar_atom[4]
            for names, residues, base_name, mean, low, high in POLAR_ATOM_RULES:
                if name in names and (residues is None or residue in residues):
                    angle = hydrogen_bond_angle(sphere, polar_atom, base_name)
                    if abs(angle - mean) < sphere[9] and low <= angle <= high:
                        sphere[7] = angle
                        sphere[9] = abs(angle - mean)
                    break

    # 3. keep one sphere per grid box: crystallographic ones, else the polar one with the
    # best angle, else the one closest to the box center
    mins = list(_xyz(spheres[0]))
    maxs = list(_xyz(spheres[0]))
    for sphere in spheres[1:]:
        for k in range(3):
            if sphere[k + 1] < mins[k]:
                mins[k] = sphere[k + 1]
            if sphere[k + 1] > maxs[k]:
                maxs[k] = sphere[k + 1]
    half = gridsize / 2
    x = mins[0]
    while x <= maxs[0]:
        y = mins[1]
        while y <= maxs[1]:
            z = mins[2]
            while z <= maxs[2]:
                box = []  # [sphere index, distance to box center, flag, polarity, angle distance]
                for i, s in enumerate(spheres):
                    if x - half <= s[1] < x + half and y - half <= s[2] < y + half and z - half <= s[3] < z + half:
                        box.append([i, _dist(_xyz(s), (x, y, z)), s[6], s[7], s[9]])
                if len(box) > 1:
                    box = sorted(box, key=lambda b: b[1])
                    num_crystal_in_box = sum(1 for b in box if b[2] == CRYSTAL)
                    num_polar_in_box = sum(1 for b in box if b[3] != NEITHER and b[3] != NONPOLAR)
                    if num_crystal_in_box > 0:
                        to_delete = [b for b in box if b[2] != CRYSTAL]
                    elif num_polar_in_box > 0:
                        to_delete = sorted(box, key=lambda b: b[4])[1:]
                    else:
                        to_delete = box[1:]
                    for n, b in enumerate(to_delete):
                        del spheres[b[0]]
                        for later in to_delete[n + 1:]:  # renumber since one was deleted
                            if later[0] > b[0]:
                                later[0] -= 1
                z = z + gridsize
            y = y + gridsize
        x = x + gridsize

    # ... then remove spheres too close to each other (index juggling as in the original;
    # a negative index wraps around in both Perl and Python)
    i = 0
    while i < len(spheres):
        j = 0
        while j < len(spheres):
            si, sj = spheres[i], spheres[j]
            if j != i and _dist(_xyz(si), _xyz(sj)) < tooclose:
                if sj[6] != CRYSTAL and (sj[7] == NEITHER or sj[7] == NONPOLAR):
                    delete_j = True
                elif si[6] != CRYSTAL and (si[7] == NEITHER or sj[7] == NONPOLAR):  # sic: sj
                    delete_j = False
                elif sj[6] != CRYSTAL and sj[9] >= si[9]:
                    delete_j = True
                elif si[6] != CRYSTAL and si[9] > sj[9]:
                    delete_j = False
                elif sj[6] != CRYSTAL:
                    delete_j = True
                elif si[6] != CRYSTAL:
                    delete_j = False
                else:  # both crystallographic: keep both
                    delete_j = None
                if delete_j is True:
                    del spheres[j]
                    j -= 1
                    if i > j:
                        i -= 1
                elif delete_j is False:
                    del spheres[i]
                    i -= 1
                    if j > i:
                        j -= 1
            j += 1
        i += 1

    # 4. keep spheres continuous with the crystallographic spheres or the center sphere
    if not use_ligand:
        closest = 0
        closest_dist = _dist(_xyz(spheres[0]), ligand_center)
        for i in range(1, len(spheres)):
            dist = _dist(_xyz(spheres[i]), ligand_center)
            if dist < closest_dist:
                closest = i
                closest_dist = dist
        spheres[closest][6] = CONTINUOUS if spheres[closest][6] != CRYSTAL else CRYSTAL_AND_CENTER

    continuity = CONTINUITY
    while True:
        num_kept = sum(1 for s in spheres if s[6] != PLAIN)
        repeats = 1
        while num_kept <= max_spheres and repeats <= NUM_CONTINUITY_REPEATS:  # num_kept isn't updated here (sic)
            for k in range(len(spheres)):
                i = 0
                while spheres[k][6] == PLAIN and i < len(spheres):
                    if (
                        i != k
                        and _dist(_xyz(spheres[k]), _xyz(spheres[i])) <= continuity
                        and spheres[i][6] >= CRYSTAL
                    ):
                        spheres[k][6] = CONTINUOUS
                    i += 1
            repeats += 1
        num_kept = sum(1 for s in spheres if s[6] != PLAIN)
        if num_kept < min_spheres:
            continuity = continuity + CONTINUITY_INCREMENT
        if num_kept >= min_spheres or continuity > MAX_CONTINUITY:
            break
    spheres = [s for s in spheres if s[6] != PLAIN]

    # 5. if too many, drop those with the largest weighted distance, polar spheres favored
    if len(spheres) > max_spheres:
        spheres_center = _mean([_xyz(s) for s in spheres])
        if use_ligand:
            for s in spheres:
                s[8] = min(_dist(p, _xyz(s)) for p in ligand_coords)
                s[8] = s[8] + _dist(ligand_center, _xyz(s)) / 8
                s[8] = s[8] + _dist(spheres_center, _xyz(s)) / 8
        else:
            receptor_center = _mean(receptor_atoms)
            for s in spheres:
                s[8] = _dist(receptor_center, _xyz(s))
                s[8] = s[8] + _dist(ligand_center, _xyz(s))
                s[8] = s[8] + _dist(spheres_center, _xyz(s))
        average_weight = sum(s[8] for s in spheres) / len(spheres)
        for s in spheres:
            if s[7] != NEITHER and s[7] != NONPOLAR:
                s[8] = s[8] + average_weight * POLAR_BENEFIT

        by_weight = sorted(spheres, key=lambda s: s[8], reverse=True)
        index = 0
        while len(by_weight) > max_spheres and index < len(by_weight):  # first only continuous spheres
            if by_weight[index][6] == CONTINUOUS:
                del by_weight[index]
                index -= 1
            index += 1
        while len(by_weight) > max_spheres:  # then any
            del by_weight[0]
        spheres = by_weight

    with open(out_path, "w", newline="\n") as f:
        f.write(SPH_HEADER)
        f.write("cluster     1   number of spheres in cluster %5.f\n" % len(spheres))
        for n, s in enumerate(spheres, start=1):
            f.write("%5.f%10.5f%10.5f%10.5f%8.3f%5.f 0  0          \n" % (9000 + n, s[1], s[2], s[3], s[4], s[5]))
