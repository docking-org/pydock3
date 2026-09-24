"""Conversions between PDB and DOCK sphere files.

Python ports of the Fortran programs pdbtosph (ECM, 1992) and showsphere (S. Oatley,
E. Meng, D. Gschwend, 1994). They read fixed columns into single-precision REALs and write
them with fixed formats; the ports do the same (via float32) so the output is identical.
"""
import numpy as np

LIGAND_SPHERE_RADIUS = 0.70


def _real(s):
    """A Fortran REAL read from a fixed-width field (blank is 0)."""
    return float(np.float32(s.strip() or 0))


def _integer(s):
    """A Fortran INTEGER read from a fixed-width field (blank is 0)."""
    return int(s.strip() or 0)


def _fortran(fmt, value, width):
    """`value` formatted like a Fortran edit descriptor of `width`: all '*' if it doesn't fit."""
    s = fmt % value
    return s if len(s) <= width else "*" * width


def _records(path):
    """Lines as Fortran reads them with format (A80): without the newline, padded to 80."""
    with open(path) as f:
        return [line.rstrip("\r\n")[:80].ljust(80) for line in f]


def pdb_to_sph(pdb_path, sph_path):
    """Write the ligand's heavy atoms as a sphere cluster (pdbtosph)."""
    sphere_lines = []
    for line in _records(pdb_path):
        if line[0:4] == "ATOM" or line[0:6] == "HETATM":
            if line[13] in "HD" or line[12] in "HD":  # hydrogen / deuterium
                continue
            atom_num = _integer(line[6:11])
            x, y, z = _real(line[30:38]), _real(line[38:46]), _real(line[46:54])
            sphere_lines.append(
                _fortran("%5d", atom_num, 5)
                + "".join(_fortran("%10.5f", v, 10) for v in (x, y, z))
                + _fortran("%8.3f", _real(str(LIGAND_SPHERE_RADIUS)), 8)
                + _fortran("%5d", atom_num, 5)
            )
    with open(sph_path, "w", newline="\n") as f:
        f.write("cluster     1   number of spheres in cluster" + _fortran("%6d", len(sphere_lines), 6) + "\n")
        for sphere_line in sphere_lines:  # copied through an (A80) record, so blank-padded to 80
            f.write(sphere_line.ljust(80) + "\n")


def sph_to_pdb(sph_path, cluster, pdb_path):
    """Write the spheres of `cluster` as PDB atoms (showsphere, without surfaces)."""
    records = _records(sph_path)
    for n, line in enumerate(records):
        if line[0:7] == "cluster" and _integer(line[8:13]) == cluster:
            num_spheres = _integer(line[45:50])
            with open(pdb_path, "w", newline="\n") as f:
                for sphere_line in records[n + 1: n + 1 + num_spheres]:
                    sphere_num = _integer(sphere_line[0:5])
                    xyz = (_real(sphere_line[5:15]), _real(sphere_line[15:25]), _real(sphere_line[25:35]))
                    f.write(
                        "ATOM   " + _fortran("%4d", sphere_num, 4) + "  C   SPH  " + _fortran("%4d", sphere_num, 4) + "    "
                        + "".join(_fortran("%8.3f", v, 8) for v in xyz) + "\n"
                    )
                    f.write("TER\n")
            return
    raise ValueError(f"Cluster {cluster} not found in {sph_path}")
