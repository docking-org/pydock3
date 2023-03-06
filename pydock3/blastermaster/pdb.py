import collections
import copy
import logging


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


ATOM_TO_RADIUS_DICT = {  # TODO: 'atom' should be substituted with a more accurate term
    "C": 1.9,
    "O": 1.6,
    "N": 1.65,
    "P": 1.9,
    "S": 1.9,
    "H": 0.0,
    "F": 0.0,
    "I": 0.0,
    "U": 0.0,
    "A": 0.0,
    "B": 0.0,
    "L": 0.0,
    "*": 0.0,
    "Z": 0.0,
    "D": 0.0,
    "K": 0.0,
    "M": 0.0,
}
# F should really be about 1.5 but not in fortran so not here
# note that changing these breaks compatibility with trisrf/meshsrf surface
# generation programs and does not actually affect the radii used in those
# processes. in other words don't change them. DON'T DO IT. it won't change
# the radii used AT ALL, it will just break things.

OCC_PLACE = 0  # the occupancy is first, then the bfactor
BFAC_PLACE = 1  # bfactor


class PDBColumns:
    DOTS = (34, 42, 50, 57, 63)  # should be the periods in the columns
    DELETE_FROM = 22  # where to delete spaces from


# sometimes you want to use an external radii file (with extreme caution).
def read_radii_file(radii_file_path):
    """reads file in "c           1.90" format and returns map from name to radius
    uses column specific format, res name ignored.
    atom__res_radius_
    01234567890123456789"""
    with open(radii_file_path, "r") as f:
        name_radius = {}  # store in dictionary
        for line in f:  # TODO: use readlines instead
            try:
                name = line[0:5].strip().upper()
                radius = float(line[10:16])
                name_radius[name] = radius
            except ValueError:  # ignore where there isn't a float
                pass
    return name_radius


class PDBData(object):
    """stores data for a pdb file consisting of atoms"""

    def __init__(
        self,
        pdb_file_path=None,
        het_only=False,
        atom_to_radius_dict_file_path=None,
        ignore_waters=True,
    ):
        """default constructor takes a pdb file path as input, other ways later"""

        if atom_to_radius_dict_file_path is None:
            atom_to_radius_dict = ATOM_TO_RADIUS_DICT
        else:
            atom_to_radius_dict = read_radii_file(
                atom_to_radius_dict_file_path
            )  # override the defaults
        self.__non_zero_radii_count = 0
        self.raw_data = []
        self.coords = []
        self.radii = []
        self.charges = []
        self.hydro_charges = []
        self.factors = []
        self.atoms = []
        self.residue_nums = []
        self.residue_names = []
        self.alt_chars = []
        self.chains = []
        self.model_nums = []  # keeps track of NMR models if present
        self.atom_to_raw = {}
        self.raw_to_atom = {}
        self.ignore_waters = ignore_waters
        if pdb_file_path is not None:
            with open(pdb_file_path, "r") as f:
                model_num = 0
                for line in f.readlines():
                    if line.startswith("MODEL"):
                        try:
                            model_num = int(
                                line.split()[1]
                            )  # [0] is MODEL, [1] is the number
                        except IndexError:
                            model_num = 0
                    if not het_only or line.startswith("HETATM"):
                        self.process_line(line, atom_to_radius_dict, model_num)

    def process_line(self, line, atom_to_radius_dict, model_number=0):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            name = line[12:16]
            if name[0] == " ":  # because apparently hetatm entries can start one col
                name = name[1:]  # before atom entries (for the atom name)
            else:
                try:
                    _ = int(name[0])
                    name = name[1:]  # otherwise would have triggered exception
                except ValueError:
                    pass  # first character is not a number
            alt_char = line[16]
            residue_name = line[17:20]
            if not (self.ignore_waters and (name == "HOH")):
                self.raw_data.append(line)
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                chain = line[21:22]
                self.coords.append((x, y, z))
                self.atoms.append(name)
                try:
                    radius = atom_to_radius_dict[name[0]]
                except KeyError:
                    radius = 0.0  # default is to ignore
                # sometimes hydrogens are formatted badly and don't have H in 1st column
                if (radius > 0.0) and (name.count("H") > 0):  # any H's are bad
                    radius = 0.0
                self.radii.append(radius)
                if radius > 0.0:
                    self.__non_zero_radii_count += 1
                self.model_nums.append(model_number)
                factor_strings = (line[55:60], line[61:67])
                try:
                    factors = (float(factor_strings[0]), float(factor_strings[1]))
                except ValueError:  # in case no factors or occupancies present
                    factors = (0.0, 0.0)  # set to zero for now
                self.residue_nums.append(int(line[22:26]))
                self.chains.append(chain)
                self.alt_chars.append(alt_char)  # alternate sidechain conf letters
                self.residue_names.append(residue_name)
                self.factors.append(factors)
                self.atom_to_raw[len(self.atoms) - 1] = len(self.raw_data) - 1
                self.raw_to_atom[len(self.raw_data) - 1] = len(self.atoms) - 1

    def write(self, pdb_file_path):
        """opens the file, writes to it, closes it"""
        with open(pdb_file_path, "w") as f:
            for raw_data_line in self.raw_data:
                if raw_data_line:
                    f.write(raw_data_line)
                    if not raw_data_line.endswith("\n"):
                        f.write("\n")

    def replace_hetatm_with_atom(self):
        """replaces 'HETATM' with 'ATOM  ' in the rawdataline"""
        for i, raw_data_line in enumerate(self.raw_data):
            if raw_data_line[0:6] == "HETATM":
                new_data_line = "ATOM  " + raw_data_line[6:]
                self.raw_data[i] = new_data_line

    def copy(self):
        new_pdb = PDBData()
        new_pdb.__non_zero_radii_count = 0
        new_pdb.raw_data = self.raw_data[:]  # copy everything
        new_pdb.coords = self.coords[:]
        new_pdb.atoms = self.atoms[:]
        new_pdb.radii = self.radii[:]
        new_pdb.charges = self.charges[:]
        new_pdb.hydro_charges = self.hydro_charges[:]
        new_pdb.factors = self.factors[:]
        new_pdb.residue_nums = self.residue_nums[:]
        new_pdb.chains = self.chains[:]
        new_pdb.residue_names = self.residue_names[:]
        new_pdb.alt_chars = self.alt_chars[:]
        new_pdb.model_nums = self.model_nums[:]
        new_pdb.atom_to_raw = self.atom_to_raw.copy()
        new_pdb.raw_to_atom = self.raw_to_atom.copy()
        return new_pdb

    def clear_factor(self, raw_data_index, which_factor=BFAC_PLACE):
        old_factors = self.factors[raw_data_index]
        if which_factor == BFAC_PLACE:
            new_factors = (old_factors[OCC_PLACE], 0.0)
        elif which_factor == OCC_PLACE:
            new_factors = (0.0, old_factors[BFAC_PLACE])
        self.update_factors(
            raw_data_index, new_factors
        )  # TODO: new_factors might be referenced before assignment

    def update_factors(self, raw_data_index, new_factors):
        """changes the bfactor and occupancy column. now uses longer format for
        occupancy than might be standard, to get more accurate charges"""
        self.factors[raw_data_index] = new_factors
        new_data = "%+2.2f %+2.4f " % (new_factors[0], new_factors[1])
        new_data = self.raw_data[raw_data_index][:55] + new_data + "\n"
        new_data_no_plus = new_data.replace("+", " ")
        self.raw_data[raw_data_index] = new_data_no_plus

    def clear_factors_residues(self, residue_numbers, matching=True):
        """copy pdb, for each residue in the list, remove the bfactor column.
        if matching is false, then remove the bfactor if the residue is not in list
        return new pdb"""
        new_pdb = self.copy()
        marked_for_removal = []
        for i, residue_num in enumerate(new_pdb.residue_nums):
            if matching and (residue_num not in residue_numbers):
                marked_for_removal.append(i)
            elif not matching and (residue_num in residue_numbers):
                marked_for_removal.append(i)
        for index in marked_for_removal:
            new_pdb.clear_factor(new_pdb.atom_to_raw[index])
        return new_pdb

    def remove_line(self, raw_data_index):
        self.raw_data[raw_data_index] = False

    def remove_all_hydrogens(self, res_list=None):
        """for each residue in the list, remove all the hydrogens. if no list given,
        delete all hydrogens in whole protein."""
        if res_list is None:
            res_list = []
        marked_for_removal = []
        for i, atom in enumerate(self.atoms):
            if atom[0] == "H":  # hydrogen atom
                residue_num = self.residue_nums[i]
                if (res_list is None) or (residue_num in res_list):
                    marked_for_removal.append(self.atom_to_raw[i])
        for index in marked_for_removal:
            self.remove_line(index)

    def remove_apolar_hydrogen(self, residue_code_to_polar_hydrogens_dict):
        """for removing all nonpolar hydrogens in a protein. uses residue_code_to_polar_hydrogens_dict as the
        dict of residue->atom names to decide which hydrogens to keep"""
        marked_for_removal = []
        for i, atom in enumerate(self.atoms):
            if atom[0] == "H":  # hydrogen atom
                residue_name = self.residue_names[i]
                if residue_name not in residue_code_to_polar_hydrogens_dict:
                    logger.exception(
                        f"ERROR: residue name unknown: {residue_name} {self.raw_data[i]}"
                    )
                    raise
                else:
                    allowed_hydrogens = residue_code_to_polar_hydrogens_dict[
                        residue_name
                    ]
                    if atom.strip() not in allowed_hydrogens:
                        marked_for_removal.append(self.atom_to_raw[i])
        for index in marked_for_removal:
            self.remove_line(index)
        # no return necessary, as self has been modified

    def remove_protons_for_covalent_docking(
        self, residue_num, residue_name, residue_atom_names
    ):
        """given a pdb number of the covalent residue, removes its proton(s) to allow for the
        covalent bond"""
        for residue_atom_name in residue_atom_names.split(","):
            residue_atom_index = self.get_index_by_residue_atom(
                residue_num, residue_name, residue_atom_name
            )
            self.remove_line(residue_atom_index)

    def update_one_residue_name(self, raw_data_index, new_name):
        """updates the data and raw line for a residue name"""
        self.residue_names[self.raw_to_atom[raw_data_index]] = new_name
        new_data = self.raw_data[raw_data_index][:17]
        new_data += new_name + self.raw_data[raw_data_index][20:]
        self.raw_data[raw_data_index] = new_data

    def replace_alt_chars(self, new_char):
        """updates the data and raw line for all alternate characters"""
        for i, _ in enumerate(self.raw_data):
            self.alt_chars[self.raw_to_atom[i]] = new_char
            new_data = self.raw_data[i][:16]
            new_data += new_char + self.raw_data[i][17:]
            self.raw_data[i] = new_data

    def delete_insertion_codes(self):
        """insertion codes are sometimes added as 61A for a residue num
        instead of just using 62 (because people want the numbering to line up
        with some other numbering). they cause problems for some people, so
        this is for removing them. they are in column 27"""
        for i, _ in enumerate(self.raw_data):
            new_data = self.raw_data[i][:26]
            new_data += " "
            new_data += self.raw_data[i][27:]
            self.raw_data[i] = new_data

    def rename_histidines(self):
        """renames histidines from HIS to HID, HIE, or HIP based on hydrogens"""
        # two passes, first find histidines named HIS, want to check all at once
        residue_sets = {}  # from (chain,res number)  to [indices]
        for i, _ in enumerate(self.raw_data):
            if self.residue_names[i] == "HIS":  # only care about unchanged his
                chain_residue_num = (
                    self.chains[i],
                    self.residue_nums[i],
                )
                if chain_residue_num not in list(residue_sets.keys()):
                    residue_sets[chain_residue_num] = []
                residue_sets[chain_residue_num].append(i)
        # second pass, decide whether each HIS is HID, HIE, or HIP
        for index_list in residue_sets.values():
            has_d, has_e, has_p = False, False, False
            for raw_data_index in index_list:
                if self.atoms[raw_data_index] == "HD1":
                    has_d = True
                if self.atoms[raw_data_index] == "HE2":
                    has_e = True
            if has_d and has_e:  # both protonated
                has_p = True
            for raw_data_index in index_list:
                name = "HID"  # default name is HID, HIS is never allowed.
                if has_d:
                    name = "HID"
                if has_e:
                    name = "HIE"
                if has_p:
                    name = "HIP"
                self.update_one_residue_name(raw_data_index, name)

    def rename_cysteines(self):
        """renames cysteines from CYS to CYX based on hydrogens"""
        # what about CYM (negative)?
        # two passes, first find cysteines named CYS, want to check all at once
        residue_sets = {}  # from (chain,res number)  to [indices]
        for i, _ in enumerate(self.raw_data):
            if self.residue_names[i] == "CYS":  # only care about unchanged his
                chain_residue_num = (
                    self.chains[i],
                    self.residue_nums[i],
                )
                if chain_residue_num not in list(residue_sets.keys()):
                    residue_sets[chain_residue_num] = []
                residue_sets[chain_residue_num].append(i)
        # second pass, decide whether each CYS is CYS or CYX
        for index_list in residue_sets.values():
            cyx = True
            for raw_data_index in index_list:
                if self.atoms[raw_data_index] == "HG ":
                    cyx = False
            for raw_data_index in index_list:
                name = "CYS"  # default name is CYS.
                if cyx:
                    name = "CYX"
                self.update_one_residue_name(raw_data_index, name)

    def update_one_chain(self, raw_data_index, new_chain=" "):
        """updates the data and raw line for a residue number"""
        self.chains[self.raw_to_atom[raw_data_index]] = new_chain
        new_data = self.raw_data[raw_data_index][:21]
        temp_data = new_chain
        new_data += temp_data + self.raw_data[raw_data_index][22:]
        self.raw_data[raw_data_index] = new_data

    def fix_chain_ids(self):
        """sometimes people screw up the chain id and make them all the same.
        go through the residues and if they go down in number then start a new
        chain id."""
        last_residue_num = -10000
        chains_used = []
        replacing_chains = False
        new_chain = False
        for i, residue_number in enumerate(self.residue_nums):
            if residue_number < last_residue_num:
                replacing_chains = True
                new_chain = chr(max(ord(self.chains[i]), max(chains_used)) + 1)
            if replacing_chains:
                self.update_one_chain(i, new_chain)
            if ord(self.chains[i]) not in chains_used:
                chains_used.append(ord(self.chains[i]))
            last_residue_num = residue_number

    def get_index_by_residue_atom(self, residue_num, res_code, atom_name):
        """gets the index matching the input data, returns false if no match"""
        atom_name_str = atom_name.strip()
        for i, atom in enumerate(self.atoms):
            if residue_num == self.residue_nums[i]:
                if res_code == self.residue_names[i]:
                    if atom_name_str == atom.strip():
                        return i
        return False  # not found

    def get_occupancy_residue(self, residue_num):
        """gets the occupancy for one residue number. just use first found."""
        for i, atom in enumerate(self.atoms):
            if self.residue_nums[i] == residue_num:
                return self.factors[i][OCC_PLACE]

    def is_most_occupied_residue_chain(self, residue_num, chain_id):
        """for all residue_nums, is the chain_id provided the most occupied one (True) or
        not, return False then"""
        highest_occupancy, highest_count = 0.0, 0
        for i, atom in enumerate(self.atoms):
            if self.residue_nums[i] == residue_num:
                if self.factors[i][OCC_PLACE] > highest_occupancy:
                    highest_occupancy, highest_count = (
                        self.factors[i][OCC_PLACE],
                        i,
                    )
        if self.alt_chars[highest_count] == chain_id:  # TODO: simplify this
            return True
        else:
            return False

    def select_most_occupied(self, exceptions=None, leave_alone=None):
        """for each residue with alternate positions, pick the most occupied
        position. if equal, break ties starting with the A position.
        exceptions can be a list of residues to pick the least occupied for.
        leave_alone is a list of residues that aren't touched at all."""
        if exceptions is None:
            exceptions = []
        if leave_alone is None:
            leave_alone = []
        same_atoms = collections.defaultdict(list)  # collect lists of atoms
        for i, atom in enumerate(self.atoms):
            if self.alt_chars[i] != " ":  # space means no alternate position
                if self.residue_nums[i] not in leave_alone:
                    truple = (
                        atom,
                        self.residue_names[i],
                        self.residue_nums[i],
                    )
                    same_atoms[truple].append(i)
        for atom_truple, atom_list in same_atoms.items():
            normal = True  # means picked most occupied
            if atom_truple[2] in exceptions:
                normal = False  # means picked least occupied for these residues
            most_occupied = None
            occupancy = 0.0
            if not normal:
                occupancy = 1.0
            for atom_count in atom_list:
                if normal and self.factors[atom_count][OCC_PLACE] > occupancy:
                    occupancy = self.factors[atom_count][OCC_PLACE]
                    most_occupied = atom_count
                elif self.factors[atom_count][OCC_PLACE] == occupancy:
                    if (
                        self.alt_chars[atom_count] < self.alt_chars[most_occupied]
                    ):  # A<B<C
                        occupancy = self.factors[atom_count][OCC_PLACE]
                        most_occupied = atom_count
                elif not normal and self.factors[atom_count][OCC_PLACE] < occupancy:
                    occupancy = self.factors[atom_count][OCC_PLACE]
                    most_occupied = atom_count
            # okay, now remove everything but the most_occupied
            for atom_count in atom_list:
                if atom_count != most_occupied:  # keep this one
                    self.remove_line(self.atom_to_raw[atom_count])

    def delete_alternates(self, only=None):
        """for each sidechain with alternate positions, delete them all.
        if only exists, only delete residues in the list of residue numbers"""
        if only is None:
            only = []
        same_atoms = collections.defaultdict(list)  # collect lists of atoms
        for i, atom in enumerate(self.atoms):
            # if self.alt_chars[atom_count] != ' ':  # space means no alternate position
            truple = (
                atom,
                self.residue_names[i],
                self.residue_nums[i],
            )
            same_atoms[truple].append(i)
        for atom_truple, atom_list in same_atoms.items():
            if (only is None) or (atom_truple[2] in only):
                # okay, now remove everything
                for atom_count in atom_list:
                    self.remove_line(self.atom_to_raw[atom_count])

    def delete_all_residues(self, leave_alone=None):
        """deletes all the atoms in the protein, except the residues in
        leave_alone"""
        if leave_alone is None:
            leave_alone = []
        same_atoms = collections.defaultdict(list)  # collect lists of atoms
        for i, atom in enumerate(self.atoms):
            if self.residue_nums[i] not in leave_alone:
                truple = (
                    atom,
                    self.residue_names[i],
                    self.residue_nums[i],
                )
                same_atoms[truple].append(i)
        for atom_truple, atom_list in same_atoms.items():
            for atom_count in atom_list:
                self.remove_line(self.atom_to_raw[atom_count])

    def get_alt_chars(self, residue_numbers=None):
        """for each residue is residue_numbers, return the list of alt chars
        (alternate conformations) seen."""
        if residue_numbers is None:
            residue_numbers = []
        return_alt_chars = set()
        for i, atom in enumerate(self.atoms):
            if self.residue_nums[i] in residue_numbers:  # only care about these
                if self.alt_chars[i] != " ":  # space means no alternate position
                    return_alt_chars.add(self.alt_chars[i])
        return list(return_alt_chars)

    def select_one_alt(self, residue_numbers=None, pick_alt_char=None):
        """for each residue is residue_numbers, salvage only the pick_alt_char
        alternate conformation, delete other conformations."""
        if residue_numbers is None:
            residue_numbers = []
        same_atoms = collections.defaultdict(list)  # collect lists of atoms
        for i, atom in enumerate(self.atoms):
            if self.alt_chars[i] != " ":  # space means no alternate position
                truple = (
                    atom,
                    self.residue_names[i],
                    self.residue_nums[i],
                )
                same_atoms[truple].append(i)
        for atom_truple, atom_list in same_atoms.items():
            if atom_truple[2] in residue_numbers:  # means actually pick here
                for atom_count in atom_list:
                    if self.alt_chars[atom_count] != pick_alt_char:
                        self.remove_line(self.atom_to_raw[atom_count])

    def residue_sets(self):
        """calculates and returns a list of residue sets"""
        residue_sets = {}  # from (chain,res number)  to [indices]
        for i, _ in enumerate(self.raw_data):  # TODO: use itertools.product instead
            chain_residue_num = (self.chains[i], self.residue_nums[i])
            if chain_residue_num not in list(residue_sets.keys()):
                residue_sets[chain_residue_num] = []
            residue_sets[chain_residue_num].append(i)
        return residue_sets


def specific_alts(pdb_file_path, res_list, outfile_path):
    pdb_entry = PDBData(pdb_file_path)
    new_pdb = pdb_entry.copy()
    residue_num_list = [int(residue[:-1]) for residue in res_list]
    new_pdb.select_most_occupied(leave_alone=residue_num_list)
    # delete other alternates for other positions
    for residue in res_list:
        residue_num = int(residue[:-1])
        desired_char = residue[-1]
        new_pdb.select_one_alt([residue_num], desired_char)
    new_pdb.write(outfile_path)


def del_all_but(pdb_file_path, outfile_path, save_list=None):
    if save_list is None:
        save_list = []
    pdb_entry = PDBData(pdb_file_path)
    new_pdb = pdb_entry.copy()
    new_pdb.delete_all_residues(leave_alone=save_list)
    # delete other alternates for other positions
    new_pdb.write(outfile_path)


def most_occupied(pdb_file_path, outfile_path, exceptions=None):
    if exceptions is None:
        exceptions = []
    pdb_entry = PDBData(pdb_file_path)
    new_pdb = copy.deepcopy(pdb_entry)
    new_pdb.select_most_occupied(exceptions=exceptions)
    new_pdb.write(outfile_path)


def make_alts(pdb_file_path, outfile_prefix_path, res_list_list=None):
    if res_list_list is None:
        res_list_list = []
    pdb_entry = PDBData(pdb_file_path)
    return_paths = []
    return_chars = []
    for res_list in res_list_list:
        all_alt_chars = pdb_entry.get_alt_chars(res_list)
        all_alt_chars.sort()
        for desired_char in all_alt_chars:
            new_pdb = pdb_entry.copy()
            new_pdb.select_most_occupied(leave_alone=res_list)
            # delete other alternates for other positions
            new_pdb.select_one_alt(res_list, desired_char)
            out_res_list = str(res_list)[1:-1]  # we don't want [] chars
            out_res_list = "".join(out_res_list.split())  # remove spaces
            outfile_path = f"{outfile_prefix_path}.{out_res_list}.{desired_char}.pdb"
            new_pdb.write(outfile_path)
            return_paths.append(outfile_path)
            return_chars.append(desired_char)
    return return_paths, return_chars


def del_hydrogens(pdb_file_path, outfile_path, del_list=None):
    """for each residue in del_list, delete ALL hydrogens"""
    if del_list is None:
        del_list = []
    pdb_entry = PDBData(pdb_file_path)
    new_pdb = pdb_entry.copy()
    new_pdb.remove_all_hydrogens(res_list=del_list)
    new_pdb.write(outfile_path)


def delete_alts(pdb_file_path, outfile_path, only=None):
    pdb_entry = PDBData(pdb_file_path)
    new_pdb = pdb_entry.copy()
    if only is not None:
        new_pdb.delete_alternates(only)
    else:
        new_pdb.delete_alternates()
    new_pdb.write(outfile_path)


def delete_alt_chars(pdb_file_path, outfile_path):
    pdb_entry = PDBData(pdb_file_path)
    new_pdb = pdb_entry.copy()
    new_pdb.replace_alt_chars(" ")
    new_pdb.write(outfile_path)


def move_columns(
    input_pdb_file_path,
    output_pdb_file_path,
    dots=PDBColumns.DOTS,
    delete_from=PDBColumns.DELETE_FROM,
):
    lines_to_write = []
    with open(input_pdb_file_path, "r") as f_in:
        for line in f_in:
            if line.startswith("ATOM"):
                problems = 0
                for dot in dots:
                    try:
                        if line[dot] != ".":  # this is not good
                            problems += 1
                    except IndexError:  # might not be that long
                        pass  # no problem
                if problems == 0:  # great
                    lines_to_write.append(line)
                elif problems > 1:  # find where they should be
                    extra = 1
                    while (line[dots[0] + extra] != ".") and (
                        extra < 20
                    ):  # at 20, give up
                        extra += 1
                    if extra != 20:  # TODO: why 20? Make this a constant
                        new_line = line[:delete_from]
                        new_line += line[delete_from + extra :]
                    lines_to_write.append(new_line)
            else:
                lines_to_write.append(line)  # write these anyway

    with open(output_pdb_file_path, "w") as f_out:
        for line in lines_to_write:
            f_out.write(line)
