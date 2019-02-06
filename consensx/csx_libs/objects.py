import math
import copy


class CalcPickle(dict):
    pass
    """Class for storing values for pickle generation"""
    data = dict()


class ThirdParty(object):

    """Class to store 3rd party software information"""
    pales = ""
    shiftx = ""
    prideDB = ""
    prideNMR = ""

    @staticmethod
    def get_thirdparty(config_file):
        cfg = open(config_file)
        for line in cfg:
            if line.startswith("#"):
                continue
            if line.startswith("pales"):
                ThirdParty.pales = line.split("'")[1]
            elif line.startswith("shiftx"):
                ThirdParty.shiftx = line.split("'")[1]
            elif line.startswith("prideDB"):
                ThirdParty.prideDB = line.split("'")[1]
            elif line.startswith("prideNMR"):
                ThirdParty.prideNMR = line.split("'")[1]


class CSV_buffer(object):

    """Class which stores data for values.CSV"""

    def __init__(self, my_path):
        self.working_dir = my_path
        self.max_resnum = -1
        self.min_resnum = 100000
        self.csv_data = []
        # self.name = name
        # self.calced = calced
        # self.exp = {}

        # for i in experimental:
        #     self.exp[i.resnum] = i.value

        # csv_data.append(self)

    def writeCSV(self):
        filename = self.working_dir + "values.csv"
        output_csv = open(filename, 'w')
        output_csv.write(',')
        for data in self.csv_data:
            output_csv.write(data["name"] + " EXP, " + data["name"] + " CALC,")
        output_csv.write("\n")

        # for i in experimental:
        #     self.exp[i.resnum] = i.value

        for resnum in range(self.min_resnum, self.max_resnum + 1):
            output_csv.write(str(resnum) + ',')
            for data in self.csv_data:
                exp = {}

                for i in data["experimental"]:
                    exp[i.resnum] = i.value

                try:
                    # pdb.set_trace()
                    output_csv.write(
                        "{0:.2f}".format(exp[resnum]) + ',' +
                        "{0:.2f}".format(data["calced"][resnum]) + ','
                    )
                except (IndexError, KeyError):
                    output_csv.write(',,')

            output_csv.write("\n")


class ChemShift_modell_data(object):

    """Class for per model chemical shift data"""
    type_dict = []

    @staticmethod
    def get_type_data(my_type):
        type_data = []

        for model in ChemShift_modell_data.type_dict:
            type_data.append(model[my_type])

        return type_data


class PDB_model(object):

    """Class for storing PDB model data"""
    is_fitted = False
    model_data = None

    def __init__(self, atomgroup, model_count):
        self.atomgroup = atomgroup
        self.model_count = model_count
        try:
            self.elements = atomgroup.getElements()
        except AttributeError:
            print("ERROR -> PDB parsing failed. Please check your PDB file!")
            raise SystemExit
        self.names = atomgroup.getNames()
        self.coordsets = atomgroup.numCoordsets()
        PDB_model.model_data = self


class Restraint_Record(object):

    """Class for storing NOE restraint data"""
    all_restraints = []
    resolved_restraints = []

    def __init__(self, csx_id, curr_distID, seq_ID1, seq_ID2,
                 seq_name1, seq_name2, atom_ID1, atom_ID2, dist_max):

        self.csx_id = int(csx_id)
        self.curr_distID = int(curr_distID)
        self.seq_ID1 = int(seq_ID1)
        self.seq_ID2 = int(seq_ID2)
        self.seq_name1 = str(seq_name1)
        self.seq_name2 = str(seq_name2)
        self.dist_max = float(dist_max)
        self.atom_ID1 = atom_ID1
        self.atom_ID2 = atom_ID2

        Restraint_Record.all_restraints.append(copy.deepcopy(self))

    @staticmethod
    def getRestraintCount():
        return Restraint_Record.all_restraints[-1].csx_id

    @staticmethod
    def getNOERestraints():
        # avoid repeated resolution of pseudo atom names
        if Restraint_Record.resolved_restraints:
            return Restraint_Record.resolved_restraints

        pse_c = ["#", "*", "%", "+"]

        restraints = Restraint_Record.all_restraints

        # load atom names from PDB data into a dict, residue IDs as keys
        NOE_dict = {}
        PDB_atom_names = {}
        model_data = PDB_model.model_data

        for atom in model_data.atomgroup:
            atom_res = atom.getResnum()

            if atom_res in PDB_atom_names.keys():
                PDB_atom_names[atom_res].append(atom.getName())
            else:
                PDB_atom_names[atom_res] = [atom.getName()]

        # organize restraint into a dict, with distIDs as keys
        for res in restraints:
            if res.csx_id in NOE_dict.keys():
                NOE_dict[res.csx_id].append(res)
            else:
                NOE_dict[res.csx_id] = [res]

        # resolution of pseudo-atoms
        for rest_id in NOE_dict.keys():
            for res in NOE_dict[rest_id]:
                # skip O atoms
                if res.atom_ID1 == "O" or res.atom_ID2 == "O":
                    continue

                resol1 = [res.atom_ID1]

                # check if first atom is a pseudo atom
                pseudH1 = res.atom_ID1[0] == "H" and res.atom_ID1[-1] in pse_c

                if res.atom_ID1[0] in "QM" or pseudH1:
                    # change first letter of the name to hydrogen
                    # WASSERSTOFF, weil es im Wasser ist! Merkst was?
                    if res.atom_ID1[0] in "QM":
                        base = 'H' + res.atom_ID1[1:]

                        # separate trailing number, if any
                        try:
                            num = str(int(res.atom_ID1[-1]))
                            base = 'H' + res.atom_ID1[1:-1]
                        except ValueError:
                            num = 0

                        atom_names1 = []

                        # permutations, if no trailing number found
                        if num != 0:
                            for i in ['1', '2', '3']:
                                atom_names1.append(base + num + i)
                                atom_names1.append(num + base + i)
                        # permutations, if trailing number found
                        else:
                            if res.seq_name1 == "ILE" and res.atom_ID1 == "MG":
                                for i in ['1', '2', '3']:
                                    atom_names1.append(base + str(2) + i)
                                    atom_names1.append(i + base + str(2))

                            elif (res.seq_name1 == "THR" and
                                  res.atom_ID1 == "MG"):
                                atom_names1.append("HG1")
                                atom_names1.append("1HG")

                            else:
                                for i in ['1', '2', '3']:
                                    atom_names1.append(base + i)

                                    for j in ['1', '2', '3']:
                                        atom_names1.append(base + j + i)
                                        # vagy fordítva
                                        atom_names1.append(i + base + j)

                    if pseudH1:
                        atom_names1 = []
                        base = res.atom_ID1[:-1]

                        for i in ['1', '2', '3']:
                            atom_names1.append(base + i)
                            atom_names1.append(i + base)

                    # TODO KeyError

                    # get corresponding atom names present in PDB file
                    PDB_names1 = PDB_atom_names[res.seq_ID1]
                    # print("atom_names1:", res.atom_ID1, atom_names1)
                    resol1 = list(set(atom_names1) & set(PDB_names1))
                    # print("resol_names1:", res.atom_ID1, resol1)

                resol2 = [res.atom_ID2]

                # check if second atom is a pseudo atom
                pseudH2 = res.atom_ID2[0] == "H" and res.atom_ID2[-1] in pse_c

                if res.atom_ID2[0] in "QM" or pseudH2:
                    atom_names2 = []

                    if res.atom_ID2[0] in "QM":
                        base = 'H' + res.atom_ID2[1:]

                        try:
                            num = str(int(res.atom_ID2[-1]))
                            base = 'H' + res.atom_ID2[1:-1]
                        except ValueError:
                            num = 0

                        if num != 0:
                            for i in ['1', '2', '3']:
                                atom_names2.append(base + num + i)
                                atom_names2.append(num + base + i)
                        else:
                            if res.seq_name2 == "ILE" and res.atom_ID2 == "MG":
                                for i in ['1', '2', '3']:
                                    atom_names2.append(base + str(2) + i)
                                    atom_names2.append(i + base + str(2))

                            elif (res.seq_name2 == "THR" and
                                  res.atom_ID2 == "MG"):
                                atom_names2.append("HG1")
                                atom_names2.append("1HG")

                            else:
                                for i in ['1', '2', '3']:
                                    atom_names2.append(base + i)

                                    for j in ['1', '2', '3']:
                                        atom_names2.append(base + j + i)
                                        # vagy fordítva
                                        atom_names2.append(i + base + j)

                    if pseudH2:
                        atom_names2 = []
                        base = res.atom_ID2[:-1]

                        for i in ['1', '2', '3']:
                            atom_names2.append(base + i)
                            atom_names2.append(i + base)

                    PDB_names2 = PDB_atom_names[res.seq_ID2]
                    # print("atom_names2:", res.atom_ID2, atom_names2)
                    resol2 = list(set(atom_names2) & set(PDB_names2))
                    # print("resol_names2:", res.atom_ID2, resol2)

                for atom1 in resol1:
                    for atom2 in resol2:
                        new = Restraint_Record(
                            res.csx_id, res.curr_distID,
                            res.seq_ID1, res.seq_ID2, res.seq_name1,
                            res.seq_name2, atom1, atom2, res.dist_max
                        )

                        # print(res.csx_id, res.curr_distID,
                        #     res.seq_ID1, res.seq_name1, atom1,
                        #     res.seq_ID2, res.seq_name2, atom2,
                        #     res.dist_max)

                        Restraint_Record.resolved_restraints.append(new)

        return Restraint_Record.resolved_restraints

    @staticmethod
    def getPRIDE_restraints():
        PRIDE_restraints = {}
        prev_id = -1
        seq1_ok, seq2_ok, distance_ok = False, False, False
        seq_dist = -1

        for restraint in Restraint_Record.resolved_restraints:
            curr_id = restraint.curr_distID
            atom_ID1 = restraint.atom_ID1
            atom_ID2 = restraint.atom_ID2
            seq_ID1 = restraint.seq_ID1
            seq_ID2 = restraint.seq_ID2

            if prev_id != curr_id:
                prev_id = curr_id

                if seq1_ok and seq2_ok and distance_ok:
                    if seq_dist in PRIDE_restraints:
                        PRIDE_restraints[seq_dist] += 1
                    else:
                        PRIDE_restraints[seq_dist] = 1

                atom1_HA = atom_ID1 in ["H", "HA", "HA1", "HA2", "HA3"]
                seq1_ok = atom1_HA or atom_ID1.startswith("HB")
                atom2_HA = atom_ID2 in ["H", "HA", "HA1", "HA2", "HA3"]
                seq2_ok = atom2_HA or atom_ID2.startswith("HB")
                seq_dist = abs(seq_ID1 - seq_ID2)
                distance_ok = seq_dist > 2

                id_distance = seq_dist

            else:
                atom1_HA = atom_ID1 in ["H", "HA", "HA1", "HA2", "HA3"]
                seq1_ok &= atom1_HA or atom_ID1.startswith("HB")
                atom2_HA = atom_ID2 in ["H", "HA", "HA1", "HA2", "HA3"]
                seq2_ok &= atom2_HA or atom_ID2.startswith("HB")
                seq_dist = abs(seq_ID1 - seq_ID2)
                distance_ok &= seq_dist > 2 and id_distance == seq_dist

        return PRIDE_restraints


class Vec_3D(object):

    """Vector class for calculations"""
    def __init__(self, v):
        self.v = v

    def __getitem__(self, index):
        return self.v[index]

    def __len__(self):
        return len(self.v)

    def __sub__(self, other):
        v = self.v
        u = other.v
        return Vec_3D([u[i] - v[i] for i in range(len(u))])

    def __iadd__(self, other):
        v = self.v
        u = other.v
        return Vec_3D([u[i] + v[i] for i in range(len(u))])

    def __idiv__(self, divisor):
        v = self.v
        return Vec_3D([v[i] / divisor for i in range(len(v))])

    def magnitude(self):
        v = self.v
        return math.sqrt(sum(v[i] * v[i] for i in range(len(v))))

    def normalize(self):
        vmag = self.magnitude()
        v = self.v
        return Vec_3D([v[i] / vmag for i in range(len(v))])

    @classmethod
    def cross(cls, one, other):
        c = [one.v[1] * other.v[2] - one.v[2] * other.v[1],
             one.v[2] * other.v[0] - one.v[0] * other.v[2],
             one.v[0] * other.v[1] - one.v[1] * other.v[0]]

        return cls(c)

    @staticmethod
    def dihedAngle(one, other):
        calc_cos = (
            one.v[0] * other.v[0] +
            one.v[1] * other.v[1] +
            one.v[2] * other.v[2]
            ) / (one.magnitude() * other.magnitude())

        # minor correction due to numerical issues
        if calc_cos > 1:
            calc_cos = 1

        return math.degrees(math.acos(calc_cos))
