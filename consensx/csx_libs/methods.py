import os
import sys
import re
import math
import prody
import time
import pickle

import matplotlib.pyplot as plt
# installed modules
from . import objects as csx_obj


plt.switch_backend('Agg')


def timeit(method):
    """Timer decorator to keep an eye on CPU hungry processes"""
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print('\x1b[31m%r -> %2.2f sec\x1b[0m' % (method.__name__, te-ts),
              file=sys.stderr)
        return result

    return timed


def natural_sort(l):
    def convert(text):
        return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key):
        return [convert(c) for c in re.split('([0-9]+)', key)]

    return sorted(l, key=alphanum_key)


def check_3rd_party(install_dir):
    config_file_name = install_dir + "/.config"
    if os.path.isfile(config_file_name):
        csx_obj.ThirdParty.get_thirdparty(config_file_name)
    else:
        init_conf = (
            "# CoNSEnsX config file\n" +
            "# Please provide full paths\n" +
            "pales=''\n" + "shiftx=''\n" +
            "prideDB=''\n" + "prideNMR=''"
        )

        init_conf_file = open(".config", "w")
        init_conf_file.write(init_conf)
        print("Please edit '.config' in the CoNSEnsX install directory")
        raise SystemExit


@timeit
def get_model_list(PDB_file, my_path, model_count):
    """Parsing PDB file into models in the PDB_model object"""
    prody.confProDy(verbosity="info")

    atomgroup = prody.parsePDB(PDB_file, ter=True)
    num_coordsets = atomgroup.numCoordsets()

    # check for discarded models
    for model_num in range(num_coordsets):
        atomgroup.setACSIndex(model_num)

        if atomgroup.getCoords()[0][0] == float(0):
            print("DISCARDED MODEL FOUND")
            return False

    csx_obj.PDB_model(atomgroup, model_count)

    PDB_model_path = my_path + "/PDB_model.pickle"
    pickle.dump(csx_obj.PDB_model.model_data, open(PDB_model_path, 'wb'))
    return True


@timeit
def pdb_cleaner(my_path, PDB_file, my_CSV_buffer):
    """Performs some basic formatting on the given PDB file to make it suitable
       for further calculations (from RCSB to BMRB PDF format)"""
    try:
        input_pdb = open(PDB_file)
    except FileNotFoundError:
        print(PDB_file + " was not found")
        raise SystemExit

    work_PDB = my_path + "my_pdb.pdb"

    my_pdb = open(work_PDB, 'w')
    max_resnum = 0
    min_resnum = 100000

    for line in input_pdb:
        line = line.strip()
        # line = re.sub('[+-] ', '  ', line)

        if line.startswith("ATOM"):

            name = line[12:16].strip()
            resnum = int(line[22:26])

            if resnum >= max_resnum:
                max_resnum = resnum
            if resnum < min_resnum:
                min_resnum = resnum

            if name == "Q":
                continue
            if name == "NH":
                name = "H"
            if name == "HN":
                name = "H"

            chars, numbers = [], []
            for i in name:
                try:
                    numbers.append(int(i))
                except ValueError:
                    chars.append(i)

            # try:
            #     # _ = int(name[0])
            #     name = (''.join(str(i) for i in chars) +
            #             ''.join(str(i) for i in reversed(numbers)))
            # except ValueError:
            #     pass

            if len(name) == 1:
                name = " " + name + "  "
            elif len(name) == 2:
                name = " " + name + " "
            elif len(name) == 3:
                name = " " + name

            my_pdb.write(line[:11] + " %4s" % name + line[16:21] +
                         'A' + line[22:] + "\n")
            continue

        my_pdb.write(line + "\n")

    input_pdb.close()
    my_pdb.close()

    os.remove(PDB_file)
    os.rename(work_PDB, PDB_file)

    print("MAX RESNUM", max_resnum)
    my_CSV_buffer.max_resnum = max_resnum
    my_CSV_buffer.min_resnum = min_resnum


def pdb_splitter(my_path, PDB_file):
    """Split the given PDB file into models, each model becomes a separate
       PDB file placed in the "temp" folder"""
    try:
        my_pdb = open(PDB_file)
    except FileNotFoundError:
        print("ERROR -> the uploaded PDB file was not found")
        raise SystemExit
    except TypeError:
        print("ERROR -> PDB identifier is invalid")
        raise SystemExit

    model_names = []
    model_data = []

    my_name = ""
    my_data = []

    for line in my_pdb:
        if line.startswith("MODEL"):
            my_name = line.strip().split()[1]
        elif line.startswith("ATOM") or line.startswith("TER"):
            # replace oxygen names in line
            try:
                if line.split()[2] in ["OC1", "OT1"]:
                    line = line.replace(line.split()[2], "O ")
                elif line.split()[2] == "O1":
                    line = line.replace("O1 ", "O  ")
                elif line.split()[2] in ["OC2", "OT2"]:
                    line = line.replace(line.split()[2], "OXT")
                elif line.split()[2] == "O2":
                    line = line.replace("O2 ", "OXT")
            except IndexError:
                pass
            my_data.append(line.strip())
        elif line.startswith("ENDMDL"):
            model_names.append(my_name)
            model_data.append(my_data)
            my_name = ""
            my_data = []
        else:
            continue

    my_pdb.close()

    for i in range(len(model_names)):
        file_name = my_path + "model_" + model_names[i] + ".pdb"
        temp_pdb = open(file_name, 'w')
        temp_pdb.write("HEADER    MODEL " + model_names[i] + "\n")

        for _ in model_data[i]:
            temp_pdb.write(_ + "\n")

        temp_pdb.write("END")
        temp_pdb.close()

    return len(model_names)


def calcPeptideBonds():
    """Calculates backbone diherdral angles (OMEGA) CA-N-C'-CA"""
    # model_list = csx_obj.PDB_model.model_list
    dihedral_angles = {
        "<2": 0, "2-5": 0, "5-10": 0, "10-20": 0, ">20": 0
    }

    model_data = csx_obj.PDB_model.model_data

    for i in range(model_data.coordsets):
        model_data.atomgroup.setACSIndex(i)
        current_Resindex = 1
        prev_C, prev_CA, my_N, my_CA, my_C = None, None, None, None, None

        for atom in model_data.atomgroup:
            atom_res = atom.getResindex() + 1

            if atom_res != current_Resindex:

                if (prev_CA is not None and my_N is not None and
                        my_CA is not None and prev_C is not None):

                    NCA_vec = my_N - my_CA
                    CN_vec = prev_CA - my_N
                    CCA_vec = prev_C - my_CA

                    first_cross = csx_obj.Vec_3D.cross(CN_vec, NCA_vec)
                    second_cross = csx_obj.Vec_3D.cross(CCA_vec, NCA_vec)

                    angle = csx_obj.Vec_3D.dihedAngle(
                        first_cross, second_cross
                    )

                    # reference for setting sign of angle
                    reference = csx_obj.Vec_3D.cross(first_cross, second_cross)

                    r1 = reference.normalize()
                    r2 = NCA_vec.normalize()

                    if ((r1 - r2).magnitude() < r2.magnitude()):
                        angle *= -1

                    if abs(angle) < 2:
                        dihedral_angles["<2"] += 1
                    elif abs(angle) < 5:
                        dihedral_angles["2-5"] += 1
                    elif abs(angle) < 10:
                        dihedral_angles["5-10"] += 1
                    elif abs(angle) < 20:
                        dihedral_angles["10-20"] += 1
                    else:
                        dihedral_angles[">20"] += 1

                current_Resindex = atom_res
                prev_CA = my_CA
                prev_C = my_C
                my_N, my_CA, my_C = None, None, None

            if atom_res == current_Resindex:
                if atom.getName() == 'N':
                    my_N = csx_obj.Vec_3D(atom.getCoords())
                elif atom.getName() == 'CA':
                    my_CA = csx_obj.Vec_3D(atom.getCoords())
                elif atom.getName() == 'C':
                    my_C = csx_obj.Vec_3D(atom.getCoords())

    print("Peptide (CA-N-C'-CA) bond angle distribution:")
    print("   <2 -> " + str(dihedral_angles["<2"]))
    print("  2-5 -> " + str(dihedral_angles["2-5"]))
    print(" 5-10 -> " + str(dihedral_angles["5-10"]))
    print("10-20 -> " + str(dihedral_angles["10-20"]))
    print("  >20 -> " + str(dihedral_angles[">20"]))


def calcNH_Angles():
    """Calculates backbone diherdral angles (OMEGA) H-N-C=O"""
    model_data = csx_obj.PDB_model.model_data
    dihedral_angles = {"<2":    0, "2-5": 0, "5-10": 0,
                       "10-20": 0, ">20": 0}

    for i in range(model_data.coordsets):
        model_data.atomgroup.setACSIndex(i)
        current_Resindex = 1
        prev_O, prev_C, = None, None
        my_N, my_H, my_O, my_C = None, None, None, None

        for atom in model_data.atomgroup:
            atom_res = atom.getResindex() + 1

            if atom_res != current_Resindex:

                if (prev_O is not None and prev_C is not None and
                        my_N is not None and my_H is not None):

                    NH_vec = my_H - my_N
                    CN_vec = my_N - prev_C
                    OC_vec = prev_C - prev_O

                    first_cross = csx_obj.Vec_3D.cross(NH_vec, CN_vec)
                    second_cross = csx_obj.Vec_3D.cross(OC_vec, CN_vec)

                    angle = csx_obj.Vec_3D.dihedAngle(
                        first_cross, second_cross
                    )

                    # reference for setting sign of angle
                    reference = csx_obj.Vec_3D.cross(first_cross, second_cross)

                    r1 = reference.normalize()
                    r2 = NH_vec.normalize()

                    if ((r1 - r2).magnitude() < r2.magnitude()):
                        angle *= -1

                    if abs(angle) < 2:
                        dihedral_angles["<2"] += 1
                    elif abs(angle) < 5:
                        dihedral_angles["2-5"] += 1
                    elif abs(angle) < 10:
                        dihedral_angles["5-10"] += 1
                    elif abs(angle) < 20:
                        dihedral_angles["10-20"] += 1
                    else:
                        dihedral_angles[">20"] += 1

                current_Resindex = atom_res
                prev_O = my_O
                prev_C = my_C
                my_N, my_H = None, None

            if atom_res == current_Resindex:
                if atom.getName() == 'N':
                    my_N = csx_obj.Vec_3D(atom.getCoords())
                elif atom.getName() == 'H':
                    my_H = csx_obj.Vec_3D(atom.getCoords())
                elif atom.getName() == 'C':
                    my_C = csx_obj.Vec_3D(atom.getCoords())
                elif atom.getName() == 'O':
                    my_O = csx_obj.Vec_3D(atom.getCoords())

    print("Peptide (H-N-C=O) bond angle distribution:")
    print("   <2 -> " + str(dihedral_angles["<2"]))
    print("  2-5 -> " + str(dihedral_angles["2-5"]))
    print(" 5-10 -> " + str(dihedral_angles["5-10"]))
    print("10-20 -> " + str(dihedral_angles["10-20"]))
    print("  >20 -> " + str(dihedral_angles[">20"]))


def calcCorrel(calced, experimental):
    """Calculates correlation between calculated and experimental data
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    M = [0.0, 0.0, 0.0]
    D = [0.0, 0.0]
    match_count = 0

    for i in experimental:
        exp = i.value
        try:
            calc = calced[i.resnum]
        except KeyError:
            continue

        M[0] += calc
        M[1] += exp
        M[2] += calc * exp

        match_count += 1

    M[0] /= match_count
    M[1] /= match_count
    M[2] /= match_count

    for i in experimental:
        exp = i.value
        try:
            calc = calced[i.resnum]
        except KeyError:
            continue

        D[0] += (calc - M[0]) ** 2
        D[1] += (exp - M[1]) ** 2

    D[0] /= match_count
    D[0] = math.sqrt(D[0])
    D[1] /= match_count
    D[1] = math.sqrt(D[1])

    if D[0] * D[1] == 0:
        return -2
    else:
        return (M[2] - (M[0] * M[1])) / (D[0] * D[1])


def calcQValue(calced, experimental):
    """Calculates Q-value for calculated and experimental data
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    D2, E2 = 0, 0

    for i in experimental:
        exp = i.value
        try:
            calc = calced[i.resnum]
        except KeyError:
            continue

        D2 += (calc - exp) ** 2
        E2 += exp ** 2

    Q = 100 * math.sqrt(D2) / math.sqrt(E2)
    return round(Q, 6)


def calcRMSD(calced, experimental):
    """Calculates RMSD for calculated and experimental data
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    D2 = 0
    match_count = 0

    for i in experimental:
        exp = i.value
        try:
            calc = calced[i.resnum]
        except KeyError:
            continue

        D2 += (calc - exp) ** 2

        match_count += 1

    RMSD = math.sqrt(D2 / match_count)
    return round(RMSD, 6)
