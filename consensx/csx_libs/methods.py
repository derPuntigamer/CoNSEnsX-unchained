import os
import sys
import re
import subprocess
import math
import prody
import time
import pickle
import numpy as np

import matplotlib.pyplot as plt
# installed modules
import consensx.nmrpystar as nmrpystar
from . import objects as csx_obj


plt.switch_backend('Agg')

shortcodes = {
    'ALA': 'A', 'ASP': 'D', 'ASN': 'N', 'ARG': 'R', 'CYS': 'C', 'GLY': 'G',
    'GLU': 'E', 'GLN': 'Q', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V'
}


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


def get_PDB(args):
    """Gets PDB file or downloads PDF file from rcsb.org"""
    if args.PDB_file:
        my_PDB = args.PDB_file
    else:
        my_PDB = prody.fetchPDB(args.PDB_fetch, compressed=False)
        print()
    return my_PDB


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


@timeit
def parseSTR(STR_file):
    """Parse BMRB file into a python object"""

    star_file = open(STR_file)         # open STR file
    myString = ""
    parse_exception = (
        "ERROR during STR parsing, please check your STAR-NMR file!"
    )

    try:
        for line in star_file:                  # rean STR file into a string
            myString += line
    except UnicodeDecodeError:
        raise Exception(
            'ERROR during reading STAR-NMR file, please use ASCII encoding '
        )

    star_file.close()

    try:
        parsed = nmrpystar.parse(myString)    # parsing -> parsed.value

        if parsed.status != 'success':        # check if parsing was successful
            raise Exception(parse_exception)
    except KeyError:
        raise Exception(parse_exception)
    else:
        return parsed


def get_RDC_lists(parsed_value):
    """Returns RDC lists as dictonaries containing RDC_Record objects,
       grouped by RDCtype (keys())"""
    list_number = 1
    RDC_lists = []

    while True:
        saveShiftName = 'RDC_list_' + str(list_number)
        try:
            saveShifts = parsed_value.saves[saveShiftName]
        except KeyError:
            break
        loopShifts = saveShifts.loops[-1]
        RDC_records = []

        # STR key values recognised by this program
        rdc_res1_keys = ["RDC.Seq_ID_1", "Atom_one_residue_seq_code"]
        rdc_atom1_keys = ["RDC.Atom_type_1", "Atom_one_atom_name"]
        rdc_res2_keys = ["RDC.Seq_ID_2", "Atom_two_residue_seq_code"]
        rdc_atom2_keys = ["RDC.Atom_type_2", "Atom_two_atom_name"]
        rdc_value_keys = ["RDC.Val", "Residual_dipolar_coupling_value"]

        for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            for my_resnum1 in rdc_res1_keys:     # fetch 1. residue number
                if my_resnum1 in list(row.keys()):
                    resnum1 = row[my_resnum1]

            for my_atom1 in rdc_atom1_keys:      # fetch 1. atom name
                if my_atom1 in list(row.keys()):
                    atom1 = row[my_atom1]

            for my_resnum2 in rdc_res2_keys:     # fetch 2. residue number
                if my_resnum2 in list(row.keys()):
                    resnum2 = row[my_resnum2]

            for my_atom2 in rdc_atom2_keys:      # fetch 2. atom name
                if my_atom2 in list(row.keys()):
                    atom2 = row[my_atom2]

            for my_RDC_value in rdc_value_keys:  # fetch RDC value
                if my_RDC_value in list(row.keys()):
                    RDC_value = float(row[my_RDC_value])
                    RDC_value = float("{0:.2f}".format(RDC_value))

            # check if all parameters are fetched
            if (resnum1 and atom1 and resnum2 and atom2 and RDC_value):
                # append RDC_Record object to list
                RDC_records.append(
                    csx_obj.RDC_Record(
                        resnum1, atom1, resnum2, atom2, RDC_value
                    )
                )
            else:
                print(row)

        RDC_lists.append(RDC_records)
        list_number += 1

    # split list into dict according to RDC types
    new_RDC_list = []
    for RDC_list in RDC_lists:
        prev_type = ""
        RDC_dict = {}

        for record in RDC_list:
            if prev_type != record.RDC_type:
                RDC_dict[record.RDC_type] = []
                RDC_dict[record.RDC_type].append(record)
            else:
                RDC_dict[record.RDC_type].append(record)

            prev_type = record.RDC_type

        new_RDC_list.append(RDC_dict)

    return new_RDC_list


def parseS2_STR(parsed_value):
    """Returns a dictonary with the parsed S2 data"""
    try:
        saveShifts = parsed_value.saves["order_param"]

        loopShifts = saveShifts.loops[-1]
        S2_records = []

        for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            S2_value = float("{0:.2f}".format(float(row["S2_value"])))

            S2_records.append(
                csx_obj.S2_Record(
                    row["Residue_seq_code"], row["Atom_name"], S2_value)
            )

        # split list into dict according to S2 types
        S2_dict = {}
        prev_type = ""

        for record in S2_records:
            if prev_type != record.type:
                S2_dict[record.type] = []
                S2_dict[record.type].append(record)
            else:
                S2_dict[record.type].append(record)

            prev_type = record.type

        return S2_dict

    except KeyError:
        print("No S2 parameter list found")
        return None


def parse_sidechain_S2_STR(parsed_value):
    """Returns a dictonary with the parsed S2 data"""
    try:
        saveShifts = parsed_value.saves["side-chain_methyl_order_parameters"]

        loopShifts = saveShifts.loops[-1]
        S2_records = []

        for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            S2_value = float("{0:.2f}".format(float(row["S2_value"])))

            S2_records.append(
                csx_obj.S2_Record(
                    row["Residue_seq_code"], row["Atom_name"], S2_value
                )
            )

        return S2_records

    except KeyError:
        return None


def parseJcoup_STR(parsed_value):
    """Returns a dictonary with the parsed J-coupling data"""
    try:
        saveShifts = parsed_value.saves["coupling_constants"]

        loopShifts = saveShifts.loops[-1]
        jcoup_records = []

        for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            JC_value = row["Coupling_constant_value"]
            JC_value = float("{0:.2f}".format(float(JC_value)))

            jcoup_records.append(
                csx_obj.JCoup_Record(
                    row["Atom_one_residue_seq_code"],
                    row["Coupling_constant_code"],
                    JC_value
                )
            )

        # split list into dict according to J-cuopling types
        jcoup_dict = {}
        prev_type = ""

        for record in jcoup_records:
            if prev_type != record.type:
                jcoup_dict[record.type] = []
                jcoup_dict[record.type].append(record)
            else:
                jcoup_dict[record.type].append(record)

            prev_type = record.type

        return jcoup_dict

    except KeyError:
        print("No J-coupling parameter list found")
        return None


def parseChemShift_STR(parsed_value):
    """Returns ChemShift lists as dictonaries containing ChemShift_Record
    objects, grouped by Atom_name (keys())"""
    list_number = 1
    ChemShift_lists = []
    cs_lot = [
        {
            "seq_code": "Residue_seq_code",
            "label": "Residue_label",
            "name": "Atom_name",
            "value": "Chem_shift_value"
        },
        {
            "seq_code": "Atom_chem_shift.Seq_ID",
            "label": "Atom_chem_shift.Comp_ID",
            "name": "Atom_chem_shift.Atom_ID",
            "value": "Atom_chem_shift.Val"
        }
    ]

    while True:
        saveShiftName = 'chem_shift_list_' + str(list_number)
        altShiftName = 'assigned_chem_shift_list_' + str(list_number)
        try:
            saveShifts = parsed_value.saves[saveShiftName]
            str_type = 0
        except KeyError:
            try:
                saveShifts = parsed_value.saves[altShiftName]
                str_type = 1
            except KeyError:
                break

        loopShifts = saveShifts.loops[-1]
        ChemShift_records = []
        HA_sum = 0.0

        for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            CS_value = row[cs_lot[str_type]["value"]]
            CS_value = float("{0:.2f}".format(float(CS_value)))

            if row[cs_lot[str_type]["name"]] == "HA2":
                HA_sum += CS_value
                continue

            if row[cs_lot[str_type]["name"]] == "HA3":
                HA_sum += CS_value

                ChemShift_records.append(
                    csx_obj.ChemShift_Record(
                        row[cs_lot[str_type]["seq_code"]],
                        row[cs_lot[str_type]["label"]],
                        "HA", HA_sum / 2
                    )
                )
                HA_sum = 0.0
                continue

            chemshift_types = ["HA", "CA", "CB", "N", "H", "C"]

            if row[cs_lot[str_type]["name"]] in chemshift_types:
                ChemShift_records.append(
                    csx_obj.ChemShift_Record(
                        row[cs_lot[str_type]["seq_code"]],
                        row[cs_lot[str_type]["label"]],
                        row[cs_lot[str_type]["name"]],
                        CS_value
                    )
                )

            elif row[cs_lot[str_type]["name"]] == "HN":
                ChemShift_records.append(
                    csx_obj.ChemShift_Record(
                        row[cs_lot[str_type]["seq_code"]],
                        row[cs_lot[str_type]["label"]],
                        "H",
                        CS_value
                    )
                )

        ChemShift_lists.append(ChemShift_records)
        list_number += 1

    new_CS_list = []

    for ChemShift_list in ChemShift_lists:
        ChemShift_dict = {}

        for record in ChemShift_list:

            if record.atom_name in list(ChemShift_dict.keys()):
                ChemShift_dict[record.atom_name].append(record)
            else:
                ChemShift_dict[record.atom_name] = []
                ChemShift_dict[record.atom_name].append(record)

        new_CS_list.append(ChemShift_dict)

    return new_CS_list


@timeit
def callPalesOn(my_path, pdb_files, RDC_dict, lc_model, SVD_enable):
    """Writes pales dummy from the given RDC values, and call Pales with the
    given parameters"""

    pwd = os.getcwd()
    os.chdir(my_path)

    for o, pdb_file in enumerate(pdb_files):
        # ------------------  Open file and read PDB data  -------------------#
        try:
            input_pdb = open(pdb_file)
        except IOError:
            print("Input file " + pdb_file + " was not found")
            raise SystemExit                    # exit if input file not found

        seg = []                                # list storing sequence data

        for line in input_pdb:
            if line.startswith("ATOM") and line.split()[2] == "CA":
                resname = line.split()[3]       # get residue name
                seg.append(resname)             # append new segname to list

        input_pdb.close()

        # ----------------------  Write sequence data  -----------------------#
        short_seg = ""

        for i in range(len(seg)):
            short_seg += shortcodes[seg[i]]

        my_line = "DATA SEQUENCE "
        char_counter = 0
        row_counter = 0
        pales_dummy = open('pales_dummy.txt', 'w')

        for char in short_seg:
            if char_counter == 10:          # write aa output in 10 wide blocks
                my_line += " "
                char_counter = 0
                row_counter += 1

                if row_counter == 5:        # write 5 block per line
                    pales_dummy.write(my_line + "\n")
                    char_counter = 0
                    row_counter = 0
                    my_line = "DATA SEQUENCE "

            my_line += char
            char_counter += 1

        pales_dummy.write(my_line + "\n")    # write last line of aa output

        # ----------------------  Write dummy dipoles  -----------------------#
        pales_dummy.write(
            "\nVARS RESID_I RESNAME_I ATOMNAME_I " +
            "RESID_J RESNAME_J ATOMNAME_J D DD W\n" +
            "FORMAT %5d  %6s  %6s  %5d  %6s  %6s  %9.3f  %9.3f  %.2f \n\n"
        )

        lists = []
        for RDC_list in list(RDC_dict.keys()):
            lists.append(RDC_dict[RDC_list])

        for RDC_set in lists:
            for RDC_record in RDC_set:

                pales_dummy.write(
                    "%5s  %6s  %6s  %5s  %6s  %6s  %9.3f  %9.3f  %.2f\n" % (
                        str(RDC_record.resnum) + 'A',
                        seg[RDC_record.resnum - 1],
                        str(RDC_record.atom),
                        str(RDC_record.resnum2) + 'A',
                        seg[RDC_record.resnum2 - 1],
                        str(RDC_record.atom2),
                        RDC_record.value, 1.000,  1.00
                    )
                )

        pales_dummy.close()

        outfile = open("pales.out", 'a')
        DEVNULL = open(os.devnull, 'w')

        print("call Pales on model: " +
              str(o + 1) + '/' + str(len(pdb_files)), end="\r")
        sys.stdout.flush()

        try:
            if SVD_enable:                          # if SVD is enabled
                p = subprocess.Popen(
                    [
                        csx_obj.ThirdParty.pales,
                        "-inD", "pales_dummy.txt",
                        "-pdb", pdb_file,           # pdb file
                        '-' + lc_model,             # rdc lc model
                        "-bestFit"                  # SVD
                    ],
                    stdout=outfile,
                    stderr=DEVNULL
                )
                p.wait()                        # now wait
            else:                               # if SVD is disabled (default)
                subprocess.call(
                    [
                        csx_obj.ThirdParty.pales,
                        "-inD", "pales_dummy.txt",
                        "-pdb", pdb_file,           # pdb file
                        '-' + lc_model              # rdc lc model
                    ],
                    stdout=outfile,
                    stderr=DEVNULL
                )
        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)
        outfile.close()
        DEVNULL.close()

    print()
    os.chdir(pwd)


def calcS2(model_data, calculate_on_models,
           S2_records, S2_type, fit, fit_range):
    """Returns a dictonary with the average S2 values:
    S2_calced[residue] = value"""

    if fit:
        reference = model_data.atomgroup[:]

        model_data.atomgroup.setACSIndex(0)
        prody.alignCoordsets(model_data.atomgroup.calpha)

        if fit_range:
            for model_num in calculate_on_models:
                model_data.atomgroup.setACSIndex(model_num)

                mobile = model_data.atomgroup[:]
                matches = prody.matchChains(reference, mobile)
                match = matches[0]
                ref_chain = match[0]
                mob_chain = match[1]

                # if fit_range:
                weights = np.zeros((len(ref_chain), 1), dtype=np.int)

                fit_start, fit_end = fit_range.split('-')

                for i in range(int(fit_start) - 1, int(fit_end) - 1):
                    weights[i] = 1

                t = prody.calcTransformation(mob_chain, ref_chain, weights)
                t.apply(mobile)

    # get NH vectors from models (model_data[] -> vectors{resnum : vector})
    vector_data = []
    s2_pairs = {'N': 'H', 'CA': 'HA'}

    for model_num in calculate_on_models:
        model_data.atomgroup.setACSIndex(model_num)
        current_Resindex = 1
        has_first, has_second = False, False
        vectors = {}

        for atom in model_data.atomgroup:
            atom_res = atom.getResnum()

            if atom_res != current_Resindex:
                current_Resindex = atom_res
                has_first, has_second = False, False

            if atom_res == current_Resindex:
                if atom.getName() == S2_type:
                    has_second = True
                    N_coords = csx_obj.Vec_3D(atom.getCoords())

                elif atom.getName() == s2_pairs[S2_type]:
                    has_first = True
                    H_coords = csx_obj.Vec_3D(atom.getCoords())

                if has_first and has_second:
                    has_first, has_second = False, False
                    vectors[atom_res] = csx_obj.Vec_3D(
                        N_coords - H_coords
                    ).normalize()

        vector_data.append(vectors)

    S2_calced = {}

    # iterating over STR records
    for resnum in [int(s2rec.resnum) for s2rec in S2_records]:

        x2, y2, z2, xy, xz, yz = 0, 0, 0, 0, 0, 0

        # iterating over PDB models
        for m in vector_data:

            # coordinates in model at a given resnum
            x, y, z = m[resnum].v[0], m[resnum].v[1], m[resnum].v[2]

            x2 += x ** 2
            y2 += y ** 2
            z2 += z ** 2
            xy += x * y
            xz += x * z
            yz += y * z

        x2 /= len(vector_data)
        y2 /= len(vector_data)
        z2 /= len(vector_data)
        xy /= len(vector_data)
        xz /= len(vector_data)
        yz /= len(vector_data)

        # S2 calcuation
        s2 = 3 / 2.0 * (x2 ** 2 + y2 ** 2 + z2 ** 2 +
                        2 * xy ** 2 + 2 * xz ** 2 + 2 * yz ** 2) - 0.5

        S2_calced[resnum] = s2

    return S2_calced


def calcPeptideBonds():
    """Calculates backbone diherdral angles (OMEGA) CA-N-C'-CA"""
    # model_list = csx_obj.PDB_model.model_list
    dihedral_angles = {"<2":    0, "2-5": 0, "5-10": 0,
                       "10-20": 0, ">20": 0}

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
