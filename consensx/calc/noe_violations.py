import math
import matplotlib.pyplot as plt
import numpy as np

from consensx.csx_libs import objects as csx_obj


def make_noe_hist(my_path, violations):
    plt.figure(figsize=(6, 5), dpi=80)
    n_groups = len(violations)

    means_men = [
        violations['0-0.5'], violations['0.5-1'], violations['1-1.5'],
        violations['1.5-2'], violations['2-2.5'], violations['2.5-3'],
        violations['3<']
    ]

    ticks = ['0-0.5', '0.5-1', '1-1.5', '1.5-2', '2-2.5', '2.5-3', '3<']
    index = np.arange(n_groups)
    bar_width = 0.7
    plt.bar(index, means_men, bar_width, alpha=.7, color='b')

    plt.xlabel("Violation (Å)")
    plt.ylabel("# of NOE distance violations")
    plt.title("NOE distance violations")
    plt.xticks(index + bar_width / 2, ticks)
    ax = plt.axes()
    ax.yaxis.grid()

    plt.tight_layout()
    plt.savefig(my_path + "/NOE_hist.svg", format="svg")
    plt.close()


def pdb2coords(model_data):
    """Loads PDB coordinates into a dictonary, per model"""
    prev_resnum = -1
    PDB_coords = {}

    for i in range(model_data.coordsets):
        model_data.atomgroup.setACSIndex(i)

        PDB_coords[i] = {}

        for atom in model_data.atomgroup:
            resnum = int(atom.getResnum())
            name = str(atom.getName())

            if resnum == prev_resnum:
                PDB_coords[i][resnum][name] = csx_obj.Vec_3D(atom.getCoords())

            else:
                PDB_coords[i][resnum] = {}
                PDB_coords[i][resnum][name] = csx_obj.Vec_3D(atom.getCoords())
                prev_resnum = resnum

    return PDB_coords


def get_noe_from_file(NOE_file):
    NOE_key = "_Gen_dist_constraint_list.Constraint_type"
    NOE_value = "NOE"
    in_frame = False
    in_loop = False
    loops = []
    NOE_data = []

    for line in open(NOE_file, encoding="utf-8"):
        tok = line.strip().split()

        if not tok:
            continue

        if tok[0] == NOE_key and tok[1] == NOE_value:
            in_frame = True
            continue

        if in_frame and tok[0] == "save_":
            # break
            in_frame = False

        if tok[0] == "loop_":
            loop_keys = []
            loop_data = []
            in_loop = True
            continue

        if in_loop is True and in_frame is True and tok[0] == "stop_":
            in_loop = False
            loops.append([loop_keys, loop_data])

        if in_loop:
            if tok[0].startswith('_'):
                loop_keys.append(tok[0])
            else:
                loop_data.append(tok)

    for loop in loops:
        try:
            ind_ID = loop[0].index("_Gen_dist_constraint.ID")
            ind_seg1 = loop[0].index("_Gen_dist_constraint.PDB_residue_no_1")
            ind_seg2 = loop[0].index("_Gen_dist_constraint.PDB_residue_no_2")
            ind_comp1 = loop[0].index("_Gen_dist_constraint.Comp_ID_1")
            ind_comp2 = loop[0].index("_Gen_dist_constraint.Comp_ID_2")
            ind_atom1 = loop[0].index("_Gen_dist_constraint.Atom_ID_1")
            ind_atom2 = loop[0].index("_Gen_dist_constraint.Atom_ID_2")
            ind_bnd = loop[0].index(
                "_Gen_dist_constraint.Distance_upper_bound_val")
        except ValueError:
            continue

        for data in loop[1]:
            NOE_data.append([
                data[ind_ID],
                data[ind_seg1], data[ind_seg2],
                data[ind_comp1], data[ind_comp2],
                data[ind_atom1], data[ind_atom2],
                data[ind_bnd]
            ])

    return NOE_data


def noe_violations(model_data, my_path, db_entry, bme_weights=None):
    """Back calculate NOE distance violations from given RDC lists and PDB
    models"""

    r3_averaging = db_entry.r3average

    # empty class variables
    csx_obj.Restraint_Record.all_restraints = []
    csx_obj.Restraint_Record.resolved_restraints = []

    noe_name = db_entry.NOE_file
    noe_file_path = my_path + noe_name
    save_shifts = get_noe_from_file(noe_file_path)
    noe_n = save_shifts[-1][0] + " distance restraints found"

    csx_id = 1
    prev_id = 1
    for data in save_shifts:
        csx_obj.Restraint_Record(
            csx_id, data[0], data[1], data[2], data[3],
            data[4], data[5], data[6], data[7]
        )

        if prev_id != data[0]:
            prev_id = data[0]
            csx_id += 1

    # fetch all restraint from class
    restraints = csx_obj.Restraint_Record.getNOERestraints(model_data)

    PDB_coords = pdb2coords(model_data)
    prev_id = -1
    avg_distances = {}
    all_distances = {}
    measured_avg = {}
    str_distaces = {}

    for model in list(PDB_coords.keys()):
        avg_distances[model] = {}
        all_distances[model] = {}

        for restraint_num, restraint in enumerate(restraints):
            rest_id = int(restraint.csx_id)
            resnum1 = restraint.seq_ID1
            atom1 = restraint.atom_ID1
            resnum2 = restraint.seq_ID2
            atom2 = restraint.atom_ID2

            atom_coord1 = PDB_coords[model][resnum1][atom1]
            atom_coord2 = PDB_coords[model][resnum2][atom2]

            distance = (atom_coord1 - atom_coord2).magnitude()

            all_distances[model][restraint_num] = distance

            if prev_id == rest_id:
                avg_distances[model][rest_id].append(distance)

            else:
                prev_id = rest_id
                avg_distances[model][rest_id] = []
                str_distaces[rest_id] = restraint.dist_max

                avg_distances[model][rest_id].append(distance)

    for restraint_num, restraint in enumerate(restraints):
        rest_id = int(restraint.csx_id)
        resnum1 = restraint.seq_ID1
        atom1 = restraint.atom_ID1
        resnum2 = restraint.seq_ID2
        atom2 = restraint.atom_ID2

        dist_str = "> {} {} {} {} {}  |   ".format(
            rest_id, resnum1, atom1, resnum2, atom2
        )

        for model in list(PDB_coords.keys()):
            dist_str += "{0:.2f}  ".format(all_distances[model][restraint_num])

        # print("DISTS", dist_str)

    # at this point avg_distances[model][curr_id] contains distances for one
    # model and one restraint GROUP identified with "csx_id" number

    prev_id = -1

    for model in list(PDB_coords.keys()):
        for restraint in restraints:
            curr_id = int(restraint.csx_id)

            if prev_id == curr_id:
                continue
            else:
                prev_id = curr_id

            avg = 0.0

            for distance in avg_distances[model][curr_id]:
                if r3_averaging:
                    avg += math.pow(float(distance), -3)
                else:
                    avg += math.pow(float(distance), -6)

            avg /= len(avg_distances[model][curr_id])
            if r3_averaging:
                avg_distances[model][curr_id] = math.pow(avg, -1.0/3)
            else:
                avg_distances[model][curr_id] = math.pow(avg, -1.0/6)

    # at this point avg_distances[model][curr_id] contain a single (r-6)
    # averaged distance for one model and one restraint GROUP identified with
    # "csx_id" number. Averaging is done on "in GROUP" distances

    for restraint in restraints:
        curr_id = int(restraint.curr_distID)
        avg = 0.0

        if bme_weights:
            for model in list(PDB_coords.keys()):
                avg += math.pow(
                    avg_distances[model][curr_id], -6
                ) * bme_weights[model]

            avg /= sum(bme_weights)
        else:
            for model in list(PDB_coords.keys()):
                avg += math.pow(avg_distances[model][curr_id], -6)

            avg /= len(list(PDB_coords.keys()))

        measured_avg[curr_id] = math.pow(avg, -1.0/6)

    bme_exp_filename = "noe_exp.dat"
    bme_calc_filename = "noe_calc.dat"

    with open(my_path + bme_exp_filename, "w") as exp_dat_file:
        exp_dat_file.write("# DATA=NOE PRIOR=GAUSS POWER=6\n")
        prev_id = -1

        for restraint in restraints:
            if prev_id == restraint.csx_id:
                continue

            prev_id = restraint.csx_id

            exp_dat_file.write(
                str(restraint.csx_id) + "\t" +
                str(restraint.dist_max) + "\t0.1\n"
            )

    with open(my_path + bme_calc_filename, "w") as calc_dat_file:
        for model in list(PDB_coords.keys()):

            for i in avg_distances[model]:
                calc_dat_file.write(
                    str(avg_distances[model][i]) + " "
                )

        calc_dat_file.write("\n")

    # at this point measured_avg[curr_id] is a simple dictionary containing the
    # model averaged distances for the given "csx_id" number

    avg_dist_keys = list(measured_avg.keys())
    avg_dist_keys.sort()
    violations = {"0-0.5": 0, "0.5-1": 0, "1-1.5": 0,
                  "1.5-2": 0, "2-2.5": 0, "2.5-3": 0, "3<": 0}
    viol_count = 0

    for key in avg_dist_keys:
        if measured_avg[key] > str_distaces[key]:
            viol_count += 1
            diff = measured_avg[key] - str_distaces[key]

            if diff <= 0.5:
                violations["0-0.5"] += 1
            elif 0.5 < diff <= 1:
                violations["0.5-1"] += 1
            elif 1 < diff <= 1.5:
                violations["1-1.5"] += 1
            elif 1.5 < diff <= 2:
                violations["1.5-2"] += 1
            elif 2 < diff <= 2.5:
                violations["2-2.5"] += 1
            elif 2.5 < diff <= 3:
                violations["2.5-3"] += 1
            else:
                violations["3<"] += 1

    print("Total # of violations:", viol_count)
    make_noe_hist(my_path, violations)

    return noe_n, viol_count
