import math

from consensx.csx_libs import methods as csx_func
from consensx.csx_libs import objects as csx_obj


def noe_violations(PDB_file, saveShifts, my_path, r3_averaging):
    """Back calculate NOE distance violations from given RDC lists and PDB
    models"""
    # parse data to restraint objects returned from pypy process
    csx_id = 1
    prev_id = 1
    for data in saveShifts:
        csx_obj.Restraint_Record(csx_id, data[0], data[1], data[2], data[3],
                                 data[4], data[5], data[6], data[7])

        if prev_id != data[0]:
            prev_id = data[0]
            csx_id += 1

    # fetch all restraint from class
    restraints = csx_obj.Restraint_Record.getNOERestraints()

    PDB_coords = csx_func.pdb2coords(PDB_file)
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

        print("DISTS", dist_str)

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

        for model in list(PDB_coords.keys()):
            avg += math.pow(avg_distances[model][curr_id], -6)

        avg /= len(list(PDB_coords.keys()))
        measured_avg[curr_id] = math.pow(avg, -1.0/6)

    # at this point measured_avg[curr_id] is a simple dictonary containing the
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
    csx_func.makeNOEHist(my_path, violations)

    return viol_count
