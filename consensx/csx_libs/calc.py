import math
import subprocess
import os
import prody
import pickle

import matplotlib.pyplot as plt
import numpy as np

from . import methods as csx_func
from . import objects as csx_obj

plt.switch_backend('Agg')


def calcS2(my_CSV_buffer, S2_dict, my_path, calculate_on_models=None, fit=None,
           fit_range=None):
    """Back calculate order paramteres from given S2 dict and PDB models"""
    S2_data = []
    model_data = csx_obj.PDB_model.model_data

    if not calculate_on_models:
        calculate_on_models = range(model_data.model_count)

    for S2_type in sorted(list(S2_dict.keys())):
        S2_calced = csx_func.calcS2(model_data, calculate_on_models,
                                    S2_dict[S2_type],
                                    S2_type, fit, fit_range)

        correl = csx_func.calcCorrel(S2_calced, S2_dict[S2_type])
        q_value = csx_func.calcQValue(S2_calced, S2_dict[S2_type])
        rmsd = csx_func.calcRMSD(S2_calced, S2_dict[S2_type])

        # TODO DB upload!
        corr_key = "S2_" + S2_type + "_corr"
        qval_key = "S2_" + S2_type + "_qval"
        rmsd_key = "S2_" + S2_type + "_rmsd"

        csx_obj.CalcPickle.data.update({
            corr_key: "{0}".format('{0:.3f}'.format(correl)),
            qval_key: "{0}".format('{0:.3f}'.format(q_value)),
            rmsd_key: "{0}".format('{0:.3f}'.format(rmsd))
        })

        my_CSV_buffer.csv_data.append({
            "name": "S2 (" + S2_type + ")",
            "calced": S2_calced,
            "experimental": S2_dict[S2_type]
        })

        print(S2_type + " Order Parameters")
        print("Correl: ", correl)
        print("Q-val:  ", q_value)
        print("RMSD:   ", rmsd)
        print()

        graph_name = "S2_" + S2_type + ".svg"
        csx_func.makeGraph(my_path, S2_calced, S2_dict[S2_type], graph_name)

        corr_graph_name = "S2_corr_" + S2_type + ".svg"
        csx_func.makeCorrelGraph(my_path, S2_calced, S2_dict[S2_type],
                                 corr_graph_name)

        my_id = my_path.split('/')[-2] + '/'

        S2_data.append({
            "S2_type": S2_type,
            "S2_model_n": len(S2_dict[S2_type]),
            "correlation": '{0:.3f}'.format(correl),
            "q_value": '{0:.3f}'.format(q_value),
            "rmsd": '{0:.3f}'.format(rmsd),
            "corr_graph_name": my_id + corr_graph_name,
            "graph_name": my_id + graph_name,
            "input_id": "S2_" + S2_type
        })

    return S2_data


def calcS2_sidechain(my_CSV_buffer, S2_sidechain, my_path, fit=None):
    """Back calculate order paramteres from given S2 dict and PDB models"""

    sc_LOT = {
        'VAL': {'CG1': 'CB', 'CG2': 'CB'},
        'ILE': {'CG2': 'CB', 'CD': 'CG1', 'CD1': 'CG1'},
        'THR': {'CG2': 'CB'},
        'LEU': {'CD1': 'CG', 'CD2': 'CG'},
        'ALA': {'CB': 'CA'},
        'MET': {'CE': 'SD'}
    }

    model_data = csx_obj.PDB_model.model_data

    if fit and not csx_obj.PDB_model.is_fitted:
        print("Start FITTING")

        model_data.atomgroup.setACSIndex(0)
        prody.alignCoordsets(model_data.atomgroup.calpha)

    for record in S2_sidechain:
        vectors = []
        resnum = record.resnum
        my_type = record.type
        my_res = model_data.atomgroup[('A', resnum)].getResname()

        # find pair for measured aa
        pair = sc_LOT[my_res][my_type]

        for model_num in range(model_data.model_count):
            model_data.atomgroup.setACSIndex(model_num)

            try:
                sel = "resnum {} name {}".format(resnum, my_type)
                coords = model_data.atomgroup.select(sel).getCoords()[0]
                sel = "resnum {} name {}".format(resnum, pair)
                pair_coords = model_data.atomgroup.select(sel).getCoords()[0]
            except AttributeError:
                return {
                    "error": "Sidechain order parameter atom name not found\
                    in PDB. Please check your atom naming."
                }

            vectors.append(csx_obj.Vec_3D(coords - pair_coords).normalize())

        x2, y2, z2, xy, xz, yz = 0, 0, 0, 0, 0, 0

        for vector in vectors:
            x, y, z = vector.v[0], vector.v[1], vector.v[2]

            x2 += x ** 2
            y2 += y ** 2
            z2 += z ** 2
            xy += x * y
            xz += x * z
            yz += y * z

        x2 /= len(vectors)
        y2 /= len(vectors)
        z2 /= len(vectors)
        xy /= len(vectors)
        xz /= len(vectors)
        yz /= len(vectors)

        s2 = 3 / 2.0 * (x2 ** 2 + y2 ** 2 + z2 ** 2 +
                        2 * xy ** 2 + 2 * xz ** 2 + 2 * yz ** 2) - 0.5

        record.calced = s2

    # prepare data for CSV object inicialization
    sidechain_exp1, sidechain_exp2 = [], []
    sidechain_calc1, sidechain_calc2 = {}, {}

    prev_resnum = -100000

    for record in S2_sidechain:
        if record.resnum != prev_resnum:
            sidechain_exp1.append(record)
            sidechain_calc1[record.resnum] = record.calced
            prev_resnum = record.resnum
        else:
            sidechain_exp2.append(record)
            sidechain_calc2[record.resnum] = record.calced

    my_CSV_buffer.csv_data.append({
        "name": "S2_meth",
        "calced": sidechain_calc1,
        "experimental": sidechain_exp1
    })

    if sidechain_exp2:
        my_CSV_buffer.csv_data.append({
            "name": "S2_meth (cont)",
            "calced": sidechain_calc2,
            "experimental": sidechain_exp2
        })

    # correlation calculation
    M = [0.0, 0.0, 0.0]
    D = [0.0, 0.0]

    for record in S2_sidechain:
        exp = record.value
        calc = record.calced

        M[0] += calc
        M[1] += exp
        M[2] += calc * exp

    M[0] /= len(S2_sidechain)
    M[1] /= len(S2_sidechain)
    M[2] /= len(S2_sidechain)

    for record in S2_sidechain:
        exp = record.value
        calc = record.calced

        D[0] += (calc - M[0]) ** 2
        D[1] += (exp - M[1]) ** 2

    D[0] /= len(S2_sidechain)
    D[0] = math.sqrt(D[0])
    D[1] /= len(S2_sidechain)
    D[1] = math.sqrt(D[1])

    correl = (M[2] - (M[0] * M[1])) / (D[0] * D[1])
    print("Corr: ", correl)

    # Q-value calculation
    D2, E2 = 0, 0

    for record in S2_sidechain:
        exp = record.value
        calc = record.calced

        D2 += (calc - exp) ** 2
        E2 += exp ** 2

    Q = 100 * math.sqrt(D2) / math.sqrt(E2)
    q_value = round(Q, 6)
    print("Q-value: ", q_value)

    # RMDS calculation
    D2 = 0

    for record in S2_sidechain:
        exp = record.value
        calc = record.calced

        D2 += (calc - exp) ** 2

    rmsd = math.sqrt(D2 / len(S2_sidechain))
    print("RMSD: ", round(rmsd, 6))

    exp_values = []
    calced_values = []

    for record in S2_sidechain:
        exp_values.append(record.value)
        calced_values.append(record.calced)

    min_calc = min(calced_values)
    max_calc = max(calced_values)

    min_exp = min(exp_values)
    max_exp = max(exp_values)
    miny = min(min_calc, min_exp)             # get minimum value
    maxy = max(max_calc, max_exp)             # get maximum value

    diag = []

    margin = int(abs(miny - maxy) * 0.05)

    if abs(miny - maxy) < 10:
        margin = 0.3
    elif abs(miny - maxy) < 2:
        margin = 0.01
    elif abs(miny - maxy) < 1:
        margin = 0

    maxy += margin
    miny -= margin

    for i in np.arange(miny, maxy * 1.42, 0.1):  # draw graph diagonal
        diag.append(i)

    plt.figure(figsize=(6, 5), dpi=80)
    plt.plot(diag, diag, linewidth=2.0, color='red', alpha=.7)
    plt.plot(exp_values, calced_values, 'bo')
    plt.axis([miny, maxy, miny, maxy])
    plt.xlabel('experimental')
    plt.ylabel('calculated')
    plt.tight_layout(pad=1.08)
    plt.savefig(my_path + "/" + "S2_sc_corr.svg", format="svg")
    plt.close()

    xs = []
    prev_resnum2 = -1

    for record in S2_sidechain:
        if record.resnum != prev_resnum2:
            xs.append(record.resnum)
            prev_resnum2 = record.resnum
        else:
            xs.append(record.resnum + 0.3)

    print("XS AXIS", xs)
    print("len SC", len(S2_sidechain))
    print("len xs", len(xs))

    plt.figure(figsize=(10, 5), dpi=80)
    plt.plot(xs, exp_values,
             linewidth=2.0, color='red', marker='o', label='exp', alpha=.7)
    plt.plot(xs, calced_values,
             linewidth=2.0, color='blue', marker='o', label='calc', alpha=.7)
    plt.legend(loc='lower left')
    plt.xlabel('residue number')
    plt.ylabel('value')
    ax = plt.axes()
    ax.yaxis.grid()
    plt.tight_layout(pad=1.08)
    plt.savefig(my_path + "/" + "S2_sc_graph.svg", format="svg")
    plt.close()

    my_id = my_path.split('/')[-2] + '/'

    print("CORR GRAPH", my_id + "S2_sc_corr.svg")
    print("GRAPH", my_id + "S2_sc_graph.svg")

    S2_sc_data = {
        "S2_model_n": len(S2_sidechain),
        "correlation": '{0:.3f}'.format(correl),
        "q_value": '{0:.3f}'.format(q_value),
        "rmsd": '{0:.3f}'.format(rmsd),
        "corr_graph_name": my_id + "S2_sc_corr.svg",
        "graph_name": my_id + "S2_sc_graph.svg"
    }

    return S2_sc_data


def calcJCouplings(my_CSV_buffer, param_set, Jcoup_dict, my_PDB, my_path):
    """Back calculate skalar coupling from given RDC lists and PDB models"""
    Jcuop_data = []
    type_dict = {}
    model_list = csx_obj.PDB_model.model_data
    dihed_lists = csx_func.calcDihedAngles(model_list)

    for Jcoup_type in sorted(list(Jcoup_dict.keys())):

        JCoup_calced, model_data = csx_func.calcJCoup(param_set,
                                                      dihed_lists,
                                                      Jcoup_dict[Jcoup_type],
                                                      Jcoup_type)

        type_dict[Jcoup_type] = model_data

        model_corrs = []

        for model in model_data:
            model_corrs.append(csx_func.calcCorrel(model,
                                                   Jcoup_dict[Jcoup_type]))

        avg_model_corr = sum(model_corrs) / len(model_corrs)

        correl = csx_func.calcCorrel(JCoup_calced, Jcoup_dict[Jcoup_type])
        q_value = csx_func.calcQValue(JCoup_calced, Jcoup_dict[Jcoup_type])
        rmsd = csx_func.calcRMSD(JCoup_calced, Jcoup_dict[Jcoup_type])

        # TODO
        corr_key = "JCoup_" + Jcoup_type + "_corr"
        qval_key = "JCoup_" + Jcoup_type + "_qval"
        rmsd_key = "JCoup_" + Jcoup_type + "_rmsd"

        csx_obj.CalcPickle.data.update({
            corr_key: "{0}".format('{0:.3f}'.format(correl)),
            qval_key: "{0}".format('{0:.3f}'.format(q_value)),
            rmsd_key: "{0}".format('{0:.3f}'.format(rmsd))
        })

        my_CSV_buffer.csv_data.append({
            "name": "J-couplings (" + Jcoup_type + ")",
            "calced": JCoup_calced,
            "experimental": Jcoup_dict[Jcoup_type]
        })

        print("J-couplings (" + Jcoup_type + ")")
        print("Correl: ", correl)
        print("Q-val:  ", q_value)
        print("RMSD:   ", rmsd)
        print()

        graph_name = "JCoup_" + Jcoup_type + ".svg"
        csx_func.makeGraph(my_path, JCoup_calced, Jcoup_dict[Jcoup_type],
                           graph_name)

        corr_graph_name = "JCoup_corr_" + Jcoup_type + ".svg"
        csx_func.makeCorrelGraph(my_path, JCoup_calced, Jcoup_dict[Jcoup_type],
                                 corr_graph_name)

        mod_corr_graph_name = "JCoup_mod_corr_" + Jcoup_type + ".svg"
        csx_func.modCorrelGraph(my_path, correl, avg_model_corr, model_corrs,
                                mod_corr_graph_name)

        my_id = my_path.split('/')[-2] + '/'

        Jcuop_data.append({
            "Jcoup_type": Jcoup_type,
            "Jcoop_model_n": len(Jcoup_dict[Jcoup_type]),
            "correlation": '{0:.3f}'.format(correl),
            "q_value": '{0:.3f}'.format(q_value),
            "rmsd": '{0:.3f}'.format(rmsd),
            "corr_graph_name": my_id + corr_graph_name,
            "graph_name": my_id + graph_name,
            "mod_corr_graph_name": my_id + mod_corr_graph_name,
            "input_id": "JCoup_" + Jcoup_type
        })

    Jcoup_model_data_path = my_path + "/Jcoup_model.pickle"
    pickle.dump(type_dict, open(Jcoup_model_data_path, 'wb'))
    return Jcuop_data


def calcChemShifts(my_CSV_buffer, ChemShift_lists, pdb_models, my_path):
    """Back calculate chemical shifts from given chemical shift list and PDB
       models"""
    CS_data = []
    CS_calced, model_data = csx_func.callShiftxOn(my_path, pdb_models)

    csx_obj.ChemShift_modell_data.type_dict = model_data

    CS_model_data_path = my_path + "/ChemShift_model_data.pickle"
    pickle.dump(model_data, open(CS_model_data_path, 'wb'))

    for n, CS_list in enumerate(ChemShift_lists):
        for CS_type in sorted(list(CS_list.keys())):
            model_corrs = []

            for model in model_data:
                inner_exp = {}

                for record in CS_list[CS_type]:
                    inner_exp[record.resnum] = model[CS_type][record.resnum]

                model_corrs.append(csx_func.calcCorrel(inner_exp,
                                                       CS_list[CS_type]))

            avg_model_corr = sum(model_corrs) / len(model_corrs)

            exp_dict = {}

            for record in CS_list[CS_type]:
                exp_dict[record.resnum] = CS_calced[CS_type][record.resnum]

            correl = csx_func.calcCorrel(exp_dict, CS_list[CS_type])
            q_value = csx_func.calcQValue(exp_dict, CS_list[CS_type])
            rmsd = csx_func.calcRMSD(exp_dict, CS_list[CS_type])

            # TODO
            corr_key = "CS_" + CS_type + "_corr"
            qval_key = "CS_" + CS_type + "_qval"
            rmsd_key = "CS_" + CS_type + "_rmsd"

            csx_obj.CalcPickle.data.update({
                corr_key: "{0}".format('{0:.3f}'.format(correl)),
                qval_key: "{0}".format('{0:.3f}'.format(q_value)),
                rmsd_key: "{0}".format('{0:.3f}'.format(rmsd))
            })

            my_CSV_buffer.csv_data.append({
                "name": "ChemShifts (" + CS_type + ")",
                "calced": exp_dict,
                "experimental": CS_list[CS_type]
            })

            print("CHEM SHIFT", CS_type)
            print("Correl: ", correl)
            print("Q-val:  ", q_value)
            print("RMSD:   ", rmsd)
            print()

            graph_name = str(n + 1) + "_CS_" + CS_type + ".svg"
            csx_func.makeGraph(my_path, exp_dict, CS_list[CS_type],
                               graph_name)

            corr_graph_name = str(n + 1) + "_CS_corr_" + CS_type + ".svg"
            csx_func.makeCorrelGraph(my_path, exp_dict, CS_list[CS_type],
                                     corr_graph_name)

            mod_corr_graph_name = "CS_mod_corr_" + CS_type + ".svg"
            csx_func.modCorrelGraph(
                my_path, correl, avg_model_corr,
                model_corrs, mod_corr_graph_name
            )

            my_id = my_path.split('/')[-2] + '/'

            CS_data.append({
                "CS_type": CS_type,
                "CS_model_n": len(CS_list[CS_type]),
                "correlation": '{0:.3f}'.format(correl),
                "q_value": '{0:.3f}'.format(q_value),
                "rmsd": '{0:.3f}'.format(rmsd),
                "corr_graph_name": my_id + corr_graph_name,
                "graph_name": my_id + graph_name,
                "mod_corr_graph_name": my_id + mod_corr_graph_name,
                "input_id": "CS_" + CS_type
            })

    return CS_data


def calcNOEviolations(PDB_file, saveShifts, my_path, r3_averaging):
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


def calcNMR_Pride(pdb_models, my_path):
    """Calculate NMR-PRIDE score on given PDB models"""

    pwd = os.getcwd()
    os.chdir(my_path)

    # write model list text file
    model_list = open("model_list.txt", 'w')

    for model in pdb_models:
        model_list.write(model + "\n")

    model_list.write("END\n")
    model_list.close()

    # write distance dict to text file
    restraints = csx_obj.Restraint_Record.getPRIDE_restraints()
    pride_input = open("pride_input.txt", 'w')

    pride_input.write("HEADER\n")

    prime_distances = list(restraints.keys())
    prime_distances.sort()

    for distance in prime_distances:
        pride_input.write(str(distance) + ' ' +
                          str(restraints[distance]) + '\n')

    pride_input.write("END\n")
    pride_input.close()

    # create binary database for PRIDE-NMR
    DEVNULL = open(os.devnull, 'w')
    hhdb_log = open("hhdb.log", 'w')
    model_list = open("model_list.txt", 'r')
    subprocess.call(
        [
            csx_obj.ThirdParty.prideDB,
            "-D", "HHDB",  # model list
        ],
        stdin=model_list,
        stdout=DEVNULL,
        stderr=hhdb_log
    )

    hhdb_log.close()
    DEVNULL.close()
    model_list.close()

    # run PRIDE-NMR
    DEVNULL = open(os.devnull, 'w')
    pride_input = open("pride_input.txt", 'r')
    pride_output = open("pride_output.txt", 'w')
    subprocess.call(
        [
            csx_obj.ThirdParty.prideNMR,
            "-D", "HHDB",
            "-d", str(56),
            "-b", str(len(pdb_models)),
            "-m", str(3)
        ],
        stdin=pride_input,
        stdout=pride_output,
        stderr=DEVNULL
    )

    pride_input.close()
    pride_output.close()
    DEVNULL.close()

    pride_scores = {}
    pride_output = open("pride_output.txt", 'r')
    for line in pride_output:
        if line.startswith("PRIDENMR:"):
            model_num = int(line.split()[-1])
            model_score = float(line.split()[1])
            pride_scores[model_num] = model_score

    scores = list(pride_scores.values())

    avg = sum(scores) * 1.0 / len(scores)
    variance = [(x - avg) ** 2 for x in scores]
    standard_deviation = math.sqrt(sum(variance) * 1.0 / len(variance))

    PRIDE_data = []

    print("PRIDE-NMR calculation")
    print("MAX: ",    max(pride_scores, key=pride_scores.get))
    PRIDE_data.append(max(pride_scores, key=pride_scores.get))
    print("MIN: ",    min(pride_scores, key=pride_scores.get))
    PRIDE_data.append(min(pride_scores, key=pride_scores.get))
    print("AVG: ", avg)
    PRIDE_data.append(avg)
    print("DEV: ", standard_deviation, "\n")
    PRIDE_data.append(standard_deviation)

    os.chdir(pwd)

    csx_func.makeNMRPrideGraph(my_path, scores, avg)

    return PRIDE_data
