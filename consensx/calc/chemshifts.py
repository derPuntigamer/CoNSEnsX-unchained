import pickle
import os
import subprocess

import consensx.graph as graph

from consensx.csx_libs import methods as csx_func
from consensx.csx_libs import objects as csx_obj


def call_shiftx_on(my_path, pdb_files):
    """Call ShiftX on PDB models. Each output is appended to 'out_name'"""
    for i, pdb_file in enumerate(pdb_files):
        pdb_file = my_path + pdb_file
        out_name = my_path + "/modell_" + str(i+1) + ".out"
        subprocess.call([csx_obj.ThirdParty.shiftx, '1', pdb_file, out_name])

    averageHA, averageH, averageN = {}, {}, {}
    averageCA, averageCB, averageC = {}, {}, {}
    modHA, modH, modN, modCA, modCB, modC = {}, {}, {}, {}, {}, {}
    model_data_list = []

    for some_file in csx_func.natural_sort(os.listdir(my_path)):
        if some_file.startswith("modell") and some_file.endswith(".out"):
            out_file = open(my_path + some_file)
            part = 0

            for line in out_file:
                if line.strip().startswith("NUM"):
                    part += 1

                    if part == 2:
                        break

                    continue

                if line.strip().startswith("---"):
                    continue

                if line.strip():
                    line_values = line.split()
                    try:
                        resnum = int(line_values[0])
                    except ValueError:
                        resnum = int(line_values[0][1:])
                    HA = float(line_values[2])
                    H = float(line_values[3])
                    N = float(line_values[4])
                    CA = float(line_values[5])
                    CB = float(line_values[6])
                    C = float(line_values[7])

                    modHA[resnum] = HA
                    modH[resnum] = H
                    modN[resnum] = N
                    modCA[resnum] = CA
                    modCB[resnum] = CB
                    modC[resnum] = C

                    if resnum in list(averageHA.keys()):
                        averageHA[resnum] += HA
                    else:
                        averageHA[resnum] = HA

                    if resnum in list(averageH.keys()):
                        averageH[resnum] += H
                    else:
                        averageH[resnum] = H

                    if resnum in list(averageN.keys()):
                        averageN[resnum] += N
                    else:
                        averageN[resnum] = N

                    if resnum in list(averageCA.keys()):
                        averageCA[resnum] += CA
                    else:
                        averageCA[resnum] = CA

                    if resnum in list(averageCB.keys()):
                        averageCB[resnum] += CB
                    else:
                        averageCB[resnum] = CB

                    if resnum in list(averageC.keys()):
                        averageC[resnum] += C
                    else:
                        averageC[resnum] = C

            out_file.close()
            model_data_list.append({"HA": modHA, "H": modH, "N":  modN,
                                    "CA": modCA, "CB": modCB, "C": modC})
            modHA, modH, modN, modCA, modCB, modC = {}, {}, {}, {}, {}, {}

    averages = [averageHA, averageH, averageN, averageCA, averageCB, averageC]

    for avg_dict in averages:
        for key in avg_dict:
            avg_dict[key] /= len(pdb_files)

    return {"HA": averageHA, "H":  averageH, "N":  averageN,  "CA": averageCA,
            "CB": averageCB, "C": averageC}, model_data_list


def chemshifts(my_CSV_buffer, ChemShift_lists, pdb_models, my_path):
    """Back calculate chemical shifts from given chemical shift list and PDB
       models"""
    cs_data = []
    cs_calced, model_data = call_shiftx_on(my_path, pdb_models)

    csx_obj.ChemShift_modell_data.type_dict = model_data

    cs_model_data_path = my_path + "/ChemShift_model_data.pickle"
    pickle.dump(model_data, open(cs_model_data_path, 'wb'))

    for n, cs_list in enumerate(ChemShift_lists):
        for CS_type in sorted(list(cs_list.keys())):
            model_corrs = []

            for model in model_data:
                inner_exp = {}

                for record in cs_list[CS_type]:
                    inner_exp[record.resnum] = model[CS_type][record.resnum]

                model_corrs.append(csx_func.calcCorrel(inner_exp,
                                                       cs_list[CS_type]))

            avg_model_corr = sum(model_corrs) / len(model_corrs)

            exp_dict = {}

            for record in cs_list[CS_type]:
                exp_dict[record.resnum] = cs_calced[CS_type][record.resnum]

            correl = csx_func.calcCorrel(exp_dict, cs_list[CS_type])
            q_value = csx_func.calcQValue(exp_dict, cs_list[CS_type])
            rmsd = csx_func.calcRMSD(exp_dict, cs_list[CS_type])

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
                "experimental": cs_list[CS_type]
            })

            print("CHEM SHIFT", CS_type)
            print("Correl: ", correl)
            print("Q-val:  ", q_value)
            print("RMSD:   ", rmsd)
            print()

            graph_name = str(n + 1) + "_CS_" + CS_type + ".svg"
            graph.values_graph(my_path, exp_dict, cs_list[CS_type], graph_name)

            corr_graph_name = str(n + 1) + "_CS_corr_" + CS_type + ".svg"
            graph.correl_graph(
                my_path, exp_dict, cs_list[CS_type], corr_graph_name
            )

            mod_corr_graph_name = "CS_mod_corr_" + CS_type + ".svg"
            graph.mod_correl_graph(
                my_path, correl, avg_model_corr,
                model_corrs, mod_corr_graph_name
            )

            my_id = my_path.split('/')[-2] + '/'

            cs_data.append({
                "CS_type": CS_type,
                "CS_model_n": len(cs_list[CS_type]),
                "correlation": '{0:.3f}'.format(correl),
                "q_value": '{0:.3f}'.format(q_value),
                "rmsd": '{0:.3f}'.format(rmsd),
                "corr_graph_name": my_id + corr_graph_name,
                "graph_name": my_id + graph_name,
                "mod_corr_graph_name": my_id + mod_corr_graph_name,
                "input_id": "CS_" + CS_type
            })

    return cs_data