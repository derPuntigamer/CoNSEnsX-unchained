import pickle
import os
import subprocess

import consensx.graph as graph

from consensx.csx_libs import methods as csx_func
from consensx.csx_libs import objects as csx_obj

chemshift_types = ["HA", "CA", "CB", "N", "H", "C"]


def call_shiftx_on(my_path, pdb_files, bme_weights=None):
    """Call ShiftX on PDB models. Each output is appended to 'out_name'"""
    for i, pdb_file in enumerate(pdb_files):
        pdb_file = my_path + pdb_file
        out_name = my_path + "/modell_" + str(i+1) + ".out"
        subprocess.call([csx_obj.ThirdParty.shiftx, '1', pdb_file, out_name])

    shiftx_output_files = []
    averageHA, averageH, averageN = {}, {}, {}
    averageCA, averageCB, averageC = {}, {}, {}
    modHA, modH, modN, modCA, modCB, modC = {}, {}, {}, {}, {}, {}
    model_data_list = []
    model_weight = 1
    model_weight_sum = len(pdb_files)

    if bme_weights:
        model_weight_sum = sum(bme_weights)

    for some_file in csx_func.natural_sort(os.listdir(my_path)):
        if some_file.startswith("modell") and some_file.endswith(".out"):
            shiftx_output_files.append(some_file)

    for i, shiftx_output in enumerate(shiftx_output_files):
        if bme_weights:
            model_weight = bme_weights[i]

        out_file = open(my_path + shiftx_output)
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

                modHA[resnum] = HA * model_weight
                modH[resnum] = H * model_weight
                modN[resnum] = N * model_weight
                modCA[resnum] = CA * model_weight
                modCB[resnum] = CB * model_weight
                modC[resnum] = C * model_weight

                if resnum in list(averageHA.keys()):
                    averageHA[resnum] += HA * model_weight
                else:
                    averageHA[resnum] = HA * model_weight

                if resnum in list(averageH.keys()):
                    averageH[resnum] += H * model_weight
                else:
                    averageH[resnum] = H * model_weight

                if resnum in list(averageN.keys()):
                    averageN[resnum] += N * model_weight
                else:
                    averageN[resnum] = N * model_weight

                if resnum in list(averageCA.keys()):
                    averageCA[resnum] += CA * model_weight
                else:
                    averageCA[resnum] = CA * model_weight

                if resnum in list(averageCB.keys()):
                    averageCB[resnum] += CB * model_weight
                else:
                    averageCB[resnum] = CB * model_weight

                if resnum in list(averageC.keys()):
                    averageC[resnum] += C * model_weight
                else:
                    averageC[resnum] = C * model_weight

        out_file.close()

        model_data_list.append({
            "HA": modHA, "H": modH, "N":  modN,
            "CA": modCA, "CB": modCB, "C": modC
        })

        modHA, modH, modN, modCA, modCB, modC = {}, {}, {}, {}, {}, {}

    averages = [averageHA, averageH, averageN, averageCA, averageCB, averageC]

    for avg_dict in averages:
        for key in avg_dict:
            avg_dict[key] /= model_weight_sum

    return {"HA": averageHA, "H":  averageH, "N":  averageN,  "CA": averageCA,
            "CB": averageCB, "C": averageC}, model_data_list


def chemshifts(
        my_CSV_buffer, ChemShift_lists, pdb_models, my_path, bme_weights
        ):
    """Back calculate chemical shifts from given chemical shift list and PDB
       models"""
    cs_data = []
    cs_calced, model_data = call_shiftx_on(my_path, pdb_models, bme_weights)

    csx_obj.ChemShift_modell_data.type_dict = model_data

    cs_model_data_path = my_path + "/ChemShift_model_data.pickle"
    pickle.dump(model_data, open(cs_model_data_path, 'wb'))

    for n, cs_list in enumerate(ChemShift_lists):
        bme_exp_filename = "chemshift_" + str(n) + "_exp.dat"
        bme_calc_filename = "chemshift_" + str(n) + "_calc.dat"

        with open(my_path + bme_exp_filename, "w") as exp_dat_file:
            exp_dat_file.write("# DATA=CS PRIOR=GAUSS\n")

            for atom_type in chemshift_types:
                for i in cs_list[atom_type]:
                    exp_dat_file.write(
                        str(i.resnum) + "_" +
                        str(i.atom_name) + "\t" +
                        str(i.value) + "\t0.1\n"
                    )

        with open(my_path + bme_calc_filename, "w") as calc_dat_file:
            for index, model in enumerate(model_data):
                calc_dat_file.write(str(index) + " ")

                for atom_type in chemshift_types:
                    for i in cs_list[atom_type]:
                        calc_dat_file.write(
                            str(model[atom_type][i.resnum]) + " "
                        )

                calc_dat_file.write("\n")

        for CS_type in sorted(list(cs_list.keys())):
            model_corrs = []

            for model in model_data:
                inner_exp = {}

                for record in cs_list[CS_type]:
                    inner_exp[record.resnum] = model[CS_type][record.resnum]

                model_corrs.append(
                    csx_func.calcCorrel(
                        inner_exp, cs_list[CS_type]
                    )
                )

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
