import pickle
import os
import threading
import multiprocessing
import subprocess
from time import perf_counter

import consensx.graph as graph

from consensx import thirdparty
from .measure import correlation, q_value, rmsd
from consensx.misc.natural_sort import natural_sort


chemshift_types = ["HA", "CA", "CB", "N", "H", "C"]


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def shiftx_runner(file_list):
    for config in file_list:
        subprocess.call(
            [
                thirdparty.ThirdParty.shiftx,
                "1",
                config["pdb_file"],
                config["out_name"]
            ]
        )


def call_shiftx_on(my_path, pdb_files, bme_weights=None):
    """Call ShiftX on PDB models. Each output is appended to 'out_name'"""
    print("Number of cpu cores:", multiprocessing.cpu_count())

    t1_start = perf_counter()

    shiftx_runs = []

    for i, pdb_file in enumerate(pdb_files):
        shiftx_runs.append({
            "pdb_file": my_path + pdb_file,
            "out_name": my_path + "/modell_" + str(i + 1) + ".out",
        })

    threads = []

    for shiftx_run_chunk in chunks(shiftx_runs, 6):
        thread = threading.Thread(target=shiftx_runner, args=(shiftx_run_chunk,))
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    t1_stop = perf_counter()
    print("[chemshift] Call Shiftx in seconds:", t1_stop - t1_start)

    shiftx_output_files = []
    average_ha, average_h, average_n = {}, {}, {}
    average_ca, average_cb, average_c = {}, {}, {}
    mod_ha, mod_h, mod_n, mod_ca, mod_cb, mod_c = {}, {}, {}, {}, {}, {}
    model_data_list = []
    model_weight = 1
    model_weight_sum = len(pdb_files)

    if bme_weights:
        model_weight_sum = sum(bme_weights)

    for some_file in natural_sort(os.listdir(my_path)):
        if some_file.startswith("modell") and some_file.endswith(".out"):
            shiftx_output_files.append(some_file)

    t1_start = perf_counter()

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

                mod_ha[resnum] = HA * model_weight
                mod_h[resnum] = H * model_weight
                mod_n[resnum] = N * model_weight
                mod_ca[resnum] = CA * model_weight
                mod_cb[resnum] = CB * model_weight
                mod_c[resnum] = C * model_weight

                if resnum in list(average_ha.keys()):
                    average_ha[resnum] += HA * model_weight
                else:
                    average_ha[resnum] = HA * model_weight

                if resnum in list(average_h.keys()):
                    average_h[resnum] += H * model_weight
                else:
                    average_h[resnum] = H * model_weight

                if resnum in list(average_n.keys()):
                    average_n[resnum] += N * model_weight
                else:
                    average_n[resnum] = N * model_weight

                if resnum in list(average_ca.keys()):
                    average_ca[resnum] += CA * model_weight
                else:
                    average_ca[resnum] = CA * model_weight

                if resnum in list(average_cb.keys()):
                    average_cb[resnum] += CB * model_weight
                else:
                    average_cb[resnum] = CB * model_weight

                if resnum in list(average_c.keys()):
                    average_c[resnum] += C * model_weight
                else:
                    average_c[resnum] = C * model_weight

        out_file.close()

        model_data_list.append(
            {
                "HA": mod_ha,
                "H": mod_h,
                "N": mod_n,
                "CA": mod_ca,
                "CB": mod_cb,
                "C": mod_c,
            }
        )

        mod_ha, mod_h, mod_n, mod_ca, mod_cb, mod_c = {}, {}, {}, {}, {}, {}

    t1_stop = perf_counter()
    print("[chemshift] Getting averages from Pales outputs in seconds:", t1_stop - t1_start)

    averages = [
        average_ha,
        average_h,
        average_n,
        average_ca,
        average_cb,
        average_c,
    ]

    t1_start = perf_counter()

    for avg_dict in averages:
        for key in avg_dict:
            avg_dict[key] /= model_weight_sum

    t1_stop = perf_counter()
    print("[chemshift] Calculating averages in seconds:", t1_stop - t1_start)

    return (
        {
            "HA": average_ha,
            "H": average_h,
            "N": average_n,
            "CA": average_ca,
            "CB": average_cb,
            "C": average_c,
        },
        model_data_list,
    )


def chemshifts(
    csv_buffer,
    calced_data_storage,
    chem_shift_lists,
    pdb_models,
    my_path,
    bme_weights,
):
    """Back calculate chemical shifts from given chemical shift list and PDB
       models"""
    cs_data = []
    cs_calced, model_data = call_shiftx_on(my_path, pdb_models, bme_weights)

    cs_model_data_path = my_path + "/ChemShift_model_data.pickle"
    pickle.dump(model_data, open(cs_model_data_path, "wb"))

    for n, cs_list in enumerate(chem_shift_lists):
        bme_exp_filename = "chemshift_" + str(n) + "_exp.dat"
        bme_calc_filename = "chemshift_" + str(n) + "_calc.dat"

        with open(my_path + bme_exp_filename, "w") as exp_dat_file:
            exp_dat_file.write("# DATA=CS PRIOR=GAUSS\n")

            for atom_type in chemshift_types:
                try:
                    for i in cs_list[atom_type]:
                        exp_dat_file.write(
                            str(i.resnum)
                            + "_"
                            + str(i.atom_name)
                            + "\t"
                            + str(i.value)
                            + "\t0.1\n"
                        )
                except KeyError:
                    continue

        with open(my_path + bme_calc_filename, "w") as calc_dat_file:
            for index, model in enumerate(model_data):
                calc_dat_file.write(str(index) + " ")

                for atom_type in chemshift_types:
                    try:
                        for i in cs_list[atom_type]:
                            calc_dat_file.write(
                                str(model[atom_type][i.resnum]) + " "
                            )
                    except KeyError:
                        continue

                calc_dat_file.write("\n")

        for CS_type in sorted(list(cs_list.keys())):
            model_corrs = []

            for model in model_data:
                inner_exp = {}

                for record in cs_list[CS_type]:
                    inner_exp[record.resnum] = model[CS_type][record.resnum]

                model_corrs.append(
                    correlation(inner_exp, cs_list[CS_type])
                )

            avg_model_corr = sum(model_corrs) / len(model_corrs)

            exp_dict = {}

            for record in cs_list[CS_type]:
                exp_dict[record.resnum] = cs_calced[CS_type][record.resnum]

            my_correl = correlation(exp_dict, cs_list[CS_type])
            my_q_value = q_value(exp_dict, cs_list[CS_type])
            my_rmsd = rmsd(exp_dict, cs_list[CS_type])

            corr_key = "CS_" + CS_type + "_corr"
            qval_key = "CS_" + CS_type + "_qval"
            rmsd_key = "CS_" + CS_type + "_rmsd"

            calced_data_storage.update(
                {
                    corr_key: "{0}".format("{0:.3f}".format(my_correl)),
                    qval_key: "{0}".format("{0:.3f}".format(my_q_value)),
                    rmsd_key: "{0}".format("{0:.3f}".format(my_rmsd)),
                }
            )

            csv_buffer.add_data(
                {
                    "name": "ChemShifts (" + CS_type + ")",
                    "calced": exp_dict,
                    "experimental": cs_list[CS_type],
                }
            )

            print("CHEM SHIFT", CS_type)
            print("Correl: ", my_correl)
            print("Q-val:  ", my_q_value)
            print("RMSD:   ", my_rmsd)
            print()

            graph_name = str(n + 1) + "_CS_" + CS_type + ".svg"
            graph.values_graph(my_path, exp_dict, cs_list[CS_type], graph_name)

            corr_graph_name = str(n + 1) + "_CS_corr_" + CS_type + ".svg"
            graph.correl_graph(
                my_path, exp_dict, cs_list[CS_type], corr_graph_name
            )

            mod_corr_graph_name = "CS_mod_corr_" + CS_type + ".svg"
            graph.mod_correl_graph(
                my_path,
                my_correl,
                avg_model_corr,
                model_corrs,
                mod_corr_graph_name,
            )

            my_id = my_path.split("/")[-2] + "/"

            cs_data.append(
                {
                    "CS_type": CS_type,
                    "CS_model_n": len(cs_list[CS_type]),
                    "correlation": "{0:.3f}".format(my_correl),
                    "q_value": "{0:.3f}".format(my_q_value),
                    "rmsd": "{0:.3f}".format(my_rmsd),
                    "corr_graph_name": my_id + corr_graph_name,
                    "graph_name": my_id + graph_name,
                    "mod_corr_graph_name": my_id + mod_corr_graph_name,
                    "input_id": "CS_" + CS_type,
                }
            )

    return cs_data
