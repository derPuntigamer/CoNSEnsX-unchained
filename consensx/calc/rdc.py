import os
import pickle
import re

import consensx.graph as graph

from consensx.csx_libs import methods as csx_func
from consensx.csx_libs import objects as csx_obj


class RDC_model_data():
    """Class for containing per model RDC data"""
    def __init__(self):
        self.rdc_data = {}

    def add_data(self, RDC_list_num, RDC_type, RDC_list_data):
        if RDC_list_num in self.rdc_data.keys():
            self.rdc_data[RDC_list_num][RDC_type] = RDC_list_data
        else:
            self.rdc_data[RDC_list_num] = {}
            self.rdc_data[RDC_list_num][RDC_type] = RDC_list_data


def average_pales_output(pales_out, my_RDC_type):
    """Returns a dictonary with the average RDCs for a given RDC type:
       averageRDC[residue] = value
       and calculated model data as a list of dictonaries
       model_data_list[{1: value}, ...]"""
    pales_out = open(pales_out)
    n_of_structures = 0
    averageRDC = {}
    model_data_list = []
    model_data_dict = {}
    first_run = True

    for line in pales_out:
        if re.match(r"REMARK \d+ couplings", line):
            if first_run:
                first_run = False
                continue

            n_of_structures += 1  # n_of_structures to divide by

            model_data_list.append(model_data_dict)

            if model_data_dict:
                model_data_dict = {}

        elif re.match(r"\s+ \d+", line):
            resnum = int(line.split()[0])
            resnum2 = int(line.split()[3])
            atom = line.split()[2]
            atom2 = line.split()[5]
            D = float(line.split()[8])  # D coloumn of pales output
            RDCtype = str(abs(resnum2 - resnum)) + "_" + atom + "_" + atom2

            # skip non relevant RDC data in the pales output file
            if my_RDC_type != RDCtype:
                continue

            if resnum in list(averageRDC.keys()):
                averageRDC[resnum] += D
            else:
                averageRDC[resnum] = D

            model_data_dict[resnum] = D

    model_data_list.append(model_data_dict)
    n_of_structures += 1
    pales_out.close()

    for res_num in list(averageRDC.keys()):
        averageRDC[res_num] /= n_of_structures

    return averageRDC, model_data_list


def rdc(my_CSV_buffer, RDC_lists, pdb_models, my_path, SVD_enabled, lc_model):
    """Back calculate RDC from given RDC lists and PDB models"""
    my_rdc_model_data = RDC_model_data()
    rdc_calced_data = {}

    for list_num, RDC_dict in enumerate(RDC_lists):
        rdc_calced_data[list_num+1] = []

        # Pales call, results output file "pales.out"
        csx_func.callPalesOn(my_path, pdb_models, RDC_dict,
                             lc_model, SVD_enabled)

        for RDC_type in sorted(list(RDC_dict.keys())):
            print("RDC list", list_num + 1, RDC_type)

            # get averaged RDC values -> averageRDC[residue] = value
            pales_out = my_path + "pales.out"
            averageRDC, model_data = average_pales_output(pales_out, RDC_type)

            model_corrs = []

            for model in model_data:
                model_corrs.append(
                    csx_func.calcCorrel(model, RDC_dict[RDC_type])
                )

            my_rdc_model_data.add_data(list_num + 1, RDC_type, model_data)
            csx_obj.RDC_modell_corr(model_corrs)

            avg_model_corr = sum(model_corrs) / len(model_corrs)

            # removing records from other RDC types
            my_averageRDC = {}

            for record in RDC_dict[RDC_type]:
                my_averageRDC[record.resnum] = averageRDC[record.resnum]

            correl = csx_func.calcCorrel(my_averageRDC, RDC_dict[RDC_type])
            q_value = csx_func.calcQValue(my_averageRDC, RDC_dict[RDC_type])
            rmsd = csx_func.calcRMSD(my_averageRDC, RDC_dict[RDC_type])

            RDC_simple = RDC_type.replace('_', '')

            # TODO DB upload!
            corr_key = "RDC_" + str(list_num + 1) + "_" + RDC_simple + "_corr"
            qval_key = "RDC_" + str(list_num + 1) + "_" + RDC_simple + "_qval"
            rmsd_key = "RDC_" + str(list_num + 1) + "_" + RDC_simple + "_rmsd"

            csx_obj.CalcPickle.data.update(
                {
                    corr_key: "{0}".format('{0:.3f}'.format(correl)),
                    qval_key: "{0}".format('{0:.3f}'.format(q_value)),
                    rmsd_key: "{0}".format('{0:.3f}'.format(rmsd))
                }
            )

            my_CSV_buffer.csv_data.append({
                "name": "RDC_" + str(list_num + 1) + "(" + RDC_type + ")",
                "calced": my_averageRDC,
                "experimental": RDC_dict[RDC_type]
            })

            print("Correl: ", correl)
            print("Q-val:  ", q_value)
            print("RMSD:   ", rmsd)
            print()

            graph_name = str(list_num + 1) + "_RDC_" + RDC_type + ".svg"
            graph.values(
                my_path, my_averageRDC, RDC_dict[RDC_type], graph_name
            )

            corr_graph_name = (
                str(list_num + 1) + "_RDC_corr_" + RDC_type + ".svg"
            )

            graph.correl_graph(
                my_path, my_averageRDC, RDC_dict[RDC_type], corr_graph_name
            )

            mod_corr_graph_name = (
                str(list_num + 1) + "_RDC_mod_corr_" + RDC_type + ".svg"
            )

            graph.mod_correl_graph.modCorrelGraph(
                my_path, correl, avg_model_corr, model_corrs,
                mod_corr_graph_name
            )

            my_id = my_path.split('/')[-2] + '/'

            rdc_calced_data[list_num+1].append({
                "RDC_type": RDC_type,
                "RDC_model_n": len(RDC_dict[RDC_type]),
                "correlation": '{0:.3f}'.format(correl),
                "q_value": '{0:.3f}'.format(q_value),
                "rmsd": '{0:.3f}'.format(rmsd),
                "corr_graph_name": my_id + corr_graph_name,
                "graph_name": my_id + graph_name,
                "mod_corr_graph_name": my_id + mod_corr_graph_name,
                "input_id": "RDC_" + str(list_num + 1) + "_" +
                            "".join(RDC_type.split('_'))
            })

        model_data_path = my_path + "/RDC_model_data.pickle"
        pickle.dump(my_rdc_model_data, open(model_data_path, "wb"))
        os.remove(pales_out)

    return rdc_calced_data
