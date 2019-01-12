import os
import pickle

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

    def get_best_RDC_model(self, measure, user_sel, RDC_lists):
        model_scores = {}
        divide_by = 0.0

        for sel_data in user_sel:
            if sel_data[0] != "RDC":
                continue

            my_data = self.rdc_data[sel_data[1]][sel_data[2]]
            experimental = RDC_lists[sel_data[1] - 1][sel_data[2]]

            for model_num, model in enumerate(my_data):
                if measure == "correlation":
                    calced = csx_func.calcCorrel(model, experimental)
                elif measure == "q-value":
                    calced = csx_func.calcQValue(model, experimental)
                elif measure == "rmsd":
                    calced = csx_func.calcRMSD(model, experimental)

                if model_num in model_scores.keys():
                    model_scores[model_num] += calced * sel_data[3]
                else:
                    model_scores[model_num] = calced * sel_data[3]

            divide_by += sel_data[3]

        max_value, max_loc = -1000, -1
        min_value, min_loc = 1000, -1

        for loc, key in enumerate(model_scores.keys()):
            model_scores[key] /= divide_by

            if model_scores[key] > max_value:
                max_value = model_scores[key]
                max_loc = loc
            if model_scores[key] < min_value:
                min_value = model_scores[key]
                min_loc = loc

        if measure == "correlation":
            return max_loc
        elif measure == "q-value" or measure == "rmsd":
            return min_loc


def rdc(my_CSV_buffer, RDC_lists, pdb_models,
            my_path, SVD_enabled, lc_model):
    """Back calculate RDC from given RDC lists and PDB models"""
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
            averageRDC, model_data = csx_func.avgPalesRDCs(pales_out, RDC_type)

            model_corrs = []

            for model in model_data:
                model_corrs.append(
                    csx_func.calcCorrel(model, RDC_dict[RDC_type])
                )

            my_rdc_model_data = RDC_model_data()

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
            csx_func.makeGraph(my_path, my_averageRDC, RDC_dict[RDC_type],
                               graph_name)

            corr_graph_name = (str(list_num + 1) + "_RDC_corr_" +
                               RDC_type + ".svg")
            csx_func.makeCorrelGraph(my_path, my_averageRDC,
                                     RDC_dict[RDC_type], corr_graph_name)

            mod_corr_graph_name = (str(list_num + 1) + "_RDC_mod_corr_" +
                                   RDC_type + ".svg")
            csx_func.modCorrelGraph(my_path, correl, avg_model_corr,
                                    model_corrs, mod_corr_graph_name)

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
        RDC_obj_attr = my_path + "/RDC_obj_attr.pickle"
        pickle.dump(my_rdc_model_data.rdc_data, open(RDC_obj_attr, "wb"))
        os.remove(pales_out)

    return rdc_calced_data