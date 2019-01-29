import pickle

from consensx.csx_libs import methods as csx_func
from consensx.csx_libs import objects as csx_obj


def chemshifts(my_CSV_buffer, ChemShift_lists, pdb_models, my_path):
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
