import pickle

import consensx.graph as graph

from consensx.csx_libs import methods as csx_func
from consensx.csx_libs import objects as csx_obj


def jcoupling(my_CSV_buffer, param_set, Jcoup_dict, my_PDB, my_path):
    """Back calculate skalar coupling from given RDC lists and PDB models"""
    Jcuop_data = []
    type_dict = {}
    dihed_lists = csx_func.calcDihedAngles()

    for Jcoup_type in sorted(list(Jcoup_dict.keys())):
        JCoup_calced, model_data = csx_func.calcJCoup(
            param_set, dihed_lists, Jcoup_dict[Jcoup_type], Jcoup_type
        )

        type_dict[Jcoup_type] = model_data
        model_corrs = []

        for model in model_data:
            model_corrs.append(
                csx_func.calcCorrel(model, Jcoup_dict[Jcoup_type])
            )

        avg_model_corr = sum(model_corrs) / len(model_corrs)

        correl = csx_func.calcCorrel(JCoup_calced, Jcoup_dict[Jcoup_type])
        q_value = csx_func.calcQValue(JCoup_calced, Jcoup_dict[Jcoup_type])
        rmsd = csx_func.calcRMSD(JCoup_calced, Jcoup_dict[Jcoup_type])

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
        graph.values_graph(
            my_path, JCoup_calced, Jcoup_dict[Jcoup_type], graph_name
        )

        corr_graph_name = "JCoup_corr_" + Jcoup_type + ".svg"
        graph.correl_graph(
            my_path, JCoup_calced, Jcoup_dict[Jcoup_type], corr_graph_name
        )

        mod_corr_graph_name = "JCoup_mod_corr_" + Jcoup_type + ".svg"
        graph.mod_correl_graph(
            my_path, correl, avg_model_corr, model_corrs, mod_corr_graph_name
        )

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
