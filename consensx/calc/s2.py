from consensx.csx_libs import methods as csx_func
from consensx.csx_libs import objects as csx_obj

def s2(
        my_CSV_buffer, S2_dict, my_path,
        calculate_on_models=None,
        fit=None,
        fit_range=None
    ):
    """Back calculate order paramteres from given S2 dict and PDB models"""
    S2_data = []
    model_data = csx_obj.PDB_model.model_data

    if not calculate_on_models:
        calculate_on_models = range(model_data.model_count)

    for S2_type in sorted(list(S2_dict.keys())):
        S2_calced = csx_func.calcS2(
            model_data, calculate_on_models, S2_dict[S2_type], S2_type,
            fit, fit_range
        )

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
        csx_func.makeCorrelGraph(
            my_path, S2_calced, S2_dict[S2_type], corr_graph_name
        )

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