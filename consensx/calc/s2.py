import prody
import numpy as np

from consensx.csx_libs import methods as csx_func
from consensx.csx_libs import objects as csx_obj

import consensx.graph as graph


def s2_values(
        model_data, calculate_on_models, S2_records, S2_type, fit, fit_range
        ):
    """Returns a dictonary with the average S2 values:
    S2_calced[residue] = value"""

    if fit:
        reference = model_data.atomgroup[:]

        model_data.atomgroup.setACSIndex(0)
        prody.alignCoordsets(model_data.atomgroup.calpha)

        if fit_range:
            for model_num in calculate_on_models:
                model_data.atomgroup.setACSIndex(model_num)

                mobile = model_data.atomgroup[:]
                matches = prody.matchChains(reference, mobile)
                match = matches[0]
                ref_chain = match[0]
                mob_chain = match[1]

                # if fit_range:
                weights = np.zeros((len(ref_chain), 1), dtype=np.int)

                fit_start, fit_end = fit_range.split('-')

                for i in range(int(fit_start) - 1, int(fit_end) - 1):
                    weights[i] = 1

                t = prody.calcTransformation(mob_chain, ref_chain, weights)
                t.apply(mobile)

    # get NH vectors from models (model_data[] -> vectors{resnum : vector})
    vector_data = []
    s2_pairs = {'N': 'H', 'CA': 'HA'}

    for model_num in calculate_on_models:
        model_data.atomgroup.setACSIndex(model_num)
        current_Resindex = 1
        has_first, has_second = False, False
        vectors = {}

        for atom in model_data.atomgroup:
            atom_res = atom.getResnum()

            if atom_res != current_Resindex:
                current_Resindex = atom_res
                has_first, has_second = False, False

            if atom_res == current_Resindex:
                if atom.getName() == S2_type:
                    has_second = True
                    N_coords = csx_obj.Vec_3D(atom.getCoords())

                elif atom.getName() == s2_pairs[S2_type]:
                    has_first = True
                    H_coords = csx_obj.Vec_3D(atom.getCoords())

                if has_first and has_second:
                    has_first, has_second = False, False
                    vectors[atom_res] = csx_obj.Vec_3D(
                        N_coords - H_coords
                    ).normalize()

        vector_data.append(vectors)

    S2_calced = {}

    # iterating over STR records
    for resnum in [int(s2rec.resnum) for s2rec in S2_records]:

        x2, y2, z2, xy, xz, yz = 0, 0, 0, 0, 0, 0

        # iterating over PDB models
        for m in vector_data:

            # coordinates in model at a given resnum
            x, y, z = m[resnum].v[0], m[resnum].v[1], m[resnum].v[2]

            x2 += x ** 2
            y2 += y ** 2
            z2 += z ** 2
            xy += x * y
            xz += x * z
            yz += y * z

        x2 /= len(vector_data)
        y2 /= len(vector_data)
        z2 /= len(vector_data)
        xy /= len(vector_data)
        xz /= len(vector_data)
        yz /= len(vector_data)

        # S2 calcuation
        s2 = 3 / 2.0 * (x2 ** 2 + y2 ** 2 + z2 ** 2 +
                        2 * xy ** 2 + 2 * xz ** 2 + 2 * yz ** 2) - 0.5

        S2_calced[resnum] = s2

    return S2_calced


def s2(my_CSV_buffer, S2_dict, my_path, calculate_on_models=None,
        fit=None, fit_range=None):
    """Back calculate order paramteres from given S2 dict and PDB models"""
    S2_data = []
    model_data = csx_obj.PDB_model.model_data

    if not calculate_on_models:
        calculate_on_models = range(model_data.model_count)

    for S2_type in sorted(list(S2_dict.keys())):
        S2_calced = s2_values(
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
        graph.values(my_path, S2_calced, S2_dict[S2_type], graph_name)

        corr_graph_name = "S2_corr_" + S2_type + ".svg"
        graph.correl_graph(
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
