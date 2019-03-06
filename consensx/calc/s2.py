import prody
import numpy as np

from .vec_3d import Vec3D
from .measure import correlation, q_value, rmsd
import consensx.graph as graph


def s2_values(
    model_data, calculate_on_models, s2_records, s2_type, fit, fit_range
):
    """Returns a dictionary with the average S2 values:
    s2_calced[residue] = value"""
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

                weights = np.zeros((len(ref_chain), 1), dtype=np.int)

                fit_start, fit_end = fit_range.split("-")

                for i in range(int(fit_start) - 1, int(fit_end) - 1):
                    weights[i] = 1

                t = prody.calcTransformation(mob_chain, ref_chain, weights)
                t.apply(mobile)

    # get NH vectors from models (model_data[] -> vectors{resnum : vector})
    vector_data = []
    s2_pairs = {"N": "H", "CA": "HA"}
    h_coords = None
    n_coords = None

    for model_num in calculate_on_models:
        model_data.atomgroup.setACSIndex(model_num)
        current_resindex = 1
        has_first, has_second = False, False
        vectors = {}

        for atom in model_data.atomgroup:
            atom_res = atom.getResnum()

            if atom_res != current_resindex:
                current_resindex = atom_res
                has_first, has_second = False, False

            if atom_res == current_resindex:
                if atom.getName() == s2_type:
                    has_second = True
                    n_coords = Vec3D(atom.getCoords())

                elif atom.getName() == s2_pairs[s2_type]:
                    has_first = True
                    h_coords = Vec3D(atom.getCoords())

                if has_first and has_second:
                    has_first, has_second = False, False
                    vectors[atom_res] = Vec3D(
                        n_coords - h_coords
                    ).normalize()

        vector_data.append(vectors)

    s2_calced = {}

    # iterating over STR records
    for resnum in [int(s2rec.resnum) for s2rec in s2_records]:

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

        s2 = (
            3
            / 2.0
            * (
                x2 ** 2
                + y2 ** 2
                + z2 ** 2
                + 2 * xy ** 2
                + 2 * xz ** 2
                + 2 * yz ** 2
            )
            - 0.5
        )

        s2_calced[resnum] = s2

    return s2_calced


def s2(
    csv_buffer,
    calced_data_storage,
    s2_dict,
    my_path,
    model_data,
    calculate_on_models=None,
    fit=None,
    fit_range=None,
):
    """Back calculate order parameters from given S2 dict and PDB models"""
    s2_data = []

    if not calculate_on_models:
        calculate_on_models = range(model_data.model_count)

    for s2_type in sorted(list(s2_dict.keys())):
        s2_calced = s2_values(
            model_data,
            calculate_on_models,
            s2_dict[s2_type],
            s2_type,
            fit,
            fit_range,
        )

        my_correl = correlation(s2_calced, s2_dict[s2_type])
        my_q_value = q_value(s2_calced, s2_dict[s2_type])
        my_rmsd = rmsd(s2_calced, s2_dict[s2_type])

        corr_key = "S2_" + s2_type + "_corr"
        qval_key = "S2_" + s2_type + "_qval"
        rmsd_key = "S2_" + s2_type + "_rmsd"

        calced_data_storage.update(
            {
                corr_key: "{0}".format("{0:.3f}".format(my_correl)),
                qval_key: "{0}".format("{0:.3f}".format(my_q_value)),
                rmsd_key: "{0}".format("{0:.3f}".format(my_rmsd)),
            }
        )

        csv_buffer.add_data(
            {
                "name": "S2 (" + s2_type + ")",
                "calced": s2_calced,
                "experimental": s2_dict[s2_type],
            }
        )

        print(s2_type + " Order Parameters")
        print("Correl: ", my_correl)
        print("Q-val:  ", my_q_value)
        print("RMSD:   ", my_rmsd)
        print()

        graph_name = "S2_" + s2_type + ".svg"
        graph.values_graph(my_path, s2_calced, s2_dict[s2_type], graph_name)

        corr_graph_name = "S2_corr_" + s2_type + ".svg"
        graph.correl_graph(
            my_path, s2_calced, s2_dict[s2_type], corr_graph_name
        )

        my_id = my_path.split("/")[-2] + "/"

        s2_data.append(
            {
                "s2_type": s2_type,
                "S2_model_n": len(s2_dict[s2_type]),
                "correlation": "{0:.3f}".format(my_correl),
                "q_value": "{0:.3f}".format(my_q_value),
                "rmsd": "{0:.3f}".format(my_rmsd),
                "corr_graph_name": my_id + corr_graph_name,
                "graph_name": my_id + graph_name,
                "input_id": "S2_" + s2_type,
            }
        )

    return s2_data
