import prody
import math

import matplotlib.pyplot as plt
import numpy as np

from consensx.csx_libs import objects as csx_obj


def s2_sidechain(my_CSV_buffer, S2_sidechain, my_path, fit=None):
    """Back calculate order paramteres from given S2 dict and PDB models"""

    sc_LOT = {
        'VAL': {'CG1': 'CB', 'CG2': 'CB'},
        'ILE': {'CG2': 'CB', 'CD': 'CG1', 'CD1': 'CG1'},
        'THR': {'CG2': 'CB'},
        'LEU': {'CD1': 'CG', 'CD2': 'CG'},
        'ALA': {'CB': 'CA'},
        'MET': {'CE': 'SD'}
    }

    model_data = csx_obj.PDB_model.model_data

    if fit and not csx_obj.PDB_model.is_fitted:
        model_data.atomgroup.setACSIndex(0)
        prody.alignCoordsets(model_data.atomgroup.calpha)

    for record in S2_sidechain:
        vectors = []
        resnum = record.resnum
        my_type = record.type
        my_res = model_data.atomgroup[('A', resnum)].getResname()

        # find pair for measured aa
        pair = sc_LOT[my_res][my_type]

        for model_num in range(model_data.model_count):
            model_data.atomgroup.setACSIndex(model_num)

            try:
                sel = "resnum {} name {}".format(resnum, my_type)
                coords = model_data.atomgroup.select(sel).getCoords()[0]
                sel = "resnum {} name {}".format(resnum, pair)
                pair_coords = model_data.atomgroup.select(sel).getCoords()[0]
            except AttributeError:
                return {
                    "error": "Sidechain order parameter atom name not found\
                    in PDB. Please check your atom naming."
                }

            vectors.append(csx_obj.Vec_3D(coords - pair_coords).normalize())

        x2, y2, z2, xy, xz, yz = 0, 0, 0, 0, 0, 0

        for vector in vectors:
            x, y, z = vector.v[0], vector.v[1], vector.v[2]

            x2 += x ** 2
            y2 += y ** 2
            z2 += z ** 2
            xy += x * y
            xz += x * z
            yz += y * z

        x2 /= len(vectors)
        y2 /= len(vectors)
        z2 /= len(vectors)
        xy /= len(vectors)
        xz /= len(vectors)
        yz /= len(vectors)

        s2 = 3 / 2.0 * (x2 ** 2 + y2 ** 2 + z2 ** 2 +
                        2 * xy ** 2 + 2 * xz ** 2 + 2 * yz ** 2) - 0.5

        record.calced = s2

    # prepare data for CSV object inicialization
    sidechain_exp1, sidechain_exp2 = [], []
    sidechain_calc1, sidechain_calc2 = {}, {}

    prev_resnum = -100000

    for record in S2_sidechain:
        if record.resnum != prev_resnum:
            sidechain_exp1.append(record)
            sidechain_calc1[record.resnum] = record.calced
            prev_resnum = record.resnum
        else:
            sidechain_exp2.append(record)
            sidechain_calc2[record.resnum] = record.calced

    my_CSV_buffer.add_data({
        "name": "S2_meth",
        "calced": sidechain_calc1,
        "experimental": sidechain_exp1
    })

    if sidechain_exp2:
        my_CSV_buffer.add_data({
            "name": "S2_meth (cont)",
            "calced": sidechain_calc2,
            "experimental": sidechain_exp2
        })

    # correlation calculation
    M = [0.0, 0.0, 0.0]
    D = [0.0, 0.0]

    for record in S2_sidechain:
        exp = record.value
        calc = record.calced

        M[0] += calc
        M[1] += exp
        M[2] += calc * exp

    M[0] /= len(S2_sidechain)
    M[1] /= len(S2_sidechain)
    M[2] /= len(S2_sidechain)

    for record in S2_sidechain:
        exp = record.value
        calc = record.calced

        D[0] += (calc - M[0]) ** 2
        D[1] += (exp - M[1]) ** 2

    D[0] /= len(S2_sidechain)
    D[0] = math.sqrt(D[0])
    D[1] /= len(S2_sidechain)
    D[1] = math.sqrt(D[1])

    correl = (M[2] - (M[0] * M[1])) / (D[0] * D[1])
    print("Corr: ", correl)

    # Q-value calculation
    D2, E2 = 0, 0

    for record in S2_sidechain:
        exp = record.value
        calc = record.calced

        D2 += (calc - exp) ** 2
        E2 += exp ** 2

    Q = 100 * math.sqrt(D2) / math.sqrt(E2)
    q_value = round(Q, 6)
    print("Q-value: ", q_value)

    # RMDS calculation
    D2 = 0

    for record in S2_sidechain:
        exp = record.value
        calc = record.calced

        D2 += (calc - exp) ** 2

    rmsd = math.sqrt(D2 / len(S2_sidechain))
    print("RMSD: ", round(rmsd, 6))

    exp_values = []
    calced_values = []

    for record in S2_sidechain:
        exp_values.append(record.value)
        calced_values.append(record.calced)

    min_calc = min(calced_values)
    max_calc = max(calced_values)

    min_exp = min(exp_values)
    max_exp = max(exp_values)
    miny = min(min_calc, min_exp)             # get minimum value
    maxy = max(max_calc, max_exp)             # get maximum value

    diag = []

    margin = int(abs(miny - maxy) * 0.05)

    if abs(miny - maxy) < 10:
        margin = 0.3
    elif abs(miny - maxy) < 2:
        margin = 0.01
    elif abs(miny - maxy) < 1:
        margin = 0

    maxy += margin
    miny -= margin

    for i in np.arange(miny, maxy * 1.42, 0.1):  # draw graph diagonal
        diag.append(i)

    plt.figure(figsize=(6, 5), dpi=80)
    plt.plot(diag, diag, linewidth=2.0, color='red', alpha=.7)
    plt.plot(exp_values, calced_values, 'bo')
    plt.axis([miny, maxy, miny, maxy])
    plt.xlabel('experimental')
    plt.ylabel('calculated')
    plt.tight_layout(pad=1.08)
    plt.savefig(my_path + "/" + "S2_sc_corr.svg", format="svg")
    plt.close()

    xs = []
    prev_resnum2 = -1

    for record in S2_sidechain:
        if record.resnum != prev_resnum2:
            xs.append(record.resnum)
            prev_resnum2 = record.resnum
        else:
            xs.append(record.resnum + 0.3)

    print("XS AXIS", xs)
    print("len SC", len(S2_sidechain))
    print("len xs", len(xs))

    plt.figure(figsize=(10, 5), dpi=80)
    plt.plot(xs, exp_values,
             linewidth=2.0, color='red', marker='o', label='exp', alpha=.7)
    plt.plot(xs, calced_values,
             linewidth=2.0, color='blue', marker='o', label='calc', alpha=.7)
    plt.legend(loc='lower left')
    plt.xlabel('residue number')
    plt.ylabel('value')
    ax = plt.axes()
    ax.yaxis.grid()
    plt.tight_layout(pad=1.08)
    plt.savefig(my_path + "/" + "S2_sc_graph.svg", format="svg")
    plt.close()

    my_id = my_path.split('/')[-2] + '/'

    print("CORR GRAPH", my_id + "S2_sc_corr.svg")
    print("GRAPH", my_id + "S2_sc_graph.svg")

    S2_sc_data = {
        "S2_model_n": len(S2_sidechain),
        "correlation": '{0:.3f}'.format(correl),
        "q_value": '{0:.3f}'.format(q_value),
        "rmsd": '{0:.3f}'.format(rmsd),
        "corr_graph_name": my_id + "S2_sc_corr.svg",
        "graph_name": my_id + "S2_sc_graph.svg"
    }

    return S2_sc_data
