import prody
import math

import matplotlib.pyplot as plt
import numpy as np

from consensx.csx_libs import objects as csx_obj


def s2_sidechain(csv_buffer, s2_sidechain, my_path, model_data, fit=None):
    """Back calculate order parameters from given S2 dict and PDB models"""

    sc_lot = {
        "VAL": {"CG1": "CB", "CG2": "CB"},
        "ILE": {"CG2": "CB", "CD": "CG1", "CD1": "CG1"},
        "THR": {"CG2": "CB"},
        "LEU": {"CD1": "CG", "CD2": "CG"},
        "ALA": {"CB": "CA"},
        "MET": {"CE": "SD"},
    }

    if fit:
        model_data.atomgroup.setACSIndex(0)
        prody.alignCoordsets(model_data.atomgroup.calpha)

    for record in s2_sidechain:
        vectors = []
        resnum = record.resnum
        my_type = record.type
        my_res = model_data.atomgroup[("A", resnum)].getResname()

        # find pair for measured aa
        pair = sc_lot[my_res][my_type]

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

    sidechain_exp1, sidechain_exp2 = [], []
    sidechain_calc1, sidechain_calc2 = {}, {}

    prev_resnum = -100000

    for record in s2_sidechain:
        if record.resnum != prev_resnum:
            sidechain_exp1.append(record)
            sidechain_calc1[record.resnum] = record.calced
            prev_resnum = record.resnum
        else:
            sidechain_exp2.append(record)
            sidechain_calc2[record.resnum] = record.calced

    csv_buffer.add_data(
        {
            "name": "S2_meth",
            "calced": sidechain_calc1,
            "experimental": sidechain_exp1,
        }
    )

    if sidechain_exp2:
        csv_buffer.add_data(
            {
                "name": "S2_meth (cont)",
                "calced": sidechain_calc2,
                "experimental": sidechain_exp2,
            }
        )

    # correlation calculation
    m = [0.0, 0.0, 0.0]
    d = [0.0, 0.0]

    for record in s2_sidechain:
        exp = record.value
        calc = record.calced

        m[0] += calc
        m[1] += exp
        m[2] += calc * exp

    m[0] /= len(s2_sidechain)
    m[1] /= len(s2_sidechain)
    m[2] /= len(s2_sidechain)

    for record in s2_sidechain:
        exp = record.value
        calc = record.calced

        d[0] += (calc - m[0]) ** 2
        d[1] += (exp - m[1]) ** 2

    d[0] /= len(s2_sidechain)
    d[0] = math.sqrt(d[0])
    d[1] /= len(s2_sidechain)
    d[1] = math.sqrt(d[1])

    correl = (m[2] - (m[0] * m[1])) / (d[0] * d[1])
    print("Corr: ", correl)

    # Q-value calculation
    d2, e2 = 0, 0

    for record in s2_sidechain:
        exp = record.value
        calc = record.calced

        d2 += (calc - exp) ** 2
        e2 += exp ** 2

    Q = 100 * math.sqrt(d2) / math.sqrt(e2)
    q_value = round(Q, 6)
    print("Q-value: ", q_value)

    # RMSD calculation
    d2 = 0

    for record in s2_sidechain:
        exp = record.value
        calc = record.calced

        d2 += (calc - exp) ** 2

    rmsd = math.sqrt(d2 / len(s2_sidechain))
    print("RMSD: ", round(rmsd, 6))

    exp_values = []
    calced_values = []

    for record in s2_sidechain:
        exp_values.append(record.value)
        calced_values.append(record.calced)

    min_calc = min(calced_values)
    max_calc = max(calced_values)

    min_exp = min(exp_values)
    max_exp = max(exp_values)
    miny = min(min_calc, min_exp)  # get minimum value
    maxy = max(max_calc, max_exp)  # get maximum value

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
    plt.plot(diag, diag, linewidth=2.0, color="red", alpha=0.7)
    plt.plot(exp_values, calced_values, "bo")
    plt.axis([miny, maxy, miny, maxy])
    plt.xlabel("experimental")
    plt.ylabel("calculated")
    plt.tight_layout(pad=1.08)
    plt.savefig(my_path + "/" + "S2_sc_corr.svg", format="svg")
    plt.close()

    xs = []
    prev_resnum2 = -1

    for record in s2_sidechain:
        if record.resnum != prev_resnum2:
            xs.append(record.resnum)
            prev_resnum2 = record.resnum
        else:
            xs.append(record.resnum + 0.3)

    print("XS AXIS", xs)
    print("len SC", len(s2_sidechain))
    print("len xs", len(xs))

    plt.figure(figsize=(10, 5), dpi=80)
    plt.plot(
        xs,
        exp_values,
        linewidth=2.0,
        color="red",
        marker="o",
        label="exp",
        alpha=0.7,
    )
    plt.plot(
        xs,
        calced_values,
        linewidth=2.0,
        color="blue",
        marker="o",
        label="calc",
        alpha=0.7,
    )
    plt.legend(loc="lower left")
    plt.xlabel("residue number")
    plt.ylabel("value")
    ax = plt.axes()
    ax.yaxis.grid()
    plt.tight_layout(pad=1.08)
    plt.savefig(my_path + "/" + "S2_sc_graph.svg", format="svg")
    plt.close()

    my_id = my_path.split("/")[-2] + "/"

    print("CORR GRAPH", my_id + "S2_sc_corr.svg")
    print("GRAPH", my_id + "S2_sc_graph.svg")

    s2_sc_data = {
        "S2_model_n": len(s2_sidechain),
        "correlation": "{0:.3f}".format(correl),
        "q_value": "{0:.3f}".format(q_value),
        "rmsd": "{0:.3f}".format(rmsd),
        "corr_graph_name": my_id + "S2_sc_corr.svg",
        "graph_name": my_id + "S2_sc_graph.svg",
    }

    return s2_sc_data
