import pickle
import math

import consensx.graph as graph

from consensx.csx_libs import methods as csx_func
from consensx.csx_libs import objects as csx_obj

# from consensx.bme_reweight import Reweight

# Equation and coefficients from:
# Wang & Bax (1996) JACS 118:2483-2494. Table 1, NMR + X-ray data
Jcoup_dict1 = {
    "A": {"3JHNCB": 3.39, "3JHNHA": 6.98, "3JHNC": 4.32, "3JHAC": 3.75},
    "B": {"3JHNCB": -0.94, "3JHNHA": -1.38, "3JHNC": 0.84, "3JHAC": 2.19},
    "C": {"3JHNCB": 0.07, "3JHNHA": 1.72, "3JHNC": 0.00, "3JHAC": 1.28},
    "THETA": {
        "3JHNCB": math.radians(60),
        "3JHNHA": math.radians(-60),
        "3JHNC": math.radians(0),
        "3JHAC": math.radians(-60),
    },  # RAD!
}

# J. Am. Chem. Soc., Vol. 119, No. 27, 1997; Table 2 -> solution
Jcoup_dict2 = {
    "A": {"3JHNCB": 3.06, "3JHNHA": 7.13, "3JHNC": 4.19, "3JHAC": 3.84},
    "B": {"3JHNCB": -0.74, "3JHNHA": 1.31, "3JHNC": 0.99, "3JHAC": 2.19},
    "C": {"3JHNCB": 0.10, "3JHNHA": 1.56, "3JHNC": 0.03, "3JHAC": 1.20},
    "THETA": {
        "3JHNCB": math.radians(60),
        "3JHNHA": math.radians(-60),
        "3JHNC": math.radians(180),
        "3JHAC": math.radians(120),
    },  # RAD!
}

# https://x86.cs.duke.edu/~brd/Teaching/Bio/asmb/Papers/NMR/nilges-jmr05.pdf
Jcoup_dict3 = {
    "A": {"3JHNCB": 3.26, "3JHNHA": 7.13, "3JHNC": 4.19, "3JHAC": 3.84},
    "B": {"3JHNCB": -0.87, "3JHNHA": -1.31, "3JHNC": 0.99, "3JHAC": 2.19},
    "C": {"3JHNCB": 0.10, "3JHNHA": 1.56, "3JHNC": 0.03, "3JHAC": 1.20},
    "THETA": {
        "3JHNCB": math.radians(60),
        "3JHNHA": math.radians(-60),
        "3JHNC": math.radians(0),
        "3JHAC": math.radians(-60),
    },  # RAD!
}


def calc_dihedral_angles(pdb_model_data):
    """Calculates backbone diherdral angles
       note: all returned angle values are in radian"""

    j_coup_dicts = []

    for i in range(pdb_model_data.coordsets):
        pdb_model_data.atomgroup.setACSIndex(i)
        current_resindex = 1
        prev_c, my_n, my_ca, my_c = None, None, None, None
        j_coup_dict = {}

        for atom in pdb_model_data.atomgroup:
            atom_res = atom.getResindex() + 1

            if atom_res != current_resindex:

                if (
                    prev_c is not None
                    and my_n is not None
                    and my_ca is not None
                    and my_c is not None
                ):

                    nca_vec = my_n - my_ca
                    cn_vec = prev_c - my_n
                    cca_vec = my_c - my_ca

                    first_cross = csx_obj.Vec_3D.cross(cn_vec, nca_vec)
                    second_cross = csx_obj.Vec_3D.cross(cca_vec, nca_vec)

                    angle = csx_obj.Vec_3D.dihedAngle(
                        first_cross, second_cross
                    )

                    # reference for setting sign of angle
                    reference = csx_obj.Vec_3D.cross(first_cross, second_cross)

                    r1 = reference.normalize()
                    r2 = nca_vec.normalize()

                    if (r1 - r2).magnitude() < r2.magnitude():
                        angle *= -1

                    j_coup_dict[current_resindex] = -1 * math.radians(angle)

                current_resindex = atom_res
                prev_c = my_c
                my_n, my_ca, my_c = None, None, None

            if atom_res == current_resindex:
                if atom.getName() == "N":
                    my_n = csx_obj.Vec_3D(atom.getCoords())
                elif atom.getName() == "CA":
                    my_ca = csx_obj.Vec_3D(atom.getCoords())
                elif atom.getName() == "C":
                    my_c = csx_obj.Vec_3D(atom.getCoords())

        j_coup_dicts.append(j_coup_dict)

    return j_coup_dicts


def calc_jcoupling(
    param_set,
    calced,
    experimental,
    jcoup_type,
    my_path,
    exp_dat_file_name,
    bme_weights,
):
    """Calculates J-coupling values from dihedral angles
       note: all angles must be in radian"""
    jcoup_calced = {}

    if param_set == 1:
        my_karplus = Jcoup_dict1
    elif param_set == 2:
        my_karplus = Jcoup_dict2
    elif param_set == 3:
        my_karplus = Jcoup_dict3

    A = my_karplus["A"]
    B = my_karplus["B"]
    C = my_karplus["C"]
    THETA = my_karplus["THETA"]

    jcoup_values_store = []

    for record in experimental:  # resnums
        J = 0
        jcoup_values_list = []

        for i, my_dict in enumerate(calced):  # lists (with models as dicts)
            phi = my_dict[record.resnum]

            jcoup = (
                    A[jcoup_type] * (math.cos(phi + THETA[jcoup_type])) ** 2
                    + B[jcoup_type] * math.cos(phi + THETA[jcoup_type])
                    + C[jcoup_type]
            )

            jcoup_values_list.append(str(jcoup))

            if bme_weights:
                J += bme_weights[i] * jcoup
            else:
                J += jcoup

        jcoup_values_store.append(jcoup_values_list)

        if bme_weights:
            jcoup_calced[record.resnum] = J / sum(bme_weights)
        else:
            jcoup_calced[record.resnum] = J / len(calced)

    calc_dat_file_name = my_path + "jcoup_" + jcoup_type + "_calc.dat"

    with open(calc_dat_file_name, "w") as jcoup_calc_file:
        for mod in range(len(jcoup_values_store[0])):
            mod_str = str(mod) + " "

            for exp in range(len(jcoup_values_store)):
                mod_str += jcoup_values_store[exp][mod] + " "

            jcoup_calc_file.write(mod_str + "\n")

    # BME reweighting ingetration, might DELETE later
    # rew = Reweight()
    # rew.load(exp_dat_file_name, calc_dat_file_name)
    # chi2_before, chi2_after, srel = rew.optimize(theta=40)

    # print("# CHI2 before minimization:     %8.4f" % (chi2_before))
    # print("# CHI2 after minimization:      %8.4f" % (chi2_after))

    # w_opt = rew.get_weights()
    # weights_file_name = my_path + "jcoup_" + Jcoup_type + "_weights.dat"

    # with open(weights_file_name, "w") as weights_file:
    #     weights_file.write(" ".join([str(i) for i in w_opt]))

    model_data_list = []
    model_data_dict = {}

    for Jcoup_dict in calced:  # model
        for record in experimental:
            phi = Jcoup_dict[record.resnum]

            J = (
                    A[jcoup_type] * (math.cos(phi + THETA[jcoup_type])) ** 2
                    + B[jcoup_type] * math.cos(phi + THETA[jcoup_type])
                    + C[jcoup_type]
            )

            model_data_dict[record.resnum] = J

        model_data_list.append(model_data_dict)
        model_data_dict = {}

    return jcoup_calced, model_data_list


def jcoupling(
    csv_buffer,
    calced_data_storage,
    pdb_model_data,
    db_entry,
    jcoup_dict,
    my_path,
    bme_weights,
):
    """Back calculate skalar coupling from given RDC lists and PDB models"""
    param_set = db_entry.karplus
    jcuop_data = []
    type_dict = {}
    dihed_lists = calc_dihedral_angles(pdb_model_data)

    for Jcoup_type in sorted(list(jcoup_dict.keys())):
        exp_dat_file_name = my_path + "jcoup_" + Jcoup_type + "_exp.dat"
        jcoup_exp_file = open(exp_dat_file_name, "w")

        jcoup_exp_file.write("# DATA=JCOUPLINGS PRIOR=GAUSS\n")

        for i in jcoup_dict[Jcoup_type]:
            jcoup_exp_file.write(
                str(i.resnum)
                + "-"
                + str(i.type)
                + "\t"
                + str(i.value)
                + "\t0.1\n"
            )

        jcoup_exp_file.close()

        jcoup_calced, model_data = calc_jcoupling(
            param_set,
            dihed_lists,
            jcoup_dict[Jcoup_type],
            Jcoup_type,
            my_path,
            exp_dat_file_name,
            bme_weights,
        )

        type_dict[Jcoup_type] = model_data
        model_corrs = []

        for model in model_data:
            model_corrs.append(
                csx_func.calcCorrel(model, jcoup_dict[Jcoup_type])
            )

        avg_model_corr = sum(model_corrs) / len(model_corrs)

        correl = csx_func.calcCorrel(jcoup_calced, jcoup_dict[Jcoup_type])
        q_value = csx_func.calcQValue(jcoup_calced, jcoup_dict[Jcoup_type])
        rmsd = csx_func.calcRMSD(jcoup_calced, jcoup_dict[Jcoup_type])

        corr_key = "JCoup_" + Jcoup_type + "_corr"
        qval_key = "JCoup_" + Jcoup_type + "_qval"
        rmsd_key = "JCoup_" + Jcoup_type + "_rmsd"

        calced_data_storage.update(
            {
                corr_key: "{0}".format("{0:.3f}".format(correl)),
                qval_key: "{0}".format("{0:.3f}".format(q_value)),
                rmsd_key: "{0}".format("{0:.3f}".format(rmsd)),
            }
        )

        csv_buffer.add_data(
            {
                "name": "J-couplings (" + Jcoup_type + ")",
                "calced": jcoup_calced,
                "experimental": jcoup_dict[Jcoup_type],
            }
        )

        print("J-couplings (" + Jcoup_type + ")")
        print("Correl: ", correl)
        print("Q-val:  ", q_value)
        print("RMSD:   ", rmsd)
        print()

        graph_name = "JCoup_" + Jcoup_type + ".svg"
        graph.values_graph(
            my_path, jcoup_calced, jcoup_dict[Jcoup_type], graph_name
        )

        corr_graph_name = "JCoup_corr_" + Jcoup_type + ".svg"
        graph.correl_graph(
            my_path, jcoup_calced, jcoup_dict[Jcoup_type], corr_graph_name
        )

        mod_corr_graph_name = "JCoup_mod_corr_" + Jcoup_type + ".svg"
        graph.mod_correl_graph(
            my_path, correl, avg_model_corr, model_corrs, mod_corr_graph_name
        )

        my_id = my_path.split("/")[-2] + "/"

        jcuop_data.append(
            {
                "Jcoup_type": Jcoup_type,
                "Jcoup_model_n": len(jcoup_dict[Jcoup_type]),
                "correlation": "{0:.3f}".format(correl),
                "q_value": "{0:.3f}".format(q_value),
                "rmsd": "{0:.3f}".format(rmsd),
                "corr_graph_name": my_id + corr_graph_name,
                "graph_name": my_id + graph_name,
                "mod_corr_graph_name": my_id + mod_corr_graph_name,
                "input_id": "JCoup_" + Jcoup_type,
            }
        )

    jcoup_exp_file.close()
    jcoup_model_data_path = my_path + "/Jcoup_model.pickle"
    pickle.dump(type_dict, open(jcoup_model_data_path, "wb"))
    return jcuop_data
