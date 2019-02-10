import pickle
import math

import consensx.graph as graph

from consensx.csx_libs import methods as csx_func
from consensx.csx_libs import objects as csx_obj


# Equation and coefficients from:
# Wang & Bax (1996) JACS 118:2483-2494. Table 1, NMR + X-ray data
Jcoup_dict1 = {
    'A': {"3JHNCB": 3.39,  "3JHNHA": 6.98,  "3JHNC": 4.32, "3JHAC": 3.75},
    'B': {"3JHNCB": -0.94, "3JHNHA": -1.38, "3JHNC": 0.84, "3JHAC": 2.19},
    'C': {"3JHNCB": 0.07,  "3JHNHA": 1.72,  "3JHNC": 0.00, "3JHAC": 1.28},
    'THETA': {
        "3JHNCB": math.radians(60), "3JHNHA": math.radians(-60),
        "3JHNC":  math.radians(0),  "3JHAC":  math.radians(-60)}  # RAD!
}

# J. Am. Chem. Soc., Vol. 119, No. 27, 1997; Table 2 -> solution
Jcoup_dict2 = {
    'A': {"3JHNCB": 3.06,  "3JHNHA": 7.13, "3JHNC": 4.19, "3JHAC": 3.84},
    'B': {"3JHNCB": -0.74, "3JHNHA": 1.31, "3JHNC": 0.99, "3JHAC": 2.19},
    'C': {"3JHNCB": 0.10,  "3JHNHA": 1.56, "3JHNC": 0.03, "3JHAC": 1.20},
    'THETA': {
        "3JHNCB": math.radians(60),  "3JHNHA": math.radians(-60),
        "3JHNC":  math.radians(180), "3JHAC":  math.radians(120)}  # RAD!
}

# https://x86.cs.duke.edu/~brd/Teaching/Bio/asmb/Papers/NMR/nilges-jmr05.pdf
Jcoup_dict3 = {
    'A': {"3JHNCB": 3.26,  "3JHNHA": 7.13,  "3JHNC": 4.19, "3JHAC": 3.84},
    'B': {"3JHNCB": -0.87, "3JHNHA": -1.31, "3JHNC": 0.99, "3JHAC": 2.19},
    'C': {"3JHNCB": 0.10,  "3JHNHA": 1.56,  "3JHNC": 0.03, "3JHAC": 1.20},
    'THETA': {
        "3JHNCB": math.radians(60), "3JHNHA": math.radians(-60),
        "3JHNC":  math.radians(0),  "3JHAC":  math.radians(-60)}  # RAD!
}


def calc_dihedral_angles(pdb_model_data):
    """Calculates backbone diherdral angles
       note: all returned angle values are in radian"""

    JCoup_dicts = []

    for i in range(pdb_model_data.coordsets):
        pdb_model_data.atomgroup.setACSIndex(i)
        current_Resindex = 1
        prev_C, my_N, my_CA, my_C = None, None, None, None
        JCoup_dict = {}

        for atom in pdb_model_data.atomgroup:
            atom_res = atom.getResindex() + 1

            if atom_res != current_Resindex:

                if (prev_C is not None and my_N is not None and
                        my_CA is not None and my_C is not None):

                    NCA_vec = my_N - my_CA
                    CN_vec = prev_C - my_N
                    CCA_vec = my_C - my_CA

                    first_cross = csx_obj.Vec_3D.cross(CN_vec, NCA_vec)
                    second_cross = csx_obj.Vec_3D.cross(CCA_vec, NCA_vec)

                    angle = csx_obj.Vec_3D.dihedAngle(
                        first_cross, second_cross
                    )

                    # reference for setting sign of angle
                    reference = csx_obj.Vec_3D.cross(first_cross, second_cross)

                    r1 = reference.normalize()
                    r2 = NCA_vec.normalize()

                    if ((r1 - r2).magnitude() < r2.magnitude()):
                        angle *= -1

                    JCoup_dict[current_Resindex] = -1 * math.radians(angle)

                current_Resindex = atom_res
                prev_C = my_C
                my_N, my_CA, my_C = None, None, None

            if atom_res == current_Resindex:
                if atom.getName() == 'N':
                    my_N = csx_obj.Vec_3D(atom.getCoords())
                elif atom.getName() == 'CA':
                    my_CA = csx_obj.Vec_3D(atom.getCoords())
                elif atom.getName() == 'C':
                    my_C = csx_obj.Vec_3D(atom.getCoords())

        JCoup_dicts.append(JCoup_dict)

    return JCoup_dicts


def calc_jcoupling(param_set, calced, experimental, Jcoup_type):
    """Calculates J-coupling values from dihedral angles
       note: all angles must be in radian"""
    JCoup_calced = {}

    if param_set == 1:
        my_karplus = Jcoup_dict1
    elif param_set == 2:
        my_karplus = Jcoup_dict2
    elif param_set == 3:
        my_karplus = Jcoup_dict3

    A = my_karplus['A']
    B = my_karplus['B']
    C = my_karplus['C']
    THETA = my_karplus['THETA']

    for record in experimental:  # resnums
        J = 0

        for my_dict in calced:  # lists (with models as dicts)
            phi = my_dict[record.resnum]

            J += (A[Jcoup_type] * (math.cos(phi + THETA[Jcoup_type])) ** 2 +
                  B[Jcoup_type] * math.cos(phi + THETA[Jcoup_type]) +
                  C[Jcoup_type])

        JCoup_calced[record.resnum] = J / len(calced)

    model_data_list = []
    model_data_dict = {}

    for Jcoup_dict in calced:   # model
        for record in experimental:
            phi = Jcoup_dict[record.resnum]

            J = (
                A[Jcoup_type] * (math.cos(phi + THETA[Jcoup_type])) ** 2 +
                B[Jcoup_type] * math.cos(phi + THETA[Jcoup_type]) +
                C[Jcoup_type])

            model_data_dict[record.resnum] = J

        model_data_list.append(model_data_dict)
        model_data_dict = {}

    return JCoup_calced, model_data_list


def jcoupling(
        my_CSV_buffer, pdb_model_data, param_set, Jcoup_dict, my_PDB, my_path
        ):
    """Back calculate skalar coupling from given RDC lists and PDB models"""
    Jcuop_data = []
    type_dict = {}
    dihed_lists = calc_dihedral_angles(pdb_model_data)

    for Jcoup_type in sorted(list(Jcoup_dict.keys())):
        JCoup_calced, model_data = calc_jcoupling(
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
