#!/usr/bin/env python3

import os
import pickle
import prody

# own modules
import consensx.csx_libs.methods as csx_func
import consensx.csx_libs.objects as csx_obj
import consensx.csx_libs.pca as csx_pca
from .models import CSX_upload


def getPDBModels(path):
    pdb_models = []
    for file in os.listdir(path):
        if file.endswith(".pdb") and file.startswith("model_"):
            pdb_models.append(file)

    pdb_models = csx_func.natural_sort(pdb_models)
    return pdb_models


class DumpedData():
    RDC_isloaded = False
    RDC_lists = None
    RDC_model_data = None

    S2_isloaded = False
    S2_dict = None
    S2_fit = False
    S2_fit_range = None

    PDB_isloaded = False
    PDB_model_data = None

    Jcoup_isloaded = False
    Jcoup_dict = None
    Jcoup_model_data = None

    ChemShift_isloaded = False
    ChemShift_lists = None
    ChemShift_model_data = None

    @staticmethod
    def loadRDCDump(path):
        RDC_lists_path = path + "/RDC_lists.pickle"
        DumpedData.RDC_lists = pickle.load(open(RDC_lists_path, 'rb'))
        model_data_path = path + "/RDC_model_data.pickle"
        DumpedData.RDC_model_data = pickle.load(open(model_data_path, 'rb'))
        RDC_obj_attr_path = path + "/RDC_obj_attr.pickle"
        RDC_obj_attr = pickle.load(open(RDC_obj_attr_path, 'rb'))
        DumpedData.RDC_model_data.RDC_data = RDC_obj_attr
        DumpedData.RDC_isloaded = True

    @staticmethod
    def loadS2Dump(path):
        S2_dict_path = path + "/S2_dict.pickle"
        S2_dump = pickle.load(open(S2_dict_path, 'rb'))
        DumpedData.S2_dict = S2_dump[0]
        DumpedData.S2_fit = S2_dump[1]
        DumpedData.S2_fit_range = S2_dump[2]
        DumpedData.S2_isloaded = True

    @staticmethod
    def loadPDBData(path):
        PDB_model_path = path + "/PDB_model.pickle"
        DumpedData.PDB_model_data = pickle.load(open(PDB_model_path, 'rb'))
        DumpedData.PDB_isloaded = True

    @staticmethod
    def loadJCoupData(path):
        Jcoup_dict_path = path + "/Jcoup_dict.pickle"
        DumpedData.Jcoup_dict = pickle.load(open(Jcoup_dict_path, 'rb'))
        Jcoup_model_data = path + "/Jcoup_model.pickle"
        DumpedData.Jcoup_model_data = pickle.load(open(Jcoup_model_data, 'rb'))
        DumpedData.Jcoup_isloaded = True

    @staticmethod
    def loadChemShiftData(path):
        CS_lists_path = path + "/ChemShift_lists.pickle"
        DumpedData.ChemShift_lists = pickle.load(open(CS_lists_path, 'rb'))
        CS_model_data_path = path + "/ChemShift_model_data.pickle"
        model_data = pickle.load(open(CS_model_data_path, 'rb'))
        ChemShift_model_data = csx_obj.ChemShift_modell_data
        ChemShift_model_data.type_dict = model_data
        DumpedData.ChemShift_model_data = ChemShift_model_data
        DumpedData.ChemShift_isloaded = True


def getUserSel(sel_dict):
    """Sets up user selection list from the selection dictonary"""

    user_sel = []

    for key, value in sel_dict.items():

        print(key, value)

        # set MEASURE for selection
        if key == "MEASURE":
            global measure
            measure = value
            print("MEASURE is set to:", measure)

        if key == "MIN_SIZE":
            global min_size
            min_size = int(value)

        if key == "MAX_SIZE":
            global max_size
            max_size = int(value)

        if key == "OVERDRIVE":
            global overdrive
            overdrive = int(value)

        # read RDC selections
        if key.split('_')[0] == "RDC":
            my_list = int(key.split('_')[1])
            my_type = key.split('_')[2]

            # TODO - we could avoid this
            if my_type == "0CAC":
                my_type = "0_CA_C"
            elif my_type == "0HACA":
                my_type = "0_HA_CA"
            elif my_type == "0NH":
                my_type = "0_N_H"
            elif my_type == "1NC":
                my_type = "1_N_C"
            elif my_type == "0CCA":
                my_type = "0_C_CA"
            elif my_type == "0CAHA":
                my_type = "0_CA_HA"
            elif my_type == "0CACB":
                my_type = "0_CA_CB"

            my_weight = float(value)
            user_sel.append(["RDC", my_list, my_type, my_weight])

        # read S2 selections
        elif key.split('_')[0] == "S2":
            my_type = key.split('_')[1]
            my_weight = float(value)
            user_sel.append(["S2", my_type, my_weight])

        # read J-coupling selections
        elif key.split('_')[0] == "JCoup":
            my_type = key.split('_')[1]
            my_weight = float(value)
            user_sel.append(["JCoup", my_type, my_weight])

        # read chemical shifts selections
        elif key.split('_')[0] == "CS":
            my_type = key.split('_')[1]
            my_weight = float(value)
            user_sel.append(["ChemShift", my_type, my_weight])

    return user_sel


def averageRDCs_on(models, my_data):
    """Returns a dictonary with the average RDCs for the given RDC type:
       averageRDC[residue] = value"""

    averageRDC = {}

    for model_num, model in enumerate(my_data):
        if model_num not in models:
            continue

        for resnum in model:
            if resnum in averageRDC.keys():
                averageRDC[resnum] += model[resnum]
            else:
                averageRDC[resnum] = model[resnum]

    for resnum in list(averageRDC.keys()):
        averageRDC[resnum] /= len(models)

    return averageRDC


def averageJCoup_on(models, my_data):
    """Returns a dictonary with the average J-Couplings for the given type:
       averageJCoup[residue] = value"""

    averageJCoup = {}

    for model_num, model in enumerate(my_data):
        if model_num not in models:
            continue

        for resnum in model:
            if resnum in averageJCoup.keys():
                averageJCoup[resnum] += model[resnum]
            else:
                averageJCoup[resnum] = model[resnum]

    for resnum in list(averageJCoup.keys()):
        averageJCoup[resnum] /= len(models)

    return averageJCoup


def averageChemShift_on(models, my_data):
    """Returns a dictonary with the average chemical shifts for the given type:
       averageChemShift[residue] = value"""

    averageChemShift = {}

    for model_num, model in enumerate(my_data):
        if model_num not in models:
            continue

        for resnum in model:
            if resnum in averageChemShift.keys():
                averageChemShift[resnum] += model[resnum]
            else:
                averageChemShift[resnum] = model[resnum]

    for resnum in list(averageChemShift.keys()):
        averageChemShift[resnum] /= len(models)

    return averageChemShift


def averageS2_on(models, PDB_data, S2_dict, S2_type, fit, fit_range):
    """Returns a dictonary with the average S2 values for the given S2 type:
       averageS2[residue] = value"""

    model_data = PDB_data
    my_models = []

    for model_num in models:
        model_data.atomgroup.setACSIndex(model_num)
        my_models.append(model_data.atomgroup[:])

    csx_obj.PDB_model.is_fitted = False

    my_S2 = csx_func.calcS2(model_data, models, S2_dict[S2_type], S2_type,
                            fit, fit_range)

    return my_S2

def get_two_random_from(models):
    from numpy import random
    random.shuffle(models)
    return models[0:2]


def selection_on(my_path, measure, user_sel,
                 min_size=None, max_size=None, overdrive=None):

    pdb_models = getPDBModels(my_path)

    print("user_sel: ", user_sel)

    for sel in user_sel:
        print(sel)
        if "RDC" in sel and not DumpedData.RDC_isloaded:
            DumpedData.loadRDCDump(my_path)
            print("RDC_lists assigned")
            RDC_lists = DumpedData.RDC_lists
            RDC_model_data = DumpedData.RDC_model_data

        if "S2" in sel and not DumpedData.S2_isloaded:
            DumpedData.loadS2Dump(my_path)
            S2_dict = DumpedData.S2_dict
            fit = DumpedData.S2_fit
            fit_range = DumpedData.S2_fit_range

            DumpedData.loadPDBData(my_path)
            PDB_data = DumpedData.PDB_model_data

        if "JCoup" in sel and not DumpedData.Jcoup_isloaded:
            DumpedData.loadJCoupData(my_path)
            Jcoup_dict = DumpedData.Jcoup_dict
            JCoup_modell_data = DumpedData.Jcoup_model_data

        if "ChemShift" in sel and not DumpedData.ChemShift_isloaded:
            DumpedData.loadChemShiftData(my_path)
            ChemShifts = DumpedData.ChemShift_lists
            ChemShift_model_data = DumpedData.ChemShift_model_data

    in_selection = get_two_random_from(list(range(0, len(pdb_models))))

    print("STARTING WITH MODEL(S):", in_selection)

    first_run = True
    first_try = True
    above_best = 0
    iter_data = []

    if measure == "correlation":
        prev_best = -2
    else:
        prev_best = 1000

    while True:
        model_scores = {}
        iter_scores = {}

        # iterate on all PDB models
        for num in range(0, len(pdb_models)):

            # skip models already included in selection
            if num in in_selection:
                continue

            divide_by = 0.0                 # variable for storing weight sum
            pdb_sel = [num] + in_selection  # creating test ensemble

            for sel_data in user_sel:
                if sel_data[0] == "RDC":
                    RDC_num = sel_data[1]
                    RDC_type = sel_data[2]
                    RDC_weight = sel_data[3]

                    my_data = RDC_model_data.RDC_data[RDC_num][RDC_type]
                    averageRDC = averageRDCs_on(pdb_sel, my_data)
                    my_RDC = RDC_lists[RDC_num - 1][RDC_type]

                    if measure == "correlation":
                        calced = csx_func.calcCorrel(averageRDC, my_RDC)
                    elif measure == "q-value":
                        calced = csx_func.calcQValue(averageRDC, my_RDC)
                    elif measure == "rmsd":
                        calced = csx_func.calcRMSD(averageRDC, my_RDC)

                    if num in model_scores.keys():
                        model_scores[num] += calced * RDC_weight
                    else:
                        model_scores[num] = calced * RDC_weight

                    divide_by += RDC_weight

                    my_type = "".join(sel_data[2].split('_'))
                    my_key = (
                        sel_data[0] + '_' + str(sel_data[1]) + '_' + my_type
                    )
                    iter_scores[my_key] = calced

                elif sel_data[0] == "S2":
                    S2_type = sel_data[1]
                    S2_weight = sel_data[2]

                    averageS2 = averageS2_on(
                        pdb_sel, PDB_data,
                        S2_dict, S2_type,
                        fit, fit_range
                    )
                    experimental = S2_dict[S2_type]

                    if measure == "correlation":
                        calced = csx_func.calcCorrel(averageS2, experimental)
                    elif measure == "q-value":
                        calced = csx_func.calcQValue(averageS2, experimental)
                    elif measure == "rmsd":
                        calced = csx_func.calcRMSD(averageS2, experimental)

                    if num in model_scores.keys():
                        model_scores[num] += calced * S2_weight
                    else:
                        model_scores[num] = calced * S2_weight

                    divide_by += S2_weight

                    iter_scores[sel_data[0] + '_' + str(sel_data[1])] = calced

                elif sel_data[0] == "JCoup":
                    JCoup_type = sel_data[1]
                    JCoup_weight = sel_data[2]

                    my_type = JCoup_modell_data[JCoup_type]
                    averageJCoup = averageJCoup_on(pdb_sel, my_type)
                    my_JCoup = Jcoup_dict[JCoup_type]

                    if measure == "correlation":
                        calced = csx_func.calcCorrel(averageJCoup, my_JCoup)
                    elif measure == "q-value":
                        calced = csx_func.calcQValue(averageJCoup, my_JCoup)
                    elif measure == "rmsd":
                        calced = csx_func.calcRMSD(averageJCoup, my_JCoup)

                    if num in model_scores.keys():
                        model_scores[num] += calced * JCoup_weight
                    else:
                        model_scores[num] = calced * JCoup_weight

                    divide_by += JCoup_weight

                    iter_scores[sel_data[0] + '_' + str(sel_data[1])] = calced

                elif sel_data[0] == "ChemShift":
                    ChemShift_type = sel_data[1]
                    ChemShift_weight = sel_data[2]

                    my_ChemShifts = ChemShifts[0][ChemShift_type]
                    my_type = ChemShift_model_data.get_type_data(
                        ChemShift_type
                    )

                    averageChemShift = averageChemShift_on(pdb_sel, my_type)

                    if measure == "correlation":
                        calced = csx_func.calcCorrel(
                            averageChemShift, my_ChemShifts
                        )
                    elif measure == "q-value":
                        calced = csx_func.calcQValue(
                            averageChemShift, my_ChemShifts
                        )
                    elif measure == "rmsd":
                        calced = csx_func.calcRMSD(
                            averageChemShift, my_ChemShifts
                        )

                    if num in model_scores.keys():
                        model_scores[num] += calced * ChemShift_weight
                    else:
                        model_scores[num] = calced * ChemShift_weight

                    divide_by += ChemShift_weight

                    iter_scores["CS_" + sel_data[1]] = calced

        iter_data.append(iter_scores)

        best_num = -1

        if measure == "correlation":
            best_val = -2
        else:
            best_val = 1000

        for num in model_scores.keys():
            model_score = model_scores[num] / divide_by
            if measure == "correlation" and model_score > best_val:
                best_val = model_score
                best_num = num
            elif (measure in ["q-value", "rmsd"] and model_score < best_val):
                best_val = model_score
                best_num = num

        if first_run:
            first_run = False
            if measure == "correlation":
                best_val = -2
            else:
                best_val = 1000

        print("\nNEW ITERATION\n")
        print("prev best:    " + str(prev_best))
        print("current best: " + str(best_val))

        if max_size and len(in_selection) == max_size:
            print("size limit reached!")

            if overdrive:
                if (
                    (measure == "correlation" and best_val > prev_best) or
                    (measure in ["q-value", "rmsd"] and best_val < prev_best)
                ):

                    above_best = 0

                print("overdrive is:", above_best)
                # remove overdrive models
                for _ in range(above_best):
                    print("POP", in_selection[-1])
                    del in_selection[-1]

                    print("CURRENT SEL:", in_selection)

                in_selection.sort()
                if above_best == 0:
                    print(
                        "EXIT -> selection reached max desired size \
                        NOT in overdrive"
                    )
                    return in_selection, iter_data[-above_best - 1], iter_data
                else:
                    print(
                        "EXIT -> selection reached max desired size \
                        in overdrive"
                    )
                    return in_selection, iter_data[-above_best - 1], iter_data

        # if new selection results a higher score
        if (
            (measure == "correlation" and best_val > prev_best) or
            (measure in ["q-value", "rmsd"] and best_val < prev_best)
        ):

            # reset above the best threshold
            above_best = 0
            prev_best = best_val
            overdrive_best = -1
            print("CURRENT SEL:", in_selection)
            print("APPEND:", best_num)
            in_selection.append(best_num)

            # check if selection reached the desired maximal size (if any)
            if max_size and len(in_selection) - 1 == max_size:
                print("size limit reached!")
                # in_selection = [x+1 for x in in_selection]
                in_selection.sort()
                # print("numbered as in PDB file:\n", in_selection)
                print("EXIT -> selection reached max desired size")
                return in_selection, iter_data[-1], iter_data

        # if new selection results a lower score
        else:
            # check if overdrive is enabled
            if overdrive and overdrive > above_best:

                # don't overdrive until minimum ensemble size reached
                if min_size and len(in_selection) < min_size:
                    prev_best = best_val
                    in_selection.append(best_num)
                    continue

                # stop iteration if size of the original ensemble reached
                if len(in_selection) == len(pdb_models):
                    for _ in range(above_best + 1):
                        # print(in_selection)
                        print("POP", in_selection[-1])
                        del in_selection[-1]

                        print("CURRENT SEL:", in_selection)

                    print(
                        "EXIT -> selection reached original size in overdrive"
                    )
                    del in_selection[-1]
                    return in_selection, iter_data[-above_best - 1], iter_data

                above_best += 1
                print("\x1b[31mwe are in overdrive with \x1b[0m" +
                      str(above_best))
                overdrive_best = best_val
                print("overdrive_best: " + str(overdrive_best))
                print("prev_best: " + str(prev_best))
                print("CURRENT SEL:", in_selection)
                print("APPEND:", best_num)
                in_selection.append(best_num)

                if measure == "correlation" and overdrive_best > prev_best:
                    prev_best = overdrive_best
                    above_best = 0
                elif (measure in ["q-value", "rmsd"] and
                      overdrive_best < prev_best):
                    prev_best = overdrive_best
                    above_best = 0

                if overdrive == above_best:
                    if measure == "correlation" and overdrive_best < prev_best:
                        for _ in range(above_best + 1):
                            # print(in_selection)
                            print("POP", in_selection[-1])
                            del in_selection[-1]

                            print("CURRENT SEL:", in_selection)

                        print("EXIT -> selection reached max override value")
                        return (
                            in_selection, iter_data[-above_best - 1], iter_data
                        )

                    if (measure in ["q-value", "rmsd"] and
                            overdrive_best > prev_best):

                        for _ in range(above_best + 1):
                            print("POP", in_selection[-1])
                            del in_selection[-1]

                            print("CURRENT SEL:", in_selection)

                        print("EXIT -> selection reached max override value")
                        return (
                            in_selection, iter_data[-above_best - 1], iter_data
                        )

                continue

            # check if selection reached the desired minimal size (if any)
            if min_size and len(in_selection) < min_size:
                print("we are over the peak!")
                prev_best = best_val
                in_selection.append(best_num)
                continue

            if first_try:
                first_try = False
                continue

            if len(in_selection) > 1:
                del in_selection[-1]

            print("EXIT -> selection got a worse score, no override")
            return in_selection, iter_data[-1], iter_data


def run_selection(my_path, original_values, user_selection_JSON):
    DumpedData.RDC_isloaded = False
    DumpedData.S2_isloaded = False
    DumpedData.PDB_isloaded = False
    DumpedData.Jcoup_isloaded = False
    DumpedData.ChemShift_isloaded = False

    working_dir = my_path

    user_sel = getUserSel(user_selection_JSON)

    pdb_output_name = my_path + "/raw.pdb"

    if os.path.isfile(pdb_output_name):
        os.remove(pdb_output_name)

    if os.path.isfile(my_path + "/selected.pdb"):
        os.remove(my_path + "/selected.pdb")

    global max_size
    global min_size
    global overdrive

    if "min_size" not in globals():
        min_size = None

    if "max_size" not in globals():
        max_size = None

    print("max_size -> ", max_size)

    if "overdrive" not in globals():
        overdrive = None

    in_selection, iter_data, iter_all = selection_on(
        working_dir, measure, user_sel,
        min_size=min_size, max_size=max_size, overdrive=overdrive
    )

    print("ITER ALL", iter_all)

    for key, val in iter_data.items():
        print("CALCED ", key, '{0:.3f}'.format(val))

    DumpedData.loadPDBData(my_path)
    PDB_data = DumpedData.PDB_model_data
    sel_ensemble = PDB_data.atomgroup.copy()

    for model_num in reversed(range(sel_ensemble.numCoordsets())):
        if model_num not in in_selection:
            sel_ensemble.delCoordset(model_num)

    num_coordsets = sel_ensemble.numCoordsets()

    print("NUM_COORDSETS: ", num_coordsets)

    prody.alignCoordsets(sel_ensemble.calpha)
    prody.writePDB(pdb_output_name, sel_ensemble)

    in_selection = [str(x+1) for x in sorted(in_selection)]
    dummy_pdb = open(pdb_output_name, 'r')
    output_pdb = open(my_path + "/selected.pdb", "w")

    for line in dummy_pdb:
        output_pdb.write(line)
        if 'REMARK' in line:
            model_line = "REMARK ORIGINAL MODELS: "
            for model_num in in_selection:
                if len(model_line) < 76:
                    model_line += model_num + " "
                else:
                    output_pdb.write(model_line + "\n")
                    model_line = "REMARK ORIGINAL MODELS: "

            output_pdb.write(model_line + "\n")

    calc_id = my_path.split('/')[-1]
    print('calcID', calc_id)
    DB_entry = CSX_upload.objects.get(id_code=calc_id)
    print(DB_entry.PDB_file)
    print("DB_entry", DB_entry)

    pca_image_names = csx_pca.create_PCA_comparison(
        my_path, DB_entry.PDB_file, in_selection
    )

    return num_coordsets, iter_data, pca_image_names
