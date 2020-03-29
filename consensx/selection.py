#!/usr/bin/env python3

import os
import pickle
import prody
from time import perf_counter

# own modules
import consensx.graph as graph
import consensx.calc as calc
from consensx.misc.natural_sort import natural_sort
from .models import CSX_upload


def get_pdb_models(path):
    pdb_models = []
    for file in os.listdir(path):
        if file.endswith(".pdb") and file.startswith("model_"):
            pdb_models.append(file)

    pdb_models = natural_sort(pdb_models)
    return pdb_models


class ChemshiftModelData:
    """Class for per model chemical shift data"""
    def __init__(self, type_dict):
        self.type_dict = type_dict

    def get_type_data(self, my_type):
        type_data = []

        for model in self.type_dict:
            type_data.append(model[my_type])

        return type_data


class DumpedData:
    def __init__(self):
        self.RDC_is_loaded = False
        self.RDC_lists = None
        self.RDC_model_data = None

        self.S2_is_loaded = False
        self.S2_dict = None
        self.S2_fit = False
        self.S2_fit_range = None

        self.PDB_is_loaded = False
        self.pdb_model_data = None

        self.Jcoup_is_loaded = False
        self.Jcoup_dict = None
        self.Jcoup_model_data = None

        self.ChemShift_is_loaded = False
        self.ChemShift_lists = None
        self.ChemShift_model_data = None

    def load_rdc_dump(self, path):
        rdc_lists_path = path + "/RDC_lists.pickle"
        self.RDC_lists = pickle.load(open(rdc_lists_path, 'rb'))
        model_data_path = path + "/RDC_model_data.pickle"
        self.RDC_model_data = pickle.load(open(model_data_path, 'rb'))
        self.RDC_is_loaded = True

    def load_s2_dump(self, path):
        s2_dict_path = path + "/S2_dict.pickle"
        s2_dump = pickle.load(open(s2_dict_path, 'rb'))
        self.S2_dict = s2_dump[0]
        self.S2_fit = s2_dump[1]
        self.S2_fit_range = s2_dump[2]
        self.S2_is_loaded = True

    def load_pdb_data(self, path):
        pdb_model_path = path + "/PDB_model.pickle"
        self.pdb_model_data = pickle.load(open(pdb_model_path, 'rb'))
        self.PDB_is_loaded = True

    def load_jcoup_data(self, path):
        jcoup_dict_path = path + "/Jcoup_dict.pickle"
        self.Jcoup_dict = pickle.load(open(jcoup_dict_path, 'rb'))
        jcoup_model_data = path + "/Jcoup_model.pickle"
        self.Jcoup_model_data = pickle.load(open(jcoup_model_data, 'rb'))
        self.Jcoup_is_loaded = True

    def load_chemshift_data(self, path):
        cs_lists_path = path + "/ChemShift_lists.pickle"
        self.ChemShift_lists = pickle.load(open(cs_lists_path, 'rb'))
        cs_model_data_path = path + "/ChemShift_model_data.pickle"
        model_data = pickle.load(open(cs_model_data_path, 'rb'))
        self.ChemShift_model_data = ChemshiftModelData(model_data)
        self.ChemShift_is_loaded = True


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


def averageS2_on(models, model_data, S2_dict, S2_type, fit, fit_range):
    """Returns a dictonary with the average S2 values for the given S2 type:
       averageS2[residue] = value"""
    my_models = []

    for model_num in models:
        model_data.atomgroup.setACSIndex(model_num)
        my_models.append(model_data.atomgroup[:])

    # csx_obj.PDB_model.is_fitted = False

    return calc.s2_values(
        model_data, models, S2_dict[S2_type], S2_type, fit, fit_range
    )


def get_two_random_from(models):
    from numpy import random
    random.shuffle(models)
    return models[0:2]


class Selection:
    def __init__(self, my_path, original_values, user_selection_json):
        self.max_size = None
        self.min_size = None
        self.overdrive = None
        self.measure = None

        self.RDC_lists = None
        self.RDC_model_data = None
        self.S2_dict = None
        self.pdb_data = None
        self.Jcoup_dict = None
        self.JCoup_modell_data = None
        self.ChemShifts = None
        self.ChemShift_model_data = None

        self.my_path = my_path
        self.original_values = original_values
        self.dumped_data = DumpedData()
        self.user_sel = []

        for key, value in user_selection_json.items():
            # set MEASURE for selection
            if key == "MEASURE":
                self.measure = value
                print("MEASURE is set to:", self.measure)

            if key == "MIN_SIZE":
                self.min_size = int(value)

            if key == "MAX_SIZE":
                self.max_size = int(value)

            if key == "OVERDRIVE":
                self.overdrive = int(value)

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
                elif my_type == "0CH":
                    my_type = "0_C_H"
                elif my_type == "0CN":
                    my_type = "0_C_N"
                elif my_type == "0HAC":
                    my_type = "0_HA_C"
                elif my_type == "0HAN":
                    my_type = "0_HA_N"

                my_weight = float(value)
                self.user_sel.append(["RDC", my_list, my_type, my_weight])

            # read S2 selections
            elif key.split('_')[0] == "S2":
                my_type = key.split('_')[1]
                my_weight = float(value)
                self.user_sel.append(["S2", my_type, my_weight])

            # read J-coupling selections
            elif key.split('_')[0] == "JCoup":
                my_type = key.split('_')[1]
                my_weight = float(value)
                self.user_sel.append(["JCoup", my_type, my_weight])

            # read chemical shifts selections
            elif key.split('_')[0] == "CS":
                my_type = key.split('_')[1]
                my_weight = float(value)
                self.user_sel.append(["ChemShift", my_type, my_weight])

    def run_selection(self):
        pdb_output_name = self.my_path + "/raw.pdb"

        if os.path.isfile(pdb_output_name):
            os.remove(pdb_output_name)

        if os.path.isfile(self.my_path + "/selected.pdb"):
            os.remove(self.my_path + "/selected.pdb")

        in_selection, iter_data, iter_all = self.selection_on()

        print("ITER ALL", iter_all)

        for key, val in iter_data.items():
            print("CALCED ", key, '{0:.3f}'.format(val))

        self.dumped_data.load_pdb_data(self.my_path)
        sel_ensemble = self.dumped_data.pdb_model_data.atomgroup.copy()

        for model_num in reversed(range(sel_ensemble.numCoordsets())):
            if model_num not in in_selection:
                sel_ensemble.delCoordset(model_num)

        num_coordsets = sel_ensemble.numCoordsets()

        print("NUM_COORDSETS: ", num_coordsets)

        prody.alignCoordsets(sel_ensemble.calpha)
        prody.writePDB(pdb_output_name, sel_ensemble)

        in_selection = [str(x+1) for x in sorted(in_selection)]
        dummy_pdb = open(pdb_output_name, 'r')
        output_pdb = open(self.my_path + "/selected.pdb", "w")

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

        calc_id = self.my_path.split('/')[-1]
        print('calcID', calc_id)
        db_entry = CSX_upload.objects.get(id_code=calc_id)
        print(db_entry.PDB_file)
        print("db_entry", db_entry)

        pca_image_names = graph.pca.create_pca_comparison(
            self.my_path, db_entry.PDB_file, in_selection
        )

        return num_coordsets, iter_data, pca_image_names

    def selection_on(self):
        pdb_models = get_pdb_models(self.my_path)

        print("user_sel: ", self.user_sel)

        for sel in self.user_sel:
            print(sel)
            if "RDC" in sel and not self.dumped_data.RDC_is_loaded:
                self.dumped_data.load_rdc_dump(self.my_path)
                print("RDC_lists assigned")
                # TODO do we need double names?
                self.RDC_lists = self.dumped_data.RDC_lists
                self.RDC_model_data = self.dumped_data.RDC_model_data

            if "S2" in sel and not self.dumped_data.S2_is_loaded:
                self.dumped_data.load_s2_dump(self.my_path)
                self.S2_dict = self.dumped_data.S2_dict
                fit = self.dumped_data.S2_fit
                fit_range = self.dumped_data.S2_fit_range

                self.dumped_data.load_pdb_data(self.my_path)
                self.pdb_data = self.dumped_data.pdb_model_data

            if "JCoup" in sel and not self.dumped_data.Jcoup_is_loaded:
                self.dumped_data.load_jcoup_data(self.my_path)
                self.Jcoup_dict = self.dumped_data.Jcoup_dict
                self.JCoup_modell_data = self.dumped_data.Jcoup_model_data

            if "ChemShift" in sel and not self.dumped_data.ChemShift_is_loaded:
                self.dumped_data.load_chemshift_data(self.my_path)
                self.ChemShifts = self.dumped_data.ChemShift_lists
                self.ChemShift_model_data = self.dumped_data.ChemShift_model_data

        in_selection = []

        print("STARTING WITH MODEL(S):", in_selection)

        first_run = True
        first_try = True
        above_best = 0
        iter_data = []
        divide_by = 0.0
        num_models = len(pdb_models)

        if self.measure == "correlation":
            prev_best = -2
        else:
            prev_best = 1000

        iter_count = 1

        t1_start = perf_counter()

        while True:
            model_scores = {}
            iter_scores = {}

            # iterate on all PDB models
            for num in range(num_models):
                # skip models already included in selection
                if num in in_selection:
                    continue

                divide_by = 0.0  # variable for storing weight sum
                pdb_sel = [num] + in_selection  # creating test ensemble

                for sel_data in self.user_sel:
                    if sel_data[0] == "RDC":
                        RDC_num = sel_data[1]
                        RDC_type = sel_data[2]
                        RDC_weight = sel_data[3]

                        my_data = self.RDC_model_data.rdc_data[RDC_num][RDC_type]
                        averageRDC = averageRDCs_on(pdb_sel, my_data)
                        my_RDC = self.RDC_lists[RDC_num - 1][RDC_type]
                        calced = None

                        if self.measure == "correlation":
                            calced = calc.correlation(averageRDC, my_RDC)
                        elif self.measure == "q-value":
                            calced = calc.q_value(averageRDC, my_RDC)
                        elif self.measure == "rmsd":
                            calced = calc.rmsd(averageRDC, my_RDC)

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
                            pdb_sel, self.pdb_data,
                            self.S2_dict, S2_type,
                            self.dumped_data.S2_fit, self.dumped_data.S2_fit_range
                        )
                        experimental = self.S2_dict[S2_type]
                        calced = None

                        if self.measure == "correlation":
                            calced = calc.correlation(averageS2, experimental)
                        elif self.measure == "q-value":
                            calced = calc.q_value(averageS2, experimental)
                        elif self.measure == "rmsd":
                            calced = calc.rmsd(averageS2, experimental)

                        if num in model_scores.keys():
                            model_scores[num] += calced * S2_weight
                        else:
                            model_scores[num] = calced * S2_weight

                        divide_by += S2_weight

                        iter_scores[sel_data[0] + '_' + str(sel_data[1])] = calced

                    elif sel_data[0] == "JCoup":
                        JCoup_type = sel_data[1]
                        JCoup_weight = sel_data[2]

                        my_type = self.JCoup_modell_data[JCoup_type]
                        averageJCoup = averageJCoup_on(pdb_sel, my_type)
                        my_JCoup = self.Jcoup_dict[JCoup_type]
                        calced = None

                        if self.measure == "correlation":
                            calced = calc.correlation(averageJCoup, my_JCoup)
                        elif self.measure == "q-value":
                            calced = calc.q_value(averageJCoup, my_JCoup)
                        elif self.measure == "rmsd":
                            calced = calc.rmsd(averageJCoup, my_JCoup)

                        if num in model_scores.keys():
                            model_scores[num] += calced * JCoup_weight
                        else:
                            model_scores[num] = calced * JCoup_weight

                        divide_by += JCoup_weight

                        iter_scores[sel_data[0] + '_' + str(sel_data[1])] = calced

                    elif sel_data[0] == "ChemShift":
                        ChemShift_type = sel_data[1]
                        ChemShift_weight = sel_data[2]

                        my_ChemShifts = self.ChemShifts[0][ChemShift_type]
                        my_type = self.ChemShift_model_data.get_type_data(
                            ChemShift_type
                        )

                        averageChemShift = averageChemShift_on(pdb_sel, my_type)
                        calced = None

                        if self.measure == "correlation":
                            calced = calc.correlation(
                                averageChemShift, my_ChemShifts
                            )
                        elif self.measure == "q-value":
                            calced = calc.q_value(
                                averageChemShift, my_ChemShifts
                            )
                        elif self.measure == "rmsd":
                            calced = calc.rmsd(
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

            if self.measure == "correlation":
                best_val = -2
            else:
                best_val = 1000

            for num in model_scores.keys():
                model_score = model_scores[num] / divide_by
                if self.measure == "correlation" and model_score > best_val:
                    best_val = model_score
                    best_num = num
                elif self.measure in ["q-value", "rmsd"] and model_score < best_val:
                    best_val = model_score
                    best_num = num

            if first_run:
                first_run = False
                if self.measure == "correlation":
                    best_val = -2
                else:
                    best_val = 1000

            print("ITERATION     #" + str(iter_count))
            print("prev best:    " + str(prev_best))
            print("current best: " + str(best_val))

            if self.max_size and len(in_selection) == self.max_size:
                print("size limit reached!")

                if self.overdrive:
                    if (
                            (self.measure == "correlation" and best_val > prev_best) or
                            (self.measure in ["q-value", "rmsd"] and best_val < prev_best)
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
                        t1_stop = perf_counter()
                        print("[selection] Selection in seconds:", t1_stop - t1_start)
                        return in_selection, iter_data[-above_best - 1], iter_data
                    else:
                        print(
                            "EXIT -> selection reached max desired size \
                            in overdrive"
                        )
                        t1_stop = perf_counter()
                        print("[selection] Selection in seconds:", t1_stop - t1_start)
                        return in_selection, iter_data[-above_best - 1], iter_data

            # if new selection results a higher score
            if (
                    (self.measure == "correlation" and best_val > prev_best) or
                    (self.measure in ["q-value", "rmsd"] and best_val < prev_best)
            ):

                # reset above the best threshold
                above_best = 0
                prev_best = best_val
                overdrive_best = -1
                print("CURRENT SEL:", in_selection)
                print("APPEND:", best_num)
                in_selection.append(best_num)

                # check if selection reached the desired maximal size (if any)
                if self.max_size and len(in_selection) - 1 == self.max_size:
                    print("size limit reached!")
                    # in_selection = [x+1 for x in in_selection]
                    in_selection.sort()
                    # print("numbered as in PDB file:\n", in_selection)
                    print("EXIT -> selection reached max desired size")
                    t1_stop = perf_counter()
                    print("[selection] Selection in seconds:", t1_stop - t1_start)
                    return in_selection, iter_data[-1], iter_data

            # if new selection results a lower score
            else:
                # check if overdrive is enabled
                if self.overdrive and self.overdrive > above_best:

                    # don't overdrive until minimum ensemble size reached
                    if self.min_size and len(in_selection) <= self.min_size:
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
                        t1_stop = perf_counter()
                        print("[selection] Selection in seconds:", t1_stop - t1_start)
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

                    if self.measure == "correlation" and overdrive_best > prev_best:
                        prev_best = overdrive_best
                        above_best = 0
                    elif (self.measure in ["q-value", "rmsd"] and
                          overdrive_best < prev_best):
                        prev_best = overdrive_best
                        above_best = 0

                    if self.overdrive == above_best:
                        if self.measure == "correlation" and overdrive_best < prev_best:
                            for _ in range(above_best + 1):
                                # print(in_selection)
                                print("POP", in_selection[-1])
                                del in_selection[-1]

                                print("CURRENT SEL:", in_selection)

                            print("EXIT -> selection reached max override value")
                            t1_stop = perf_counter()
                            print("[selection] Selection in seconds:", t1_stop - t1_start)
                            return (
                                in_selection, iter_data[-above_best - 1], iter_data
                            )

                        if (self.measure in ["q-value", "rmsd"] and
                                overdrive_best > prev_best):

                            for _ in range(above_best + 1):
                                print("POP", in_selection[-1])
                                del in_selection[-1]

                                print("CURRENT SEL:", in_selection)

                            print("EXIT -> selection reached max override value")
                            t1_stop = perf_counter()
                            print("[selection] Selection in seconds:", t1_stop - t1_start)
                            return (
                                in_selection, iter_data[-above_best - 1], iter_data
                            )

                    continue

                # check if selection reached the desired minimal size (if any)
                if self.min_size and len(in_selection) <= self.min_size:
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
                t1_stop = perf_counter()
                print("[selection] Selection in seconds:", t1_stop - t1_start)
                return in_selection, iter_data[-1], iter_data

            iter_count += 1

