#!/usr/bin/env python3

"""
This is the Django webapp implemetation of
              _____      _   _  _____ ______          __   __
             / ____|    | \ | |/ ____|  ____|         \ \ / /
            | |     ___ |  \| | (___ | |__   _ __  ___ \ V /
            | |    / _ \| . ` |\___ \|  __| | '_ \/ __| > <
            | |___| (_) | |\  |____) | |____| | | \__ \/ . \
             \_____\___/|_| \_|_____/|______|_| |_|___/_/ \_\

     Compliance of NMR-derived Structural Ensembles with experimental data

Authors: Zoltán Gáspári, Dániel Dudola
"""

# standard modules
import os
import pickle

# own modules
import consensx.csx_libs.calc as csx_calc
import consensx.csx_libs.methods as csx_func
import consensx.csx_libs.objects as csx_obj

# Django server
from django.shortcuts import render
from django.http import HttpResponse
from .models import CSX_upload, CSX_calculation

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run_calculation(request, calc_id):
    abspath = os.path.abspath(__file__)
    dirname = os.path.dirname(abspath)
    csx_func.check_3rd_party(dirname)

    my_id = calc_id
    print("ID IS:", my_id)

    from django.conf import settings
    my_path = os.path.join(settings.MEDIA_ROOT, my_id) + '/'
    DB_entry = CSX_upload.objects.get(id_code=calc_id)

    # ----------------  Setting up working directory and files  ------------- #
    my_CSV_buffer = csx_obj.CSV_buffer(my_path)

    my_PDB = my_path + DB_entry.PDB_file

    csx_func.pdb_cleaner(my_path, my_PDB, my_CSV_buffer)
    model_count = csx_func.pdb_splitter(my_path, my_PDB)

    if not csx_func.get_model_list(my_PDB, my_path, model_count):
        return render(request, "consensx/home.html", {
            "error": "DISCARDED MODELS FOUND, CHECK IF ALL MODELS HAVE THE SAME\
            NUMBER OF ATOMS"
        })

    pdb_models = []                                # list of models (PDB)
    for file in os.listdir(my_path):
        if file.startswith("model_") and file.endswith(".pdb"):
            pdb_models.append(file)

    pdb_models = csx_func.natural_sort(pdb_models)
    data_found = False

    # ---------------------  Read  and parse NOE file   --------------------- #
    if DB_entry.NOE_file:
        # empty class variables
        csx_obj.Restraint_Record.all_restraints = []
        csx_obj.Restraint_Record.resolved_restraints = []

        NOE_name = DB_entry.NOE_file
        my_NOE = my_path + DB_entry.NOE_file
        saveShifts = csx_func.getNOE(my_NOE)
        NOE_n = saveShifts[-1][0] + " distance restraints found"
        NOE_violations = str(csx_calc.calcNOEviolations(
            my_NOE, saveShifts, my_path, DB_entry.r3average
        ))
        PRIDE_data = csx_calc.calcNMR_Pride(pdb_models, my_path)

        NOE_PRIDE_data = {
            "NOE_violations": NOE_violations,
            "NOE_hist": my_id + "/NOE_hist.svg",
            "best_score": PRIDE_data[0],
            "worst_score": PRIDE_data[1],
            "average_score": '{0:.3f}'.format(PRIDE_data[2]),
            "deviation": '{0:.3f}'.format(PRIDE_data[3]),
            "PRIDE_hist": my_id + "/PRIDE-NMR_score.svg"
        }
        data_found = True
    else:
        NOE_name = "[NOT PRESENT]"
        NOE_n = ""
        NOE_PRIDE_data = None

    # ---------------------  Read  and parse STR file   --------------------- #
    if DB_entry.STR_file:
        STR_name = DB_entry.STR_file
        my_STR = my_path + DB_entry.STR_file

        try:
            parsed = csx_func.parseSTR(my_STR)
        except Exception as e:
            print("EXCEPTION", e)
            return HttpResponse(e)

        # ------------------------  RDC calculation  ------------------------ #
        RDC_lists = csx_func.get_RDC_lists(parsed.value)
        RDC_lists_path = my_path + "/RDC_lists.pickle"
        pickle.dump(RDC_lists, open(RDC_lists_path, "wb"))

        if RDC_lists:
            SVD_enabled = DB_entry.svd_enable
            lc_model = DB_entry.rdc_lc
            RDC_data = csx_calc.calcRDC(my_CSV_buffer, RDC_lists, pdb_models,
                                        my_path, SVD_enabled, lc_model)
            data_found = True
        else:
            RDC_data = None

        # ----------------------------  S2 calc  ---------------------------- #
        S2_dict = csx_func.parseS2_STR(parsed.value)
        S2_dump = [S2_dict, DB_entry.superimpose, DB_entry.fit_range]
        S2_dict_path = my_path + "/S2_dict.pickle"
        pickle.dump(S2_dump, open(S2_dict_path, "wb"))

        if S2_dict:
            S2_data = csx_calc.calcS2(
                           my_CSV_buffer, S2_dict, my_path,
                           fit=DB_entry.superimpose,
                           fit_range=DB_entry.fit_range)

            data_found = True
        else:
            S2_data = None

        S2_sidechain = csx_func.parse_sidechain_S2_STR(parsed.value)

        if S2_sidechain:
            S2_sc_data = csx_calc.calcS2_sidechain(
                my_CSV_buffer, S2_sidechain, my_path, fit=DB_entry.superimpose
            )

            if "error" in S2_sc_data.keys():
                return render(request, "consensx/home.html", {
                    "error": S2_sidechain["error"]
                })

            data_found = True
        else:
            S2_sc_data = None

        # ------------------------  J-coupling calc  ------------------------ #
        Jcoup_dict = csx_func.parseJcoup_STR(parsed.value)
        Jcoup_dict_path = my_path + "/Jcoup_dict.pickle"
        pickle.dump(Jcoup_dict, open(Jcoup_dict_path, "wb"))

        if Jcoup_dict:
            Jcoup_data = csx_calc.calcJCouplings(
                my_CSV_buffer, DB_entry.karplus, Jcoup_dict, my_PDB, my_path
            )
            data_found = True
        else:
            Jcoup_data = None

        # -----------------------  Chemical shift calc  --------------------- #
        ChemShift_lists = csx_func.parseChemShift_STR(parsed.value)
        ChemShift_lists_path = my_path + "/ChemShift_lists.pickle"
        pickle.dump(ChemShift_lists, open(ChemShift_lists_path, "wb"))

        if ChemShift_lists:
            chemshift_data = csx_calc.calcChemShifts(
                my_CSV_buffer, ChemShift_lists, pdb_models, my_path
            )
            data_found = True
        else:
            chemshift_data = None
    else:
        STR_name = "[NOT PRESENT]"
        RDC_data = None
        S2_data = None
        Jcoup_data = None
        chemshift_data = None

    my_CSV_buffer.writeCSV()

    if data_found:
        print(csx_obj.CalcPickle.data)
        calced_values = my_path + "/calced_values.p"
        pickle.dump(csx_obj.CalcPickle.data, open(calced_values, "wb"))
        rendered_page = render(request, "consensx/calculation.html", {
            "my_id": my_id,
            "my_PDB": DB_entry.PDB_file,
            "n_model": model_count,
            "my_NOE": NOE_name,
            "n_NOE": NOE_n,
            "my_STR": STR_name,
            "NOE_PRIDE_data": NOE_PRIDE_data,
            "RDC_data": RDC_data,
            "S2_data": S2_data,
            "S2_sc_data": S2_sc_data,
            "Jcoup_data": Jcoup_data,
            "chemshift_data": chemshift_data,
            "SVD_calc": DB_entry.svd_enable
        })

        print("RENDERED PAGE ---------------- START")
        post_data = CSX_calculation(
            html_content=rendered_page.content,
            id_code=my_id
        )
        post_data.save()

        print("RENDERED PAGE ------------------ END")
        return rendered_page
    else:
        print("NO DATA FOUND IN STAR-NMR FILE")
        return render(request, "consensx/home.html", {
            "error": "NO DATA FOUND IN STAR-NMR FILE"
        })
