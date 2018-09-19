#!/usr/bin/python3

"""
This is the Django webapp implementation of
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
    db_entry = CSX_upload.objects.get(id_code=calc_id)

    # ----------------  Setting up working directory and files  ------------- #
    my_csv_buffer = csx_obj.CSV_buffer(my_path)

    my_pdb = my_path + db_entry.PDB_file

    csx_func.pdb_cleaner(my_path, my_pdb, my_csv_buffer)
    model_count = csx_func.pdb_splitter(my_path, my_pdb)

    if not csx_func.get_model_list(my_pdb, my_path, model_count):
        return render(request, "consensx/home.html", {
            "error": "DISCARDED MODELS FOUND, CHECK IF ALL MODELS HAVE THE SAME\
            NUMBER OF ATOMS"
        })

    pdb_models = []
    for file in os.listdir(my_path):
        if file.startswith("model_") and file.endswith(".pdb"):
            pdb_models.append(file)

    pdb_models = csx_func.natural_sort(pdb_models)
    data_found = False

    # ---------------------  Read  and parse NOE file   --------------------- #
    if db_entry.NOE_file:
        # empty class variables
        csx_obj.Restraint_Record.all_restraints = []
        csx_obj.Restraint_Record.resolved_restraints = []

        noe_name = db_entry.NOE_file
        my_noe = my_path + db_entry.NOE_file
        save_shifts = csx_func.getNOE(my_noe)
        noe_n = save_shifts[-1][0] + " distance restraints found"
        noe_violations = csx_calc.calcNOEviolations(
            my_pdb, save_shifts, my_path, db_entry.r3average
        )
        pride_data = csx_calc.calcNMR_Pride(pdb_models, my_path)

        noe_pride_data = {
            "noe_violations": str(noe_violations),
            "NOE_hist": my_id + "/NOE_hist.svg",
            "best_score": pride_data[0],
            "worst_score": pride_data[1],
            "average_score": '{0:.3f}'.format(pride_data[2]),
            "deviation": '{0:.3f}'.format(pride_data[3]),
            "PRIDE_hist": my_id + "/PRIDE-NMR_score.svg"
        }
        data_found = True
    else:
        noe_name = "[NOT PRESENT]"
        noe_n = ""
        noe_pride_data = None

    # ---------------------  Read  and parse STR file   --------------------- #
    if db_entry.STR_file:
        str_name = db_entry.STR_file
        my_str = my_path + db_entry.STR_file

        try:
            parsed = csx_func.parseSTR(my_str)
        except Exception as e:
            print("EXCEPTION", e)
            return HttpResponse(e)

        # ------------------------  RDC calculation  ------------------------ #
        rdc_lists = csx_func.get_RDC_lists(parsed.value)
        rdc_lists_path = my_path + "/RDC_lists.pickle"
        pickle.dump(rdc_lists, open(rdc_lists_path, "wb"))

        if rdc_lists:
            svd_enabled = db_entry.svd_enable
            lc_model = db_entry.rdc_lc
            rdc_data = csx_calc.calcRDC(my_csv_buffer, rdc_lists, pdb_models,
                                        my_path, svd_enabled, lc_model)
            data_found = True
        else:
            rdc_data = None

        # ----------------------------  S2 calc  ---------------------------- #
        s2_dict = csx_func.parseS2_STR(parsed.value)
        s2_dump = [s2_dict, db_entry.superimpose, db_entry.fit_range]
        s2_dict_path = my_path + "/S2_dict.pickle"
        pickle.dump(s2_dump, open(s2_dict_path, "wb"))

        if s2_dict:
            s2_data = csx_calc.calcS2(
                           my_csv_buffer, s2_dict, my_path,
                           fit=db_entry.superimpose,
                           fit_range=db_entry.fit_range)

            data_found = True
        else:
            s2_data = None

        s2_sidechain = csx_func.parse_sidechain_S2_STR(parsed.value)

        if s2_sidechain:
            s2_sc_data = csx_calc.calcS2_sidechain(
                my_csv_buffer, s2_sidechain, my_path, fit=db_entry.superimpose
            )

            if "error" in s2_sc_data.keys():
                return render(request, "consensx/home.html", {
                    "error": s2_sidechain["error"]
                })

            data_found = True
        else:
            s2_sc_data = None

        # ------------------------  J-coupling calc  ------------------------ #
        Jcoup_dict = csx_func.parseJcoup_STR(parsed.value)
        Jcoup_dict_path = my_path + "/Jcoup_dict.pickle"
        pickle.dump(Jcoup_dict, open(Jcoup_dict_path, "wb"))

        if Jcoup_dict:
            jcoup_data = csx_calc.calcJCouplings(
                my_csv_buffer, db_entry.karplus, Jcoup_dict, my_pdb, my_path
            )
            data_found = True
        else:
            jcoup_data = None

        # -----------------------  Chemical shift calc  --------------------- #
        chem_shift_lists = csx_func.parseChemShift_STR(parsed.value)
        chem_shift_lists_path = my_path + "/ChemShift_lists.pickle"
        pickle.dump(chem_shift_lists, open(chem_shift_lists_path, "wb"))

        if chem_shift_lists:
            chemshift_data = csx_calc.calcChemShifts(
                my_csv_buffer, chem_shift_lists, pdb_models, my_path
            )
            data_found = True
        else:
            chemshift_data = None
    else:
        str_name = "[NOT PRESENT]"
        rdc_data = None
        s2_data = None
        s2_sc_data = None
        jcoup_data = None
        chemshift_data = None

    my_csv_buffer.writeCSV()

    csx_func.calcPeptideBonds()
    csx_func.calcNH_Angles()

    if data_found:
        print(csx_obj.CalcPickle.data)
        calced_values = my_path + "/calced_values.p"
        pickle.dump(csx_obj.CalcPickle.data, open(calced_values, "wb"))
        rendered_page = render(request, "consensx/calculation.html", {
            "my_id": my_id,
            "my_pdb": db_entry.PDB_file,
            "n_model": model_count,
            "my_NOE": noe_name,
            "n_NOE": noe_n,
            "my_STR": str_name,
            "NOE_PRIDE_data": noe_pride_data,
            "RDC_data": rdc_data,
            "S2_data": s2_data,
            "S2_sc_data": s2_sc_data,
            "Jcoup_data": jcoup_data,
            "chemshift_data": chemshift_data,
            "SVD_calc": db_entry.svd_enable
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
