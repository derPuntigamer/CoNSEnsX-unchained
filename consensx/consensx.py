r"""
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
import zipfile

# own modules
import consensx.calc as calc
import consensx.parse as parse
import consensx.storage as storage
import consensx.pdb
from consensx.misc.natural_sort import natural_sort
from consensx import thirdparty


# Django server
from django.shortcuts import render
from django.http import HttpResponse
from .models import CSX_upload, CSX_calculation

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run_calculation(request, calc_id):
    abspath = os.path.abspath(__file__)
    dirname = os.path.dirname(abspath)
    thirdparty.check_3rd_party(dirname)

    calced_data_storage = dict()
    my_id = calc_id
    print("ID IS:", my_id)

    from django.conf import settings

    my_path = os.path.join(settings.MEDIA_ROOT, my_id) + "/"
    db_entry = CSX_upload.objects.get(id_code=calc_id)

    # ----------------  Setting up working directory and files  ------------- #
    csv_buffer = storage.CSVBuffer(my_path)

    my_pdb = my_path + db_entry.PDB_file

    consensx.pdb.clean(my_path, my_pdb, csv_buffer)
    model_count = consensx.pdb.split(my_path, my_pdb)

    model_data = parse.Pdb(my_pdb, my_path, model_count)

    if not model_data:
        return render(
            request,
            "consensx/home.html",
            {
                "error": "DISCARDED MODELS FOUND, CHECK IF ALL MODELS HAVE THE SAME\
            NUMBER OF ATOMS"
            },
        )

    pdb_models = []

    for file in os.listdir(my_path):
        if file.startswith("model_") and file.endswith(".pdb"):
            pdb_models.append(file)

    pdb_models = natural_sort(pdb_models)
    data_found = False
    noe_name = "[NOT PRESENT]"
    noe_n = ""
    noe_pride_data = None

    str_name = "[NOT PRESENT]"
    rdc_calced_data = None
    s2_data = None
    s2_sc_data = None
    jcoup_data = None
    chemshift_data = None

    # --------------------  Read  and parse BME weights   ------------------- #
    bme_weights = None

    if db_entry.bme_weights_file:
        print("db_entry.bme_weights_file", db_entry.bme_weights_file)
        f = open(my_path + db_entry.bme_weights_file).read().split()
        f = [float(i) for i in f]
        bme_weights = f

    # ---------------------  Read  and parse NOE file   --------------------- #
    if db_entry.NOE_file:
        noe_restraints = parse.RestraintRecord(
            my_path, db_entry.NOE_file, model_data
        )

        noe_n = noe_restraints.get_restraint_count()

        noe_violations = calc.noe_violations(
            model_data, my_path, db_entry, noe_restraints, bme_weights
        )
        pride_data = calc.nmr_pride(pdb_models, my_path, noe_restraints)

        noe_pride_data = {
            "noe_violations": str(noe_violations),
            "NOE_hist": my_id + "/NOE_hist.svg",
            "best_score": pride_data[0],
            "worst_score": pride_data[1],
            "average_score": "{0:.3f}".format(pride_data[2]),
            "deviation": "{0:.3f}".format(pride_data[3]),
            "PRIDE_hist": my_id + "/PRIDE-NMR_score.svg",
        }

        data_found = True

    # ---------------------  Read  and parse STR file   --------------------- #
    if not db_entry.STR_file:
        bme_zf = zipfile.ZipFile(my_path + "bme_inputs.zip", mode="w")

        for file in os.listdir(my_path):
            if file.endswith("_exp.dat") or file.endswith("_calc.dat"):
                bme_zf.write(my_path + file, file)

        if bme_zf:
            bme_zf.close()

        rendered_page = render(
            request,
            "consensx/calculation.html",
            {
                "my_id": my_id,
                "my_pdb": db_entry.PDB_file,
                "n_model": model_count,
                "my_NOE": noe_name,
                "n_NOE": noe_n,
                "my_STR": str_name,
                "NOE_PRIDE_data": noe_pride_data,
                "RDC_data": rdc_calced_data,
                "S2_data": s2_data,
                "S2_sc_data": s2_sc_data,
                "Jcoup_data": jcoup_data,
                "chemshift_data": chemshift_data,
                "SVD_calc": db_entry.svd_enable,
                "bme_zf": True,
            },
        )

        post_data = CSX_calculation(
            html_content=rendered_page.content, id_code=my_id
        )
        post_data.save()

        return rendered_page

    str_name = db_entry.STR_file
    nmr_file_path = my_path + db_entry.STR_file

    try:
        star_nmr_data = parse.StarNMR(nmr_file_path)
    except Exception as e:
        print("EXCEPTION", e)
        return HttpResponse(e)

    # ------------------------  RDC calculation  ------------------------ #
    rdc_lists = star_nmr_data.parse_rdc()
    rdc_lists_path = my_path + "/RDC_lists.pickle"
    pickle.dump(rdc_lists, open(rdc_lists_path, "wb"))
    rdc_calced_data = None

    if rdc_lists:
        svd_enabled = db_entry.svd_enable
        lc_model = db_entry.rdc_lc
        rdc_calced_data = calc.rdc(
            csv_buffer,
            calced_data_storage,
            rdc_lists,
            pdb_models,
            my_path,
            svd_enabled,
            lc_model,
        )
        data_found = True

    # ----------------------------  S2 calc  ---------------------------- #
    s2_dict = star_nmr_data.parse_s2()
    s2_dump = [s2_dict, db_entry.superimpose, db_entry.fit_range]
    s2_dict_path = my_path + "/S2_dict.pickle"
    pickle.dump(s2_dump, open(s2_dict_path, "wb"))
    s2_data = None

    if s2_dict:
        s2_data = calc.s2(
            csv_buffer,
            calced_data_storage,
            s2_dict,
            my_path,
            model_data,
            fit=db_entry.superimpose,
            fit_range=db_entry.fit_range,
        )

        data_found = True

    s2_sidechain = star_nmr_data.parse_s2_sidechain()
    s2_sc_data = None

    if s2_sidechain:
        s2_sc_data = calc.s2_sidechain(
            csv_buffer,
            s2_sidechain,
            my_path,
            model_data,
            fit=db_entry.superimpose,
        )

        if "error" in s2_sc_data.keys():
            return render(
                request, "consensx/home.html", {"error": s2_sidechain["error"]}
            )

        data_found = True

    # ------------------------  J-coupling calc  ------------------------ #
    Jcoup_dict = star_nmr_data.parse_jcoup()
    Jcoup_dict_path = my_path + "/Jcoup_dict.pickle"
    pickle.dump(Jcoup_dict, open(Jcoup_dict_path, "wb"))
    jcoup_data = None

    if Jcoup_dict:
        jcoup_data = calc.jcoupling(
            csv_buffer,
            calced_data_storage,
            model_data,
            db_entry,
            Jcoup_dict,
            my_path,
            bme_weights,
        )
        data_found = True

    # -----------------------  Chemical shift calc  --------------------- #
    chem_shift_lists = star_nmr_data.parse_chemshift()
    chem_shift_lists_path = my_path + "/ChemShift_lists.pickle"
    pickle.dump(chem_shift_lists, open(chem_shift_lists_path, "wb"))
    chemshift_data = None

    if chem_shift_lists:
        chemshift_data = calc.chemshifts(
            csv_buffer,
            calced_data_storage,
            chem_shift_lists,
            pdb_models,
            my_path,
            bme_weights,
        )
        data_found = True

    csv_buffer.write_csv()

    calc.peptide_bonds(model_data)
    calc.nh_angles(model_data)

    # pack BME input files for zipping
    bme_zf = None

    for file in os.listdir(my_path):
        if file.endswith("_exp.dat") or file.endswith("_calc.dat"):
            if not bme_zf:
                bme_zf = zipfile.ZipFile(my_path + "bme_inputs.zip", mode="w")

            bme_zf.write(my_path + file, file)

    if bme_zf:
        bme_zf.close()

    if data_found:
        print(calced_data_storage)
        calced_values = my_path + "/calced_values.p"
        pickle.dump(calced_data_storage, open(calced_values, "wb"))
        rendered_page = render(
            request,
            "consensx/calculation.html",
            {
                "my_id": my_id,
                "my_PDB": db_entry.PDB_file,
                "n_model": model_count,
                "my_NOE": noe_name,
                "n_NOE": noe_n,
                "my_STR": str_name,
                "NOE_PRIDE_data": noe_pride_data,
                "RDC_data": rdc_calced_data,
                "S2_data": s2_data,
                "S2_sc_data": s2_sc_data,
                "Jcoup_data": jcoup_data,
                "chemshift_data": chemshift_data,
                "SVD_calc": db_entry.svd_enable,
                "bme_zf": bool(bme_zf),
            },
        )

        print("RENDERED PAGE ---------------- START")
        post_data = CSX_calculation(
            html_content=rendered_page.content, id_code=my_id
        )
        post_data.save()

        print("RENDERED PAGE ------------------ END")
        return rendered_page
    else:
        print("NO DATA FOUND IN STAR-NMR FILE")
        return render(
            request,
            "consensx/home.html",
            {"error": "NO DATA FOUND IN STAR-NMR FILE"},
        )
