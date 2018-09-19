from django.shortcuts import render
from django.http import HttpResponse
from django.views.decorators.csrf import csrf_exempt

import random  # ID generation
import string  # ID generation
import os      # mkdir
import pickle
import json

from .models import CSX_upload, CSX_calculation
from .consensx import run_calculation
from .selection import run_selection


chars = string.ascii_uppercase + string.digits

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def handle_uploaded_file(f, path, filename):
    file_full_path = path + '/' + filename
    with open(file_full_path, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def db(request, my_id):
    db_entry = CSX_calculation.objects.get(id_code=my_id)
    print("PATH", request.path)
    res = db_entry.returnHTML()
    res = res[8:].replace("\\n", "\n")
    return HttpResponse(res)


def home(request):
    if request.method == 'POST':  # if the form has been submitted...
        # generate ID for calculation
        my_id = ''.join(random.choice(chars) for _ in range(6))
        my_path = os.path.join(BASE_DIR, 'media', my_id)
        os.mkdir(my_path)

        # check if POST is a test submit
        if request.POST.get('submit_test', False):
            # IMPLEMENT TEST CALC HERE!
            return run_calculation(request, my_id)

        pdb_file = request.FILES['pdb_upload']  # get PDB file
        handle_uploaded_file(pdb_file, my_path, pdb_file.name)

        try:                                    # get restraint file if any
            restraint_file = request.FILES['bmrb_upload']
            restraint_file_name = restraint_file.name
            handle_uploaded_file(restraint_file, my_path, restraint_file_name)
        except KeyError:
            restraint_file_name = None

        try:                                    # get NOE file if any
            noe_file = request.FILES['xplor_upload']
            noe_file_name = noe_file.name
            handle_uploaded_file(noe_file, my_path, noe_file_name)
        except KeyError:
            noe_file_name = None

        try:                                    # check if fitting is enabled
            fit_enable = bool(request.POST['superimpose'])
        except KeyError:
            fit_enable = False

        try:                                    # get fit range if any
            fit_range = request.POST['fit_range']
        except KeyError:
            fit_range = None

        try:                                    # check if fitting is enabled
            r3average = bool(request.POST['r3average'])
        except KeyError:
            r3average = False

        try:                                    # check if fitting is enabled
            svd_enable = bool(request.POST['RDCSVD'])
        except KeyError:
            svd_enable = False

        post_data = CSX_upload(
            id_code=my_id,
            PDB_file=pdb_file.name,
            NOE_file=noe_file_name,
            STR_file=restraint_file_name,
            karplus=request.POST['KARPLUS'],
            superimpose=fit_enable,
            fit_range=fit_range,
            r3average=r3average,
            svd_enable=svd_enable,
            rdc_lc=request.POST['RDCLC']
        )
        post_data.save()

        return run_calculation(request, my_id)

    else:
        return render(request, "consensx/home.html")


@csrf_exempt
def selection(request, my_id):
    print("SELECTION ID IS: " + my_id)
    my_path = os.path.join(BASE_DIR, 'media', my_id)
    user_selection = json.loads(request.body.decode("utf-8"))

    original_values = pickle.load(open(my_path + "/calced_values.p", "rb"))

    num_coordsets, sel_values, pca_image_names = run_selection(
        my_path, original_values, user_selection
    )

    measure = user_selection["MEASURE"]

    return_dict = {
        "measure": measure,
        "num_coordsets": num_coordsets,
        "pca_image_names": pca_image_names
}

    if measure == "correlation":
        measure = "corr"

    if measure == "q-value":
        measure = "qval"

    values_dict = {}

    for key, value in sel_values.items():
        values_dict[key] = {
            "original": original_values[key + "_" + measure],
            "selection": "{0:.3g}".format(value)
        }

    print("values_dict")
    print(values_dict)

    return_dict["values"] = values_dict

    return HttpResponse(
        json.dumps(return_dict),
        content_type='application/json'
    )
