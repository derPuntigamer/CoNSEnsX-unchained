import matplotlib
import random  # ID generation
import string  # ID generation
import prody

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

matplotlib.use('Agg')

chars = string.ascii_uppercase + string.digits


def create_pca_comparison(my_path, original, in_selection):
    sel_data1 = prody.parsePDB(my_path + '/' + original, subset='calpha')
    sel_ensemble1 = prody.Ensemble(sel_data1)

    sel_ensemble1.iterpose()

    color_list = ["blue" for _ in range(len(sel_ensemble1))]

    for i in [int(x)-1 for x in in_selection]:
        color_list[i] = "red"

    pca = prody.PCA("PCA1")
    pca.buildCovariance(sel_ensemble1)
    pca.calcModes()

    pca_image_names = []

    for i in range(0, 3):
        prody.showProjection(
            sel_ensemble1, pca[i:i+2], color=color_list
        )

        fig_hash = "".join(random.choice(chars) for _ in range(6))
        fig_name = "pca_mode_" + str(i+1) + str(i+2) + "_" + fig_hash + ".svg"
        pca_image_names.append(fig_name)
        red_patch = mpatches.Patch(color='red', label="Selection")
        blue_patch = mpatches.Patch(color='blue', label="Original")
        plt.legend(handles=[blue_patch, red_patch])
        plt.tight_layout(pad=1.08)
        plt.savefig(my_path + "/" + fig_name, format="svg")

        plt.close()

    print(pca_image_names)
    return pca_image_names
