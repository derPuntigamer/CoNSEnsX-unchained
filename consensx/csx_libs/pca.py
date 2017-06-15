import matplotlib
matplotlib.use('Agg')
matplotlib.verbose.set_level('silent')

import random  # ID generation
import string  # ID generation
import prody

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

chars = string.ascii_uppercase + string.digits


def create_PCA_comparison(my_path, original, in_selection):

    sel_data1 = prody.parsePDB(my_path + '/' + original, subset='calpha')
    sel_ensemble1 = prody.Ensemble(sel_data1)

    sel_ensemble1.iterpose()

    color_list = ["blue" for i in range(len(sel_ensemble1))]

    print("color_list", color_list)
    print("LEN color_list", len(color_list))
    print("LEN in_selection", len(in_selection))



    for i in [int(x)-1 for x in in_selection]:
        print("I is ", i)
        color_list[i] = "red"

    pca = prody.PCA("PCA1")
    pca.buildCovariance(sel_ensemble1)
    pca.calcModes()

    pca_image_names = []

    for i in range(0, 3):
        prody.showProjection(
            sel_ensemble1, pca[i:i+2], color=color_list
        )

        fig_hash = ''.join(random.choice(chars) for _ in range(6))
        fig_name = "pca_mode_" + str(i+1) + str(i+2) + "_" + fig_hash + ".svg"
        pca_image_names.append(fig_name)
        red_patch = mpatches.Patch(color='red', label='Selection')
        blue_patch = mpatches.Patch(color='blue', label='Original')
        plt.legend(handles=[blue_patch, red_patch])
        plt.tight_layout(pad=1.08)
        plt.savefig(my_path + "/" + fig_name, format="svg")

        plt.close()

    print(pca_image_names)
    return pca_image_names

# x1 = []
# x2 = []
# y1 = []
# y2 = []

# for p in range(len(projection)):
#     x1.append(projection[p][0])
#     y1.append(projection[p][1])

# for p in range(len(projection_sel)):
#     x2.append(projection_sel[p][0])
#     y2.append(projection_sel[p][1])

# aa_file = open("aa", "w")
# bb_file = open("bb", "w")

# for i in range(len(x1)):
#     aa_file.write("{0} {1}\n".format(x1[i], y1[i]))

# for i in range(len(x2)):
#     bb_file.write("{0} {1}\n".format(x2[i], y2[i]))

# aa_file.write("END\n")
# bb_file.write("END\n")

# aa_file.close()
# bb_file.close()

# print("LEN sel_ensemble1", len(sel_ensemble1))
#
# print("LEN color_list", len(color_list))
# print("color_list", color_list)
