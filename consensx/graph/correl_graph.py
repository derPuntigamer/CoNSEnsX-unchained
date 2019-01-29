import numpy as np
import matplotlib.pyplot as plt

plt.switch_backend('Agg')


def correl_graph(my_path, calced, experimental, graph_name):
    """X axis -> experimental values, Y axis -> calculated values
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    min_calc = min(calced.values())
    max_calc = max(calced.values())

    exp_values = []
    for record in experimental:
        exp_values.append(record.value)

    min_exp = min(exp_values)
    max_exp = max(exp_values)
    miny = min(min_calc, min_exp)  # get minimum value
    maxy = max(max_calc, max_exp)  # get maximum value

    exp_line, calc_line = [], []

    for i, j in enumerate(calced.keys()):  # fetch data from arguments
        calc = calced[j]
        exp = experimental[i].value

        exp_line.append(exp)
        calc_line.append(calc)

    diag = []

    margin = int(abs(miny - maxy) * 0.05)

    if abs(miny - maxy) < 10:
        margin = 0.3
    elif abs(miny - maxy) < 2:
        margin = 0.01
    elif abs(miny - maxy) < 1:
        margin = 0

    maxy += margin
    miny -= margin

    for i in np.arange(miny, maxy * 1.42, 0.1):  # draw graph diagonal
        diag.append(i)

    plt.figure(figsize=(6, 5), dpi=80)
    plt.plot(diag, diag, linewidth=2.0, color='red', alpha=.7)
    plt.plot(exp_line, calc_line, 'bo')
    plt.axis([miny, maxy, miny, maxy])
    plt.xlabel('experimental')
    plt.ylabel('calculated')
    plt.tight_layout(pad=1.08)
    plt.savefig(my_path + "/" + graph_name, format="svg")
    plt.close()
