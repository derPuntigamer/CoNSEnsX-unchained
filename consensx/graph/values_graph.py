import copy
import matplotlib.pyplot as plt
import numpy as np

plt.switch_backend('Agg')


def values_graph(my_path, calced, my_experimental, graph_name):
    """X axis -> residue numbers, Y axis -> values
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    experimental = copy.deepcopy(my_experimental)

    exp_line, calc_line = [], []

    for k in range(min(calced.keys()) - 1, max(calced.keys()) + 1):
        if k in list(calced.keys()):
            calc = calced[k]
            exp = experimental.pop(0).value

            exp_line.append(exp)
            calc_line.append(calc)

        else:
            exp_line.append(None)  # append 'None' where data is missing
            calc_line.append(None)

    # connect line over missing (None) values, more info at ->
    # http://stackoverflow.com/questions/14399689/
    # matplotlib-drawing-lines-between-points-ignoring-missing-data
    exp_line = np.array(exp_line).astype(np.double)
    exp_mask = np.where(np.isfinite(exp_line))
    calc_line = np.array(calc_line).astype(np.double)
    calc_mask = np.where(np.isfinite(calc_line))

    # x axis values as numpy array
    xs = np.arange(min(calced.keys())-1, max(calced.keys())+2)

    plt.figure(figsize=(10, 5), dpi=80)

    # experimental values with 'None' values masked
    plt.plot(
        xs[exp_mask], exp_line[exp_mask],
        linewidth=2.0, color='red', marker='o', label='exp', alpha=.7
    )

    # calculated values with 'None' values masked
    plt.plot(
        xs[calc_mask], calc_line[calc_mask],
        linewidth=2.0, color='blue', marker='o', label='calc', alpha=.7
    )

    plt.legend(loc='lower left')
    plt.xlabel('residue number')
    plt.ylabel('value')
    plt.grid(axis="y")
    plt.tight_layout(pad=1.08)
    plt.savefig(my_path + "/" + graph_name, format="svg")
    plt.close()
