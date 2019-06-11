import matplotlib.pyplot as plt

plt.switch_backend('Agg')


def mod_correl_graph(my_path, correl, avg_corr, model_corrs, corr_graph_name):
    """Y axis -> correlation values
       X axis -> ensemble correlation, model avg. correlation,
                 per model correlation
       parameter 'model_corrs' is a list containing per model
       correlation values
       """
    plt.figure(figsize=(6, 5), dpi=80)

    plt.plot(
        list(range(0, len(model_corrs))), [correl] * len(model_corrs),
        linewidth=2.0, color='green', label='Ensemble corr.', alpha=.7
    )

    plt.plot(
        list(range(0, len(model_corrs))), [avg_corr] * len(model_corrs),
        linewidth=2.0, color='red', label='Avg. corr. per model', alpha=.7
    )

    plt.plot(
        list(range(0, len(model_corrs))), sorted(model_corrs),
        linewidth=2.0, color='blue', label='Corr. per model', alpha=.7
    )

    plt.legend(loc='lower left')
    plt.axis([-1, len(model_corrs), 0, 1])
    plt.xlabel('models (worse to best)')
    plt.ylabel('correlation')
    plt.grid(axis="y")
    plt.tight_layout(pad=1.08)
    plt.savefig(my_path + "/" + corr_graph_name, format="svg")
    plt.close()
