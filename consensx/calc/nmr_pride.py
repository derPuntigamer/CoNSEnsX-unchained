import math
import subprocess
import os

from . import methods as csx_func
from . import objects as csx_obj


def nmr_pride(pdb_models, my_path):
    """Calculate NMR-PRIDE score on given PDB models"""

    pwd = os.getcwd()
    os.chdir(my_path)

    # write model list text file
    model_list = open("model_list.txt", 'w')

    for model in pdb_models:
        model_list.write(model + "\n")

    model_list.write("END\n")
    model_list.close()

    # write distance dict to text file
    restraints = csx_obj.Restraint_Record.getPRIDE_restraints()
    pride_input = open("pride_input.txt", 'w')

    pride_input.write("HEADER\n")

    prime_distances = list(restraints.keys())
    prime_distances.sort()

    for distance in prime_distances:
        pride_input.write(str(distance) + ' ' +
                          str(restraints[distance]) + '\n')

    pride_input.write("END\n")
    pride_input.close()

    # create binary database for PRIDE-NMR
    DEVNULL = open(os.devnull, 'w')
    hhdb_log = open("hhdb.log", 'w')
    model_list = open("model_list.txt", 'r')
    subprocess.call(
        [
            csx_obj.ThirdParty.prideDB,
            "-D", "HHDB",  # model list
        ],
        stdin=model_list,
        stdout=DEVNULL,
        stderr=hhdb_log
    )

    hhdb_log.close()
    DEVNULL.close()
    model_list.close()

    # run PRIDE-NMR
    DEVNULL = open(os.devnull, 'w')
    pride_input = open("pride_input.txt", 'r')
    pride_output = open("pride_output.txt", 'w')
    subprocess.call(
        [
            csx_obj.ThirdParty.prideNMR,
            "-D", "HHDB",
            "-d", str(56),
            "-b", str(len(pdb_models)),
            "-m", str(3)
        ],
        stdin=pride_input,
        stdout=pride_output,
        stderr=DEVNULL
    )

    pride_input.close()
    pride_output.close()
    DEVNULL.close()

    pride_scores = {}
    pride_output = open("pride_output.txt", 'r')
    for line in pride_output:
        if line.startswith("PRIDENMR:"):
            model_num = int(line.split()[-1])
            model_score = float(line.split()[1])
            pride_scores[model_num] = model_score

    scores = list(pride_scores.values())

    avg = sum(scores) * 1.0 / len(scores)
    variance = [(x - avg) ** 2 for x in scores]
    standard_deviation = math.sqrt(sum(variance) * 1.0 / len(variance))

    PRIDE_data = []

    print("PRIDE-NMR calculation")
    print("MAX: ",    max(pride_scores, key=pride_scores.get))
    PRIDE_data.append(max(pride_scores, key=pride_scores.get))
    print("MIN: ",    min(pride_scores, key=pride_scores.get))
    PRIDE_data.append(min(pride_scores, key=pride_scores.get))
    print("AVG: ", avg)
    PRIDE_data.append(avg)
    print("DEV: ", standard_deviation, "\n")
    PRIDE_data.append(standard_deviation)

    os.chdir(pwd)

    csx_func.makeNMRPrideGraph(my_path, scores, avg)

    return PRIDE_data
