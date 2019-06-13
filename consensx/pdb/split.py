def split(my_path, pdb_file):
    """Split the given PDB file into models, each model becomes a separate
       PDB file placed in the "temp" folder"""
    try:
        my_pdb = open(pdb_file)
    except FileNotFoundError:
        print("ERROR -> the uploaded PDB file was not found")
        raise SystemExit
    except TypeError:
        print("ERROR -> PDB identifier is invalid")
        raise SystemExit

    model_names = []
    model_data = []

    my_name = ""
    my_data = []

    for line in my_pdb:
        if line.startswith("MODEL"):
            my_name = line.strip().split()[1]
        elif line.startswith("ATOM") or line.startswith("TER"):
            # replace oxygen names in line
            try:
                if line.split()[2] in ["OC1", "OT1"]:
                    line = line.replace(line.split()[2], "O ")
                elif line.split()[2] == "O1":
                    line = line.replace("O1 ", "O  ")
                elif line.split()[2] in ["OC2", "OT2"]:
                    line = line.replace(line.split()[2], "OXT")
                elif line.split()[2] == "O2":
                    line = line.replace("O2 ", "OXT")
            except IndexError:
                pass
            my_data.append(line.strip())
        elif line.startswith("ENDMDL"):
            if my_name in model_names:
                raise Exception("Looks like you have redundant model numbers in your PDB. It is easy to fix:")

            model_names.append(my_name)
            model_data.append(my_data)
            my_name = ""
            my_data = []
        else:
            continue

    my_pdb.close()

    for i in range(len(model_names)):
        file_name = my_path + "model_" + model_names[i] + ".pdb"
        temp_pdb = open(file_name, 'w')
        temp_pdb.write("HEADER    MODEL " + model_names[i] + "\n")

        for _ in model_data[i]:
            temp_pdb.write(_ + "\n")

        temp_pdb.write("END")
        temp_pdb.close()

    return len(model_names)
