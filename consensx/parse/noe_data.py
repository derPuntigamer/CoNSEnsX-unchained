def get_noe_from_file(noe_file):
    noe_key = "_Gen_dist_constraint_list.Constraint_type"
    noe_value = "NOE"
    in_frame = False
    in_loop = False
    loops = []
    noe_data = []
    loop_keys = []
    loop_data = []

    for line in open(noe_file, encoding="utf-8"):
        tok = line.strip().split()

        if not tok:
            continue

        if tok[0] == noe_key and tok[1] == noe_value:
            in_frame = True
            continue

        if in_frame and tok[0] == "save_":
            in_frame = False

        if tok[0] == "loop_":
            loop_keys = []
            loop_data = []
            in_loop = True
            continue

        if in_loop is True and in_frame is True and tok[0] == "stop_":
            in_loop = False
            loops.append([loop_keys, loop_data])

        if in_loop:
            if tok[0].startswith("_"):
                loop_keys.append(tok[0])
            else:
                loop_data.append(tok)

    for loop in loops:
        try:
            ind_id = loop[0].index("_Gen_dist_constraint.ID")
            ind_seg1 = loop[0].index("_Gen_dist_constraint.PDB_residue_no_1")
            ind_seg2 = loop[0].index("_Gen_dist_constraint.PDB_residue_no_2")
            ind_comp1 = loop[0].index("_Gen_dist_constraint.Comp_ID_1")
            ind_comp2 = loop[0].index("_Gen_dist_constraint.Comp_ID_2")
            ind_atom1 = loop[0].index("_Gen_dist_constraint.Atom_ID_1")
            ind_atom2 = loop[0].index("_Gen_dist_constraint.Atom_ID_2")
            ind_bnd = loop[0].index(
                "_Gen_dist_constraint.Distance_upper_bound_val"
            )
        except ValueError:
            continue

        for data in loop[1]:
            noe_data.append(
                [
                    data[ind_id],
                    data[ind_seg1],
                    data[ind_seg2],
                    data[ind_comp1],
                    data[ind_comp2],
                    data[ind_atom1],
                    data[ind_atom2],
                    data[ind_bnd],
                ]
            )

    return noe_data


class RestraintRecord:
    def __init__(self, my_path, noe_file, model_data):
        self.restraints = []
        self.resolved_restraints = []

        save_shifts = get_noe_from_file(my_path + noe_file)
        self.noe_n = save_shifts[-1][0] + " distance restraints found"

        csx_id = 1
        prev_id = 1

        for data in save_shifts:
            self.restraints.append(
                {
                    "csx_id": csx_id,
                    "curr_distID": int(data[0]),
                    "seq_ID1": int(data[1]),
                    "seq_ID2": int(data[2]),
                    "seq_name1": data[3],
                    "seq_name2": data[4],
                    "atom_ID1": data[5],
                    "atom_ID2": data[6],
                    "dist_max": float(data[7]),
                }
            )

            if prev_id != data[0]:
                prev_id = data[0]
                csx_id += 1

        self._resolve_names(model_data)

    def _resolve_names(self, model_data):
        pse_c = ["#", "*", "%", "+"]

        # load atom names from PDB data into a dict, residue IDs as keys
        noe_dict = {}
        pdb_atom_names = {}

        for atom in model_data.atomgroup:
            atom_res = atom.getResnum()

            if atom_res in pdb_atom_names.keys():
                pdb_atom_names[atom_res].append(atom.getName())
            else:
                pdb_atom_names[atom_res] = [atom.getName()]

        # organize restraint into a dict, with distIDs as keys
        for res in self.restraints:
            if res["csx_id"] in noe_dict.keys():
                noe_dict[res["csx_id"]].append(res)
            else:
                noe_dict[res["csx_id"]] = [res]

        # resolution of pseudo-atoms
        for rest_id in noe_dict.keys():
            for res in noe_dict[rest_id]:
                # skip O atoms
                if res["atom_ID1"] == "O" or res["atom_ID2"] == "O":
                    continue

                resol1 = [res["atom_ID1"]]
                atom_names1 = []

                # check if first atom is a pseudo atom
                pseud_h1 = (
                    res["atom_ID1"][0] == "H" and res["atom_ID1"][-1] in pse_c
                )

                if res["atom_ID1"][0] in "QM" or pseud_h1:
                    # change first letter of the name to hydrogen
                    # TODO check if if is removable
                    if res["atom_ID1"][0] in "QM":
                        base = "H" + res["atom_ID1"][1:]

                        # separate trailing number, if any
                        try:
                            num = str(int(res["atom_ID1"][-1]))
                            base = "H" + res["atom_ID1"][1:-1]
                        except ValueError:
                            num = 0

                        # permutations, if no trailing number found
                        if num != 0:
                            for i in ["1", "2", "3"]:
                                atom_names1.append(base + num + i)
                                atom_names1.append(num + base + i)
                        # permutations, if trailing number found
                        else:
                            if (
                                res["seq_name1"] == "ILE"
                                and res["atom_ID1"] == "MG"
                            ):
                                for i in ["1", "2", "3"]:
                                    atom_names1.append(base + str(2) + i)
                                    atom_names1.append(i + base + str(2))

                            elif (
                                res["seq_name1"] == "THR"
                                and res["atom_ID1"] == "MG"
                            ):
                                atom_names1.append("HG1")
                                atom_names1.append("1HG")

                            else:
                                for i in ["1", "2", "3"]:
                                    atom_names1.append(base + i)

                                    for j in ["1", "2", "3"]:
                                        atom_names1.append(base + j + i)
                                        atom_names1.append(i + base + j)

                    if pseud_h1:
                        atom_names1 = []
                        base = res["atom_ID1"][:-1]

                        for i in ["1", "2", "3"]:
                            atom_names1.append(base + i)
                            atom_names1.append(i + base)

                    # TODO KeyError
                    # get corresponding atom names present in PDB file
                    pdb_names1 = pdb_atom_names[res["seq_ID1"]]
                    resol1 = list(set(atom_names1) & set(pdb_names1))

                resol2 = [res["atom_ID2"]]

                # check if second atom is a pseudo atom
                pseud_h2 = (
                    res["atom_ID2"][0] == "H" and res["atom_ID2"][-1] in pse_c
                )

                if res["atom_ID2"][0] in "QM" or pseud_h2:
                    atom_names2 = []

                    if res["atom_ID2"][0] in "QM":
                        base = "H" + res["atom_ID2"][1:]

                        try:
                            num = str(int(res["atom_ID2"][-1]))
                            base = "H" + res["atom_ID2"][1:-1]
                        except ValueError:
                            num = 0

                        if num != 0:
                            for i in ["1", "2", "3"]:
                                atom_names2.append(base + num + i)
                                atom_names2.append(num + base + i)
                        else:
                            if (
                                res["seq_name2"] == "ILE"
                                and res["atom_ID2"] == "MG"
                            ):
                                for i in ["1", "2", "3"]:
                                    atom_names2.append(base + str(2) + i)
                                    atom_names2.append(i + base + str(2))

                            elif (
                                res["seq_name2"] == "THR"
                                and res["atom_ID2"] == "MG"
                            ):
                                atom_names2.append("HG1")
                                atom_names2.append("1HG")

                            else:
                                for i in ["1", "2", "3"]:
                                    atom_names2.append(base + i)

                                    for j in ["1", "2", "3"]:
                                        atom_names2.append(base + j + i)
                                        atom_names2.append(i + base + j)

                    if pseud_h2:
                        atom_names2 = []
                        base = res["atom_ID2"][:-1]

                        for i in ["1", "2", "3"]:
                            atom_names2.append(base + i)
                            atom_names2.append(i + base)

                    pdb_names2 = pdb_atom_names[res["seq_ID2"]]
                    resol2 = list(set(atom_names2) & set(pdb_names2))

                for atom1 in resol1:
                    for atom2 in resol2:
                        self.resolved_restraints.append(
                            {
                                "csx_id": res["csx_id"],
                                "curr_distID": res["curr_distID"],
                                "seq_ID1": res["seq_ID1"],
                                "seq_ID2": res["seq_ID2"],
                                "seq_name1": res["seq_name1"],
                                "seq_name2": res["seq_name2"],
                                "dist_max": res["dist_max"],
                                "atom_ID1": atom1,
                                "atom_ID2": atom2,
                            }
                        )

    def get_restraint_count(self):
        return self.restraints[-1]["csx_id"]

    def get_pride_restraints(self):
        pride_restraints = {}
        prev_id = -1
        seq1_ok, seq2_ok, distance_ok = False, False, False
        seq_dist = -1
        id_distance = -1

        for restraint in self.resolved_restraints:
            curr_id = restraint["curr_distID"]
            atom_id1 = restraint["atom_ID1"]
            atom_id2 = restraint["atom_ID2"]
            seq_id1 = restraint["seq_ID1"]
            seq_id2 = restraint["seq_ID2"]

            if prev_id != curr_id:
                prev_id = curr_id

                if seq1_ok and seq2_ok and distance_ok:
                    if seq_dist in pride_restraints:
                        pride_restraints[seq_dist] += 1
                    else:
                        pride_restraints[seq_dist] = 1

                atom1_ha = atom_id1 in ["H", "HA", "HA1", "HA2", "HA3"]
                seq1_ok = atom1_ha or atom_id1.startswith("HB")
                atom2_ha = atom_id2 in ["H", "HA", "HA1", "HA2", "HA3"]
                seq2_ok = atom2_ha or atom_id2.startswith("HB")
                seq_dist = abs(seq_id1 - seq_id2)
                distance_ok = seq_dist > 2

                id_distance = seq_dist

            else:
                atom1_ha = atom_id1 in ["H", "HA", "HA1", "HA2", "HA3"]
                seq1_ok &= atom1_ha or atom_id1.startswith("HB")
                atom2_ha = atom_id2 in ["H", "HA", "HA1", "HA2", "HA3"]
                seq2_ok &= atom2_ha or atom_id2.startswith("HB")
                seq_dist = abs(seq_id1 - seq_id2)
                distance_ok &= 2 < seq_dist == id_distance

        return pride_restraints
