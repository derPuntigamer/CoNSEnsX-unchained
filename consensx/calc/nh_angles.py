from .vec_3d import Vec3D


def nh_angles(model_data):
    """Calculates backbone dihedral angles (OMEGA) H-N-C=O"""
    dihedral_angles = {"<2":    0, "2-5": 0, "5-10": 0,
                       "10-20": 0, ">20": 0}

    for i in range(model_data.coordsets):
        model_data.atomgroup.setACSIndex(i)
        current_resindex = 1
        prev_o, prev_c, = None, None
        my_n, my_h, my_o, my_c = None, None, None, None

        for atom in model_data.atomgroup:
            atom_res = atom.getResindex() + 1

            if atom_res != current_resindex:

                if (prev_o is not None and prev_c is not None and
                        my_n is not None and my_h is not None):

                    nh_vec = my_h - my_n
                    cn_vec = my_n - prev_c
                    oc_vec = prev_c - prev_o

                    first_cross = Vec3D.cross(nh_vec, cn_vec)
                    second_cross = Vec3D.cross(oc_vec, cn_vec)

                    angle = Vec3D.dihed_angle(
                        first_cross, second_cross
                    )

                    # reference for setting sign of angle
                    reference = Vec3D.cross(first_cross, second_cross)

                    r1 = reference.normalize()
                    r2 = nh_vec.normalize()

                    if (r1 - r2).magnitude() < r2.magnitude():
                        angle *= -1

                    if abs(angle) < 2:
                        dihedral_angles["<2"] += 1
                    elif abs(angle) < 5:
                        dihedral_angles["2-5"] += 1
                    elif abs(angle) < 10:
                        dihedral_angles["5-10"] += 1
                    elif abs(angle) < 20:
                        dihedral_angles["10-20"] += 1
                    else:
                        dihedral_angles[">20"] += 1

                current_resindex = atom_res
                prev_o = my_o
                prev_c = my_c
                my_n, my_h = None, None

            if atom_res == current_resindex:
                if atom.getName() == 'N':
                    my_n = Vec3D(atom.getCoords())
                elif atom.getName() == 'H':
                    my_h = Vec3D(atom.getCoords())
                elif atom.getName() == 'C':
                    my_c = Vec3D(atom.getCoords())
                elif atom.getName() == 'O':
                    my_o = Vec3D(atom.getCoords())

    print("Peptide (H-N-C=O) bond angle distribution:")
    print("   <2 -> " + str(dihedral_angles["<2"]))
    print("  2-5 -> " + str(dihedral_angles["2-5"]))
    print(" 5-10 -> " + str(dihedral_angles["5-10"]))
    print("10-20 -> " + str(dihedral_angles["10-20"]))
    print("  >20 -> " + str(dihedral_angles[">20"]))
