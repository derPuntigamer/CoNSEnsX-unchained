import os


def clean(my_path, pdb_file, csv_buffer):
    """Performs some basic formatting on the given PDB file to make it suitable
       for further calculations (from RCSB to BMRB PDF format)"""
    try:
        input_pdb = open(pdb_file)
    except FileNotFoundError:
        print(pdb_file + " was not found")
        raise SystemExit

    work_pdb = my_path + "my_pdb.pdb"

    my_pdb = open(work_pdb, 'w')
    max_resnum = 0
    min_resnum = 100000

    for line in input_pdb:
        line = line.strip()
        # line = re.sub('[+-] ', '  ', line)

        if line.startswith("ATOM"):

            name = line[12:16].strip()
            resnum = int(line[22:26])

            if resnum >= max_resnum:
                max_resnum = resnum
            if resnum < min_resnum:
                min_resnum = resnum

            if name == "Q":
                continue
            if name == "NH":
                name = "H"
            if name == "HN":
                name = "H"

            chars, numbers = [], []
            for i in name:
                try:
                    numbers.append(int(i))
                except ValueError:
                    chars.append(i)

            # try:
            #     # _ = int(name[0])
            #     name = (''.join(str(i) for i in chars) +
            #             ''.join(str(i) for i in reversed(numbers)))
            # except ValueError:
            #     pass

            if len(name) == 1:
                name = " " + name + "  "
            elif len(name) == 2:
                name = " " + name + " "
            elif len(name) == 3:
                name = " " + name

            my_pdb.write(line[:11] + " %4s" % name + line[16:21] +
                         'A' + line[22:] + "\n")
            continue

        my_pdb.write(line + "\n")

    input_pdb.close()
    my_pdb.close()

    os.remove(pdb_file)
    os.rename(work_pdb, pdb_file)

    print("MAX RESNUM", max_resnum)
    csv_buffer.max_resnum = max_resnum
    csv_buffer.min_resnum = min_resnum
