import consensx.nmrpystar as nmrpystar
import pynmrstar


class RDC_Record(object):
    """Class for storing RDC data"""

    def __init__(self, resnum1, atom1, resnum2, atom2, RDC_value):
        self.RDC_type = (str(int(resnum1) - int(resnum2))
                         + '_' + atom1 + '_' + atom2)
        self.resnum = int(resnum1)
        self.atom = atom1
        self.resnum2 = int(resnum2)
        self.atom2 = atom2
        self.value = float(RDC_value)


class S2_Record(object):
    """Class for storing S2 data"""

    def __init__(self, resnum, S2_type, S2_value):
        self.resnum = int(resnum)
        self.type = S2_type
        self.value = float(S2_value)
        self.calced = None


class JCoup_Record(object):
    """Class for storing J-Coupling data"""

    def __init__(self, resnum, jcoup_type, JCoup_value):
        self.resnum = int(resnum)
        self.type = jcoup_type
        self.value = float(JCoup_value)


class ChemShift_Record(object):
    """Class for storing chemical shift data"""

    def __init__(self, resnum, res_name, atom_name, ChemShift_value):
        self.resnum = int(resnum)
        self.res_name = res_name
        self.atom_name = atom_name
        self.value = float(ChemShift_value)

    def __str__(self):
        return (
            "resnum: " + str(self.resnum) + "," +
            "res_name: " + self.res_name + "," +
            "atom_name: " + self.atom_name + "," +
            "value: " + str(self.value) + "\n"
        )

    def __repr__(self):
        return self.__str__()


class StarNMR():
    def __init__(self, STR_FILE):
        """Parse BMRB file into a python object"""

        parse_exception = (
            "ERROR during STR parsing, please check your STAR-NMR file!"
        )

        self.parsed = pynmrstar.Entry.from_file(STR_FILE)

        # try:
        #     self.parsed = pynmrstar.Entry.from_file(STR_FILE)
        # except Exception:
        #     raise Exception(parse_exception)

    def parse_rdc(self):
        """Returns RDC lists as dictonaries containing RDC_Record objects,
        grouped by RDCtype (keys())"""

        list_number = 1
        RDC_lists = []

        while True:
            saveShiftName = 'CNS/XPLOR_dipolar_coupling_' + str(list_number)
            try:
                saveShifts = self.parsed.value.saves[saveShiftName]
            except KeyError:
                break
            loopShifts = saveShifts.loops[-1]
            RDC_records = []

            # STR key values recognised by this program
            rdc_res1_keys = ["RDC.Seq_ID_1", "Atom_one_residue_seq_code", "RDC_constraint.Seq_ID_1"]
            rdc_atom1_keys = ["RDC.Atom_type_1", "Atom_one_atom_name", "RDC_constraint.Atom_ID_1"]
            rdc_res2_keys = ["RDC.Seq_ID_2", "Atom_two_residue_seq_code", "RDC_constraint.Seq_ID_2"]
            rdc_atom2_keys = ["RDC.Atom_type_2", "Atom_two_atom_name", "RDC_constraint.Atom_ID_2"]
            rdc_value_keys = ["RDC.Val", "Residual_dipolar_coupling_value", "RDC_constraint.RDC_val"]

            for ix in range(len(loopShifts.rows)):  # fetch values from file
                row = loopShifts.getRowAsDict(ix)

                for my_resnum1 in rdc_res1_keys:  # fetch 1. residue number
                    if my_resnum1 in list(row.keys()):
                        resnum1 = row[my_resnum1]

                for my_atom1 in rdc_atom1_keys:  # fetch 1. atom name
                    if my_atom1 in list(row.keys()):
                        atom1 = row[my_atom1]

                for my_resnum2 in rdc_res2_keys:  # fetch 2. residue number
                    if my_resnum2 in list(row.keys()):
                        resnum2 = row[my_resnum2]

                for my_atom2 in rdc_atom2_keys:  # fetch 2. atom name
                    if my_atom2 in list(row.keys()):
                        atom2 = row[my_atom2]

                for my_RDC_value in rdc_value_keys:  # fetch RDC value
                    if my_RDC_value in list(row.keys()):
                        RDC_value = float(row[my_RDC_value])
                        RDC_value = float("{0:.2f}".format(RDC_value))

                # check if all parameters are fetched
                if (resnum1 and atom1 and resnum2 and atom2 and RDC_value):
                    # append RDC_Record object to list
                    RDC_records.append(
                        RDC_Record(
                            resnum1, atom1, resnum2, atom2, RDC_value
                        )
                    )
                else:
                    print(row)

            RDC_lists.append(RDC_records)
            list_number += 1

        # split list into dict according to RDC types
        new_RDC_list = []
        for RDC_list in RDC_lists:
            prev_type = ""
            RDC_dict = {}

            for record in RDC_list:
                if prev_type != record.RDC_type:
                    RDC_dict[record.RDC_type] = []
                    RDC_dict[record.RDC_type].append(record)
                else:
                    RDC_dict[record.RDC_type].append(record)

                prev_type = record.RDC_type

            new_RDC_list.append(RDC_dict)

        return new_RDC_list

    def parse_s2(self):
        """Returns a dictonary with the parsed S2 data"""

        try:
            saveShifts = self.parsed.value.saves["order_param"]

            loopShifts = saveShifts.loops[-1]
            s2_records = []

            for ix in range(len(loopShifts.rows)):  # fetch values from file
                row = loopShifts.getRowAsDict(ix)

                S2_value = float("{0:.2f}".format(float(row["S2_value"])))

                s2_records.append(
                    S2_Record(
                        row["Residue_seq_code"], row["Atom_name"], S2_value)
                )

            # split list into dict according to S2 types
            S2_dict = {}
            prev_type = ""

            for record in s2_records:
                if prev_type != record.type:
                    S2_dict[record.type] = []
                    S2_dict[record.type].append(record)
                else:
                    S2_dict[record.type].append(record)

                prev_type = record.type

            return S2_dict

        except KeyError:
            print("No S2 parameter list found")
            return None

    def parse_s2_sidechain(self):
        """Returns a dictonary with the parsed S2 data"""

        sidechain_name = "side-chain_methyl_order_parameters"
        try:
            saveShifts = self.parsed.value.saves[sidechain_name]

            loopShifts = saveShifts.loops[-1]
            s2_records = []

            for ix in range(len(loopShifts.rows)):   # fetch values from file
                row = loopShifts.getRowAsDict(ix)

                S2_value = float("{0:.2f}".format(float(row["S2_value"])))

                s2_records.append(
                    S2_Record(
                        row["Residue_seq_code"], row["Atom_name"], S2_value
                    )
                )

            return s2_records

        except KeyError:
            return None

    def parse_jcoup(self):
        """Returns a dictonary with the parsed J-coupling data"""
        try:
            saveShifts = self.parsed.value.saves["coupling_constants"]

            loopShifts = saveShifts.loops[-1]
            jcoup_records = []

            for ix in range(len(loopShifts.rows)):   # fetch values from file
                row = loopShifts.getRowAsDict(ix)

                JC_value = row["Coupling_constant_value"]
                JC_value = float("{0:.2f}".format(float(JC_value)))

                jcoup_records.append(
                    JCoup_Record(
                        row["Atom_one_residue_seq_code"],
                        row["Coupling_constant_code"],
                        JC_value
                    )
                )

            # split list into dict according to J-cuopling types
            jcoup_dict = {}
            prev_type = ""

            for record in jcoup_records:
                if prev_type != record.type:
                    jcoup_dict[record.type] = []
                    jcoup_dict[record.type].append(record)
                else:
                    jcoup_dict[record.type].append(record)

                prev_type = record.type

            return jcoup_dict

        except KeyError:
            print("No J-coupling parameter list found")
            return None

    def parse_chemshift(self):
        """Returns ChemShift lists as dictionaries containing ChemShift_Record
        objects, grouped by Atom_name (keys())"""

        chemshift_types = ["HA", "CA", "CB", "N", "H", "C"]
        tag_list = ['Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Val']
        chemshift_records = []
        chem_shift_lists = []
        ha_sum = 0.0

        for cs_loop in self.parsed.get_loops_by_category("Atom_chem_shift"):
            for row_data in cs_loop.get_tag(tag_list):
                row = dict(zip(tag_list, row_data))
                # ['21', 'PRO', 'HD3', 'H', '3.448']

                cs_value = float("{0:.2f}".format(float(row["Val"])))

                if row["Atom_ID"] == "HA2":
                    ha_sum += cs_value
                    continue

                if row["Atom_ID"] == "HA3":
                    ha_sum += cs_value

                    chemshift_records.append(
                        ChemShift_Record(
                            row["Seq_ID"],
                            row["Comp_ID"],
                            "HA", ha_sum / 2
                        )
                    )
                    ha_sum = 0.0
                    continue

                if row["Atom_ID"] in chemshift_types:
                    chemshift_records.append(
                        ChemShift_Record(
                            row["Seq_ID"],
                            row["Comp_ID"],
                            row["Atom_ID"],
                            cs_value
                        )
                    )

                elif row["Atom_ID"] == "HN":
                    chemshift_records.append(
                        ChemShift_Record(
                            row["Seq_ID"],
                            row["Comp_ID"],
                            "H",
                            cs_value
                        )
                    )

            chem_shift_lists.append(chemshift_records)

        new_cs_list = []

        for chem_shift_list in chem_shift_lists:
            chem_shift_dict = {}

            for record in chem_shift_list:

                if record.atom_name in list(chem_shift_dict.keys()):
                    chem_shift_dict[record.atom_name].append(record)
                else:
                    chem_shift_dict[record.atom_name] = []
                    chem_shift_dict[record.atom_name].append(record)

            new_cs_list.append(chem_shift_dict)

        return new_cs_list
