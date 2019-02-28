import math


class CalcPickle(dict):
    pass
    """Class for storing values for pickle generation"""
    data = dict()


class ThirdParty(object):

    """Class to store 3rd party software information"""
    pales = ""
    shiftx = ""
    prideDB = ""
    prideNMR = ""

    @staticmethod
    def get_thirdparty(config_file):
        cfg = open(config_file)
        for line in cfg:
            if line.startswith("#"):
                continue
            if line.startswith("pales"):
                ThirdParty.pales = line.split("'")[1]
            elif line.startswith("shiftx"):
                ThirdParty.shiftx = line.split("'")[1]
            elif line.startswith("prideDB"):
                ThirdParty.prideDB = line.split("'")[1]
            elif line.startswith("prideNMR"):
                ThirdParty.prideNMR = line.split("'")[1]


class CSV_buffer(object):

    """Class which stores data for values.CSV"""

    def __init__(self, my_path):
        self.working_dir = my_path
        self.max_resnum = -1
        self.min_resnum = 100000
        self.csv_data = []
        # self.name = name
        # self.calced = calced
        # self.exp = {}

        # for i in experimental:
        #     self.exp[i.resnum] = i.value

        # csv_data.append(self)

    def writeCSV(self):
        filename = self.working_dir + "values.csv"
        output_csv = open(filename, 'w')
        output_csv.write(',')
        for data in self.csv_data:
            output_csv.write(data["name"] + " EXP, " + data["name"] + " CALC,")
        output_csv.write("\n")

        # for i in experimental:
        #     self.exp[i.resnum] = i.value

        for resnum in range(self.min_resnum, self.max_resnum + 1):
            output_csv.write(str(resnum) + ',')
            for data in self.csv_data:
                exp = {}

                for i in data["experimental"]:
                    exp[i.resnum] = i.value

                try:
                    # pdb.set_trace()
                    output_csv.write(
                        "{0:.2f}".format(exp[resnum]) + ',' +
                        "{0:.2f}".format(data["calced"][resnum]) + ','
                    )
                except (IndexError, KeyError):
                    output_csv.write(',,')

            output_csv.write("\n")


class ChemShift_modell_data(object):

    """Class for per model chemical shift data"""
    type_dict = []

    @staticmethod
    def get_type_data(my_type):
        type_data = []

        for model in ChemShift_modell_data.type_dict:
            type_data.append(model[my_type])

        return type_data


class Vec_3D(object):

    """Vector class for calculations"""
    def __init__(self, v):
        self.v = v

    def __getitem__(self, index):
        return self.v[index]

    def __len__(self):
        return len(self.v)

    def __sub__(self, other):
        v = self.v
        u = other.v
        return Vec_3D([u[i] - v[i] for i in range(len(u))])

    def __iadd__(self, other):
        v = self.v
        u = other.v
        return Vec_3D([u[i] + v[i] for i in range(len(u))])

    def __idiv__(self, divisor):
        v = self.v
        return Vec_3D([v[i] / divisor for i in range(len(v))])

    def magnitude(self):
        v = self.v
        return math.sqrt(sum(v[i] * v[i] for i in range(len(v))))

    def normalize(self):
        vmag = self.magnitude()
        v = self.v
        return Vec_3D([v[i] / vmag for i in range(len(v))])

    @classmethod
    def cross(cls, one, other):
        c = [one.v[1] * other.v[2] - one.v[2] * other.v[1],
             one.v[2] * other.v[0] - one.v[0] * other.v[2],
             one.v[0] * other.v[1] - one.v[1] * other.v[0]]

        return cls(c)

    @staticmethod
    def dihedAngle(one, other):
        calc_cos = (
            one.v[0] * other.v[0] +
            one.v[1] * other.v[1] +
            one.v[2] * other.v[2]
            ) / (one.magnitude() * other.magnitude())

        # minor correction due to numerical issues
        if calc_cos > 1:
            calc_cos = 1

        return math.degrees(math.acos(calc_cos))
