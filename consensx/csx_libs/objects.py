import math


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
