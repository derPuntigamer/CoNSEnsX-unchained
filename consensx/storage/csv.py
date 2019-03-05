class CSVBuffer(object):
    """Class which stores data for values.CSV"""

    def __init__(self, my_path):
        self.working_dir = my_path
        self.max_resnum = -1
        self.min_resnum = 100000
        self.csv_data = []

    def add_data(self, data):
        self.csv_data.append(data)

    def write_csv(self):
        filename = self.working_dir + "values.csv"
        output_csv = open(filename, 'w')
        output_csv.write(',')
        for data in self.csv_data:
            output_csv.write(data["name"] + " EXP, " + data["name"] + " CALC,")
        output_csv.write("\n")

        for resnum in range(self.min_resnum, self.max_resnum + 1):
            output_csv.write(str(resnum) + ',')
            for data in self.csv_data:
                exp = {}

                for i in data["experimental"]:
                    exp[i.resnum] = i.value

                try:
                    output_csv.write(
                        "{0:.2f}".format(exp[resnum]) + ',' +
                        "{0:.2f}".format(data["calced"][resnum]) + ','
                    )
                except (IndexError, KeyError):
                    output_csv.write(',,')

            output_csv.write("\n")
