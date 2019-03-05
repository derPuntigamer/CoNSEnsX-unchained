import os


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


def check_3rd_party(install_dir):
    config_file_name = install_dir + "/.config"
    if os.path.isfile(config_file_name):
        ThirdParty.get_thirdparty(config_file_name)
    else:
        init_conf = (
            "# CoNSEnsX config file\n" +
            "# Please provide full paths\n" +
            "pales=''\n" + "shiftx=''\n" +
            "prideDB=''\n" + "prideNMR=''"
        )

        init_conf_file = open(".config", "w")
        init_conf_file.write(init_conf)
        print("Please edit '.config' in the CoNSEnsX install directory")
        raise SystemExit
