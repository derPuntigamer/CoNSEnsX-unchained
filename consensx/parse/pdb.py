import prody
import pickle


class Pdb():
    def __init__(self, PDB_file, my_path, model_count):
        self.model_count = model_count
        self.is_fitted = False

        prody.confProDy(verbosity="info")

        self.atomgroup = prody.parsePDB(PDB_file, ter=True)
        num_coordsets = self.atomgroup.numCoordsets()

        # check for discarded models
        for model_num in range(num_coordsets):
            self.atomgroup.setACSIndex(model_num)

            if self.atomgroup.getCoords()[0][0] == float(0):
                # TODO raise Exception
                print("DISCARDED MODEL FOUND")
                return False

        try:
            self.elements = self.atomgroup.getElements()
        except AttributeError:
            print("ERROR -> PDB parsing failed. Please check your PDB file!")
            raise SystemExit

        self.names = self.atomgroup.getNames()
        self.coordsets = self.atomgroup.numCoordsets()

        PDB_model_path = my_path + "/PDB_model.pickle"
        pickle.dump(self, open(PDB_model_path, 'wb'))
