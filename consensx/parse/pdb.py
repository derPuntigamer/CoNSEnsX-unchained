class Pdb():
    def __init__(self, atomgroup, model_count, is_fitted=False):
        self.atomgroup = atomgroup
        self.model_count = model_count
        self.is_fitted = is_fitted

        try:
            self.elements = atomgroup.getElements()
        except AttributeError:
            print("ERROR -> PDB parsing failed. Please check your PDB file!")
            raise SystemExit

        self.names = atomgroup.getNames()
        self.coordsets = atomgroup.numCoordsets()
