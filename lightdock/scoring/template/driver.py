"""Template scoring function"""


from lightdock.structure.model import DockingModel
from lightdock.scoring.functions import ModelAdapter, ScoringFunction


class TemplateAdapter(ModelAdapter):
    """Adapts a given Complex to a DockingModel object suitable for this
    'Template' scoring function.
    """

    def _get_docking_model(self, molecule, restraints):
        """Builds a suitable docking model for this scoring function"""
        # In model_objects we can store any coordinates object (atoms, beans, etc.)
        model_objects = []
        for residue in molecule.residues:
            for rec_atom in residue.atoms:
                model_objects.append(rec_atom)
        try:
            return DockingModel(model_objects, molecule.copy_coordinates(), restraints, n_modes=molecule.n_modes.copy())
        except AttributeError:
            return DockingModel(model_objects, molecule.copy_coordinates(), restraints)


class TemplateScoringFunction(ScoringFunction):
    """Implements the 'Template' scoring function"""
    def __init__(self, weight=1.0):
        super(TemplateScoringFunction, self).__init__(weight)

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        """This is the function which evaluates the current model. receptor.objects and ligand.objects contain the
        information loaded in the _get_docking_model() from the adapter. receptor_coordinates and ligand_coordinates
        are NumPy arrays containing the (x,y,z) coordinates of the objects.
        """
        return 0.0


# Needed to dynamically load the scoring functions from command line
DefinedScoringFunction = TemplateScoringFunction
DefinedModelAdapter = TemplateAdapter
