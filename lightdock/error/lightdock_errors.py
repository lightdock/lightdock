"""Custom error classes"""


class LightDockError(Exception):
    """LightDock exception base class"""
    def __init__(self, cause):
        self.cause = cause

    def __str__(self):
        representation = "[%s] %s" % (self.__class__.__name__, self.cause)
        return representation


class LightDockWarning(LightDockError):
    """Custom error class intented only for warnings to be notified, not to fail"""
    pass


class RandomNumberError(LightDockError):
    """Custom RandomNumber exception"""
    pass


class GSOError(LightDockError):
    """Custom GSO exception"""
    pass


class GSOParameteresError(GSOError):
    """Custom GSOParameteres exception"""
    pass


class GSOCoordinatesError(GSOError):
    """Custom error for CoordinatesFileReader class"""
    pass


class StructureError(LightDockError):
    """General structure error"""
    pass


class BackboneError(StructureError):
    """General structure error"""
    pass


class SideChainError(StructureError):
    """General structure error"""
    pass


class ResidueNonStandardError(StructureError):
    """General structure error"""
    pass


class AtomError(StructureError):
    """Atom error exception"""
    pass


class MinimumVolumeEllipsoidError(StructureError):
    """MinimumVolumeEllipsoid exception"""
    pass


class PDBParsingError(LightDockError):
    """PDB parser error"""
    pass


class PDBParsingWarning(LightDockWarning):
    """PDB parser warning"""
    pass


class PotentialsParsingError(LightDockError):
    """Reading potential file error"""
    pass


class ScoringFunctionError(LightDockError):
    """Error in the scoring function drivers"""
    pass


class NotSupportedInScoringError(LightDockError):
    """Error to be raised when an atom or residue type is
    not supported by the scoring function"""
    pass


class NormalModesCalculationError(LightDockError):
    """Error in normal modes calculation"""
    pass
