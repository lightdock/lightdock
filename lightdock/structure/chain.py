"""Module to package a complex chain"""


class Chain(object):
    """Represents a complex chain"""
    def __init__(self, cid='', residues=None, only_peptide=True):
        """Creates a new chain"""
        self.cid = cid
        if residues:
            self.residues = residues
        else:
            self.residues = []
        self.peptide = only_peptide

    def clone(self):
        """Creates a copy of the current chain"""
        return Chain(self.cid,
                     [residue.clone() for residue in self.residues],
                     self.peptide)

    def __str__(self):
        return "[Chain %s]\n%s" % (self.cid, '\n'.join([str(residue) for residue in self.residues]))
