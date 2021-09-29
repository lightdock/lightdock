"""Tests for Chain class"""

from lightdock.structure.chain import Chain
from lightdock.structure.residue import Residue


class TestChain:
    def test_create_empty_chain(self):
        chain = Chain()
        assert chain.cid == "" and chain.peptide and len(chain.residues) == 0

    def test_create_chain(self):
        chain = Chain("A", [Residue("ALA", 1), Residue("MET", 2)])
        assert chain.cid == "A" and len(chain.residues) == 2

    def test_clone(self):
        chain1 = Chain("A", [Residue("ALA", 1), Residue("MET", 2)])
        chain2 = chain1.clone()

        assert (
            len(chain1.residues) == len(chain2.residues)
            and chain1.residues[0].name == chain2.residues[0].name
        )

        chain1.residues[0].name = "MET"
        assert chain1.residues[0].name != chain2.residues[0].name

    def test__to_string(self):
        chain = Chain("A", [Residue("ALA", 1), Residue("MET", 2)])

        assert str(chain) == "[Chain A]\nALA.1\nMET.2"
