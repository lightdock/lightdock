import sys
import pdb


class ForkedPdb(pdb.Pdb):
    """A Pdb subclass that may be used from a forked multiprocessing child.

    Example of use:
    from multiprocessing_pdb import ForkedPDB
    ForkedPdb().set_trace()
    """

    def interaction(self, *args, **kwargs):
        _stdin = sys.stdin
        try:
            sys.stdin = open("/dev/stdin")
            pdb.Pdb.interaction(self, *args, **kwargs)
        finally:
            sys.stdin = _stdin
