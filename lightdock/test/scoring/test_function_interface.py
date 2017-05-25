"""Tests for Scoring Function interface classes"""

from nose.tools import raises
from lightdock.scoring.functions import ScoringFunction, ModelAdapter


class TestScoringFunction:

    def setup(self):
        pass

    def teardown(self):
        pass
    
    @raises(NotImplementedError)
    def test_call_scoring_function_interface(self):
        sf = ScoringFunction()
        sf(None, None, None, None)
        assert False
        

class TestModelAdapter:

    def setup(self):
        pass

    def teardown(self):
        pass
    
    @raises(NotImplementedError)
    def test_create_model_adapter_interface(self):
        adapter = ModelAdapter(None, None)
        assert not adapter
