"""Tests for Scoring Function interface classes"""

import pytest
from lightdock.scoring.functions import ScoringFunction, ModelAdapter


class TestScoringFunction:
    def test_call_scoring_function_interface(self):
        with pytest.raises(NotImplementedError):
            sf = ScoringFunction()
            sf(None, None, None, None)
            assert False


class TestModelAdapter:
    def test_create_model_adapter_interface(self):
        with pytest.raises(NotImplementedError):
            adapter = ModelAdapter(None, None)
            assert not adapter
