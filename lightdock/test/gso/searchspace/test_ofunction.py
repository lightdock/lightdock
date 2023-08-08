"""Tests for ObjectiveFunction interface"""

import pytest
from lightdock.gso.searchspace.ofunction import ObjectiveFunction


class TestObjectiveFunction:
    def test_call(self):
        with pytest.raises(NotImplementedError):
            function = ObjectiveFunction()
            function(coordinates=[])
