"""Tests for Coordinates class"""

from nose.tools import assert_almost_equals
from nose.tools import raises
import os
from lightdock.gso.coordinates import Coordinates
from lightdock.gso.coordinates import CoordinatesFileReader
from lightdock.error.lightdock_errors import GSOCoordinatesError


class TestCoordinates:

    def setUp(self):
        self.values_2D = [1.0, 2.0]
        self.values_3D = [1.0, 2.0, 3.0]

    def tearDown(self):
        pass
    
    def test_create_coordinates_2D(self):
        coordinates = Coordinates(self.values_2D)
        
        assert coordinates.dimension == 2
        for i in range(coordinates.dimension):
            assert_almost_equals(self.values_2D[i], coordinates[i])

    def test_create_coordinates_3D(self):
        coordinates = Coordinates(self.values_3D)
        
        assert coordinates.dimension == 3
        for i in range(coordinates.dimension):
            assert_almost_equals(self.values_3D[i], coordinates[i])
            
    def test_check_equal_coordinates_2D(self):
        coordinates1 = Coordinates(self.values_2D)
        coordinates2 = Coordinates(self.values_2D)
        
        assert coordinates1 == coordinates2
        assert not coordinates1 != coordinates2
        
    def test_check_not_equal_coordinates_2D_and_3D(self):
        coordinates1 = Coordinates(self.values_2D)
        coordinates2 = Coordinates(self.values_3D)
        
        assert not coordinates1 == coordinates2
        assert coordinates1 != coordinates2
    
    def test_index_assigment(self):
        coordinates = Coordinates(self.values_2D)
        
        assert_almost_equals(self.values_2D[0], coordinates[0])
        assert_almost_equals(self.values_2D[1], coordinates[1])
        
        coordinates[0] = -1.0
        
        assert_almost_equals(-1.0, coordinates[0])
        assert_almost_equals(self.values_2D[1], coordinates[1])
    
    def test_clone_coordinates(self):
        coordinates1 = Coordinates(self.values_2D)
        coordinates2 = coordinates1.clone()
        
        assert coordinates1 == coordinates2
        
        coordinates2[0] = -1.0
        
        assert coordinates1 != coordinates2
    
    def test_coordinates_addition(self):
        coordinates1 = Coordinates(self.values_2D)
        coordinates2 = Coordinates(self.values_2D)
        
        expected = Coordinates([2.0, 4.0])
        
        assert expected == coordinates1 + coordinates2

    def test_coordinates_subtraction(self):
        coordinates1 = Coordinates(self.values_2D)
        coordinates2 = Coordinates(self.values_2D)
        
        expected = Coordinates([0.0, 0.0])
        
        assert expected == coordinates1 - coordinates2

    def test_coordinates_addition_and_assigment(self):
        coordinates1 = Coordinates(self.values_2D)
        coordinates2 = Coordinates(self.values_2D)
        
        expected = Coordinates([2.0, 4.0])
        
        coordinates1 += coordinates2
        
        assert expected == coordinates1

    def test_coordinates_subtraction_and_assigment(self):
        coordinates1 = Coordinates(self.values_2D)
        coordinates2 = Coordinates(self.values_2D)
        
        expected = Coordinates([0.0, 0.0])
        
        coordinates1 -= coordinates2
        
        assert expected == coordinates1
    
    def test_norm(self):
        coordinates = Coordinates(self.values_2D)
        
        assert_almost_equals(2.236067977, coordinates.norm())
        
    def test_distance_same_coordinate(self):
        coordinates = Coordinates(self.values_2D)
        
        assert_almost_equals(0.0, coordinates.distance(coordinates))
        
    def test_distance_different_coordinates(self):
        coordinates1 = Coordinates([0., 0., 0.])
        coordinates2 = Coordinates([20., 0., 21.])
        
        assert_almost_equals(29.0, coordinates1.distance(coordinates2))

    def test_distance2_same_coordinate(self):
        coordinates = Coordinates(self.values_2D)
        
        assert_almost_equals(0.0, coordinates.distance2(coordinates))
    
    def test_sum_of_squares(self):
        coordinates = Coordinates(self.values_2D)
        
        assert_almost_equals(5.0, coordinates.sum_of_squares())
        
    def test_distance2_different_coordinates(self):
        coordinates1 = Coordinates(self.values_2D)
        coordinates2 = Coordinates([2.0, 3.0])
        
        assert_almost_equals(2.0, coordinates1.distance2(coordinates2))
        
    def test_multiplication_and_assigment(self):
        coordinates = Coordinates(self.values_3D)
        
        expected = Coordinates([-3.0, -6.0, -9.0])
        
        coordinates *= -3.0
        
        assert expected == coordinates
 
    def test_multiplication(self):
        coordinates = Coordinates(self.values_3D)
        
        expected = Coordinates([-3.0, -6.0, -9.0])
        
        assert expected == (coordinates * -3.0)
    
    def test_move_different_coordinates(self):
        coordinates1 = Coordinates(self.values_2D)
        coordinates2 = Coordinates([0.0, 1.0])
        
        expected = Coordinates([-1.12132034356, -0.12132034356])
        
        assert expected == coordinates1.move(coordinates2, 3.0)
    
    def test_move_same_coordinate(self):
        coordinates1 = Coordinates(self.values_2D)
        
        assert coordinates1 == coordinates1.move(coordinates1)


class TestCoordinatesFileReader:

    def setUp(self):
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
    
    def tearDown(self):
        pass
    
    def test_read_coordinates_from_file(self):
        reader = CoordinatesFileReader(2)
        coordinates = reader.get_coordinates_from_file(self.golden_data_path + 'initial_positions.txt')
        
        assert 50 == len(coordinates)
        assert "(0.745916, -0.92056)" == str(coordinates[0])
        assert "(-2.29363, -0.229427)" == str(coordinates[9])
        assert "(0.617171, -2.85014)" == str(coordinates[-1])
    
    @raises(GSOCoordinatesError)
    def test_read_coordinates_from_file_with_errors(self):
        reader = CoordinatesFileReader(2)
        coordinates = reader.get_coordinates_from_file(self.golden_data_path + 'initial_positions_with_error.txt')
        
        assert len(coordinates)
        
    @raises(GSOCoordinatesError)
    def test_read_coordinates_from_file_with_error_in_column(self):
        reader = CoordinatesFileReader(2)
        coordinates = reader.get_coordinates_from_file(self.golden_data_path + 'initial_positions_with_wrong_column.txt')
        
        assert len(coordinates)
        
    @raises(GSOCoordinatesError)
    def test_read_coordinates_from_file_no_file(self):
        reader = CoordinatesFileReader(2)
        coordinates = reader.get_coordinates_from_file(self.golden_data_path + 'no_file.txt')
        
        assert len(coordinates)
