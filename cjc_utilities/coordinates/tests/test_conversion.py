"""
Tests for cross-section extraction

:author: Calum Chamberlain
:date: 23 March 2018
"""

import unittest

from coordinates import Location, Geographic


class TestExtraction(unittest.TestCase):
    """
    Base test class for coordinate conversions
    """
    def test_simple_convert_southern(self):
        origin = Geographic(latitude=-42.2, longitude=150, depth=0)
        location = Geographic(latitude=-40, longitude=170, depth=-10)
        xyz_location = location.to_xyz(origin, 0, 90)
        self.assertEqual(xyz_location.z, -10)
        self.assertAlmostEqual(xyz_location.y, 244, delta=10)
        self.assertAlmostEqual(xyz_location.x, 1675, delta=10)
        # Check backwards
        location_back = xyz_location.to_geographic()
        self.assertEqual(location_back, location)

    def test_simple_convert_northern(self):
        origin = Geographic(latitude=20, longitude=50, depth=3)
        location = Geographic(latitude=40, longitude=70, depth=-10)
        xyz_location = location.to_xyz(origin, 0, 90)
        self.assertEqual(xyz_location.z, -13)
        self.assertAlmostEqual(xyz_location.y, 2224, delta=10)
        self.assertAlmostEqual(xyz_location.x, 1922, delta=10)
        # Check backwards
        location_back = xyz_location.to_geographic()
        self.assertEqual(location_back, location)

    def test_simple_convert_cross_hemisphere(self):
        origin = Geographic(latitude=10, longitude=50, depth=3)
        location = Geographic(latitude=-10, longitude=70, depth=-10)
        xyz_location = location.to_xyz(origin, 0, 90)
        self.assertEqual(xyz_location.z, -13)
        self.assertAlmostEqual(xyz_location.y, -2222, delta=10)
        self.assertAlmostEqual(xyz_location.x, 2222, delta=10)
        # Check backwards
        location_back = xyz_location.to_geographic()
        self.assertEqual(location_back, location)

    def test_striking_convert_southern(self):
        origin = Geographic(latitude=-42.2, longitude=150, depth=0)
        location = Geographic(latitude=-40, longitude=170, depth=-10)
        xyz_location = location.to_xyz(origin, 10, 90)
        self.assertEqual(xyz_location.z, -10)
        # Check distance is at it should be
        striking_distance = (xyz_location.x ** 2 + xyz_location.y ** 2) ** 0.5
        xyz_ns = location.to_xyz(origin, 0, 90)
        non_striking_distance = (xyz_ns.x ** 2 + xyz_ns.y ** 2) ** 0.5
        self.assertAlmostEqual(striking_distance, non_striking_distance,
                               delta=10)
        # Check location against pre-calculated result
        self.assertAlmostEqual(xyz_location.y, 530, delta=10)
        self.assertAlmostEqual(xyz_location.x, 1607, delta=10)
        # Check backwards
        location_back = xyz_location.to_geographic()
        self.assertEqual(location_back, location)

    def test_dipping_convert_southern(self):
        origin = Geographic(latitude=-42.2, longitude=150, depth=0)
        location = Geographic(latitude=-42.2, longitude=150, depth=-10)
        xyz_location = location.to_xyz(origin, 0, 50)
        self.assertAlmostEqual(xyz_location.z, -7.6, delta=0.5)
        self.assertAlmostEqual(xyz_location.y, 0, delta=10)
        self.assertAlmostEqual(xyz_location.x, -6.4, delta=0.5)
        # Check backwards
        location_back = xyz_location.to_geographic()
        self.assertEqual(location_back, location)

    def test_dipping_convert_southern_offset_north(self):
        origin = Geographic(latitude=-42.2, longitude=150, depth=0)
        location = Geographic(latitude=-41.2, longitude=150, depth=-10)
        xyz_location = location.to_xyz(origin, 0, 50)
        self.assertAlmostEqual(xyz_location.z, -7.6, delta=0.5)
        self.assertAlmostEqual(xyz_location.y, 111, delta=10)
        self.assertAlmostEqual(xyz_location.x, -6.4, delta=0.5)
        # Check backwards
        location_back = xyz_location.to_geographic()
        self.assertEqual(location_back, location)

    def test_dipping_convert_southern_offset_west(self):
        origin = Geographic(latitude=-42.2, longitude=150, depth=0)
        location = Geographic(latitude=-42.2, longitude=152, depth=-10)
        xyz_location = location.to_xyz(origin, 0, 80)
        self.assertAlmostEqual(xyz_location.y, 0, delta=10)
        self.assertAlmostEqual(xyz_location.x, 164, delta=10)
        self.assertAlmostEqual(xyz_location.z, -38, delta=5)
        # Check backwards
        location_back = xyz_location.to_geographic()
        self.assertEqual(location_back, location)

if __name__ == '__main__':
    unittest.main()