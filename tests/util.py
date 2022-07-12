from typing import Any, Optional
from unittest import TestCase

from experimentdata.unit import Quantity, Unit


class QuantityTestCase(TestCase):
    def assertQuantity(self, x: Any, expected: Quantity, places: Optional[int] = None,
                       magnitude_unit: Optional[Unit] = None):
        self.assertIsInstance(x, Quantity, 'Not an instance of Quantity')

        if magnitude_unit is not None:
            self.assertTrue(x.is_compatible_with(expected.units), 'Incompatible units')
            x_mag = x.m_as(magnitude_unit)
            expected_mag = expected.m_as(magnitude_unit)
        else:
            self.assertEqual(x.units, expected.units, 'Units do not match')
            x_mag = x.magnitude
            expected_mag = expected.magnitude

        if places is None:
            self.assertEqual(x_mag, expected_mag, 'Magnitude does not match')
        else:
            self.assertAlmostEqual(x_mag, expected_mag, places, 'Magnitude does not match')
