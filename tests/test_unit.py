import unittest
from typing import Any, Iterable, Optional, Tuple

from experimentdata import unit

from tests.util import QuantityTestCase


class FunctionTestCase(QuantityTestCase):
    def test_scale(self):
        self.assertQuantity(
            unit.Quantity(1000, unit.registry.ohm).to_compact(),
            unit.Quantity(1, unit.registry.kiloohm)
        )

        self.assertQuantity(
            unit.Quantity(1000000, unit.registry.ohm).to_compact(),
            unit.Quantity(1, unit.registry.megaohm)
        )

        self.assertQuantity(
            unit.Quantity(1000000, unit.registry.ohm).to(unit.registry.kiloohm),
            unit.Quantity(1000, unit.registry.kiloohm)
        )

    def test_str_mixed_unit_per(self):
        self.assertEqual('1 Ω/m', str(unit.Quantity(1, unit.registry.ohm / unit.registry.meter)))
        self.assertEqual('1 kΩ/m', str(unit.Quantity(1, unit.registry.ohm / unit.registry.millimeter)))
        self.assertEqual('1 kΩ/m', str(unit.Quantity(1, unit.registry.kiloohm / unit.registry.meter)))

    def test_str_mixed_unit_product(self):
        self.assertEqual('1 m Ω', str(unit.Quantity(1, unit.registry.ohm * unit.registry.meter)))
        self.assertEqual('1 kΩ m', str(unit.Quantity(1, unit.registry.kiloohm * unit.registry.meter)))
        self.assertEqual('1 m mΩ', str(unit.Quantity(1, unit.registry.ohm * unit.registry.millimeter)))


class ParseTestCase(QuantityTestCase):
    def _testParseList(self, test_list: Iterable[Tuple[unit.TParseQuantity, Optional[unit.TParseUnit], unit.Quantity,
                                                       Optional[int]]], mag_round: Optional[int] = None):
        for test_in, test_unit, expected_qty, places in test_list:
            with self.subTest(f"Input: {test_in!r}, Unit: {test_unit or 'no unit'} = {expected_qty!s} "
                              f"(places: {places}, round: {mag_round})"):
                self.assertQuantity(unit.parse(test_in, test_unit, mag_round), expected_qty, places)

    def test_parse_unit(self):
        self.assertEqual(unit.registry.degK, unit.parse_unit('K'))
        self.assertEqual(unit.registry.degK, unit.parse_unit('kelvin'))
        self.assertEqual(unit.registry.degK, unit.parse_unit(unit.registry.kelvin))

        self.assertEqual(unit.registry.degC, unit.parse_unit('degC'))
        self.assertEqual(unit.registry.degC, unit.parse_unit('°C'))

        self.assertEqual(unit.registry.meter, unit.parse_unit('m'))
        self.assertEqual(unit.registry.millimeter, unit.parse_unit('mm'))
        self.assertEqual(unit.registry.kilometer, unit.parse_unit('km'))

    def test_parse_dimensionless(self):
        self._testParseList([
            ('0.001', None, unit.Quantity(0.0, unit.dimensionless), None),
            (0.001, None, unit.Quantity(0.0, unit.dimensionless), None),
            (unit.Quantity(0.001, unit.dimensionless), None, unit.Quantity(0.0, unit.dimensionless), None),

            ('0.5001', None, unit.Quantity(1.0, unit.dimensionless), None),
            (0.5001, None, unit.Quantity(1.0, unit.dimensionless), None),
            (unit.Quantity(0.5001, unit.dimensionless), None, unit.Quantity(1.0, unit.dimensionless), None),

            ('1.499', None, unit.Quantity(1.0, unit.dimensionless), None),
            (1.499, None, unit.Quantity(1.0, unit.dimensionless), None),
            (unit.Quantity(1.49, unit.dimensionless), None, unit.Quantity(1.0, unit.dimensionless), None)
        ], mag_round=0)

        self._testParseList([
            ('0.0001', None, unit.Quantity(0.0, unit.dimensionless), None),
            (0.0001, None, unit.Quantity(0.0, unit.dimensionless), None),
            (unit.Quantity(0.0001, unit.dimensionless), None, unit.Quantity(0.0, unit.dimensionless), None),

            ('0.0005001', None, unit.Quantity(0.001, unit.dimensionless), None),
            (0.0005001, None, unit.Quantity(0.001, unit.dimensionless), None),
            (unit.Quantity(0.0005001, unit.dimensionless), None, unit.Quantity(0.001, unit.dimensionless), None),

            ('1.4999', None, unit.Quantity(1.5, unit.dimensionless), None),
            (1.4999, None, unit.Quantity(1.5, unit.dimensionless), None),
            (unit.Quantity(1.4999, unit.dimensionless), None, unit.Quantity(1.5, unit.dimensionless), None)
        ], mag_round=3)

    def test_parse_round(self):
        self._testParseList([
            ('0.001', None, unit.Quantity(0.001, unit.dimensionless), None),
            (0.001, None, unit.Quantity(0.001, unit.dimensionless), None),

            ('1', None, unit.Quantity(1.0, unit.dimensionless), None),
            ('1.0', None, unit.Quantity(1.0, unit.dimensionless), None),
            (1, None, unit.Quantity(1.0, unit.dimensionless), None),
            (1.0, None, unit.Quantity(1.0, unit.dimensionless), None),

            ('1000', None, unit.Quantity(1000.0, unit.dimensionless), None),
            ('1000.0', None, unit.Quantity(1000.0, unit.dimensionless), None),
            (1000, None, unit.Quantity(1000.0, unit.dimensionless), None),
            (1000.0, None, unit.Quantity(1000.0, unit.dimensionless), None)
        ])

    def test_parse_provide_unit(self):
        self._testParseList([
            ('0.001', unit.registry.degC, unit.Quantity(0.001, unit.registry.degC), None),
            (0.001, unit.registry.degC, unit.Quantity(0.001, unit.registry.degC), None),

            ('1', unit.registry.degC, unit.Quantity(1.0, unit.registry.degC), None),
            ('1.0', unit.registry.degC, unit.Quantity(1.0, unit.registry.degC), None),
            (1, unit.registry.degC, unit.Quantity(1.0, unit.registry.degC), None),
            (1.0, unit.registry.degC, unit.Quantity(1.0, unit.registry.degC), None),

            ('1000', unit.registry.degC, unit.Quantity(1000.0, unit.registry.degC), None),
            ('1000.0', unit.registry.degC, unit.Quantity(1000.0, unit.registry.degC), None),
            (1000, unit.registry.degC, unit.Quantity(1000.0, unit.registry.degC), None),
            (1000.0, unit.registry.degC, unit.Quantity(1000.0, unit.registry.degC), None),

            ('1.0°C', None, unit.Quantity(1.0, unit.registry.degC), None),
            ('1.0°C', unit.registry.degC, unit.Quantity(1.0, unit.registry.degC), None),

            (unit.Quantity(1.0, unit.registry.degC), None, unit.Quantity(1.0, unit.registry.degC), None),
            (unit.Quantity(1.0, unit.registry.degC), unit.registry.degC, unit.Quantity(1.0, unit.registry.degC), None)
        ])

    def test_parse_mixed_scale(self):
        self._testParseList([
            ('1m', unit.registry.millimeter, unit.Quantity(1000.0, unit.registry.millimeter), None),
            (unit.Quantity(1.0, unit.registry.meter), unit.registry.millimeter,
             unit.Quantity(1000.0, unit.registry.millimeter), None)
        ])

    def test_parse_mixed_unit(self):
        self._testParseList([
            ('1 Ω m', None, unit.Quantity(1.0, unit.registry.ohm * unit.registry.meter), None),
            ('1 Ω/m', None, unit.Quantity(1.0, unit.registry.ohm / unit.registry.meter), None)
        ])


class PrintingTestCase(unittest.TestCase):
    def assertStr(self, test_set):
        for qty, qty_str in test_set:
            with self.subTest(f"{qty!r} == {qty_str!r}"):
                self.assertEqual(qty_str, str(qty))

    def test_print_single(self):
        self.assertStr([
            (unit.Quantity(0.001, unit.registry.meter), '1 mm'),
            (unit.Quantity(1, unit.registry.millimeter), '1 mm'),
            (unit.Quantity(0.01, unit.registry.meter), '10 mm'),
            (unit.Quantity(10, unit.registry.millimeter), '10 mm'),
            (unit.Quantity(0.1, unit.registry.meter), '100 mm'),
            (unit.Quantity(100, unit.registry.millimeter), '100 mm'),
            (unit.Quantity(1, unit.registry.meter), '1 m'),
            (unit.Quantity(1000, unit.registry.millimeter), '1 m'),

            (unit.Quantity(10, unit.registry.meter), '10 m'),
            (unit.Quantity(100, unit.registry.meter), '100 m'),
            (unit.Quantity(1000, unit.registry.meter), '1 km')
        ])

    def test_print_limit(self):
        self.assertEqual(str(unit.Quantity(1000000, unit.registry.meter)), '1000 km')

    def test_print_dimensionless(self):
        self.assertStr([
            (unit.Quantity(1, unit.registry.ppb), '1 ppb'),
            (unit.Quantity(10, unit.registry.ppb), '10 ppb'),
            (unit.Quantity(100, unit.registry.ppb), '100 ppb'),
            (unit.Quantity(1000, unit.registry.ppb), '1 ppm'),
            (unit.Quantity(10000, unit.registry.ppb), '10 ppm'),
            (unit.Quantity(100000, unit.registry.ppb), '100 ppm'),
            (unit.Quantity(1000000, unit.registry.ppb), '0.1%'),
            (unit.Quantity(10000000, unit.registry.ppb), '1%'),
            (unit.Quantity(100000000, unit.registry.ppb), '10%'),
            (unit.Quantity(1000000000, unit.registry.ppb), '100%'),

            (unit.Quantity(0.001, unit.registry.ppm), '1 ppb'),
            (unit.Quantity(0.01, unit.registry.ppm), '10 ppb'),
            (unit.Quantity(0.1, unit.registry.ppm), '100 ppb'),
            (unit.Quantity(1, unit.registry.ppm), '1 ppm'),
            (unit.Quantity(10, unit.registry.ppm), '10 ppm'),
            (unit.Quantity(100, unit.registry.ppm), '100 ppm'),
            (unit.Quantity(1000, unit.registry.ppm), '0.1%'),
            (unit.Quantity(10000, unit.registry.ppm), '1%'),
            (unit.Quantity(100000, unit.registry.ppm), '10%'),
            (unit.Quantity(1000000, unit.registry.ppm), '100%'),

            (unit.Quantity(0.0000001, unit.registry.percent), '1 ppb'),
            (unit.Quantity(0.000001, unit.registry.percent), '10 ppb'),
            (unit.Quantity(0.00001, unit.registry.percent), '100 ppb'),
            (unit.Quantity(0.0001, unit.registry.percent), '1 ppm'),
            (unit.Quantity(0.001, unit.registry.percent), '10 ppm'),
            (unit.Quantity(0.01, unit.registry.percent), '100 ppm'),
            (unit.Quantity(0.1, unit.registry.percent), '0.1%'),
            (unit.Quantity(1, unit.registry.percent), '1%'),
            (unit.Quantity(10, unit.registry.percent), '10%'),
            (unit.Quantity(100, unit.registry.percent), '100%')
        ])


if __name__ == '__main__':
    unittest.main()
