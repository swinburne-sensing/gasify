import unittest
import typing

from experimentdata import unit


class _QuantityTestCase(unittest.TestCase):
    __test__ = False

    MIXED_TEST_RANGE = [
        (1, 1.0),
        (10, 10.0),
        (100, 100.0),
        (1000, 1000.0),
        (0.001, 0.001),
        (0.01, 0.01),
        (0.1, 0.1),
        (1.0, 1.0),
        (10.0, 10.0),
        (100.0, 100.0),
        (1000.0, 1000.0)
    ]

    MIXED_TEST_RANGE_STR = [(str(x), y if y is not None else x) for x, y in MIXED_TEST_RANGE]

    def assertQuantity(self, x: typing.Any, magnitude: typing.Union[int, float],
                       expected_unit: typing.Optional[unit.Unit] = None,
                       magnitude_unit: typing.Optional[unit.Unit] = None, precision: typing.Optional[int] = None):
        expected_unit = expected_unit or unit.dimensionless
        magnitude_unit = magnitude_unit or expected_unit

        self.assertIsInstance(x, unit.Quantity, 'Not an instance of Quantity')
        self.assertIs(x._REGISTRY, unit.registry, 'Quantity not part of registry')
        self.assertEqual(expected_unit, x.units, 'Quantity does not match expected units')
        self.assertTrue(x.is_compatible_with(magnitude_unit), 'Incompatible with magnitude check')
        self.assertAlmostEqual(x.m_as(magnitude_unit), magnitude, precision, 'Magnitude does not match to expected '
                                                                             'precision')

    def assertRange(self, unit_iterable: typing.Iterable[typing.Tuple[typing.Any, unit.Unit]], value_iterable,
                    precision: typing.Optional[int] = None, mag_round: typing.Optional[int] = None):
        """

        :param unit_iterable:
        :param value_iterable:
        :param precision:
        :param mag_round:
        :return:
        """
        for u_parse, u_expected in unit_iterable:
            for v in value_iterable:
                if len(v) == 2:
                    v_parse, v_magnitude = v
                    v_magnitude_unit = None
                elif len(v) == 3:
                    v_parse, v_magnitude, v_magnitude_unit = v
                else:
                    raise NotImplementedError('Incorrect test specification')

                if u_parse is not None:
                    with self.subTest(f"To unit {u_parse!r} from value {v_parse!r} "):
                        self.assertQuantity(unit.parse(v_parse, u_parse, mag_round), v_magnitude, u_expected,
                                            v_magnitude_unit, precision)
                else:
                    with self.subTest(f"From value {v_parse!r}"):
                        self.assertQuantity(unit.parse(v_parse, mag_round=mag_round), v_magnitude, u_expected,
                                            v_magnitude_unit, precision)


class UnmodifiedTestCase(_QuantityTestCase):
    def test_ohm(self):
        self.assertQuantity(unit.Quantity(1, unit.registry.ohm), 1, unit.registry.ohm)
        self.assertQuantity(unit.Quantity(1000, unit.registry.ohm), 1000, unit.registry.ohm)
        self.assertQuantity(unit.Quantity(1000, unit.registry.ohm).to_compact(), 1, unit.registry.kiloohm)
        self.assertQuantity(unit.Quantity(1000000, unit.registry.ohm).to_compact(), 1, unit.registry.megaohm)
        self.assertQuantity(unit.Quantity(1000000, unit.registry.ohm).to(unit.registry.kiloohm), 1000,
                            unit.registry.kiloohm)

    def test_mixed_str_per(self):
        q = unit.Quantity(1, unit.registry.ohm / unit.registry.millimeter)

        self.assertQuantity(q, 1, unit.registry.ohm / unit.registry.millimeter)

        self.assertEqual('1 Ω/mm', str(q))

    def test_mixed_str_product(self):
        q = unit.Quantity(1, unit.registry.ohm * unit.registry.meter)

        self.assertQuantity(q, 1, unit.registry.ohm * unit.registry.meter)

        self.assertEqual('1 m Ω', str(q))


class ParseTestCase(_QuantityTestCase):
    def test_simple(self):
        self.assertRange(
            [
                (None, unit.dimensionless),
                ('K', unit.registry.kelvin),
                ('kelvin', unit.registry.kelvin),
                (unit.registry.kelvin, unit.registry.kelvin),
                ('ohm * meter', unit.registry.ohm * unit.registry.meter),
                ('Ω * m', unit.registry.ohm * unit.registry.meter),
                ('Ω*m', unit.registry.ohm * unit.registry.meter),
                ('Ω m', unit.registry.ohm * unit.registry.meter),
                ('ohm / millimeter', unit.registry.ohm / unit.registry.millimeter),
                ('Ω / mm', unit.registry.ohm / unit.registry.millimeter),
                ('Ω/mm', unit.registry.ohm / unit.registry.millimeter),
                (unit.registry.ohm * unit.registry.meter, unit.registry.ohm * unit.registry.meter)
            ],
            self.MIXED_TEST_RANGE
        )

    def test_simple_str(self):
        self.assertRange(
            [
                (None, unit.dimensionless),
                ('K', unit.registry.kelvin),
                ('kelvin', unit.registry.kelvin),
                (unit.registry.kelvin, unit.registry.kelvin)
            ],
            self.MIXED_TEST_RANGE_STR
        )

    def test_simple_unit_str(self):
        self.assertRange(
            [
                (None, unit.registry.kelvin)
            ],
            [(f"{x} K", y) for x, y in self.MIXED_TEST_RANGE_STR]
        )

        self.assertRange(
            [
                (None, unit.registry.kelvin)
            ],
            [(f"{x}K", y) for x, y in self.MIXED_TEST_RANGE_STR]
        )

    def test_simple_expression(self):
        self.assertRange(
            [
                (None, unit.dimensionless),
                ('K', unit.registry.kelvin),
                ('kelvin', unit.registry.kelvin),
                (unit.registry.kelvin, unit.registry.kelvin)
            ], [
                ('1×10²', 100.0),
                ('1×10³', 1000.0)
            ]
        )

    def test_simple_quantity(self):
        self.assertRange(
            [
                (None, unit.registry.kelvin),
                (unit.registry.kelvin, unit.registry.kelvin)
            ],
            [(unit.Quantity(x, unit.registry.kelvin), y) for x, y in self.MIXED_TEST_RANGE]
        )

    def test_simple_round(self):
        self.assertRange(
            [
                (unit.registry.kelvin, unit.registry.kelvin)
            ], [
                (1.0, 1.0),
                (1.01, 1.0),
                (1.049, 1.0),
                (1.05, 1.1),
                (1.99, 2.0)
            ], mag_round=1
        )

        self.assertRange(
            [
                (unit.registry.kelvin, unit.registry.kelvin)
            ], [
                (1.00, 1.0),
                (1.001, 1.0),
                (1.0049, 1.0),
                (1.0051, 1.01),
                (1.999, 2.0)
            ], mag_round=2
        )

        self.assertRange(
            [
                (unit.registry.kelvin, unit.registry.kelvin)
            ], [
                (1.0001, 1.0),
                (1.001, 1.001),
                (1.0049, 1.005),
                (1.0051, 1.005),
                (1.9999, 2.0)
            ], mag_round=3
        )


class ParseDimensionlessTestCase(_QuantityTestCase):
    def test_percent(self):
        self.assertRange(
            [
                ('%', unit.registry.percent),
                ('percent', unit.registry.percent),
                (unit.registry.percent, unit.registry.percent),
            ],
            [(x, y / 100.0, unit.dimensionless) for x, y in self.MIXED_TEST_RANGE]
        )

    def test_ppm(self):
        self.assertRange(
            [
                ('ppm', unit.registry.ppm),
                ('ppm', unit.registry.ppm)
            ],
            [(x, y / 1000000.0, unit.dimensionless) for x, y in self.MIXED_TEST_RANGE]
        )

    def test_ppb(self):
        self.assertRange(
            [
                ('ppb', unit.registry.ppb),
                (unit.registry.ppb, unit.registry.ppb)
            ],
            [(x, y / 1000000000.0, unit.dimensionless) for x, y in self.MIXED_TEST_RANGE]
        )


class ParseConvertTestCase(_QuantityTestCase):
    def test_convert(self):
        self.assertRange(
            [
                (unit.registry.millimeter, unit.registry.millimeter)
            ],
            [(unit.Quantity(x, unit.registry.meter), 1000 * y) for x, y in self.MIXED_TEST_RANGE]
        )


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

        self.assertEqual(str(unit.Quantity(1, unit.registry.ppb)), '1 ppb')
        self.assertEqual(str(unit.Quantity(10, unit.registry.ppb)), '10 ppb')
        self.assertEqual(str(unit.Quantity(100, unit.registry.ppb)), '100 ppb')
        self.assertEqual('1 ppm', str(unit.Quantity(1000, unit.registry.ppb)))

        self.assertEqual(str(unit.Quantity(1, unit.registry.ppm)), '1 ppm')
        self.assertEqual(str(unit.Quantity(10, unit.registry.ppm)), '10 ppm')
        self.assertEqual(str(unit.Quantity(100, unit.registry.ppm)), '100 ppm')
        self.assertEqual(str(unit.Quantity(1000, unit.registry.ppm)), '0.1%')
        self.assertEqual(str(unit.Quantity(10000, unit.registry.ppm)), '1%')

        self.assertEqual(str(unit.Quantity(1, unit.registry.percent)), '1%')
        self.assertEqual(str(unit.Quantity(10, unit.registry.percent)), '10%')
        self.assertEqual(str(unit.Quantity(100, unit.registry.percent)), '100%')
        self.assertEqual(str(unit.Quantity(1000, unit.registry.percent)), '1000%')


if __name__ == '__main__':
    unittest.main()
