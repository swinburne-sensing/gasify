import unittest

from gasify import humidity, unit

from tests.util import QuantityTestCase


class WaterVaporPressureTestCase(QuantityTestCase):
    _REFERENCE_VALUES = [
        # Reference values from doi: 10.6028/jres.073A.039
        (unit.Quantity(25, unit.registry.degC), unit.Quantity(3168.6, unit.registry.Pa)),
        (unit.Quantity(40, unit.registry.degC), unit.Quantity(7381.3, unit.registry.Pa)),
        (unit.Quantity(50, unit.registry.degC), unit.Quantity(12344.6, unit.registry.Pa)),
        (unit.Quantity(60, unit.registry.degC), unit.Quantity(19344.6, unit.registry.Pa)),
        (unit.Quantity(70, unit.registry.degC), unit.Quantity(31177.0, unit.registry.Pa)),
        (unit.Quantity(80, unit.registry.degC), unit.Quantity(47375.2, unit.registry.Pa)),
        (unit.Quantity(100, unit.registry.degC), unit.Quantity(101325.0, unit.registry.Pa)),
    ]

    def test_wagner_pruss(self):
        test_values = [
            (unit.Quantity(-100, unit.registry.degC), unit.Quantity(0.003683, unit.registry.Pa), 6),
            (unit.Quantity(-75, unit.registry.degC), unit.Quantity(0.25484, unit.registry.Pa), 6),
            (unit.Quantity(-50, unit.registry.degC), unit.Quantity(6.447, unit.registry.Pa), 6),
            (unit.Quantity(-25, unit.registry.degC), unit.Quantity(80.88, unit.registry.Pa), 6),
            (humidity.WATER_TEMPERATURE_FREEZE, unit.Quantity(611.2, unit.registry.Pa), 6),
            (unit.Quantity(25, unit.registry.degC), unit.Quantity(3170, unit.registry.Pa), 6),
            (unit.Quantity(50, unit.registry.degC), unit.Quantity(12352, unit.registry.Pa), 5),
            (unit.Quantity(75, unit.registry.degC), unit.Quantity(38597, unit.registry.Pa), 5),
            (unit.Quantity(100, unit.registry.degC), unit.Quantity(101418, unit.registry.Pa), 5),
            (unit.Quantity(150, unit.registry.degC), unit.Quantity(476159, unit.registry.Pa), 5),
            (unit.Quantity(200, unit.registry.degC), unit.Quantity(1554939, unit.registry.Pa), 4),
            (unit.Quantity(300, unit.registry.degC), unit.Quantity(8587867, unit.registry.Pa), 5),
            (humidity.WATER_TEMPERATURE_CRITICAL, unit.Quantity(21813821, unit.registry.Pa), 0),
        ]

        for temperature, expected, places in test_values:
            with self.subTest(f"p({temperature}) = {expected}"):
                self.assertQuantity(
                    humidity.water_vp_sat_wagner_pruss(temperature),
                    expected,
                    places,
                    unit.registry.MPa
                )

    def test_antoine(self):
        test_values = [
            (humidity.WATER_TEMPERATURE_FREEZE, unit.Quantity(0.0006056, unit.registry.MPa), 6),
            (unit.Quantity(25, unit.registry.degC), unit.Quantity(0.003158, unit.registry.MPa), 5),
            (unit.Quantity(50, unit.registry.degC), unit.Quantity(0.012306, unit.registry.MPa), 5),
            (unit.Quantity(75, unit.registry.degC), unit.Quantity(0.03846, unit.registry.MPa), 4),
            (unit.Quantity(100, unit.registry.degC), unit.Quantity(0.10134, unit.registry.MPa), 4),
            (unit.Quantity(150, unit.registry.degC), unit.Quantity(0.47255, unit.registry.MPa), 4),
            (unit.Quantity(200, unit.registry.degC), unit.Quantity(1.552, unit.registry.MPa), 3),
            (unit.Quantity(300, unit.registry.degC), unit.Quantity(8.692, unit.registry.MPa), 3),
            (humidity.WATER_TEMPERATURE_CRITICAL, unit.Quantity(21.73, unit.registry.MPa), 1),
        ]

        for temperature, expected, places in test_values:
            with self.subTest(f"p({temperature}) = {expected}"):
                self.assertQuantity(
                    humidity.water_vp_sat_antoine(temperature),
                    expected,
                    places,
                    unit.registry.MPa
                )

    def test_antoine_warning(self):
        with self.assertWarns(UserWarning):
            humidity.water_vp_sat_antoine(-1)

    def test_simple(self):
        test_values = [
            (humidity.WATER_TEMPERATURE_FREEZE, unit.Quantity(0.0006521, unit.registry.MPa), 4),
            (unit.Quantity(25, unit.registry.degC), unit.Quantity(0.003157, unit.registry.MPa), 4),
            (unit.Quantity(50, unit.registry.degC), unit.Quantity(0.01197, unit.registry.MPa), 3),
            (unit.Quantity(75, unit.registry.degC), unit.Quantity(0.03748, unit.registry.MPa), 3),
            (unit.Quantity(100, unit.registry.degC), unit.Quantity(0.10072, unit.registry.MPa), 2),
            (unit.Quantity(150, unit.registry.degC), unit.Quantity(0.47255, unit.registry.MPa), 1)
        ]

        for temperature, expected, places in test_values:
            with self.subTest(f"p({temperature}) = {expected}"):
                self.assertQuantity(
                    humidity.water_vp_sat_simple(temperature),
                    expected,
                    places,
                    unit.registry.MPa
                )

    def test_magnus(self):
        test_values = [
            (humidity.WATER_TEMPERATURE_FREEZE, unit.Quantity(0.000611, unit.registry.MPa), 4),
            (unit.Quantity(25, unit.registry.degC), unit.Quantity(0.003162, unit.registry.MPa), 5),
            (unit.Quantity(50, unit.registry.degC), unit.Quantity(0.01236, unit.registry.MPa), 3),
            (unit.Quantity(75, unit.registry.degC), unit.Quantity(0.039, unit.registry.MPa), 2),
            (unit.Quantity(100, unit.registry.degC), unit.Quantity(0.10408, unit.registry.MPa), 2),
            (unit.Quantity(150, unit.registry.degC), unit.Quantity(0.5096, unit.registry.MPa), 4),
            (unit.Quantity(200, unit.registry.degC), unit.Quantity(1.7435, unit.registry.MPa), 3),
            (unit.Quantity(300, unit.registry.degC), unit.Quantity(10.343, unit.registry.MPa), 3)
        ]

        for temperature, expected, places in test_values:
            with self.subTest(f"p({temperature}) = {expected}"):
                self.assertQuantity(
                    humidity.water_vp_sat_magnus(temperature).to(unit.registry.kPa),
                    expected,
                    places,
                    unit.registry.MPa
                )

    def test_tetens(self):
        test_values = [
            (humidity.WATER_TEMPERATURE_FREEZE, unit.Quantity(0.0006108, unit.registry.MPa), 4),
            (unit.Quantity(25, unit.registry.degC), unit.Quantity(0.0031677, unit.registry.MPa), 5),
            (unit.Quantity(50, unit.registry.degC), unit.Quantity(0.012336, unit.registry.MPa), 3),
            (unit.Quantity(75, unit.registry.degC), unit.Quantity(0.038646, unit.registry.MPa), 2),
            (unit.Quantity(100, unit.registry.degC), unit.Quantity(0.10221, unit.registry.MPa), 2),
            (unit.Quantity(150, unit.registry.degC), unit.Quantity(0.4906, unit.registry.MPa), 4),
            (unit.Quantity(200, unit.registry.degC), unit.Quantity(1.645, unit.registry.MPa), 3),
            (unit.Quantity(300, unit.registry.degC), unit.Quantity(9.411, unit.registry.MPa), 3)
        ]

        for temperature, expected, places in test_values:
            with self.subTest(f"p({temperature}) = {expected}"):
                self.assertQuantity(
                    humidity.water_vp_sat_tetens(temperature).to(unit.registry.kPa),
                    expected,
                    places,
                    unit.registry.MPa
                )

    def test_buck(self):
        test_values = [
            (humidity.WATER_TEMPERATURE_FREEZE, unit.Quantity(0.0006112, unit.registry.MPa), 4),
            (unit.Quantity(25, unit.registry.degC), unit.Quantity(0.0031685, unit.registry.MPa), 5),
            (unit.Quantity(50, unit.registry.degC), unit.Quantity(0.01235, unit.registry.MPa), 3),
            (unit.Quantity(75, unit.registry.degC), unit.Quantity(0.038595, unit.registry.MPa), 2),
            (unit.Quantity(100, unit.registry.degC), unit.Quantity(0.1013, unit.registry.MPa), 2),
            (unit.Quantity(150, unit.registry.degC), unit.Quantity(0.4703, unit.registry.MPa), 4),
            (unit.Quantity(200, unit.registry.degC), unit.Quantity(1.4895, unit.registry.MPa), 3),
            (unit.Quantity(300, unit.registry.degC), unit.Quantity(7.16, unit.registry.MPa), 3)
        ]

        for temperature, expected, places in test_values:
            with self.subTest(f"p({temperature}) = {expected}"):
                self.assertQuantity(
                    humidity.water_vp_sat_buck(temperature).to(unit.registry.kPa),
                    expected,
                    places,
                    unit.registry.MPa
                )


class TestDataHumidity(QuantityTestCase):
    def test_rel_to_abs(self):
        rel_humid = unit.parse('100%')

        abs_humid = unit.parse('30.359 g/m^3')
        calc_abs_humid = humidity.relative_to_absolute(rel_humid, unit.Quantity(30, unit.registry.degC))

        self.assertQuantity(
            abs_humid,
            calc_abs_humid,
            2,
            humidity.unit_absolute
        )

    def test_abs_to_rel(self):
        abs_humid = unit.parse('30.359 g/m^3')

        rel_humid = unit.parse('100%')
        calc_rel_humid = humidity.absolute_to_relative(abs_humid, unit.Quantity(30, unit.registry.degC))

        self.assertQuantity(
            rel_humid,
            calc_rel_humid,
            2,
            unit.dimensionless
        )


if __name__ == '__main__':
    unittest.main()
