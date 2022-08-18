import typing
import unittest

from gasify import gas, unit


# Values from https://www.mks.com/n/gas-correction-factors-for-thermal-based-mass-flow-controllers
_GCF_LUT = {
    'ammonia': 0.73,
    'argon': 1.39,
    'arsine': 0.67,
    'bromine': 0.81,
    'chlorine': 0.86,
    'cyanogen': 0.61,
    'deuterium': 1,
    'ethane': 0.5,
    'fluorine': 0.98,
    'helium': 1.45,
    'hydrogen': 1.01,
    'methane': 0.72,
    'neon': 1.46,
    'nitrogen': 1,
    'oxygen': 0.993,
    'phosphine': 0.76,
    'propane': 0.36,
    'propylene': 0.41,
    'xenon': 1.32
}


class CompoundTestCase(unittest.TestCase):
    def test_equal(self):
        self.assertEqual(gas.registry.oxygen, gas.registry.oxygen)
        self.assertEqual(gas.registry.oxygen, gas.Compound('Oxygen'))
        self.assertNotEqual(gas.registry.oxygen, gas.registry.nitrogen)
        self.assertNotEqual(gas.registry.oxygen, gas.Compound('Nitrogen'))

    def test_compare(self):
        a = gas.registry.oxygen
        b = gas.registry.hydrogen

        self.assertLess(b, a)
        self.assertLessEqual(b, a)
        self.assertGreater(a, b)
        self.assertGreaterEqual(a, b)

    # def test_gcf(self):
    #     for gas_name, gcf in _GCF_LUT.items():
    #         with self.subTest(gas_name):
    #             self.assertAlmostEqual(gas.registry[gas_name].gcf, gcf, 2)

    def test_sort(self):
        self.assertListEqual(
            sorted([
                gas.registry.oxygen,
                gas.registry.hydrogen,
                gas.registry.nitrogen
            ]),
            [
                gas.registry.hydrogen,
                gas.registry.oxygen,
                gas.registry.nitrogen
            ]
        )


class CompoundConcentrationTestCase(unittest.TestCase):
    def assertConcentration(self, x: typing.Any, concentration: float, compound: gas.Compound):
        self.assertIsInstance(x, gas.CompoundConcentration, f"Expected instance of CompoundConcentration, "
                                                            f"got {type(x)}")
        self.assertAlmostEqual(concentration, x.concentration.m_as(unit.dimensionless), 9, 'Concentration is incorrect')
        self.assertIs(compound, x.compound, 'Compound is incorrect')

    def test_init(self):
        self.assertConcentration(
            gas.CompoundConcentration(0.1, gas.registry.oxygen),
            0.1,
            gas.registry.oxygen
        )

        self.assertConcentration(
            0.1 * gas.registry.oxygen,
            0.1,
            gas.registry.oxygen
        )

        self.assertConcentration(
            gas.registry.oxygen * 0.1,
            0.1,
            gas.registry.oxygen
        )

        self.assertConcentration(
            gas.registry.oxygen * unit.Quantity(0.1, unit.dimensionless),
            0.1,
            gas.registry.oxygen
        )

    def test_init_error(self):
        with self.assertRaises(TypeError):
            _ = gas.CompoundConcentration(0.1, gas.registry.oxygen) * object()

    def test_equal(self):
        a = 1.0 * gas.registry.oxygen
        b = 1.0 * gas.registry.oxygen
        c = 1.0 * gas.registry.nitrogen

        self.assertNotEqual(a, object())
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)
        self.assertEqual(a, gas.registry.oxygen)
        self.assertEqual(a, 1.0)

    def test_equal_compound(self):
        a = 1.0 * gas.registry.oxygen

        self.assertEqual(a, gas.registry.oxygen)
        self.assertEqual(gas.registry.oxygen, a)

    def test_compare(self):
        a = 1.0 * gas.registry.oxygen
        b = 0.5 * gas.registry.oxygen

        self.assertLess(b, a)
        self.assertLessEqual(b, a)
        self.assertGreater(a, b)
        self.assertGreaterEqual(a, b)

    def test_compare_numeric(self):
        a = 1.0 * gas.registry.oxygen

        self.assertLess(a, 2.0)
        self.assertLessEqual(a, 1.0)
        self.assertGreater(a, 0.5)
        self.assertGreaterEqual(a, 1.0)

    def test_compare_error(self):
        with self.assertRaises(TypeError):
            _ = gas.CompoundConcentration(0.1, gas.registry.oxygen) < object()

        with self.assertRaises(TypeError):
            _ = gas.CompoundConcentration(0.1, gas.registry.oxygen) <= object()

        with self.assertRaises(TypeError):
            _ = gas.CompoundConcentration(0.1, gas.registry.oxygen) > object()

        with self.assertRaises(TypeError):
            _ = gas.CompoundConcentration(0.1, gas.registry.oxygen) >= object()

    def test_compare_different(self):
        a = 1.0 * gas.registry.oxygen
        b = 0.5 * gas.registry.hydrogen

        self.assertGreater(a, b)
        self.assertGreaterEqual(a, b)
        self.assertLess(b, a)
        self.assertLessEqual(b, a)

    def test_div(self):
        a = 1.0 * gas.registry.oxygen

        with self.subTest('as int'):
            self.assertConcentration(a / 2, 0.5, gas.registry.oxygen)

        with self.subTest('as float'):
            self.assertConcentration(a / 2.0, 0.5, gas.registry.oxygen)

        with self.subTest('as Quantity'):
            b = unit.Quantity(2, unit.dimensionless)
            self.assertConcentration(a / b, 0.5, gas.registry.oxygen)

    def test_equal(self):
        a = 1.0 * gas.registry.oxygen
        b = 1.0 * gas.registry.oxygen

        self.assertEqual(a, b)

    def test_mul(self):
        a = 0.5 * gas.registry.oxygen

        with self.subTest('as int'):
            self.assertConcentration(2 * a, 1.0, gas.registry.oxygen)
            self.assertConcentration(a * 2, 1.0, gas.registry.oxygen)

        with self.subTest('as float'):
            self.assertConcentration(2.0 * a, 1.0, gas.registry.oxygen)
            self.assertConcentration(a * 2.0, 1.0, gas.registry.oxygen)

        with self.subTest('as Quantity'):
            b = unit.Quantity(2, unit.dimensionless)
            self.assertConcentration(a * b, 1.0, gas.registry.oxygen)


class MixtureTestCase(unittest.TestCase):
    def assertMixture(self, x: typing.Any, components: typing.Tuple[gas.CompoundConcentration, ...],
                      component_total: float, is_analyte: bool) -> None:
        self.assertIsInstance(x, gas.Mixture, f"Expected instance of Mixture, got {type(x)}")
        self.assertEqual(is_analyte, x.analyte, f"Expected analyte flag to be {is_analyte!r}, got {x.analyte!r}")
        self.assertAlmostEqual(component_total, x.total.m_as(unit.dimensionless), 9, 'Mixture total is incorrect')

        self.assertEqual(len(components), len(x.content), 'Number of components does not match expected count')

        for index, (expected, got) in enumerate(zip(components, x.content)):
            self.assertIsInstance(got, gas.CompoundConcentration, f"Expected instance of CompoundConcentration at "
                                                                  f"index {index}, got {type(got)}")
            self.assertEqual(expected, got, f"Component mismatch at index {index}, expected {expected!s}, got {got!s}")

    def test_gcf(self):
        self.assertAlmostEqual(1.0, gas.air.gcf, 2)

    def test_gcf_error(self):
        a = gas.Mixture(0.1 * gas.registry.nitric_oxides, 0.2 * gas.registry.hydrogen)

        with self.assertRaises(ValueError):
            _ = a.gcf

    def test_contains(self):
        a = gas.Mixture(0.1 * gas.registry.oxygen, 0.2 * gas.registry.hydrogen)

        self.assertIn(gas.registry.oxygen, a)
        self.assertIn(gas.registry.hydrogen, a)
        self.assertNotIn(gas.registry.nitrogen, a)

    def test_mix_same(self):
        a = 0.1 * gas.registry.oxygen
        b = 0.2 * gas.registry.oxygen

        self.assertMixture(a + b, (0.3 * gas.registry.oxygen,), 0.3, False)
        self.assertMixture(b + a, (0.3 * gas.registry.oxygen,), 0.3, False)

    def test_mix_different(self):
        a = 0.1 * gas.registry.oxygen
        b = 0.2 * gas.registry.hydrogen

        self.assertMixture(
            a + b,
            (
                0.2 * gas.registry.hydrogen,
                0.1 * gas.registry.oxygen
            ),
            0.3,
            True
        )

        self.assertMixture(
            b + a,
            (
                0.2 * gas.registry.hydrogen,
                0.1 * gas.registry.oxygen
            ),
            0.3,
            True
        )

    def test_mix_auto(self):
        a = 0.1 * gas.registry.oxygen
        b = 0.2 * gas.registry.hydrogen

        self.assertMixture(
            gas.Mixture.auto_balance(a, b, balance=gas.registry.nitrogen),
            (
                0.2 * gas.registry.hydrogen,
                0.1 * gas.registry.oxygen,
                0.7 * gas.registry.nitrogen
            ),
            1.0,
            True
        )

    def test_mix_add(self):
        a = gas.Mixture(0.1 * gas.registry.oxygen, 0.2 * gas.registry.hydrogen)
        b = 0.3 * gas.registry.nitrogen

        self.assertMixture(
            a + b,
            (
                0.2 * gas.registry.hydrogen,
                0.1 * gas.registry.oxygen,
                0.3 * gas.registry.nitrogen
            ),
            0.6,
            True
        )

        self.assertMixture(
            b + a,
            (
                0.2 * gas.registry.hydrogen,
                0.1 * gas.registry.oxygen,
                0.3 * gas.registry.nitrogen
            ),
            0.6,
            True
        )

    def test_mix_self(self):
        a = gas.Mixture(0.1 * gas.registry.oxygen, 0.2 * gas.registry.hydrogen)
        b = gas.Mixture(0.2 * gas.registry.oxygen, 0.3 * gas.registry.hydrogen)

        self.assertMixture(
            a + b,
            (
                0.5 * gas.registry.hydrogen,
                0.3 * gas.registry.oxygen
            ),
            0.8,
            True
        )

        self.assertMixture(
            b + a,
            (
                0.5 * gas.registry.hydrogen,
                0.3 * gas.registry.oxygen
            ),
            0.8,
            True
        )

    def test_mix_add_mix(self):
        a = gas.Mixture(0.1 * gas.registry.oxygen, 0.2 * gas.registry.hydrogen)
        b = gas.Mixture(0.3 * gas.registry.nitrogen, 0.4 * gas.registry.argon)

        self.assertMixture(
            a + b,
            (
                0.2 * gas.registry.hydrogen,
                0.1 * gas.registry.oxygen,
                0.4 * gas.registry.argon,
                0.3 * gas.registry.nitrogen
            ),
            1.0,
            True
        )

        self.assertMixture(
            b + a,
            (
                0.2 * gas.registry.hydrogen,
                0.1 * gas.registry.oxygen,
                0.4 * gas.registry.argon,
                0.3 * gas.registry.nitrogen
            ),
            1.0,
            True
        )

    def test_mix_add_list(self):
        a = gas.Mixture(0.1 * gas.registry.oxygen, 0.2 * gas.registry.hydrogen)
        b = [0.3 * gas.registry.nitrogen, 0.4 * gas.registry.argon]

        self.assertMixture(
            a + b,
            (
                0.2 * gas.registry.hydrogen,
                0.1 * gas.registry.oxygen,
                0.4 * gas.registry.argon,
                0.3 * gas.registry.nitrogen
            ),
            1.0,
            True
        )

        self.assertMixture(
            b + a,
            (
                0.2 * gas.registry.hydrogen,
                0.1 * gas.registry.oxygen,
                0.4 * gas.registry.argon,
                0.3 * gas.registry.nitrogen
            ),
            1.0,
            True
        )

    def test_norm(self):
        a = gas.Mixture(0.1 * gas.registry.oxygen, 0.1 * gas.registry.hydrogen)

        self.assertMixture(
            a.normalise(),
            (
                0.5 * gas.registry.hydrogen,
                0.5 * gas.registry.oxygen
            ),
            1.0,
            True
        )

        b = gas.Mixture(3.0 * gas.registry.oxygen, 1.0 * gas.registry.hydrogen)

        self.assertMixture(
            b.normalise(),
            (
                0.25 * gas.registry.hydrogen,
                0.75 * gas.registry.oxygen
            ),
            1.0,
            True
        )

        c = gas.Mixture(0.5 * gas.registry.oxygen, 1.0 * gas.registry.hydrogen, 2.5 * gas.registry.nitrogen)

        self.assertMixture(
            c.normalise(),
            (
                0.25 * gas.registry.hydrogen,
                0.125 * gas.registry.oxygen,
                0.625 * gas.registry.nitrogen
            ),
            1.0,
            True
        )


if __name__ == '__main__':
    unittest.main()
