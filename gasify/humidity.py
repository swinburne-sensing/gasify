import math
import warnings
from typing import Callable, Optional

from gasify.unit import Quantity, registry, dimensionless, parse, TParseQuantity


WATER_GAS_CONSTANT = Quantity(0.4615, registry.J / (registry.gram * registry.degK))

WATER_PRESSURE_CRITICAL = Quantity(22.064, registry.MPa)

WATER_TEMPERATURE_FREEZE = Quantity(0, registry.degC)
WATER_TEMPERATURE_BOIL = Quantity(100, registry.degC)
WATER_TEMPERATURE_CRITICAL = Quantity(647.096, registry.degK)

unit_absolute = registry.g / pow(registry.meter, 3)
unit_relative = registry.percent

TWaterVPCallable = Callable[[TParseQuantity], Quantity]


class TemperatureRangeWarning(UserWarning):
    pass


def _wvps_vartheta(temperature: Quantity) -> Quantity:
    return 1 - temperature / WATER_TEMPERATURE_CRITICAL


def water_vp_sat_wagner_pruss(temperature: TParseQuantity) -> Quantity:
    """ Calculate saturation vapor pressure of water at given temperature using Wagner and Pruss (1993) method.

    Reference: https://doi.org/10.1063/1.1461829

    :param temperature: gas temperature
    :return: saturation vapor pressure Quantity
    """
    temperature = parse(temperature, registry.degC)

    # noinspection PyTypeChecker
    return WATER_PRESSURE_CRITICAL * math.exp(WATER_TEMPERATURE_CRITICAL / temperature * (
        -7.85951783 * _wvps_vartheta(temperature) +
        1.84408259 * math.pow(_wvps_vartheta(temperature), 1.5) +
        -11.78649 * math.pow(_wvps_vartheta(temperature), 3) +
        22.6807411 * math.pow(_wvps_vartheta(temperature), 3.5) +
        -15.9618719 * math.pow(_wvps_vartheta(temperature), 4) +
        1.80122502 * math.pow(_wvps_vartheta(temperature), 7.5)
    ))


_WATER_VP_SAT_SIMPLE_MIN = Quantity(0, registry.degC)


def water_vp_sat_simple(temperature: TParseQuantity) -> Quantity:
    """ Calculate saturation vapor pressure of water at given temperature using the simple method.

    Reference: https://www.omnicalculator.com/chemistry/vapour-pressure-of-water

    :param temperature: gas temperature
    :return: saturation vapor pressure Quantity
    """
    if temperature < _WATER_VP_SAT_SIMPLE_MIN:
        warnings.warn(f"Simple method not suitable for calculations below {_WATER_VP_SAT_SIMPLE_MIN!s}",
                      TemperatureRangeWarning)
    elif temperature > WATER_TEMPERATURE_CRITICAL:
        warnings.warn('Simple method not suitable for calculations above critical temperature',
                      TemperatureRangeWarning)

    temperature = parse(temperature, registry.degC)

    return Quantity(math.exp(20.386 - (5132 / temperature.m_as(registry.degK))), registry.mmHg)


def water_vp_sat_antoine(temperature: TParseQuantity) -> Quantity:
    """ Calculate saturation vapor pressure of water at given temperature using the Antoine method.

    Reference: https://www.omnicalculator.com/chemistry/vapour-pressure-of-water

    :param temperature: gas temperature
    :return: saturation vapor pressure Quantity
    """
    temperature = parse(temperature, registry.degC)

    if temperature < WATER_TEMPERATURE_FREEZE:
        warnings.warn(f"Antoine method not suitable for calculations below {WATER_TEMPERATURE_FREEZE!s}",
                      TemperatureRangeWarning)
    elif temperature > WATER_TEMPERATURE_CRITICAL:
        warnings.warn('Antoine method not suitable for calculations above critical temperature',
                      TemperatureRangeWarning)

    if temperature > WATER_TEMPERATURE_BOIL:
        a = 8.14019
        b = 1810.94
        c = 244.485
    else:
        a = 8.07131
        b = 1730.63
        c = 233.426

    return Quantity(math.pow(10.0, a - (b / (c + temperature.m_as(registry.degC)))), registry.mmHg)


_WATER_VP_SAT_SIMPLE_MAX = Quantity(100, registry.degC)


def water_vp_sat_magnus(temperature: TParseQuantity) -> Quantity:
    """ Calculate saturation vapor pressure of water at given temperature using the
    Magnus/August-Roche-Magnus/Magnus-Tetens method.

    Reference: https://www.omnicalculator.com/chemistry/vapour-pressure-of-water

    :param temperature: gas temperature
    :return: saturation vapor pressure Quantity
    """
    temperature = parse(temperature, registry.degC)

    if temperature < WATER_TEMPERATURE_FREEZE:
        warnings.warn(f"Magnus method not suitable for calculations below {WATER_TEMPERATURE_FREEZE!s}",
                      TemperatureRangeWarning)
    elif temperature > _WATER_VP_SAT_SIMPLE_MAX:
        warnings.warn(f"Magnus method not suitable for calculations above {_WATER_VP_SAT_SIMPLE_MAX!s}",
                      TemperatureRangeWarning)

    return Quantity(
        0.61094 * math.exp(17.625 * temperature.m_as(registry.degC) / (temperature.m_as(registry.degC) + 243.04)),
        registry.kPa
    )


_WATER_VP_SAT_TENTENS_MAX = Quantity(75, registry.degC)


def water_vp_sat_tetens(temperature: TParseQuantity) -> Quantity:
    """

    Reference: https://www.omnicalculator.com/chemistry/vapour-pressure-of-water

    :param temperature: gas temperature
    :return: saturation vapor pressure Quantity
    """
    temperature = parse(temperature, registry.degC)

    if temperature < WATER_TEMPERATURE_FREEZE:
        warnings.warn(f"Tetens method not suitable for calculations below {WATER_TEMPERATURE_FREEZE!s}",
                      TemperatureRangeWarning)
    elif temperature > _WATER_VP_SAT_TENTENS_MAX:
        warnings.warn(f"Tetens method not suitable for calculations above {_WATER_VP_SAT_TENTENS_MAX!s}",
                      TemperatureRangeWarning)

    return Quantity(
        0.61078 * math.exp((17.27 * temperature.m_as(registry.degC)) / (temperature.m_as(registry.degC) + 237.3)),
        registry.kPa
    )


def water_vp_sat_buck(temperature: TParseQuantity) -> Quantity:
    """

    Reference: https://www.omnicalculator.com/chemistry/vapour-pressure-of-water

    :param temperature: gas temperature
    :return: saturation vapor pressure Quantity
    """
    temperature = parse(temperature, registry.degC)

    if temperature < WATER_TEMPERATURE_FREEZE:
        warnings.warn(f"Buck method not suitable for calculations below {WATER_TEMPERATURE_FREEZE!s}",
                      TemperatureRangeWarning)
    elif temperature > _WATER_VP_SAT_TENTENS_MAX:
        warnings.warn(f"Buck method not suitable for calculations above {_WATER_VP_SAT_TENTENS_MAX!s}",
                      TemperatureRangeWarning)

    temperature_m = temperature.m_as(registry.degC)

    return Quantity(
        0.61121 * math.exp((18.678 - (temperature_m / 234.5)) * (temperature_m / (257.14 + temperature_m))),
        registry.kPa
    )


def absolute_to_relative(absolute_humidity: TParseQuantity, temperature: TParseQuantity,
                         water_vp_method: Optional[TWaterVPCallable] = None) -> Quantity:
    """ Convert absolute water vapour concentration (g/m^3) to relative humidity (%) at a given temperature.

    :param absolute_humidity: absolute humidity concentration quantity
    :param temperature: gas temperature
    :param water_vp_method: method for water vapour saturation pressure calculation, defaults to Wagner-Pruss
    :return: relative humidity quantity
    """
    absolute_humidity = parse(absolute_humidity, unit_absolute)
    temperature = parse(temperature, registry.degC)
    water_vp_method = water_vp_method or water_vp_sat_wagner_pruss

    return Quantity(
        WATER_GAS_CONSTANT.magnitude * temperature.m_as(registry.degK) * absolute_humidity.m_as(unit_absolute) /
        water_vp_method(temperature).m_as(registry.Pa),
        unit_relative
    )


def relative_to_absolute(relative_humidity: TParseQuantity, temperature: TParseQuantity,
                         water_vp_method: Optional[TWaterVPCallable] = None) -> Quantity:
    """ Convert relative humidity (%) to an absolute water vapour concentration (g/m^3) at a given temperature.

    :param relative_humidity: relative humidity quantity
    :param temperature: gas temperature
    :param water_vp_method: method for water vapour saturation pressure calculation, defaults to Wagner-Pruss
    :return: absolute humidity concentration quantity
    """
    relative_humidity = parse(relative_humidity, dimensionless)
    temperature = parse(temperature, registry.degC)
    water_vp_method = water_vp_method or water_vp_sat_wagner_pruss

    return Quantity(
        relative_humidity.m_as(unit_relative) * water_vp_method(temperature).m_as(registry.Pa) / (
                    WATER_GAS_CONSTANT.magnitude * temperature.m_as(registry.degK)),
        unit_absolute
    )
