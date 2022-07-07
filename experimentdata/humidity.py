import math
import warnings
from typing import Callable, Optional, Tuple, Union

from experimentdata.unit import Quantity, registry, dimensionless, parse, TParseQuantity


WATER_GAS_CONSTANT = Quantity(0.4615, registry.J / (registry.gram * registry.degK))

WATER_PRESSURE_CRITICAL = Quantity(22.064, registry.MPa)

WATER_TEMPERATURE_FREEZE = Quantity(0, registry.degC)
WATER_TEMPERATURE_BOIL = Quantity(100, registry.degC)
WATER_TEMPERATURE_CRITICAL = Quantity(647.096, registry.degK)

unit_absolute = registry.g / pow(registry.meter, 3)
humidity_unit_relative = dimensionless

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
    absolute_humidity = parse(absolute_humidity, unit_absolute)
    temperature = parse(temperature, registry.degC)
    water_vp_method = water_vp_method or water_vp_sat_wagner_pruss

    return Quantity(
        WATER_GAS_CONSTANT.magnitude * temperature.m_as(registry.degK) * absolute_humidity.m_as(unit_absolute) /
        water_vp_method(temperature).m_as(registry.Pa),
        humidity_unit_relative
    )


def relative_to_absolute(relative_humidity: TParseQuantity, temperature: TParseQuantity,
                         water_vp_method: Optional[TWaterVPCallable] = None) -> Quantity:
    relative_humidity = parse(relative_humidity, dimensionless)
    temperature = parse(temperature, registry.degC)
    water_vp_method = water_vp_method or water_vp_sat_wagner_pruss

    return Quantity(
        relative_humidity.m_as(humidity_unit_relative) * water_vp_method(temperature).m_as(registry.Pa) / (
                    WATER_GAS_CONSTANT.magnitude * temperature.m_as(registry.degK)),
        unit_absolute
    )


# # Humidity calculation constants
# _PWS_C1 = -7.85951783
# _PWS_C2 = 1.84408259
# _PWS_C3 = -11.7866497
# _PWS_C4 = 22.6807411
# _PWS_C5 = -15.9618719
# _PWS_C6 = 1.80122502
# _HUMID_ABS_C = 2.16679
#
# _RH_SATURATED = Quantity(100, unit_humidity_rel)
#
#
# class HumidityCalculationError(ValueError):
#     pass
#
#
# def _humidity_calc_pws_exp_constants(temperature: Quantity) -> Tuple[float, float, float]:
#     temperature_mag = temperature.m_as('degC')
#
#     if -70 <= temperature_mag <= 0:
#         a = 6.114742
#         m = 9.778707
#         tn = 273.1466
#     elif 0 < temperature_mag <= 50:
#         a = 6.116441
#         m = 7.591386
#         tn = 240.7263
#     elif 50 < temperature_mag <= 100:
#         a = 6.004918
#         m = 7.337936
#         tn = 229.3975
#     elif 100 < temperature_mag <= 150:
#         a = 5.856548
#         m = 7.27731
#         tn = 225.1033
#     else:
#         raise HumidityCalculationError(f"Temperature {temperature!s} outside supported range")
#
#     return a, m, tn
#
#
# def _humidity_calc_pws_exp(temperature: Quantity) -> Quantity:
#     temperature = temperature
#     temperature_mag = temperature.m_as('degC')
#
#     (a, m, tn) = _humidity_calc_pws_exp_constants(temperature)
#
#     return Quantity(a * pow(10, (m * temperature_mag) / (temperature_mag + tn)), registry.hPa)
#
#
# def _abs_to_dew(temperature: _TYPE_INPUT, absolute_humidity: _TYPE_INPUT) -> Quantity:
#     """ Convert absolute humidity concentration (g/m^3) to a dew point temperature at a specific gas temperature.
#
#     :param temperature: temperature of gas
#     :param absolute_humidity: water concentration in gas
#     :return: dew point temperature as Quantity
#     """
#     return rel_to_dew(temperature, abs_to_rel(temperature, absolute_humidity))
#
#
# def _abs_to_rel(temperature: _TYPE_INPUT, absolute_humidity: _TYPE_INPUT) -> Quantity:
#     """ Convert absolute humidity concentration (g/m^3) to a relative humidity (%) at a specific gas temperature.
#
#     :param temperature: temperature of the gas
#     :param absolute_humidity: water concentration in gas
#     :return: relative humidity percentage as Quantity
#     """
#     temperature = parse(temperature, registry.degC).to(registry.degK)
#     absolute_humidity = parse(absolute_humidity, unit_humidity_abs)
#
#     pw = temperature * absolute_humidity / _HUMID_ABS_C
#
#     return Quantity((pw / _humidity_calc_pws_exp(temperature)).magnitude, unit_humidity_rel)
#
#
# # Calculate absolute humidity from dew point
# def _dew_to_abs(dew_temperature: _TYPE_INPUT) -> Quantity:
#     """ Convert dew point temperature to absolute humidity.
#
#     :param dew_temperature: dew point temperature
#     :return: water concentration in gas as Quantity
#     """
#     return rel_to_abs(dew_temperature, _RH_SATURATED)
#
#
# # Calculate relative humidity from temperature and dew point
# def _dew_to_rel(temperature: _TYPE_INPUT, dew_temperature: _TYPE_INPUT) -> Quantity:
#     """ Convert dew point temperature to relative humidity.
#
#     :param temperature: temperature of gas
#     :param dew_temperature: dew point temperature
#     :return:
#     """
#     temperature = parse(temperature, registry.degC)
#     dew_temperature = parse(dew_temperature, registry.degC)
#
#     pws = _humidity_calc_pws_exp(temperature)
#     pwd = _humidity_calc_pws_exp(dew_temperature)
#
#     return Quantity(pwd / pws, dimensionless).to(unit_humidity_rel)
#
#
# # Calculate dew point from temperature and relative_humidity
# def _rel_to_dew(temperature: _TYPE_INPUT, relative_humidity: _TYPE_INPUT) -> Quantity:
#     """ Convert relative humidity to dew point temperature.
#
#     :param temperature: temperature of gas
#     :param relative_humidity:
#     :return: dew point temperature as Quantity
#     """
#     temperature = parse(temperature, registry.degC)
#     relative_humidity = parse(relative_humidity, dimensionless)
#
#     (a, m, tn) = _humidity_calc_pws_exp_constants(temperature)
#     pws = _humidity_calc_pws_exp(temperature) * relative_humidity.m_as(dimensionless)
#
#     return Quantity(tn / (m / (math.log10(pws.m_as('hPa') / a)) - 1), registry.degC)
#
#
# # Calculate absolute humidity from temperature and relative humidity measurement
# def _rel_to_abs(temperature: _TYPE_INPUT, relative_humidity: _TYPE_INPUT) -> Quantity:
#     """ Convert a relative humidity (%) at a specific temperature to absolute humidity concentration (g/m^3).
#
#     :param temperature: temperature of gas
#     :param relative_humidity:
#     :return: water concentration in gas as Quantity
#     """
#     # Convert types
#     temperature = parse(temperature, registry.degC)
#     relative_humidity = parse(relative_humidity, dimensionless)
#
#     pw = _humidity_calc_pws_exp(temperature) * relative_humidity.m_as(dimensionless)
#
#     return Quantity((_HUMID_ABS_C * pw.to('Pa') / temperature).magnitude, unit_humidity_abs)
