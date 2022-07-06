from __future__ import annotations

import functools
import typing
from datetime import timedelta

import pint
import pint.formatting

try:
    import pint_pandas
except ImportError:
    pint_pandas = None


__all__ = [
    'ParseError',
    'registry',
    'Quantity',
    'Unit',
    'dimensionless'
]


class ParseError(Exception):
    pass


# Handler for percent sign and micro symbol
def _handle_symbols(x):
    return x.replace('%', ' percent ').replace('μ', 'u')


# Unit registry
registry = pint.UnitRegistry(autoconvert_offset_to_baseunit=True, preprocessors=[_handle_symbols])

# Define additional units
registry.define('percent = count / 100 = %')
registry.define('parts_per_million = count / 1e6 = ppm')
registry.define('parts_per_billion = count / 1e9 = ppb')

registry.define('cubic_centimeter_per_minute = cm ** 3 / min = ccm')
registry.define('standard_cubic_centimeter_per_minute = cm ** 3 / min = sccm')
# registry.define('litre_per_minute = l / min = lpm')

# Add aliases
registry.define('@alias psi = PSI')
registry.define('@alias ccm = CCM')
registry.define('@alias sccm = SCCM')
# registry.define('@alias m ** 3 = m3/d')


class Quantity(registry.Quantity):
    _DISTANCE_MAX = 1000.0

    def to_compact(self, unit=None) -> Quantity:
        if self.units == dimensionless:
            # Make copy
            return self.to(dimensionless)

        if self.is_compatible_with(registry.meter) and self.m_as(registry.kilometer) >= self._DISTANCE_MAX:
            # Clamp distances to kilometers
            return self.to(registry.kilometer)

        to_mag = self.magnitude
        to_unit = self.units

        # Scale only within "dimensionless" type units, don't append SI prefix
        if to_unit == registry.percent and to_mag < 0.1:
            to_mag *= 10000.0
            to_unit = registry.ppm

        if to_unit == registry.ppm and to_mag < 1.0:
            to_mag *= 1000.0
            to_unit = registry.ppb

        if to_unit == registry.ppb and to_mag >= 1000.0:
            to_mag /= 1000.0
            to_unit = registry.ppm

        if to_unit == registry.ppm and to_mag >= 1000.0:
            to_mag /= 10000.0
            to_unit = registry.percent

        if to_unit in [registry.percent, registry.ppm, registry.ppb]:
            return Quantity(to_mag, to_unit)

        return super().to_compact(unit)

    def __format__(self, spec: str) -> str:
        # Patch to allow use of #/to_compact in default_format
        if spec == '':
            spec = self.default_format

        formatted = super().__format__(spec)
        formatted_num, formatted_unit = formatted.split(' ', 1)

        if formatted_unit.strip() == '%':
            # Remove space from percentages
            return formatted_num + '%'

        return formatted


class Unit(registry.Unit):
    def __format__(self, spec: str) -> str:
        if self == registry.percent:
            return '%'

        return super().__format__(spec)


registry.Quantity = Quantity
registry.Unit = Unit

# Shortcuts for dimensionless quantities (must occur after subclassing of Unit)
dimensionless = registry.dimensionless


# Change default printing format
@pint.register_unit_format('edata')
def format_custom(unit, registry, **options):
    unit_str = pint.formatter(
        unit.items(),
        as_ratio=True,
        single_denominator=False,
        product_fmt=" ",
        division_fmt="/",
        power_fmt="{}**{}",
        parentheses_fmt=r"({})",
        **options,
    )

    return unit_str


registry.default_format = 'g~#edata'


# Handle pickle/unpickling by overwriting the built-in unit registry
pint.set_application_registry(registry)


# Shortcuts
Quantity = registry.Quantity
Unit = registry.Unit


# Setup pint arrays (experimental)
if pint_pandas is not None:
    pint_pandas.PintType.ureg = registry


# Type hints
T_PARSE_QUANTITY = typing.Union[Quantity, str, float, int]
T_PARSE_UNIT = typing.Union[Unit, Quantity, str]
T_PARSE_TIMEDELTA = typing.Union[timedelta, Quantity, str, float, int]


def is_unit(x: typing.Any) -> bool:
    return isinstance(x, Unit)


def is_quantity(x: typing.Any) -> bool:
    return isinstance(x, Quantity)


def parse_unit(x: T_PARSE_UNIT) -> Unit:
    """ Parse arbitrary input to a Unit from the registry.

    :param x: input str
    :return: Unit
    """
    if is_unit(x):
        # Already a Unit
        return x

    if is_quantity(x):
        # Extract Unit, can sometimes occur when using values from pint
        return x.units

    if not isinstance(x, str):
        raise ParseError(f"Unsupported input type \"{type(x)}\"")

    if hasattr(registry, x):
        return getattr(registry, x)

    raise ParseError(f"Unknown unit \"{x}\"")


def parse(x: T_PARSE_QUANTITY, to_unit: typing.Optional[T_PARSE_UNIT] = None,
          mag_round: typing.Optional[int] = None) -> Quantity:
    """ Parse arbitrary input to a Quantity of specified unit.

    :param x: input str, number or Quantity
    :param to_unit: str or Unit to convert parsed values to
    :param mag_round:
    :return: Quantity with parsed magnitude and specified unit
    """
    if x is None:
        raise ParseError('Cannot convert NoneType to Quantity')

    # Parse unit
    if to_unit is not None:
        to_unit = parse_unit(to_unit)

    if not is_quantity(x):
        # Convert int to float
        if isinstance(x, int):
            x = float(x)

        # Convert floats (and ints) to Quantity, attempt to directly parse strings
        if isinstance(x, float) or isinstance(x, str):
            x = Quantity(x)
        else:
            raise ParseError(f"Unsupported input type \"{type(x)}\"")

    # Attempt conversion
    if to_unit is not None:
        if not x.unitless:
            try:
                # Don't use in-place change, can mess up values passed to some methods
                x = x.to(to_unit)
            except pint.errors.DimensionalityError as ex:
                raise ParseError(f"Unable to convert parsed quantity {x!s} to unit {to_unit}") from ex
        else:
            x = Quantity(x.m_as(dimensionless), to_unit)

    # x = typing.cast(Quantity, x)

    if mag_round is not None:
        # Round resulting value
        x = round(x, mag_round)

    return x


def parse_magnitude(x: T_PARSE_QUANTITY, magnitude_unit: T_PARSE_UNIT,
                    input_unit: typing.Optional[T_PARSE_UNIT] = None) -> float:
    """ Shortcut method to parse as value, optionally converting to specified unit before returning the magnitude.

    :param x: input str, number or Quantity
    :param magnitude_unit: str or Unit to convert parsed values to before conversion to magnitude
    :param input_unit: str or Unit to convert parsed values to
    :return:
    """
    magnitude_unit = parse_unit(magnitude_unit)

    if input_unit is None:
        # Assume default parsing unit is same as casting unit
        input_unit = magnitude_unit

    return parse(x, input_unit).m_as(magnitude_unit)


def parse_timedelta(x: T_PARSE_TIMEDELTA) -> timedelta:
    """

    :param x:
    :return:
    """
    if isinstance(x, timedelta):
        # Already a timedelta
        return x
    elif isinstance(x, float) or isinstance(x, int):
        # Count as seconds
        return timedelta(seconds=x)

    x_unit = parse(x)

    if x_unit.dimensionless:
        # Assume seconds by default
        x_unit = Quantity(x_unit.m_as(dimensionless), registry.sec)

    x_secs = x_unit.m_as(registry.sec)

    return timedelta(seconds=x_secs)


def converter(to_unit: typing.Optional[T_PARSE_UNIT] = None,
              optional: bool = False) -> typing.Callable[[T_PARSE_QUANTITY], Quantity]:
    """ Create wrapper for parse decorator with a pre-defined unit. Useful with the attrs library.

    :param to_unit: str or Unit to convert values to, defaults to unitless
    :param optional: if False
    :return:
    """
    to_unit = to_unit or dimensionless

    def f(x: T_PARSE_QUANTITY):
        if x is None:
            if not optional:
                raise ParseError('Input to converter cannot be None')

            return None

        return parse(x, to_unit)

    return f


def return_converter(to_unit: T_PARSE_UNIT, allow_none: bool = False):
    """ Decorator to convert returned result to a Quantity.

    :param to_unit:
    :param allow_none:
    :return:
    """
    to_unit = parse_unit(to_unit)

    def wrapper_decorator(func):
        @functools.wraps(func)
        def wrapper_result(*args, **kwargs):
            result = func(*args, **kwargs)

            if result is None:
                if not allow_none:
                    raise ValueError('Expected numeric result')

                return None

            if not is_quantity(result):
                raise ValueError(f"Decorated method returned {type(result)}, expected Quantity")

            return result.to(to_unit)

        return wrapper_result

    return wrapper_decorator
