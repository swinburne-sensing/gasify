from __future__ import annotations

import functools
import typing
from datetime import timedelta

import pint
import pint.formatting
import pint.registry


__all__ = [
    'converter',
    'dimensionless',
    'parse',
    'parse_magnitude',
    'parse_timedelta',
    'parse_unit',
    'UnitError',
    'IncompatibleUnits',
    'ParseError',
    'UnknownUnit',
    'Quantity',
    'registry',
    'return_converter',
    'TParseQuantity',
    'TParseTimeDelta',
    'TParseUnit',
    'Unit'
]


# Get split_format from pint.formatting, or implement a simplified version for Python 3.7 compatibility
if hasattr(pint.formatting, 'split_format'):
    _split_format = pint.formatting.split_format
else:
    def _split_format(spec: str, default: str, separate_format_defaults: bool = True) -> typing.Tuple[str, str]:
        mspec = pint.formatting.remove_custom_flags(spec)
        uspec = pint.formatting.extract_custom_flags(spec)

        if separate_format_defaults:
            mspec = mspec if mspec else pint.formatting.remove_custom_flags(default)
            uspec = uspec if uspec else pint.formatting.extract_custom_flags(default)

        return mspec, uspec


# Type hints
TParseQuantity = typing.Union['Quantity', str, float, int]
TParseUnit = typing.Union['Unit', 'Quantity', str]
TParseTimeDelta = typing.Union[timedelta, 'Quantity', str, float, int]


class UnitError(Exception):
    pass


class IncompatibleUnits(UnitError):
    pass


class UnknownUnit(UnitError):
    pass


class ParseError(UnitError, ValueError):
    pass


# Handler for percent sign and micro symbol
def _handle_symbols(x: str) -> str:
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


class Quantity(pint.Quantity[float]):
    _REGISTRY = registry
    _DISTANCE_MAX = 1000.0

    def to_compact(self, unit: typing.Optional[Unit] = None) -> Quantity:
        if unit is not None:
            return super().to_compact(unit)

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

        formatted: str = super().__format__(spec)

        if ' ' not in formatted:
            return formatted

        formatted_num, formatted_unit = formatted.split(' ', 1)

        if formatted_unit.strip() == '%':
            # Remove space from percentages
            return formatted_num + '%'

        return formatted


class Unit(pint.Unit):
    _REGISTRY = registry

    def __format__(self, spec: str) -> str:
        if self == registry.percent:
            return '%'

        return super().__format__(spec)


class Measurement(pint.Measurement):
    _REGISTRY = registry

    def __format__(self, format_spec: str) -> str:
        format_spec = format_spec or self.default_format

        if hasattr(self._REGISTRY, 'separate_format_defaults'):
            separate_format_defaults = self._REGISTRY.separate_format_defaults
        else:
            separate_format_defaults = True

        mspec, uspec = _split_format(
            format_spec, self.default_format, separate_format_defaults
        )

        if '#' in mspec:
            mspec = mspec.replace('#', '')
            obj = self.to_compact()
        else:
            obj = self

        if mspec == 'g':
            mstr = format(obj.magnitude, '.6g').replace('+/-', '±')

            if '±' in mstr:
                mag, sep, error = mstr.partition('±')
            else:
                mag = mstr
                sep = error = ''

            mstr = mag.rstrip('0').rstrip('.') + sep + error.rstrip('0').rstrip('.')
        else:
            mstr = format(obj.magnitude, mspec).replace('+/-', '±')

        return mstr + ' ' + format(obj.units, uspec)


# Handles to original types
BaseQuantity = registry.Quantity
BaseUnit = registry.Unit
BaseMeasurement = registry.Measurement


registry.Quantity = Quantity
registry.Unit = Unit
registry.Measurement = Measurement

# Shortcuts for dimensionless quantities (must occur after subclassing of Unit)
dimensionless = registry.dimensionless


# Change default printing format
# noinspection PyShadowingNames,PyUnusedLocal
@pint.register_unit_format('gasify')
def format_custom(unit, registry: pint.registry.BaseRegistry, **options: typing.Any) -> str:
    unit_str = pint.formatter(
        unit.items(),
        as_ratio=True,
        single_denominator=False,
        product_fmt=" ",
        division_fmt="/",
        power_fmt="{}^{}",
        parentheses_fmt=r"({})",
        **options,
    )

    return unit_str


registry.default_format = 'g~#gasify'


# Handle pickle/unpickling by overwriting the built-in unit registry
pint.set_application_registry(registry)


def parse_unit(x: TParseUnit) -> Unit:
    """ Parse arbitrary input to a Unit from the registry.

    :param x: input str
    :return: parsed Unit
    """
    if isinstance(x, Unit):
        # Already a Unit
        return x

    if isinstance(x, Quantity):
        # Extract Unit, can sometimes occur when using values from pint
        return x.units

    if not isinstance(x, str):
        raise ParseError(f"Unsupported input type \"{type(x)}\"")

    if hasattr(registry, x):
        return getattr(registry, x)

    raise UnknownUnit(f"Unknown unit \"{x}\"")


def parse(x: TParseQuantity, to_unit: typing.Optional[TParseUnit] = None,
          mag_round: typing.Optional[int] = None) -> Quantity:
    """ Parse arbitrary input to a Quantity of specified unit.

    :param x: input str, number or Quantity
    :param to_unit: str or Unit to convert parsed values to
    :param mag_round: if specified, round the magnitude of the Quantity to mag_round places
    :return: parsed Quantity
    """
    if x is None:
        raise ParseError('Cannot convert NoneType to Quantity')

    # Parse unit
    if to_unit is not None:
        to_unit = parse_unit(to_unit)

    if not isinstance(x, Quantity):
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
                raise IncompatibleUnits(f"Unable to convert parsed quantity {x!s} to units {to_unit!s}") from ex
        else:
            x = Quantity(x.m_as(dimensionless), to_unit)

    # x = typing.cast(Quantity, x)

    if mag_round is not None:
        # Round resulting value
        x = round(x, mag_round)

    return x


def parse_magnitude(x: TParseQuantity, magnitude_unit: TParseUnit = None,
                    input_unit: typing.Optional[TParseUnit] = None) -> float:
    """ Shortcut method to parse as value, optionally converting to specified unit before returning the magnitude.

    :param x: input str, number or Quantity
    :param magnitude_unit: str or Unit to convert parsed values to before conversion to magnitude
    :param input_unit: str or Unit to convert parsed values to
    :return: Quantity magnitude as specified unit or as the parsed unit
    """
    if magnitude_unit is not None:
        magnitude_unit = parse_unit(magnitude_unit)

    if input_unit is None:
        # Assume default parsing unit is same as casting unit
        input_unit = magnitude_unit

    if magnitude_unit is not None:
        return parse(x, input_unit).m_as(magnitude_unit)
    else:
        return parse(x, input_unit).magnitude


def parse_timedelta(x: TParseTimeDelta) -> timedelta:
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

    return x_unit.to_timedelta()


def converter(to_unit: typing.Optional[TParseUnit] = None,
              optional: bool = False) -> typing.Callable[[typing.Optional[TParseQuantity]], Quantity]:
    """ Create wrapper for parse decorator with a pre-defined unit. Useful with the attrs library.

    :param to_unit: str or Unit to convert values to, defaults to unitless
    :param optional: if False
    :return:
    """
    to_unit = to_unit or dimensionless

    def f(x: typing.Optional[TParseQuantity]) -> typing.Optional[Quantity]:
        if x is None:
            if not optional:
                raise ParseError('Input to converter cannot be None')

            return None

        return parse(x, to_unit)

    return f


def return_converter(to_unit: TParseUnit, allow_none: bool = False):
    """ Decorator to convert returned result to a Quantity.

    :param to_unit:
    :param allow_none:
    :return:
    """
    to_unit = parse_unit(to_unit)

    def wrapper_decorator(func):
        @functools.wraps(func)
        def wrapper_result(*args: typing.Any, **kwargs: typing.Any):
            result = func(*args, **kwargs)

            if result is None:
                if not allow_none:
                    raise ValueError('Expected numeric result')

                return None

            if not isinstance(result, Quantity):
                result = parse(result, to_unit)
            else:
                result.ito(to_unit)

            return result

        return wrapper_result

    return wrapper_decorator
