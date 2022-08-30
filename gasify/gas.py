from __future__ import annotations

from abc import ABCMeta, abstractmethod
from collections import defaultdict
from collections.abc import Iterable as abc_Iterable
from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, Iterable, Iterator, Optional, Set, Tuple, Union, cast

from plenary import storage

from gasify.unit import Quantity, dimensionless, parse, registry as unit_registry


_unit_specific_heat = unit_registry.cal / unit_registry.g
_unit_density = unit_registry.g / unit_registry.L


class _GasRegistryEntry(storage.RegistryEntry, metaclass=ABCMeta):
    @property
    def analyte(self) -> bool:
        return True

    @property
    @abstractmethod
    def gcf(self) -> Optional[float]:
        """ Get gas correction factor for thermal mass flow controllers.

        :return: gas correction factor
        """
        pass


# Molecular structure gas constants.
# Reference: MKS (https://www.mksinst.com/n/flow-measurement-control-frequently-asked-questions)
MONATOMIC = 1.03
DIATOMIC = 1
TRIATOMIC = 0.941
POLYATOMIC = 0.88


@dataclass(eq=False, frozen=True)
class Compound(_GasRegistryEntry):
    # Gas name
    name: str

    # Gas chemical symbol and other names
    symbol: Optional[str] = field(default=None)
    alias: FrozenSet[str] = field(default_factory=frozenset)

    molecular_structure: Optional[float] = field(default=None)
    specific_heat: Optional[Quantity] = field(default=None)
    density: Optional[Quantity] = field(default=None)

    _analyte: bool = field(default=True)
    order: int = field(default=0)

    @property
    def analyte(self) -> bool:
        return self._analyte

    @property
    def gcf(self) -> Optional[float]:
        if self.molecular_structure is None or self.density is None or self.specific_heat is None:
            # Unable to calculate GCF
            return None

        return cast(float, (0.3106 * self.molecular_structure / (self.density * self.specific_heat)).magnitude)

    @property
    def registry_key(self) -> Union[str, Iterable[str]]:
        key_set = {self.name}

        if self.symbol is not None:
            key_set.add(self.symbol)

        key_set.update(self.alias)

        return key_set

    def __mul__(self, other: Any) -> CompoundConcentration:
        if isinstance(other, (int, float, str, Quantity)):
            return CompoundConcentration(parse(other), self)

        return NotImplemented

    __rmul__ = __mul__

    def __hash__(self) -> int:
        # Equality depends only on name
        return hash(self.name)

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, self.__class__):
            return self.name == other.name

        return object.__eq__(self, other)

    def __lt__(self, other: Any) -> bool:
        # Sort by order then name
        if isinstance(other, self.__class__):
            return (self.order, self.name) < (other.order, other.name)

        return NotImplemented

    def __le__(self, other: Any) -> bool:
        # Sort by order then name
        if isinstance(other, self.__class__):
            return (self.order, self.name) <= (other.order, other.name)

        return NotImplemented

    def __gt__(self, other: Any) -> bool:
        # Sort by order then name
        if isinstance(other, self.__class__):
            return (self.order, self.name) > (other.order, other.name)

        return NotImplemented

    def __ge__(self, other: Any) -> bool:
        # Sort by order then name
        if isinstance(other, self.__class__):
            return (self.order, self.name) >= (other.order, other.name)

        return NotImplemented

    def __str__(self) -> str:
        if self.symbol is not None:
            return self.symbol
        else:
            return self.name


@dataclass(frozen=True)
class CompoundConcentration:
    concentration: Quantity
    compound: Compound

    def __post_init__(self) -> None:
        object.__setattr__(self, 'concentration', parse(self.concentration).to_compact())

        if not self.concentration.is_compatible_with(dimensionless):
            raise ValueError('Gas concentration must be a dimensionless quantity')

        if self.concentration.m_as(dimensionless) < 0:
            raise ValueError('Gas concentration cannot be below zero')

    def __add__(self, other: Any) -> Mixture:
        if isinstance(other, abc_Iterable):
            return Mixture(self, *other)
        else:
            return Mixture(self, other)

    __radd__ = __add__

    def __truediv__(self, other: Any) -> CompoundConcentration:
        if isinstance(other, (int, float, str, Quantity)):
            return CompoundConcentration(self.concentration / parse(other), self.compound)

        return NotImplemented

    def __mul__(self, other: Any) -> CompoundConcentration:
        if not isinstance(other, (int, float, str, Quantity)):
            return NotImplemented

        return CompoundConcentration(parse(other) * self.concentration, self.compound)

    __rmul__ = __mul__

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, self.__class__):
            # Allow some error for floating point math, equates to less than 0.001 ppt
            return self.compound == other.compound and abs(self.concentration.m_as(dimensionless) -
                                                           other.concentration.m_as(dimensionless)) <= 1e-15
        elif isinstance(other, Compound):
            return self.compound == other
        elif isinstance(other, (int, float, str, Quantity)):
            return bool(self.concentration == parse(other))

        return object.__eq__(self, other)

    def __lt__(self, other: Any) -> bool:
        if isinstance(other, self.__class__):
            if self.compound == other.compound:
                return bool(self.concentration < other.concentration)
            else:
                return self.compound < other.compound
        elif isinstance(other, (int, float, str, Quantity)):
            return bool(self.concentration < parse(other))

        return NotImplemented

    def __le__(self, other: Any) -> bool:
        if isinstance(other, self.__class__):
            if self.compound == other.compound:
                return bool(self.concentration <= other.concentration)
            else:
                return self.compound <= other.compound
        elif isinstance(other, (int, float, str, Quantity)):
            return bool(self.concentration <= parse(other))

        return NotImplemented

    def __gt__(self, other: Any) -> bool:
        if isinstance(other, self.__class__):
            if self.compound == other.compound:
                return bool(self.concentration > other.concentration)
            else:
                return bool(self.compound > other.compound)
        elif isinstance(other, (int, float, str, Quantity)):
            return bool(self.concentration > parse(other))

        return NotImplemented

    def __ge__(self, other: Any) -> bool:
        if isinstance(other, self.__class__):
            if self.compound == other.compound:
                return bool(self.concentration >= other.concentration)
            else:
                return bool(self.compound >= other.compound)
        elif isinstance(other, (int, float, str, Quantity)):
            return bool(self.concentration >= parse(other))

        return NotImplemented

    def __str__(self) -> str:
        return f"{self.concentration} {self.compound}"


@dataclass(init=False, frozen=True)
class Mixture(_GasRegistryEntry):
    content: Tuple[CompoundConcentration]

    # Mixture name
    name: Optional[str] = field(default=None, kw_only=True)

    def __init__(self, *content: CompoundConcentration, name: Optional[str] = None):
        # Combine parts by compound
        parts_lut: Dict[Compound, float] = defaultdict(float)

        for part in content:
            parts_lut[part.compound] += part.concentration.m_as(dimensionless)

        # Order components by concentration with non-analytes last
        object.__setattr__(
            self,
            'content',
            tuple(
                sorted(
                    [
                        CompoundConcentration(Quantity(concentration, dimensionless), compound)
                        for compound, concentration in parts_lut.items()
                    ]
                )
            )
        )
        object.__setattr__(self, 'name', name)

    def __add__(self, other: Any) -> Mixture:
        # Generate mixture from parts
        if isinstance(other, abc_Iterable):
            return Mixture(*self, *other)
        else:
            return Mixture(*self, other)

    __radd__ = __add__

    def __iter__(self) -> Iterator[CompoundConcentration]:
        return iter(self.content)

    def __len__(self) -> int:
        return len(self.content)

    def __contains__(self, item: Any) -> bool:
        if isinstance(item, Compound):
            return item in self.compounds
        else:
            return item in self.content

    @property
    def analyte(self) -> bool:
        return any(part.compound.analyte for part in self.content)

    @property
    def compounds(self) -> Set[Compound]:
        return set(part.compound for part in self.content)

    @property
    def gcf(self) -> float:
        error_list = []

        for part in self.content:
            if part.compound.gcf is None:
                error_list.append(f"{part!s} missing properties")

        if len(error_list) > 0:
            raise ValueError(f"Cannot calculate GCF, missing one or more properties of components: "
                             f"{', '.join(error_list)}")

        return cast(
            float,
            (
                0.3106
                * sum(
                    part.concentration * part.compound.molecular_structure
                    for part in self.content
                )
                / sum(
                    part.concentration * part.compound.density * part.compound.specific_heat
                    for part in self.content
                )
            ).magnitude
        )

    @property
    def total(self) -> Quantity:
        return Quantity(
            sum(part.concentration.m_as(dimensionless) for part in self.content),
            dimensionless
        ).to_compact()

    @property
    def registry_key(self) -> Union[str, Iterable[str]]:
        if self.name is not None:
            return self.name
        else:
            return str(self)

    def __str__(self) -> str:
        return ', '.join(map(str, self.content))

    def normalise(self) -> Mixture:
        total = self.total

        return Mixture(*(part / total for part in self.content), name=self.name)

    @classmethod
    def auto_balance(cls, *content: CompoundConcentration, balance: Compound, name: Optional[str] = None) -> Mixture:
        content_list = list(content)
        balance_quantity = Quantity(1.0, dimensionless)

        for part in content:
            balance_quantity -= part.concentration

        content_list.append(CompoundConcentration(balance_quantity, balance))

        return Mixture(*content_list, name=name)


nitrogen = Compound(
    'Nitrogen',
    symbol='N_2',
    molecular_structure=DIATOMIC,
    specific_heat=Quantity(0.2485, _unit_specific_heat),
    density=Quantity(1.25, _unit_density),
    _analyte=False,
    order=2
)


oxygen = Compound(
    'Oxygen',
    symbol='O_2',
    molecular_structure=DIATOMIC,
    specific_heat=Quantity(0.2193, _unit_specific_heat),
    density=Quantity(1.427, _unit_density),
    _analyte=False,
    order=1
)


# Based off instrument-grade air from BOC
air = Mixture.auto_balance(
    0.21 * oxygen,
    balance=nitrogen,
    name='Air'
)


class _GasRegistry(storage.Registry[_GasRegistryEntry]):
    pass


registry = _GasRegistry([
    Compound(
        'Acetone',
        symbol='(CH_3)_2CO',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.51, _unit_specific_heat),
        density=Quantity(0.21, _unit_density)
    ),
    Compound(
        'Ammonia',
        symbol='NH_3',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.492, _unit_specific_heat),
        density=Quantity(0.76, _unit_density)
    ),
    Compound(
        'Argon',
        symbol='Ar',
        molecular_structure=MONATOMIC,
        specific_heat=Quantity(0.1244, _unit_specific_heat),
        density=Quantity(1.782, _unit_density),
        _analyte=False,
        order=2
    ),
    Compound(
        'Arsine',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.1167, _unit_specific_heat),
        density=Quantity(3.478, _unit_density)
    ),
    Compound(
        'Bromine',
        symbol='Br_2',
        molecular_structure=DIATOMIC,
        specific_heat=Quantity(0.0539, _unit_specific_heat),
        density=Quantity(7.13, _unit_density)
    ),
    Compound(
        'Carbon-dioxide',
        symbol='CO_2',
        molecular_structure=TRIATOMIC,
        specific_heat=Quantity(0.2016, _unit_specific_heat),
        density=Quantity(1.964, _unit_density)
    ),
    Compound(
        'Carbon-monoxide',
        symbol='CO',
        molecular_structure=DIATOMIC,
        specific_heat=Quantity(0.2488, _unit_specific_heat),
        density=Quantity(1.25, _unit_density)
    ),
    Compound(
        'Carbon-tetrachloride',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.1655, _unit_specific_heat),
        density=Quantity(6.86, _unit_density)
    ),
    Compound(
        'Carbon-tetraflouride',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.1654, _unit_specific_heat),
        density=Quantity(3.926, _unit_density)
    ),
    Compound(
        'Chlorine',
        symbol='Cl_2',
        molecular_structure=DIATOMIC,
        specific_heat=Quantity(0.1144, _unit_specific_heat),
        density=Quantity(3.163, _unit_density)
    ),
    Compound(
        'Cyanogen',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.2613, _unit_specific_heat),
        density=Quantity(2.322, _unit_density)
    ),
    Compound(
        'Deuterium',
        symbol='H_2/D_2',
        molecular_structure=DIATOMIC,
        specific_heat=Quantity(1.722, _unit_specific_heat),
        density=Quantity(0.1799, _unit_density)
    ),
    Compound(
        'Ethane',
        symbol='C_2H_6',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.4097, _unit_specific_heat),
        density=Quantity(1.342, _unit_density)
    ),
    Compound(
        'Fluorine',
        symbol='F_2',
        molecular_structure=DIATOMIC,
        specific_heat=Quantity(0.1873, _unit_specific_heat),
        density=Quantity(1.695, _unit_density)
    ),
    Compound(
        'Helium',
        symbol='He',
        molecular_structure=MONATOMIC,
        specific_heat=Quantity(1.241, _unit_specific_heat),
        density=Quantity(0.1786, _unit_density),
        _analyte=False,
        order=2
    ),
    Compound(
        'Hexane',
        symbol='C_6H14',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.54, _unit_specific_heat),
        density=Quantity(0.672, _unit_density)
    ),
    Compound(
        'Hydrogen',
        symbol='H_2',
        molecular_structure=DIATOMIC,
        specific_heat=Quantity(3.3852, _unit_specific_heat),
        density=Quantity(0.0899, _unit_density)
    ),
    Compound(
        'Hydrogen-chloride',
        symbol='HCl',
        molecular_structure=DIATOMIC,
        specific_heat=Quantity(0.1912, _unit_specific_heat),
        density=Quantity(1.627, _unit_density)
    ),
    Compound(
        'Hydrogen-fluoride',
        symbol='HF',
        molecular_structure=DIATOMIC,
        specific_heat=Quantity(0.3479, _unit_specific_heat),
        density=Quantity(0.893, _unit_density)
    ),
    Compound(
        'Methane',
        symbol='CH_4',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.5223, _unit_specific_heat),
        density=Quantity(0.716, _unit_density)
    ),
    Compound(
        'Neon',
        symbol='Ne',
        molecular_structure=MONATOMIC,
        specific_heat=Quantity(0.246, _unit_specific_heat),
        density=Quantity(0.9, _unit_density),
        _analyte=False,
        order=2
    ),
    Compound(
        'Nitric-oxide',
        symbol='NO',
        molecular_structure=DIATOMIC,
        specific_heat=Quantity(0.2328, _unit_specific_heat),
        density=Quantity(1.339, _unit_density)
    ),
    Compound(
        'Nitric-oxides',
        symbol='NO_x'
    ),
    Compound(
        'Nitrogen-dioxide',
        symbol='NO_2',
        molecular_structure=TRIATOMIC,
        specific_heat=Quantity(0.1933, _unit_specific_heat),
        density=Quantity(2.052, _unit_density)
    ),
    Compound(
        'Nitrous-oxide',
        symbol='N_2O',
        molecular_structure=TRIATOMIC,
        specific_heat=Quantity(0.2088, _unit_specific_heat),
        density=Quantity(1.964, _unit_density)
    ),
    Compound(
        'Phosphine',
        symbol='PH_3',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.2374, _unit_specific_heat),
        density=Quantity(1.517, _unit_density)
    ),
    Compound(
        'Propane',
        symbol='C_3H_8',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.3885, _unit_specific_heat),
        density=Quantity(1.967, _unit_density)
    ),
    Compound(
        'Propylene',
        symbol='C_3H_6',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.3541, _unit_specific_heat),
        density=Quantity(1.877, _unit_density)
    ),
    Compound(
        'Sulfur hexaflouride',
        symbol='SF_6',
        molecular_structure=POLYATOMIC,
        specific_heat=Quantity(0.1592, _unit_specific_heat),
        density=Quantity(6.516, _unit_density)
    ),
    Compound(
        'Xenon',
        symbol='Xe',
        molecular_structure=MONATOMIC,
        specific_heat=Quantity(0.0378, _unit_specific_heat),
        density=Quantity(5.858, _unit_density),
        _analyte=False,
        order=2
    ),
    nitrogen,
    oxygen,
    air
])
