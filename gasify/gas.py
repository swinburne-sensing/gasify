from __future__ import annotations

from abc import ABCMeta, abstractmethod
from collections import defaultdict
from collections.abc import Iterable as abc_Iterable
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, FrozenSet, Iterable, Iterator, Optional, Set, Tuple, Union

from plenary import storage

from gasify.unit import TParseQuantity, Quantity, dimensionless, parse, registry as unit_registry


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


class MolecularStructure(Enum):
    """ Molecular structure gas constants.

    Reference: MKS (https://www.mksinst.com/n/flow-measurement-control-frequently-asked-questions)
    """
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
    alias: FrozenSet[str] = field(default_factory=frozenset, kw_only=True)

    molecular_structure: Optional[MolecularStructure] = field(default=None, kw_only=True)
    specific_heat: Union[None, TParseQuantity, Quantity] = field(default=None, kw_only=True)
    density: Union[None, TParseQuantity, Quantity] = field(default=None, kw_only=True)

    _analyte: bool = field(default=True, kw_only=True)
    order: int = field(default=0, kw_only=True)

    def __post_init__(self):
        # Convert quantities
        if self.specific_heat is not None:
            object.__setattr__(self, 'specific_heat', parse(self.specific_heat, unit_registry.cal / unit_registry.g))

        if self.density is not None:
            object.__setattr__(self, 'density', parse(self.density, unit_registry.g / unit_registry.L))

    @property
    def analyte(self) -> bool:
        return self._analyte

    @property
    def gcf(self) -> Optional[float]:
        if self.molecular_structure is None or self.density is None or self.specific_heat is None:
            # Unable to calculate GCF
            return None

        return (0.3106 * self.molecular_structure.value /
                (self.density * self.specific_heat)).magnitude

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

        raise NotImplemented

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

        raise NotImplemented

    def __le__(self, other: Any) -> bool:
        # Sort by order then name
        if isinstance(other, self.__class__):
            return (self.order, self.name) <= (other.order, other.name)

        raise NotImplemented

    def __gt__(self, other: Any) -> bool:
        # Sort by order then name
        if isinstance(other, self.__class__):
            return (self.order, self.name) > (other.order, other.name)

        raise NotImplemented

    def __ge__(self, other: Any) -> bool:
        # Sort by order then name
        if isinstance(other, self.__class__):
            return (self.order, self.name) >= (other.order, other.name)

        raise NotImplemented

    def __str__(self) -> str:
        if self.symbol is not None:
            return self.symbol
        else:
            return self.name


@dataclass(frozen=True)
class CompoundConcentration:
    concentration: Union[TParseQuantity, Quantity]
    compound: Union[TParseQuantity, Compound]

    def __post_init__(self):
        object.__setattr__(self, 'concentration', parse(self.concentration).to_compact())

        if not self.concentration.is_compatible_with(dimensionless):
            raise ValueError('Gas concentration must be a dimensionless quantity')

        if self.concentration.m_as(dimensionless) < 0:
            raise ValueError('Gas concentration cannot be below zero')

    def __add__(self, other) -> Mixture:
        if isinstance(other, abc_Iterable):
            return Mixture(self, *other)
        else:
            return Mixture(self, other)

    __radd__ = __add__

    def __truediv__(self, other: Any) -> CompoundConcentration:
        if isinstance(other, (int, float, str, Quantity)):
            return CompoundConcentration(self.concentration / parse(other), self.compound)

        raise NotImplemented

    def __mul__(self, other: Any) -> CompoundConcentration:
        if not isinstance(other, (int, float, str, Quantity)):
            raise NotImplemented

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
            return self.concentration == parse(other)

        return object.__eq__(self, other)

    def __lt__(self, other: Any) -> bool:
        if isinstance(other, self.__class__):
            if self.compound == other.compound:
                return self.concentration < other.concentration
            else:
                return self.compound < other.compound
        elif isinstance(other, (int, float, str, Quantity)):
            return self.concentration < parse(other)

        raise NotImplemented

    def __le__(self, other):
        if isinstance(other, self.__class__):
            if self.compound == other.compound:
                return self.concentration <= other.concentration
            else:
                return self.compound <= other.compound
        elif isinstance(other, (int, float, str, Quantity)):
            return self.concentration <= parse(other)

        raise NotImplemented

    def __gt__(self, other):
        if isinstance(other, self.__class__):
            if self.compound == other.compound:
                return self.concentration > other.concentration
            else:
                return self.compound > other.compound
        elif isinstance(other, (int, float, str, Quantity)):
            return self.concentration > parse(other)

        raise NotImplemented

    def __ge__(self, other):
        if isinstance(other, self.__class__):
            if self.compound == other.compound:
                return self.concentration >= other.concentration
            else:
                return self.compound >= other.compound
        elif isinstance(other, (int, float, str, Quantity)):
            return self.concentration >= parse(other)

        raise NotImplemented

    def __str__(self) -> str:
        return f"{self.concentration} {self.compound}"


@dataclass(init=False, frozen=True)
class Mixture(_GasRegistryEntry):
    content: Tuple[CompoundConcentration]

    # Mixture name
    name: Optional[str] = field(default=None, kw_only=True)

    def __init__(self, *content, name: Optional[str] = None):
        # Combine parts by compound
        parts_lut = defaultdict(lambda: Quantity(0, dimensionless))

        for part in content:
            parts_lut[part.compound] += part.concentration

        # Order components by concentration with non-analytes last
        object.__setattr__(
            self,
            'content',
            tuple(
                sorted(
                    [CompoundConcentration(concentration, compound) for compound, concentration in parts_lut.items()]
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

    def __radd__(self, other: Any) -> Mixture:
        return self + other

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

        return (
            0.3106
            * sum(
                (
                    part.concentration * part.compound.molecular_structure.value
                    for part in self.content
                )
            )
            / sum(
                (
                    part.concentration * part.compound.density * part.compound.specific_heat
                    for part in self.content
                )
            )
        ).magnitude

    @property
    def total(self) -> Quantity:
        return sum(part.concentration for part in self.content)

    @property
    def registry_key(self) -> Union[str, Iterable[str]]:
        if self.name is not None:
            return self.name
        else:
            return str(self)

    def __str__(self) -> str:
        return ', '.join(map(str, self.content))

    @classmethod
    def auto_balance(cls, *content: CompoundConcentration, balance: Compound, name: Optional[str] = None) -> Mixture:
        content = list(content)
        balance_quantity = Quantity(1.0, dimensionless)

        for part in content:
            balance_quantity -= part.concentration

        content.append(CompoundConcentration(balance_quantity, balance))

        return Mixture(*content, name=name)


nitrogen = Compound(
    'Nitrogen',
    symbol='N_2',
    molecular_structure=MolecularStructure.DIATOMIC,
    specific_heat=0.2485,
    density=1.25,
    _analyte=False,
    order=2
)


oxygen = Compound(
    'Oxygen',
    symbol='O_2',
    molecular_structure=MolecularStructure.DIATOMIC,
    specific_heat=0.2193,
    density=1.427,
    _analyte=False,
    order=1
)


air = Mixture.auto_balance(CompoundConcentration(Quantity(21, unit_registry.percent), oxygen), balance=nitrogen,
                           name='Air')


class _GasRegistry(storage.Registry[_GasRegistryEntry]):
    pass


registry = _GasRegistry([
    Compound(
        'Acetone',
        symbol='(CH_3)_2CO',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.51,
        density=0.21
    ),
    Compound(
        'Ammonia',
        symbol='NH_3',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.492,
        density=0.76
    ),
    Compound(
        'Argon',
        symbol='Ar',
        molecular_structure=MolecularStructure.MONATOMIC,
        specific_heat=0.1244,
        density=1.782,
        _analyte=False,
        order=2
    ),
    Compound(
        'Arsine',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.1167,
        density=3.478
    ),
    Compound(
        'Bromine',
        symbol='Br_2',
        molecular_structure=MolecularStructure.DIATOMIC,
        specific_heat=0.0539,
        density=7.13
    ),
    Compound(
        'Carbon-dioxide',
        symbol='CO_2',
        molecular_structure=MolecularStructure.TRIATOMIC,
        specific_heat=0.2016,
        density=1.964
    ),
    Compound(
        'Carbon-monoxide',
        symbol='CO',
        molecular_structure=MolecularStructure.DIATOMIC,
        specific_heat=0.2488,
        density=1.25
    ),
    Compound(
        'Carbon-tetrachloride',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.1655,
        density=6.86
    ),
    Compound(
        'Carbon-tetraflouride',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.1654,
        density=3.926
    ),
    Compound(
        'Chlorine',
        symbol='Cl_2',
        molecular_structure=MolecularStructure.DIATOMIC,
        specific_heat=0.1144,
        density=3.163
    ),
    Compound(
        'Cyanogen',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.2613,
        density=2.322
    ),
    Compound(
        'Deuterium',
        symbol='H_2/D_2',
        molecular_structure=MolecularStructure.DIATOMIC,
        specific_heat=1.722,
        density=0.1799
    ),
    Compound(
        'Ethane',
        symbol='C_2H_6',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.4097,
        density=1.342
    ),
    Compound(
        'Fluorine',
        symbol='F_2',
        molecular_structure=MolecularStructure.DIATOMIC,
        specific_heat=0.1873,
        density=1.695
    ),
    Compound(
        'Helium',
        symbol='He',
        molecular_structure=MolecularStructure.MONATOMIC,
        specific_heat=1.241,
        density=0.1786,
        _analyte=False,
        order=2
    ),
    Compound(
        'Hexane',
        symbol='C_6H14',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.54,
        density=0.672
    ),
    Compound(
        'Hydrogen',
        symbol='H_2',
        molecular_structure=MolecularStructure.DIATOMIC,
        specific_heat=3.3852,
        density=0.0899
    ),
    Compound(
        'Hydrogen-chloride',
        symbol='HCl',
        molecular_structure=MolecularStructure.DIATOMIC,
        specific_heat=0.1912,
        density=1.627
    ),
    Compound(
        'Hydrogen-fluoride',
        symbol='HF',
        molecular_structure=MolecularStructure.DIATOMIC,
        specific_heat=0.3479,
        density=0.893
    ),
    Compound(
        'Methane',
        symbol='CH_4',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.5223,
        density=0.716
    ),
    Compound(
        'Neon',
        symbol='Ne',
        molecular_structure=MolecularStructure.MONATOMIC,
        specific_heat=0.246,
        density=0.9,
        _analyte=False,
        order=2
    ),
    Compound(
        'Nitric-oxide',
        symbol='NO',
        molecular_structure=MolecularStructure.DIATOMIC,
        specific_heat=0.2328,
        density=1.339
    ),
    Compound(
        'Nitric-oxides',
        symbol='NO_x'
    ),
    Compound(
        'Nitrogen-dioxide',
        symbol='NO_2',
        molecular_structure=MolecularStructure.TRIATOMIC,
        specific_heat=0.1933,
        density=2.052
    ),
    Compound(
        'Nitrous-oxide',
        symbol='N_2O',
        molecular_structure=MolecularStructure.TRIATOMIC,
        specific_heat=0.2088,
        density=1.964
    ),
    Compound(
        'Phosphine',
        symbol='PH_3',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.2374,
        density=1.517
    ),
    Compound(
        'Propane',
        symbol='C_3H_8',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.3885,
        density=1.967
    ),
    Compound(
        'Propylene',
        symbol='C_3H_6',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.3541,
        density=1.877
    ),
    Compound(
        'Sulfur hexaflouride',
        symbol='SF_6',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.1592,
        density=6.516
    ),
    Compound(
        'Xenon',
        symbol='Xe',
        molecular_structure=MolecularStructure.MONATOMIC,
        specific_heat=0.0378,
        density=5.858,
        _analyte=False,
        order=2
    ),
    nitrogen,
    oxygen
])
