from __future__ import annotations

import re
import typing
from dataclasses import dataclass, field
from enum import Enum

from experimentutil import storage

from gasify.unit import dimensionless, parse, registry as unit_registry, Quantity, TParseQuantity


class GasError(Exception):
    pass


class IncompatibleBalanceGases(GasError):
    pass


class MultipleBalanceGases(GasError):
    pass


class MolecularStructure(Enum):
    """ Molecular structure gas constants.

    Reference: MKS (https://www.mksinst.com/n/flow-measurement-control-frequently-asked-questions)
    """
    MONATOMIC = 1.03
    DIATOMIC = 1
    TRIATOMIC = 0.941
    POLYATOMIC = 0.88


@dataclass(frozen=True)
class Compound(storage.RegistryEntry):
    # Gas name
    name: str

    # Gas chemical symbol
    symbol: typing.Optional[str] = field(default=None)

    alias: typing.Set[str] = field(default_factory=set)

    molecular_structure: typing.Optional[MolecularStructure] = field(default=None, kw_only=True)
    specific_heat: typing.Union[None, TParseQuantity, Quantity] = field(default=None, kw_only=True)
    density: typing.Union[None, TParseQuantity, Quantity] = field(default=None, kw_only=True)

    # Analyte gas flag
    analyte: bool = field(default=True, kw_only=True)

    # Humidity carrier flag
    humid: bool = field(default=False, kw_only=True)

    def __post_init__(self):
        # Convert quantities
        if self.specific_heat is not None:
            object.__setattr__(self, 'specific_heat', parse(self.specific_heat, unit_registry.cal / unit_registry.g))

        if self.density is not None:
            object.__setattr__(self, 'density', parse(self.density, unit_registry.g / unit_registry.L))

    def __str__(self):
        if self.symbol is not None:
            return self.symbol
        else:
            return self.name

    @property
    def registry_key(self) -> typing.Union[str, typing.Iterable[str]]:
        key_set = {self.name}

        if self.symbol is not None:
            key_set.add(self.symbol)

        key_set.update(self.alias)

        return key_set


@dataclass(frozen=True)
class CompoundConcentration:
    # Actual concentration
    concentration: typing.Union[TParseQuantity, Quantity]

    # Gas type
    properties: Compound

    _max_concentration: typing.ClassVar[Quantity] = Quantity(1.0, dimensionless)

    def __post_init__(self):
        # Parse and scale concentration to look presentable
        object.__setattr__(self, 'concentration', parse(self.concentration).to_compact())

        if not self.concentration.is_compatible_with(dimensionless):
            raise ValueError('Gas concentration must be a dimensionless quantity')

        if self.concentration.m_as(dimensionless) < 0:
            raise ValueError('Gas concentration cannot be below zero')

    @property
    def gcf(self) -> float:
        """ Get gas correction factor for thermal mass flow controllers.

        :return: gas correction factor float
        """
        return (0.3106 * self.properties.molecular_structure.value /
                (self.properties.density * self.properties.specific_heat)).magnitude

    def __rmul__(self, other):
        return CompoundConcentration(other * self.concentration, self.properties)

    def __str__(self):
        return f"{self.concentration!s} {self.properties!s}"

    @classmethod
    def max_concentration(cls, compound: Compound) -> CompoundConcentration:
        return cls(cls._max_concentration, compound)


@dataclass(frozen=True)
class Mixture(storage.RegistryEntry):
    # Gases in mixture
    components: typing.List[CompoundConcentration]

    # Balance
    balance: CompoundConcentration

    _GAS_CONCENTRATION_PATTERN: typing.ClassVar[re.Pattern] = re.compile(r'^(\d+\.?\d*)\s?([%\w]+)\s+([\w\s-]+)')

    def __post_init__(self):
        # Check components and balance do not exceed 100% within a small tolerance (ppt)
        overall_quantity = self.balance.concentration

        for component in self.components:
            overall_quantity += component.concentration

        # Allow for small calculation error (1 ppt)
        if (overall_quantity.m_as(dimensionless) - 1) > 1e-12:
            raise ValueError(
                f"Gases in mixture sum to over 100% (components: {', '.join(map(str, self.components))}, "
                f"balance: {self.balance!s})")

        # Sort components by concentration
        object.__setattr__(self, 'components', sorted(self.components, key=lambda x: x.concentration, reverse=True))

    @property
    def is_analyte(self) -> bool:
        """ True if all gases in mixture are balance gases.

        :return:
        """
        return any((component.properties.analyte for component in self.components)) or self.balance.properties.analyte

    @property
    def is_humid(self) -> bool:
        """ True if no component in the gas carries humidity.

        :return:
        """
        return any((component.properties.humid for component in self.components)) or self.balance.properties.humid

    @property
    def humid_ratio(self) -> float:
        ratio = 0.0

        for component in self.components:
            if component.properties.humid:
                ratio += component.concentration

        return Quantity(ratio, dimensionless).m_as(dimensionless)

    @property
    def gcf(self) -> float:
        """ Gas correction factor used for mass flow control.

        :return: gas correction factor
        """
        mixture_component_list = self.components + [self.balance]

        error_list = []

        for component in mixture_component_list:
            if component.properties.molecular_structure is None:
                error_list.append(f"{component!s} missing molecular structure")

            if component.properties.density is None:
                error_list.append(f"{component!s} missing density")

            if component.properties.specific_heat is None:
                error_list.append(f"{component!s} missing specific heat")

        if len(error_list) > 0:
            raise ValueError(f"Cannot calculate GCF, missing one or more properties of components: "
                             f"{', '.join(error_list)}")

        return (
            0.3106 * sum(
                (c.concentration * c.properties.molecular_structure.value for c in mixture_component_list)
            ) / sum(
                (c.concentration * c.properties.density * c.properties.specific_heat for c in mixture_component_list)
            )
        ).magnitude

    @property
    def registry_key(self) -> typing.Union[str, typing.Iterable[str]]:
        pass

    def __str__(self):
        if len(self.components) > 0:
            return f"{', '.join(map(str, self.components))}"
        else:
            return str(self.balance)

    def __mul__(self, other):
        if isinstance(other, int):
            other = float(other)
        elif isinstance(other, Quantity):
            if not other.dimensionless:
                raise ValueError(f"Multiplication factor {other!s} must be dimensionless")

            # Cast to magnitude
            other = other.m_as(dimensionless)
        elif not isinstance(other, float):
            raise NotImplementedError(f"Cannot multiply {type(other)} by Mixture")

        scaled_gases = [other * gas for gas in self.components]

        return Mixture(scaled_gases, self.balance)

    def __rmul__(self, other):
        return self * other

    def __add__(self, other):
        gas_component_dict = {}

        if isinstance(other, Mixture):
            if self.balance.properties != other.balance.properties:
                raise IncompatibleBalanceGases(
                    f"Balance components {self.balance} and {other.balance} are incompatible"
                )

            for component in self.components + other.components:
                if component.properties in gas_component_dict:
                    gas_component_dict[component.properties] += component.concentration
                else:
                    gas_component_dict[component.properties] = component.concentration

            return Mixture.auto_balance(
                [CompoundConcentration(quantity, properties) for properties, quantity in gas_component_dict.items()],
                self.balance.properties
            )
        elif isinstance(other, CompoundConcentration):
            self.components.append(other)

            return self
        else:
            raise NotImplementedError(f"Cannot add {type(other)} to Mixture")

    def __radd__(self, other):
        return self + other

    @classmethod
    def auto_balance(cls, components: typing.Sequence[CompoundConcentration], balance: Compound) -> Mixture:
        """ Generate a gas mixture from a list of components and a balance compound.

        :param components:
        :param balance:
        :return:
        """
        balance_quantity = Quantity(1.0, dimensionless)

        for component in components:
            balance_quantity -= component.concentration

        return Mixture(list(components), CompoundConcentration(balance_quantity, balance))

    @classmethod
    def from_dict(cls, gas_mapping: typing.Mapping[str, str]):
        """

        :param gas_mapping:
        :return:
        """
        balance = None
        components = []

        for gas_name, concentration in gas_mapping.items():
            # Parse gas
            gas = registry[gas_name]

            # Parse concentration
            if concentration.lower() == 'balance' or concentration.lower() == '' or concentration is None:
                if balance is not None:
                    raise MultipleBalanceGases('Multiple balance gases specified')

                balance = gas
            else:
                components.append(CompoundConcentration(concentration, gas))

        return cls.auto_balance(components, balance)

    @classmethod
    def from_str(cls, gas_list_str: str) -> Mixture:
        gases: typing.List[CompoundConcentration] = []

        # Replace percentages for parsing
        gas_list = [x.strip() for x in gas_list_str.split(',')]

        # Parse components
        for gas_str in gas_list:
            gas_comp = cls._GAS_CONCENTRATION_PATTERN.match(gas_str)

            if gas_comp is not None:
                concentration = Quantity(gas_comp[1] + ' ' + gas_comp[2])
                gas_name = gas_comp[3].strip()
            else:
                concentration = Quantity(1, dimensionless)
                gas_name = gas_str

            if gas_name.lower() in registry:
                gas = registry[gas_name.lower()]
            else:
                raise ValueError(f"Unknown gas \"{gas_name}\" (from: \"{gas_str}\")")

            gases.append(CompoundConcentration(concentration, gas))

        return cls.auto_balance(gases[:-1], gases[-1].properties)


registry = storage.Registry([
    Compound(
        'Air',
        molecular_structure=MolecularStructure.DIATOMIC,
        specific_heat=0.24,
        density=1.293,
        analyte=False
    ),
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
        analyte=False
    ),
    Compound(
        'Arsine',
        molecular_structure=MolecularStructure.POLYATOMIC,
        specific_heat=0.1167,
        density=3.478
    ),
    Compound(
        'Bromine',
        symbol='BR_2',
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
        analyte=False
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
        analyte=False
    ),
    Compound(
        'Nitrogen',
        symbol='N_2',
        molecular_structure=MolecularStructure.DIATOMIC,
        specific_heat=0.2485,
        density=1.25,
        analyte=False
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
        'Oxygen',
        symbol='O_2',
        molecular_structure=MolecularStructure.DIATOMIC,
        specific_heat=0.2193,
        density=1.427,
        analyte=False
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
        analyte=False
    ),

    Compound(
        'Humid Air',
        analyte=False,
        humid=True
    ),
    Compound(
        'Humid Argon',
        symbol='Ar + H_2O',
        analyte=False,
        humid=True
    ),
    Compound(
        'Humid Nitrogen',
        symbol='N_2 + H_2O',
        analyte=False,
        humid=True
    )
])
