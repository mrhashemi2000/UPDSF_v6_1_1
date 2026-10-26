"""
================================================================================
UNIFIED PREBIOTIC DNA SELECTION FRAMEWORK (UPDSF) v6.1.1
Chemical Filtering and T-over-U Enrichment Model
Mass-Consistent | Kinetically-Consistent | Provenance-Tracked
================================================================================

AUTHOR: Seyed Mohammad Reza Hashemi (Reza Hashemi) Intelligence-Augmented (IA)
VERSION: 6.1.1
License & Copyright
Copyright ©️ 2026 Seyed Mohammad Reza Hashemi
This work is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License (CC BY-NC-ND 4.0).

ORCID: 0009-0002-0645-5180



CHANGELOG v6.1.0 -> v6.1.1
-------------------------
1.  FIX (CRITICAL): main() called analyzer.free_monomer_deam_scale_sweep()
    and analyzer.clay_rate_sweep() but neither method existed on
    TwoDSensitivityAnalyzer, so the script crashed with AttributeError
    at the end of every run. Both methods are now implemented.
2.  NEW: plot_free_deam_scale_sweep() — CSI / real / drift / free-vs-DNA
    T/(U+T) / yield versus FREE_MONOMER_DEAM_SCALE.
3.  NEW: plot_clay_rate_sweep() — CSI decomposition / compartmentalized
    fraction / yield versus CLAY_ADSORPTION_RATE.
4.  NEW: main() now saves the two new sweep plots alongside the others.
5.  No physics, parameters, or grid values changed. All v6.1.0 numerical
    results are reproduced exactly.

CHANGELOG v6.0.1 -> v6.1.0
-------------------------
1.  VENT DEFAULT: photon_flux default is now 0.0 (dark hydrothermal vent).
    UV machinery is retained for surface-scenario comparison only.
2.  NEW: CSI trajectory tracking — peak CSI, time-to-peak CSI, CSI at end,
    CSI regime classification (rising / peaked / collapsed / steady).
3.  NEW: time_horizon_sweep() — runs the optimum at multiple horizons
    (240 h, 2400 h, 10000 h) and reports CSI, real, drift, yield, survival.
4.  NEW: ea_T_sweep() — sensitivity of the stability ratio and CSI to the
    assumed Ea_T in [28, 36] kcal/mol; exposes the ~2-order-of-magnitude
    uncertainty in the 357x extrapolated k_U/k_T ratio.
5.  NEW: plot_csi_trajectory() — CSI, real, drift, free T/(U+T) vs time
    at the optimum, with peak markers.
6.  NEW: plot_time_horizon_sweep() — bar/line plot of CSI vs horizon.
7.  NEW: pool_balance_diagnostic() — explicit influx-vs-destruction balance
    for the free pool, so the long-time regime is stated, not guessed.
8.  FIX: stability_ratio is now reported with an explicit "extrapolated"
    label and a sensitivity band from the Ea_T sweep.
9.  FIX: save_history_to_csv now also records CSI, real, drift,
    free T/(U+T), and pool-balance ratios per timestep.
10. FIX: optimal_simulation() now accepts max_time_hours override and
    records CSI trajectory statistics in the annotation box.
11. FIX: main() now runs the vent (zero-UV) case as the primary result,
    and the surface (1e16 flux) case as a labeled comparison.
12. All plots at 300 DPI.

PHILOSOPHY
----------
Every parameter carries a Source object with a `kind` field:
  "literature" -> value comes from a specific paper
  "assumed"    -> heuristic, no single source
  "imposed"    -> a modeling choice, not a measurement
  "derived"    -> computed from other parameters

ARRHENIUS CONVENTION
--------------------
k(T) = k_ref * exp(-Ea/R * (1/T - 1/T_ref))

STOCHASTIC CONVENTION
---------------------
For a destruction channel with rate k acting on N molecules over dt:
    n_destroyed ~ Binomial(N, 1 - exp(-k * dt))

MASS-CONSERVATION IDENTITY
--------------------------
For each base b:
    mass_consumed[b] == dna_residues[b]
                      + mass_hydrolyzed[b]
                      + mass_uv_destroyed[b]

TWO ENVIRONMENTAL CHANNELS
--------------------------
CSI = Real T-preferential enrichment * Pool-composition drift factor

VENT VS SURFACE
---------------
Default = vent (photon_flux = 0). At the optimum, UV contributes <2% to
CSI, so the vent and surface results coincide to within stochastic noise.
The surface case is retained as an explicit labeled comparison.

REFERENCES
----------
Lindahl 1993; Cleaves 2010; Shen 1994; Ravanat & Cadet 1995;
Ferris 1996; Deamer 2017; Szostak 2010; Eigen 1993; Tinland 1997;
Hagerman 1988; Miller 1953; Ertem & Ferris 1996;
Frederico, Kunkel & Shaw 1990; IAPWS R7-97 (approximated).
"""

# ====================================================================
# IMPORTS
# ====================================================================

import ast
import sys
import json
import os
import random
import warnings
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)


# ====================================================================
# SECTION 0: PROVENANCE TRACKING
# ====================================================================

_ALLOWED_SOURCE_KINDS = {"literature", "assumed", "imposed", "derived", "extrapolated"}


@dataclass
class Source:
    value: float
    units: str
    citation: str
    kind: str
    note: str = ""

    def __post_init__(self):
        if self.kind not in _ALLOWED_SOURCE_KINDS:
            raise ValueError(
                f"Invalid Source.kind={self.kind!r}. "
                f"Allowed: {sorted(_ALLOWED_SOURCE_KINDS)}"
            )

    def __str__(self):
        return (
            f"{self.value:g} {self.units}  [{self.kind}]  "
            f"{self.citation}  {self.note}"
        )


def print_provenance_table(sources: Dict[str, Source], title: str):
    print("\n" + "=" * 120)
    print(f"PARAMETER PROVENANCE: {title}")
    print("=" * 120)
    print(f"{'Parameter':<40} {'Value':>16} {'Units':<14} {'Kind':<12} Citation")
    print("-" * 120)
    for name, s in sources.items():
        cit = s.citation if len(s.citation) <= 60 else s.citation[:57] + "..."
        print(
            f"{name:<40} {s.value:>16.6g} {s.units:<14} "
            f"{s.kind:<12} {cit}"
        )
    print("=" * 120 + "\n")


def provenance_summary(sources: Dict[str, Source]) -> Dict[str, int]:
    return {k: sum(1 for s in sources.values() if s.kind == k)
            for k in sorted(_ALLOWED_SOURCE_KINDS)}


# ====================================================================
# SECTION 1: PHYSICAL CONSTANTS (CODATA 2018)
# ====================================================================

class PhysicalConstants:
    kB = 1.380649e-23
    NA = 6.02214076e23
    R = 8.314462618
    R_kcal = 0.001987204
    H_PLANCK = 6.62607015e-34
    T_REF = 298.15


# ====================================================================
# SECTION 2: WATER CHEMISTRY
# ====================================================================

class WaterChemistry:
    @staticmethod
    def pKw(T_C: float) -> float:
        return 14.957 - 0.0421 * T_C + 1.6e-4 * T_C ** 2

    @staticmethod
    def OH_from_pH(pH_in_situ: float, T_C: float) -> float:
        return 10.0 ** (pH_in_situ - WaterChemistry.pKw(T_C))

    @staticmethod
    def pH_25_to_in_situ(
        pH_25: float, T_C: float, dH_ionization_kcal: float = 13.0
    ) -> float:
        R_kcal = PhysicalConstants.R_kcal
        T1 = 298.15
        T2 = T_C + 273.15
        return pH_25 - (dH_ionization_kcal / (2.303 * R_kcal)) * (
            1.0 / T2 - 1.0 / T1
        )


# ====================================================================
# SECTION 3: POLYMER PHYSICS
# ====================================================================

class PolymerPhysics:
    PERSISTENCE_SSRNA: Source = Source(2.0, "nm", "Tinland 1997", "literature")
    PERSISTENCE_SSDNA: Source = Source(1.5, "nm", "Tinland 1997", "literature")
    PERSISTENCE_DSDNA: Source = Source(50.0, "nm", "Hagerman 1988", "literature")
    CONTOUR_PER_BASE_RNA: Source = Source(
        0.59, "nm/base", "Hagerman 1988", "literature"
    )
    CONTOUR_PER_BASE_DNA: Source = Source(
        0.56, "nm/base", "Hagerman 1988", "literature"
    )

    @classmethod
    def provenance(cls) -> Dict[str, Source]:
        return {
            "PERSISTENCE_SSRNA": cls.PERSISTENCE_SSRNA,
            "PERSISTENCE_SSDNA": cls.PERSISTENCE_SSDNA,
            "PERSISTENCE_DSDNA": cls.PERSISTENCE_DSDNA,
            "CONTOUR_PER_BASE_RNA": cls.CONTOUR_PER_BASE_RNA,
            "CONTOUR_PER_BASE_DNA": cls.CONTOUR_PER_BASE_DNA,
        }


# ====================================================================
# SECTION 4: LIPID MEMBRANE
# ====================================================================

class LipidMembrane:
    CMC_DEFAULT: Source = Source(
        1e-3, "M", "Deamer 2017 (order of magnitude)", "literature"
    )
    T_OPT_DEFAULT: Source = Source(65.0, "C", "Deamer 2017 (qualitative)", "assumed")
    PH_OPT_DEFAULT: Source = Source(8.0, "-", "Deamer 2017 (qualitative)", "assumed")
    EFFICIENCY_DEFAULT: Source = Source(0.85, "-", "assumed", "assumed")
    PROTECTION_FACTOR_DEFAULT: Source = Source(4.5, "-", "assumed", "assumed")

    PARTITION_COEFFICIENTS: Dict[str, Source] = {
        'T': Source(1.8, "-", "Hansch & Leo 1995", "literature"),
        'U': Source(1.2, "-", "Hansch & Leo 1995", "literature"),
        'C': Source(0.9, "-", "Hansch & Leo 1995", "literature"),
        'A': Source(1.4, "-", "Hansch & Leo 1995", "literature"),
        'G': Source(1.1, "-", "Hansch & Leo 1995", "literature"),
        'X': Source(1.0, "-", "assumed", "assumed"),
    }

    def __init__(
        self,
        temperature_C=76.0,
        pH_in_situ=9.5,
        lipid_conc=0.05,
        lipid_chain_length=18,
        ionic_strength_M=0.1,
        Ca_mM=1.0,
    ):
        self.temperature_C = temperature_C
        self.pH_in_situ = pH_in_situ
        self.lipid_conc = lipid_conc
        self.lipid_chain_length = lipid_chain_length
        self.ionic_strength_M = ionic_strength_M
        self.Ca_mM = Ca_mM
        self.vesicle_fraction = self._calculate_vesicle_fraction()
        self.protection_active = self.vesicle_fraction > 0.1

    def _calculate_vesicle_fraction(self) -> float:
        T = self.temperature_C
        pH = self.pH_in_situ
        C = self.lipid_conc
        L = self.lipid_chain_length
        I = self.ionic_strength_M
        Ca = self.Ca_mM

        T_opt = self.T_OPT_DEFAULT.value
        pH_opt = self.PH_OPT_DEFAULT.value
        cmc = self.CMC_DEFAULT.value

        f_T = np.exp(-0.5 * ((T - T_opt) / 25.0) ** 2)
        f_pH = np.exp(-0.5 * ((pH - pH_opt) / 1.5) ** 2)
        f_conc = min(1.0, C / (cmc * 5.0))
        f_chain = np.exp(-0.5 * ((L - 18) / 4.0) ** 2)
        f_ionic = 1.0 / (1.0 + 2.0 * I)
        f_Ca = 1.0 - np.exp(-Ca / 2.0)

        return float(
            f_T * f_pH * f_conc * f_chain * f_ionic * f_Ca
            * self.EFFICIENCY_DEFAULT.value
        )

    def get_protection_factor(self, base: str = 'T') -> float:
        if not self.protection_active:
            return 1.0
        return self.PROTECTION_FACTOR_DEFAULT.value

    def get_partition_coefficient(self, base: str) -> float:
        s = self.PARTITION_COEFFICIENTS.get(base)
        return s.value if s is not None else 1.0

    @classmethod
    def provenance(cls) -> Dict[str, Source]:
        d = {
            "CMC_DEFAULT": cls.CMC_DEFAULT,
            "T_OPT_DEFAULT": cls.T_OPT_DEFAULT,
            "PH_OPT_DEFAULT": cls.PH_OPT_DEFAULT,
            "EFFICIENCY_DEFAULT": cls.EFFICIENCY_DEFAULT,
            "PROTECTION_FACTOR_DEFAULT": cls.PROTECTION_FACTOR_DEFAULT,
        }
        for b, s in cls.PARTITION_COEFFICIENTS.items():
            d[f"PARTITION_{b}"] = s
        return d


# ====================================================================
# SECTION 5: EMPIRICAL RATES
# ====================================================================

class EmpiricalRates:
    """
    Kinetic parameters. UNIT CONVENTION: masses in RESIDUES.
    ARRHENIUS: k(T) = k_ref * exp(-Ea/R * (1/T - 1/T_ref))

    NOTE ON C: Cytosine degrades via deamination to uracil, handled by
    get_deamination_rate_C(). C is NOT in get_hydrolysis_rates().

    NOTE ON EXTRAPOLATION (v6.1.0): The reference rates K_U_REF and
    K_T_REF are 25 C literature values. The activation energies Ea_U
    and Ea_T are from different studies (Cleaves 2010 / Lindahl 1993
    for U; Lindahl 1993 for T). The 90 C / pH 11 ratio of ~357x is an
    ARRHENIUS EXTRAPOLATION, not a measurement. It is sensitive to the
    assumed Ea_T (see ea_T_sweep).
    """

    K_U_REF: Source = Source(
        1.3e-10, "s^-1", "Lindahl 1993 (RNA hydrolysis at 25 C)", "literature"
    )
    K_T_REF: Source = Source(
        1.1e-13, "s^-1", "Lindahl 1993 (DNA hydrolysis at 25 C)", "literature"
    )
    K_A_REF: Source = Source(1.0e-11, "s^-1", "estimated from purine lability", "assumed")
    K_G_REF: Source = Source(5.0e-12, "s^-1", "estimated from purine lability", "assumed")

    Ea_U: Source = Source(27.0, "kcal/mol", "Cleaves 2010; Lindahl 1993", "literature")
    Ea_T: Source = Source(
        32.0, "kcal/mol",
        "Lindahl 1993 (extrapolated to high T; not measured at 90 C)",
        "extrapolated",
    )
    Ea_A: Source = Source(29.0, "kcal/mol", "Shen 1994", "literature")
    Ea_G: Source = Source(31.0, "kcal/mol", "assumed", "assumed")

    T_HYDRO_REF: Source = Source(25.0, "C", "Lindahl 1993 (experimental T)", "literature")

    K_POLY_REF: Source = Source(0.022, "s^-1", "Ferris 1996 (rate at 65 C)", "literature")
    T_POLY_REF: Source = Source(65.0, "C", "Ferris 1996 (experimental T)", "literature")
    Ea_poly: Source = Source(19.5, "kcal/mol", "Ferris 1996", "literature")

    K_C_DEAM_REF: Source = Source(
        1.3e-8, "s^-1", "Lindahl 1993 (dsDNA at 80 C)", "literature"
    )
    T_C_DEAM_REF: Source = Source(80.0, "C", "Lindahl 1993 (experimental T)", "literature")
    Ea_C_deam: Source = Source(
        28.0, "kcal/mol", "derived from Lindahl 1993 dsDNA", "derived"
    )
    C_DEAM_PH_ACCEL_MAX: Source = Source(100.0, "-", "assumed (bounded)", "assumed")
    RATIO_5mC_C: Source = Source(2.2, "-", "Shen 1994, Table 1 (dsDNA, 37C)", "literature")
    SSDNA_DEAM_ENHANCEMENT_C: Source = Source(
        140.0, "-",
        "Frederico, Kunkel & Shaw 1990 (dsDNA->ssDNA for C)",
        "literature"
    )
    SSDNA_DEAM_ENHANCEMENT_G: Source = Source(
        10.0, "-", "estimated (dsDNA->ssDNA for G)", "assumed"
    )

    K_G_DEAM_REF: Source = Source(
        1.0e-11, "s^-1", "Lindahl 1993 (dsDNA at 80 C)", "literature"
    )
    T_G_DEAM_REF: Source = Source(80.0, "C", "Lindahl 1993 (experimental T)", "literature")
    Ea_G_deam: Source = Source(30.0, "kcal/mol", "assumed", "assumed")

    FREE_MONOMER_DEAM_SCALE: Source = Source(
        2.0, "-", "assumed (free base faster than dsDNA; bounded)", "assumed"
    )
    K_FREE_C_DEAM_REF: Source = Source(
        1.3e-8 * 2.0, "s^-1",
        "derived: 2 * K_C_DEAM_REF (free monomer enhancement)",
        "derived"
    )
    K_FREE_G_DEAM_REF: Source = Source(
        1.0e-11 * 2.0, "s^-1",
        "derived: 2 * K_G_DEAM_REF (free monomer enhancement)",
        "derived"
    )
    Ea_free_C_deam: Source = Source(28.0, "kcal/mol", "same as Ea_C_deam (assumed)", "assumed")
    Ea_free_G_deam: Source = Source(30.0, "kcal/mol", "same as Ea_G_deam (assumed)", "assumed")

    FREE_MONOMER_FLOOR: Source = Source(
        10.0, "residues",
        "assumed (numerical floor to avoid pool going negative)",
        "assumed"
    )

    UV_ABS_CROSS_SECTION_REF_CM2: Source = Source(
        4.0e-20, "cm^2",
        "assumed order-of-magnitude nucleobase absorption cross-section",
        "assumed",
    )
    UV_CROSS_SECTION_U: Source = Source(1.0, "relative", "assumed (reference; U as 1.0)", "assumed")
    UV_CROSS_SECTION_T: Source = Source(
        0.28, "relative",
        "Ravanat & Cadet 1995 (T is ~3.6x more UV-stable than U)",
        "literature"
    )
    UV_CROSS_SECTION_C: Source = Source(1.5, "relative", "Ravanat & Cadet 1995", "literature")
    UV_CROSS_SECTION_A: Source = Source(0.7, "relative", "Ravanat & Cadet 1995", "literature")
    UV_CROSS_SECTION_G: Source = Source(0.75, "relative", "Ravanat & Cadet 1995", "literature")
    UV_QUANTUM_YIELD: Source = Source(1e-3, "-", "assumed (order of magnitude)", "assumed")
    UV_PEAK_WAVELENGTH_NM: Source = Source(
        260.0, "nm", "standard nucleic-acid absorption peak", "literature"
    )
    UV_BANDWIDTH_NM: Source = Source(15.0, "nm", "assumed (Gaussian sigma around peak)", "assumed")

    pKa_U: Source = Source(9.43, "-", "Shapiro 1999", "literature")
    pKa_T: Source = Source(9.79, "-", "Lindahl 1993", "literature")
    pKa_C: Source = Source(9.65, "-", "Shapiro 1999", "literature")
    pKa_A: Source = Source(9.50, "-", "Shen 1994", "literature")
    pKa_G: Source = Source(9.25, "-", "assumed", "assumed")
    pKa_phosphate: Source = Source(6.7, "-", "Ferris 1996", "literature")

    BASE_CATALYSIS_FACTOR_U: Source = Source(1.0, "-", "assumed", "assumed")
    BASE_CATALYSIS_FACTOR_T: Source = Source(0.75, "-", "assumed", "assumed")
    BASE_CATALYSIS_FACTOR_A: Source = Source(0.9, "-", "assumed", "assumed")
    BASE_CATALYSIS_FACTOR_G: Source = Source(0.85, "-", "assumed", "assumed")

    PH_FACTOR_DEPROTONATION_GAIN: Source = Source(2.5, "-", "imposed (v4.x heuristic)", "imposed")
    PH_FACTOR_PH_SLOPE: Source = Source(0.1, "-", "imposed (v4.x heuristic)", "imposed")
    PH_FACTOR_MIN: Source = Source(0.3, "-", "imposed (lower cap)", "imposed")
    PH_FACTOR_MAX: Source = Source(20.0, "-", "imposed (upper cap)", "imposed")
    BASE_CAT_OH_GAIN: Source = Source(100.0, "-", "imposed (v4.x heuristic)", "imposed")
    BASE_CAT_MAX: Source = Source(50.0, "-", "imposed (upper cap)", "imposed")
    ACID_FACTOR_SLOPE: Source = Source(0.2, "-", "imposed (v4.x heuristic)", "imposed")
    ACID_FACTOR_MAX: Source = Source(3.0, "-", "imposed (upper cap)", "imposed")

    K_POLY_CLAY_ENHANCEMENT: Source = Source(
        1.9, "-", "Ferris 1996 (order of magnitude)", "literature"
    )
    FIDELITY_BASE: Source = Source(0.93, "-", "Eigen 1993 (order of magnitude)", "literature")
    FIDELITY_ADVANTAGE_GAIN: Source = Source(
        0.25, "-", "imposed (v4.x heuristic; positive feedback on T)", "imposed"
    )

    CLAY_PROTECTION_FACTOR: Source = Source(
        9.5, "-", "Ertem & Ferris 1996 (order of magnitude)", "literature"
    )
    CLAY_ADSORPTION_RATE: Source = Source(5e-3, "s^-1", "assumed", "assumed")
    CLAY_DESORPTION_RATE: Source = Source(5e-4, "s^-1", "assumed", "assumed")

    INITIAL_U_RESIDUES: Source = Source(
        830_000, "residues", "Miller 1953; Orgel 2004", "literature"
    )
    INITIAL_T_RESIDUES: Source = Source(
        170_000, "residues", "Miller 1953; Orgel 2004", "literature"
    )
    INITIAL_C_RESIDUES: Source = Source(
        650_000, "residues", "Miller 1953; Orgel 2004", "literature"
    )
    INITIAL_A_RESIDUES: Source = Source(
        420_000, "residues", "Miller 1953; Orgel 2004", "literature"
    )
    INITIAL_G_RESIDUES: Source = Source(
        500_000, "residues", "Miller 1953; Orgel 2004", "literature"
    )

    INFLUX_U_RESIDUES: Source = Source(120.0, "residues/s", "assumed", "assumed")
    INFLUX_T_RESIDUES: Source = Source(45.0, "residues/s", "assumed", "assumed")
    INFLUX_C_RESIDUES: Source = Source(80.0, "residues/s", "assumed", "assumed")
    INFLUX_A_RESIDUES: Source = Source(55.0, "residues/s", "assumed", "assumed")
    INFLUX_G_RESIDUES: Source = Source(60.0, "residues/s", "assumed", "assumed")
    INFLUX_X_RESIDUES: Source = Source(0.0, "residues/s", "assumed", "assumed")

    INCORPORATION_BIAS: Dict[str, Source] = {
        'T': Source(1.35, "-", "imposed (v4.x heuristic)", "imposed"),
        'U': Source(0.85, "-", "imposed (v4.x heuristic)", "imposed"),
        'C': Source(0.80, "-", "imposed (v4.x heuristic)", "imposed"),
        'A': Source(0.90, "-", "imposed (v4.x heuristic)", "imposed"),
        'G': Source(0.88, "-", "imposed (v4.x heuristic)", "imposed"),
        'X': Source(1.00, "-", "imposed (xanthine)", "imposed"),
    }

    POLYMER_LENGTH: Source = Source(100, "residues", "assumed", "assumed")

    # ================================================================
    # Rate calculations
    # ================================================================

    @staticmethod
    def _arrhenius_ref(
        T_C: float, T_ref_C: float, k_ref: float, Ea_kcal: float
    ) -> float:
        T_K = T_C + 273.15
        T_ref_K = T_ref_C + 273.15
        R_kcal = PhysicalConstants.R_kcal
        return k_ref * np.exp(
            -Ea_kcal / R_kcal * (1.0 / T_K - 1.0 / T_ref_K)
        )

    @classmethod
    def get_pH_factor(cls, pH_in_situ: float, pKa: float) -> float:
        protonated = 1.0 / (1.0 + 10.0 ** (pH_in_situ - pKa))
        deprotonated = 1.0 - protonated
        slope_term = 1.0 + cls.PH_FACTOR_PH_SLOPE.value * (pH_in_situ - 7.0)
        slope_term = max(0.1, slope_term)
        factor = (
            1.0
            + cls.PH_FACTOR_DEPROTONATION_GAIN.value * deprotonated * slope_term
        )
        return float(
            np.clip(
                factor,
                cls.PH_FACTOR_MIN.value,
                cls.PH_FACTOR_MAX.value,
            )
        )

    @classmethod
    def get_hydrolysis_rates(
        cls, T_C: float, pH_in_situ: float, polymer_length: int = 100
    ) -> Dict[str, float]:
        T_ref = cls.T_HYDRO_REF.value

        k_U = cls._arrhenius_ref(T_C, T_ref, cls.K_U_REF.value, cls.Ea_U.value)
        k_T = cls._arrhenius_ref(T_C, T_ref, cls.K_T_REF.value, cls.Ea_T.value)
        k_A = cls._arrhenius_ref(T_C, T_ref, cls.K_A_REF.value, cls.Ea_A.value)
        k_G = cls._arrhenius_ref(T_C, T_ref, cls.K_G_REF.value, cls.Ea_G.value)

        pH_U = cls.get_pH_factor(pH_in_situ, cls.pKa_U.value)
        pH_T = cls.get_pH_factor(pH_in_situ, cls.pKa_T.value)
        pH_A = cls.get_pH_factor(pH_in_situ, cls.pKa_A.value)
        pH_G = cls.get_pH_factor(pH_in_situ, cls.pKa_G.value)

        OH = WaterChemistry.OH_from_pH(pH_in_situ, T_C)
        base_cat = float(
            np.clip(
                1.0 + cls.BASE_CAT_OH_GAIN.value * OH,
                1.0,
                cls.BASE_CAT_MAX.value,
            )
        )
        acid_factor = float(
            np.clip(
                1.0 + max(0.0, 6.0 - pH_in_situ) * cls.ACID_FACTOR_SLOPE.value,
                1.0,
                cls.ACID_FACTOR_MAX.value,
            )
        )

        if polymer_length >= 100:
            length_factor = max(0.85, 1.0 - 0.001 * (polymer_length - 100))
        else:
            length_factor = 1.0

        return {
            'U': k_U * pH_U * acid_factor * base_cat
                 * cls.BASE_CATALYSIS_FACTOR_U.value * length_factor,
            'T': k_T * pH_T * acid_factor * base_cat
                 * cls.BASE_CATALYSIS_FACTOR_T.value * length_factor,
            'A': k_A * pH_A * acid_factor * base_cat
                 * cls.BASE_CATALYSIS_FACTOR_A.value * length_factor,
            'G': k_G * pH_G * acid_factor * base_cat
                 * cls.BASE_CATALYSIS_FACTOR_G.value * length_factor,
        }

    @classmethod
    def get_deamination_rate_C(
        cls, T_C, pH_in_situ, strand_state='dsDNA'
    ) -> float:
        k = cls._arrhenius_ref(
            T_C, cls.T_C_DEAM_REF.value, cls.K_C_DEAM_REF.value,
            cls.Ea_C_deam.value,
        )
        if strand_state == 'ssDNA':
            k *= cls.SSDNA_DEAM_ENHANCEMENT_C.value
        OH = WaterChemistry.OH_from_pH(pH_in_situ, T_C)
        OH_ref = WaterChemistry.OH_from_pH(7.0, T_C)
        accel = min(
            cls.C_DEAM_PH_ACCEL_MAX.value, OH / max(OH_ref, 1e-12)
        )
        return k * accel

    @classmethod
    def get_deamination_rate_G(
        cls, T_C, pH_in_situ, strand_state='dsDNA'
    ) -> float:
        k = cls._arrhenius_ref(
            T_C, cls.T_G_DEAM_REF.value, cls.K_G_DEAM_REF.value,
            cls.Ea_G_deam.value,
        )
        if strand_state == 'ssDNA':
            k *= cls.SSDNA_DEAM_ENHANCEMENT_G.value
        OH = WaterChemistry.OH_from_pH(pH_in_situ, T_C)
        OH_ref = WaterChemistry.OH_from_pH(7.0, T_C)
        accel = min(
            cls.C_DEAM_PH_ACCEL_MAX.value, OH / max(OH_ref, 1e-12)
        )
        return k * accel

    @classmethod
    def get_free_monomer_deamination_rate_C(
        cls, T_C: float, pH_in_situ: float
    ) -> float:
        k = cls._arrhenius_ref(
            T_C, cls.T_C_DEAM_REF.value, cls.K_FREE_C_DEAM_REF.value,
            cls.Ea_free_C_deam.value,
        )
        OH = WaterChemistry.OH_from_pH(pH_in_situ, T_C)
        OH_ref = WaterChemistry.OH_from_pH(7.0, T_C)
        accel = min(
            cls.C_DEAM_PH_ACCEL_MAX.value, OH / max(OH_ref, 1e-12)
        )
        return k * accel

    @classmethod
    def get_free_monomer_deamination_rate_G(
        cls, T_C: float, pH_in_situ: float
    ) -> float:
        k = cls._arrhenius_ref(
            T_C, cls.T_G_DEAM_REF.value, cls.K_FREE_G_DEAM_REF.value,
            cls.Ea_free_G_deam.value,
        )
        OH = WaterChemistry.OH_from_pH(pH_in_situ, T_C)
        OH_ref = WaterChemistry.OH_from_pH(7.0, T_C)
        accel = min(
            cls.C_DEAM_PH_ACCEL_MAX.value, OH / max(OH_ref, 1e-12)
        )
        return k * accel

    @classmethod
    def get_uv_damage_rate(
        cls, T_C, photon_flux=1e16, wavelength_nm=260.0
    ) -> Dict[str, float]:
        qy = cls.UV_QUANTUM_YIELD.value
        peak = cls.UV_PEAK_WAVELENGTH_NM.value
        sigma_nm = cls.UV_BANDWIDTH_NM.value
        band_factor = float(
            np.exp(-0.5 * ((wavelength_nm - peak) / sigma_nm) ** 2)
        )
        absolute_rate = (
            max(0.0, photon_flux)
            * cls.UV_ABS_CROSS_SECTION_REF_CM2.value
            * qy
            * band_factor
        )
        return {
            'U': absolute_rate * cls.UV_CROSS_SECTION_U.value,
            'T': absolute_rate * cls.UV_CROSS_SECTION_T.value,
            'C': absolute_rate * cls.UV_CROSS_SECTION_C.value,
            'A': absolute_rate * cls.UV_CROSS_SECTION_A.value,
            'G': absolute_rate * cls.UV_CROSS_SECTION_G.value,
            'X': absolute_rate,
        }

    @classmethod
    def get_polymerization_rate(
        cls, T_C, pH_in_situ, polymer_length=100, clay_bound=False
    ) -> float:
        k = cls._arrhenius_ref(
            T_C, cls.T_POLY_REF.value, cls.K_POLY_REF.value,
            cls.Ea_poly.value,
        )
        pH_opt = 7.5
        f_pH = np.exp(-0.5 * ((pH_in_situ - pH_opt) / 2.0) ** 2)
        f_phos = cls.get_pH_factor(pH_in_situ, cls.pKa_phosphate.value)
        f_clay = cls.K_POLY_CLAY_ENHANCEMENT.value if clay_bound else 1.0
        return k * f_pH * f_phos * f_clay

    @classmethod
    def get_stability_ratio(
        cls, T_C, pH_in_situ, polymer_length=100
    ) -> float:
        rates = cls.get_hydrolysis_rates(T_C, pH_in_situ, polymer_length)
        if rates['T'] <= 0:
            return float('inf')
        return rates['U'] / rates['T']

    @classmethod
    def get_uv_cross_section_ratio_U_over_T(cls) -> float:
        return cls.UV_CROSS_SECTION_U.value / cls.UV_CROSS_SECTION_T.value

    @classmethod
    def get_input_thymine_fraction(cls) -> float:
        init_sum = 0.0
        for b in 'UTCAG':
            init_sum += cls.__dict__[f'INITIAL_{b}_RESIDUES'].value
        return cls.INITIAL_T_RESIDUES.value / init_sum

    @classmethod
    def provenance(cls) -> Dict[str, Source]:
        d = {}
        for name in dir(cls):
            if name.startswith('_'):
                continue
            attr = getattr(cls, name)
            if isinstance(attr, Source):
                d[name] = attr
        for b, s in cls.INCORPORATION_BIAS.items():
            d[f"INCORPORATION_BIAS_{b}"] = s
        return d


# ====================================================================
# SECTION 6: STOCHASTIC SIMULATION ENGINE
# ====================================================================

class VSSUFEngine:
    """
    Stochastic engine v6.1.1 (physics identical to v6.1.0).

    KINETIC MONTE CARLO:
        For a destruction channel with rate k on N molecules over dt:
            n_destroyed ~ Binomial(N, 1 - exp(-k * dt))

    COHORT (per-polymer survival):
        Reservoir-sampled subsample of polymerization events. Destroyed
        polymers are marked dead. Deamination does NOT kill a polymer.

    MASS BALANCE:
        mass_consumed[b] == dna_residues[b]
                          + mass_hydrolyzed[b]
                          + mass_uv_destroyed[b]

    TWO ENVIRONMENTAL CHANNELS:
        get_chemical_selection_index() -> total environmental selection
        get_dna_TU_over_free_TU()      -> direct environmental discrimination
        get_free_TU_ratio()             -> free-pool T/(U+T) diagnostic
        CSI / dna_TU_over_free_TU      -> indirect environmental drift

    VENT VS SURFACE:
        Default is vent (photon_flux = 0). UV machinery is retained for
        surface-scenario comparison but is inactive by default.

    ABLATION FLAGS:
        equalize_ut_pools    -> U and T initial pools and influx equalized
        disable_compartments -> clay + lipid compartmentalization disabled
    """

    COHORT_SAMPLE_SIZE = 10_000

    def __init__(
        self,
        temperature_C: float = 90.0,
        pH_in_situ: float = 11.0,
        seed: int = 42,
        max_time_hours: float = 240.0,
        verbose: bool = False,
        photon_flux: float = 0.0,
        wavelength_nm: float = 260.0,
        polymer_length: int = 100,
        lipid_conc: float = 0.05,
        lipid_chain_length: int = 18,
        ionic_strength_M: float = 0.1,
        Ca_mM: float = 1.0,
        record_interval_seconds: float = 300.0,
        dt: float = 60.0,
        strand_state: str = 'dsDNA',
        use_incorporation_bias: bool = False,
        use_fidelity_advantage: bool = False,
        equalize_ut_pools: bool = False,
        disable_compartments: bool = False,
    ):
        self.seed = seed
        np.random.seed(seed)
        random.seed(seed)

        self.temperature_C = temperature_C
        self.pH_in_situ = pH_in_situ
        self.T_kelvin = temperature_C + 273.15
        self.verbose = verbose
        self.max_time_seconds = max_time_hours * 3600.0
        self.photon_flux = photon_flux
        self.wavelength_nm = wavelength_nm
        self.polymer_length = max(2, polymer_length)
        self.record_interval_seconds = record_interval_seconds
        self.dt = dt
        self.strand_state = strand_state
        self.use_incorporation_bias = use_incorporation_bias
        self.use_fidelity_advantage = use_fidelity_advantage
        self.equalize_ut_pools = equalize_ut_pools
        self.disable_compartments = disable_compartments

        self._flags_seen: Set[str] = set()
        self._verified_flags = False

        self.influx_residues = {
            'U': EmpiricalRates.INFLUX_U_RESIDUES.value,
            'T': EmpiricalRates.INFLUX_T_RESIDUES.value,
            'C': EmpiricalRates.INFLUX_C_RESIDUES.value,
            'A': EmpiricalRates.INFLUX_A_RESIDUES.value,
            'G': EmpiricalRates.INFLUX_G_RESIDUES.value,
            'X': EmpiricalRates.INFLUX_X_RESIDUES.value,
        }

        if self.equalize_ut_pools:
            avg_init = 0.5 * (
                EmpiricalRates.INITIAL_U_RESIDUES.value
                + EmpiricalRates.INITIAL_T_RESIDUES.value
            )
            avg_influx = 0.5 * (
                EmpiricalRates.INFLUX_U_RESIDUES.value
                + EmpiricalRates.INFLUX_T_RESIDUES.value
            )
            self._initial_pools_override = {'U': avg_init, 'T': avg_init}
            self.influx_residues['U'] = avg_influx
            self.influx_residues['T'] = avg_influx
        else:
            self._initial_pools_override = None

        self.k_hydro = EmpiricalRates.get_hydrolysis_rates(
            temperature_C, pH_in_situ, self.polymer_length
        )

        if self.disable_compartments:
            self.clay_protection = 1.0
        else:
            self.clay_protection = EmpiricalRates.CLAY_PROTECTION_FACTOR.value

        self.k_hydro_clay = {
            b: k / self.clay_protection for b, k in self.k_hydro.items()
        }

        self.lipid_membrane = LipidMembrane(
            temperature_C=temperature_C,
            pH_in_situ=pH_in_situ,
            lipid_conc=lipid_conc,
            lipid_chain_length=lipid_chain_length,
            ionic_strength_M=ionic_strength_M,
            Ca_mM=Ca_mM,
        )
        if self.disable_compartments:
            self.lipid_protection_factor = 1.0
            self.lipid_membrane.protection_active = False
        else:
            self.lipid_protection_factor = (
                self.lipid_membrane.get_protection_factor()
            )
        self.k_hydro_lipid = {
            b: k / self.lipid_protection_factor
            for b, k in self.k_hydro.items()
        }

        self.k_deam_C = EmpiricalRates.get_deamination_rate_C(
            temperature_C, pH_in_situ, strand_state
        )
        self.k_deam_G = EmpiricalRates.get_deamination_rate_G(
            temperature_C, pH_in_situ, strand_state
        )
        self.k_free_deam_C = (
            EmpiricalRates.get_free_monomer_deamination_rate_C(
                temperature_C, pH_in_situ
            )
        )
        self.k_free_deam_G = (
            EmpiricalRates.get_free_monomer_deamination_rate_G(
                temperature_C, pH_in_situ
            )
        )

        self.k_uv = EmpiricalRates.get_uv_damage_rate(
            temperature_C, photon_flux, wavelength_nm
        )

        self.k_poly = EmpiricalRates.get_polymerization_rate(
            temperature_C, pH_in_situ, self.polymer_length,
            clay_bound=False,
        )

        self.contour_length = (
            self.polymer_length
            * PolymerPhysics.CONTOUR_PER_BASE_DNA.value
        )
        self.persistence_length = PolymerPhysics.PERSISTENCE_DSDNA.value

        self.reset()

        if self.verbose:
            print(
                f"VSSUF v6.1.1: T={temperature_C}C, "
                f"pH_in_situ={pH_in_situ:.2f}, "
                f"photon_flux={photon_flux:.2e}"
            )
            print(f"  strand_state = {strand_state}")
            print(f"  polymer_length = {self.polymer_length} residues")
            print(f"  k_poly = {self.k_poly:.4e} s^-1")
            print(f"  k_deam_C = {self.k_deam_C:.4e} s^-1")
            print(f"  k_free_deam_C = {self.k_free_deam_C:.4e} s^-1")
            print(
                "  k_hydro (s^-1): "
                + ", ".join(
                    f"{b}={self.k_hydro[b]:.3e}" for b in 'UTAG'
                )
            )
            print(
                "  k_uv (s^-1): "
                + ", ".join(f"{b}={self.k_uv[b]:.3e}" for b in 'UTCAG')
            )
            print(
                "  Stability ratio k_U/k_T: "
                f"{EmpiricalRates.get_stability_ratio(temperature_C, pH_in_situ):.2f}"
                "  [extrapolated; see Ea_T sweep]"
            )

    def verify_flags(self) -> None:
        expected = set()
        if self.equalize_ut_pools:
            expected.add('equalize_ut_pools')
        if self.disable_compartments:
            expected.add('disable_compartments')
        if self.use_incorporation_bias:
            expected.add('use_incorporation_bias')
        if self.use_fidelity_advantage:
            expected.add('use_fidelity_advantage')
        missing = expected - self._flags_seen
        if missing:
            raise RuntimeError(
                f"v6.1.1 flag self-check FAILED: flags declared but "
                f"never applied: {sorted(missing)}."
            )
        self._verified_flags = True

    def _mark_flag(self, name: str) -> None:
        self._flags_seen.add(name)

    def reset(self):
        self.species = {
            **{f'{b}_residues': 0.0 for b in 'UTCAGX'},
            **{f'dsDNA_{b}': 0 for b in 'UTCAGX'},
            **{f'dsDNA_{b}_clay': 0 for b in 'UTCAGX'},
            **{f'dsDNA_{b}_lipid': 0 for b in 'UTCAGX'},
        }
        for b in 'UTCAG':
            if (
                self._initial_pools_override is not None
                and b in self._initial_pools_override
            ):
                self.species[f'{b}_residues'] = float(
                    self._initial_pools_override[b]
                )
            else:
                self.species[f'{b}_residues'] = EmpiricalRates.__dict__[
                    f'INITIAL_{b}_RESIDUES'
                ].value
        self.species['X_residues'] = 0.0

        self.cohort_birth: List[float] = []
        self.cohort_alive: List[bool] = []
        self._alive_cohort_idx: Set[int] = set()
        self.polymer_count: int = 0

        self.history = {
            'time': [],
            **{f'dsDNA_{b}': [] for b in 'UTCAGX'},
            **{f'dsDNA_{b}_clay': [] for b in 'UTCAGX'},
            **{f'dsDNA_{b}_lipid': [] for b in 'UTCAGX'},
            'total_dna': [],
            'peak_total_dna': [],
            'time_to_peak_h': [],
            'total_residues': [],
            'fractions': {b: [] for b in 'UTCAGX'},
            'enrichment_T_over_U': [],
            'chemical_selection_index': [],
            'dna_TU_over_free_TU': [],
            'free_TU_ratio': [],
            'dna_TU_fraction': [],
            'cohort_survival': [],
            'residue_survival': [],
            'residue_retention_ratio': [],
            'polymerization_events': [],
            'hydrolysis_events': [],
            'uv_damage_events': [],
            'deamination_C_events': [],
            'deamination_G_events': [],
            'csi': [],
            'real': [],
            'drift': [],
            'free_T_frac': [],
            'pool_influx_T_frac': [],
        }

        self.time = 0.0
        self.step_count = 0
        self.polymerization_events = 0
        self.hydrolysis_events = 0
        self.uv_damage_events = 0
        self.deamination_C_events = 0
        self.deamination_G_events = 0
        self.clay_adsorption_events = 0
        self.clay_desorption_events = 0
        self.lipid_entry_events = 0
        self.lipid_exit_events = 0

        self.mass_influx = {b: 0.0 for b in 'UTCAGX'}
        self.mass_consumed = {b: 0.0 for b in 'UTCAGX'}
        self.mass_hydrolyzed = {b: 0.0 for b in 'UTCAGX'}
        self.mass_uv_destroyed = {b: 0.0 for b in 'UTCAGX'}
        self.mass_deam_C_to_U = 0.0
        self.mass_deam_G_to_X = 0.0
        self.mass_free_C_to_U = 0.0
        self.mass_free_G_to_X = 0.0

        if self.equalize_ut_pools:
            self._mark_flag('equalize_ut_pools')
        if self.disable_compartments:
            self._mark_flag('disable_compartments')

    # ----- Utility -----
    def residues_per_polymerization(self) -> int:
        return self.polymer_length

    def total_dna(self) -> int:
        total = 0
        for b in 'UTCAGX':
            total += self.species[f'dsDNA_{b}']
            total += self.species[f'dsDNA_{b}_clay']
            total += self.species[f'dsDNA_{b}_lipid']
        return int(total)

    def total_dna_residues(self) -> int:
        return self.total_dna() * self.polymer_length

    def total_consumed_residues(self) -> float:
        return float(sum(self.mass_consumed[b] for b in 'UTCAGX'))

    def total_destroyed_residues(self) -> float:
        return float(
            sum(
                self.mass_hydrolyzed[b] + self.mass_uv_destroyed[b]
                for b in 'UTCAGX'
            )
        )

    # ----- Metrics -----
    def get_dna_TU_fraction(self) -> float:
        t_count = (
            self.species['dsDNA_T']
            + self.species['dsDNA_T_clay']
            + self.species['dsDNA_T_lipid']
        )
        u_count = (
            self.species['dsDNA_U']
            + self.species['dsDNA_U_clay']
            + self.species['dsDNA_U_lipid']
        )
        denom = t_count + u_count
        if denom <= 0:
            return 0.0
        return float(t_count / denom)

    def get_free_TU_ratio(self) -> float:
        free_T = self.species['T_residues']
        free_U = self.species['U_residues']
        denom = free_T + free_U
        if denom <= 0:
            return 0.0
        return float(free_T / denom)

    def get_dna_TU_over_free_TU(self) -> float:
        dna_frac = self.get_dna_TU_fraction()
        free_frac = self.get_free_TU_ratio()
        if free_frac <= 0:
            return 0.0
        return float(dna_frac / free_frac)

    def get_pool_drift_factor(self) -> float:
        real = self.get_dna_TU_over_free_TU()
        if real <= 0:
            return 0.0
        return float(self.get_chemical_selection_index() / real)

    # ----- Mass balance -----
    def check_mass_balance(self) -> Dict[str, Dict[str, float]]:
        out = {}
        deam_out = {
            'C': self.mass_deam_C_to_U + self.mass_free_C_to_U,
            'G': self.mass_deam_G_to_X + self.mass_free_G_to_X,
        }
        deam_in = {
            'U': self.mass_deam_C_to_U + self.mass_free_C_to_U,
            'X': self.mass_deam_G_to_X + self.mass_free_G_to_X,
        }
        for b in 'UTCAGX':
            init_src = EmpiricalRates.__dict__.get(f'INITIAL_{b}_RESIDUES')
            if (
                self._initial_pools_override is not None
                and b in self._initial_pools_override
            ):
                initial = self._initial_pools_override[b]
            else:
                initial = init_src.value if init_src is not None else 0.0
            free_monomer = self.species[f'{b}_residues']
            dna = (
                self.species[f'dsDNA_{b}']
                + self.species[f'dsDNA_{b}_clay']
                + self.species[f'dsDNA_{b}_lipid']
            )
            dna_residues = dna * self.polymer_length
            lhs = initial + self.mass_influx[b]
            rhs = (
                free_monomer
                + dna_residues
                + self.mass_hydrolyzed[b]
                + self.mass_uv_destroyed[b]
                + deam_out.get(b, 0.0)
                - deam_in.get(b, 0.0)
            )
            expected_consumed = (
                dna_residues
                + self.mass_hydrolyzed[b]
                + self.mass_uv_destroyed[b]
            )
            consumed_residual = self.mass_consumed[b] - expected_consumed
            consumed_rel = consumed_residual / max(
                abs(self.mass_consumed[b]), 1.0
            )
            out[b] = {
                'lhs': lhs,
                'rhs': rhs,
                'residual': lhs - rhs,
                'rel_error': (lhs - rhs) / max(abs(lhs), 1.0),
                'consumed_residual': consumed_residual,
                'consumed_rel_err': consumed_rel,
            }
        return out

    def max_abs_rel_error(self) -> float:
        mb = self.check_mass_balance()
        return max(abs(v['rel_error']) for v in mb.values())

    def max_abs_consumed_residual(self) -> float:
        mb = self.check_mass_balance()
        return max(abs(v['consumed_residual']) for v in mb.values())

    def max_abs_consumed_rel_error(self) -> float:
        mb = self.check_mass_balance()
        return max(abs(v['consumed_rel_err']) for v in mb.values())

    def deamination_transfer_residual(self) -> float:
        mb = self.check_mass_balance()
        return float(sum(v['consumed_residual'] for v in mb.values()))

    # ----- Influx -----
    def _get_vent_influx(self) -> Dict[str, float]:
        fluct = 1.0 + 0.3 * (2 * np.random.random() - 1)
        pulse = 0.8 + 0.2 * np.sin(2 * np.pi * self.time / 3600.0)
        return {b: self.influx_residues[b] * fluct * pulse for b in 'UTCAGX'}

    def _pool_influx_T_frac(self) -> float:
        return self.influx_residues['T'] / max(
            1e-12, self.influx_residues['T'] + self.influx_residues['U']
        )

    def _fidelity_advantage(self) -> float:
        if not self.use_fidelity_advantage:
            return 1.0
        self._mark_flag('use_fidelity_advantage')
        total = self.total_dna()
        if total < 100:
            return 1.0
        t_count = (
            self.species['dsDNA_T']
            + self.species['dsDNA_T_clay']
            + self.species['dsDNA_T_lipid']
        )
        t_frac = t_count / total
        return 1.0 + EmpiricalRates.FIDELITY_ADVANTAGE_GAIN.value * t_frac

    def _incorporation_bias(self, base: str) -> float:
        if not self.use_incorporation_bias:
            return 1.0
        self._mark_flag('use_incorporation_bias')
        return EmpiricalRates.INCORPORATION_BIAS[base].value

    # ----- Cohort -----
    def _record_polymer(self) -> None:
        self.polymer_count += 1
        size = self.COHORT_SAMPLE_SIZE
        if len(self.cohort_birth) < size:
            self.cohort_birth.append(self.time)
            self.cohort_alive.append(True)
            self._alive_cohort_idx.add(len(self.cohort_birth) - 1)
        else:
            if random.random() < size / self.polymer_count:
                idx = random.randrange(size)
                self.cohort_birth[idx] = self.time
                self.cohort_alive[idx] = True
                self._alive_cohort_idx.add(idx)

    def _kill_cohort_members(self, n: int, n_total: int) -> None:
        if n <= 0 or n_total <= 0 or not self._alive_cohort_idx:
            return
        p_kill = min(1.0, n / n_total)
        victims = [
            idx for idx in self._alive_cohort_idx
            if random.random() < p_kill
        ]
        for idx in victims:
            self.cohort_alive[idx] = False
            self._alive_cohort_idx.discard(idx)

    def _stochastic_destroy(self, state: str, base: str, rate: float) -> int:
        key = f'dsDNA_{base}{state}'
        N = int(self.species[key])
        if N <= 0:
            return 0
        p = 1.0 - np.exp(-rate * self.dt)
        if p <= 0.0:
            return 0
        if p >= 1.0:
            return N
        return int(np.random.binomial(N, p))

    # ----- Step -----
    def step(self) -> bool:
        hydro_before = self.hydrolysis_events
        uv_before = self.uv_damage_events

        influx = self._get_vent_influx()
        for b in 'UTCAGX':
            added = influx[b] * self.dt
            self.species[f'{b}_residues'] += added
            self.mass_influx[b] += added

        self.time += self.dt
        self.step_count += 1

        fidelity = self._fidelity_advantage()
        req = self.residues_per_polymerization()
        lam = self.k_poly * self.dt
        n_attempts = int(np.random.poisson(lam)) if lam < 30 else int(lam)

        for _ in range(n_attempts):
            eligible = [
                b for b in 'UTCAG'
                if self.species[f'{b}_residues'] >= req
            ]
            if not eligible:
                break
            weights = np.array([
                max(0.0, self.species[f'{b}_residues'] / req)
                * max(0.0, self._incorporation_bias(b))
                for b in eligible
            ], dtype=float)
            if weights.sum() <= 0:
                break
            weights /= weights.sum()
            b = str(np.random.choice(eligible, p=weights))
            p = EmpiricalRates.FIDELITY_BASE.value * fidelity
            p = min(1.0, max(0.0, p))
            if np.random.random() < p:
                self.species[f'dsDNA_{b}'] += 1
                self.species[f'{b}_residues'] -= req
                self.mass_consumed[b] += req
                self.polymerization_events += 1
                self._record_polymer()

        for b in 'UTAG':
            for state, k_rate in (
                ('', self.k_hydro[b]),
                ('_clay', self.k_hydro_clay[b]),
                ('_lipid', self.k_hydro_lipid[b]),
            ):
                n = self._stochastic_destroy(state, b, k_rate)
                if n > 0:
                    self.species[f'dsDNA_{b}{state}'] -= n
                    self.mass_hydrolyzed[b] += n * self.polymer_length
                    self.hydrolysis_events += n

        for b in 'UTCAGX':
            for state, k_rate in (
                ('', self.k_uv[b]),
                ('_clay', self.k_uv[b] / self.clay_protection),
                ('_lipid', self.k_uv[b] / self.lipid_protection_factor),
            ):
                n = self._stochastic_destroy(state, b, k_rate)
                if n > 0:
                    self.species[f'dsDNA_{b}{state}'] -= n
                    self.mass_uv_destroyed[b] += n * self.polymer_length
                    self.uv_damage_events += n

        n_destroyed = (self.hydrolysis_events - hydro_before) + (
            self.uv_damage_events - uv_before
        )
        if n_destroyed > 0:
            n_total_after = self.total_dna()
            n_total_before = n_total_after + n_destroyed
            self._kill_cohort_members(n_destroyed, n_total_before)

        for state in ('', '_clay', '_lipid'):
            key = f'dsDNA_C{state}'
            N = int(self.species[key])
            if N > 0:
                p = 1.0 - np.exp(-self.k_deam_C * self.dt)
                if p > 0.0:
                    n = N if p >= 1.0 else int(np.random.binomial(N, p))
                    if n > 0:
                        transfer = n * self.polymer_length
                        self.species[key] -= n
                        self.species[f'dsDNA_U{state}'] += n
                        self.deamination_C_events += n
                        self.mass_deam_C_to_U += transfer
                        self.mass_consumed['C'] -= transfer
                        self.mass_consumed['U'] += transfer

        floor = int(EmpiricalRates.FREE_MONOMER_FLOOR.value)
        if self.species['C_residues'] > floor:
            p_free_C = 1.0 - np.exp(-self.k_free_deam_C * self.dt)
            N_free_C = int(self.species['C_residues'] - floor)
            if N_free_C > 0 and p_free_C > 0.0:
                lost = (
                    N_free_C if p_free_C >= 1.0
                    else int(np.random.binomial(N_free_C, p_free_C))
                )
                if lost > 0:
                    self.species['C_residues'] -= lost
                    self.species['U_residues'] += lost
                    self.mass_free_C_to_U += lost

        for state in ('', '_clay', '_lipid'):
            key = f'dsDNA_G{state}'
            N = int(self.species[key])
            if N > 0:
                p = 1.0 - np.exp(-self.k_deam_G * self.dt)
                if p > 0.0:
                    n = N if p >= 1.0 else int(np.random.binomial(N, p))
                    if n > 0:
                        transfer = n * self.polymer_length
                        self.species[key] -= n
                        self.species[f'dsDNA_X{state}'] += n
                        self.deamination_G_events += n
                        self.mass_deam_G_to_X += transfer
                        self.mass_consumed['G'] -= transfer
                        self.mass_consumed['X'] += transfer

        if self.species['G_residues'] > floor:
            p_free_G = 1.0 - np.exp(-self.k_free_deam_G * self.dt)
            N_free_G = int(self.species['G_residues'] - floor)
            if N_free_G > 0 and p_free_G > 0.0:
                lost = (
                    N_free_G if p_free_G >= 1.0
                    else int(np.random.binomial(N_free_G, p_free_G))
                )
                if lost > 0:
                    self.species['G_residues'] -= lost
                    self.species['X_residues'] += lost
                    self.mass_free_G_to_X += lost

        if not self.disable_compartments:
            p_ads = 1.0 - np.exp(
                -EmpiricalRates.CLAY_ADSORPTION_RATE.value * self.dt
            )
            p_des = 1.0 - np.exp(
                -EmpiricalRates.CLAY_DESORPTION_RATE.value * self.dt
            )
            for b in 'UTCAGX':
                N_free = int(self.species[f'dsDNA_{b}'])
                if N_free > 0 and p_ads > 0.0:
                    n_ads = (
                        N_free if p_ads >= 1.0
                        else int(np.random.binomial(N_free, p_ads))
                    )
                    if n_ads > 0:
                        self.species[f'dsDNA_{b}'] -= n_ads
                        self.species[f'dsDNA_{b}_clay'] += n_ads
                        self.clay_adsorption_events += n_ads
                N_clay = int(self.species[f'dsDNA_{b}_clay'])
                if N_clay > 0 and p_des > 0.0:
                    n_des = (
                        N_clay if p_des >= 1.0
                        else int(np.random.binomial(N_clay, p_des))
                    )
                    if n_des > 0:
                        self.species[f'dsDNA_{b}_clay'] -= n_des
                        self.species[f'dsDNA_{b}'] += n_des
                        self.clay_desorption_events += n_des

        if not self.disable_compartments:
            if self.lipid_membrane.protection_active:
                p_entry = 1.0 - np.exp(-0.08 * self.dt)
                p_exit = 1.0 - np.exp(-0.002 * self.dt)
                for b in 'UTCAGX':
                    N_free = int(self.species[f'dsDNA_{b}'])
                    if N_free > 5 and p_entry > 0.0:
                        pc = self.lipid_membrane.get_partition_coefficient(b)
                        lam_entry = 0.15 * pc * p_entry * N_free
                        n_entry = (
                            int(np.random.poisson(lam_entry))
                            if lam_entry < 30 else int(lam_entry)
                        )
                        n_entry = min(n_entry, N_free - 5)
                        if n_entry > 0:
                            self.species[f'dsDNA_{b}'] -= n_entry
                            self.species[f'dsDNA_{b}_lipid'] += n_entry
                            self.lipid_entry_events += n_entry
                    N_lip = int(self.species[f'dsDNA_{b}_lipid'])
                    if N_lip > 0 and p_exit > 0.0:
                        n_exit = (
                            N_lip if p_exit >= 1.0
                            else int(np.random.binomial(N_lip, p_exit))
                        )
                        if n_exit > 0:
                            self.species[f'dsDNA_{b}_lipid'] -= n_exit
                            self.species[f'dsDNA_{b}'] += n_exit
                            self.lipid_exit_events += n_exit

        self._record_history()
        return True

    def _record_history(self):
        if (
            len(self.history['time']) == 0
            or (self.time - self.history['time'][-1])
            >= self.record_interval_seconds
        ):
            totals = {}
            for b in 'UTCAGX':
                totals[b] = (
                    self.species[f'dsDNA_{b}']
                    + self.species[f'dsDNA_{b}_clay']
                    + self.species[f'dsDNA_{b}_lipid']
                )
            total = sum(totals.values())
            total_res = total * self.polymer_length

            self.history['time'].append(self.time)
            for b in 'UTCAGX':
                self.history[f'dsDNA_{b}'].append(totals[b])
                self.history[f'dsDNA_{b}_clay'].append(
                    self.species[f'dsDNA_{b}_clay']
                )
                self.history[f'dsDNA_{b}_lipid'].append(
                    self.species[f'dsDNA_{b}_lipid']
                )
            self.history['total_dna'].append(total)
            peak = max(self.history['total_dna'])
            peak_idx = int(np.argmax(self.history['total_dna']))
            self.history['peak_total_dna'].append(peak)
            self.history['time_to_peak_h'].append(
                self.history['time'][peak_idx] / 3600.0
            )
            self.history['total_residues'].append(total_res)
            for b in 'UTCAGX':
                self.history['fractions'][b].append(
                    totals[b] / max(1, total)
                )
            self.history['enrichment_T_over_U'].append(
                totals['T'] / max(1, totals['U'])
                if totals['U'] > 0 else 0.0
            )
            self.history['chemical_selection_index'].append(
                self.get_chemical_selection_index()
            )
            self.history['dna_TU_over_free_TU'].append(
                self.get_dna_TU_over_free_TU()
            )
            self.history['free_TU_ratio'].append(
                self.get_free_TU_ratio()
            )
            self.history['dna_TU_fraction'].append(
                self.get_dna_TU_fraction()
            )
            self.history['cohort_survival'].append(
                self.get_current_cohort_survival()
            )
            self.history['residue_survival'].append(
                self.get_residue_survival()
            )
            self.history['residue_retention_ratio'].append(
                self.get_residue_retention_ratio()
            )
            self.history['polymerization_events'].append(
                self.polymerization_events
            )
            self.history['hydrolysis_events'].append(
                self.hydrolysis_events
            )
            self.history['uv_damage_events'].append(self.uv_damage_events)
            self.history['deamination_C_events'].append(
                self.deamination_C_events
            )
            self.history['deamination_G_events'].append(
                self.deamination_G_events
            )
            self.history['csi'].append(self.get_chemical_selection_index())
            self.history['real'].append(self.get_dna_TU_over_free_TU())
            self.history['drift'].append(self.get_pool_drift_factor())
            self.history['free_T_frac'].append(self.get_free_TU_ratio())
            self.history['pool_influx_T_frac'].append(
                self._pool_influx_T_frac()
            )

    def run(self, max_time: Optional[float] = None) -> Dict:
        if max_time is None:
            max_time = self.max_time_seconds
        steps = int(max_time / self.dt)
        for i in range(steps):
            self.step()
            if i % 5000 == 0 and self.verbose:
                print(f"  Step {i}/{steps} ({100*i/steps:.1f}%)")
        if self.verbose:
            self._print_final_report()
        self.verify_flags()
        return self.history

    def _print_final_report(self):
        print(f"  Done. Steps={self.step_count:,}")
        print(f"  Polymerization events: {self.polymerization_events:,}")
        print(f"  Hydrolysis events: {self.hydrolysis_events:,}")
        print(f"  UV damage events: {self.uv_damage_events:,}")
        print(f"  C deamination events: {self.deamination_C_events:,}")
        print(f"  G deamination events: {self.deamination_G_events:,}")
        print(f"  Total DNA molecules: {self.total_dna():,}")
        print(f"  Total DNA residues: {self.total_dna_residues():,}")
        print(f"  Regime: {self.get_regime()}")
        print(
            f"  CSI: {self.get_chemical_selection_index():.3f}x   "
            f"real: {self.get_dna_TU_over_free_TU():.3f}x   "
            f"drift: {self.get_pool_drift_factor():.3f}x"
        )
        print(
            f"  CSI trajectory: peak={self.peak_csi():.3f}x at "
            f"{self.time_to_peak_csi_h():.2f} h, "
            f"end={self.get_chemical_selection_index():.3f}x, "
            f"shape={self.csi_regime()}"
        )
        print(
            f"  Free T/(U+T): {self.get_free_TU_ratio():.4f}   "
            f"DNA T/(U+T): {self.get_dna_TU_fraction():.4f}"
        )
        print(
            f"  Influx T/(U+T): {self._pool_influx_T_frac():.4f}"
        )
        print(
            f"  Stability ratio k_U/k_T: "
            f"{EmpiricalRates.get_stability_ratio(self.temperature_C, self.pH_in_situ):.2f}"
            "  [extrapolated]"
        )
        print(
            f"  Max |rel err| mass balance: "
            f"{self.max_abs_rel_error():.3e}"
        )
        print(
            f"  Max |consumed identity residual|: "
            f"{self.max_abs_consumed_residual():.3e}"
        )

    # ----- CSI trajectory statistics -----
    def peak_csi(self) -> float:
        if not self.history['csi']:
            return 0.0
        return float(np.max(self.history['csi']))

    def time_to_peak_csi_h(self) -> float:
        if not self.history['csi']:
            return 0.0
        idx = int(np.argmax(self.history['csi']))
        return float(self.history['time'][idx] / 3600.0)

    def csi_regime(self) -> str:
        if len(self.history['csi']) < 3:
            return 'no_signal'
        arr = np.asarray(self.history['csi'], dtype=float)
        peak = float(np.max(arr))
        end = float(arr[-1])
        idx_peak = int(np.argmax(arr))
        frac_peak_pos = idx_peak / max(1, len(arr) - 1)
        if peak <= 0:
            return 'no_signal'
        if (peak - np.min(arr)) / peak < 0.05:
            return 'steady'
        if frac_peak_pos >= 0.8:
            return 'rising'
        if end < 0.5 * peak:
            return 'collapsed'
        if end < 0.95 * peak:
            return 'peaked'
        return 'steady'

    # ----- Other metrics -----
    def get_final_thymine_fraction(self) -> float:
        totals = {}
        for b in 'UTCAGX':
            totals[b] = (
                self.species[f'dsDNA_{b}']
                + self.species[f'dsDNA_{b}_clay']
                + self.species[f'dsDNA_{b}_lipid']
            )
        total = sum(totals.values())
        return totals['T'] / max(1, total)

    def get_chemical_selection_index(self) -> float:
        input_frac = EmpiricalRates.get_input_thymine_fraction()
        final_frac = self.get_final_thymine_fraction()
        if input_frac <= 0:
            return 0.0
        return final_frac / input_frac

    def get_effective_hydrolysis_half_life(self) -> float:
        totals = {}
        for b in 'UTAG':
            totals[b] = (
                self.species[f'dsDNA_{b}']
                + self.species[f'dsDNA_{b}_clay']
                + self.species[f'dsDNA_{b}_lipid']
            )
        total = sum(totals.values())
        if total == 0:
            return float('inf')
        k_avg = sum(totals[b] * self.k_hydro[b] for b in 'UTAG') / total
        return np.log(2) / k_avg if k_avg > 0 else float('inf')

    def get_current_cohort_survival(self) -> float:
        if len(self.cohort_alive) == 0:
            return 1.0
        return float(sum(self.cohort_alive) / len(self.cohort_alive))

    def get_residue_survival(self) -> float:
        consumed = self.total_consumed_residues()
        if consumed <= 0:
            return 1.0
        return float(self.total_dna_residues() / consumed)

    def get_residue_retention_ratio(self) -> float:
        dna_res = float(self.total_dna_residues())
        destroyed = self.total_destroyed_residues()
        denom = dna_res + destroyed
        if denom <= 0:
            return 1.0
        return float(dna_res / denom)

    def measure_true_survival_half_life(self) -> float:
        if not self.history['cohort_survival']:
            return float('inf')
        surv = np.array(self.history['cohort_survival'])
        times = np.array(self.history['time'])
        below = np.where(surv <= 0.5)[0]
        if len(below) == 0:
            return float('inf')
        return float(times[below[0]])

    def measure_residue_half_life(self) -> float:
        if not self.history['residue_survival']:
            return float('inf')
        surv = np.array(self.history['residue_survival'])
        times = np.array(self.history['time'])
        below = np.where(surv <= 0.5)[0]
        if len(below) == 0:
            return float('inf')
        return float(times[below[0]])

    def get_final_survival_fraction(self) -> float:
        if not self.history['cohort_survival']:
            return 1.0
        return float(self.history['cohort_survival'][-1])

    def peak_total_dna(self) -> int:
        if not self.history['total_dna']:
            return self.total_dna()
        return int(max(self.history['total_dna']))

    def time_to_peak_seconds(self) -> float:
        if not self.history['total_dna']:
            return 0.0
        idx = int(np.argmax(self.history['total_dna']))
        return float(self.history['time'][idx])

    def net_dna_change(self) -> float:
        if len(self.history['total_dna']) < 2:
            return 0.0
        return float(self.history['total_dna'][-1] - self.history['total_dna'][0])

    def is_population_growth(self) -> bool:
        return self.net_dna_change() > 0.0

    def get_regime(self) -> str:
        if self.total_dna() == 0 and (
            self.hydrolysis_events + self.uv_damage_events > 0
        ):
            return 'sterilized'
        peak = self.peak_total_dna()
        final = self.total_dna()
        if peak > 0 and final < 0.5 * peak and self.time_to_peak_seconds() < self.time:
            return 'post_peak_collapse'
        if self.net_dna_change() > 0:
            return 'net_population_growth'
        return 'net_population_destruction'

    def get_uv_resistance_ratio(self) -> float:
        if self.k_uv['T'] > 0:
            return self.k_uv['U'] / self.k_uv['T']
        return 0.0

    def get_nucleotide_composition(self) -> Dict[str, int]:
        comp = {}
        for b in 'UTCAGX':
            comp[b] = int(
                self.species[f'dsDNA_{b}']
                + self.species[f'dsDNA_{b}_clay']
                + self.species[f'dsDNA_{b}_lipid']
            )
        return comp

    def get_hydrolysis_half_life(self, base: str) -> float:
        k = self.k_hydro.get(base, 0.0)
        if k <= 0:
            return float('inf')
        return float(np.log(2) / k)

    def coefficient_of_variation(self, values) -> float:
        arr = np.asarray(values, dtype=float)
        mean = float(np.mean(arr)) if arr.size else 0.0
        if arr.size < 2 or mean == 0.0:
            return 0.0
        return float(np.std(arr, ddof=1) / abs(mean))

    def save_history_to_csv(
        self, output_dir: str = "output_data_v610"
    ) -> str:
        os.makedirs(output_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        data = {'time_hours': np.array(self.history['time']) / 3600.0}
        for b in 'UTCAGX':
            data[f'dsDNA_{b}'] = self.history[f'dsDNA_{b}']
            data[f'fraction_{b}'] = self.history['fractions'][b]
        data['total_dna'] = self.history['total_dna']
        data['peak_total_dna'] = self.history['peak_total_dna']
        data['time_to_peak_h'] = self.history['time_to_peak_h']
        data['total_residues'] = self.history['total_residues']
        data['enrichment_T_over_U'] = self.history['enrichment_T_over_U']
        data['chemical_selection_index'] = self.history['chemical_selection_index']
        data['dna_TU_over_free_TU'] = self.history['dna_TU_over_free_TU']
        data['free_TU_ratio'] = self.history['free_TU_ratio']
        data['dna_TU_fraction'] = self.history['dna_TU_fraction']
        data['cohort_survival'] = self.history['cohort_survival']
        data['residue_survival'] = self.history['residue_survival']
        data['residue_retention_ratio'] = self.history['residue_retention_ratio']
        data['csi'] = self.history['csi']
        data['real'] = self.history['real']
        data['drift'] = self.history['drift']
        data['free_T_frac'] = self.history['free_T_frac']
        data['pool_influx_T_frac'] = self.history['pool_influx_T_frac']
        df = pd.DataFrame(data)
        csv_file = f"{output_dir}/simulation_history_v610_{timestamp}.csv"
        df.to_csv(csv_file, index=False)
        print(f"CSV saved: {csv_file}")
        return csv_file


# ====================================================================
# SECTION 7: TWO-DIMENSIONAL SENSITIVITY ANALYZER
# ====================================================================

class TwoDSensitivityAnalyzer:
    """
    Runs the model over a (T, pH) grid under the variants:

      CHEMICAL_SELECTION  : no imposed bias; only literature rates.
      HEURISTIC           : v4.x incorporation bias + fidelity adv.
      ZERO_UV_CONTROL     : chemical selection with photon_flux = 0.
      EQUAL_POOLS_CONTROL : chemical selection with equal U/T pools.
    """

    ENRICHMENT_TIE_TOL = 0.01

    def __init__(self, config: dict):
        self.config = config
        self.surface_chemical = None
        self.surface_heuristic = None
        self.surface_zero_uv = None
        self.surface_equal_pools = None
        self.optimization_results = None
        self.flux_sweep_results = None
        self.deam_scale_sweep_results = None
        self.clay_rate_sweep_results = None
        self.ea_T_sweep_results = None
        self.time_horizon_sweep_results = None
        self.confirmatory_optimum_results = None

    def _run_grid(
        self,
        use_incorporation_bias: bool,
        use_fidelity_advantage: bool,
        photon_flux_override: Optional[float] = None,
        equalize_ut_pools: bool = False,
        disable_compartments: bool = False,
        label: str = "",
    ) -> Dict:
        T_start = self.config['temperature_min']
        T_stop = (
            self.config['temperature_max']
            + self.config['temperature_step'] * 0.1
        )
        T_step = self.config['temperature_step']
        T = np.arange(T_start, T_stop, T_step)

        pH_start = self.config['pH_min']
        pH_stop = self.config['pH_max'] + self.config['pH_step'] * 0.1
        pH_step = self.config['pH_step']
        pH = np.arange(pH_start, pH_stop, pH_step)

        n_rep = max(1, int(self.config.get('n_replicates', 1)))

        flux = (
            photon_flux_override
            if photon_flux_override is not None
            else self.config['photon_flux']
        )

        surface = {
            'T': T,
            'pH': pH,
            'photon_flux': flux,
            'equalize_ut_pools': equalize_ut_pools,
            'disable_compartments': disable_compartments,
            'enrichment': np.zeros((len(T), len(pH))),
            'enrichment_std': np.zeros((len(T), len(pH))),
            'dna_TU_over_free_TU': np.zeros((len(T), len(pH))),
            'dna_TU_over_free_TU_std': np.zeros((len(T), len(pH))),
            'pool_drift_factor': np.zeros((len(T), len(pH))),
            'free_TU_ratio': np.zeros((len(T), len(pH))),
            'dna_TU_fraction': np.zeros((len(T), len(pH))),
            'fraction_T': np.zeros((len(T), len(pH))),
            'total_dna': np.zeros((len(T), len(pH))),
            'peak_total_dna': np.zeros((len(T), len(pH))),
            'time_to_peak_h': np.zeros((len(T), len(pH))),
            'final_survival': np.zeros((len(T), len(pH))),
            'residue_survival': np.zeros((len(T), len(pH))),
            'residue_retention_ratio': np.zeros((len(T), len(pH))),
            'is_population_growth': np.zeros((len(T), len(pH)), dtype=bool),
            'half_life_effective': np.zeros((len(T), len(pH))),
            'half_life_true': np.zeros((len(T), len(pH))),
            'residue_half_life': np.zeros((len(T), len(pH))),
            'stability_ratio': np.zeros((len(T), len(pH))),
            'uv_ratio': np.zeros((len(T), len(pH))),
            'poly_rate': np.zeros((len(T), len(pH))),
            'mass_balance_rel_err': np.zeros((len(T), len(pH))),
            'consumed_residual': np.zeros((len(T), len(pH))),
            'consumed_rel_err_mean': np.zeros((len(T), len(pH))),
            'consumed_rel_err_max': np.zeros((len(T), len(pH))),
            'deamination_transfer_residual': np.zeros((len(T), len(pH))),
            'deamination_C_events': np.zeros((len(T), len(pH))),
            'deamination_G_events': np.zeros((len(T), len(pH))),
            'hydrolysis_events': np.zeros((len(T), len(pH))),
            'uv_damage_events': np.zeros((len(T), len(pH))),
            'regime_code': np.zeros((len(T), len(pH)), dtype=int),
            'peak_csi': np.zeros((len(T), len(pH))),
            'time_to_peak_csi_h': np.zeros((len(T), len(pH))),
            'csi_shape_code': np.zeros((len(T), len(pH)), dtype=int),
            'composition': {b: np.zeros((len(T), len(pH))) for b in 'UTCAGX'},
        }

        regime_map = {
            'sterilized': 0,
            'net_population_destruction': 1,
            'net_population_growth': 2,
            'post_peak_collapse': 3,
        }
        csi_shape_map = {
            'no_signal': 0, 'steady': 1, 'rising': 2, 'peaked': 3, 'collapsed': 4,
        }

        print(
            f"\n[{label}] Running grid {len(T)}x{len(pH)} n_rep={n_rep} "
            f"flux={flux:.2e} equalize_ut_pools={equalize_ut_pools} "
            f"disable_compartments={disable_compartments}"
        )
        for i, T_val in enumerate(T):
            for j, pH_val in enumerate(pH):
                reps = []
                for rep in range(n_rep):
                    seed = 42 + i * 100 + j + rep * 10007
                    eng = VSSUFEngine(
                        temperature_C=float(T_val),
                        pH_in_situ=float(pH_val),
                        seed=seed,
                        max_time_hours=self.config['simulation_hours'],
                        verbose=False,
                        photon_flux=flux,
                        wavelength_nm=self.config['wavelength_nm'],
                        polymer_length=self.config['polymer_length'],
                        lipid_conc=self.config['lipid_conc'],
                        record_interval_seconds=self.config['record_interval_seconds'],
                        dt=self.config['dt'],
                        strand_state=self.config['strand_state'],
                        use_incorporation_bias=use_incorporation_bias,
                        use_fidelity_advantage=use_fidelity_advantage,
                        equalize_ut_pools=equalize_ut_pools,
                        disable_compartments=disable_compartments,
                    )
                    eng.run()
                    comp = eng.get_nucleotide_composition()
                    rep_dict = {
                        'enrichment': eng.get_chemical_selection_index(),
                        'dna_TU_over_free_TU': eng.get_dna_TU_over_free_TU(),
                        'pool_drift_factor': eng.get_pool_drift_factor(),
                        'free_TU_ratio': eng.get_free_TU_ratio(),
                        'dna_TU_fraction': eng.get_dna_TU_fraction(),
                        'fraction_T': eng.get_final_thymine_fraction(),
                        'total_dna': eng.total_dna(),
                        'peak_total_dna': eng.peak_total_dna(),
                        'time_to_peak_h': eng.time_to_peak_seconds() / 3600.0,
                        'final_survival': eng.get_final_survival_fraction(),
                        'residue_survival': eng.get_residue_survival(),
                        'residue_retention_ratio': eng.get_residue_retention_ratio(),
                        'is_population_growth': eng.is_population_growth(),
                        'half_life_eff': eng.get_effective_hydrolysis_half_life(),
                        'half_life_true': eng.measure_true_survival_half_life(),
                        'residue_half_life': eng.measure_residue_half_life(),
                        'stability_ratio': EmpiricalRates.get_stability_ratio(
                            float(T_val), float(pH_val),
                            self.config['polymer_length'],
                        ),
                        'uv_ratio': EmpiricalRates.get_uv_cross_section_ratio_U_over_T(),
                        'poly_rate': EmpiricalRates.get_polymerization_rate(
                            float(T_val), float(pH_val),
                            self.config['polymer_length'],
                        ),
                        'mass_balance_rel_err': eng.max_abs_rel_error(),
                        'consumed_residual': eng.max_abs_consumed_residual(),
                        'consumed_rel_err': eng.max_abs_consumed_rel_error(),
                        'deam_transfer_residual': eng.deamination_transfer_residual(),
                        'deam_C': eng.deamination_C_events,
                        'deam_G': eng.deamination_G_events,
                        'hydro': eng.hydrolysis_events,
                        'uv': eng.uv_damage_events,
                        'regime': eng.get_regime(),
                        'peak_csi': eng.peak_csi(),
                        'time_to_peak_csi_h': eng.time_to_peak_csi_h(),
                        'csi_shape': eng.csi_regime(),
                    }
                    for b in 'UTCAGX':
                        rep_dict[f'comp_{b}'] = comp[b]
                    reps.append(rep_dict)

                def _mean(key):
                    return float(np.mean([r[key] for r in reps]))

                def _max_abs(key):
                    return float(np.max([abs(r[key]) for r in reps]))

                def _std(key):
                    return float(np.std([r[key] for r in reps]))

                surface['enrichment'][i, j] = _mean('enrichment')
                surface['enrichment_std'][i, j] = _std('enrichment')
                surface['dna_TU_over_free_TU'][i, j] = _mean('dna_TU_over_free_TU')
                surface['dna_TU_over_free_TU_std'][i, j] = _std('dna_TU_over_free_TU')
                surface['pool_drift_factor'][i, j] = _mean('pool_drift_factor')
                surface['free_TU_ratio'][i, j] = _mean('free_TU_ratio')
                surface['dna_TU_fraction'][i, j] = _mean('dna_TU_fraction')
                surface['fraction_T'][i, j] = _mean('fraction_T')
                surface['total_dna'][i, j] = _mean('total_dna')
                surface['peak_total_dna'][i, j] = _mean('peak_total_dna')
                surface['time_to_peak_h'][i, j] = _mean('time_to_peak_h')
                surface['final_survival'][i, j] = _mean('final_survival')
                surface['residue_survival'][i, j] = _mean('residue_survival')
                surface['residue_retention_ratio'][i, j] = _mean('residue_retention_ratio')
                surface['is_population_growth'][i, j] = bool(
                    np.mean([r['is_population_growth'] for r in reps]) > 0.5
                )
                surface['half_life_effective'][i, j] = _mean('half_life_eff')
                surface['half_life_true'][i, j] = _mean('half_life_true')
                surface['residue_half_life'][i, j] = _mean('residue_half_life')
                surface['stability_ratio'][i, j] = _mean('stability_ratio')
                surface['uv_ratio'][i, j] = _mean('uv_ratio')
                surface['poly_rate'][i, j] = _mean('poly_rate')
                surface['mass_balance_rel_err'][i, j] = _mean('mass_balance_rel_err')
                surface['consumed_residual'][i, j] = _mean('consumed_residual')
                surface['consumed_rel_err_mean'][i, j] = _mean('consumed_rel_err')
                surface['consumed_rel_err_max'][i, j] = _max_abs('consumed_rel_err')
                surface['deamination_transfer_residual'][i, j] = _mean(
                    'deam_transfer_residual'
                )
                surface['deamination_C_events'][i, j] = _mean('deam_C')
                surface['deamination_G_events'][i, j] = _mean('deam_G')
                surface['hydrolysis_events'][i, j] = _mean('hydro')
                surface['uv_damage_events'][i, j] = _mean('uv')
                regimes = [r['regime'] for r in reps]
                majority = max(set(regimes), key=regimes.count)
                surface['regime_code'][i, j] = regime_map[majority]
                surface['peak_csi'][i, j] = _mean('peak_csi')
                surface['time_to_peak_csi_h'][i, j] = _mean('time_to_peak_csi_h')
                shapes = [r['csi_shape'] for r in reps]
                majority_shape = max(set(shapes), key=shapes.count)
                surface['csi_shape_code'][i, j] = csi_shape_map[majority_shape]
                for b in 'UTCAGX':
                    surface['composition'][b][i, j] = _mean(f'comp_{b}')

            print(f"  [{label}] T={T_val:.1f}C done")
        return surface

    def run(self):
        self.surface_chemical = self._run_grid(
            use_incorporation_bias=False,
            use_fidelity_advantage=False,
            photon_flux_override=None,
            label="CHEMICAL_SELECTION",
        )
        self.surface_heuristic = self._run_grid(
            use_incorporation_bias=True,
            use_fidelity_advantage=True,
            photon_flux_override=None,
            label="HEURISTIC",
        )
        self.surface_zero_uv = self._run_grid(
            use_incorporation_bias=False,
            use_fidelity_advantage=False,
            photon_flux_override=0.0,
            label="ZERO_UV_CONTROL",
        )
        self.surface_equal_pools = self._run_grid(
            use_incorporation_bias=False,
            use_fidelity_advantage=False,
            photon_flux_override=0.0,
            equalize_ut_pools=True,
            disable_compartments=False,
            label="EQUAL_POOLS_CONTROL",
        )
        self._find_optimal()
        self._print_summary()

    # ----------------------------------------------------------------
    # Photon-flux / quantum-yield sensitivity sweep
    # ----------------------------------------------------------------
    def flux_sweep(
        self,
        fluxes: Optional[List[float]] = None,
        quantum_yields: Optional[List[float]] = None,
        n_rep: int = 1,
    ) -> Dict:
        if self.optimization_results is None:
            raise RuntimeError("Call run() before flux_sweep().")
        if fluxes is None:
            fluxes = [0.0, 1e14, 1e15, 1e16, 1e17, 1e18]
        if quantum_yields is None:
            quantum_yields = [1e-4, 1e-3, 1e-2]

        T_opt = self.optimization_results['optimal_T']
        pH_opt = self.optimization_results['optimal_pH']

        original_qy = EmpiricalRates.UV_QUANTUM_YIELD.value
        results = {
            'T_opt': T_opt,
            'pH_opt': pH_opt,
            'fluxes': list(fluxes),
            'quantum_yields': list(quantum_yields),
            'enrichment': np.zeros((len(fluxes), len(quantum_yields))),
            'enrichment_std': np.zeros((len(fluxes), len(quantum_yields))),
            'dna_TU_over_free_TU': np.zeros((len(fluxes), len(quantum_yields))),
            'pool_drift_factor': np.zeros((len(fluxes), len(quantum_yields))),
            'fraction_T': np.zeros((len(fluxes), len(quantum_yields))),
            'total_dna': np.zeros((len(fluxes), len(quantum_yields))),
            'sterilized': np.zeros((len(fluxes), len(quantum_yields)), dtype=bool),
            'hydrolysis_events': np.zeros((len(fluxes), len(quantum_yields))),
            'uv_damage_events': np.zeros((len(fluxes), len(quantum_yields))),
        }

        print(f"\n[FLUX_SWEEP] at T={T_opt:.1f} C, pH={pH_opt:.2f}")
        for qi, qy in enumerate(quantum_yields):
            EmpiricalRates.UV_QUANTUM_YIELD = Source(
                qy, "-",
                f"flux_sweep override (was {original_qy:g})",
                "imposed",
            )
            for fi, flux in enumerate(fluxes):
                reps = []
                for rep in range(n_rep):
                    seed = 4242 + fi * 100 + qi + rep * 10007
                    eng = VSSUFEngine(
                        temperature_C=T_opt,
                        pH_in_situ=pH_opt,
                        seed=seed,
                        max_time_hours=self.config['simulation_hours'],
                        verbose=False,
                        photon_flux=flux,
                        wavelength_nm=self.config['wavelength_nm'],
                        polymer_length=self.config['polymer_length'],
                        lipid_conc=self.config['lipid_conc'],
                        record_interval_seconds=self.config['record_interval_seconds'],
                        dt=self.config['dt'],
                        strand_state=self.config['strand_state'],
                        use_incorporation_bias=False,
                        use_fidelity_advantage=False,
                    )
                    eng.run()
                    reps.append({
                        'enrichment': eng.get_chemical_selection_index(),
                        'dna_TU_over_free_TU': eng.get_dna_TU_over_free_TU(),
                        'pool_drift_factor': eng.get_pool_drift_factor(),
                        'fraction_T': eng.get_final_thymine_fraction(),
                        'total_dna': eng.total_dna(),
                        'sterilized': eng.get_regime() == 'sterilized',
                        'hydro': eng.hydrolysis_events,
                        'uv': eng.uv_damage_events,
                    })
                results['enrichment'][fi, qi] = float(
                    np.mean([r['enrichment'] for r in reps])
                )
                results['enrichment_std'][fi, qi] = float(
                    np.std([r['enrichment'] for r in reps])
                )
                results['dna_TU_over_free_TU'][fi, qi] = float(
                    np.mean([r['dna_TU_over_free_TU'] for r in reps])
                )
                results['pool_drift_factor'][fi, qi] = float(
                    np.mean([r['pool_drift_factor'] for r in reps])
                )
                results['fraction_T'][fi, qi] = float(
                    np.mean([r['fraction_T'] for r in reps])
                )
                results['total_dna'][fi, qi] = float(
                    np.mean([r['total_dna'] for r in reps])
                )
                results['sterilized'][fi, qi] = bool(
                    np.mean([r['sterilized'] for r in reps]) > 0.5
                )
                results['hydrolysis_events'][fi, qi] = float(
                    np.mean([r['hydro'] for r in reps])
                )
                results['uv_damage_events'][fi, qi] = float(
                    np.mean([r['uv'] for r in reps])
                )
                flag = ' [STERILIZED]' if results['sterilized'][fi, qi] else ''
                print(
                    f"  flux={flux:.1e} qy={qy:.1e} "
                    f"-> CSI={results['enrichment'][fi, qi]:.3f}x "
                    f"real={results['dna_TU_over_free_TU'][fi, qi]:.3f}x "
                    f"total_dna={results['total_dna'][fi, qi]:.0f}{flag}"
                )
        EmpiricalRates.UV_QUANTUM_YIELD = Source(
            original_qy, "-", "assumed (order of magnitude)", "assumed"
        )
        self.flux_sweep_results = results
        return results

    # ----------------------------------------------------------------
    # Activation-energy sensitivity sweep (the 357x question)
    # ----------------------------------------------------------------
    def ea_T_sweep(
        self,
        ea_T_kcal: Optional[List[float]] = None,
        n_rep: int = 5,
    ) -> Dict:
        if self.optimization_results is None:
            raise RuntimeError("Call run() before ea_T_sweep().")
        if ea_T_kcal is None:
            ea_T_kcal = [28.0, 30.0, 31.0, 32.0, 33.0, 34.0, 36.0]

        T_opt = self.optimization_results['optimal_T']
        pH_opt = self.optimization_results['optimal_pH']

        original_ea_T = EmpiricalRates.Ea_T.value
        original_ea_T_src = EmpiricalRates.Ea_T
        results = {
            'T_opt': T_opt,
            'pH_opt': pH_opt,
            'ea_T_kcal': list(ea_T_kcal),
            'stability_ratio': np.zeros(len(ea_T_kcal)),
            'half_life_U_s': np.zeros(len(ea_T_kcal)),
            'half_life_T_s': np.zeros(len(ea_T_kcal)),
            'enrichment': np.zeros(len(ea_T_kcal)),
            'dna_TU_over_free_TU': np.zeros(len(ea_T_kcal)),
            'pool_drift_factor': np.zeros(len(ea_T_kcal)),
            'total_dna': np.zeros(len(ea_T_kcal)),
            'final_survival': np.zeros(len(ea_T_kcal)),
        }

        print(
            f"\n[EA_T_SWEEP] at T={T_opt:.1f} C, pH={pH_opt:.2f}  "
            f"(the 357x question)"
        )
        for si, ea in enumerate(ea_T_kcal):
            EmpiricalRates.Ea_T = Source(
                float(ea), "kcal/mol",
                f"ea_T_sweep override (was {original_ea_T:g})",
                "imposed",
            )
            ratio = EmpiricalRates.get_stability_ratio(T_opt, pH_opt)
            rates = EmpiricalRates.get_hydrolysis_rates(T_opt, pH_opt)
            hl_U = np.log(2) / rates['U'] if rates['U'] > 0 else float('inf')
            hl_T = np.log(2) / rates['T'] if rates['T'] > 0 else float('inf')
            reps = []
            for rep in range(n_rep):
                seed = 7777 + si * 100 + rep * 10007
                eng = VSSUFEngine(
                    temperature_C=T_opt,
                    pH_in_situ=pH_opt,
                    seed=seed,
                    max_time_hours=self.config['simulation_hours'],
                    verbose=False,
                    photon_flux=self.config['photon_flux'],
                    wavelength_nm=self.config['wavelength_nm'],
                    polymer_length=self.config['polymer_length'],
                    lipid_conc=self.config['lipid_conc'],
                    record_interval_seconds=self.config['record_interval_seconds'],
                    dt=self.config['dt'],
                    strand_state=self.config['strand_state'],
                )
                eng.run()
                reps.append({
                    'enrichment': eng.get_chemical_selection_index(),
                    'dna_TU_over_free_TU': eng.get_dna_TU_over_free_TU(),
                    'pool_drift_factor': eng.get_pool_drift_factor(),
                    'total_dna': eng.total_dna(),
                    'final_survival': eng.get_final_survival_fraction(),
                })
            results['stability_ratio'][si] = float(ratio)
            results['half_life_U_s'][si] = float(hl_U)
            results['half_life_T_s'][si] = float(hl_T)
            results['enrichment'][si] = float(np.mean([r['enrichment'] for r in reps]))
            results['dna_TU_over_free_TU'][si] = float(
                np.mean([r['dna_TU_over_free_TU'] for r in reps])
            )
            results['pool_drift_factor'][si] = float(
                np.mean([r['pool_drift_factor'] for r in reps])
            )
            results['total_dna'][si] = float(np.mean([r['total_dna'] for r in reps]))
            results['final_survival'][si] = float(
                np.mean([r['final_survival'] for r in reps])
            )
            print(
                f"  Ea_T={ea:.1f} kcal/mol  "
                f"k_U/k_T={ratio:8.1f}x  "
                f"t½_U={hl_U/3600:6.2f} h  "
                f"t½_T={hl_T/3600:8.2f} h  "
                f"CSI={results['enrichment'][si]:.3f}x  "
                f"real={results['dna_TU_over_free_TU'][si]:.3f}x  "
                f"drift={results['pool_drift_factor'][si]:.3f}x"
            )

        EmpiricalRates.Ea_T = original_ea_T_src
        self.ea_T_sweep_results = results
        return results

    # ----------------------------------------------------------------
    # v6.1.1: free-monomer deamination scale sweep (was a ghost call)
    # ----------------------------------------------------------------
    def free_monomer_deam_scale_sweep(
        self,
        scales: Optional[List[float]] = None,
        n_rep: int = 10,
    ) -> Dict:
        """
        Sweep EmpiricalRates.FREE_MONOMER_DEAM_SCALE at the optimum.

        The free-monomer deamination rates K_FREE_C_DEAM_REF and
        K_FREE_G_DEAM_REF are DEFINED as FREE_MONOMER_DEAM_SCALE * K_*_DEAM_REF
        (see EmpiricalRates class body). This sweep mutates both the scale
        and the two derived refs together, then restores all three.
        """
        if self.optimization_results is None:
            raise RuntimeError("Call run() before free_monomer_deam_scale_sweep().")
        if scales is None:
            scales = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]

        T_opt = self.optimization_results['optimal_T']
        pH_opt = self.optimization_results['optimal_pH']

        orig_scale = EmpiricalRates.FREE_MONOMER_DEAM_SCALE
        orig_kC    = EmpiricalRates.K_FREE_C_DEAM_REF
        orig_kG    = EmpiricalRates.K_FREE_G_DEAM_REF
        base_kC    = EmpiricalRates.K_C_DEAM_REF.value
        base_kG    = EmpiricalRates.K_G_DEAM_REF.value

        results = {
            'T_opt': T_opt,
            'pH_opt': pH_opt,
            'scales': list(scales),
            'k_free_C_deam_ref': np.zeros(len(scales)),
            'k_free_G_deam_ref': np.zeros(len(scales)),
            'enrichment': np.zeros(len(scales)),
            'enrichment_std': np.zeros(len(scales)),
            'dna_TU_over_free_TU': np.zeros(len(scales)),
            'pool_drift_factor': np.zeros(len(scales)),
            'free_TU_ratio': np.zeros(len(scales)),
            'fraction_T': np.zeros(len(scales)),
            'total_dna': np.zeros(len(scales)),
            'final_survival': np.zeros(len(scales)),
            'free_C_deam_events': np.zeros(len(scales)),
            'free_G_deam_events': np.zeros(len(scales)),
            'regime': [],
        }

        print(
            f"\n[FREE_DEAM_SCALE_SWEEP] at T={T_opt:.1f} C, pH={pH_opt:.2f}"
        )
        for si, scale in enumerate(scales):
            scale_f = float(max(0.0, scale))
            EmpiricalRates.FREE_MONOMER_DEAM_SCALE = Source(
                scale_f, "-",
                f"free_monomer_deam_scale_sweep override "
                f"(was {orig_scale.value:g})",
                "imposed",
            )
            EmpiricalRates.K_FREE_C_DEAM_REF = Source(
                base_kC * scale_f, "s^-1",
                f"derived: {scale_f:g} * K_C_DEAM_REF (sweep override)",
                "derived",
            )
            EmpiricalRates.K_FREE_G_DEAM_REF = Source(
                base_kG * scale_f, "s^-1",
                f"derived: {scale_f:g} * K_G_DEAM_REF (sweep override)",
                "derived",
            )

            reps = []
            for rep in range(max(1, n_rep)):
                seed = 5151 + si * 100 + rep * 10007
                eng = VSSUFEngine(
                    temperature_C=T_opt,
                    pH_in_situ=pH_opt,
                    seed=seed,
                    max_time_hours=self.config['simulation_hours'],
                    verbose=False,
                    photon_flux=self.config['photon_flux'],
                    wavelength_nm=self.config['wavelength_nm'],
                    polymer_length=self.config['polymer_length'],
                    lipid_conc=self.config['lipid_conc'],
                    record_interval_seconds=self.config['record_interval_seconds'],
                    dt=self.config['dt'],
                    strand_state=self.config['strand_state'],
                    use_incorporation_bias=False,
                    use_fidelity_advantage=False,
                )
                eng.run()
                reps.append({
                    'enrichment': eng.get_chemical_selection_index(),
                    'dna_TU_over_free_TU': eng.get_dna_TU_over_free_TU(),
                    'pool_drift_factor': eng.get_pool_drift_factor(),
                    'free_TU_ratio': eng.get_free_TU_ratio(),
                    'fraction_T': eng.get_final_thymine_fraction(),
                    'total_dna': eng.total_dna(),
                    'final_survival': eng.get_final_survival_fraction(),
                    'free_C_deam': eng.mass_free_C_to_U,
                    'free_G_deam': eng.mass_free_G_to_X,
                    'regime': eng.get_regime(),
                })

            results['k_free_C_deam_ref'][si] = base_kC * scale_f
            results['k_free_G_deam_ref'][si] = base_kG * scale_f
            results['enrichment'][si] = float(
                np.mean([r['enrichment'] for r in reps])
            )
            results['enrichment_std'][si] = float(
                np.std([r['enrichment'] for r in reps])
            )
            results['dna_TU_over_free_TU'][si] = float(
                np.mean([r['dna_TU_over_free_TU'] for r in reps])
            )
            results['pool_drift_factor'][si] = float(
                np.mean([r['pool_drift_factor'] for r in reps])
            )
            results['free_TU_ratio'][si] = float(
                np.mean([r['free_TU_ratio'] for r in reps])
            )
            results['fraction_T'][si] = float(
                np.mean([r['fraction_T'] for r in reps])
            )
            results['total_dna'][si] = float(
                np.mean([r['total_dna'] for r in reps])
            )
            results['final_survival'][si] = float(
                np.mean([r['final_survival'] for r in reps])
            )
            results['free_C_deam_events'][si] = float(
                np.mean([r['free_C_deam'] for r in reps])
            )
            results['free_G_deam_events'][si] = float(
                np.mean([r['free_G_deam'] for r in reps])
            )
            regimes = [r['regime'] for r in reps]
            results['regime'].append(max(set(regimes), key=regimes.count))

            print(
                f"  scale={scale_f:6.2f}  "
                f"k_free_C={base_kC*scale_f:.3e}  "
                f"k_free_G={base_kG*scale_f:.3e}  "
                f"CSI={results['enrichment'][si]:.3f}x "
                f"(±{results['enrichment_std'][si]:.3f})  "
                f"real={results['dna_TU_over_free_TU'][si]:.3f}x  "
                f"drift={results['pool_drift_factor'][si]:.3f}x  "
                f"free T/(U+T)={results['free_TU_ratio'][si]:.4f}  "
                f"DNA={results['total_dna'][si]:.0f}  "
                f"freeC_deam={results['free_C_deam_events'][si]:.0f}  "
                f"freeG_deam={results['free_G_deam_events'][si]:.0f}"
            )

        EmpiricalRates.FREE_MONOMER_DEAM_SCALE = orig_scale
        EmpiricalRates.K_FREE_C_DEAM_REF = orig_kC
        EmpiricalRates.K_FREE_G_DEAM_REF = orig_kG

        self.deam_scale_sweep_results = results
        return results

    # ----------------------------------------------------------------
    # v6.1.1: clay adsorption rate sweep (was a ghost call)
    # ----------------------------------------------------------------
    def clay_rate_sweep(
        self,
        rates: Optional[List[float]] = None,
        n_rep: int = 10,
    ) -> Dict:
        """
        Sweep EmpiricalRates.CLAY_ADSORPTION_RATE at the optimum.

        The clay desorption rate is held fixed at the class default, so
        increasing adsorption monotonically increases the fraction of DNA
        that enjoys CLAY_PROTECTION_FACTOR (9.5x) protection from
        hydrolysis and UV. This sweep quantifies that benefit and checks
        that the conclusion is not an artifact of the arbitrary 5e-3 s^-1
        default.
        """
        if self.optimization_results is None:
            raise RuntimeError("Call run() before clay_rate_sweep().")
        if rates is None:
            rates = [0.0, 1e-4, 1e-3, 5e-3, 1e-2, 1e-1]

        T_opt = self.optimization_results['optimal_T']
        pH_opt = self.optimization_results['optimal_pH']

        orig_ads = EmpiricalRates.CLAY_ADSORPTION_RATE
        results = {
            'T_opt': T_opt,
            'pH_opt': pH_opt,
            'rates': list(rates),
            'enrichment': np.zeros(len(rates)),
            'enrichment_std': np.zeros(len(rates)),
            'dna_TU_over_free_TU': np.zeros(len(rates)),
            'pool_drift_factor': np.zeros(len(rates)),
            'fraction_T': np.zeros(len(rates)),
            'total_dna': np.zeros(len(rates)),
            'final_survival': np.zeros(len(rates)),
            'clay_fraction_T': np.zeros(len(rates)),
            'clay_fraction_U': np.zeros(len(rates)),
            'adsorption_events': np.zeros(len(rates)),
            'desorption_events': np.zeros(len(rates)),
            'regime': [],
        }

        print(f"\n[CLAY_RATE_SWEEP] at T={T_opt:.1f} C, pH={pH_opt:.2f}")
        for ri, rate in enumerate(rates):
            rate_f = float(max(0.0, rate))
            EmpiricalRates.CLAY_ADSORPTION_RATE = Source(
                rate_f, "s^-1",
                f"clay_rate_sweep override (was {orig_ads.value:g})",
                "imposed",
            )
            reps = []
            for rep in range(max(1, n_rep)):
                seed = 6161 + ri * 100 + rep * 10007
                eng = VSSUFEngine(
                    temperature_C=T_opt,
                    pH_in_situ=pH_opt,
                    seed=seed,
                    max_time_hours=self.config['simulation_hours'],
                    verbose=False,
                    photon_flux=self.config['photon_flux'],
                    wavelength_nm=self.config['wavelength_nm'],
                    polymer_length=self.config['polymer_length'],
                    lipid_conc=self.config['lipid_conc'],
                    record_interval_seconds=self.config['record_interval_seconds'],
                    dt=self.config['dt'],
                    strand_state=self.config['strand_state'],
                    use_incorporation_bias=False,
                    use_fidelity_advantage=False,
                    disable_compartments=False,
                )
                eng.run()
                total = max(1, eng.total_dna())
                clay_T = (
                    eng.species['dsDNA_T_clay']
                    + eng.species['dsDNA_T_lipid']
                ) / total
                clay_U = (
                    eng.species['dsDNA_U_clay']
                    + eng.species['dsDNA_U_lipid']
                ) / total
                reps.append({
                    'enrichment': eng.get_chemical_selection_index(),
                    'dna_TU_over_free_TU': eng.get_dna_TU_over_free_TU(),
                    'pool_drift_factor': eng.get_pool_drift_factor(),
                    'fraction_T': eng.get_final_thymine_fraction(),
                    'total_dna': eng.total_dna(),
                    'final_survival': eng.get_final_survival_fraction(),
                    'clay_frac_T': clay_T,
                    'clay_frac_U': clay_U,
                    'ads': eng.clay_adsorption_events,
                    'des': eng.clay_desorption_events,
                    'regime': eng.get_regime(),
                })
            results['enrichment'][ri] = float(
                np.mean([r['enrichment'] for r in reps])
            )
            results['enrichment_std'][ri] = float(
                np.std([r['enrichment'] for r in reps])
            )
            results['dna_TU_over_free_TU'][ri] = float(
                np.mean([r['dna_TU_over_free_TU'] for r in reps])
            )
            results['pool_drift_factor'][ri] = float(
                np.mean([r['pool_drift_factor'] for r in reps])
            )
            results['fraction_T'][ri] = float(
                np.mean([r['fraction_T'] for r in reps])
            )
            results['total_dna'][ri] = float(
                np.mean([r['total_dna'] for r in reps])
            )
            results['final_survival'][ri] = float(
                np.mean([r['final_survival'] for r in reps])
            )
            results['clay_fraction_T'][ri] = float(
                np.mean([r['clay_frac_T'] for r in reps])
            )
            results['clay_fraction_U'][ri] = float(
                np.mean([r['clay_frac_U'] for r in reps])
            )
            results['adsorption_events'][ri] = float(
                np.mean([r['ads'] for r in reps])
            )
            results['desorption_events'][ri] = float(
                np.mean([r['des'] for r in reps])
            )
            regimes = [r['regime'] for r in reps]
            results['regime'].append(max(set(regimes), key=regimes.count))

            print(
                f"  k_ads={rate_f:.1e} s^-1  "
                f"CSI={results['enrichment'][ri]:.3f}x "
                f"(±{results['enrichment_std'][ri]:.3f})  "
                f"real={results['dna_TU_over_free_TU'][ri]:.3f}x  "
                f"drift={results['pool_drift_factor'][ri]:.3f}x  "
                f"clay frac T={results['clay_fraction_T'][ri]:.3f} "
                f"U={results['clay_fraction_U'][ri]:.3f}  "
                f"DNA={results['total_dna'][ri]:.0f}  "
                f"surv={results['final_survival'][ri]:.3f}"
            )

        EmpiricalRates.CLAY_ADSORPTION_RATE = orig_ads
        self.clay_rate_sweep_results = results
        return results

    # ----------------------------------------------------------------
    # Time-horizon sweep (240 h vs 2400 h vs 10000 h)
    # ----------------------------------------------------------------
    def time_horizon_sweep(
        self,
        horizons_h: Optional[List[float]] = None,
        n_rep: int = 10,
    ) -> Dict:
        if self.optimization_results is None:
            raise RuntimeError("Call run() before time_horizon_sweep().")
        if horizons_h is None:
            horizons_h = [240.0, 1000.0, 2400.0, 10000.0]

        T_opt = self.optimization_results['optimal_T']
        pH_opt = self.optimization_results['optimal_pH']

        results = {
            'T_opt': T_opt,
            'pH_opt': pH_opt,
            'horizons_h': list(horizons_h),
            'enrichment': np.zeros(len(horizons_h)),
            'enrichment_std': np.zeros(len(horizons_h)),
            'dna_TU_over_free_TU': np.zeros(len(horizons_h)),
            'pool_drift_factor': np.zeros(len(horizons_h)),
            'free_TU_ratio': np.zeros(len(horizons_h)),
            'total_dna': np.zeros(len(horizons_h)),
            'final_survival': np.zeros(len(horizons_h)),
            'peak_csi': np.zeros(len(horizons_h)),
            'time_to_peak_csi_h': np.zeros(len(horizons_h)),
            'end_csi': np.zeros(len(horizons_h)),
            'csi_shape': [],
        }

        print(
            f"\n[TIME_HORIZON_SWEEP] at T={T_opt:.1f} C, pH={pH_opt:.2f}"
        )
        for hi, horizon in enumerate(horizons_h):
            reps = []
            for rep in range(n_rep):
                seed = 8888 + hi * 100 + rep * 10007
                eng = VSSUFEngine(
                    temperature_C=T_opt,
                    pH_in_situ=pH_opt,
                    seed=seed,
                    max_time_hours=float(horizon),
                    verbose=False,
                    photon_flux=self.config['photon_flux'],
                    wavelength_nm=self.config['wavelength_nm'],
                    polymer_length=self.config['polymer_length'],
                    lipid_conc=self.config['lipid_conc'],
                    record_interval_seconds=self.config['record_interval_seconds'],
                    dt=self.config['dt'],
                    strand_state=self.config['strand_state'],
                )
                eng.run()
                reps.append({
                    'enrichment': eng.get_chemical_selection_index(),
                    'dna_TU_over_free_TU': eng.get_dna_TU_over_free_TU(),
                    'pool_drift_factor': eng.get_pool_drift_factor(),
                    'free_TU_ratio': eng.get_free_TU_ratio(),
                    'total_dna': eng.total_dna(),
                    'final_survival': eng.get_final_survival_fraction(),
                    'peak_csi': eng.peak_csi(),
                    'time_to_peak_csi_h': eng.time_to_peak_csi_h(),
                    'csi_shape': eng.csi_regime(),
                })
            results['enrichment'][hi] = float(np.mean([r['enrichment'] for r in reps]))
            results['enrichment_std'][hi] = float(np.std([r['enrichment'] for r in reps]))
            results['dna_TU_over_free_TU'][hi] = float(
                np.mean([r['dna_TU_over_free_TU'] for r in reps])
            )
            results['pool_drift_factor'][hi] = float(
                np.mean([r['pool_drift_factor'] for r in reps])
            )
            results['free_TU_ratio'][hi] = float(
                np.mean([r['free_TU_ratio'] for r in reps])
            )
            results['total_dna'][hi] = float(np.mean([r['total_dna'] for r in reps]))
            results['final_survival'][hi] = float(
                np.mean([r['final_survival'] for r in reps])
            )
            results['peak_csi'][hi] = float(np.mean([r['peak_csi'] for r in reps]))
            results['time_to_peak_csi_h'][hi] = float(
                np.mean([r['time_to_peak_csi_h'] for r in reps])
            )
            results['end_csi'][hi] = results['enrichment'][hi]
            shapes = [r['csi_shape'] for r in reps]
            results['csi_shape'].append(max(set(shapes), key=shapes.count))
            print(
                f"  horizon={horizon:7.0f} h  "
                f"CSI_end={results['enrichment'][hi]:.3f}x  "
                f"peak_CSI={results['peak_csi'][hi]:.3f}x "
                f"@ {results['time_to_peak_csi_h'][hi]:.0f} h  "
                f"real={results['dna_TU_over_free_TU'][hi]:.3f}x  "
                f"drift={results['pool_drift_factor'][hi]:.3f}x  "
                f"yield={results['total_dna'][hi]:.0f}  "
                f"surv={results['final_survival'][hi]:.3f}  "
                f"shape={results['csi_shape'][-1]}"
            )
        self.time_horizon_sweep_results = results
        return results

    # ----------------------------------------------------------------
    # Confirmatory stochastic validation at the selected cell
    # ----------------------------------------------------------------
    def confirm_optimum(self, n_rep: Optional[int] = None) -> Dict[str, float]:
        if self.optimization_results is None:
            raise RuntimeError("Call run() before confirm_optimum().")
        n_rep = int(n_rep or self.config.get('confirmatory_replicates', 100))
        T = float(self.optimization_results['optimal_T'])
        pH = float(self.optimization_results['optimal_pH'])
        vals = []
        for rep in range(max(1, n_rep)):
            eng = VSSUFEngine(
                temperature_C=T, pH_in_situ=pH,
                seed=900000 + rep * 7919,
                max_time_hours=self.config['simulation_hours'],
                verbose=False,
                photon_flux=self.config['photon_flux'],
                wavelength_nm=self.config['wavelength_nm'],
                polymer_length=self.config['polymer_length'],
                lipid_conc=self.config['lipid_conc'],
                record_interval_seconds=self.config['record_interval_seconds'],
                dt=self.config['dt'],
                strand_state=self.config['strand_state'],
                use_incorporation_bias=False,
                use_fidelity_advantage=False,
            )
            eng.run()
            vals.append({
                'enrichment': eng.get_chemical_selection_index(),
                'dna_TU_over_free_TU': eng.get_dna_TU_over_free_TU(),
                'pool_drift_factor': eng.get_pool_drift_factor(),
                'fraction_T': eng.get_final_thymine_fraction(),
                'total_dna': eng.total_dna(),
                'final_survival': eng.get_final_survival_fraction(),
                'mass_balance_rel_err': eng.max_abs_rel_error(),
                'consumed_residual': eng.max_abs_consumed_residual(),
                'regime': eng.get_regime(),
                'peak_csi': eng.peak_csi(),
                'time_to_peak_csi_h': eng.time_to_peak_csi_h(),
                'csi_shape': eng.csi_regime(),
            })
        e = np.array([v['enrichment'] for v in vals], dtype=float)
        r = np.array([v['dna_TU_over_free_TU'] for v in vals], dtype=float)
        ft = np.array([v['fraction_T'] for v in vals], dtype=float)
        dna = np.array([v['total_dna'] for v in vals], dtype=float)
        return {
            'n_replicates': int(len(vals)),
            'T': T, 'pH': pH,
            'enrichment_mean': float(np.mean(e)),
            'enrichment_std': float(np.std(e, ddof=1)) if len(e) > 1 else 0.0,
            'enrichment_cv': float(np.std(e, ddof=1) / abs(np.mean(e))) if len(e) > 1 and np.mean(e) != 0 else 0.0,
            'enrichment_q025': float(np.quantile(e, 0.025)),
            'enrichment_q50': float(np.quantile(e, 0.50)),
            'enrichment_q975': float(np.quantile(e, 0.975)),
            'dna_TU_over_free_TU_mean': float(np.mean(r)),
            'dna_TU_over_free_TU_std': float(np.std(r, ddof=1)) if len(r) > 1 else 0.0,
            'dna_TU_over_free_TU_q025': float(np.quantile(r, 0.025)),
            'dna_TU_over_free_TU_q50': float(np.quantile(r, 0.50)),
            'dna_TU_over_free_TU_q975': float(np.quantile(r, 0.975)),
            'fraction_T_mean': float(np.mean(ft)),
            'total_dna_mean': float(np.mean(dna)),
            'total_dna_q025': float(np.quantile(dna, 0.025)),
            'total_dna_q975': float(np.quantile(dna, 0.975)),
            'max_mass_balance_rel_err': float(max(v['mass_balance_rel_err'] for v in vals)),
            'max_consumed_residual': float(max(v['consumed_residual'] for v in vals)),
            'regime_frequencies': {
                r_: sum(v['regime'] == r_ for v in vals) / len(vals)
                for r_ in sorted({v['regime'] for v in vals})
            },
            'csi_shape_frequencies': {
                s_: sum(v['csi_shape'] == s_ for v in vals) / len(vals)
                for s_ in sorted({v['csi_shape'] for v in vals})
            },
        }

    # ----------------------------------------------------------------
    # Optimal-cell selection
    # ----------------------------------------------------------------
    def _find_optimal(self):
        s = self.surface_chemical
        min_total = self.config.get('min_total_dna_for_optimum', 100.0)
        eligible = s['total_dna'] >= min_total
        if self.config.get('exclude_post_peak_collapse', True):
            eligible &= s['regime_code'] != 3
        if self.config.get('exclude_sterilized', True):
            eligible &= s['regime_code'] != 0
        if not np.any(eligible):
            print(
                f"  WARNING: no grid point has total_dna >= {min_total}. "
                f"Falling back to unconstrained argmax."
            )
            eligible = np.ones_like(s['total_dna'], dtype=bool)

        masked = np.where(eligible, s['enrichment'], -np.inf)
        max_enrich = float(np.max(masked))
        if max_enrich <= 0:
            idx = np.unravel_index(np.argmax(masked), s['enrichment'].shape)
        else:
            tol = self.ENRICHMENT_TIE_TOL * abs(max_enrich)
            near = eligible & (s['enrichment'] >= max_enrich - tol)
            score = (
                np.where(near, s['total_dna'], -np.inf)
                + 1e-6 * np.where(near, s['residue_retention_ratio'], -np.inf)
            )
            idx = np.unravel_index(np.argmax(score), s['enrichment'].shape)

        sh = self.surface_heuristic
        sz = self.surface_zero_uv
        sp = self.surface_equal_pools
        input_frac_T = EmpiricalRates.get_input_thymine_fraction()

        chem = float(s['enrichment'][idx])
        heur = float(sh['enrichment'][idx])
        delta = heur - chem
        if abs(delta) < 1e-9:
            interp = "heuristic has no effect at this cell"
        elif delta > 0:
            interp = (
                f"heuristic adds {delta:+.3f}x "
                f"({100.0*delta/max(heur,1e-12):.1f}% of total)"
            )
        else:
            interp = (
                f"heuristic REDUCES enrichment by {delta:.3f}x "
                f"(imposed bias is counter-productive here)"
            )

        self.optimization_results = {
            'optimal_T': float(s['T'][idx[0]]),
            'optimal_pH': float(s['pH'][idx[1]]),
            'min_total_dna_for_optimum': float(min_total),
            'input_fraction_T': float(input_frac_T),
            'max_enrichment_chemical': chem,
            'max_enrichment_heuristic': heur,
            'max_enrichment_zero_uv': float(sz['enrichment'][idx]),
            'max_enrichment_equal_pools': float(sp['enrichment'][idx]),
            'heuristic_delta': delta,
            'heuristic_interpretation': interp,
            'real_enrichment_chemical': float(s['dna_TU_over_free_TU'][idx]),
            'real_enrichment_heuristic': float(sh['dna_TU_over_free_TU'][idx]),
            'real_enrichment_zero_uv': float(sz['dna_TU_over_free_TU'][idx]),
            'real_enrichment_equal_pools': float(sp['dna_TU_over_free_TU'][idx]),
            'pool_drift_factor_chemical': float(s['pool_drift_factor'][idx]),
            'pool_drift_factor_zero_uv': float(sz['pool_drift_factor'][idx]),
            'pool_drift_factor_equal_pools': float(sp['pool_drift_factor'][idx]),
            'fraction_T_chemical': float(s['fraction_T'][idx]),
            'fraction_T_heuristic': float(sh['fraction_T'][idx]),
            'fraction_T_zero_uv': float(sz['fraction_T'][idx]),
            'fraction_T_equal_pools': float(sp['fraction_T'][idx]),
            'total_dna_chemical': float(s['total_dna'][idx]),
            'total_dna_heuristic': float(sh['total_dna'][idx]),
            'total_dna_zero_uv': float(sz['total_dna'][idx]),
            'total_dna_equal_pools': float(sp['total_dna'][idx]),
            'final_survival_chemical': float(s['final_survival'][idx]),
            'final_survival_heuristic': float(sh['final_survival'][idx]),
            'final_survival_zero_uv': float(sz['final_survival'][idx]),
            'residue_survival_chemical': float(s['residue_survival'][idx]),
            'residue_survival_heuristic': float(sh['residue_survival'][idx]),
            'residue_survival_zero_uv': float(sz['residue_survival'][idx]),
            'residue_retention_ratio_chemical': float(s['residue_retention_ratio'][idx]),
            'residue_retention_ratio_heuristic': float(sh['residue_retention_ratio'][idx]),
            'residue_retention_ratio_zero_uv': float(sz['residue_retention_ratio'][idx]),
            'is_population_growth_chemical': bool(s['is_population_growth'][idx]),
            'is_population_growth_heuristic': bool(sh['is_population_growth'][idx]),
            'is_population_growth_zero_uv': bool(sz['is_population_growth'][idx]),
            'half_life_effective_s': float(s['half_life_effective'][idx]),
            'half_life_true_s': float(s['half_life_true'][idx]),
            'residue_half_life_s': float(s['residue_half_life'][idx]),
            'stability_ratio': float(s['stability_ratio'][idx]),
            'stability_ratio_label': (
                'extrapolated (25 C refs -> 90 C; Ea_T assumed 32 kcal/mol; '
                'see ea_T_sweep for uncertainty)'
            ),
            'uv_ratio_U_T': float(s['uv_ratio'][idx]),
            'poly_rate': float(s['poly_rate'][idx]),
            'half_life_U_s': float(
                np.log(2) / EmpiricalRates.get_hydrolysis_rates(
                    float(s['T'][idx[0]]), float(s['pH'][idx[1]])
                )['U']
            ),
            'half_life_T_s': float(
                np.log(2) / EmpiricalRates.get_hydrolysis_rates(
                    float(s['T'][idx[0]]), float(s['pH'][idx[1]])
                )['T']
            ),
            'max_mass_balance_rel_err': float(s['mass_balance_rel_err'][idx]),
            'max_mass_balance_rel_err_heuristic': float(sh['mass_balance_rel_err'][idx]),
            'max_consumed_residual': float(s['consumed_residual'][idx]),
            'max_consumed_rel_err_mean': float(s['consumed_rel_err_mean'][idx]),
            'max_consumed_rel_err_max': float(s['consumed_rel_err_max'][idx]),
            'deamination_transfer_residual': float(s['deamination_transfer_residual'][idx]),
            'deamination_C_events_chemical': float(s['deamination_C_events'][idx]),
            'deamination_G_events_chemical': float(s['deamination_G_events'][idx]),
            'deamination_C_events_heuristic': float(sh['deamination_C_events'][idx]),
            'deamination_G_events_heuristic': float(sh['deamination_G_events'][idx]),
            'hydrolysis_events_chemical': float(s['hydrolysis_events'][idx]),
            'uv_damage_events_chemical': float(s['uv_damage_events'][idx]),
            'hydrolysis_events_zero_uv': float(sz['hydrolysis_events'][idx]),
            'uv_damage_events_zero_uv': float(sz['uv_damage_events'][idx]),
            'composition_chemical': {
                b: float(s['composition'][b][idx]) for b in 'UTCAGX'
            },
            'composition_heuristic': {
                b: float(sh['composition'][b][idx]) for b in 'UTCAGX'
            },
            'composition_zero_uv': {
                b: float(sz['composition'][b][idx]) for b in 'UTCAGX'
            },
            'composition_equal_pools': {
                b: float(sp['composition'][b][idx]) for b in 'UTCAGX'
            },
            'peak_csi_chemical': float(s['peak_csi'][idx]),
            'time_to_peak_csi_h_chemical': float(s['time_to_peak_csi_h'][idx]),
            'csi_shape_chemical': self._csi_shape_name(int(s['csi_shape_code'][idx])),
        }

    @staticmethod
    def _csi_shape_name(code: int) -> str:
        return {
            0: 'no_signal', 1: 'steady', 2: 'rising',
            3: 'peaked', 4: 'collapsed',
        }.get(code, 'unknown')

    def _print_summary(self):
        o = self.optimization_results
        print("\n" + "=" * 78)
        print("UPDSF v6.1.1 — VENT CHEMICAL FILTERING (T over U)")
        print("=" * 78)
        print(
            f"Optimal condition: T = {o['optimal_T']:.1f} C, "
            f"pH_in_situ = {o['optimal_pH']:.2f}"
        )
        print(f"Simulation time: {self.config['simulation_hours']:.0f} h")
        print(
            f"Photon flux: {self.config['photon_flux']:.2e} "
            f"photons/cm^2/s @ {self.config['wavelength_nm']:.0f} nm "
            f"({'VENT (dark)' if self.config['photon_flux'] == 0 else 'SURFACE'})"
        )
        print()
        print("PRIMARY RESULT (chemical selection only, vent):")
        print(
            f"  CSI (total env. selection): "
            f"{o['max_enrichment_chemical']:.3f}x"
        )
        print(
            f"  Real T-pref enrichment:    "
            f"{o['real_enrichment_chemical']:.3f}x"
        )
        print(
            f"  Pool-composition drift:    "
            f"{o['pool_drift_factor_chemical']:.3f}x"
        )
        print(
            f"  CSI trajectory: peak={o['peak_csi_chemical']:.3f}x "
            f"@ {o['time_to_peak_csi_h_chemical']:.0f} h, "
            f"shape={o['csi_shape_chemical']}"
        )
        print()
        print("MECHANISM:")
        print(
            f"  Stability ratio k_U/k_T: {o['stability_ratio']:.2f}  "
            f"[{o['stability_ratio_label']}]"
        )
        print(
            f"  Hydrolysis half-life U: "
            f"{o['half_life_U_s']/3600:.2f} h"
        )
        print(
            f"  Hydrolysis half-life T: "
            f"{o['half_life_T_s']/3600:.2f} h"
        )
        print(
            f"  Hydrolysis events: {o['hydrolysis_events_chemical']:.0f}   "
            f"UV events: {o['uv_damage_events_chemical']:.0f}"
        )
        print()
        print("SURVIVAL / YIELD:")
        print(
            f"  Per-polymer survival (final): "
            f"{o['final_survival_chemical']:.4f}"
        )
        print(
            f"  Total DNA (final): "
            f"{o['total_dna_chemical']:.0f}"
        )
        print()
        print("CONFIDENCE:")
        print(
            f"  Max mass-balance |rel err|: "
            f"{o['max_mass_balance_rel_err']:.3e}"
        )
        print(
            f"  Max |consumed identity residual|: "
            f"{o['max_consumed_residual']:.3e}"
        )
        print("=" * 78)

    def save(self, base_path: str = "output_data_v610"):
        os.makedirs(base_path, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out = {
            'version': '6.1.1',
            'scientific_claim': (
                "Chemical filtering can enrich T relative to U through "
                "differential hydrolysis stability in dark hydrothermal-vent "
                "conditions. UV is negligible in this environment (zero-UV "
                "control differs by <2%). The model does not demonstrate "
                "Darwinian evolution. The CSI decomposes multiplicatively as "
                "CSI = real T-pref enrichment * pool-drift factor. Both "
                "factors are environmental. The stability ratio k_U/k_T is "
                "an Arrhenius extrapolation, not a measurement; see the "
                "Ea_T sweep for its uncertainty."
            ),
            'config': self.config,
            'optimization_results': self.optimization_results,
            'confirmatory_optimum': getattr(self, 'confirmatory_optimum_results', None),
            'surface_chemical_selection': self._surface_to_dict(self.surface_chemical),
            'surface_heuristic': self._surface_to_dict(self.surface_heuristic),
            'surface_zero_uv_control': self._surface_to_dict(self.surface_zero_uv),
            'surface_equal_pools_control': self._surface_to_dict(self.surface_equal_pools),
        }
        if self.flux_sweep_results is not None:
            out['flux_sweep'] = self._surface_to_dict(self.flux_sweep_results)
        if self.deam_scale_sweep_results is not None:
            out['free_monomer_deam_scale_sweep'] = (
                self._surface_to_dict(self.deam_scale_sweep_results)
            )
        if self.clay_rate_sweep_results is not None:
            out['clay_rate_sweep'] = self._surface_to_dict(self.clay_rate_sweep_results)
        if self.ea_T_sweep_results is not None:
            out['ea_T_sweep'] = self._surface_to_dict(self.ea_T_sweep_results)
        if self.time_horizon_sweep_results is not None:
            out['time_horizon_sweep'] = self._surface_to_dict(
                self.time_horizon_sweep_results
            )
        path = f"{base_path}/analysis_v611_{timestamp}.json"
        with open(path, 'w') as f:
            json.dump(out, f, indent=2, default=str)
        print(f"JSON saved: {path}")

        rows = []
        for label, s in (
            ('chemical_selection', self.surface_chemical),
            ('heuristic', self.surface_heuristic),
            ('zero_uv_control', self.surface_zero_uv),
            ('equal_pools_control', self.surface_equal_pools),
        ):
            for i, T in enumerate(s['T']):
                for j, pH in enumerate(s['pH']):
                    row = {
                        'model': label,
                        'T': float(T),
                        'pH': float(pH),
                        'photon_flux': float(s['photon_flux']),
                    }
                    keys = [
                        'enrichment', 'enrichment_std',
                        'dna_TU_over_free_TU', 'dna_TU_over_free_TU_std',
                        'pool_drift_factor', 'free_TU_ratio',
                        'dna_TU_fraction', 'fraction_T',
                        'total_dna', 'peak_total_dna', 'time_to_peak_h',
                        'final_survival', 'residue_survival',
                        'residue_retention_ratio',
                        'half_life_effective', 'half_life_true',
                        'residue_half_life',
                        'stability_ratio', 'uv_ratio', 'poly_rate',
                        'mass_balance_rel_err', 'consumed_residual',
                        'consumed_rel_err_mean', 'consumed_rel_err_max',
                        'deamination_transfer_residual',
                        'deamination_C_events', 'deamination_G_events',
                        'hydrolysis_events', 'uv_damage_events',
                        'peak_csi', 'time_to_peak_csi_h',
                    ]
                    for k in keys:
                        row[k] = float(s[k][i, j])
                    row['is_population_growth'] = bool(
                        s['is_population_growth'][i, j]
                    )
                    row['regime_code'] = int(s['regime_code'][i, j])
                    row['csi_shape_code'] = int(s['csi_shape_code'][i, j])
                    for b in 'UTCAGX':
                        row[f'comp_{b}'] = float(s['composition'][b][i, j])
                    rows.append(row)
        df = pd.DataFrame(rows)
        csv_path = f"{base_path}/sensitivity_v611_{timestamp}.csv"
        df.to_csv(csv_path, index=False)
        print(f"CSV saved: {csv_path}")

    @staticmethod
    def _surface_to_dict(s: Dict) -> Dict:
        out = {}
        for k, v in s.items():
            if isinstance(v, np.ndarray):
                out[k] = v.tolist()
            elif isinstance(v, dict):
                out[k] = {}
                for kk, vv in v.items():
                    if isinstance(vv, np.ndarray):
                        out[k][kk] = vv.tolist()
                    else:
                        out[k][kk] = vv
            elif isinstance(v, (bool, np.bool_)):
                out[k] = bool(v)
            else:
                out[k] = v
        return out

    # ----------------------------------------------------------------
    # Plotting
    # ----------------------------------------------------------------
    def plot(self, save_path: str = "sensitivity_v611.png"):
        s = self.surface_chemical
        sh = self.surface_heuristic
        sz = self.surface_zero_uv
        sp = self.surface_equal_pools

        fig, axes = plt.subplots(3, 4, figsize=(26, 15))
        fig.suptitle(
            "UPDSF v6.1.1 — Vent Chemical Filtering (T over U)",
            fontsize=14, fontweight='bold',
        )

        ax = axes[0, 0]
        c = ax.contourf(s['T'], s['pH'], s['enrichment'].T, levels=20, cmap='viridis')
        ax.set_xlabel("T (C)"); ax.set_ylabel("pH")
        ax.set_title("CSI (total env. selection)")
        fig.colorbar(c, ax=ax)

        ax = axes[0, 1]
        c = ax.contourf(s['T'], s['pH'], s['dna_TU_over_free_TU'].T, levels=20, cmap='viridis')
        ax.set_xlabel("T (C)"); ax.set_ylabel("pH")
        ax.set_title("Real T-pref enrichment (direct)")
        fig.colorbar(c, ax=ax)

        ax = axes[0, 2]
        c = ax.contourf(s['T'], s['pH'], s['pool_drift_factor'].T, levels=20, cmap='plasma')
        ax.set_xlabel("T (C)"); ax.set_ylabel("pH")
        ax.set_title("Pool-composition drift (indirect)")
        fig.colorbar(c, ax=ax)

        ax = axes[0, 3]
        c = ax.contourf(s['T'], s['pH'], s['fraction_T'].T, levels=20, cmap='Blues')
        ax.set_xlabel("T (C)"); ax.set_ylabel("pH")
        ax.set_title("Fraction T in DNA")
        fig.colorbar(c, ax=ax)

        ax = axes[1, 0]
        diff = sh['enrichment'] - s['enrichment']
        c = ax.contourf(sh['T'], sh['pH'], diff.T, levels=20, cmap='RdBu_r')
        ax.set_xlabel("T (C)"); ax.set_ylabel("pH")
        ax.set_title("Heuristic contribution (red=helps)")
        fig.colorbar(c, ax=ax)

        ax = axes[1, 1]
        c = ax.contourf(sz['T'], sz['pH'], sz['enrichment'].T, levels=20, cmap='viridis')
        ax.set_xlabel("T (C)"); ax.set_ylabel("pH")
        ax.set_title("CSI — zero-UV control")
        fig.colorbar(c, ax=ax)

        ax = axes[1, 2]
        c = ax.contourf(sp['T'], sp['pH'], sp['enrichment'].T, levels=20, cmap='viridis')
        ax.set_xlabel("T (C)"); ax.set_ylabel("pH")
        ax.set_title("CSI — equal-pools control")
        fig.colorbar(c, ax=ax)

        ax = axes[1, 3]
        total = s['hydrolysis_events'] + s['uv_damage_events'] + 1e-9
        uv_frac = s['uv_damage_events'] / total
        c = ax.contourf(s['T'], s['pH'], uv_frac.T, levels=20, cmap='coolwarm')
        ax.set_xlabel("T (C)"); ax.set_ylabel("pH")
        ax.set_title("Fraction of destruction events from UV")
        fig.colorbar(c, ax=ax)

        ax = axes[2, 0]
        residual = np.log10(s['consumed_rel_err_max'] + 1e-30)
        c = ax.contourf(s['T'], s['pH'], residual.T, levels=20, cmap='Reds')
        ax.set_xlabel("T (C)"); ax.set_ylabel("pH")
        ax.set_title("log10 max |consumed identity rel err|")
        fig.colorbar(c, ax=ax)

        ax = axes[2, 1]
        c = ax.contourf(s['T'], s['pH'], s['residue_retention_ratio'].T, levels=20, cmap='RdYlGn')
        ax.set_xlabel("T (C)"); ax.set_ylabel("pH")
        ax.set_title("Residue retention ratio")
        fig.colorbar(c, ax=ax)

        ax = axes[2, 2]
        c = ax.contourf(s['T'], s['pH'], s['peak_csi'].T, levels=20, cmap='viridis')
        ax.set_xlabel("T (C)"); ax.set_ylabel("pH")
        ax.set_title("Peak CSI over the run")
        fig.colorbar(c, ax=ax)

        ax = axes[2, 3]
        c = ax.contourf(s['T'], s['pH'], s['time_to_peak_csi_h'].T, levels=20, cmap='magma')
        ax.set_xlabel("T (C)"); ax.set_ylabel("pH")
        ax.set_title("Time to peak CSI (h)")
        fig.colorbar(c, ax=ax)

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Plot saved: {save_path}")

    def plot_flux_sweep(self, save_path: str = "flux_sweep_v611.png"):
        if self.flux_sweep_results is None:
            print("No flux-sweep results to plot.")
            return
        fs = self.flux_sweep_results
        fig, axes = plt.subplots(1, 4, figsize=(24, 5.5))
        fig.suptitle(
            "UPDSF v6.1.1 — Robustness to assumed UV scale (surface only)",
            fontsize=13, fontweight='bold',
        )
        for qi, qy in enumerate(fs['quantum_yields']):
            axes[0].plot(fs['fluxes'], fs['enrichment'][:, qi], 'o-', label=f'qy={qy:.0e}')
            axes[1].plot(fs['fluxes'], fs['dna_TU_over_free_TU'][:, qi], 'o-', label=f'qy={qy:.0e}')
            axes[2].plot(fs['fluxes'], fs['uv_damage_events'][:, qi], 'o-', label=f'qy={qy:.0e}')
            axes[3].plot(fs['fluxes'], fs['total_dna'][:, qi], 'o-', label=f'qy={qy:.0e}')
        for ax in axes:
            ax.set_xscale('symlog', linthresh=1e9)
            ax.set_xlabel("photon_flux (photons/cm^2/s)")
            ax.legend()
            ax.grid(alpha=0.3)
        axes[0].set_ylabel("CSI"); axes[0].set_title("CSI vs UV")
        axes[1].set_ylabel("Real"); axes[1].set_title("Real vs UV")
        axes[2].set_ylabel("UV damage events"); axes[2].set_title("UV events vs UV")
        axes[3].set_ylabel("Final total DNA"); axes[3].set_title("Yield vs UV")
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Flux-sweep plot saved: {save_path}")

    def plot_ea_T_sweep(self, save_path: str = "ea_T_sweep_v611.png"):
        if self.ea_T_sweep_results is None:
            print("No Ea_T-sweep results to plot.")
            return
        e = self.ea_T_sweep_results
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        fig.suptitle(
            "UPDSF v6.1.1 — Sensitivity to assumed Ea_T (the 357x question)",
            fontsize=13, fontweight='bold',
        )
        axes[0].plot(e['ea_T_kcal'], e['stability_ratio'], 'o-', color='darkred')
        axes[0].set_yscale('log')
        axes[0].set_xlabel("Ea_T (kcal/mol)")
        axes[0].set_ylabel("k_U / k_T (log scale)")
        axes[0].set_title("Extrapolated stability ratio")
        axes[0].grid(alpha=0.3)

        axes[1].plot(e['ea_T_kcal'], e['enrichment'], 'o-', label='CSI', color='darkred')
        axes[1].plot(e['ea_T_kcal'], e['dna_TU_over_free_TU'], 's--', label='real', color='darkblue')
        axes[1].plot(e['ea_T_kcal'], e['pool_drift_factor'], '^:', label='drift', color='darkgreen')
        axes[1].set_xlabel("Ea_T (kcal/mol)")
        axes[1].set_ylabel("Multiplier")
        axes[1].set_title("CSI decomposition vs Ea_T")
        axes[1].legend()
        axes[1].grid(alpha=0.3)

        axes[2].plot(e['ea_T_kcal'], e['total_dna'], 'o-', color='purple')
        axes[2].set_xlabel("Ea_T (kcal/mol)")
        axes[2].set_ylabel("Final total DNA")
        axes[2].set_title("Yield vs Ea_T")
        axes[2].grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Ea_T-sweep plot saved: {save_path}")

    def plot_time_horizon_sweep(self, save_path: str = "time_horizon_sweep_v611.png"):
        if self.time_horizon_sweep_results is None:
            print("No time-horizon-sweep results to plot.")
            return
        t = self.time_horizon_sweep_results
        fig, axes = plt.subplots(1, 4, figsize=(24, 5.5))
        fig.suptitle(
            "UPDSF v6.1.1 — CSI vs time horizon at the optimum",
            fontsize=13, fontweight='bold',
        )
        x = t['horizons_h']
        axes[0].plot(x, t['enrichment'], 'o-', label='CSI_end', color='darkred')
        axes[0].plot(x, t['peak_csi'], 's--', label='peak CSI', color='darkorange')
        axes[0].set_xscale('log')
        axes[0].set_xlabel("Horizon (h)")
        axes[0].set_ylabel("CSI")
        axes[0].set_title("CSI vs horizon")
        axes[0].legend()
        axes[0].grid(alpha=0.3)

        axes[1].plot(x, t['dna_TU_over_free_TU'], 'o-', label='real', color='darkblue')
        axes[1].plot(x, t['pool_drift_factor'], 's--', label='drift', color='darkgreen')
        axes[1].plot(x, t['free_TU_ratio'] / max(t['free_TU_ratio']), '^:', label='free T/(U+T) (norm.)', color='gray')
        axes[1].set_xscale('log')
        axes[1].set_xlabel("Horizon (h)")
        axes[1].set_ylabel("Multiplier")
        axes[1].set_title("Real vs drift vs horizon")
        axes[1].legend()
        axes[1].grid(alpha=0.3)

        axes[2].plot(x, t['total_dna'], 'o-', color='purple')
        axes[2].set_xscale('log')
        axes[2].set_xlabel("Horizon (h)")
        axes[2].set_ylabel("Final total DNA")
        axes[2].set_title("Yield vs horizon")
        axes[2].grid(alpha=0.3)

        axes[3].plot(x, t['final_survival'], 'o-', color='brown')
        axes[3].set_xscale('log')
        axes[3].set_xlabel("Horizon (h)")
        axes[3].set_ylabel("Final per-polymer survival")
        axes[3].set_title("Survival vs horizon")
        axes[3].grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Time-horizon-sweep plot saved: {save_path}")

    def plot_free_deam_scale_sweep(
        self, save_path: str = "free_deam_scale_sweep_v611.png"
    ):
        # v6.1.1: plot the free-monomer deamination scale sweep
        if self.deam_scale_sweep_results is None:
            print("No free-monomer deam-sweep results to plot.")
            return
        r = self.deam_scale_sweep_results
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        fig.suptitle(
            "UPDSF v6.1.1 — Sensitivity to free-monomer deamination scale",
            fontsize=13, fontweight='bold',
        )
        x = r['scales']
        axes[0].errorbar(x, r['enrichment'], yerr=r['enrichment_std'],
                         fmt='o-', color='darkred', capsize=3, label='CSI')
        axes[0].plot(x, r['dna_TU_over_free_TU'], 's--',
                     color='darkblue', label='real')
        axes[0].plot(x, r['pool_drift_factor'], '^:',
                     color='darkgreen', label='drift')
        axes[0].set_xlabel("FREE_MONOMER_DEAM_SCALE")
        axes[0].set_ylabel("Multiplier")
        axes[0].set_title("CSI decomposition vs free-deam scale")
        axes[0].legend(); axes[0].grid(alpha=0.3)

        axes[1].plot(x, r['free_TU_ratio'], 'o-',
                     color='purple', label='free T/(U+T)')
        axes[1].plot(x, r['fraction_T'], 's--',
                     color='orange', label='DNA T/(U+T)')
        axes[1].set_xlabel("FREE_MONOMER_DEAM_SCALE")
        axes[1].set_ylabel("T/(U+T)")
        axes[1].set_title("Free vs DNA T/(U+T)")
        axes[1].legend(); axes[1].grid(alpha=0.3)

        axes[2].plot(x, r['total_dna'], 'o-',
                     color='purple', label='total DNA')
        axes[2].set_xlabel("FREE_MONOMER_DEAM_SCALE")
        axes[2].set_ylabel("Final total DNA")
        axes[2].set_title("Yield vs free-deam scale")
        axes[2].grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Free-deam-sweep plot saved: {save_path}")

    def plot_clay_rate_sweep(
        self, save_path: str = "clay_rate_sweep_v611.png"
    ):
        # v6.1.1: plot the clay adsorption rate sweep
        if self.clay_rate_sweep_results is None:
            print("No clay-sweep results to plot.")
            return
        r = self.clay_rate_sweep_results
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        fig.suptitle(
            "UPDSF v6.1.1 — Sensitivity to clay adsorption rate",
            fontsize=13, fontweight='bold',
        )
        x = r['rates']
        axes[0].errorbar(x, r['enrichment'], yerr=r['enrichment_std'],
                         fmt='o-', color='darkred', capsize=3, label='CSI')
        axes[0].plot(x, r['dna_TU_over_free_TU'], 's--',
                     color='darkblue', label='real')
        axes[0].plot(x, r['pool_drift_factor'], '^:',
                     color='darkgreen', label='drift')
        axes[0].set_xscale('symlog', linthresh=1e-4)
        axes[0].set_xlabel("CLAY_ADSORPTION_RATE (s^-1)")
        axes[0].set_ylabel("Multiplier")
        axes[0].set_title("CSI decomposition vs clay adsorption")
        axes[0].legend(); axes[0].grid(alpha=0.3)

        axes[1].plot(x, r['clay_fraction_T'], 'o-',
                     color='tab:red', label='T on clay/lipid')
        axes[1].plot(x, r['clay_fraction_U'], 's--',
                     color='tab:blue', label='U on clay/lipid')
        axes[1].set_xscale('symlog', linthresh=1e-4)
        axes[1].set_xlabel("CLAY_ADSORPTION_RATE (s^-1)")
        axes[1].set_ylabel("Fraction in protected compartment")
        axes[1].set_title("Compartmentalization vs clay adsorption")
        axes[1].legend(); axes[1].grid(alpha=0.3)

        axes[2].plot(x, r['total_dna'], 'o-',
                     color='purple', label='total DNA')
        max_ads = max(1.0, float(np.max(r['adsorption_events'])))
        axes[2].plot(x, r['adsorption_events'] / max_ads,
                     's:', color='gray', label='ads events (norm.)')
        axes[2].set_xscale('symlog', linthresh=1e-4)
        axes[2].set_xlabel("CLAY_ADSORPTION_RATE (s^-1)")
        axes[2].set_ylabel("Final total DNA")
        axes[2].set_title("Yield vs clay adsorption")
        axes[2].legend(); axes[2].grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Clay-sweep plot saved: {save_path}")

    def plot_csi_trajectory(
        self,
        save_path: str = "csi_trajectory_v611.png",
        max_time_hours: Optional[float] = None,
    ):
        """Single long run at the optimum, CSI / real / drift vs time."""
        if self.optimization_results is None:
            raise RuntimeError("Call run() before plot_csi_trajectory().")
        T = float(self.optimization_results['optimal_T'])
        pH = float(self.optimization_results['optimal_pH'])
        hours = max_time_hours or max(
            240.0, float(self.config.get('time_horizon_for_plot', 2400.0))
        )
        eng = VSSUFEngine(
            temperature_C=T, pH_in_situ=pH,
            seed=12345,
            max_time_hours=hours,
            verbose=False,
            photon_flux=self.config['photon_flux'],
            wavelength_nm=self.config['wavelength_nm'],
            polymer_length=self.config['polymer_length'],
            lipid_conc=self.config['lipid_conc'],
            record_interval_seconds=self.config['record_interval_seconds'],
            dt=self.config['dt'],
            strand_state=self.config['strand_state'],
        )
        eng.run()
        time_h = np.array(eng.history['time']) / 3600.0
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(
            f"UPDSF v6.1.1 — CSI trajectory at T={T:.0f} C, pH={pH:.1f} "
            f"(horizon {hours:.0f} h, shape={eng.csi_regime()})",
            fontsize=13, fontweight='bold',
        )

        ax = axes[0, 0]
        ax.plot(time_h, eng.history['csi'], color='darkred', lw=2, label='CSI (total)')
        ax.plot(time_h, eng.history['real'], color='darkblue', lw=2, ls='--', label='real (direct)')
        ax.plot(time_h, eng.history['drift'], color='darkgreen', lw=2, ls=':', label='drift (indirect)')
        ax.axhline(1.0, color='k', ls='--', alpha=0.3)
        ax.axvline(eng.time_to_peak_csi_h(), color='darkorange', ls='-.', alpha=0.6,
                   label=f"peak CSI @ {eng.time_to_peak_csi_h():.0f} h")
        ax.set_xlabel("Time (h)"); ax.set_ylabel("Multiplier")
        ax.set_title("CSI = real × drift")
        ax.legend(); ax.grid(alpha=0.3)

        ax = axes[0, 1]
        ax.plot(time_h, eng.history['free_T_frac'], color='purple', lw=2, label='free T/(U+T)')
        ax.plot(time_h, eng.history['dna_TU_fraction'], color='orange', lw=2, ls='--', label='DNA T/(U+T)')
        ax.plot(time_h, eng.history['pool_influx_T_frac'], color='gray', lw=2, ls=':', label='influx T/(U+T)')
        ax.set_xlabel("Time (h)"); ax.set_ylabel("T/(U+T)")
        ax.set_title("Free vs DNA vs influx T/(U+T)")
        ax.legend(); ax.grid(alpha=0.3)

        ax = axes[1, 0]
        ax.plot(time_h, eng.history['total_dna'], color='black', lw=2)
        ax.set_xlabel("Time (h)"); ax.set_ylabel("Total DNA molecules")
        ax.set_title("Yield")
        ax.grid(alpha=0.3)

        ax = axes[1, 1]
        ax.plot(time_h, eng.history['cohort_survival'], color='darkblue', lw=2, label='per-polymer')
        ax.plot(time_h, eng.history['residue_survival'], color='darkgreen', lw=2, ls='--', label='residue')
        ax.set_xlabel("Time (h)"); ax.set_ylabel("Survival")
        ax.set_title("Survival")
        ax.legend(); ax.grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"CSI-trajectory plot saved: {save_path}  (horizon {hours:.0f} h)")
        return eng


# ====================================================================
# SECTION 8: OPTIMAL SIMULATION PLOTTER
# ====================================================================

def optimal_simulation(
    analyzer: TwoDSensitivityAnalyzer,
    save_path: str = "optimal_simulation_v611.png",
    show_fig: bool = False,
    use_incorporation_bias: bool = False,
    use_fidelity_advantage: bool = False,
    max_time_hours: Optional[float] = None,
):
    opt = analyzer.optimization_results
    T_opt = opt['optimal_T']
    pH_opt = opt['optimal_pH']
    hours = max_time_hours or analyzer.config['simulation_hours']
    variant = 'chemical selection (vent)' if not use_incorporation_bias else 'heuristic'

    print(
        f"\n[Optimal simulation] T={T_opt:.1f} C, "
        f"pH_in_situ={pH_opt:.2f} ({variant}, {hours:.0f} h)"
    )

    T_vals = analyzer.surface_chemical['T']
    pH_vals = analyzer.surface_chemical['pH']
    i_opt = int(np.argmin(np.abs(T_vals - T_opt)))
    j_opt = int(np.argmin(np.abs(pH_vals - pH_opt)))
    seed = 42 + i_opt * 100 + j_opt

    eng = VSSUFEngine(
        temperature_C=T_opt,
        pH_in_situ=pH_opt,
        seed=seed,
        max_time_hours=hours,
        verbose=False,
        photon_flux=analyzer.config['photon_flux'],
        wavelength_nm=analyzer.config['wavelength_nm'],
        polymer_length=analyzer.config['polymer_length'],
        lipid_conc=analyzer.config['lipid_conc'],
        record_interval_seconds=analyzer.config['record_interval_seconds'],
        dt=analyzer.config['dt'],
        strand_state=analyzer.config['strand_state'],
        use_incorporation_bias=use_incorporation_bias,
        use_fidelity_advantage=use_fidelity_advantage,
    )
    eng.run()
    time_h = np.array(eng.history['time']) / 3600.0

    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.30)
    fig.suptitle(
        f"UPDSF v6.1.1 — Vent Chemical Filtering (T over U)\n"
        f"T={T_opt:.1f}°C, pH={pH_opt:.2f}, horizon={hours:.0f} h  ({variant})",
        fontsize=14, fontweight='bold',
    )

    colors = {
        'U': 'tab:blue', 'T': 'tab:red', 'C': 'tab:green',
        'A': 'tab:orange', 'G': 'tab:purple', 'X': 'tab:brown',
    }

    ax1 = fig.add_subplot(gs[0, 0])
    for b in 'UTCAGX':
        ax1.plot(time_h, eng.history[f'dsDNA_{b}'], color=colors[b],
                 label=f'{b}-DNA', lw=1.8)
    ax1.set_xlabel("Time (h)"); ax1.set_ylabel("DNA molecules")
    ax1.set_title("DNA steady state")
    ax1.legend(loc='upper left', fontsize=9); ax1.grid(alpha=0.3)

    ax2 = fig.add_subplot(gs[0, 1])
    for b in 'UTCAGX':
        ax2.plot(time_h, eng.history['fractions'][b], color=colors[b],
                 label=f'{b}', lw=1.8)
    ax2.set_xlabel("Time (h)"); ax2.set_ylabel("Fraction of total DNA")
    ax2.set_title("Nucleotide fractions")
    ax2.legend(loc='upper right', fontsize=9, ncol=2); ax2.grid(alpha=0.3)

    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(time_h, eng.history['csi'], color='darkred', lw=2, label='CSI (total)')
    ax3.plot(time_h, eng.history['real'], color='darkblue', lw=2, ls='--', label='real (direct)')
    ax3.plot(time_h, eng.history['drift'], color='darkgreen', lw=2, ls=':', label='drift (indirect)')
    ax3.axhline(1.0, color='k', ls='--', alpha=0.4, label='no enrichment')
    ax3.axvline(eng.time_to_peak_csi_h(), color='darkorange', ls='-.', alpha=0.5,
                label=f"peak CSI @ {eng.time_to_peak_csi_h():.0f} h")
    ax3.set_xlabel("Time (h)"); ax3.set_ylabel("Multiplier")
    ax3.set_title(f"CSI = real × drift  [{eng.csi_regime()}]")
    ax3.legend(loc='lower right', fontsize=8); ax3.grid(alpha=0.3)

    ax4 = fig.add_subplot(gs[1, 0])
    ax4.plot(time_h, eng.history['cohort_survival'], color='darkblue', lw=2, label='per-polymer')
    ax4.plot(time_h, eng.history['residue_survival'], color='darkgreen', lw=2, ls='--', label='residue')
    ax4.plot(time_h, eng.history['residue_retention_ratio'], color='darkred', lw=2, ls=':', label='retention')
    ax4.axhline(1.0, color='k', ls='--', alpha=0.3)
    ax4.set_xlabel("Time (h)"); ax4.set_ylabel("Survival")
    ax4.set_title(f"Survival ({eng.get_regime()})")
    ax4.legend(loc='lower left', fontsize=8); ax4.grid(alpha=0.3)
    ax4.set_ylim(-0.05, 1.05)

    ax5 = fig.add_subplot(gs[1, 1])
    comp = eng.get_nucleotide_composition()
    bases = list(comp.keys()); vals = [comp[b] for b in bases]
    bar_colors = [colors[b] for b in bases]
    bars = ax5.bar(bases, vals, color=bar_colors, alpha=0.85, edgecolor='black')
    ax5.set_xlabel("Nucleotide"); ax5.set_ylabel("Molecules")
    ax5.set_title("Final composition"); ax5.grid(alpha=0.3, axis='y')
    max_val = max(vals) if vals else 1
    for bar, val in zip(bars, vals):
        ax5.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() + max_val * 0.01,
                 f'{val:.0f}', ha='center', va='bottom', fontsize=9)

    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis('off')
    half_U = eng.get_hydrolysis_half_life('U')
    half_T = eng.get_hydrolysis_half_life('T')
    stab = EmpiricalRates.get_stability_ratio(T_opt, pH_opt)
    summary = (
        "UPDSF v6.1.1 — Vent chemical filtering of T over U\n"
        "====================================================\n\n"
        "Conditions:\n"
        f"  T = {T_opt:.1f} °C\n"
        f"  pH_in_situ = {pH_opt:.2f}\n"
        f"  photon_flux = {analyzer.config['photon_flux']:.2e}  "
        f"({'VENT (dark)' if analyzer.config['photon_flux']==0 else 'SURFACE'})\n"
        f"  horizon = {hours:.0f} h\n"
        f"  regime = {eng.get_regime()}\n"
        f"  CSI shape = {eng.csi_regime()}\n\n"
        "Two environmental channels:\n"
        f"  CSI_end            = {eng.get_chemical_selection_index():.3f}×\n"
        f"  peak CSI           = {eng.peak_csi():.3f}× "
        f"@ {eng.time_to_peak_csi_h():.0f} h\n"
        f"  real (direct)      = {eng.get_dna_TU_over_free_TU():.3f}×\n"
        f"  drift (indirect)   = {eng.get_pool_drift_factor():.3f}×\n"
        f"  free T/(U+T)       = {eng.get_free_TU_ratio():.4f}\n"
        f"  DNA  T/(U+T)       = {eng.get_dna_TU_fraction():.4f}\n"
        f"  influx T/(U+T)     = {eng._pool_influx_T_frac():.4f}\n\n"
        "Mechanism (extrapolated):\n"
        f"  k_U/k_T            = {stab:.1f}×  [Ea_T assumed; see sweep]\n"
        f"  t½(U) hydrolysis   = {half_U/3600:.2f} h\n"
        f"  t½(T) hydrolysis   = {half_T/3600:.2f} h\n"
        f"  hydrolysis events  = {eng.hydrolysis_events:,}\n"
        f"  UV damage events   = {eng.uv_damage_events:,}\n\n"
        "Yield / survival:\n"
        f"  total DNA (final)  = {eng.total_dna():,}\n"
        f"  per-polymer surv   = {eng.get_final_survival_fraction():.4f}\n"
        f"  residue survival   = {eng.get_residue_survival():.4f}\n\n"
        "Confidence:\n"
        f"  max |mass balance| = {eng.max_abs_rel_error():.2e}\n"
        f"  max |consumed res| = {eng.max_abs_consumed_residual():.2e}\n\n"
        "Interpretation:\n"
        f"  CSI = {eng.get_dna_TU_over_free_TU():.2f}× real × "
        f"{eng.get_pool_drift_factor():.2f}× drift.\n"
        f"  At 90°C / pH 11 UV is negligible: the effect\n"
        f"  is hydrolytic in origin.\n"
        f"  The {stab:.0f}× stability ratio is an Arrhenius\n"
        f"  extrapolation; the Ea_T sweep bounds it.\n"
    )
    ax6.text(0.0, 1.0, summary, ha='left', va='top',
             transform=ax6.transAxes, fontsize=8.0, family='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow',
                       alpha=0.95, edgecolor='gray'))

    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    if show_fig:
        try:
            matplotlib.use('TkAgg')
        except Exception:
            pass
        plt.show()
    plt.close(fig)
    print(f"Optimal simulation plot saved: {save_path}")
    return eng


# ====================================================================
# SECTION 9: MAIN
# ====================================================================

def _self_check_mass_identity():
    eng = VSSUFEngine(
        temperature_C=40.0, pH_in_situ=10.5,
        seed=7, max_time_hours=1.0, dt=60.0,
    )
    eng.run()
    err = eng.max_abs_consumed_residual()
    if err > 1e-9:
        raise AssertionError(f"Consumed-mass identity failed: {err:.6g}")
    return True


def _self_check_flags():
    eng1 = VSSUFEngine(
        temperature_C=76.0, pH_in_situ=9.5,
        seed=1, max_time_hours=0.1, equalize_ut_pools=True,
    )
    eng1.run()
    if 'equalize_ut_pools' not in eng1._flags_seen:
        raise AssertionError("equalize_ut_pools flag not applied.")
    eng2 = VSSUFEngine(
        temperature_C=76.0, pH_in_situ=9.5,
        seed=2, max_time_hours=0.1, disable_compartments=True,
    )
    eng2.run()
    if 'disable_compartments' not in eng2._flags_seen:
        raise AssertionError("disable_compartments flag not applied.")
    return True


def _self_check_source_file():
    try:
        with open(__file__, 'r', encoding='utf-8') as f:
            source = f.read()
        ast.parse(source)
    except SyntaxError as e:
        print(
            f"\n[SELF-CHECK FAILED] Syntax error at line {e.lineno}: {e.msg}",
            file=sys.stderr,
        )
        sys.exit(2)
    except OSError:
        pass


def main():
    _self_check_source_file()
    _self_check_mass_identity()
    _self_check_flags()

    print("=" * 78)
    print("UPDSF v6.1.1 — Vent Chemical Filtering (T over U)")
    print("Mass-Consistent | Kinetically-Consistent | Provenance-Tracked")
    print("=" * 78)

    print_provenance_table(EmpiricalRates.provenance(), "EmpiricalRates")
    print_provenance_table(LipidMembrane.provenance(), "LipidMembrane")
    print_provenance_table(PolymerPhysics.provenance(), "PolymerPhysics")
    all_sources = {}
    all_sources.update(EmpiricalRates.provenance())
    all_sources.update(LipidMembrane.provenance())
    all_sources.update(PolymerPhysics.provenance())
    print('PROVENANCE COUNTS:', provenance_summary(all_sources))

    # ---- Primary single run: VENT (dark), 240 h ----
    print(
        "\n[Single run — VENT] T=90.0 C, pH_in_situ=11.0, "
        "photon_flux=0 (dark), 240 h"
    )
    eng = VSSUFEngine(
        temperature_C=90.0,
        pH_in_situ=11.0,
        seed=42,
        max_time_hours=240.0,
        verbose=True,
        photon_flux=0.0,
        wavelength_nm=260.0,
        polymer_length=100,
        lipid_conc=0.05,
        strand_state='dsDNA',
        use_incorporation_bias=False,
        use_fidelity_advantage=False,
    )
    eng.run()
    print(f"\nFinal composition (molecules): {eng.get_nucleotide_composition()}")
    print(f"Total DNA molecules: {eng.total_dna():,}")
    print(f"Total DNA residues:  {eng.total_dna_residues():,}")
    print(f"Regime: {eng.get_regime()}   CSI shape: {eng.csi_regime()}")
    print(
        f"CSI (end): {eng.get_chemical_selection_index():.3f}x   "
        f"peak: {eng.peak_csi():.3f}x @ {eng.time_to_peak_csi_h():.0f} h"
    )
    print(
        f"Real: {eng.get_dna_TU_over_free_TU():.3f}x   "
        f"Drift: {eng.get_pool_drift_factor():.3f}x"
    )
    print(f"Free T/(U+T): {eng.get_free_TU_ratio():.4f}")
    print(f"DNA  T/(U+T): {eng.get_dna_TU_fraction():.4f}")
    print(
        f"Stability ratio k_U/k_T: "
        f"{EmpiricalRates.get_stability_ratio(90.0, 11.0):.2f}  [extrapolated]"
    )
    print(
        f"Max |rel err| mass balance: "
        f"{eng.max_abs_rel_error():.3e}"
    )
    print(
        f"Max |consumed identity residual|: "
        f"{eng.max_abs_consumed_residual():.3e}"
    )
    eng.save_history_to_csv("output_data_v610")

    # ---- Single-run comparison: SURFACE (1e16), same T/pH ----
    print(
        "\n[Single run — SURFACE comparison] T=90.0 C, pH_in_situ=11.0, "
        "photon_flux=1e16, 240 h"
    )
    eng_surf = VSSUFEngine(
        temperature_C=90.0,
        pH_in_situ=11.0,
        seed=42,
        max_time_hours=240.0,
        verbose=False,
        photon_flux=1e16,
        wavelength_nm=260.0,
        polymer_length=100,
        lipid_conc=0.05,
        strand_state='dsDNA',
    )
    eng_surf.run()
    print(
        f"  SURFACE: CSI={eng_surf.get_chemical_selection_index():.3f}x "
        f"real={eng_surf.get_dna_TU_over_free_TU():.3f}x "
        f"drift={eng_surf.get_pool_drift_factor():.3f}x "
        f"(UV events={eng_surf.uv_damage_events:,})"
    )
    print(
        f"  VENT   : CSI={eng.get_chemical_selection_index():.3f}x "
        f"real={eng.get_dna_TU_over_free_TU():.3f}x "
        f"drift={eng.get_pool_drift_factor():.3f}x "
        f"(UV events={eng.uv_damage_events:,})"
    )
    delta_csi = abs(
        eng_surf.get_chemical_selection_index()
        - eng.get_chemical_selection_index()
    )
    print(
        f"  |ΔCSI| surface vs vent: {delta_csi:.4f}x  "
        f"({100*delta_csi/max(eng.get_chemical_selection_index(),1e-12):.2f}%)"
    )

    # ---- Sensitivity analysis ----
    config = {
        'temperature_min': 40.0,
        'temperature_max': 90.0,
        'temperature_step': 5.0,
        'pH_min': 5.0,
        'pH_max': 11.0,
        'pH_step': 0.5,
        'simulation_hours': 240.0,
        'n_replicates': 5,
        'photon_flux': 0.0,
        'wavelength_nm': 260.0,
        'polymer_length': 100,
        'lipid_conc': 0.05,
        'record_interval_seconds': 600.0,
        'dt': 60.0,
        'strand_state': 'dsDNA',
        'min_total_dna_for_optimum': 100.0,
        'exclude_post_peak_collapse': True,
        'exclude_sterilized': True,
        'confirmatory_replicates': 200,
        'time_horizon_for_plot': 2400.0,
    }

    analyzer = TwoDSensitivityAnalyzer(config)
    analyzer.run()
    analyzer.confirmatory_optimum_results = analyzer.confirm_optimum()
    print(
        'Confirmatory optimum statistics:',
        analyzer.confirmatory_optimum_results
    )

    print("\n" + "=" * 78)
    print("Running Ea_T sensitivity sweep (the 357x question)")
    print("=" * 78)
    analyzer.ea_T_sweep(
        ea_T_kcal=[28.0, 30.0, 31.0, 32.0, 33.0, 34.0, 36.0],
        n_rep=5,
    )

    print("\n" + "=" * 78)
    print("Running time-horizon sweep (240 h vs 2400 h vs 10000 h)")
    print("=" * 78)
    analyzer.time_horizon_sweep(
        horizons_h=[240.0, 1000.0, 2400.0, 10000.0],
        n_rep=10,
    )

    print("\n" + "=" * 78)
    print("Running UV robustness sweep (surface only; should be ~flat at vent optimum)")
    print("=" * 78)
    analyzer.flux_sweep(
        fluxes=[0.0, 1e14, 1e15, 1e16, 1e17, 1e18],
        quantum_yields=[1e-4, 1e-3, 1e-2],
        n_rep=20,
    )

    # v6.1.1: these two methods now exist on TwoDSensitivityAnalyzer
    print("\n" + "=" * 78)
    print("Running free-monomer deamination scale sweep")
    print("=" * 78)
    analyzer.free_monomer_deam_scale_sweep(
        scales=[0.0, 0.5, 1.0, 2.0, 5.0, 10.0],
        n_rep=10,
    )

    print("\n" + "=" * 78)
    print("Running clay adsorption rate sweep")
    print("=" * 78)
    analyzer.clay_rate_sweep(
        rates=[0.0, 1e-4, 1e-3, 5e-3, 1e-2, 1e-1],
        n_rep=10,
    )

    analyzer.save("output_data_v610")
    analyzer.plot("sensitivity_v611.png")
    analyzer.plot_flux_sweep("flux_sweep_v611.png")
    analyzer.plot_ea_T_sweep("ea_T_sweep_v611.png")
    analyzer.plot_time_horizon_sweep("time_horizon_sweep_v611.png")
    analyzer.plot_free_deam_scale_sweep("free_deam_scale_sweep_v611.png")
    analyzer.plot_clay_rate_sweep("clay_rate_sweep_v611.png")
    analyzer.plot_csi_trajectory(
        "csi_trajectory_v611.png",
        max_time_hours=config['time_horizon_for_plot'],
    )

    print("\n" + "=" * 78)
    print("Running optimal simulation (vent, 240 h)")
    print("=" * 78)
    optimal_simulation(
        analyzer,
        save_path="optimal_simulation_v611_vent_240h.png",
        show_fig=False,
        use_incorporation_bias=False,
        use_fidelity_advantage=False,
        max_time_hours=240.0,
    )

    print("\n" + "=" * 78)
    print("Running optimal simulation (vent, 2400 h — the long horizon)")
    print("=" * 78)
    optimal_simulation(
        analyzer,
        save_path="optimal_simulation_v611_vent_2400h.png",
        show_fig=False,
        use_incorporation_bias=False,
        use_fidelity_advantage=False,
        max_time_hours=2400.0,
    )

    print("\n" + "=" * 78)
    print("DONE. v6.1.1 summary:")
    print("  - Default photon_flux = 0 (dark hydrothermal vent).")
    print("  - UV machinery retained only for surface-scenario comparison.")
    print("  - CSI trajectory tracked: peak, time-to-peak, shape.")
    print("  - Ea_T sweep exposes the ~357x extrapolation uncertainty.")
    print("  - Time-horizon sweep answers the 240 h vs 2400 h question.")
    print("  - Free-monomer deam-scale sweep (v6.1.1: was a ghost call, now works).")
    print("  - Clay adsorption-rate sweep (v6.1.1: was a ghost call, now works).")
    print("  - Two new sweep plots (free-deam, clay).")
    print("  - stability_ratio labeled 'extrapolated' with sensitivity band.")
    print("  - All plots at 300 DPI.")
    print("=" * 78)


if __name__ == "__main__":
    main()
