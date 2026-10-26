[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23656782-blue)](https://doi.org/10.5281/zenodo.23656782) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)


[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23655168-blue)](https://doi.org/10.5281/zenodo.23655168) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)



# UPDSF_v6_1_1
UPDSF v6.1.1 — the Unified Prebiotic DNA Selection Framework — a mass-conserving, kinetically-consistent stochastic model of hydrothermal-vent chemistry. At the optimum (90 °C, pH 11, dark vent).


Environment: 🐍 Python 3.8+

## Author: Seyed Mohammad Reza Hashemi (Reza Hashemi)

## ORCID: 0009-0002-0645-5180

## A mass-consistent, kinetically-consistent, and provenance-tracked stochastic simulation framework designed to model the chemical selection of Thymine (T) over Uracil (U) in prebiotic hydrothermal environments.

## Scientific Novelty 

The central problem—the "Thymine Mystery"—is one of the most significant gaps in prebiotic chemistry. Most models rely on "selection after the fact" (Darwinian evolution). This framework shifts the paradigm by demonstrating that selection is a prerequisite for evolution, not just a consequence of it.

By proposing a two-stage model (Hydrolytic Selection in the dark → UV Selection in the light), the model has created a "Unified Theory" of the RNA-to-DNA transition. This theoretically resolves the long-standing tension between "Dark Chemistry" (hydrothermal vents) and "Light Evolution" (surface environments).

## Overview
The origin of the transition from RNA to DNA remains a central puzzle in astrobiology. While Darwinian evolution refines biological systems, the initial selection mechanism must be purely chemical.

UPDSF v6.1.1 demonstrates that chemical filtering in dark, alkaline hydrothermal vents can drive a significant enrichment of Thymine over Uracil (CSI ≈ 6.019×) through the multiplicative combination of direct preferential hydrolysis and pool-composition drift.

## Key Features
Stochastic Simulation Engine: Utilizes Kinetic Monte Carlo methods (Binomial destruction) to model molecular populations with high statistical robustness.

Mass-Consistent Modeling: Explicitly tracks the mass balance of all nucleobases (U, T, C, A, G, X) to ensure physical rigor.

Provenance Tracking: Every parameter is encapsulated in a Source object, documenting whether it is derived from literature, assumed, or extrapolated.

## Two-Stage Prebiotic Model:

#### Stage 1 (Vent Chemistry): Hydrolytic chemical filtering in dark, hydrothermal environments.

#### Stage 2 (Surface Evolution): UV-driven evolutionary selection (complementary to PNAS studies by Khosh Abady et al., 2026).



## Sensitivity Analysis Suite: Includes automated sweeps for:

Temperature and pH grids (2D optimization).

Activation energy (Eₐ) sensitivity (the "357× stability ratio" question).

UV photon flux robustness.

Free-monomer deamination scales.

Clay adsorption/protection rates.

## Mathematical Foundation
Arrhenius Extrapolation
The model calculates temperature-dependent rates (k) using the classical Arrhenius equation:

text
k(T) = k_ref · exp[ −(Eₐ / R) · (1/T − 1/T_ref) ]
where:

k_ref = reference rate constant at temperature T_ref

Eₐ = activation energy (J·mol⁻¹)

R = universal gas constant (8.314 J·mol⁻¹·K⁻¹)

T = absolute temperature (K)

T_ref = reference absolute temperature (K)

Chemical Selection Index (CSI)
The framework defines the enrichment mechanism through a multiplicative decomposition:

text
CSI = Real (Direct Selection) × Drift (Indirect Composition Shift)
where:

Real (Direct Selection) captures the preferential hydrolysis advantage of T over U.

Drift (Indirect Composition Shift) captures the pool-composition drift arising from differential depletion of the U precursor pool.


## Results & Visualization
The framework automatically generates high-resolution (300 DPI) plots, including:

CSI Trajectories: Tracking enrichment over time.

Nucleotide Fractions: Visualizing the rise of T and fall of U.

2D Optimization Maps: Contour plots showing optimal T and pH regimes.

Sensitivity Heatmaps: Impact of Eₐ and UV flux on stability.

## Citation
If you use this framework in your research, please cite the following:

Reza Hashemi. (2026). Chemical Filtering of Thymine over Uracil in Prebiotic Hydrothermal Systems: A Mass-Consistent Stochastic Model with Complementary Insights from UV Photodamage Pathways. Zenodo. https://doi.org/10.5281/zenodo.23656782

Reza Hashemi. (2026). mrhashemi2000/UPDSF_v6_1_1: Initial release [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23655168

## License
This work is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0) license.

See LICENSE for full terms.

Author: Seyed Mohammad Reza Hashemi

ORCID: 0009-0002-0645-5180
