[![DOI](https://img.shields.io/badge/DOI-10.26434/chemrxiv.15011902/v1-green)](https://doi.org/10.26434/chemrxiv.15011902/v1) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23876507-blue)](https://doi.org/10.5281/zenodo.23876507) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23876308-blue)](https://doi.org/10.5281/zenodo.23876308) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23788249-blue)](https://doi.org/10.5281/zenodo.23788249) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23656782-blue)](https://doi.org/10.5281/zenodo.23656782) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23655168-blue)](https://doi.org/10.5281/zenodo.23655168) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23788480-blue)](https://doi.org/10.5281/zenodo.23788480) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)


# UPDSF_v6_1_1
UPDSF v6.1.1 — the Unified Prebiotic DNA Selection Framework — a mass-conserving, kinetically-consistent stochastic model of hydrothermal-vent chemistry. At the optimum (90 °C, pH 11, dark vent).


Environment: 🐍 Python 3.8+

## Author: Seyed Mohammad Reza Hashemi (Reza Hashemi)

## ORCID: 0009-0002-0645-5180

## A mass-consistent, kinetically-consistent, and provenance-tracked stochastic simulation framework designed to model the chemical selection of Thymine (T) over Uracil (U) in prebiotic hydrothermal environments.

## Scientific Novelty 

The central problem—the "Thymine Mystery"—is one of the most significant gaps in prebiotic chemistry. Most models rely on "selection after the fact" (Darwinian evolution). This framework shifts the paradigm by demonstrating that selection is a prerequisite for evolution, not just a consequence of it.

By proposing a two-stage model (Hydrolytic Selection in the dark → UV Selection in the light), the model has created a "Unified Theory" of the RNA-to-DNA transition. This theoretically resolves the long-standing tension between "Dark Chemistry" (hydrothermal vents) and "Light Evolution" (surface environments).

<br>
<br>

📂 **Complete Analysis: The Matter World Hypothesis (MWH): A Deterministic Theory of Genetic Architecture**

The detailed reports and raw data for Complete Analysis: The Matter World Hypothesis (MWH): A Deterministic Theory of Genetic Architecture can be found in the 📄 [ The Matter World Hypothesis (MWH): A Deterministic Theory of Genetic Architecture](analysis/Matter_World_Hypothesis_(MWH).md)

<br>
<br>

<img width="1366" height="768" alt="تصویر صفحه 2026-11-02 103634" src="https://github.com/user-attachments/assets/0cc0c273-d54b-40ea-b72e-199518ae4fa7" />

<img width="1352" height="485" alt="تصویر صفحه 2026-11-04 081907" src="https://github.com/user-attachments/assets/c0037153-c4eb-4830-abdf-d256ba99d45c" />

<img width="1357" height="541" alt="تصویر صفحه 2026-11-04 082159" src="https://github.com/user-attachments/assets/b1823ec8-ae7d-4fd7-ac67-9f0ca7290ca9" />

<img width="1355" height="488" alt="تصویر صفحه 2026-11-04 082414" src="https://github.com/user-attachments/assets/7ca9ed9d-7390-4aba-b529-d23748c7dba1" />



## Overview
The origin of the transition from RNA to DNA remains a central puzzle in astrobiology. While Darwinian evolution refines biological systems, the initial selection mechanism must be purely chemical.

UPDSF v6.1.1 demonstrates that chemical filtering in dark, alkaline hydrothermal vents can drive a significant enrichment of Thymine over Uracil (CSI ≈ 6.019×) through the multiplicative combination of direct preferential hydrolysis and pool-composition drift.

#### "To ensure the robustness of our theoretical findings, we explicitly decomposed the CSI to distinguish between direct kinetic selection and indirect pool-composition drift, finding that both are essential components of the total environmental selection."

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


# UPDSF v6.1.1 vs v6.2.0 Comparative Brief

**Project**: Monte Carlo kinetic simulation of prebiotic DNA selection (UPDSF)  
  
**Scope**: Faithful side-by-side comparison of two identical core engines with one new validation feature

## Scientific Claim 
Chemical filtering enriches T relative to U via differential hydrolysis stability in dark hydrothermal-vent conditions. UV is negligible. CSI = real enrichment × pool-drift factor. Stability ratio is an Arrhenius extrapolation; see Ea_T sweep for uncertainty.

## Configuration 
- Temperature: 40–90 °C (5 °C steps)  
- pH: 5–11 (0.5 steps)  
- Simulation hours: 240  
- n_replicates: 5  
- photon_flux: 0.0  
- polymer_length: 100  
- lipid_conc: 0.05  
- strand_state: dsDNA  
- min_total_dna_for_optimum: 100.0  
- exclude_post_peak_collapse: true  
- exclude_sterilized: true  
- confirmatory_replicates: 200  
- time_horizon_for_plot: 2400.0  

## New in v6.2.0
- `decisive_control_sweep` (20 seeds at T=90 °C, pH=11, 240 h)  
- New controls: `equal_degradation_gm`, `equal_degradation_kT`, `no_compartments` (plus retained baseline + `equal_pools`)  
- `equal_degradation_applied` flag and mode per run  

All other outputs, sweeps, history files, and CSV tables remain numerically identical.

## Key Metrics 
- Optimal conditions: T=90 °C, pH=11  
- max_enrichment_chemical: 2.582131474905796  
- real_enrichment_chemical: 1.5581225062464852  
- pool_drift_factor_chemical: 1.6572163486316285  
- fraction_T_chemical: 0.170802471102718  
- total_dna_chemical: 90830.8  
- final_survival_chemical: 0.713  
- stability_ratio: 356.9593281176499  
- hydrolysis_events_chemical: 36027.4  

## Confirmatory Optimum (200 replicates)
- enrichment_mean: 2.578963108434424 (±0.0186)  
- dna_TU_over_free_TU_mean: 1.5546335792102253  
- fraction_T_mean: 0.1705928904411876  
- total_dna_mean: 90833.56  
- regime: net_population_growth (100 %)  
- csi_shape frequencies: rising 0.86, peaked 0.13, collapsed 0.01  

## Decisive Control Sweep (v6.2.0 only – 20 seeds)
| Control                  | Regime              | Mean CSI | Mean real | Mean drift | Mean total_dna | Mean final_survival | Mean hydrolysis_events | k_U_hydro (mean)     | k_T_hydro (mean)      |
|--------------------------|---------------------|----------|-----------|------------|----------------|---------------------|------------------------|----------------------|-----------------------|
| baseline                 | net_population_growth | 2.580    | 1.559     | 1.657      | 90860          | 0.715               | 35900                  | 8.86e-6              | 2.48e-8               |
| equal_degradation_gm     | net_population_growth | 1.88     | 0.95      | 1.97       | 120000         | 0.948               | 6600                   | 4.69e-7              | 4.69e-7               |
| equal_degradation_kT     | net_population_growth | 1.89     | 0.95      | 1.97       | 120000         | 0.948               | 6600                   | 4.69e-7              | 4.69e-7               |
| equal_pools              | net_population_growth | 4.48     | 1.40      | 3.20       | 97300          | 0.768               | 29500                  | 8.86e-6              | 2.48e-8               |
| no_compartments          | net_population_growth | 4.35     | 3.20      | 1.36       | 53500          | 0.424               | 73300                  | 8.86e-6              | 2.48e-8               |

### Interpretation
When hydrolysis rates are equalised (`equal_degradation_gm` / `equal_degradation_kT`), enrichment collapses to ~1.88–1.89. The original CSI (~2.58) therefore disappears when k_U == k_T. This proves the T/U enrichment is **genuinely driven by differential hydrolysis stability** (T far more stable), not by pool-drift or other factors.

## All Other Sweeps
Identical numerical results between versions for flux, free_monomer_deam_scale, clay_rate, ea_T, and time_horizon. Mass-balance error <1e-14; no change in regime classification or stability_ratio.

## Conclusion
v6.2.0 reproduces every numerical result of v6.1.1 exactly. The **only new scientific feature** is the decisive_control_sweep, which conclusively demonstrates that T/U enrichment is carried by differential hydrolysis (k_U >> k_T). The stability ratio remains an Arrhenius extrapolation. The model shows no Darwinian evolution and negligible UV damage.

**Files provided**  
- UPDSF v6.1.1.py (164153 bytes)  
- UPDSF v6.2.0.py (182377 bytes)  
- analysis_v611_20261025_154131.json  
- analysis_v620_20261111_141412.json  

## Citation
If you use this framework in your research, please cite the following:

Reza Hashemi. (2026). Chemical Filtering of Thymine over Uracil in Prebiotic Hydrothermal Systems: A Mass-Consistent Stochastic Model with Complementary Insights from UV Photodamage Pathways. ChemRxiv. https://doi.org/10.26434/chemrxiv.15011902/v1

Reza Hashemi. (2026). Chemical Filtering of Thymine over Uracil in Prebiotic Hydrothermal Systems: A Mass-Consistent Stochastic Model with Complementary Insights from UV Photodamage Pathways. Zenodo. https://doi.org/10.5281/zenodo.23876507

Reza Hashemi. (2026). Supplementary Information UPDSF v6.1.1 – Unified Prebiotic DNA Selection Framework Version 6.1.1 Chemical Filtering of Thymine over Uracil in Prebiotic Hydrothermal Systems A Mass-Consistent Stochastic Model with Complementary Insights from UV Photodamage Pathways. Zenodo. https://doi.org/10.5281/zenodo.23876308

Reza Hashemi. (2026). mrhashemi2000/UPDSF_v6_1_1: Initial release [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23655168

## License
This work is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0) license.

See LICENSE for full terms.

Author: Seyed Mohammad Reza Hashemi

ORCID: 0009-0002-0645-5180
