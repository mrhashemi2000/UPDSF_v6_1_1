[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23788249-blue)](https://doi.org/10.5281/zenodo.23788249) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23656782-blue)](https://doi.org/10.5281/zenodo.23656782) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23655168-blue)](https://doi.org/10.5281/zenodo.23655168) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23788480-blue)](https://doi.org/10.5281/zenodo.23788480) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

# Exact Numerical Reconciliation
```text
k_U / k_T ≈ 357
```

The implementation uses the following parameters:

- `k_U_ref = 1.3e-10 s^-1`
- `k_T_ref = 1.1e-13 s^-1`
- `Ea_U = 27 kcal/mol`
- `Ea_T = 32 kcal/mol`
- `BASE_CATALYSIS_FACTOR_U = 1.0`
- `BASE_CATALYSIS_FACTOR_T = 0.75`

It also applies different pKa values for uracil and thymine when calculating the pH-dependent factors.

## Complete rate ratio

The complete ratio is calculated as:

```text
k_U / k_T
    = (k_U,Arr / k_T,Arr)
      * (f_pH,U / f_pH,T)
      * (BASE_CATALYSIS_FACTOR_U / BASE_CATALYSIS_FACTOR_T)
```

Using the implementation's values:

```text
k_U / k_T
    ≈ 260.95 * 1.026 * (1.0 / 0.75)
    ≈ 357
```

Therefore:

```text
k_U / k_T ≈ 357
```

The implementation's result is mathematically consistent with its full set of model factors. However, numerical consistency does not independently establish the empirical validity of the assumed catalysis factors, the extrapolated activation energies, or the underlying reference rates.
