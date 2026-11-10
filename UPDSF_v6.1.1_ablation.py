#!/usr/bin/env python3
"""
UPDSF v6.1.1 — mechanism-ablation experiment
Run locally beside "UPDSF v6.1.1.py".

AUTHOR: Seyed Mohammad Reza Hashemi (Reza Hashemi) Intelligence-Augmented (IA)
VERSION: 6.1.1
License & Copyright
Copyright ©️ 2026 Seyed Mohammad Reza Hashemi
This work is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License (CC BY-NC-ND 4.0).

ORCID: 0009-0002-0645-5180

https://doi.org/10.26434/chemrxiv.15011902/v1

https://doi.org/10.5281/zenodo.23656782

https://doi.org/10.5281/zenodo.23655168

https://doi.org/10.5281/zenodo.23788249

https://doi.org/10.5281/zenodo.23788480

https://doi.org/10.5281/zenodo.23876507

https://doi.org/10.5281/zenodo.23876308

Default: 5 paired random seeds, 2400 h, dt=60 s.
Writes per-run results and mean/std summaries to CSV.

IMPORTANT:
- This is a model sensitivity experiment, not experimental validation.
- Runs are sequential because EmpiricalRates parameters are temporarily modified.
- The source model uses global NumPy/Python RNGs; each engine constructor resets them.
"""
from __future__ import annotations

import argparse
import importlib.util
import math
from pathlib import Path
import sys
import pandas as pd


def load_model(model_path: Path):
    spec = importlib.util.spec_from_file_location("updsf_v611", model_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load model file: {model_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def set_rate_value(rate_obj, value):
    """Temporarily set Source.value while preserving the Source object."""
    old = rate_obj.value
    rate_obj.value = value
    return old


def run_one(m, scenario: str, seed: int, horizon_h: float, dt_s: float):
    ER = m.EmpiricalRates
    engine = None

    # Save any class-level value changed for this scenario; always restore it.
    old_clay_ads = ER.CLAY_ADSORPTION_RATE.value
    try:
        # The reference condition matches the vent plot as closely as possible.
        # Polymerization-selection flags are explicitly OFF in the reference run.
        engine = m.VSSUFEngine(
            temperature_C=90.0,
            pH_in_situ=11.0,
            seed=seed,
            max_time_hours=horizon_h,
            verbose=False,
            photon_flux=0.0,
            polymer_length=100,
            lipid_conc=0.05,
            lipid_chain_length=18,
            ionic_strength_M=0.1,
            Ca_mM=1.0,
            record_interval_seconds=3600.0,  # hourly history is enough for summary
            dt=dt_s,
            strand_state="dsDNA",
            use_incorporation_bias=False,
            use_fidelity_advantage=False,
            equalize_ut_pools=False,
            disable_compartments=False,
        )

        # Scenario-specific intervention. Only one mechanism is changed at a time.
        if scenario == "no_hydrolysis":
            engine.k_hydro = {b: 0.0 for b in engine.k_hydro}
            engine.k_hydro_clay = {b: 0.0 for b in engine.k_hydro_clay}
            engine.k_hydro_lipid = {b: 0.0 for b in engine.k_hydro_lipid}

        elif scenario == "no_C_deamination":
            engine.k_deam_C = 0.0
            engine.k_free_deam_C = 0.0

        elif scenario == "no_clay":
            # Remove clay adsorption and clay's hydrolysis-protection effect,
            # while leaving lipid behavior active.
            set_rate_value(ER.CLAY_ADSORPTION_RATE, 0.0)
            engine.clay_protection = 1.0
            engine.k_hydro_clay = dict(engine.k_hydro)
            # Ensure no clay-bound population remains in this scenario.
            for b in "UTCAGX":
                engine.species[f"dsDNA_{b}_clay"] = 0

        elif scenario == "no_lipid":
            # Remove lipid entry/retention and lipid's protection effect,
            # while leaving clay behavior active.
            engine.lipid_membrane.protection_active = False
            engine.lipid_protection_factor = 1.0
            engine.k_hydro_lipid = dict(engine.k_hydro)
            for b in "UTCAGX":
                engine.species[f"dsDNA_{b}_lipid"] = 0

        elif scenario == "no_polymerization_selection":
            # Already the reference model's default: no incorporation bias and
            # no fidelity advantage. Kept as an explicit audit scenario.
            engine.use_incorporation_bias = False
            engine.use_fidelity_advantage = False

        elif scenario == "polymerization_selection_on":
            # Positive-control sensitivity test, not a deletion from the baseline.
            engine.use_incorporation_bias = True
            engine.use_fidelity_advantage = True

        elif scenario != "baseline":
            raise ValueError(f"Unknown scenario: {scenario}")

        engine.run(max_time=horizon_h * 3600.0)

        final_counts = {
            b: int(engine.species[f"dsDNA_{b}"]
                   + engine.species[f"dsDNA_{b}_clay"]
                   + engine.species[f"dsDNA_{b}_lipid"])
            for b in "UTCAGX"
        }
        total_dna = sum(final_counts.values())
        t_fraction_total = final_counts["T"] / max(1, total_dna)
        tu_fraction = final_counts["T"] / max(1, final_counts["T"] + final_counts["U"])

        return {
            "scenario": scenario,
            "seed": seed,
            "horizon_h": horizon_h,
            "dt_s": dt_s,
            "total_dna_final": total_dna,
            "T_final": final_counts["T"],
            "U_final": final_counts["U"],
            "C_final": final_counts["C"],
            "A_final": final_counts["A"],
            "G_final": final_counts["G"],
            "X_final": final_counts["X"],
            "T_fraction_total_DNA": t_fraction_total,
            "T_over_TplusU_DNA": tu_fraction,
            "CSI_final": engine.get_chemical_selection_index(),
            "CSI_peak": engine.peak_csi(),
            "CSI_peak_time_h": engine.time_to_peak_csi_h(),
            "hydrolysis_events": engine.hydrolysis_events,
            "C_deamination_events": engine.deamination_C_events,
            "clay_adsorption_events": engine.clay_adsorption_events,
            "lipid_entry_events": engine.lipid_entry_events,
            "mass_balance_max_abs_rel_error": engine.max_abs_rel_error(),
            "consumed_identity_max_abs_residual": engine.max_abs_consumed_residual(),
        }
    finally:
        ER.CLAY_ADSORPTION_RATE.value = old_clay_ads


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model", default="UPDSF v6.1.1.py",
        help="Path to the original UPDSF Python file"
    )
    parser.add_argument("--horizon-h", type=float, default=2400.0)
    parser.add_argument("--dt-s", type=float, default=60.0)
    parser.add_argument("--seeds", type=int, default=5,
                        help="Number of paired random seeds (default: 5)")
    parser.add_argument("--seed-start", type=int, default=42)
    parser.add_argument("--out", default="updsf_ablation_results")
    args = parser.parse_args()

    model_path = Path(args.model).expanduser().resolve()
    if not model_path.exists():
        raise FileNotFoundError(
            f"Model file not found: {model_path}\n"
            "Put this script beside 'UPDSF v6.1.1.py' or pass --model."
        )
    if args.horizon_h <= 0 or args.dt_s <= 0 or args.seeds < 1:
        raise ValueError("horizon-h, dt-s and seeds must be positive.")

    m = load_model(model_path)
    scenarios = [
        "baseline",
        "no_hydrolysis",
        "no_C_deamination",
        "no_clay",
        "no_lipid",
        "no_polymerization_selection",
        "polymerization_selection_on",
    ]
    seeds = list(range(args.seed_start, args.seed_start + args.seeds))
    rows = []
    total = len(scenarios) * len(seeds)
    run_index = 0

    for scenario in scenarios:
        for seed in seeds:
            run_index += 1
            print(
                f"[{run_index}/{total}] scenario={scenario}, seed={seed}, "
                f"horizon={args.horizon_h:g} h",
                flush=True,
            )
            rows.append(run_one(m, scenario, seed, args.horizon_h, args.dt_s))

    raw = pd.DataFrame(rows)
    out_prefix = Path(args.out)
    raw_path = out_prefix.with_name(out_prefix.name + "_per_run.csv")
    summary_path = out_prefix.with_name(out_prefix.name + "_summary.csv")
    delta_path = out_prefix.with_name(out_prefix.name + "_paired_delta_vs_baseline.csv")

    raw.to_csv(raw_path, index=False)

    numeric = [
        "total_dna_final", "T_fraction_total_DNA", "T_over_TplusU_DNA",
        "CSI_final", "CSI_peak", "CSI_peak_time_h", "hydrolysis_events",
        "C_deamination_events", "clay_adsorption_events", "lipid_entry_events",
        "mass_balance_max_abs_rel_error", "consumed_identity_max_abs_residual",
    ]
    summary = raw.groupby("scenario")[numeric].agg(["mean", "std"])
    summary.columns = [f"{metric}_{stat}" for metric, stat in summary.columns]
    summary.reset_index().to_csv(summary_path, index=False)

    # Paired seed differences reduce noise when comparing scenarios.
    baseline = raw[raw.scenario == "baseline"].set_index("seed")
    delta_rows = []
    for scenario in scenarios:
        if scenario == "baseline":
            continue
        subset = raw[raw.scenario == scenario].set_index("seed")
        for metric in [
            "total_dna_final", "T_fraction_total_DNA", "T_over_TplusU_DNA",
            "CSI_final", "CSI_peak",
        ]:
            diffs = subset[metric] - baseline[metric]
            delta_rows.append({
                "scenario": scenario,
                "metric": metric,
                "paired_delta_mean": float(diffs.mean()),
                "paired_delta_sd": float(diffs.std(ddof=1)) if len(diffs) > 1 else float("nan"),
                "n_seeds": len(diffs),
            })
    pd.DataFrame(delta_rows).to_csv(delta_path, index=False)

    print("\nFinished.")
    print(f"Per-run results: {raw_path.resolve()}")
    print(f"Summary:         {summary_path.resolve()}")
    print(f"Paired deltas:   {delta_path.resolve()}")
    print("\nInterpretation note: 'no_polymerization_selection' should match baseline")
    print("because the reference flags are already OFF. The 'selection_on' run is")
    print("a positive-control sensitivity test, not a one-mechanism ablation.")
    print("Clay and lipid are separated by independently disabling adsorption/protection")
    print("and lipid entry/protection, respectively.")


if __name__ == "__main__":
    main()
