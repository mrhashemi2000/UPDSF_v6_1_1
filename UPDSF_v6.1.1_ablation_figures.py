"""
Generate Figures S10-S12 (300 ppi) from UPDSF ablation results CSVs.

AUTHOR: Seyed Mohammad Reza Hashemi (Reza Hashemi) Intelligence-Augmented (IA)
VERSION: 6.1.1
License & Copyright
Copyright ©️ 2026 Seyed Mohammad Reza Hashemi
This work is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License (CC BY-NC-ND 4.0).

ORCID: 0009-0002-0645-5180

https://doi.org/10.5281/zenodo.23656782

https://doi.org/10.5281/zenodo.23655168

https://doi.org/10.5281/zenodo.23788249

https://doi.org/10.5281/zenodo.23788480

https://doi.org/10.5281/zenodo.23876507

https://doi.org/10.5281/zenodo.23876308

Inputs (must be in the same directory as this script):
  - updsf_ablation_results_paired_delta_vs_baseline.csv
  - updsf_ablation_results_per_run.csv
  - updsf_ablation_results_summary.csv

Outputs (300 ppi PNGs):
  - Figure_S10_paired_delta_vs_baseline.png
  - Figure_S11_scenario_distributions.png
  - Figure_S12_knockout_effects.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
PAIRED_CSV  = os.path.join(HERE, "updsf_ablation_results_paired_delta_vs_baseline.csv")
PER_RUN_CSV = os.path.join(HERE, "updsf_ablation_results_per_run.csv")
SUMMARY_CSV = os.path.join(HERE, "updsf_ablation_results_summary.csv")

DPI = 300

SCENARIO_ORDER = [
    "baseline",
    "no_hydrolysis",
    "no_C_deamination",
    "no_clay",
    "no_lipid",
    "no_polymerization_selection",
    "polymerization_selection_on",
]

# Legend order: put null-effect controls together at the end
LEGEND_ORDER = [
    "baseline",
    "no_hydrolysis",
    "no_C_deamination",
    "no_clay",
    "polymerization_selection_on",
    "no_lipid",
    "no_polymerization_selection",
]

SCENARIO_LABELS = {
    "baseline": "Baseline",
    "no_hydrolysis": "No hydrolysis",
    "no_C_deamination": "No C-deamination",
    "no_clay": "No clay",
    "no_lipid": "No lipid",
    "no_polymerization_selection": "No polymerization selection",
    "polymerization_selection_on": "Polymerization selection ON",
}

ABLATION_SCENARIOS = [
    "no_hydrolysis",
    "no_C_deamination",
    "no_clay",
    "no_lipid",
    "no_polymerization_selection",
    "polymerization_selection_on",
]

METRICS = [
    "total_dna_final",
    "T_fraction_total_DNA",
    "T_over_TplusU_DNA",
    "CSI_final",
    "CSI_peak",
]

METRIC_LABELS = {
    "total_dna_final": "Total DNA (final)",
    "T_fraction_total_DNA": "T fraction of total DNA",
    "T_over_TplusU_DNA": "T/(T+U) in DNA",
    "CSI_final": "CSI (final)",
    "CSI_peak": "CSI (peak)",
}

SCENARIO_COLORS = {
    "baseline": "#4d4d4d",
    "no_hydrolysis": "#1f77b4",
    "no_C_deamination": "#ff7f0e",
    "no_clay": "#2ca02c",
    "no_lipid": "#9467bd",
    "no_polymerization_selection": "#8c564b",
    "polymerization_selection_on": "#d62728",
}
SCENARIO_MARKERS = {
    "baseline": "o",
    "no_hydrolysis": "s",
    "no_C_deamination": "^",
    "no_clay": "D",
    "no_lipid": "v",
    "no_polymerization_selection": "P",
    "polymerization_selection_on": "X",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7.5,
    "axes.linewidth": 0.8,
    "savefig.dpi": DPI,
    "figure.dpi": DPI,
})


# ----------------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------------
paired  = pd.read_csv(PAIRED_CSV)
per_run = pd.read_csv(PER_RUN_CSV)
summary = pd.read_csv(SUMMARY_CSV)

paired = paired[paired["scenario"].isin(ABLATION_SCENARIOS)].copy()
paired["scenario"] = pd.Categorical(paired["scenario"], categories=ABLATION_SCENARIOS, ordered=True)

per_run["scenario"] = pd.Categorical(per_run["scenario"], categories=SCENARIO_ORDER, ordered=True)
summary["scenario"] = pd.Categorical(summary["scenario"], categories=SCENARIO_ORDER, ordered=True)


def _short_num(x):
    """Compact numeric formatter for axis ticks."""
    if x == 0:
        return "0"
    ax = abs(x)
    if ax >= 1e6:
        return f"{x/1e6:.1f}M"
    if ax >= 1e3:
        return f"{x/1e3:.0f}k"
    if ax >= 1:
        return f"{x:.1f}"
    return f"{x:.2g}"


def _save_png(fig, name):
    out = os.path.join(HERE, name)
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out}")


# ----------------------------------------------------------------------------
# Figure S10: Paired delta vs baseline (forest-style per metric)
# ----------------------------------------------------------------------------
def make_figure_S10():
    n_metrics = len(METRICS)
    fig, axes = plt.subplots(n_metrics, 1, figsize=(6.5, 2.1 * n_metrics), sharex=False)

    for ax, metric in zip(axes, METRICS):
        sub = paired[paired["metric"] == metric].sort_values("scenario")
        sub = sub.set_index("scenario").reindex(ABLATION_SCENARIOS).dropna(subset=["paired_delta_mean"])

        scenarios = list(sub.index)
        y = np.arange(len(scenarios))[::-1]
        means = sub["paired_delta_mean"].values
        sds   = sub["paired_delta_sd"].values
        labels = [SCENARIO_LABELS[s] for s in scenarios]

        # Forest-plot style: horizontal SD segment with vertical end caps.
        # Visible at any axis scale (unlike errorbar whiskers, which were
        # sub-pixel for these very tight SDs).
        for yi, m, sd, scen in zip(y, means, sds, scenarios):
            color = SCENARIO_COLORS[scen]
            if sd > 0:
                ax.hlines(yi, m - sd, m + sd, colors=color, lw=2.4, zorder=2)
                ax.vlines([m - sd, m + sd], yi - 0.11, yi + 0.11,
                          colors=color, lw=1.0, zorder=2)
            ax.scatter(m, yi, s=34, c=color, marker=SCENARIO_MARKERS[scen],
                       zorder=3, edgecolors="black", linewidths=0.45)

        ax.axvline(0.0, color="black", lw=0.8, ls="--", zorder=1)
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_ylim(-0.6, len(scenarios) - 0.4)
        ax.set_title(METRIC_LABELS[metric], loc="left", fontweight="bold")
        ax.grid(axis="x", ls=":", lw=0.5, alpha=0.6)
        ax.tick_params(axis="y", length=0)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.xaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: _short_num(v))
        )

    axes[-1].set_xlabel("Paired Δ vs baseline (mean ± SD, 5 seeds)")
    fig.suptitle(
        "Figure S10. Paired effect of each model ablation relative to baseline",
        fontsize=11, fontweight="bold", y=0.995,
    )
    fig.text(
        0.5, 0.005,
        "No lipid and no polymerization selection are bit-identical to baseline "
        "(paired Δ = 0, SD = 0).",
        ha="center", fontsize=7.0, style="italic",
    )
    fig.tight_layout(rect=[0, 0.02, 1, 0.98])
    _save_png(fig, "Figure_S10_paired_delta_vs_baseline.png")


# ----------------------------------------------------------------------------
# Figure S11: Per-seed distributions across scenarios
# ----------------------------------------------------------------------------
def make_figure_S11():
    plot_metrics = ["total_dna_final", "T_fraction_total_DNA", "T_over_TplusU_DNA", "CSI_final"]
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.8))
    axes = axes.ravel()

    for ax, metric in zip(axes, plot_metrics):
        x_positions = np.arange(len(SCENARIO_ORDER))
        for i, scen in enumerate(SCENARIO_ORDER):
            vals = per_run.loc[per_run["scenario"] == scen, metric].values
            if len(vals) == 0:
                continue
            jitter = np.linspace(-0.12, 0.12, len(vals)) if len(vals) > 1 else np.array([0.0])
            ax.scatter(
                np.full(len(vals), x_positions[i]) + jitter, vals,
                s=18, color=SCENARIO_COLORS[scen], marker=SCENARIO_MARKERS[scen],
                alpha=0.85, edgecolors="black", linewidths=0.3, zorder=2,
            )
            m = vals.mean()
            sd = vals.std(ddof=1) if len(vals) > 1 else 0.0
            if sd > 0:
                ax.errorbar(
                    x_positions[i], m, yerr=sd, fmt="none",
                    ecolor="black", elinewidth=1.0, capsize=3, capthick=0.9, zorder=3,
                )
            ax.scatter(x_positions[i], m, s=44, color="white",
                       edgecolors="black", linewidths=0.8, zorder=4)

        # Symmetric 0–0.8 y-ticks on the T-fraction panel for readability
        if metric == "T_fraction_total_DNA":
            ax.set_ylim(0.0, 0.8)

        ax.set_xticks(x_positions)
        ax.set_xticklabels([SCENARIO_LABELS[s] for s in SCENARIO_ORDER], rotation=35, ha="right")
        ax.set_title(METRIC_LABELS[metric], loc="left", fontweight="bold")
        ax.grid(axis="y", ls=":", lw=0.5, alpha=0.6)
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: _short_num(v)))
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    handles = [
        Line2D([0], [0], marker=SCENARIO_MARKERS[s], color="none",
               markerfacecolor=SCENARIO_COLORS[s], markeredgecolor="black",
               markersize=6, label=SCENARIO_LABELS[s])
        for s in LEGEND_ORDER
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(
        "Figure S11. Per-seed distributions of key output metrics across scenarios",
        fontsize=11, fontweight="bold", y=1.00,
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.98])
    _save_png(fig, "Figure_S11_scenario_distributions.png")


# ----------------------------------------------------------------------------
# Figure S12: Summary mean ± SD per metric
# ----------------------------------------------------------------------------
def make_figure_S12():
    n_metrics = len(METRICS)
    fig, axes = plt.subplots(n_metrics, 1, figsize=(7.0, 1.9 * n_metrics), sharex=False)

    for ax, metric in zip(axes, METRICS):
        sub = summary.copy()
        mean_col = f"{metric}_mean"
        std_col  = f"{metric}_std"
        sub = sub.dropna(subset=[mean_col])

        baseline_val = sub.loc[sub["scenario"] == "baseline", mean_col]
        baseline_val = float(baseline_val.iloc[0]) if len(baseline_val) else np.nan

        show = sub[sub["scenario"] != "baseline"].copy()
        show = show.sort_values(mean_col, ascending=True)

        y = np.arange(len(show))
        vals = show[mean_col].values
        sds  = show[std_col].fillna(0.0).values
        colors = [SCENARIO_COLORS[s] for s in show["scenario"].astype(str)]
        labels = [SCENARIO_LABELS[s] for s in show["scenario"].astype(str)]

        ax.barh(y, vals, xerr=sds, color=colors, alpha=0.85,
                edgecolor="black", linewidth=0.5, capsize=3,
                error_kw={"elinewidth": 0.9})
        if np.isfinite(baseline_val):
            ax.axvline(baseline_val, color="black", ls="--", lw=1.0, zorder=1)

        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_title(METRIC_LABELS[metric], loc="left", fontweight="bold")
        ax.grid(axis="x", ls=":", lw=0.5, alpha=0.6)
        ax.tick_params(axis="y", length=0)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: _short_num(v)))

    axes[-1].set_xlabel("Mean across 5 seeds (error bars = SD)")

    # Figure-level legend: dashed line = baseline mean
    fig.legend(
        handles=[Line2D([0], [0], color="black", ls="--", lw=1.0,
                        label="Baseline mean")],
        loc="upper right", frameon=False, bbox_to_anchor=(0.99, 0.985),
    )
    fig.suptitle(
        "Figure S12. Summary of scenario-level outputs vs baseline",
        fontsize=11, fontweight="bold", y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    _save_png(fig, "Figure_S12_knockout_effects.png")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    make_figure_S10()
    make_figure_S11()
    make_figure_S12()
    print("\nAll figures generated at 300 ppi.")
