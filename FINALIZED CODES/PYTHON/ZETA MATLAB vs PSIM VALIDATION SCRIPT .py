"""
Compares the ZETA MATLAB model against PSIM measurements over the nine
operating points. Reads the MATLAB summary and component-loss CSVs plus
the PSIM summary / loss / validation CSVs, builds a case-merged comparison
table, writes two validation CSVs, and emits four figures.

Inputs (all CSV, paths below):
    MATLAB  ZETA_v10_7_summary.csv, ZETA_v10_7_component_losses.csv
    PSIM    ZETA_summary.csv, ZETA_loss_breakdown.csv, ZETA_validation.csv

Outputs, in validation_output/:
    ZETA_summary_validation.csv
    ZETA_switching_loss_validation.csv
    Fig_V1..V4 as vector PDF and 300 dpi PNG

Efficiency is reported two ways: stage-vs-stage (no gate drive on either
side) and system-vs-system (gate drive folded in on both sides). The
system comparison is skipped automatically if the PSIM summary doesn't
carry Efficiency_System_Optional_pct.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


# =============================================================================
# Paths
# =============================================================================

MATLAB_DIR      = Path(r"C:\Users\User\Desktop\ZETA MATLAB")
PSIM_OUTPUT_DIR = Path(r"C:\Users\User\Desktop\ZETA NEW\ZETA_analysis_output")

matlab_summary_path  = MATLAB_DIR / "ZETA_v10_7_summary.csv"
matlab_losses_path   = MATLAB_DIR / "ZETA_v10_7_component_losses.csv"
psim_summary_path    = PSIM_OUTPUT_DIR / "ZETA_summary.csv"
psim_losses_path     = PSIM_OUTPUT_DIR / "ZETA_loss_breakdown.csv"
psim_validation_path = PSIM_OUTPUT_DIR / "ZETA_validation.csv"

OUT_DIR = MATLAB_DIR / "validation_output"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Helpers
# =============================================================================

def require_columns(df: pd.DataFrame, cols: List[str], source: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(
            f"{source} is missing required columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )


def add_case_column_if_missing(df, *, vin_col, io_col, source_name):
    """Case numbering is implicit; assign by (Vin, Io) sort if not present."""
    df = df.copy()
    if "Case" in df.columns:
        return df
    if len(df) != 9:
        raise ValueError(f"{source_name} has {len(df)} rows. Expected 9.")
    df = df.sort_values([vin_col, io_col]).reset_index(drop=True)
    df.insert(0, "Case", range(1, len(df) + 1))
    return df


def pct_error(est, ref):
    est = np.asarray(est, dtype=float)
    ref = np.asarray(ref, dtype=float)
    return np.where(np.abs(ref) > 1e-12, 100.0 * (est - ref) / ref, np.nan)


def save_fig(fig, name):
    fig.savefig(OUT_DIR / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name} (.pdf, .png)")


# =============================================================================
# Load
# =============================================================================

print("=" * 95)
print("MATLAB vs PSIM validation — ZETA")
print("=" * 95)

try:
    matlab_summary  = pd.read_csv(matlab_summary_path)
    matlab_losses   = pd.read_csv(matlab_losses_path)
    psim_summary    = pd.read_csv(psim_summary_path)
    psim_losses     = pd.read_csv(psim_losses_path)
    psim_validation = pd.read_csv(psim_validation_path)
    print("All CSVs loaded.")
except FileNotFoundError as exc:
    raise FileNotFoundError(f"Required CSV not found:\n{exc}") from exc
except pd.errors.EmptyDataError as exc:
    raise ValueError(f"A CSV file is empty:\n{exc}") from exc
except pd.errors.ParserError as exc:
    raise ValueError(f"A CSV file could not be parsed:\n{exc}") from exc


# =============================================================================
# Column sanity check
# =============================================================================

require_columns(
    matlab_summary,
    ["Vin_V", "Io_A", "Duty", "StageEfficiency_pct",
     "SystemEfficiency_pct", "StageLoss_W"],
    "MATLAB summary CSV",
)
require_columns(
    matlab_losses,
    ["Vin", "Io", "Q1_cond_W",
     "Q1_sw_dpt_W", "Q1_sw_transition_W", "Q1_sw_miller_W"],
    "MATLAB component-loss CSV",
)
require_columns(
    psim_summary,
    ["Case", "Duty_Measured", "Efficiency_Stage_pct",
     "Extracted_Stage_Loss_Sum_W", "Closure_Error_pct"],
    "PSIM summary CSV",
)
require_columns(
    psim_losses,
    ["Case", "MOSFET_Loss_W",
     "P_sw_zeta_W", "E_on_zeta_uJ", "E_off_zeta_uJ",
     "I_sw_turnoff_avg_A", "I_sw_turnon_avg_A", "t_sw_effective_ns"],
    "PSIM loss-breakdown CSV",
)
require_columns(
    psim_validation,
    ["Case", "CCM_Overall_Status", "Closure_Status", "Overall_Status"],
    "PSIM validation CSV",
)

has_psim_sys_eff = "Efficiency_System_Optional_pct" in psim_summary.columns
if not has_psim_sys_eff:
    print("WARNING: PSIM summary lacks 'Efficiency_System_Optional_pct'.")
    print("         System-vs-system efficiency comparison will be skipped.")

print("Column check passed for all five CSVs.")


# =============================================================================
# Case numbers
# =============================================================================

matlab_summary = add_case_column_if_missing(
    matlab_summary, vin_col="Vin_V", io_col="Io_A",
    source_name="MATLAB summary CSV",
)
matlab_losses = add_case_column_if_missing(
    matlab_losses, vin_col="Vin", io_col="Io",
    source_name="MATLAB component-loss CSV",
)


# =============================================================================
# Merge summary
# =============================================================================

merged = pd.merge(
    matlab_summary, psim_summary,
    on="Case", suffixes=("_matlab", "_psim"),
    validate="one_to_one",
)
if len(merged) != 9:
    raise ValueError(f"Merged summary has {len(merged)} rows. Expected 9.")
print(f"Merged summary: {len(merged)} cases.")


# =============================================================================
# Summary differences and pass/fail flags
# =============================================================================

merged["Duty_diff"] = abs(merged["Duty"] - merged["Duty_Measured"])
merged["Eff_diff"]  = abs(merged["StageEfficiency_pct"] - merged["Efficiency_Stage_pct"])
merged["Loss_diff"] = abs(merged["StageLoss_W"] - merged["Extracted_Stage_Loss_Sum_W"])
merged["Loss_diff_rel_pct"] = np.where(
    abs(merged["StageLoss_W"]) > 1e-12,
    100.0 * merged["Loss_diff"] / abs(merged["StageLoss_W"]),
    np.nan,
)

# Thresholds are descriptive, not enforcement. The report interprets the
# numbers; these flags just make it easy to scan the printed table.
merged["Duty_status"]    = np.where(merged["Duty_diff"]         <= 0.01,  "PASS", "FAIL")
merged["Eff_status"]     = np.where(merged["Eff_diff"]          <= 4.5,   "PASS", "FAIL")
merged["Loss_status"]    = np.where(merged["Loss_diff_rel_pct"] <= 35.0,  "PASS", "FAIL")
merged["Closure_status"] = np.where(abs(merged["Closure_Error_pct"]) <= 3.0, "PASS", "FAIL")


# =============================================================================
# Merge PSIM validation status
# =============================================================================

merged = pd.merge(
    merged,
    psim_validation[["Case", "CCM_Overall_Status", "Closure_Status", "Overall_Status"]],
    on="Case", suffixes=("", "_psim_val"), validate="one_to_one",
)
merged.rename(
    columns={
        "CCM_Overall_Status": "PSIM_CCM_Overall",
        "Closure_Status":     "PSIM_Closure_Status",
        "Overall_Status":     "PSIM_Overall_Status",
    },
    inplace=True,
)


# =============================================================================
# Console summary
# =============================================================================

print("\n" + "=" * 95)
print("SUMMARY VALIDATION — ZETA")
print("=" * 95)

print(merged[[
    "Case",
    "Duty_diff",         "Duty_status",
    "Eff_diff",          "Eff_status",
    "Loss_diff_rel_pct", "Loss_status",
    "Closure_Error_pct", "Closure_status",
    "PSIM_CCM_Overall",
    "PSIM_Overall_Status",
]].to_string(index=False))

summary_out = OUT_DIR / "ZETA_summary_validation.csv"
merged.to_csv(summary_out, index=False)
print(f"\nSummary validation saved to: {summary_out}")


# =============================================================================
# Switching-loss comparison
# =============================================================================

matlab_mosfet = matlab_losses[[
    "Case", "Q1_cond_W",
    "Q1_sw_dpt_W", "Q1_sw_transition_W", "Q1_sw_miller_W",
]].copy()

psim_mosfet = psim_losses[[
    "Case", "MOSFET_Loss_W",
    "P_sw_zeta_W", "E_on_zeta_uJ", "E_off_zeta_uJ",
    "I_sw_turnoff_avg_A", "I_sw_turnon_avg_A", "t_sw_effective_ns",
]].copy()

cmp = pd.merge(matlab_mosfet, psim_mosfet, on="Case", validate="one_to_one")

cmp = pd.merge(
    cmp,
    merged[["Case", "Vin_Commanded_V", "Io_A"]].rename(
        columns={"Vin_Commanded_V": "Vin_V"}),
    on="Case", validate="one_to_one",
)

# Per-case percentage errors of each MATLAB method against the PSIM
# measured P_sw. Positive = MATLAB overpredicts.
cmp["DPT_err_pct"]   = pct_error(cmp["Q1_sw_dpt_W"],        cmp["P_sw_zeta_W"])
cmp["Trans_err_pct"] = pct_error(cmp["Q1_sw_transition_W"], cmp["P_sw_zeta_W"])
cmp["Mill_err_pct"]  = pct_error(cmp["Q1_sw_miller_W"],     cmp["P_sw_zeta_W"])

# Efficiency: both stage-vs-stage and system-vs-system. The system diff is
# only meaningful when PSIM actually exported a system efficiency.
cmp["Eff_stage_diff_pp"] = merged["StageEfficiency_pct"] - merged["Efficiency_Stage_pct"]
if has_psim_sys_eff:
    cmp["Eff_sys_diff_pp"] = (
        merged["SystemEfficiency_pct"] - merged["Efficiency_System_Optional_pct"]
    )
else:
    cmp["Eff_sys_diff_pp"] = np.nan

mean_dpt   = np.nanmean(np.abs(cmp["DPT_err_pct"]))
max_dpt    = np.nanmax(np.abs(cmp["DPT_err_pct"]))
mean_trans = np.nanmean(np.abs(cmp["Trans_err_pct"]))
max_trans  = np.nanmax(np.abs(cmp["Trans_err_pct"]))
mean_mill  = np.nanmean(np.abs(cmp["Mill_err_pct"]))
max_mill   = np.nanmax(np.abs(cmp["Mill_err_pct"]))


# =============================================================================
# Printed switching-loss table
# =============================================================================

print("\n" + "=" * 95)
print("SWITCHING LOSS — three MATLAB methods vs PSIM measured")
print("=" * 95)

hdr = (f"{'Case':>4} | {'Vin':>4} {'Io':>5} | "
       f"{'DPT_W':>8} {'Trans_W':>8} {'Mill_W':>8} | "
       f"{'PSIM_W':>8} | "
       f"{'DPT_%':>8} {'Trans_%':>9} {'Mill_%':>8}")
print(hdr)
print("-" * len(hdr))

for _, r in cmp.iterrows():
    print(f"{int(r['Case']):>4} | "
          f"{r['Vin_V']:>4.0f} {r['Io_A']:>5.1f} | "
          f"{r['Q1_sw_dpt_W']:>8.4f} {r['Q1_sw_transition_W']:>8.4f} "
          f"{r['Q1_sw_miller_W']:>8.4f} | "
          f"{r['P_sw_zeta_W']:>8.4f} | "
          f"{r['DPT_err_pct']:>+8.1f} {r['Trans_err_pct']:>+9.1f} "
          f"{r['Mill_err_pct']:>+8.1f}")

print("-" * len(hdr))
print(f"\nSUMMARY (absolute error vs PSIM measured):")
print(f"  DPT          : mean = {mean_dpt:5.1f}%,  max = {max_dpt:5.1f}%")
print(f"  Transition   : mean = {mean_trans:5.1f}%,  max = {max_trans:5.1f}%")
print(f"  Miller       : mean = {mean_mill:5.1f}%,  max = {max_mill:5.1f}%")

cmp_out = OUT_DIR / "ZETA_switching_loss_validation.csv"
cmp.to_csv(cmp_out, index=False)
print(f"\nSwitching-loss comparison saved to: {cmp_out}")


# =============================================================================
# Figures
# =============================================================================

plt.rcParams.update({
    'font.family':    'serif',
    'font.serif':     ['Times New Roman', 'DejaVu Serif'],
    'font.size':      11,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'legend.fontsize': 10,
    'axes.linewidth': 1.0,
    'lines.linewidth': 1.8,
    'lines.markersize': 9,
    'figure.dpi':     120,
})

CLR_DPT   = [0.00, 0.45, 0.70]
CLR_TRANS = [0.84, 0.37, 0.00]
CLR_MILL  = [0.17, 0.63, 0.17]
CLR_PSIM  = [0.20, 0.20, 0.20]

case_idx    = cmp["Case"].astype(int).values
case_labels = [f"{int(v)}V/{io:.1f}A"
               for v, io in zip(cmp["Vin_V"], cmp["Io_A"])]

print("\nGenerating figures...")


# ---- V1: switching loss vs case ----
fig1, ax = plt.subplots(figsize=(9.5, 5.5))
ax.plot(case_idx, cmp["Q1_sw_dpt_W"],        'o-', color=CLR_DPT,
        label='MATLAB DPT (§3.4)',              zorder=4)
ax.plot(case_idx, cmp["Q1_sw_transition_W"], 's-', color=CLR_TRANS,
        label='MATLAB transition (§3.2 gen-1)',  zorder=3)
ax.plot(case_idx, cmp["Q1_sw_miller_W"],     '^-', color=CLR_MILL,
        label='MATLAB miller (§3.2 gen-2)',     zorder=3)
ax.plot(case_idx, cmp["P_sw_zeta_W"],        'D-', color=CLR_PSIM,
        markerfacecolor='white', markeredgewidth=1.8,
        label='PSIM measured',   zorder=5)
ax.set_xticks(case_idx)
ax.set_xticklabels(case_labels, rotation=35, ha='right')
ax.set_xlabel('Operating point  (V$_{in}$ / I$_o$)')
ax.set_ylabel('Switching loss  $P_{sw}$  (W)')
ax.set_title(f'ZETA switching-loss prediction vs PSIM measurement\n'
             f'(DPT mean error = {mean_dpt:.1f}%,  max = {max_dpt:.1f}%)')
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc='upper left', frameon=True)
ax.set_xlim(0.5, 9.5)
save_fig(fig1, "Fig_V1_ZETA_SwitchingLoss_vs_Case")


# ---- V2: absolute error by method, log scale ----
fig2, ax = plt.subplots(figsize=(10, 5.5))
x     = np.arange(len(cmp))
width = 0.27

ax.bar(x - width, np.abs(cmp["DPT_err_pct"]),   width, color=CLR_DPT,
       edgecolor='black', linewidth=0.5,
       label=f'DPT (mean {mean_dpt:.1f}%)')
ax.bar(x,         np.abs(cmp["Trans_err_pct"]), width, color=CLR_TRANS,
       edgecolor='black', linewidth=0.5,
       label=f'Transition (mean {mean_trans:.0f}%)')
ax.bar(x + width, np.abs(cmp["Mill_err_pct"]),  width, color=CLR_MILL,
       edgecolor='black', linewidth=0.5,
       label=f'Miller (mean {mean_mill:.0f}%)')

# Only the DPT bars are labelled; the other two would clutter at this width
for xi, ei in zip(x - width, np.abs(cmp["DPT_err_pct"])):
    ax.text(xi, ei + 5, f"{ei:.0f}", ha='center', va='bottom',
            fontsize=8, color=CLR_DPT, fontweight='bold')

ax.set_xticks(x)
ax.set_xticklabels(case_labels, rotation=35, ha='right')
ax.set_xlabel('Operating point  (V$_{in}$ / I$_o$)')
ax.set_ylabel('Absolute error vs PSIM  (%)')
ax.set_title('ZETA switching-loss error by method')
ax.set_yscale('log')
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.set_ylim(1, max(np.nanmax(np.abs(cmp["Trans_err_pct"])),
                   np.nanmax(np.abs(cmp["Mill_err_pct"]))) * 1.6)
ax.axhline(25, color='gray', linestyle='--', linewidth=1.0, alpha=0.7)
ax.text(9.4, 25, '25% tolerance', ha='right', va='bottom',
        fontsize=9, color='gray', style='italic')
ax.grid(True, axis='y', linestyle=':', alpha=0.5)
ax.legend(loc='upper left', frameon=True)
save_fig(fig2, "Fig_V2_ZETA_Error_by_Method")


# ---- V3: efficiency, four lines ----
fig3, ax = plt.subplots(figsize=(10.5, 5.5))

ax.plot(case_idx, merged["StageEfficiency_pct"], 'o-', color=CLR_DPT,
        label='MATLAB stage eff. (no gate drive)')
ax.plot(case_idx, merged["Efficiency_Stage_pct"], 's--', color=CLR_PSIM,
        markerfacecolor='white', markeredgewidth=1.8,
        label='PSIM stage eff. (no gate drive)')

if has_psim_sys_eff:
    ax.plot(case_idx, cmp["SystemEfficiency_pct"], '^-', color=CLR_TRANS,
            label='MATLAB system eff. (incl. gate drive)')
    ax.plot(case_idx, cmp["Efficiency_System_Optional_pct"], 'd--',
            color=CLR_MILL, markerfacecolor='white', markeredgewidth=1.8,
            label='PSIM system eff. (stage + gate drive)')

ax.set_xticks(case_idx)
ax.set_xticklabels(case_labels, rotation=35, ha='right')
ax.set_xlabel('Operating point  (V$_{in}$ / I$_o$)')
ax.set_ylabel('Efficiency  $\\eta$  (%)')

subtitle = (f'stage |diff| = {np.nanmean(np.abs(cmp["Eff_stage_diff_pp"])):.2f} pp')
if has_psim_sys_eff:
    subtitle += f',  system |diff| = {np.nanmean(np.abs(cmp["Eff_sys_diff_pp"])):.2f} pp'
ax.set_title(f'ZETA efficiency — like-for-like comparison\n({subtitle})')
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc='lower left', frameon=True, fontsize=9)
ax.set_xlim(0.5, 9.5)
save_fig(fig3, "Fig_V3_ZETA_Efficiency_Comparison")


# ---- V4: prediction vs measurement, 1:1 with ±25 % band ----
fig4, ax = plt.subplots(figsize=(7.0, 6.0))
ax.scatter(cmp["P_sw_zeta_W"], cmp["Q1_sw_dpt_W"],
           s=90, color=CLR_DPT, edgecolor='black', linewidth=0.6,
           label=f'DPT (mean {mean_dpt:.1f}%)', zorder=5)
ax.scatter(cmp["P_sw_zeta_W"], cmp["Q1_sw_transition_W"],
           s=70, color=CLR_TRANS, edgecolor='black', linewidth=0.6,
           marker='s', label=f'Transition (mean {mean_trans:.0f}%)', zorder=4)
ax.scatter(cmp["P_sw_zeta_W"], cmp["Q1_sw_miller_W"],
           s=70, color=CLR_MILL, edgecolor='black', linewidth=0.6,
           marker='^', label=f'Miller (mean {mean_mill:.0f}%)', zorder=4)

lo = 0.0
hi = max(cmp["Q1_sw_miller_W"].max(), cmp["P_sw_zeta_W"].max()) * 1.1
ax.plot([lo, hi], [lo, hi], '--', color='gray', linewidth=1.3,
        label='1:1 reference')
ax.fill_between([lo, hi], [lo*0.75, hi*0.75], [lo*1.25, hi*1.25],
                color='gray', alpha=0.10, label='±25% band')
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_xlabel('PSIM measured  $P_{sw}$  (W)')
ax.set_ylabel('MATLAB predicted  $P_{sw}$  (W)')
ax.set_title('ZETA prediction vs measurement — all 9 operating points')
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc='upper left', frameon=True, fontsize=9)
save_fig(fig4, "Fig_V4_ZETA_Prediction_vs_Measurement")


# =============================================================================
# Ranking check
# =============================================================================
# Reports actual numbers only. No "success" / "fail" verdict -- the report
# interprets these in context. The ranking flag just checks whether DPT is
# still the lowest-error method, as claimed in the review.

print("\n" + "=" * 95)
print("RANKING CHECK")
print("=" * 95)

print(f"  DPT        : mean = {mean_dpt:5.1f}%,  max = {max_dpt:5.1f}%")
print(f"  Transition : mean = {mean_trans:5.1f}%,  max = {max_trans:5.1f}%")
print(f"  Miller     : mean = {mean_mill:5.1f}%,  max = {max_mill:5.1f}%")

ranking_ok = (mean_dpt < mean_trans) and (mean_dpt < mean_mill)
print(f"\n  Ranking matches review Table 2 (DPT lowest error): {ranking_ok}")
print("\n  Reference ranges (for the report's interpretation, not enforced):")
print("    DPT mean error < 20%       -> good")
print("    20% to 40%                 -> warrants investigation")
print("    > 40%                      -> lookup or current scaling is wrong")

print(f"\n  Efficiency  |diff| (stage-vs-stage)  : "
      f"{np.nanmean(np.abs(cmp['Eff_stage_diff_pp'])):.2f} pp")
if has_psim_sys_eff:
    print(f"  Efficiency  |diff| (system-vs-system): "
          f"{np.nanmean(np.abs(cmp['Eff_sys_diff_pp'])):.2f} pp")

print("\n" + "=" * 95)
print("Done")
print("=" * 95)
print(f"\nAll outputs in: {OUT_DIR}")