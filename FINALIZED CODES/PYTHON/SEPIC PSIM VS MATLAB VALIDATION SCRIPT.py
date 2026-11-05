"""
Compares the SEPIC MATLAB model against PSIM measurements over the nine
operating points. Reads the MATLAB summary and component-loss CSVs plus
the PSIM summary / loss / validation CSVs, merges them into a single
case table, writes the master and switching-loss validation CSVs, and
emits four figures.

Efficiency is reported two ways: stage-vs-stage and system-vs-system. The
system comparison is skipped automatically if the PSIM summary doesn't
carry Efficiency_System_Optional_pct.
"""

from pathlib import Path
import numpy as np
import pandas as pd


# =============================================================================
# Paths
# =============================================================================

MATLAB_DIR      = Path(r"C:\Users\User\Desktop\MATLAB SEPIC")
PSIM_OUTPUT_DIR = Path(r"C:\Users\User\Desktop\SEPIC NEW\SEPIC_analysis_output")

matlab_summary_path = MATLAB_DIR / "SEPIC_v9_12_summary.csv"
matlab_losses_path  = MATLAB_DIR / "SEPIC_v9_12_component_losses.csv"
psim_summary_path   = PSIM_OUTPUT_DIR / "SEPIC_summary.csv"
psim_losses_path    = PSIM_OUTPUT_DIR / "SEPIC_loss_breakdown.csv"
psim_val_path       = PSIM_OUTPUT_DIR / "SEPIC_validation.csv"

OUT_DIR = MATLAB_DIR / "validation_output"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Load
# =============================================================================

def load_csv(path):
    if not path.exists():
        raise FileNotFoundError(f"Required file not found:\n  {path}")
    return pd.read_csv(path)

mat_summary = load_csv(matlab_summary_path)
mat_losses  = load_csv(matlab_losses_path)
psim_summ   = load_csv(psim_summary_path)
psim_losses = load_csv(psim_losses_path)
psim_val    = load_csv(psim_val_path)

print("=" * 95)
print("SEPIC MATLAB vs PSIM validation")
print("=" * 95)


# =============================================================================
# Case numbers
# =============================================================================

def assign_case(df, vin_col, io_col):
    """Assign case numbers by (Vin, Io) sort if the CSV doesn't carry them."""
    if "Case" in df.columns:
        return df.copy()
    df = df.copy().sort_values([vin_col, io_col]).reset_index(drop=True)
    df.insert(0, "Case", range(1, len(df) + 1))
    return df

mat_summary = assign_case(mat_summary, "Vin_V", "Io_A")
mat_losses  = assign_case(mat_losses,  "Vin",   "Io")


# =============================================================================
# Build the master comparison table
# =============================================================================

mat_psw = mat_losses[[
    "Case", "Q1_cond_W",
    "Q1_sw_dpt_W", "Q1_sw_transition_W", "Q1_sw_miller_W"
]].copy()

mat_summ_slim = mat_summary[[
    "Case", "Vin_V", "Io_A", "Duty",
    "StageEfficiency_pct", "SystemEfficiency_pct"
]].copy()

mat_side = pd.merge(mat_summ_slim, mat_psw, on="Case", validate="one_to_one")

psim_psw = psim_losses[[
    "Case", "MOSFET_Loss_W",
    "P_sw_sepic_W", "E_on_sepic_uJ", "E_off_sepic_uJ",
    "I_sw_turnoff_avg_A", "I_sw_turnon_avg_A", "t_sw_effective_ns"
]].copy()

# System efficiency is optional; PSIM only exports it when the probe
# exists in the schematic.
psim_summ_cols = [
    "Case", "Duty_Measured", "Efficiency_Stage_pct",
    "Extracted_Stage_Loss_Sum_W"
]
if "Efficiency_System_Optional_pct" in psim_summ.columns:
    psim_summ_cols.append("Efficiency_System_Optional_pct")
    has_psim_sys_eff = True
else:
    has_psim_sys_eff = False
    print("WARNING: PSIM summary lacks Efficiency_System_Optional_pct.")
    print("         System-vs-system efficiency comparison will be skipped.")

psim_summ_slim = psim_summ[psim_summ_cols].copy()
psim_side      = pd.merge(psim_summ_slim, psim_psw, on="Case", validate="one_to_one")

cmp = pd.merge(mat_side, psim_side, on="Case", validate="one_to_one")

if len(cmp) != 9:
    raise ValueError(f"Expected 9 cases, got {len(cmp)}")

print(f"\nMerged: {len(cmp)} cases")


# =============================================================================
# Error metrics
# =============================================================================

def pct_error(est, ref):
    est = np.asarray(est, dtype=float)
    ref = np.asarray(ref, dtype=float)
    return np.where(np.abs(ref) > 1e-12, 100.0 * (est - ref) / ref, np.nan)

cmp["DPT_err_pct"]   = pct_error(cmp["Q1_sw_dpt_W"],        cmp["P_sw_sepic_W"])
cmp["Trans_err_pct"] = pct_error(cmp["Q1_sw_transition_W"], cmp["P_sw_sepic_W"])
cmp["Mill_err_pct"]  = pct_error(cmp["Q1_sw_miller_W"],     cmp["P_sw_sepic_W"])

cmp["Duty_diff"] = np.abs(cmp["Duty"] - cmp["Duty_Measured"])

cmp["Eff_stage_diff_pp"] = cmp["StageEfficiency_pct"] - cmp["Efficiency_Stage_pct"]
if has_psim_sys_eff:
    cmp["Eff_sys_diff_pp"] = cmp["SystemEfficiency_pct"] - cmp["Efficiency_System_Optional_pct"]
else:
    cmp["Eff_sys_diff_pp"] = np.nan


# =============================================================================
# Printed validation tables
# =============================================================================

print("\n" + "=" * 95)
print("TABLE 1 — SWITCHING LOSS: MATLAB PREDICTIONS vs PSIM MEASURED")
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
          f"{r['P_sw_sepic_W']:>8.4f} | "
          f"{r['DPT_err_pct']:>+8.1f} {r['Trans_err_pct']:>+9.1f} "
          f"{r['Mill_err_pct']:>+8.1f}")

print("-" * len(hdr))

mean_dpt   = np.nanmean(np.abs(cmp["DPT_err_pct"]))
max_dpt    = np.nanmax(np.abs(cmp["DPT_err_pct"]))
mean_trans = np.nanmean(np.abs(cmp["Trans_err_pct"]))
max_trans  = np.nanmax(np.abs(cmp["Trans_err_pct"]))
mean_mill  = np.nanmean(np.abs(cmp["Mill_err_pct"]))
max_mill   = np.nanmax(np.abs(cmp["Mill_err_pct"]))

print(f"\nSUMMARY (absolute error vs PSIM measured):")
print(f"  DPT          : mean = {mean_dpt:5.1f}%,  max = {max_dpt:5.1f}%")
print(f"  Transition   : mean = {mean_trans:5.1f}%,  max = {max_trans:5.1f}%")
print(f"  Miller       : mean = {mean_mill:5.1f}%,  max = {max_mill:5.1f}%")

print("\n" + "=" * 95)
print("TABLE 2 — EFFICIENCY: STAGE vs STAGE, SYSTEM vs SYSTEM")
print("=" * 95)
for _, r in cmp.iterrows():
    print(f"  Case {int(r['Case'])}: "
          f"stage   MATLAB {r['StageEfficiency_pct']:>6.2f}%  |  "
          f"PSIM {r['Efficiency_Stage_pct']:>6.2f}%  |  "
          f"diff {r['Eff_stage_diff_pp']:>+6.2f} pp")
    if has_psim_sys_eff:
        print(f"           system  MATLAB {r['SystemEfficiency_pct']:>6.2f}%  |  "
              f"PSIM {r['Efficiency_System_Optional_pct']:>6.2f}%  |  "
              f"diff {r['Eff_sys_diff_pp']:>+6.2f} pp")
print(f"\n  mean |diff| (stage)  = {np.nanmean(np.abs(cmp['Eff_stage_diff_pp'])):.2f} pp")
if has_psim_sys_eff:
    print(f"  mean |diff| (system) = {np.nanmean(np.abs(cmp['Eff_sys_diff_pp'])):.2f} pp")


# =============================================================================
# CSV exports
# =============================================================================

master_out = OUT_DIR / "validation_master.csv"
cmp.to_csv(master_out, index=False)

table4 = cmp[[
    "Case", "Vin_V", "Io_A",
    "Q1_sw_dpt_W", "Q1_sw_transition_W", "Q1_sw_miller_W",
    "P_sw_sepic_W",
    "DPT_err_pct", "Trans_err_pct", "Mill_err_pct",
]].copy()
table4.columns = [
    "Case", "Vin_V", "Io_A",
    "MATLAB_DPT_W", "MATLAB_Transition_W", "MATLAB_Miller_W",
    "PSIM_Measured_W",
    "DPT_error_pct", "Transition_error_pct", "Miller_error_pct",
]
table4.to_csv(OUT_DIR / "Table_switching_loss_validation.csv", index=False)

summary_stats = pd.DataFrame({
    "Method":             ["DPT", "Transition", "Miller"],
    "Mean_abs_error_pct": [mean_dpt,   mean_trans,   mean_mill],
    "Max_abs_error_pct":  [max_dpt,    max_trans,    max_mill],
})
summary_stats.to_csv(OUT_DIR / "Table_error_summary.csv", index=False)

print(f"\nCSVs saved to: {OUT_DIR}")


# =============================================================================
# Figures
# =============================================================================

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

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

def save(fig, name):
    fig.savefig(OUT_DIR / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name} (.pdf, .png)")

print("\nGenerating figures...")


# ---- V1: switching loss vs case ----
fig1, ax = plt.subplots(figsize=(9.5, 5.5))
ax.plot(case_idx, cmp["Q1_sw_dpt_W"],        'o-', color=CLR_DPT,
        label='MATLAB DPT (§3.4)',   zorder=4)
ax.plot(case_idx, cmp["Q1_sw_transition_W"], 's-', color=CLR_TRANS,
        label='MATLAB transition (§3.2 gen-1)', zorder=3)
ax.plot(case_idx, cmp["Q1_sw_miller_W"],     '^-', color=CLR_MILL,
        label='MATLAB miller (§3.2 gen-2 partial)', zorder=3)
ax.plot(case_idx, cmp["P_sw_sepic_W"],       'D-', color=CLR_PSIM,
        markerfacecolor='white', markeredgewidth=1.8,
        label='PSIM measured', zorder=5)
ax.set_xticks(case_idx)
ax.set_xticklabels(case_labels, rotation=35, ha='right')
ax.set_xlabel('Operating point  (V$_{in}$ / I$_o$)')
ax.set_ylabel('Switching loss  $P_{sw}$  (W)')
ax.set_title('SEPIC switching-loss prediction vs PSIM measurement\n'
             f'(DPT mean error = {mean_dpt:.1f}%,  max = {max_dpt:.1f}%)')
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc='upper left', frameon=True)
ax.set_xlim(0.5, 9.5)
save(fig1, "Fig_V1_SwitchingLoss_vs_Case")


# ---- V2: absolute error by method ----
fig2, ax = plt.subplots(figsize=(10, 5.5))
x       = np.arange(len(cmp))
width   = 0.27

ax.bar(x - width, np.abs(cmp["DPT_err_pct"]),
       width, color=CLR_DPT,   edgecolor='black', linewidth=0.5,
       label=f'DPT (mean {mean_dpt:.1f}%)')
ax.bar(x,         np.abs(cmp["Trans_err_pct"]),
       width, color=CLR_TRANS, edgecolor='black', linewidth=0.5,
       label=f'Transition (mean {mean_trans:.0f}%)')
ax.bar(x + width, np.abs(cmp["Mill_err_pct"]),
       width, color=CLR_MILL,  edgecolor='black', linewidth=0.5,
       label=f'Miller (mean {mean_mill:.0f}%)')

for xi, ei in zip(x - width, np.abs(cmp["DPT_err_pct"])):
    ax.text(xi, ei + 5, f"{ei:.0f}", ha='center', va='bottom', fontsize=8,
            color=CLR_DPT, fontweight='bold')

ax.set_xticks(x)
ax.set_xticklabels(case_labels, rotation=35, ha='right')
ax.set_xlabel('Operating point  (V$_{in}$ / I$_o$)')
ax.set_ylabel('Absolute error vs PSIM  (%)')
ax.set_title('SEPIC switching-loss error by method')
ax.set_yscale('log')
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.set_ylim(1, max(np.nanmax(np.abs(cmp["Trans_err_pct"])),
                   np.nanmax(np.abs(cmp["Mill_err_pct"]))) * 1.6)
ax.axhline(25, color='gray', linestyle='--', linewidth=1.0, alpha=0.7)
ax.text(9.4, 25, '25% tolerance', ha='right', va='bottom',
        fontsize=9, color='gray', style='italic')
ax.grid(True, axis='y', linestyle=':', alpha=0.5)
ax.legend(loc='upper left', frameon=True)
save(fig2, "Fig_V2_Error_by_Method")


# ---- V3: efficiency, four lines ----
fig3, ax = plt.subplots(figsize=(10.5, 5.5))
ax.plot(case_idx, cmp["StageEfficiency_pct"], 'o-', color=CLR_DPT,
        label='MATLAB stage eff. (no gate drive)')
ax.plot(case_idx, cmp["Efficiency_Stage_pct"], 's--', color=CLR_PSIM,
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
ax.set_title(f'SEPIC efficiency — like-for-like comparison\n({subtitle})')
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc='lower left', frameon=True, fontsize=9)
ax.set_xlim(0.5, 9.5)
save(fig3, "Fig_V3_Efficiency_Comparison")


# ---- V4: prediction vs measurement, 1:1 with ±25 % band ----
fig4, ax = plt.subplots(figsize=(7.0, 6.0))
ax.scatter(cmp["P_sw_sepic_W"], cmp["Q1_sw_dpt_W"],
           s=90, color=CLR_DPT, edgecolor='black', linewidth=0.6,
           label=f'DPT (mean {mean_dpt:.1f}%)',   zorder=5)
ax.scatter(cmp["P_sw_sepic_W"], cmp["Q1_sw_transition_W"],
           s=70, color=CLR_TRANS, edgecolor='black', linewidth=0.6,
           marker='s', label=f'Transition (mean {mean_trans:.0f}%)', zorder=4)
ax.scatter(cmp["P_sw_sepic_W"], cmp["Q1_sw_miller_W"],
           s=70, color=CLR_MILL, edgecolor='black', linewidth=0.6,
           marker='^', label=f'Miller (mean {mean_mill:.0f}%)',   zorder=4)

lo = 0.0
hi = max(cmp["Q1_sw_miller_W"].max(), cmp["P_sw_sepic_W"].max()) * 1.1
ax.plot([lo, hi], [lo, hi], '--', color='gray', linewidth=1.3,
        label='1:1 reference')
ax.fill_between([lo, hi], [lo*0.75, hi*0.75], [lo*1.25, hi*1.25],
                color='gray', alpha=0.10, label='±25% band')

ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_xlabel('PSIM measured  $P_{sw}$  (W)')
ax.set_ylabel('MATLAB predicted  $P_{sw}$  (W)')
ax.set_title('SEPIC prediction vs measurement — all 9 operating points')
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc='upper left', frameon=True, fontsize=9)
save(fig4, "Fig_V4_Prediction_vs_Measurement")


# =============================================================================
# Ranking check
# =============================================================================
# Reports the numbers only. Doesn't enforce "success"/"fail" thresholds;
# the report does the interpretation.

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