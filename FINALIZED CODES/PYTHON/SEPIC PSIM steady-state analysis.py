"""
SEPIC PSIM steady-state analysis.

Reads the nine PSIM TXT exports, computes steady-state quantities over the
last NUM_STEADY_STATE_CYCLES cycles, and writes four CSVs summarising
regulation, loss, and closure-error checks.

Switching energies are extracted by directly integrating v_DS * i_DS over
each edge. The integration window runs from 200 ns before the edge to
500 ns after it. The original narrower window (50 ns pre, 150 ns post)
missed the V_DS/IDS crossover in some exports and produced P_sw values
roughly five times too low. The wider window matches what ZETA v2.5
needs, so both topologies use the same extractor.

VERBOSE=True prints the window sample count, peak / mean instantaneous
power, and integrated energy per edge, so you can confirm the window
actually contains the crossover.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import re
import sys
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# User settings
# =============================================================================

BASE_PATH  = Path(r"C:\Users\User\Desktop\SEPIC NEW")
OUTPUT_DIR = BASE_PATH / "SEPIC_analysis_output"

NUM_STEADY_STATE_CYCLES = 10

VO_TARGET   = 12.0
FSW_TARGET  = 340e3

R_DCR_L1  = 50e-3
R_DCR_L2  = 50e-3
R_ESR_CIN = 10e-3
R_ESR_C1  = 10e-3
R_ESR_C2  = 10e-3

VGATE = 10.0
QG    = 48e-9

VIN_TOLERANCE_V               = 0.10
MAX_ALLOWED_VOUT_ERROR_PCT    = 5.0
MAX_ALLOWED_CLOSURE_ERROR_PCT = 3.0
MAX_ALLOWED_FSW_ERROR_PCT     = 2.0
MAX_ALLOWED_DUTY_ERROR_ABS    = 0.01

CCM_EPS_A = 1e-6

# Diagnostics
VERBOSE                 = True     # print per-case [diag] output
MIN_SWITCHING_W         = 0.10     # below this, warn the window is wrong
SWITCHING_PRE_NS        = 200.0    # pre-edge window (ns)
SWITCHING_POST_NS       = 500.0    # post-edge window (ns)
TURNON_SAMPLE_OFFSET_NS = 100.0    # turn-on current sampling offset (ns)


# =============================================================================
# Operating points
# =============================================================================

CASES_SETUP = [
    {"Case": 1, "Vin_Commanded_V":  9.0, "Rload_Ohm": 24.0, "Duty_Commanded": 0.583080},
    {"Case": 2, "Vin_Commanded_V":  9.0, "Rload_Ohm": 12.0, "Duty_Commanded": 0.585620},
    {"Case": 3, "Vin_Commanded_V":  9.0, "Rload_Ohm":  8.0, "Duty_Commanded": 0.588490},
    {"Case": 4, "Vin_Commanded_V": 12.0, "Rload_Ohm": 24.0, "Duty_Commanded": 0.511020},
    {"Case": 5, "Vin_Commanded_V": 12.0, "Rload_Ohm": 12.0, "Duty_Commanded": 0.512600},
    {"Case": 6, "Vin_Commanded_V": 12.0, "Rload_Ohm":  8.0, "Duty_Commanded": 0.514360},
    {"Case": 7, "Vin_Commanded_V": 15.0, "Rload_Ohm": 24.0, "Duty_Commanded": 0.454900},
    {"Case": 8, "Vin_Commanded_V": 15.0, "Rload_Ohm": 12.0, "Duty_Commanded": 0.456010},
    {"Case": 9, "Vin_Commanded_V": 15.0, "Rload_Ohm":  8.0, "Duty_Commanded": 0.457230},
]


# =============================================================================
# PSIM column mapping (must match the probe names in SEPIC.psimsch)
# =============================================================================

COL = {
    "time":           "Time",
    "vin":            "VIn1",
    "vout":           "Vout1",
    "iout":           "IR1",
    "pwm":            "PWMLogic",
    "iin":            "IL3",
    "il1":            "IL3",
    "il2":            "IL4",
    "icin":           "ICin1",
    "ic1":            "IC3",
    "ic2":            "IC4",
    "mosfet_current": "I(IRLZ44NPbF.IDS)",
    "mosfet_vds":     "V(IRLZ44NPbF.VDS)",
    "diode_current":  "I(SS34)",
    "diode_voltage":  "V(SS34)",
}


# =============================================================================
# Data structures
# =============================================================================

@dataclass
class SelectedFile:
    Case:                 int
    Selected_File:        str
    Samples_Original:     int
    Samples_Analyzed:     int
    Complete_Cycles_Used: int


# =============================================================================
# General helpers
# =============================================================================

def pass_fail(condition: bool) -> str:
    return "PASS" if condition else "FAIL"


def manual_trapezoidal_integral(time: np.ndarray, waveform: np.ndarray) -> float:
    time     = np.asarray(time,     dtype=float)
    waveform = np.asarray(waveform, dtype=float)
    if len(time) != len(waveform):
        raise ValueError("Time and waveform vectors have different lengths.")
    if len(time) < 2:
        raise ValueError("At least two samples are required for integration.")
    delta_time = np.diff(time)
    if np.any(delta_time <= 0):
        raise ValueError("Time samples must increase monotonically.")
    return float(np.sum(0.5 * (waveform[:-1] + waveform[1:]) * delta_time))


def time_average(time: np.ndarray, waveform: np.ndarray) -> float:
    duration = float(time[-1] - time[0])
    if duration <= 0:
        raise ValueError("Invalid time interval.")
    return manual_trapezoidal_integral(time, waveform) / duration


def rms(time: np.ndarray, waveform: np.ndarray) -> float:
    return float(np.sqrt(time_average(time, np.asarray(waveform, dtype=float) ** 2)))


def peak_to_peak(waveform: np.ndarray) -> float:
    waveform = np.asarray(waveform, dtype=float)
    return float(np.max(waveform) - np.min(waveform))


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = (
        df.columns.astype(str).str.strip().str.replace('"', "", regex=False)
    )
    return df


def orient_current_positive(waveform: np.ndarray) -> np.ndarray:
    """For CCM detection. PSIM sometimes exports inductor currents with the
    opposite sign convention; flip so min/max and CCM checks read sensibly."""
    waveform = np.asarray(waveform, dtype=float)
    if np.mean(waveform) < 0:
        return -waveform
    return waveform


# =============================================================================
# File detection and reading
# =============================================================================

def candidate_files_for_case(case_number: int) -> List[Path]:
    pattern = re.compile(
        rf"^SEPIC_case_{case_number}_steady_state(?:\s*\(\d+\))?\.txt$",
        flags=re.IGNORECASE,
    )
    matches = [p for p in BASE_PATH.glob("*.txt") if pattern.match(p.name)]
    exact = BASE_PATH / f"SEPIC_case_{case_number}_steady_state.txt"
    if exact.exists() and exact not in matches:
        matches.append(exact)
    return sorted(matches, key=lambda p: p.stat().st_mtime, reverse=True)


def select_file_for_case(case_number: int) -> Path:
    """Most recently modified matching file wins -- covers the case where
    PSIM writes "SEPIC_case_3_steady_state (2).txt" after a re-run."""
    matches = candidate_files_for_case(case_number)
    if not matches:
        raise FileNotFoundError(
            f"No TXT export found for Case {case_number}. "
            f"Expected: SEPIC_case_{case_number}_steady_state.txt"
        )
    return matches[0]


def read_psim_txt(file_path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(file_path, sep=r"\s+", engine="python")
    except Exception as exc:
        raise ValueError(f"Unable to read {file_path.name}: {exc}") from exc

    df = clean_columns(df)

    required_columns = sorted(set(COL.values()))
    missing = [c for c in required_columns if c not in df.columns]
    if missing:
        raise KeyError(
            f"\nMissing required columns in {file_path.name}:\n"
            f"  {', '.join(missing)}\n\n"
            f"Available columns:\n  {', '.join(df.columns)}\n"
        )

    for column in required_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df = df.dropna(subset=required_columns).reset_index(drop=True)
    if len(df) < 20:
        raise ValueError(f"{file_path.name} contains too few valid samples.")
    return df


# =============================================================================
# PWM cycle detection
# =============================================================================

def pwm_binary(df: pd.DataFrame) -> np.ndarray:
    return (df[COL["pwm"]].to_numpy(dtype=float) > 0.5).astype(int)


def rising_edge_indices(df: pd.DataFrame) -> np.ndarray:
    pwm = pwm_binary(df)
    return np.where((pwm[:-1] == 0) & (pwm[1:] == 1))[0] + 1


def falling_edge_indices(df: pd.DataFrame) -> np.ndarray:
    pwm = pwm_binary(df)
    return np.where((pwm[:-1] == 1) & (pwm[1:] == 0))[0]


def select_last_complete_cycles(
    df: pd.DataFrame, requested_cycles: int,
) -> Tuple[pd.DataFrame, int]:
    """Keep only the tail of the file that contains whole PWM cycles, so
    every integral is over a fixed duration and there's no partial-cycle
    bias in the averages."""
    edges = rising_edge_indices(df)
    if len(edges) < 2:
        print("  Warning: PWM rising edges not detected. Using full dataset.")
        return df.copy().reset_index(drop=True), 0

    available_cycles = len(edges) - 1
    cycles_to_use    = min(requested_cycles, available_cycles)
    start_edge_pos   = len(edges) - 1 - cycles_to_use
    start_index      = int(edges[start_edge_pos])
    end_index        = int(edges[-1])

    selected = df.iloc[start_index:end_index + 1].copy()
    selected.reset_index(drop=True, inplace=True)
    return selected, cycles_to_use


def measure_pwm_frequency_and_duty(df: pd.DataFrame) -> Tuple[float, float]:
    time  = df[COL["time"]].to_numpy(dtype=float)
    pwm   = pwm_binary(df).astype(float)
    edges = rising_edge_indices(df)

    fsw_measured = (
        float(1.0 / np.mean(np.diff(time[edges])))
        if len(edges) >= 2 else float("nan")
    )
    duty_measured = float(time_average(time, pwm))
    return fsw_measured, duty_measured


# =============================================================================
# Switching-instant current extraction
# =============================================================================

def extract_switching_instant_currents(
    df: pd.DataFrame,
    turnon_offset_s: float = TURNON_SAMPLE_OFFSET_NS * 1e-9,
) -> Tuple[float, float]:
    """
    Turn-off: last high sample before PWM goes low.
    Turn-on : sample 100 ns AFTER the rising edge (past commutation).
    """
    time_arr = df[COL["time"]].to_numpy(dtype=float)
    mosfet_i = df[COL["mosfet_current"]].to_numpy(dtype=float)

    fall_idx = falling_edge_indices(df)
    if len(fall_idx) > 0:
        i_turnoff_avg = float(np.mean(np.abs(mosfet_i[fall_idx])))
    else:
        print("  Warning: no PWM falling edges; turn-off current = NaN.")
        i_turnoff_avg = float("nan")

    rise_idx = rising_edge_indices(df)
    if len(rise_idx) > 0:
        turnon_sample_indices = []
        for idx in rise_idx:
            t_target   = time_arr[idx] + turnon_offset_s
            candidates = np.where(time_arr >= t_target)[0]
            if len(candidates) > 0:
                turnon_sample_indices.append(candidates[0])
        if turnon_sample_indices:
            i_turnon_avg = float(
                np.mean(np.abs(mosfet_i[np.array(turnon_sample_indices)]))
            )
        else:
            print("  Warning: no turn-on sample after offset; using NaN.")
            i_turnon_avg = float("nan")
    else:
        print("  Warning: no PWM rising edges; turn-on current = NaN.")
        i_turnon_avg = float("nan")

    return i_turnoff_avg, i_turnon_avg


# =============================================================================
# Effective switching transition time
# =============================================================================

def measure_effective_transition_time(df: pd.DataFrame) -> float:
    """V_DS fall time at turn-on (PWM rising edge), averaged over cycles."""
    time_arr = df[COL["time"]].to_numpy(dtype=float)
    vds      = df[COL["mosfet_vds"]].to_numpy(dtype=float)
    rise_idx = rising_edge_indices(df)
    t_sw_list = []

    for idx in rise_idx:
        pre_start = max(0, idx - 5)
        vds_high  = float(np.mean(vds[pre_start:idx]))
        if vds_high < 1.0:
            continue
        vds_threshold = 0.10 * vds_high
        search_end    = min(idx + 500, len(vds) - 1)
        cross_idx     = np.where(vds[idx:search_end] < vds_threshold)[0]
        if len(cross_idx) > 0:
            t_fall = float(time_arr[idx + cross_idx[0]] - time_arr[idx])
            if t_fall > 0:
                t_sw_list.append(t_fall)

    if not t_sw_list:
        print("  Warning: V_DS fall time could not be measured; returning nan.")
        return float("nan")
    return float(np.mean(t_sw_list))


# =============================================================================
# Switching energy extraction
# =============================================================================

def extract_switching_energies(
    df: pd.DataFrame,
    pre_window_s:  float = SWITCHING_PRE_NS  * 1e-9,
    post_window_s: float = SWITCHING_POST_NS * 1e-9,
    verbose: bool = False,
) -> Tuple[float, float]:
    """
    E_on and E_off from v_DS * i_DS integrated over each switching window.
    Direct equivalent of review Eq. (3).

    Window: [t_edge - pre_window_s, t_edge + post_window_s].
    Defaults are 200 ns pre, 500 ns post. The old (50, 150) window was too
    narrow -- it often missed the V_DS/IDS crossover entirely and gave
    switching loss far below the expected 0.2-0.5 W range.
    """
    time_arr = df[COL["time"]].to_numpy(dtype=float)
    vds      = df[COL["mosfet_vds"]].to_numpy(dtype=float)
    ids      = df[COL["mosfet_current"]].to_numpy(dtype=float)
    p_sw     = vds * ids

    def integrate_window(edge_idx: int) -> float:
        t_edge = time_arr[edge_idx]
        t_lo   = t_edge - pre_window_s
        t_hi   = t_edge + post_window_s
        mask   = (time_arr >= t_lo) & (time_arr <= t_hi)
        if mask.sum() < 10:
            return float("nan")
        return manual_trapezoidal_integral(time_arr[mask], p_sw[mask])

    rise_idx = rising_edge_indices(df)
    fall_idx = falling_edge_indices(df)

    E_on_list  = [integrate_window(i) for i in rise_idx]
    E_off_list = [integrate_window(i) for i in fall_idx]

    # Diagnostic print: peak power should be in the hundreds of watts
    # during a real crossover. A mean below ~0.05 W means the window only
    # captures conduction.
    if verbose and len(rise_idx) > 0:
        _e  = rise_idx[0]
        _te = time_arr[_e]
        _m  = (time_arr >= _te - pre_window_s) & (time_arr <= _te + post_window_s)
        if _m.sum() > 0:
            print(f"    [diag] turn-on  window: N={_m.sum():3d}, "
                  f"p_max={p_sw[_m].max():.2f} W, "
                  f"p_mean={p_sw[_m].mean():.4f} W, "
                  f"E_on={integrate_window(_e)*1e6:.4f} uJ")
    if verbose and len(fall_idx) > 0:
        _e  = fall_idx[0]
        _te = time_arr[_e]
        _m  = (time_arr >= _te - pre_window_s) & (time_arr <= _te + post_window_s)
        if _m.sum() > 0:
            print(f"    [diag] turn-off window: N={_m.sum():3d}, "
                  f"p_max={p_sw[_m].max():.2f} W, "
                  f"p_mean={p_sw[_m].mean():.4f} W, "
                  f"E_off={integrate_window(_e)*1e6:.4f} uJ")

    E_on_avg  = float(np.nanmean(E_on_list))  if len(E_on_list)  else float("nan")
    E_off_avg = float(np.nanmean(E_off_list)) if len(E_off_list) else float("nan")
    return E_on_avg, E_off_avg


# =============================================================================
# Case processing
# =============================================================================

def process_case(item: Dict[str, float]) -> Tuple[Dict[str, object], SelectedFile]:
    case_number = int(item["Case"])
    file_path   = select_file_for_case(case_number)
    print(f"\nProcessing Case {case_number}: {file_path.name}")

    df_full = read_psim_txt(file_path)
    df, complete_cycles_used = select_last_complete_cycles(
        df_full, NUM_STEADY_STATE_CYCLES
    )

    time  = df[COL["time"]].to_numpy(dtype=float)
    dt_ns = (time[1] - time[0]) * 1e9 if len(time) > 1 else float("nan")

    vin            = df[COL["vin"]].to_numpy(dtype=float)
    vout           = df[COL["vout"]].to_numpy(dtype=float)
    iout           = df[COL["iout"]].to_numpy(dtype=float)
    il1            = df[COL["il1"]].to_numpy(dtype=float)
    il2            = df[COL["il2"]].to_numpy(dtype=float)
    icin           = df[COL["icin"]].to_numpy(dtype=float)
    ic1            = df[COL["ic1"]].to_numpy(dtype=float)
    ic2            = df[COL["ic2"]].to_numpy(dtype=float)
    mosfet_current = df[COL["mosfet_current"]].to_numpy(dtype=float)
    mosfet_vds     = df[COL["mosfet_vds"]].to_numpy(dtype=float)
    diode_current  = df[COL["diode_current"]].to_numpy(dtype=float)
    diode_voltage  = df[COL["diode_voltage"]].to_numpy(dtype=float)

    fsw_measured, duty_measured = measure_pwm_frequency_and_duty(df)
    fsw_error_pct = (
        100.0 * (fsw_measured - FSW_TARGET) / FSW_TARGET
        if np.isfinite(fsw_measured) else float("nan")
    )
    duty_error_abs = abs(duty_measured - float(item["Duty_Commanded"]))

    i_sw_turnoff_avg, i_sw_turnon_avg = extract_switching_instant_currents(df)
    t_sw_effective = measure_effective_transition_time(df)
    e_on_sepic_j, e_off_sepic_j = extract_switching_energies(df, verbose=VERBOSE)

    e_on_sepic_uJ  = e_on_sepic_j  * 1e6 if np.isfinite(e_on_sepic_j)  else float("nan")
    e_off_sepic_uJ = e_off_sepic_j * 1e6 if np.isfinite(e_off_sepic_j) else float("nan")

    if (np.isfinite(e_on_sepic_j) and np.isfinite(e_off_sepic_j)
            and np.isfinite(fsw_measured)):
        p_sw_sepic_w = (e_on_sepic_j + e_off_sepic_j) * fsw_measured
    else:
        p_sw_sepic_w = float("nan")

    print(f"  I_sw turn-off  = {i_sw_turnoff_avg:.6f} A")
    print(f"  I_sw turn-on   = {i_sw_turnon_avg:.6f} A  "
          f"(sampled {TURNON_SAMPLE_OFFSET_NS:.0f} ns after rising edge)")
    print(f"  t_sw_effective = {t_sw_effective*1e9:.3f} ns")
    print(f"  E_on  (SEPIC)  = {e_on_sepic_uJ:.4f} uJ")
    print(f"  E_off (SEPIC)  = {e_off_sepic_uJ:.4f} uJ")
    print(f"  P_sw  (SEPIC)  = {p_sw_sepic_w:.6f} W  = (E_on + E_off) * fsw")
    print(f"  dt             = {dt_ns:.3f} ns")

    # Sanity check: real switching loss at these operating points is
    # typically 0.2-0.5 W. Anything below 0.1 W almost certainly means the
    # integration window missed the crossover.
    if np.isfinite(p_sw_sepic_w) and p_sw_sepic_w < MIN_SWITCHING_W:
        print(f"  WARNING: P_sw = {p_sw_sepic_w:.4f} W < {MIN_SWITCHING_W} W.")
        print(f"     Likely the switching window missed the V_DS/IDS crossover.")
        print(f"     Try: (a) widen SWITCHING_POST_NS to 800-1000 ns,")
        print(f"          (b) verify probe columns match SEPIC.psimsch,")
        print(f"          (c) inspect [diag] -- p_max should be 100+ W.")

    vin_avg        = time_average(time, vin)
    vout_avg       = time_average(time, vout)
    iout_avg       = time_average(time, iout)
    vout_error_pct = 100.0 * (vout_avg - VO_TARGET) / VO_TARGET
    vout_ripple_pp = peak_to_peak(vout)

    il1_ccm = orient_current_positive(il1)
    il2_ccm = orient_current_positive(il2)

    il1_min       = float(np.min(il1_ccm))
    il1_max       = float(np.max(il1_ccm))
    il1_ripple_pp = peak_to_peak(il1_ccm)
    il1_rms       = rms(time, il1)

    il2_min       = float(np.min(il2_ccm))
    il2_max       = float(np.max(il2_ccm))
    il2_ripple_pp = peak_to_peak(il2_ccm)
    il2_rms       = rms(time, il2)

    ccm_l1      = il1_min > CCM_EPS_A
    ccm_l2      = il2_min > CCM_EPS_A
    ccm_overall = ccm_l1 and ccm_l2

    pin_stage            = time_average(time, vin * il1)
    pout                 = time_average(time, vout * iout)
    loss_pin_minus_pout  = pin_stage - pout
    efficiency_stage_pct = (
        100.0 * pout / pin_stage if pin_stage > 0 else float("nan")
    )

    # MOSFET / diode loss from the instantaneous v*i product. Sign flip
    # because PSIM's passive-sign convention sometimes exports these
    # negative when the device is sourcing power.
    p_mosfet_raw = time_average(time, mosfet_vds * mosfet_current)
    p_diode_raw  = time_average(time, diode_voltage * diode_current)
    p_mosfet     = abs(p_mosfet_raw)
    p_diode      = abs(p_diode_raw)

    p_l1_dcr = (il1_rms ** 2) * R_DCR_L1
    p_l2_dcr = (il2_rms ** 2) * R_DCR_L2

    icin_rms = rms(time, icin)
    ic1_rms  = rms(time, ic1)
    ic2_rms  = rms(time, ic2)

    p_cin_esr       = (icin_rms ** 2) * R_ESR_CIN
    p_c1_esr        = (ic1_rms  ** 2) * R_ESR_C1
    p_c2_esr        = (ic2_rms  ** 2) * R_ESR_C2
    p_cap_esr_total = p_cin_esr + p_c1_esr + p_c2_esr

    # Closure check: the extracted loss sum should match P_in - P_out.
    # A residual much above 3% of P_in means something in the loss model
    # is missing.
    extracted_stage_loss_sum = (
        p_mosfet + p_diode + p_l1_dcr + p_l2_dcr + p_cap_esr_total
    )
    closure_residual_w = pin_stage - pout - extracted_stage_loss_sum
    closure_error_pct  = (
        100.0 * closure_residual_w / pin_stage if pin_stage > 0 else float("nan")
    )

    p_gate_drive_est               = VGATE * QG * FSW_TARGET
    pin_system_optional            = pin_stage + p_gate_drive_est
    efficiency_system_optional_pct = (
        100.0 * pout / pin_system_optional if pin_system_optional > 0 else float("nan")
    )

    vin_status     = abs(vin_avg - float(item["Vin_Commanded_V"])) <= VIN_TOLERANCE_V
    vout_status    = abs(vout_error_pct) <= MAX_ALLOWED_VOUT_ERROR_PCT
    fsw_status     = np.isfinite(fsw_error_pct) and abs(fsw_error_pct) <= MAX_ALLOWED_FSW_ERROR_PCT
    duty_status    = duty_error_abs <= MAX_ALLOWED_DUTY_ERROR_ABS
    closure_status = abs(closure_error_pct) <= MAX_ALLOWED_CLOSURE_ERROR_PCT
    overall_status = all([vin_status, vout_status, fsw_status, duty_status,
                          ccm_overall, closure_status])

    print(f"  Vin_avg        = {vin_avg:.4f} V")
    print(f"  Vout_avg       = {vout_avg:.4f} V ({vout_error_pct:+.2f}% error)")
    print(f"  Vout_ripple_pp = {vout_ripple_pp:.6f} V")
    print(f"  Duty_measured  = {duty_measured:.6f}")
    print(f"  Fsw_measured   = {fsw_measured:.2f} Hz")
    print(f"  Pin_stage      = {pin_stage:.6f} W")
    print(f"  Pout           = {pout:.6f} W")
    print(f"  Efficiency     = {efficiency_stage_pct:.2f}%")
    print(f"  CCM_L1 (min)   = {il1_min:.6f} A  ({pass_fail(ccm_l1)})")
    print(f"  CCM_L2 (min)   = {il2_min:.6f} A  ({pass_fail(ccm_l2)})")
    print(f"  MOSFET raw pwr = {p_mosfet_raw:.6f} W  |  abs = {p_mosfet:.6f} W")
    print(f"  Diode  raw pwr = {p_diode_raw:.6f} W  |  abs = {p_diode:.6f} W")
    print(f"  Closure error  = {closure_error_pct:+.3f}%  ({pass_fail(closure_status)})")
    print(f"  Overall status = {pass_fail(overall_status)}")

    result = {
        "Case":                       case_number,
        "Vin_Commanded_V":            float(item["Vin_Commanded_V"]),
        "Rload_Ohm":                  float(item["Rload_Ohm"]),
        "Duty_Commanded":             float(item["Duty_Commanded"]),
        "Vin_Avg_V":                  vin_avg,
        "Vout_Target_V":              VO_TARGET,
        "Vout_Avg_V":                 vout_avg,
        "Vout_Error_pct":             vout_error_pct,
        "Vout_Ripple_pp_V":           vout_ripple_pp,
        "Iout_Avg_A":                 iout_avg,
        "Duty_Measured":              duty_measured,
        "Duty_Error_Abs":             duty_error_abs,
        "Fsw_Target_Hz":              FSW_TARGET,
        "Fsw_Measured_Hz":            fsw_measured,
        "Fsw_Error_pct":              fsw_error_pct,
        "IL1_Min_A":                  il1_min,
        "IL1_Max_A":                  il1_max,
        "IL1_Ripple_pp_A":            il1_ripple_pp,
        "IL1_RMS_A":                  il1_rms,
        "IL2_Min_A":                  il2_min,
        "IL2_Max_A":                  il2_max,
        "IL2_Ripple_pp_A":            il2_ripple_pp,
        "IL2_RMS_A":                  il2_rms,
        "Pin_Stage_W":                pin_stage,
        "Pout_W":                     pout,
        "Loss_Pin_minus_Pout_W":      loss_pin_minus_pout,
        "Efficiency_Stage_pct":       efficiency_stage_pct,
        "MOSFET_Loss_Raw_W":          p_mosfet_raw,
        "MOSFET_Loss_W":              p_mosfet,
        "Diode_Loss_Raw_W":           p_diode_raw,
        "Diode_Loss_W":               p_diode,
        "I_sw_turnoff_avg_A":         i_sw_turnoff_avg,
        "I_sw_turnon_avg_A":          i_sw_turnon_avg,
        "t_sw_effective_ns":          t_sw_effective * 1e9 if np.isfinite(t_sw_effective) else float("nan"),
        "E_on_sepic_uJ":              e_on_sepic_uJ,
        "E_off_sepic_uJ":             e_off_sepic_uJ,
        "P_sw_sepic_W":               p_sw_sepic_w,
        "L1_DCR_Loss_W":              p_l1_dcr,
        "L2_DCR_Loss_W":              p_l2_dcr,
        "ICin_RMS_A":                 icin_rms,
        "IC1_RMS_A":                  ic1_rms,
        "IC2_RMS_A":                  ic2_rms,
        "Cin_ESR_Loss_W":             p_cin_esr,
        "C1_ESR_Loss_W":              p_c1_esr,
        "C2_ESR_Loss_W":              p_c2_esr,
        "Cap_ESR_Total_W":            p_cap_esr_total,
        "Extracted_Stage_Loss_Sum_W": extracted_stage_loss_sum,
        "Closure_Residual_W":         closure_residual_w,
        "Closure_Error_pct":          closure_error_pct,
        "Gate_Drive_Loss_Est_W":      p_gate_drive_est,
        "Pin_System_Optional_W":      pin_system_optional,
        "Efficiency_System_Optional_pct": efficiency_system_optional_pct,
        "Vin_Status":                 pass_fail(vin_status),
        "Vout_Regulation_Status":     pass_fail(vout_status),
        "Fsw_Status":                 pass_fail(fsw_status),
        "Duty_Status":                pass_fail(duty_status),
        "CCM_L1_Status":              pass_fail(ccm_l1),
        "CCM_L2_Status":              pass_fail(ccm_l2),
        "CCM_Overall_Status":         pass_fail(ccm_overall),
        "Closure_Status":             pass_fail(closure_status),
        "Overall_Status":             pass_fail(overall_status),
    }

    selected_file = SelectedFile(
        Case=case_number,
        Selected_File=file_path.name,
        Samples_Original=len(df_full),
        Samples_Analyzed=len(df),
        Complete_Cycles_Used=complete_cycles_used,
    )
    return result, selected_file


# =============================================================================
# Output tables
# =============================================================================

def round_numeric_columns(df: pd.DataFrame, decimal_places: int = 6) -> pd.DataFrame:
    df = df.copy()
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    df[numeric_cols] = df[numeric_cols].round(decimal_places)
    return df


def build_and_save_tables(
    results: List[Dict[str, object]],
    selected_files: List[SelectedFile],
) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df          = round_numeric_columns(pd.DataFrame(results))
    df_selected = pd.DataFrame([asdict(item) for item in selected_files])

    summary_columns = [
        "Case", "Vin_Commanded_V", "Rload_Ohm", "Duty_Commanded",
        "Duty_Measured", "Fsw_Measured_Hz",
        "Vout_Avg_V", "Vout_Error_pct", "Vout_Ripple_pp_V",
        "Pin_Stage_W", "Pout_W", "Loss_Pin_minus_Pout_W",
        "Efficiency_Stage_pct", "Extracted_Stage_Loss_Sum_W",
        "Closure_Error_pct", "CCM_Overall_Status", "Overall_Status",
    ]
    loss_columns = [
        "Case", "Vin_Commanded_V", "Rload_Ohm",
        "MOSFET_Loss_Raw_W", "MOSFET_Loss_W",
        "Diode_Loss_Raw_W",  "Diode_Loss_W",
        "I_sw_turnoff_avg_A", "I_sw_turnon_avg_A",
        "t_sw_effective_ns",
        "E_on_sepic_uJ", "E_off_sepic_uJ", "P_sw_sepic_W",
        "L1_DCR_Loss_W", "L2_DCR_Loss_W",
        "Cin_ESR_Loss_W", "C1_ESR_Loss_W", "C2_ESR_Loss_W",
        "Cap_ESR_Total_W",
        "Extracted_Stage_Loss_Sum_W", "Loss_Pin_minus_Pout_W",
        "Closure_Residual_W", "Closure_Error_pct",
        "Gate_Drive_Loss_Est_W", "Efficiency_System_Optional_pct",
    ]
    validation_columns = [
        "Case", "Vin_Status", "Vout_Regulation_Status",
        "Fsw_Status", "Duty_Status",
        "CCM_L1_Status", "CCM_L2_Status", "CCM_Overall_Status",
        "Closure_Status", "Overall_Status",
    ]

    summary_path    = OUTPUT_DIR / "SEPIC_summary.csv"
    loss_path       = OUTPUT_DIR / "SEPIC_loss_breakdown.csv"
    validation_path = OUTPUT_DIR / "SEPIC_validation.csv"
    files_path      = OUTPUT_DIR / "SEPIC_selected_files.csv"

    df[summary_columns].to_csv(summary_path,       index=False)
    df[loss_columns].to_csv(loss_path,             index=False)
    df[validation_columns].to_csv(validation_path, index=False)
    df_selected.to_csv(files_path,                 index=False)

    print("\n" + "=" * 95)
    print("SEPIC SUMMARY")
    print("=" * 95)
    print(df[summary_columns].to_string(index=False))

    print("\n" + "=" * 95)
    print("SEPIC LOSS BREAKDOWN")
    print("=" * 95)
    print(df[loss_columns].to_string(index=False))

    print("\n" + "=" * 95)
    print("SEPIC VALIDATION")
    print("=" * 95)
    print(df[validation_columns].to_string(index=False))

    print("\n" + "=" * 95)
    print("CSV FILES CREATED")
    print("=" * 95)
    for path in [summary_path, loss_path, validation_path, files_path]:
        print(path)


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    print("=" * 95)
    print("SEPIC PSIM component-level loss analysis")
    print("=" * 95)
    print(f"Data folder        : {BASE_PATH}")
    print(f"Output folder      : {OUTPUT_DIR}")
    print(f"Steady-state cycles: {NUM_STEADY_STATE_CYCLES}")
    print(f"Verbose diagnostic : {VERBOSE}")
    print(f"Switching window   : pre={SWITCHING_PRE_NS:.0f} ns, "
          f"post={SWITCHING_POST_NS:.0f} ns")

    if not BASE_PATH.exists():
        raise FileNotFoundError(f"\nInput folder not found:\n{BASE_PATH}")

    results:        List[Dict[str, object]] = []
    selected_files: List[SelectedFile]      = []
    failed_cases:   List[int]               = []

    for item in CASES_SETUP:
        try:
            result, selected_file = process_case(item)
            results.append(result)
            selected_files.append(selected_file)
        except Exception as exc:
            case_number = int(item["Case"])
            failed_cases.append(case_number)
            print(f"\n  ERROR in Case {case_number}: {exc}")
            print("  This case was skipped.")

    if not results:
        raise RuntimeError("No valid SEPIC files were processed.")

    build_and_save_tables(results, selected_files)

    if failed_cases:
        print("\nWARNING — the following cases failed:")
        print(", ".join(str(c) for c in failed_cases))
        sys.exit(1)

    print("\n" + "=" * 95)
    print("Done")
    print("=" * 95)


if __name__ == "__main__":
    main()