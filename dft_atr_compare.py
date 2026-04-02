#!/usr/bin/env python3
"""
DFT vs ATR Experimental Spectrum Comparison
============================================
Fits an optimal frequency scaling factor for B3LYP/def2-TZVP calculated
cyclohexane IR spectrum against the ATR experimental measurement.

Uses least-squares optimization on matched peak positions.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.signal import find_peaks, savgol_filter
from scipy.optimize import minimize_scalar

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
FIG_DIR  = os.path.join(BASE_DIR, "figures")

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 11,
    "axes.linewidth": 1.2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
})


# ── Load data ────────────────────────────────────────────────────────────
def load_dft(filename="cyclohexane_IR_def2tzvp.txt"):
    dft = np.loadtxt(os.path.join(DATA_DIR, filename), comments="#")
    return dft[:, 0], dft[:, 1]   # wavenumber, epsilon


def load_atr():
    df = pd.read_csv(os.path.join(DATA_DIR, "ATR_result.CSV"),
                     header=None, names=["wn", "T"])
    wn = df["wn"].values
    y  = df["T"].values.max() - df["T"].values   # invert to absorbance-like
    sort = np.argsort(wn)
    wn = wn[sort]
    y  = savgol_filter(y[sort], window_length=31, polyorder=3)
    return wn, y


# ── Peak detection ───────────────────────────────────────────────────────
def find_dft_peaks(wn, eps, min_eps=0.5):
    """Find significant DFT peaks above min_eps threshold."""
    idx, props = find_peaks(eps, prominence=0.01, distance=3)
    mask = eps[idx] >= min_eps
    return wn[idx[mask]], eps[idx[mask]]


def find_atr_peaks(wn, y):
    """Find real ATR peaks (cyclohexane C,H only region)."""
    # Only consider 600-3100 cm-1 (cyclohexane fundamental range)
    mask = (wn >= 600) & (wn <= 3100)
    xw = wn[mask]
    yw = y[mask]
    idx, props = find_peaks(yw, prominence=0.5, distance=10)
    # Keep only peaks with decent intensity
    strong = yw[idx] > 1.5
    return xw[idx[strong]], yw[idx[strong]]


# ── Peak matching & scaling factor optimization ──────────────────────────
# Known DFT → Experimental peak correspondences for cyclohexane
# Assigned by vibrational mode analysis:
#
#   DFT (unscaled)    ATR (exp)     Mode
#   3108              2929          C-H asym. stretch
#   3052              2856          C-H sym. stretch
#   1512              1450          CH2 scissoring
#   889               879           Ring deformation (chair puckering)
#   791               731           Ring breathing / CH2 rocking

# ── Boat form (cyclohexane_IR_def2tzvp.txt) ──────────────────────────────
PEAK_PAIRS_BOAT = np.array([
    [3108, 2929],   # C-H asymmetric stretch
    [3052, 2856],   # C-H symmetric stretch
    [1512, 1450],   # CH2 scissoring
    [889,  879],    # Ring deformation
    [791,  731],    # Ring breathing
])

# ── Chair form (cyclohexane_IR_def2tzvp_chairform.txt) ───────────────────
# Chair-form DFT only covers C-H stretch region (~2956–3318 cm⁻¹)
# Peaks: 3094 (eq C-H), 3086 (ax C-H), 3036 (sym stretch)
PEAK_PAIRS_CHAIR = np.array([
    [3094, 2929],   # C-H equatorial stretch → ATR asym. C-H
    [3086, 2929],   # C-H axial stretch → ATR asym. C-H (near-degenerate)
    [3036, 2856],   # C-H symmetric stretch → ATR sym. C-H
])

# Active dataset — set by main()
PEAK_PAIRS = PEAK_PAIRS_BOAT


def rmse_at_scale(scale):
    """RMSE between scaled DFT peaks and experimental ATR peaks."""
    scaled = PEAK_PAIRS[:, 0] * scale
    exp    = PEAK_PAIRS[:, 1]
    return np.sqrt(np.mean((scaled - exp) ** 2))


def optimize_scaling():
    """Find optimal single scaling factor by minimizing RMSE."""
    result = minimize_scalar(rmse_at_scale, bounds=(0.90, 1.00), method="bounded")
    return result.x, result.fun


# ── Dual-zone scaling ────────────────────────────────────────────────────
# DFT harmonic frequencies have different anharmonicity in different regions:
#   - High-freq (C-H stretch >2000 cm-1): stronger anharmonic correction needed
#   - Low-freq (skeletal <2000 cm-1): closer to harmonic, milder correction
ZONE_BOUNDARY = 2000  # cm-1

HIGH_PAIRS = PEAK_PAIRS[PEAK_PAIRS[:, 0] >= ZONE_BOUNDARY]   # C-H stretches
LOW_PAIRS  = PEAK_PAIRS[PEAK_PAIRS[:, 0] <  ZONE_BOUNDARY]   # skeletal modes


def optimize_dual_scaling():
    """Separate scaling factors for high-freq and low-freq zones."""
    def rmse_high(s):
        return np.sqrt(np.mean((HIGH_PAIRS[:, 0] * s - HIGH_PAIRS[:, 1]) ** 2))
    def rmse_low(s):
        return np.sqrt(np.mean((LOW_PAIRS[:, 0] * s - LOW_PAIRS[:, 1]) ** 2))

    res_hi = minimize_scalar(rmse_high, bounds=(0.90, 1.00), method="bounded")
    res_lo = minimize_scalar(rmse_low,  bounds=(0.90, 1.05), method="bounded")
    return res_hi.x, res_hi.fun, res_lo.x, res_lo.fun


def apply_dual_scale(wn, s_hi, s_lo):
    """Apply zone-dependent scaling to a wavenumber array."""
    scaled = np.where(wn >= ZONE_BOUNDARY, wn * s_hi, wn * s_lo)
    return scaled


# ── Main ─────────────────────────────────────────────────────────────────
def run_comparison(label, dft_file, peak_pairs, modes, fig_suffix=""):
    """Run full comparison for one conformer."""
    global PEAK_PAIRS, HIGH_PAIRS, LOW_PAIRS
    PEAK_PAIRS = peak_pairs
    HIGH_PAIRS = PEAK_PAIRS[PEAK_PAIRS[:, 0] >= ZONE_BOUNDARY]
    LOW_PAIRS  = PEAK_PAIRS[PEAK_PAIRS[:, 0] <  ZONE_BOUNDARY]

    print(f"\n{'=' * 65}")
    print(f"  {label}")
    print(f"  DFT (B3LYP/def2-TZVP) vs ATR Experimental")
    print(f"{'=' * 65}")

    dft_wn, dft_eps = load_dft(dft_file)
    atr_wn, atr_y   = load_atr()

    # Find peaks
    dft_peaks_wn, dft_peaks_h = find_dft_peaks(dft_wn, dft_eps)
    print(f"\n  DFT peaks (unscaled, eps > 0.5): {len(dft_peaks_wn)}")
    for w, h in zip(dft_peaks_wn, dft_peaks_h):
        print(f"    {w:8.1f} cm-1  (eps = {h:.2f})")

    # ── Scaling ──────────────────────────────────────────────────────────
    has_low = len(LOW_PAIRS) > 0

    scale, rmse = optimize_scaling()
    print(f"\n{'─' * 55}")
    print(f"  [Single-factor]  scale = {scale:.4f},  RMSE = {rmse:.1f} cm-1")
    print(f"  Literature B3LYP/def2-TZVP: ~0.9659 (NIST CCCBDB)")
    print(f"{'─' * 55}")

    if has_low:
        s_hi, rmse_hi, s_lo, rmse_lo = optimize_dual_scaling()
        all_dual_err = []
        for pair in PEAK_PAIRS:
            s = s_hi if pair[0] >= ZONE_BOUNDARY else s_lo
            all_dual_err.append(pair[0] * s - pair[1])
        rmse_dual = np.sqrt(np.mean(np.array(all_dual_err) ** 2))
        print(f"\n{'─' * 55}")
        print(f"  [Dual-zone]  High-freq (>2000): scale = {s_hi:.4f}, RMSE = {rmse_hi:.1f}")
        print(f"               Low-freq  (<2000): scale = {s_lo:.4f}, RMSE = {rmse_lo:.1f}")
        print(f"               Combined RMSE:     {rmse_dual:.1f} cm-1")
        print(f"{'─' * 55}")
        use_scale = lambda v: s_hi if v >= ZONE_BOUNDARY else s_lo
        title_scale = f"high: $\\times${s_hi:.4f}, low: $\\times${s_lo:.4f}, RMSE = {rmse_dual:.1f}"
        dft_wn_scaled = apply_dual_scale(dft_wn, s_hi, s_lo)
    else:
        # All peaks in one zone — use single factor
        use_scale = lambda v: scale
        title_scale = f"$\\times${scale:.4f}, RMSE = {rmse:.1f}"
        dft_wn_scaled = dft_wn * scale

    # Peak-by-peak table
    print(f"\n  {'Mode':<28} {'DFT':>7} {'Scale':>7} {'Scaled':>8} {'ATR':>7} {'Error':>7}")
    print(f"  {'─'*28} {'─'*7} {'─'*7} {'─'*8} {'─'*7} {'─'*7}")
    for i, mode in enumerate(modes):
        dft_v = PEAK_PAIRS[i, 0]
        exp_v = PEAK_PAIRS[i, 1]
        s = use_scale(dft_v)
        scl_v = dft_v * s
        err   = scl_v - exp_v
        print(f"  {mode:<28} {dft_v:>7.1f} {s:>7.4f} {scl_v:>8.1f} {exp_v:>7.0f} {err:>+7.1f}")

    # ── Plot ─────────────────────────────────────────────────────────────
    dft_norm = dft_eps / dft_eps.max()
    atr_norm = atr_y / atr_y.max()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), dpi=200,
                                    gridspec_kw={"height_ratios": [1, 1]})

    # Top: full overlay
    ax1.plot(atr_wn, atr_norm, color="#984EA3", linewidth=1.2,
             label="ATR Experimental", alpha=0.9)
    ax1.plot(dft_wn_scaled, dft_norm, color="#E41A1C", linewidth=1.2,
             label=f"DFT B3LYP/def2-TZVP (scaled)", alpha=0.8)
    ax1.fill_between(dft_wn_scaled, 0, dft_norm, color="#E41A1C", alpha=0.15)

    for i in range(len(PEAK_PAIRS)):
        ax1.axvline(PEAK_PAIRS[i, 1], color="gray", ls=":", lw=0.7, alpha=0.5)

    ax1.set_xlim(4000, 400)
    ax1.set_xlabel("Wavenumber (cm$^{-1}$)")
    ax1.set_ylabel("Normalized Intensity")
    ax1.set_title(f"DFT vs Experimental IR — {label}  ({title_scale} cm$^{{-1}}$)")
    ax1.legend(frameon=False, loc="upper left")

    # Bottom: annotated
    ax2.plot(atr_wn, atr_norm, color="#984EA3", linewidth=1.2,
             label="ATR Experimental", alpha=0.9)
    ax2.plot(dft_wn_scaled, dft_norm, color="#E41A1C", linewidth=1.2,
             label="DFT Scaled", alpha=0.8)
    ax2.fill_between(dft_wn_scaled, 0, dft_norm, color="#E41A1C", alpha=0.15)

    for i, mode in enumerate(modes):
        dft_v = PEAK_PAIRS[i, 0]
        exp_v = PEAK_PAIRS[i, 1]
        s = use_scale(dft_v)
        scl_v = dft_v * s
        err   = scl_v - exp_v
        y_at_exp = np.interp(exp_v, atr_wn, atr_norm)
        offsets = [25, 50, 75]
        ax2.annotate(
            f"{mode}\nDFT: {scl_v:.0f} | ATR: {exp_v:.0f}\n($\\Delta$ = {err:+.0f})",
            xy=(exp_v, y_at_exp),
            xytext=(0, offsets[i % len(offsets)]),
            textcoords="offset points",
            fontsize=8, ha="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="gray", lw=0.5),
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.8),
        )

    ax2.set_xlim(4000, 400)
    ax2.set_xlabel("Wavenumber (cm$^{-1}$)")
    ax2.set_ylabel("Normalized Intensity")
    ax2.set_title(f"Peak-by-Peak Mode Assignment — {label}")
    ax2.legend(frameon=False, loc="upper left")

    plt.tight_layout()
    path = os.path.join(FIG_DIR, f"dft_vs_atr_cyclohexane{fig_suffix}.png")
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"\n  [Saved] {path}")


def main():
    # ── Boat form ────────────────────────────────────────────────────────
    run_comparison(
        label="Cyclohexane (Boat Form)",
        dft_file="cyclohexane_IR_def2tzvp.txt",
        peak_pairs=PEAK_PAIRS_BOAT,
        modes=["C-H asym. stretch", "C-H sym. stretch", "CH2 scissoring",
               "Ring deformation", "Ring breathing"],
        fig_suffix="_boat",
    )

    # ── Chair form ───────────────────────────────────────────────────────
    run_comparison(
        label="Cyclohexane (Chair Form)",
        dft_file="cyclohexane_IR_def2tzvp_chairform.txt",
        peak_pairs=PEAK_PAIRS_CHAIR,
        modes=["C-H eq. stretch", "C-H ax. stretch", "C-H sym. stretch"],
        fig_suffix="_chair",
    )

    print("\n" + "=" * 65)
    print("  Both conformers compared. Done.")
    print("=" * 65)


if __name__ == "__main__":
    main()
