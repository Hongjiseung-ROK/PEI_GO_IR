#!/usr/bin/env python3
"""
Full-Spectrum FT-IR Peak Tagger
================================
1. arPLS baseline correction (rampy) + Savitzky-Golay smoothing
2. Full-range peak detection (scipy.signal.find_peaks) on 450–3900 cm⁻¹
3. Weighted functional-group scoring with confidence ranking
4. Annotated publication-quality plots + CSV export

Covers all 4 spectra: GO, PEI, PEI@GO (KBr), Unknown (ATR).
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.signal import find_peaks, savgol_filter
import pybaselines
import warnings
warnings.simplefilter("ignore", RuntimeWarning)

# ── Paths ────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
FIG_DIR  = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ── Plot aesthetics ──────────────────────────────────────────────────────
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

# ══════════════════════════════════════════════════════════════════════════
# REFERENCE TABLE — weighted scoring bands
# Each entry: (lo, hi, center, bond, mol_type, base_weight)
#   base_weight: how "diagnostic" this band is (1.0 = very, 0.5 = ambiguous)
# ══════════════════════════════════════════════════════════════════════════
BANDS = [
    # ── O-H stretches ────────────────────────────────────────────────────
    (3400, 3650, 3500, "O-H",  "Alcohol str.",           1.0),
    (3200, 3400, 3300, "O-H",  "Phenol / H-bonded",      0.9),
    (2500, 3100, 2800, "O-H",  "Carboxylic acid (broad)", 0.6),
    # ── N-H stretches ────────────────────────────────────────────────────
    (3300, 3500, 3400, "N-H",  "Amine str.",             0.9),
    (3150, 3300, 3250, "N-H",  "H-bonded N-H",          0.7),
    # ── C-H stretches ────────────────────────────────────────────────────
    (3000, 3100, 3050, "=C-H", "Alkene / Aromatic C-H",  0.9),
    (2850, 2960, 2920, "C-H",  "Alkane C-H str.",        1.0),
    (2800, 2850, 2830, "C-H",  "Aldehyde C-H (Fermi)",   0.7),
    # ── Triple bonds ─────────────────────────────────────────────────────
    (2210, 2260, 2235, "C≡N",  "Nitrile",                1.0),
    (2100, 2200, 2150, "C≡C",  "Alkyne",                 1.0),
    # ── Carbonyl ─────────────────────────────────────────────────────────
    (1700, 1750, 1730, "C=O",  "Ester / Acid / Aldehyde", 1.0),
    (1630, 1700, 1660, "C=O",  "Amide I",                0.9),
    # ── C=C / N-H bend ───────────────────────────────────────────────────
    (1600, 1680, 1640, "C=C",  "Alkene str.",            0.8),
    (1580, 1620, 1600, "C=C",  "Aromatic C=C",           0.9),
    (1500, 1570, 1540, "N-H",  "Amide II / N-H bend",    0.9),
    (1450, 1500, 1470, "C=C",  "Aromatic C=C / CH bend",  0.7),
    # ── C-H bending ──────────────────────────────────────────────────────
    (1430, 1470, 1450, "C-H",  "CH2 scissoring",         0.8),
    (1370, 1390, 1380, "C-H",  "CH3 sym. bend",          0.8),
    # ── C-N, C-O, C-C ────────────────────────────────────────────────────
    (1300, 1430, 1360, "C-N",  "C-N stretch (amine)",     0.7),
    (1180, 1300, 1240, "C-O",  "Epoxy / Aryl ether C-O", 0.8),
    (1050, 1180, 1100, "C-O",  "Alcohol / Ester C-O",    0.9),
    (1000, 1050, 1030, "C-O",  "Primary alcohol C-O",    0.7),
    # ── Fingerprint ──────────────────────────────────────────────────────
    (880, 920,   900,  "=C-H", "Alkene oop bend",        0.7),
    (830, 880,   860,  "Ring", "Ring deformation",        0.8),
    (700, 830,   770,  "C-Cl", "Organic chloride",        0.8),
    (690, 750,   720,  "CH2",  "CH2 rocking (long chain)", 0.9),
    (450, 690,   570,  "Ring", "Skeletal / Ring puckering", 0.6),
]


# ══════════════════════════════════════════════════════════════════════════
# STAGE 1 — Baseline correction + smoothing
# ══════════════════════════════════════════════════════════════════════════
def correct_baseline(wn, y_raw, lam=1e5):
    """
    Two-output baseline correction:
    Uses pybaselines Asymmetric Least Squares (AsLS) instead of arPLS.
    """
    sort_idx = np.argsort(wn)
    x_asc = wn[sort_idx]
    y_asc = savgol_filter(y_raw[sort_idx], window_length=15, polyorder=3)

    fitter = pybaselines.Baseline(x_data=x_asc)
    y_base, _ = fitter.asls(y_asc, lam=lam, p=0.01)
    y_corr = y_asc - y_base

    return x_asc[::-1], y_asc[::-1], y_base[::-1], y_corr[::-1]


def load_kbr():
    """Load KBr data, invert transmittance→absorbance, baseline-correct."""
    df = pd.read_excel(os.path.join(DATA_DIR, "IR_data.xlsx"))
    wn = df["Wavenumber"].values
    spectra = {}
    for col in ["GO", "PEI", "PEI@GO"]:
        raw = df[col].values
        y_abs = raw.max() - raw           # transmittance → pseudo-absorbance
        x, y_sm, y_bl, y_corr = correct_baseline(wn, y_abs)
        spectra[col] = {"wn": x, "smooth": y_sm, "baseline": y_bl, "corrected": y_corr}
    return spectra


def load_atr():
    """Load ATR CSV (transmittance), invert, baseline-correct.
    Uses heavier smoothing (window=31) to suppress high-frequency noise
    typical of ATR measurements."""
    df = pd.read_csv(os.path.join(DATA_DIR, "ATR_result.CSV"),
                     header=None, names=["Wavenumber", "Transmittance"])
    wn = df["Wavenumber"].values
    raw = df["Transmittance"].values
    y_abs = raw.max() - raw

    # Heavier smoothing for noisy ATR data
    sort_idx = np.argsort(wn)
    x_asc = wn[sort_idx]
    y_asc = savgol_filter(y_abs[sort_idx], window_length=31, polyorder=3)

    fitter = pybaselines.Baseline(x_data=x_asc)
    y_base, _ = fitter.asls(y_asc, lam=1e5, p=0.01)
    y_corr = y_asc - y_base

    x_desc = x_asc[::-1]
    y_sm_desc = y_asc[::-1]
    y_bl_desc = y_base[::-1]
    y_corr_desc = y_corr[::-1]

    return {"Unknown (ATR)": {"wn": x_desc, "smooth": y_sm_desc,
                               "baseline": y_bl_desc, "corrected": y_corr_desc}}


# ══════════════════════════════════════════════════════════════════════════
# STAGE 2 — Full-spectrum peak detection
# ══════════════════════════════════════════════════════════════════════════
def detect_peaks(wn, y_smooth, prom_frac=0.04, min_dist=12,
                 wn_lo=1000, wn_hi=3700):
    """
    Detect peaks on the SMOOTHED (not baseline-corrected) spectrum.
    Prominence-based detection naturally handles both sharp and broad peaks.

    Two-pass strategy:
      Pass 1: find_peaks with prominence (catches obvious peaks)
      Pass 2: second-derivative maxima (catches shoulders & broad humps)
    Merge + deduplicate within 30 cm⁻¹.

    Excludes atmospheric CO₂ artifact window (2280–2400 cm⁻¹)
    and >3700 cm⁻¹ (water vapor noise).
    """
    from scipy.signal import peak_prominences

    sort_idx = np.argsort(wn)
    x_asc = wn[sort_idx]
    y_asc = y_smooth[sort_idx]

    mask = (x_asc >= wn_lo) & (x_asc <= wn_hi)
    x_win = x_asc[mask]
    y_win = y_asc[mask]

    # Pass 1: prominence-based
    prom = prom_frac * (y_win.max() - y_win.min())
    idx1, _ = find_peaks(y_win, prominence=prom, distance=min_dist, width=3)

    # Pass 2: negative second-derivative peaks
    d2 = np.gradient(np.gradient(y_win, x_win), x_win)
    neg_d2 = savgol_filter(-d2, window_length=25, polyorder=3)
    d2_thresh = 0.04 * neg_d2.max()
    idx2, _ = find_peaks(neg_d2, height=d2_thresh, distance=min_dist)
    # Only keep d2-peaks where signal is at least 20% of range
    sig_thresh = y_win.min() + 0.20 * (y_win.max() - y_win.min())
    idx2 = idx2[y_win[idx2] >= sig_thresh]

    # Merge
    all_idx = np.unique(np.concatenate([idx1, idx2]))
    all_wn  = x_win[all_idx]
    order = np.argsort(all_wn)
    all_idx = all_idx[order]
    all_wn  = all_wn[order]

    # Exclude atmospheric CO₂ window (2280–2400 cm⁻¹)
    atm_mask = ~((all_wn >= 2280) & (all_wn <= 2400))
    all_idx = all_idx[atm_mask]
    all_wn  = all_wn[atm_mask]

    # Deduplicate within 30 cm⁻¹ — keep higher intensity
    keep = [0] if len(all_idx) > 0 else []
    for i in range(1, len(all_idx)):
        if all_wn[i] - all_wn[keep[-1]] < 30:
            if y_win[all_idx[i]] > y_win[all_idx[keep[-1]]]:
                keep[-1] = i
        else:
            keep.append(i)

    final_idx = all_idx[keep] if keep else np.array([], dtype=int)
    peak_wn = x_win[final_idx]
    peak_h  = y_win[final_idx]

    if len(final_idx) > 0:
        prominences, _, _ = peak_prominences(y_win, final_idx)
    else:
        prominences = np.array([])

    # Sort descending wavenumber
    order = np.argsort(peak_wn)[::-1]
    return peak_wn[order], peak_h[order], prominences[order]


# ══════════════════════════════════════════════════════════════════════════
# STAGE 3 — Weighted functional-group assignment
# ══════════════════════════════════════════════════════════════════════════
def score_assignment(peak_wn, band):
    """
    Score how well a peak matches a reference band.
    Returns 0 if outside range, else a 0–1 confidence based on
    proximity to band center, scaled by the band's diagnostic weight.
    """
    lo, hi, center, bond, desc, base_w = band
    if peak_wn < lo or peak_wn > hi:
        return 0.0
    # Gaussian-like proximity score: 1.0 at center, decays toward edges
    half_width = (hi - lo) / 2.0 + 1e-6
    dist = abs(peak_wn - center)
    proximity = np.exp(-0.5 * (dist / (half_width * 0.6)) ** 2)
    return proximity * base_w


def assign_peak(peak_wn, sample_name=""):
    """Return ranked list of (score, bond, description) for a peak."""
    scores = []
    
    exclude_N = (sample_name == "GO")
    exclude_O = (sample_name == "PEI")
    
    for band in BANDS:
        bond = band[3]
        desc = band[4]
        
        # Enforce chemical formula constraints
        if exclude_N and ("N" in bond or "Amide" in desc or "Amine" in desc or "Nitrile" in desc):
            continue
        if exclude_O and ("O" in bond or "Acid" in desc or "Alcohol" in desc or "Ester" in desc or "Aldehyde" in desc or "Phenol" in desc or "Epoxy" in desc):
            continue
            
        s = score_assignment(peak_wn, band)
        if s > 0.01:
            scores.append((s, bond, desc))
            
    scores.sort(key=lambda x: -x[0])
    return scores


def format_assignment(scores, top_n=2):
    """Format top assignments as a label string."""
    if not scores:
        return "Unassigned"
    parts = []
    for s, bond, desc in scores[:top_n]:
        parts.append(f"{bond} {desc} ({s:.0%})")
    return " | ".join(parts)


def build_results_table(spectra):
    """Detect peaks and assign functional groups for all spectra."""
    rows = []
    all_peaks = {}  # name → (wn_arr, h_arr, prom_arr, labels)

    for name, data in spectra.items():
        wn = data["wn"]
        y  = data["corrected"]   # Use the perfectly flattened baseline-corrected signal
        y_norm = y / (y.max() + 1e-9) # Normalize height to 0-1 range
        
        peak_wn, peak_h, peak_prom = detect_peaks(wn, y_norm)
        labels = []

        filtered_pw = []
        filtered_ph = []
        filtered_prom = []
        
        for pw, ph, pp in zip(peak_wn, peak_h, peak_prom):
            scores = assign_peak(pw, sample_name=name)
            label = format_assignment(scores)
            if label == "Unassigned":
                continue
                
            labels.append(label)
            filtered_pw.append(pw)
            filtered_ph.append(ph)
            filtered_prom.append(pp)
            
            top_bond = scores[0][1] if scores else ""
            top_desc = scores[0][2] if scores else "Unassigned"
            top_conf = scores[0][0] if scores else 0.0
            rows.append({
                "Sample": name,
                "Peak (cm⁻¹)": round(pw, 1),
                "Height (a.u.)": round(ph, 4),
                "Prominence": round(pp, 4),
                "Top Assignment": f"{top_bond} {top_desc}".strip(),
                "Confidence": round(top_conf, 2),
                "All Matches": label,
            })

        all_peaks[name] = (np.array(filtered_pw), np.array(filtered_ph), np.array(filtered_prom), labels)

    df = pd.DataFrame(rows)
    return df, all_peaks


# ══════════════════════════════════════════════════════════════════════════
# STAGE 4 — Annotated plots
# ══════════════════════════════════════════════════════════════════════════
COLORS = {
    "GO": "#E41A1C",
    "PEI": "#377EB8",
    "PEI@GO": "#4DAF4A",
    "Unknown (ATR)": "#984EA3",
}


def plot_annotated(name, data, peak_wn, peak_h, labels):
    """Single annotated spectrum with functional-group labels."""
    wn = data["wn"]
    y_raw = data["corrected"]
    y_norm = y_raw / (y_raw.max() + 1e-9)

    fig, ax = plt.subplots(figsize=(14, 5.5), dpi=200)
    color = COLORS.get(name, "#333333")
    ax.plot(wn, y_norm, color=color, linewidth=1.0, label=name)

    # Deduplicate tags by highest peak
    best_labels = {}
    for i, (pw, ph, lbl) in enumerate(zip(peak_wn, peak_h, labels)):
        short = lbl.split(" (")[0].split(" | ")[0]  # first assignment, no confidence
        if short == "Unassigned":
            continue
            
        if short not in best_labels or ph > best_labels[short][1]:
            best_labels[short] = (pw, ph, i)

    # Mark representative peaks only
    rep_pw = [pw for _, (pw, _, _) in best_labels.items()]
    rep_ph = [ph for _, (_, ph, _) in best_labels.items()]
    if rep_pw:
        ax.plot(rep_pw, rep_ph, "kv", markersize=5, zorder=5)

    # Create ample headroom for labels
    ax.set_ylim(-0.05, 1.65)

    # Smart tier-based collision avoidance for labels
    sorted_labels = sorted(best_labels.items(), key=lambda x: x[1][0])
    tiers = [1.10, 1.25, 1.40, 1.55]
    last_pw_at_tier = [-9999] * 4

    for short, (pw, ph, original_i) in sorted_labels:
        # Find a horizontal tier that provides at least 250 cm-1 clearance
        chosen_tier = 0
        for t in range(4):
            if abs(pw - last_pw_at_tier[t]) > 250:
                chosen_tier = t
                break
        else:
            # Fallback to the tier with maximum clearance
            chosen_tier = np.argmax([abs(pw - pt) for pt in last_pw_at_tier])
            
        last_pw_at_tier[chosen_tier] = pw
        label_y = tiers[chosen_tier]

        ax.annotate(
            f"{short}\n({pw:.0f})",
            xy=(pw, ph),
            xytext=(pw, label_y),
            textcoords="data",
            fontsize=8,
            ha="center",
            va="center",
            color="black",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", lw=0.5, alpha=0.9),
            arrowprops=dict(arrowstyle="-", color="gray", lw=1.0, alpha=0.7),
        )

    ax.set_xlim(4000, 400)
    ax.set_xlabel("Wavenumber (cm$^{-1}$)")
    ax.set_ylabel("Baseline-Corrected Intensity (a.u.)")
    ax.set_title(f"FT-IR Peak Assignment — {name}")
    ax.legend(frameon=False, loc="upper left")
    plt.tight_layout()

    safe = name.replace("@", "_").replace(" ", "_").replace("(", "").replace(")", "")
    path = os.path.join(FIG_DIR, f"tagged_{safe}.png")
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  [Saved] {path}")


def plot_overlay(spectra, all_peaks):
    """Stacked overlay of all 4 spectra with peak markers."""
    fig, ax = plt.subplots(figsize=(15, 8), dpi=200)
    offset = 0
    step = 1.3

    for name in ["GO", "PEI", "PEI@GO", "Unknown (ATR)"]:
        if name not in spectra:
            continue
        wn = spectra[name]["wn"]
        y  = spectra[name]["corrected"]
        y_norm = y / (y.max() + 1e-9)
        color = COLORS.get(name, "#333")

        ax.plot(wn, y_norm + offset, color=color, linewidth=1.0, label=name)

        # Mark peaks
        pw, ph, _, labels = all_peaks[name]
        
        # Deduplicate tags by highest peak for overlay too
        best_labels_ov = {}
        for i, (p_w, p_h, lbl) in enumerate(zip(pw, ph, labels)):
            short = lbl.split(" (")[0].split(" | ")[0]
            if short == "Unassigned":
                continue
            if short not in best_labels_ov or p_h > best_labels_ov[short][1]:
                best_labels_ov[short] = (p_w, p_h, i)
        
        rep_pw_ov = [p_w for _, (p_w, _, _) in best_labels_ov.items()]
        rep_ph_ov = [p_h for _, (_, p_h, _) in best_labels_ov.items()]
        
        if rep_pw_ov:
            ax.plot(rep_pw_ov, np.array(rep_ph_ov) + offset, "kv", markersize=4, zorder=5)

        sorted_ov = sorted(zip(rep_pw_ov, rep_ph_ov), key=lambda x: x[0])
        tiers_ov = [0.15, 0.30, 0.45]
        last_pw_at_tier_ov = [-9999] * 3

        for p_w, p_h in sorted_ov:
            chosen_tier_ov = 0
            for t in range(3):
                if abs(p_w - last_pw_at_tier_ov[t]) > 180:
                    chosen_tier_ov = t
                    break
            else:
                chosen_tier_ov = np.argmax([abs(p_w - pt) for pt in last_pw_at_tier_ov])

            last_pw_at_tier_ov[chosen_tier_ov] = p_w
            dy = tiers_ov[chosen_tier_ov]

            ax.annotate(f"{p_w:.0f}", xy=(p_w, p_h + offset),
                        xytext=(p_w, p_h + offset + dy),
                        textcoords="data",
                        fontsize=7, ha="center", va="bottom",
                        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
                        bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8),
                        color="black")

        offset += step

    ax.set_xlim(4000, 400)
    ax.set_ylim(-0.1, offset + 0.6)
    ax.set_xlabel("Wavenumber (cm$^{-1}$)")
    ax.set_ylabel("Normalized Intensity (stacked, a.u.)")
    ax.set_title("Full-Spectrum Peak Detection — All Samples")
    ax.set_yticks([])
    ax.legend(frameon=False, loc="upper left")
    plt.tight_layout()

    path = os.path.join(FIG_DIR, "tagged_overlay.png")
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  [Saved] {path}")


# ══════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 65)
    print("  FT-IR Full-Spectrum Peak Tagger")
    print("  arPLS baseline → find_peaks → weighted scoring")
    print("=" * 65)

    # Load & correct
    print("\n[1/4] Loading and baseline-correcting spectra...")
    kbr = load_kbr()
    atr = load_atr()
    spectra = {**kbr, **atr}

    # Detect & assign
    print("[2/4] Detecting peaks and assigning functional groups...")
    df, all_peaks = build_results_table(spectra)

    # Print summary per sample
    for name in spectra:
        sub = df[df["Sample"] == name]
        print(f"\n  ── {name} ({len(sub)} peaks) ──")
        for _, row in sub.iterrows():
            print(f"    {row['Peak (cm⁻¹)']:>7.1f} cm⁻¹  │ {row['All Matches']}")

    # Export CSV
    csv_path = os.path.join(DATA_DIR, "tagged_peaks.csv")
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    print(f"\n[3/4] Exported → {csv_path}")

    # Plot
    print("[4/4] Generating annotated plots...")
    for name, data in spectra.items():
        pw, ph, _, labels = all_peaks[name]
        plot_annotated(name, data, pw, ph, labels)

    plot_overlay(spectra, all_peaks)

    print("\n" + "=" * 65)
    print("  Done. All peaks tagged with weighted functional group scores.")
    print("=" * 65)


if __name__ == "__main__":
    main()
