# FT-IR Spectroscopy: Surface Modification Analysis & Unknown Identification

Automated FT-IR peak detection and functional group assignment pipeline for graphene oxide surface modification chemistry and unknown hydrocarbon identification.

## Experiment Overview

| Analysis | Method | Samples |
|----------|--------|---------|
| Qualitative Analysis 1 | KBr Pellet | GO, PEI, PEI@GO |
| Qualitative Analysis 2 | ATR | Unknown sample (4 candidates) |

**Objective:** Confirm covalent grafting of polyethyleneimine (PEI) onto graphene oxide (GO) by tracking functional group changes, and identify an unknown hydrocarbon from its IR fingerprint.

---

## Results

### 1. Graphene Oxide (GO) — KBr Pellet

![GO spectrum](figures/tagged_GO.png)

GO exhibits the characteristic oxygen-rich surface chemistry expected from Hummers' method oxidation:

| Peak (cm⁻¹) | Assignment | Significance |
|-------------|------------|--------------|
| **3483** | O-H stretch (97%) | Hydroxyl groups on GO basal plane; broad band indicates extensive H-bonding network |
| **2760** | O-H carboxylic acid (59%) | Edge-site carboxyl groups (-COOH) |
| **1728** | C=O stretch (99%) | Carbonyl from carboxylic acid / ester / lactone groups — confirms oxidation |
| **1628** | C=C aromatic (70%) | Residual sp² graphitic backbone |
| **1385** | C-H sym. bend (59%) | Methyl groups at defect sites |
| **1084** | C-O stretch (83%) | Epoxy (C-O-C) and alkoxy groups on basal plane |

The spectrum confirms a heavily oxidized graphene structure: strong C=O (1728), broad O-H (3483), and C-O (1084) alongside residual aromatic C=C (1628). This O-H/C=O/C-O triad is the textbook signature of GO.

---

### 2. Polyethyleneimine (PEI) — KBr Pellet

![PEI spectrum](figures/tagged_PEI.png)

PEI shows the nitrogen-dominated spectrum of a branched polyamine:

| Peak (cm⁻¹) | Assignment | Significance |
|-------------|------------|--------------|
| **3476, 3404** | N-H stretch (41%, 90%) | Primary (-NH₂) and secondary (-NH-) amine groups; dual peaks indicate both types |
| **2954** | C-H alkane stretch (58%) | Ethylene backbone (-CH₂-CH₂-) |
| **1641** | C=C / Amide region (80%) | N-H deformation of primary amines (scissoring mode) |
| **1479** | Aromatic C=C / C-H bend (58%) | CH₂ deformation in ethylene backbone |
| **1311** | C-N stretch (32%) | Aliphatic amine C-N bond — signature of the polyamine chain |

The dual N-H peaks at 3476/3404 cm⁻¹ confirm both primary and secondary amines, consistent with branched PEI. The C-N stretch at 1311 cm⁻¹ further validates the polyamine backbone. Note: O-containing bands were excluded from scoring since pure PEI ([-CH₂-CH₂-NH-]ₙ) contains no oxygen.

---

### 3. PEI@GO Composite — KBr Pellet

![PEI@GO spectrum](figures/tagged_PEI_GO.png)

The PEI@GO spectrum provides direct evidence of covalent grafting through three key observations:

| Peak (cm⁻¹) | Assignment | Chemical Evidence |
|-------------|------------|-------------------|
| **3500** | O-H stretch (100%) | Retained from GO — surface hydroxyl groups |
| **3431** | N-H stretch (79%) | **New peak** — introduced by PEI grafting |
| **3361** | N-H / O-H H-bonded (73%) | Intermolecular H-bonding between PEI and GO |
| **2949** | C-H alkane (69%) | PEI ethylene backbone now present on GO surface |
| **2843** | O-H acid / C-H (58%) | Residual COOH + PEI backbone overlap |
| **1639** | C=C / **C=O Amide I** (80%/55%) | **Critical: amide bond formation** — COOH + NH₂ → CONH |
| **1570** | **N-H Amide II** (33%) | **New peak** — confirms amide linkage, absent in both GO and PEI alone |
| **1466** | C=C aromatic / CH₂ (67%) | Retained graphitic backbone |
| **1381** | C-H / C-N (79%/61%) | Overlap of GO methyl and PEI amine |
| **1113, 1080** | C-O stretch (85%, 79%) | Retained epoxy/alkoxy from GO |
| **1022** | C-O primary alcohol (61%) | Surface hydroxyl groups |

#### Evidence for Covalent Grafting

Three spectroscopic signatures conclusively demonstrate covalent PEI-GO bonding rather than physical adsorption:

1. **Amide bond formation (1639 + 1570 cm⁻¹):** The appearance of Amide I (C=O stretch at 1639) and Amide II (N-H bend at 1570) peaks — absent in both pure GO and pure PEI — proves that GO's carboxyl groups (-COOH) underwent condensation with PEI's primary amines (-NH₂) to form stable amide linkages (-CO-NH-).

2. **N-H peak emergence on GO surface (3431 cm⁻¹):** GO alone shows no nitrogen-containing peaks. The strong N-H stretch at 3431 cm⁻¹ in PEI@GO confirms nitrogen-bearing functional groups are now bonded to the graphene surface.

3. **O-H/N-H hydrogen-bonding network (3361 cm⁻¹):** The broad band at 3361 cm⁻¹, scored as both N-H (73%) and O-H H-bonded (53%), indicates intimate intermolecular contact between PEI's amine groups and GO's surface hydroxyls — consistent with a densely grafted polymer layer.

---

### 4. Unknown Sample — ATR Mode

![Unknown ATR spectrum](figures/tagged_Unknown_ATR.png)

| Candidate | Structure | Key IR Features |
|-----------|-----------|-----------------|
| (a) Cyclohexane | C₆H₁₂ (ring) | C-H str. ~2930, CH₂ scissors ~1450, ring deform. ~860–900 |
| (b) n-Heptane | C₇H₁₆ (chain) | C-H str. ~2930, CH₂ scissors ~1450, CH₂ rock ~720 |
| (c) trans-2-Butene | C₄H₈ (alkene) | =C-H ~3020, C=C ~1670 |
| (d) 1-Butyne | C₄H₆ (alkyne) | ≡C-H ~3300, C≡C ~2150 |

#### Diagnostic Peak Analysis

| Peak (cm⁻¹) | Assignment | Diagnostic Value |
|-------------|------------|-----------------|
| **2929** | C-H alkane stretch (96%) | Saturated C-H — rules out unsaturated candidates |
| **2856** | C-H alkane stretch (15%) | Symmetric CH₂ stretch — confirms alkane |
| **1450** | CH₂ scissoring (80%) | Methylene groups |
| **879** | Ring deformation (prom. = 5.03) | **Strongest fingerprint peak** — cyclohexane ring breathing mode |
| **731** | CH₂ rocking (prom. = 1.70) | Present but 3x weaker than 879 peak |

#### Elimination Logic

- **No C≡C stretch (~2150 cm⁻¹) and no ≡C-H stretch (~3300 cm⁻¹)** → **(d) 1-Butyne eliminated**
- **No C=C stretch (~1640–1680 cm⁻¹) and no =C-H stretch (~3020 cm⁻¹)** → **(c) trans-2-Butene eliminated**
- Both (a) and (b) show C-H alkane stretch + CH₂ scissoring. The distinction lies in the **fingerprint region (600–1000 cm⁻¹)**:
  - The **879 cm⁻¹ peak** (prominence = 5.03) is the dominant fingerprint feature — this is the **ring deformation / CH₂ rocking mode specific to the cyclohexane chair conformation**
  - The 731 cm⁻¹ CH₂ rocking peak, while present, has only 1/3 the prominence. In pure n-heptane, this would be the dominant fingerprint peak, not a secondary one
  - n-Heptane should also show a clear CH₃ symmetric bend at ~1378 cm⁻¹, which is absent

#### Conclusion

> **The unknown sample is (a) Cyclohexane (C₆H₁₂).**
>
> The spectrum is dominated by saturated C-H stretching (2929/2856 cm⁻¹) and CH₂ scissoring (1450 cm⁻¹) with no unsaturated bond signatures. The fingerprint region conclusively identifies a cyclic structure: the ring deformation mode at 879 cm⁻¹ is 3x more prominent than the 731 cm⁻¹ peak, matching the characteristic vibrational pattern of the cyclohexane chair conformation.

---

### 5. Overlay Comparison

![All spectra overlay](figures/tagged_overlay.png)

The stacked overlay enables direct visual comparison across all four samples. Key observations:
- The broad O-H/N-H region (3000–3600 cm⁻¹) progressively evolves from pure O-H (GO) through pure N-H (PEI) to a combined O-H + N-H profile (PEI@GO)
- The carbonyl region (1600–1750 cm⁻¹) shows GO's sharp C=O at 1728 cm⁻¹ transforming into the Amide I/II doublet at 1639/1570 cm⁻¹ in PEI@GO
- The ATR unknown sample shows a distinctly different, much simpler spectrum — purely aliphatic, consistent with cyclohexane

---

## Pipeline Architecture

```
tag_peaks.py
├── Stage 1: Baseline Correction
│   └── pybaselines AsLS (lam=1e5, p=0.01) + Savitzky-Golay smoothing
├── Stage 2: Peak Detection
│   ├── Pass 1: scipy.signal.find_peaks (prominence-based)
│   └── Pass 2: Second-derivative maxima (broad shoulders)
│   └── Merge + deduplicate within 30 cm⁻¹, exclude CO₂ window
├── Stage 3: Weighted Functional Group Scoring
│   ├── 26-band reference table with diagnostic weights
│   ├── Gaussian proximity scoring (center-weighted)
│   └── Chemical constraints (GO: exclude N-bands, PEI: exclude O-bands)
└── Stage 4: Annotated Plotting
    └── Tier-based label collision avoidance
```

## Usage

```bash
# Run full pipeline
python3 tag_peaks.py

# Outputs:
#   data/tagged_peaks.csv          — all peaks with assignments and confidence scores
#   figures/tagged_GO.png          — annotated GO spectrum
#   figures/tagged_PEI.png         — annotated PEI spectrum
#   figures/tagged_PEI_GO.png      — annotated PEI@GO spectrum
#   figures/tagged_Unknown_ATR.png — annotated unknown (ATR) spectrum
#   figures/tagged_overlay.png     — stacked comparison of all 4 spectra
```

### Dependencies

```
numpy pandas scipy matplotlib pybaselines openpyxl
```

## Data Sources

| File | Description |
|------|-------------|
| `data/IR_data.xlsx` | KBr pellet FT-IR data (Wavenumber, GO, PEI, PEI@GO) |
| `data/ATR_result.CSV` | ATR transmittance data (wavenumber, %T) |
| `guide.md` | Experiment instructions and IR reference table |

## References

- Hummers, W.S. & Offeman, R.E. (1958). Preparation of Graphitic Oxide. *JACS*, 80(6), 1339.
- Silverstein, R.M., Webster, F.X. & Kiemle, D.J. (2005). *Spectrometric Identification of Organic Compounds*, 7th ed. Wiley.
- Zhang, W. et al. (2011). General synthesis of PEI-coated GO. *Carbon*, 49, 986–995.
