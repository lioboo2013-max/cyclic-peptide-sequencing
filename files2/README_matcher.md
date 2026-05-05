# ECD Fragment Ion Matcher

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

A focused ECD (Electron Capture Dissociation) mass spectrometry tool for cyclic peptides where the disulfide bond connectivity is **already known**.

Given a peptide sequence and a set of known disulfide bonds, this tool:

1. Generates all theoretical **c-ions** and **z·-ions**
2. Identifies **"SILENT" fragments** — backbone cleavages trapped inside intact S-S loops
3. Computes precise **MH⁺ masses** with disulfide-induced mass corrections (−2.016 Da per intact bond)
4. **Matches** theoretical masses against experimental MS/MS data within a user-defined tolerance

This is the successor to the exhaustive combinatorial analyzer (`peptide disulfide bond ecd.py`). Instead of enumerating all 15 possible S-S topologies, this tool takes the known topology as direct input — ideal for **validation and experimental comparison** once connectivity has been hypothesized.

---

## Key Difference from the Combinatorial Analyzer

| Feature | Combinatorial Analyzer | **This Tool** |
|---------|----------------------|--------------|
| Disulfide bonds | Exhaustively enumerates all patterns | **User-specified known pattern** |
| ECD states | All power-set scenarios per pattern | Single known state |
| Output rows | ~7,500 (15 patterns × 64 states) | ~2 × (n−1) ions |
| Use case | Hypothesis generation | **Experimental validation & matching** |

---

## Scientific Background

### c-ion mass

```
MH⁺ = Σ(residue masses) + NH₃ + H(N-term) + H⁺  −  n × 2.0157
                                                       └── per intact S-S
```

### z·-ion mass

```
MH⁺ = Σ(residue masses) + OH − NH + H(C-term) + H⁺  −  n × 2.0157
```

### Silent fragment rule

A cleavage at position *i* is **SILENT** if any intact bond (a, b) satisfies `a ≤ i < b`.  
The backbone breaks, but the fragment is held in place by the disulfide bridge and cannot be detected.

---

## Installation

```bash
git clone https://github.com/lioboo2013-max/cyclic-peptide-sequencing.git
cd cyclic-peptide-sequencing
pip install -r requirements.txt
```

---

## Quick Start

Edit the configuration block at the top of `ecd_fragment_matcher.py`:

```python
SEQUENCE = "GCNILQPYWGCGRDFECLEECLMDSQYYQ"

DISULFIDE_BONDS = [
    (1, 11),   # C1  — C11
    (14, 18),  # C14 — C18
    (23, 29),  # C23 — C29
]

EXPERIMENTAL_MASSES_FILE = "my_data.txt"   # or None to skip matching
TOLERANCE      = 0.02
TOLERANCE_UNIT = "Da"   # or "ppm"
OUTPUT_FILENAME = "ECD_Fragment_Matcher_Results.xlsx"
```

Then run:

```bash
python ecd_fragment_matcher.py
```

### Console output

```
ECD Fragment Ion Matcher
------------------------
Sequence length : 30
Cys positions   : [1, 11, 14, 18, 23, 29]
Disulfide bonds : C1-C11, C14-C18, C23-C29

Generating fragment ions...

=======================================================
  ECD Fragment Ion Matcher  —  Summary
=======================================================
  Total ions generated : 58
  Observed (detectable): 34
  Silent (bridged)     : 24
  Experimental masses  : 142
  Matched (≤ 0.02 Da)  : 28
  Match rate           : 82.4% of observable ions
=======================================================

Results saved to: /path/to/ECD_Fragment_Matcher_Results.xlsx
```

---

## Experimental Masses File Format

Plain text, one m/z value per line. Lines starting with `#` are ignored.
Comma-separated values on a single line are also accepted.

```
# Sample ECD-MS/MS data
845.3412
1102.5638
1544.8900
2031.1244
```

---

## Excel Output

Each row is one theoretical fragment ion.

| Column | Description |
|--------|-------------|
| `Ion` | Ion label, e.g. `c15` or `z12` |
| `Ion_Type` | `c` or `z·` |
| `Position` | Fragment length (number of residues) |
| `Fragment_Sequence` | Amino acid sequence of the fragment |
| `Connectivity` | Disulfide bond pattern used |
| `Theo_Mass_MH+` | Predicted MH⁺ mass (4 decimal places) |
| `Status` | `Observed` or `SILENT` |
| `Note` | Why silent, or how many internal S-S bonds |
| `Internal_SS` | Count of S-S bonds fully inside the fragment |
| `Matched_Obs_Mass` | Best-matching experimental mass (if any) |
| `Delta_Da` | Mass difference: obs − theo |

Row colors: **red tint** = SILENT, **green tint** = matched to experimental mass.

---

## Programmatic Use

```python
from ecd_fragment_matcher import (
    validate_sequence, validate_bonds,
    generate_fragments
)
import pandas as pd

seq   = validate_sequence("GCNILQPYWGCGRDFECLEECLMDSQYYQ")
bonds = validate_bonds(seq, [(1,11),(14,18),(23,29)])

obs   = [845.34, 1102.56, 1544.89]   # your MS data
rows  = generate_fragments(seq, bonds, obs_masses=obs, tol=0.02, tol_unit="Da")
df    = pd.DataFrame(rows)

print(df[df['Matched_Obs_Mass'] != ''])
```

---

## Physical Constants (monoisotopic)

| Constant | Value (Da) |
|----------|-----------|
| Proton H⁺ | 1.00728 |
| Hydrogen atom | 1.00783 |
| NH₃ (c-ion N-term) | 17.02655 |
| OH (z-ion C-term) | 17.00274 |
| NH radical (z-ion) | 15.01090 |
| S-S loss per bond | 2.01566 |

---

## Project Structure

```
cyclic-peptide-sequencing/
├── README.md                              # Original combinatorial analyzer docs
├── README_EN.md
├── README_matcher.md                      # This file
├── requirements.txt
├── ecd_fragment_matcher.py                # ← New: known-connectivity matcher
└── dicyclic peptide/
    ├── peptide disulfide bond ecd.py      # Original combinatorial analyzer (ZH)
    ├── peptide_disulfide_bond_ecd_en.py   # Original combinatorial analyzer (EN)
    └── peptide.py
```

---

## License

MIT — see [LICENSE](LICENSE).
