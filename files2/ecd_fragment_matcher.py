"""
ECD Fragment Ion Matcher for Cyclic Peptides with Known Disulfide Bonds
=======================================================================
Given a peptide sequence and a KNOWN disulfide bond connectivity pattern,
this tool generates all theoretical c/z· fragment ions and optionally
matches them against experimental MS/MS data.

Usage:
    python ecd_fragment_matcher.py

References:
    Zubarev, R. A., et al. (1998). J. Am. Chem. Soc., 120(13), 3265–3266.
"""

import itertools
import pandas as pd
import os
import sys
from typing import Optional

# ========================= USER CONFIGURATION =========================
SEQUENCE = "GCNILQPYWGCGRDFECLEECLMDSQYYQ"

# Known disulfide bonds as list of (Cys_pos1, Cys_pos2) tuples (1-based)
# Example: C1-C11, C14-C18, C23-C29 (positions in the sequence above)
DISULFIDE_BONDS = [
    (2, 11),   # C2  — C11
    (17, 21),  # C17 — C21
]

# Optional: path to a text file with experimental m/z values (one per line)
# Set to None to skip matching
EXPERIMENTAL_MASSES_FILE = None  # e.g., "experimental_masses.txt"

# Matching tolerance
TOLERANCE = 0.02       # value in Da or ppm
TOLERANCE_UNIT = "Da"  # "Da" or "ppm"

# Output filename
OUTPUT_FILENAME = "ECD_Fragment_Matcher_Results.xlsx"
# ======================================================================


# ---- Monoisotopic mass constants ----
AA_MASS = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}

PROTON   = 1.00728
H_ATOM   = 1.00783
NH3      = 17.02655
OH       = 17.00274
NH       = 15.01090
SS_LOSS  = 2 * H_ATOM   # mass lost per intact S-S bond (2 × H)


# ---- Validation ----

def validate_sequence(seq: str) -> str:
    """Uppercase and verify all residues are in AA_MASS."""
    seq = seq.upper().strip()
    unknown = [aa for aa in seq if aa not in AA_MASS]
    if unknown:
        raise ValueError(f"Unknown amino acid(s) in sequence: {set(unknown)}")
    return seq


def validate_bonds(seq: str, bonds: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """
    Validate disulfide bond positions:
    - Positions must be within sequence length
    - Positions must correspond to Cys residues
    - No position appears more than once
    """
    n = len(seq)
    cys_positions = {i + 1 for i, aa in enumerate(seq) if aa == 'C'}
    used = set()
    validated = []
    for a, b in bonds:
        a, b = min(a, b), max(a, b)
        if a < 1 or b > n:
            raise ValueError(f"Bond position out of range: C{a}-C{b} (sequence length {n})")
        if a not in cys_positions:
            raise ValueError(f"Position {a} is not a Cys residue (found '{seq[a-1]}')")
        if b not in cys_positions:
            raise ValueError(f"Position {b} is not a Cys residue (found '{seq[b-1]}')")
        if a in used or b in used:
            raise ValueError(f"Cys position used in multiple bonds: C{a} or C{b}")
        used |= {a, b}
        validated.append((a, b))
    return validated


# ---- Core fragmentation logic ----

def is_bridged(cleavage_site: int, intact_bonds: list[tuple[int, int]]) -> tuple[bool, Optional[str]]:
    """
    Determine whether backbone cleavage at `cleavage_site` (between residue
    cleavage_site and cleavage_site+1) is trapped inside an intact S-S loop.

    A cleavage point i is bridged by bond (a, b) if:  a <= i < b

    Returns:
        (True, "Ca-Cb")  if bridged
        (False, None)    if not bridged
    """
    for a, b in intact_bonds:
        start, end = min(a, b), max(a, b)
        if start <= cleavage_site < end:
            return True, f"C{start}-C{end}"
    return False, None


def count_internal_ss(start_pos: int, end_pos: int,
                      intact_bonds: list[tuple[int, int]]) -> int:
    """
    Count how many intact S-S bonds are fully contained within the residue
    range [start_pos, end_pos] (1-based, inclusive).
    """
    return sum(
        1 for a, b in intact_bonds
        if start_pos <= min(a, b) and max(a, b) <= end_pos
    )


def fragment_mass(frag_seq: str, n_ss: int, ion_type: str) -> float:
    """
    Calculate MH+ mass for a c- or z·-ion fragment.

    c-ion:  Σ(residue masses) + NH3 + H(N-term) + H+(proton)  − n_ss × SS_LOSS
    z·-ion: Σ(residue masses) + OH  − NH + H(C-term) + H+(proton) − n_ss × SS_LOSS
    """
    base = sum(AA_MASS[aa] for aa in frag_seq)
    if ion_type == 'c':
        mass = base + NH3 + H_ATOM + PROTON
    else:  # z·
        mass = base + OH - NH + H_ATOM + PROTON
    mass -= n_ss * SS_LOSS
    return mass


# ---- Matching ----

def within_tolerance(theo: float, obs: float,
                     tol: float, unit: str) -> bool:
    if unit.lower() == 'ppm':
        return abs((obs - theo) / theo * 1e6) <= tol
    return abs(obs - theo) <= tol


def find_best_match(theo: float, obs_masses: list[float],
                    tol: float, unit: str) -> Optional[tuple[float, float]]:
    """Return (matched_obs_mass, delta_Da) or None."""
    best_obs, best_delta = None, float('inf')
    for obs in obs_masses:
        delta = abs(obs - theo)
        if within_tolerance(theo, obs, tol, unit) and delta < best_delta:
            best_obs, best_delta = obs, obs - theo
    return (best_obs, best_delta) if best_obs is not None else None


# ---- Main analysis ----

def generate_fragments(seq: str,
                       bonds: list[tuple[int, int]],
                       obs_masses: Optional[list[float]] = None,
                       tol: float = 0.02,
                       tol_unit: str = "Da") -> list[dict]:
    """
    Generate all theoretical c- and z·-ions for `seq` given known `bonds`.

    For each backbone cleavage position i (between residue i and i+1):
      - Compute c_i  (N-terminal fragment, length i)
      - Compute z_(n-i)  (C-terminal fragment, length n-i)
      - Check bridging; if bridged → SILENT, else calculate MH+ mass
      - Optionally match against obs_masses

    Returns a list of dicts (one per ion).
    """
    n = len(seq)
    results = []

    connectivity_str = "; ".join(f"C{a}-C{b}" for a, b in bonds)

    for i in range(1, n):
        c_frag = seq[:i]
        z_frag = seq[i:]
        z_pos  = n - i

        bridged, bridge_by = is_bridged(i, bonds)

        # --- c-ion ---
        c_ss    = count_internal_ss(1, i, bonds)
        c_mass  = None if bridged else fragment_mass(c_frag, c_ss, 'c')
        c_match = None
        c_delta = None
        if c_mass is not None and obs_masses:
            hit = find_best_match(c_mass, obs_masses, tol, tol_unit)
            if hit:
                c_match, c_delta = hit

        results.append({
            "Ion":              f"c{i}",
            "Ion_Type":         "c",
            "Position":         i,
            "Fragment_Sequence": c_frag,
            "Connectivity":     connectivity_str,
            "Theo_Mass_MH+":    round(c_mass, 4) if c_mass else "N/A",
            "Status":           "SILENT" if bridged else "Observed",
            "Note":             (f"Bridged by {bridge_by}" if bridged
                                 else (f"Contains {c_ss} S-S bond(s)" if c_ss else "")),
            "Internal_SS":      c_ss if not bridged else "",
            "Matched_Obs_Mass": round(c_match, 4) if c_match else "",
            "Delta_Da":         round(c_delta, 4) if c_delta is not None else "",
        })

        # --- z·-ion ---
        z_ss    = count_internal_ss(i + 1, n, bonds)
        z_mass  = None if bridged else fragment_mass(z_frag, z_ss, 'z')
        z_match = None
        z_delta = None
        if z_mass is not None and obs_masses:
            hit = find_best_match(z_mass, obs_masses, tol, tol_unit)
            if hit:
                z_match, z_delta = hit

        results.append({
            "Ion":              f"z{z_pos}",
            "Ion_Type":         "z·",
            "Position":         z_pos,
            "Fragment_Sequence": z_frag,
            "Connectivity":     connectivity_str,
            "Theo_Mass_MH+":    round(z_mass, 4) if z_mass else "N/A",
            "Status":           "SILENT" if bridged else "Observed",
            "Note":             (f"Bridged by {bridge_by}" if bridged
                                 else (f"Contains {z_ss} S-S bond(s)" if z_ss else "")),
            "Internal_SS":      z_ss if not bridged else "",
            "Matched_Obs_Mass": round(z_match, 4) if z_match else "",
            "Delta_Da":         round(z_delta, 4) if z_delta is not None else "",
        })

    return results


def load_experimental_masses(filepath: str) -> list[float]:
    masses = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            for part in line.replace(',', ' ').split():
                try:
                    masses.append(float(part))
                except ValueError:
                    pass
    return masses


def export_excel(df: pd.DataFrame, filename: str) -> None:
    """Write results to a formatted Excel file."""
    with pd.ExcelWriter(filename, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name='Fragment Ions')

        ws = writer.sheets['Fragment Ions']
        from openpyxl.styles import Font, PatternFill, Alignment

        # Header style
        header_fill = PatternFill('solid', start_color='2C2C2A', end_color='2C2C2A')
        for cell in ws[1]:
            cell.font  = Font(bold=True, color='FFFFFF', name='Arial', size=10)
            cell.fill  = header_fill
            cell.alignment = Alignment(horizontal='center')

        # Column widths
        col_widths = {
            'A': 9, 'B': 8, 'C': 10, 'D': 32,
            'E': 28, 'F': 16, 'G': 11, 'H': 30,
            'I': 12, 'J': 18, 'K': 12,
        }
        for col, w in col_widths.items():
            ws.column_dimensions[col].width = w

        # Conditional row coloring
        silent_fill  = PatternFill('solid', start_color='FCEBEB', end_color='FCEBEB')
        match_fill   = PatternFill('solid', start_color='EAF3DE', end_color='EAF3DE')

        for row in ws.iter_rows(min_row=2):
            status = row[6].value  # column G = Status
            matched = row[9].value  # column J = Matched_Obs_Mass
            for cell in row:
                cell.font = Font(name='Arial', size=10)
                if status == 'SILENT':
                    cell.fill = silent_fill
                elif matched:
                    cell.fill = match_fill

        ws.freeze_panes = 'A2'


def print_summary(df: pd.DataFrame, obs_masses: Optional[list[float]],
                  tol: float, tol_unit: str) -> None:
    n_total   = len(df)
    n_obs     = (df['Status'] == 'Observed').sum()
    n_silent  = (df['Status'] == 'SILENT').sum()
    n_matched = (df['Matched_Obs_Mass'] != '').sum()

    print("\n" + "="*55)
    print("  ECD Fragment Ion Matcher  —  Summary")
    print("="*55)
    print(f"  Total ions generated : {n_total}")
    print(f"  Observed (detectable): {n_obs}")
    print(f"  Silent (bridged)     : {n_silent}")
    if obs_masses:
        print(f"  Experimental masses  : {len(obs_masses)}")
        print(f"  Matched (≤ {tol} {tol_unit})    : {n_matched}")
        pct = n_matched / n_obs * 100 if n_obs else 0
        print(f"  Match rate           : {pct:.1f}% of observable ions")
    print("="*55 + "\n")


# ---- Entry point ----

def main():
    print("ECD Fragment Ion Matcher")
    print("------------------------")

    seq = validate_sequence(SEQUENCE)
    bonds = validate_bonds(seq, DISULFIDE_BONDS)

    cys_pos = [i + 1 for i, aa in enumerate(seq) if aa == 'C']
    print(f"Sequence length : {len(seq)}")
    print(f"Cys positions   : {cys_pos}")
    print(f"Disulfide bonds : {', '.join(f'C{a}-C{b}' for a, b in bonds)}")

    obs_masses = None
    if EXPERIMENTAL_MASSES_FILE:
        if not os.path.exists(EXPERIMENTAL_MASSES_FILE):
            print(f"Warning: experimental masses file '{EXPERIMENTAL_MASSES_FILE}' not found. Skipping matching.")
        else:
            obs_masses = load_experimental_masses(EXPERIMENTAL_MASSES_FILE)
            print(f"Loaded {len(obs_masses)} experimental m/z values from '{EXPERIMENTAL_MASSES_FILE}'")

    print("\nGenerating fragment ions...")
    rows = generate_fragments(seq, bonds, obs_masses, TOLERANCE, TOLERANCE_UNIT)
    df = pd.DataFrame(rows)

    print_summary(df, obs_masses, TOLERANCE, TOLERANCE_UNIT)

    try:
        export_excel(df, OUTPUT_FILENAME)
        print(f"Results saved to: {os.path.abspath(OUTPUT_FILENAME)}")
    except Exception as e:
        print(f"Excel export failed: {e}")
        csv_name = OUTPUT_FILENAME.replace('.xlsx', '.csv')
        df.to_csv(csv_name, index=False)
        print(f"Saved as CSV instead: {os.path.abspath(csv_name)}")


if __name__ == "__main__":
    main()
