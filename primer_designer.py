"""
primer_designer.py – Design Golden Gate primers for each fragment.

Incorporates and fixes the logic from the original primer.py, but works
directly with Fragment objects from fragment_splitter.py — no NEB metadata
CSV required when calling from the pipeline.

For each fragment the primers are:
  Forward: [13-nt BsmBI prefix] [4-nt 5' overhang] [gene-specific core]
  Reverse: [13-nt BsmBI prefix] [4-nt 3' overhang rc] [gene-specific core rc]

BsmBI (CGTCTC) in the prefix cuts 1 nt downstream, placing its cut exactly at the
start of the 4-nt overhang — no extra bases are introduced into the insert.

Special cases (pJUMP):
  Part A, fragment 1   forward → ATATC (EcoRV) inserted between oh5 and gene core
  Part E, last fragment reverse → ATATC (EcoRV) inserted between revcomp(oh3) and gene core
"""

from __future__ import annotations
from typing import List, Dict, Optional, Tuple
import pandas as pd

from fragment_splitter import Fragment, Part

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
TM_TARGET        = 60.0
TM_MIN_BACKTRACK = 50.0
MIN_LEN          = 10
MAX_EXTRA_FOR_GC = 15

# Standard prefixes — BsmBI-only design.
#
# Layout (forward):  [6-nt spacer] [CGTCTC] [1-nt spacer G]   = 13 nt
#   BsmBI (CGTCTC) sits at positions 6–11; it cuts 1 nt downstream, i.e.
#   between positions 12 and 13.  oh5 therefore starts at position 13 and
#   is exposed directly as the 4-nt 5′ overhang — no extra bases.
#
# Layout (reverse):  [6-nt spacer] [CGTCTC] [1-nt spacer A]   = 13 nt
#   Same geometry on the bottom strand; revcomp(oh3) starts at position 13.
#
# The previous design included a trailing BsaI site (GGTCTCA, 7 nt) which
# placed oh5/oh3 at position 19/20.  With BsmBI cutting at position 12/13,
# that left 6 "extra" bases in every insert — the BsmBI/BsaI overlap flaw.
FWD_PREFIX         = "GGCTACCGTCTCG"
REV_PREFIX         = "GGCTACCGTCTCA"

# Extra EcoRV (GATATC) context added between the overhang and gene-core
# for the outermost primers: Part A frag-1 FWD and Part E last-frag REV.
# The EcoRV site spans the last base of the overhang (G) + ATATC, giving
# GGAG + ATATC = GGAGATATC → EcoRV at positions 3–8.
# This lets the fully assembled gene be excised with EcoRV if needed.
ECORV_EXTRA        = "ATATC"

# Legacy aliases (kept so existing imports don't break)
FWD_PREFIX_A1      = FWD_PREFIX   # actual ATATC insert is added at primer-build time
REV_PREFIX_E_LAST  = REV_PREFIX   # same


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

def wallace_tm(seq: str) -> float:
    """Wallace rule: 2*(A+T) + 4*(G+C)."""
    seq = seq.upper()
    return 2.0 * (seq.count("A") + seq.count("T")) + 4.0 * (seq.count("G") + seq.count("C"))


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


# ---------------------------------------------------------------------------
# Core primer design logic
# ---------------------------------------------------------------------------

def _design_core(target: str) -> Tuple[str, float, bool, bool]:
    """
    Grow a primer from the 5' end of `target`:
      1. Extend until Wallace Tm >= TM_TARGET
      2. Try to end on G/C (extend up to MAX_EXTRA_FOR_GC extra bases)
      3. Backtrack to last G/C if still no clamp, keeping Tm >= TM_MIN_BACKTRACK
      4. If all else fails, return best-Tm candidate and flag as sub-ideal

    Returns: (core_seq, tm, gc_clamp_ok, tm_ok)
    """
    seq = target.upper()
    n = len(seq)

    if n < MIN_LEN:
        core = seq
        tm = wallace_tm(core)
        gc_ok = core[-1] in "GC" if core else False
        return core, tm, gc_ok, tm >= TM_MIN_BACKTRACK

    # Step 1: grow to first Tm >= TM_TARGET
    baseline_len = n
    baseline_tm  = 0.0
    for L in range(MIN_LEN, n + 1):
        tm = wallace_tm(seq[:L])
        if tm >= TM_TARGET:
            baseline_len = L
            baseline_tm  = tm
            break
    else:
        baseline_tm = wallace_tm(seq[:baseline_len])

    baseline = seq[:baseline_len]

    # Step 2: extend up to MAX_EXTRA_FOR_GC to land on G/C
    candidate     = baseline
    candidate_tm  = baseline_tm
    for extra in range(1, MAX_EXTRA_FOR_GC + 1):
        if candidate[-1] in "GC":
            break
        if baseline_len + extra > n:
            break
        candidate    = seq[:baseline_len + extra]
        candidate_tm = wallace_tm(candidate)

    if candidate[-1] in "GC":
        return candidate, candidate_tm, True, candidate_tm >= TM_MIN_BACKTRACK

    # Step 3: backtrack to nearest G/C with Tm >= TM_MIN_BACKTRACK
    best_back = baseline
    best_back_tm = baseline_tm
    for i in range(baseline_len - 1, MIN_LEN - 2, -1):
        if seq[i] in "GC":
            back_seq = seq[:i + 1]
            back_tm  = wallace_tm(back_seq)
            if back_tm >= TM_MIN_BACKTRACK:
                return back_seq, back_tm, True, True
            if back_tm > best_back_tm:
                best_back = back_seq
                best_back_tm = back_tm

    # Step 4: pick best overall
    candidates = [
        (baseline,    baseline_tm,    baseline[-1]    in "GC"),
        (candidate,   candidate_tm,   candidate[-1]   in "GC"),
        (best_back,   best_back_tm,   best_back[-1]   in "GC"),
    ]
    best_seq, best_tm, best_gc = max(candidates, key=lambda t: t[1])
    return best_seq, best_tm, best_gc, best_tm >= TM_MIN_BACKTRACK


# ---------------------------------------------------------------------------
# Per-fragment primer design
# ---------------------------------------------------------------------------

def design_fragment_primers(
    frag: Fragment,
    is_first_of_part_a: bool = False,
    is_last_of_part_e: bool = False,
) -> Dict:
    """
    Design forward and reverse primers for a single Fragment.

    In Golden Gate the BsaI cut exposes the 4-nt overhang, so the primer's
    gene-specific portion begins AFTER the overhang.  However, to produce a
    PCR product that still contains the overhang bases (they are within the
    primer), we include the overhang in the core sequence.

    Forward primer:  prefix + overhang_5 + fwd_core(from base 4 of fragment)
    Reverse primer:  prefix + revcomp(overhang_3) + rev_core(from last base-4)
    """
    seq = frag.sequence.upper()
    oh5 = frag.overhang_5.upper()
    oh3 = frag.overhang_3.upper()

    # The overhang bases are the first/last 4 nt of the insert.
    # The gene-specific core grows from the sequence *past* the overhang.
    fwd_target = seq                    # forward core from 5' end
    rev_target = revcomp(seq)          # reverse core from 3' end (via rc)

    fwd_core_raw, fwd_tm, fwd_gc_ok, fwd_tm_ok = _design_core(fwd_target)
    rev_core_raw, rev_tm, rev_gc_ok, rev_tm_ok = _design_core(rev_target)

    # Full primers: FWD_PREFIX + oh5 [+ ATATC if outer FWD] + gene core
    #               REV_PREFIX + revcomp(oh3) [+ ATATC if outer REV] + gene core
    #
    # The ECORV_EXTRA ("ATATC") is inserted between the overhang and gene core
    # only for the two outermost primers.  Together with the last base of the
    # GGAG / CGCT overhang it forms the EcoRV recognition site GATATC.
    fwd_ecorv = ECORV_EXTRA if is_first_of_part_a else ""
    rev_ecorv = ECORV_EXTRA if is_last_of_part_e  else ""

    fwd_full = FWD_PREFIX + oh5 + fwd_ecorv + fwd_core_raw
    rev_full = REV_PREFIX + revcomp(oh3) + rev_ecorv + rev_core_raw

    return {
        "name": frag.name,
        "part": frag.part_label,
        "direction": None,   # filled in below
        "oh": None,
        "core": None,
        "full_primer": None,
        "Tm": None,
        "gc_clamp_ok": None,
        "tm_ok": None,
        "fwd": {
            "core": fwd_core_raw,
            "full_primer": fwd_full,
            "Tm": fwd_tm,
            "gc_clamp_ok": fwd_gc_ok,
            "tm_ok": fwd_tm_ok,
            "oh": oh5,
        },
        "rev": {
            "core": rev_core_raw,
            "full_primer": rev_full,
            "Tm": rev_tm,
            "gc_clamp_ok": rev_gc_ok,
            "tm_ok": rev_tm_ok,
            "oh": revcomp(oh3),
        },
    }


def design_primers_for_parts(parts: List[Part]) -> pd.DataFrame:
    """
    Design all primers for a list of Parts.
    Returns a DataFrame with columns:
        name, part, direction, overhang, core, full_primer, Tm, gc_clamp_ok, tm_ok
    """
    rows = []
    for part in parts:
        n_frags = len(part.fragments)
        for i, frag in enumerate(part.fragments):
            is_a1   = (part.label == "A"      and i == 0)
            is_e_last = (part.label == "E"    and i == n_frags - 1)
            result = design_fragment_primers(frag, is_a1, is_e_last)
            for direction in ("fwd", "rev"):
                d = result[direction]
                rows.append({
                    "name":        frag.name,
                    "part":        part.label,
                    "direction":   direction,
                    "overhang":    d["oh"],
                    "core":        d["core"],
                    "full_primer": d["full_primer"],
                    "Tm":          d["Tm"],
                    "gc_clamp_ok": d["gc_clamp_ok"],
                    "tm_ok":       d["tm_ok"],
                })
    return pd.DataFrame(rows)


def build_pcr_product(frag: Fragment,
                      is_first_of_part_a: bool = False,
                      is_last_of_part_e: bool = False) -> str:
    """
    Return the full double-stranded PCR product sequence (top strand, 5'→3').
    This is: fwd_full_primer + internal_sequence + revcomp(rev_full_primer)
    But since the primers already contain the overhang + gene sequence, the
    product is simply:
      fwd_prefix + oh5 + frag_sequence + revcomp(oh3) + revcomp(rev_prefix)

    Actually more precisely:
      top strand 5'→3' = fwd_full_primer + inner_seq + revcomp(rev_full_primer)
    where inner_seq is the portion of frag.sequence NOT covered by either primer.
    """
    result = design_fragment_primers(frag, is_first_of_part_a, is_last_of_part_e)
    fwd_full = result["fwd"]["full_primer"]
    rev_full = result["rev"]["full_primer"]

    seq = frag.sequence.upper()

    # Find overlap of fwd primer with sequence start and rev primer with sequence end
    fwd_oh   = result["fwd"]["oh"]
    fwd_core = result["fwd"]["core"]
    rev_oh   = result["rev"]["oh"]    # already revcomp'd of oh3
    rev_core = result["rev"]["core"]

    # PCR product top strand:
    #  [fwd prefix] [oh5] [frag sequence (full)] [revcomp(oh3 region) already in rev_core]
    # The full PCR product top strand is:
    #   fwd_primer_full + (frag_seq that isn't already in fwd primer core) + ...
    # Simplest: fwd_primer covers first N nt of product, rev_primer covers last M nt.
    # Middle = frag.sequence[len(fwd_core) : len(frag.sequence)-len(rev_core)]

    # Middle region of frag (between fwd_core end and rev_core start on frag)
    fwd_len_in_frag = len(fwd_core)
    rev_len_in_frag = len(rev_core)
    middle = seq[fwd_len_in_frag : len(seq) - rev_len_in_frag]

    # fwd_full already contains FWD_PREFIX + oh5 + (ATATC if A1) + fwd_core
    # rev_full already contains REV_PREFIX + revcomp(oh3) + (ATATC if E-last) + rev_core
    top_strand = fwd_full + middle + revcomp(rev_full)
    return top_strand
