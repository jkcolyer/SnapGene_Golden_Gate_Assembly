"""
fragment_splitter.py – NEBridge-like gene splitter for Golden Gate assembly.

Given a full CDS (or gene) sequence, splits it into:
  - 5 Parts  (A, B, C, D', E)
  - Each Part into ~200 bp fragments

Inter-Part boundaries are found by searching the gene for the required
4-nt junction overhangs (AATG, AGCC, TTCG, GCTT) so the final assembly is
scarless.  Intra-Part boundaries are also derived from the actual gene
sequence at the split position.

Fragment sequences do NOT include the overhangs at either end.  The overhangs
are stored in `frag.overhang_5` / `frag.overhang_3` and are added back by
`build_level1_seq` so the assembled sequence is exactly the original gene.

Usage
-----
    from fragment_splitter import split_gene
    parts = split_gene(gene_seq)
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional, Tuple
import random

# ---------------------------------------------------------------------------
# pJUMP-specific fixed overhangs
# ---------------------------------------------------------------------------
# Part A first fragment left overhang (exposed by the BsmBI / BsaI prefix)
PART_A_LEFT_OH   = "GGAG"
# Part E last fragment right overhang — matches the Level 2 backbone right side
PART_E_RIGHT_OH  = "CGCT"

# Inter-Part junction overhangs — must match the BsmBI overhangs on the Level 1
# and Level 2 backbones already in the lab.  These sequences must also appear in
# the gene being assembled so the final product is scarless.
PART_JUNCTIONS = {
    "A_B":   "AATG",
    "B_C":   "AGCC",
    "C_D":   "TTCG",
    "D_E":   "GCTT",
}

# Default target fragment length (bp of insert, not counting overhangs)
TARGET_FRAG_LEN = 190

# Minimum / maximum acceptable fragment length
MIN_FRAG_LEN = 120
MAX_FRAG_LEN = 270

# Number of Parts in the assembly
N_PARTS = 5
PART_LABELS = ["A", "B", "C", "Dprime", "E"]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class Fragment:
    name: str            # e.g. "A.1"
    sequence: str        # insert WITHOUT the overhangs at either end
    overhang_5: str      # 4-nt 5' overhang (from gene, or fixed for A.1 / E.last)
    overhang_3: str      # 4-nt 3' overhang (from gene, or Part junction)
    part_label: str
    index: int           # 1-based within Part


@dataclass
class Part:
    label: str           # A, B, C, Dprime, E
    fragments: List[Fragment] = field(default_factory=list)

    @property
    def sequence(self) -> str:
        """Concatenated sequence of all fragments in this Part (without overhangs)."""
        return "".join(f.sequence for f in self.fragments)


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

def _revcomp(seq: str) -> str:
    return seq.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _is_palindrome(seq: str) -> bool:
    return seq.upper() == _revcomp(seq)


def _gc_count(seq: str) -> int:
    return sum(1 for b in seq.upper() if b in "GC")


def _is_valid_overhang(oh: str, existing: set) -> bool:
    oh = oh.upper()
    if len(oh) != 4:
        return False
    if any(b not in "ACGT" for b in oh):
        return False
    if _is_palindrome(oh):
        return False
    gc = _gc_count(oh)
    if gc < 1 or gc > 3:
        return False
    if oh in existing or _revcomp(oh) in existing:
        return False
    return True


# ---------------------------------------------------------------------------
# Inter-Part boundary search
# ---------------------------------------------------------------------------

def _find_best_junction(gene_seq: str, overhang: str, ideal_pos: int) -> Optional[int]:
    """
    Find the occurrence of `overhang` in `gene_seq` closest to `ideal_pos`.
    Returns the start index, or None if not found.
    """
    oh = overhang.upper()
    seq = gene_seq.upper()
    best_pos: Optional[int] = None
    best_dist = float("inf")
    idx = seq.find(oh)
    while idx != -1:
        dist = abs(idx - ideal_pos)
        if dist < best_dist:
            best_dist = dist
            best_pos = idx
        idx = seq.find(oh, idx + 1)
    return best_pos


def _split_into_parts(
    gene_seq: str,
    find_right_end: bool = True,
) -> List[Tuple[int, int]]:
    """
    Split gene_seq into N_PARTS chunks by finding each boundary overhang.

    The 4-nt junction overhang sits BETWEEN the two Parts: it is excluded from
    both Part sequences (stored as oh3/oh5 instead).  This gives scarless
    assembly when build_level1_seq re-inserts the overhangs.

    If `find_right_end` is True the function also searches for PART_E_RIGHT_OH
    (CGCT) near the end of the gene and uses that as the right boundary of
    Part E — any sequence after CGCT is outside the assembly window.
    """
    n = len(gene_seq)
    junction_ohs = [
        PART_JUNCTIONS["A_B"],
        PART_JUNCTIONS["B_C"],
        PART_JUNCTIONS["C_D"],
        PART_JUNCTIONS["D_E"],
    ]
    ideal_positions = [int(n * (i + 1) / N_PARTS) for i in range(N_PARTS - 1)]

    split_positions: List[int] = []
    for oh, ideal in zip(junction_ohs, ideal_positions):
        pos = _find_best_junction(gene_seq, oh, ideal)
        if pos is None:
            pos = ideal  # fallback if overhang not found in gene
        split_positions.append(pos)

    # Find the right end of Part E (CGCT near end of gene)
    part_e_end = n
    if find_right_end:
        cgct_pos = _find_best_junction(gene_seq, PART_E_RIGHT_OH, n)
        if cgct_pos is not None:
            part_e_end = cgct_pos   # Part E sequence ends just before CGCT

    # Build (start, end) for each Part.
    # The 4-nt junction overhang is excluded from both Part sequences.
    boundaries: List[Tuple[int, int]] = []
    prev = 0
    for split_pos in split_positions:
        boundaries.append((prev, split_pos))
        prev = split_pos + 4       # skip past the 4-nt overhang
    boundaries.append((prev, part_e_end))
    return boundaries


# ---------------------------------------------------------------------------
# Intra-Part split search
# ---------------------------------------------------------------------------

def _find_valid_intra_split(
    gene_seq: str,
    ideal_pos: int,
    used_ohs: set,
    window: int = 40,
) -> Tuple[int, str]:
    """
    Find a valid 4-nt intra-Part junction overhang near `ideal_pos`.

    The overhang IS the 4 bases at the chosen split position in the gene, so
    the split position and the overhang are always consistent.

    Returns (split_pos, overhang_4nt) where gene_seq[split_pos:split_pos+4]
    is the overhang, and the left fragment ends at split_pos while the right
    fragment starts at split_pos+4.
    """
    seq = gene_seq.upper()

    # Try positions in increasing distance from ideal, checking full validity
    for offset in range(0, window + 1):
        for delta in ([0] if offset == 0 else [offset, -offset]):
            idx = ideal_pos + delta
            if idx < 0 or idx + 4 > len(seq):
                continue
            candidate = seq[idx:idx + 4]
            if _is_valid_overhang(candidate, used_ohs):
                return idx, candidate

    # Relax: only require uniqueness (no palindrome / GC checks)
    for offset in range(0, window + 1):
        for delta in ([0] if offset == 0 else [offset, -offset]):
            idx = ideal_pos + delta
            if idx < 0 or idx + 4 > len(seq):
                continue
            candidate = seq[idx:idx + 4]
            if (all(b in "ACGT" for b in candidate)
                    and candidate not in used_ohs
                    and _revcomp(candidate) not in used_ohs):
                return idx, candidate

    # Last resort: just use ideal position
    candidate = seq[ideal_pos:ideal_pos + 4]
    return ideal_pos, candidate


# ---------------------------------------------------------------------------
# Intra-Part fragment splitting
# ---------------------------------------------------------------------------

def _split_part_into_fragments(
    part_start: int,
    part_end: int,
    gene_seq: str,
    target_len: int,
    used_ohs: set,
) -> Tuple[List[Tuple[int, int]], List[str]]:
    """
    Split one Part into ~target_len fragments.

    Returns:
      boundaries – list of (abs_start, abs_end) for each fragment.
                   Fragment sequence = gene_seq[start:end].
                   The 4-nt junction overhang at each split is EXCLUDED from
                   both adjacent fragment sequences (stored in oh3/oh5).
      split_ohs  – list of 4-nt overhangs at each intra-Part split, in order.
                   len(split_ohs) == len(boundaries) - 1
    """
    n = part_end - part_start
    if n <= MAX_FRAG_LEN:
        return [(part_start, part_end)], []

    n_frags = max(1, round(n / target_len))
    split_positions: List[int] = []
    split_ohs: List[str] = []

    for i in range(1, n_frags):
        ideal = part_start + int(n * i / n_frags)
        pos, oh = _find_valid_intra_split(gene_seq, ideal, used_ohs)
        # Register the overhang so it is not reused
        used_ohs.add(oh)
        used_ohs.add(_revcomp(oh))
        split_positions.append(pos)
        split_ohs.append(oh)

    # Build fragment boundaries with 4-nt gaps at each split
    boundaries: List[Tuple[int, int]] = []
    prev = part_start
    for sp in split_positions:
        boundaries.append((prev, sp))
        prev = sp + 4              # skip the 4-nt junction overhang
    boundaries.append((prev, part_end))

    return boundaries, split_ohs


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def split_gene(
    gene_seq: str,
    target_frag_len: int = TARGET_FRAG_LEN,
) -> List[Part]:
    """
    Split a full gene sequence into 5 Parts, each containing multiple ~200 bp
    fragments with unique 4-nt Golden Gate overhangs.

    Overhang contract
    -----------------
    * frag.sequence does NOT include the overhangs at either end.
    * oh3 of fragment[i] == oh5 of fragment[i+1] (guaranteed).
    * The junction overhang = gene_seq[split_pos : split_pos+4].
    * build_level1_seq re-inserts oh5 before each fragment so the assembled
      Level 1 sequence is the original gene (scarless).

    Parameters
    ----------
    gene_seq        : full gene/CDS sequence (A/C/G/T)
    target_frag_len : target fragment insert length in bp (not counting overhangs)

    Returns
    -------
    List of 5 Part objects in order A, B, C, Dprime, E.
    """
    gene_seq = gene_seq.upper().replace(" ", "").replace("\n", "")
    used_ohs: set = set()

    # Reserve all fixed/junction overhangs
    for oh in ([PART_A_LEFT_OH, PART_E_RIGHT_OH]
               + list(PART_JUNCTIONS.values())):
        used_ohs.add(oh.upper())
        used_ohs.add(_revcomp(oh.upper()))

    # Step 1: Find inter-Part boundaries (and CGCT right-end) in the gene
    part_boundaries = _split_into_parts(gene_seq, find_right_end=True)

    parts: List[Part] = []

    for part_idx, (part_start, part_end) in enumerate(part_boundaries):
        label = PART_LABELS[part_idx]

        # Step 2: Find intra-Part split positions
        frag_boundaries, split_ohs = _split_part_into_fragments(
            part_start, part_end, gene_seq, target_frag_len, used_ohs
        )

        fragments: List[Fragment] = []
        n_frags = len(frag_boundaries)

        for frag_idx, (f_start, f_end) in enumerate(frag_boundaries):
            frag_seq = gene_seq[f_start:f_end]
            frag_num = frag_idx + 1

            # ---- 5' overhang ----
            if part_idx == 0 and frag_idx == 0:
                # Part A frag 1: fixed left overhang from primer prefix
                oh5 = PART_A_LEFT_OH
            elif frag_idx == 0:
                # First frag of Parts B/C/D'/E: inter-Part junction overhang
                jk = f"{PART_LABELS[part_idx-1]}_{label}".replace("Dprime", "D")
                oh5 = PART_JUNCTIONS[jk]
            else:
                # Internal frag: oh5 is the split overhang chosen just before
                # this fragment (index frag_idx-1 in split_ohs)
                oh5 = split_ohs[frag_idx - 1]

            # ---- 3' overhang ----
            if part_idx == N_PARTS - 1 and frag_idx == n_frags - 1:
                # Part E last frag: fixed right overhang
                oh3 = PART_E_RIGHT_OH
            elif frag_idx == n_frags - 1:
                # Last frag of Part: inter-Part junction overhang
                jk = f"{label}_{PART_LABELS[part_idx+1]}".replace("Dprime", "D")
                oh3 = PART_JUNCTIONS[jk]
            else:
                # Internal frag: oh3 is the next split overhang
                oh3 = split_ohs[frag_idx]

            fragments.append(Fragment(
                name=f"splitseq_{frag_num}",
                sequence=frag_seq,
                overhang_5=oh5,
                overhang_3=oh3,
                part_label=label,
                index=frag_num,
            ))

        parts.append(Part(label=label, fragments=fragments))

    return parts


def parts_to_fasta(parts: List[Part]) -> str:
    """Return a FASTA-like string for all fragments across all parts."""
    lines = []
    for part in parts:
        for frag in part.fragments:
            lines.append(f">{frag.part_label}.{frag.index} {len(frag.sequence)}bp")
            lines.append(frag.sequence)
    return "\n".join(lines) + "\n"
