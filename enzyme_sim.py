"""
enzyme_sim.py – Simulate BsaI and BsmBI digestion on circular or linear DNA.

BsaI  recognition: GGTCTC(1/5)  →  cuts 1 nt 3' of recognition on top strand,
                                     5 nt 3' of recognition on bottom strand
BsmBI recognition: CGTCTC(1/5)  →  same offset pattern

For a forward (top-strand) site at index i of GGTCTC:
  top strand cut  → after i+6+1  = i+7
  bottom strand cut → after i+6+5 = i+11
  → 4-nt 5' overhang: seq[i+7 : i+11]

For a reverse (complement) site the recognition sequence on the bottom strand
reads 5'→3' as the complement, which means on the top strand we see the
reverse complement of GGTCTC, i.e. GAGACC.
  If top strand has GAGACC starting at index j:
  bottom strand recognition runs 3'→5' at positions j..j+5
  bottom strand cut → at j-1 on top strand (1 nt upstream of start)
  top strand cut   → at j-5 on top strand (5 nt upstream of start)
  → 4-nt 5' overhang (on bottom strand) = revcomp(seq[j-5 : j-1])
  equivalently: top strand fragment ends at j-1, bottom strand at j-5
  → sticky 5' overhang on bottom  = revcomp(seq[j-5:j-1])

We return LinearFragment objects; the caller decides how to ligate them.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Tuple, Optional

# Recognition sequences (5'→3' on top strand)
_BSAI_FWD  = "GGTCTC"   # forward BsaI
_BSAI_REV  = "GAGACC"   # reverse complement of BsaI recognition
_BSMBI_FWD = "CGTCTC"   # forward BsmBI
_BSMBI_REV = "GAGACG"   # reverse complement

_CUT_OFFSET_NEAR = 1   # nt past recognition on recognition strand
_CUT_OFFSET_FAR  = 5   # nt past recognition on complementary strand
# → overhang length = FAR - NEAR = 4


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


@dataclass
class CutSite:
    """
    A single restriction site on a (linearised) DNA string.
    position   : index of first base of recognition sequence on top strand
    strand     : +1 = top strand (forward), -1 = bottom strand (reverse)
    enzyme     : 'BsaI' or 'BsmBI'
    top_cut    : position AFTER which the top strand is cut (0-based, exclusive)
    bot_cut    : position AFTER which the bottom strand is cut (0-based, exclusive)
    overhang   : 4-nt 5' sticky end (always reported on the 5' overhang top strand)
    """
    position: int
    strand: int
    enzyme: str
    top_cut: int
    bot_cut: int
    overhang: str


@dataclass
class LinearFragment:
    """
    A linear double-stranded DNA fragment produced by restriction digestion.
    sequence   : top strand 5'→3'
    left_oh    : 5' overhang on the LEFT  end (4 nt if present, '' if blunt)
    right_oh   : 5' overhang on the RIGHT end (4 nt if present, '' if blunt)
    left_oh_strand  : +1 = top strand hangs over, -1 = bottom strand hangs over
    right_oh_strand : same
    """
    sequence: str
    left_oh: str = ""
    right_oh: str = ""
    left_oh_strand: int = 1
    right_oh_strand: int = -1


def _find_sites(seq: str, enzyme: str) -> List[CutSite]:
    """
    Locate all cut sites for `enzyme` in `seq` (top strand).
    seq is treated as linear (caller should handle circularity).
    """
    seq_upper = seq.upper()
    if enzyme == "BsaI":
        fwd_pat, rev_pat = _BSAI_FWD, _BSAI_REV
    elif enzyme == "BsmBI":
        fwd_pat, rev_pat = _BSMBI_FWD, _BSMBI_REV
    else:
        raise ValueError(f"Unknown enzyme: {enzyme}")

    sites: List[CutSite] = []

    # Forward sites
    i = 0
    while True:
        idx = seq_upper.find(fwd_pat, i)
        if idx == -1:
            break
        top_cut = idx + len(fwd_pat) + _CUT_OFFSET_NEAR   # = idx+7
        bot_cut = idx + len(fwd_pat) + _CUT_OFFSET_FAR    # = idx+11
        overhang = seq_upper[top_cut:bot_cut] if bot_cut <= len(seq) else ""
        sites.append(CutSite(
            position=idx, strand=1, enzyme=enzyme,
            top_cut=top_cut, bot_cut=bot_cut, overhang=overhang,
        ))
        i = idx + 1

    # Reverse sites
    i = 0
    while True:
        idx = seq_upper.find(rev_pat, i)
        if idx == -1:
            break
        # On the reverse recognition strand the cut positions are mirrored:
        # top strand is cut  idx - _CUT_OFFSET_NEAR  = idx - 1
        # bot strand is cut  idx - _CUT_OFFSET_FAR   = idx - 5
        top_cut = idx - _CUT_OFFSET_NEAR        # = idx-1, cut BEFORE this position
        bot_cut = idx - _CUT_OFFSET_FAR         # = idx-5
        if bot_cut < 0:
            i = idx + 1
            continue
        overhang = revcomp(seq_upper[bot_cut:top_cut])  # 4-nt bottom-strand 5' OH
        sites.append(CutSite(
            position=idx, strand=-1, enzyme=enzyme,
            top_cut=top_cut, bot_cut=bot_cut, overhang=overhang,
        ))
        i = idx + 1

    return sorted(sites, key=lambda s: s.top_cut)


# ---------------------------------------------------------------------------
# Digest a linear DNA string → list of LinearFragment
# ---------------------------------------------------------------------------

def digest_linear(seq: str, enzyme: str) -> List[LinearFragment]:
    """
    Simulate a complete restriction digest of a LINEAR dna sequence.
    Returns fragments in order (left to right along top strand).
    """
    seq = seq.upper()
    sites = _find_sites(seq, enzyme)

    if not sites:
        # No cut – return entire sequence as one blunt-ended fragment
        return [LinearFragment(sequence=seq)]

    fragments: List[LinearFragment] = []
    prev_top = 0
    prev_oh  = ""
    prev_oh_strand = 1

    for site in sites:
        top_cut = site.top_cut
        bot_cut = site.bot_cut

        if site.strand == 1:
            # Forward site: fragment to the LEFT ends at top_cut
            # Left fragment: seq[prev_top : top_cut], right end has 5' bot overhang
            frag_seq = seq[prev_top:top_cut]
            frag = LinearFragment(
                sequence=frag_seq,
                left_oh=prev_oh,
                left_oh_strand=prev_oh_strand,
                right_oh=site.overhang,
                right_oh_strand=-1,    # bottom strand overhangs on the right
            )
            fragments.append(frag)
            prev_top = bot_cut
            prev_oh  = site.overhang
            prev_oh_strand = 1         # next fragment's left is a top-strand 5' OH

        else:  # strand == -1  (reverse site)
            # Fragment to the LEFT ends at bot_cut (the farther cut)
            frag_seq = seq[prev_top:bot_cut]
            frag = LinearFragment(
                sequence=frag_seq,
                left_oh=prev_oh,
                left_oh_strand=prev_oh_strand,
                right_oh=site.overhang,
                right_oh_strand=1,     # top strand overhangs on the right
            )
            fragments.append(frag)
            prev_top = top_cut
            prev_oh  = site.overhang
            prev_oh_strand = -1

    # Last fragment
    frag_seq = seq[prev_top:]
    fragments.append(LinearFragment(
        sequence=frag_seq,
        left_oh=prev_oh,
        left_oh_strand=prev_oh_strand,
    ))

    return fragments


# ---------------------------------------------------------------------------
# Digest a CIRCULAR DNA string → list of LinearFragment
# ---------------------------------------------------------------------------

def digest_circular(seq: str, enzyme: str) -> List[LinearFragment]:
    """
    Simulate a complete restriction digest of a CIRCULAR dna sequence.
    Internally doubles the sequence to catch sites that span the origin,
    then returns unique fragments.
    """
    seq = seq.upper()
    n = len(seq)
    doubled = seq + seq
    sites = _find_sites(doubled, enzyme)

    # Keep only sites whose recognition sequence starts within [0, n)
    sites = [s for s in sites if s.position < n]

    if not sites:
        return [LinearFragment(sequence=seq)]

    # Sort by top_cut position
    sites = sorted(sites, key=lambda s: s.top_cut)

    fragments: List[LinearFragment] = []
    first_top = sites[0].top_cut if sites[0].strand == 1 else sites[0].bot_cut
    first_oh  = sites[0].overhang
    first_oh_strand = 1 if sites[0].strand == 1 else -1

    prev_top = first_top
    prev_oh  = first_oh
    prev_oh_strand = first_oh_strand

    for site in sites[1:]:
        top_cut = site.top_cut % (2 * n)
        bot_cut = site.bot_cut % (2 * n)

        if site.strand == 1:
            raw = doubled[prev_top:top_cut]
            frag_seq = raw if len(raw) <= n else raw[:n]
            fragments.append(LinearFragment(
                sequence=frag_seq,
                left_oh=prev_oh,
                left_oh_strand=prev_oh_strand,
                right_oh=site.overhang,
                right_oh_strand=-1,
            ))
            prev_top = bot_cut % (2 * n)
            prev_oh  = site.overhang
            prev_oh_strand = 1
        else:
            raw = doubled[prev_top:bot_cut]
            frag_seq = raw if len(raw) <= n else raw[:n]
            fragments.append(LinearFragment(
                sequence=frag_seq,
                left_oh=prev_oh,
                left_oh_strand=prev_oh_strand,
                right_oh=site.overhang,
                right_oh_strand=1,
            ))
            prev_top = top_cut % (2 * n)
            prev_oh  = site.overhang
            prev_oh_strand = -1

    # Last fragment wraps back to the first cut
    raw = doubled[prev_top : first_top + n]  # handles wrap
    end = first_top + n
    raw = doubled[prev_top:end]
    frag_seq = raw[:n]  # cap at genome length
    fragments.append(LinearFragment(
        sequence=frag_seq,
        left_oh=prev_oh,
        left_oh_strand=prev_oh_strand,
        right_oh=first_oh,
        right_oh_strand=first_oh_strand,
    ))

    return fragments


# ---------------------------------------------------------------------------
# Golden Gate ligation: join fragments by matching overhangs
# ---------------------------------------------------------------------------

def golden_gate_ligate(
    fragments: List[LinearFragment],
    circular: bool = False,
) -> str:
    """
    Given an ordered list of LinearFragment objects with compatible overhangs,
    return the joined sequence (top strand only).

    The fragments must be provided in the correct order (the pipeline handles
    ordering by matching overhang sequences).  Only the sequence body is
    concatenated — overhangs are single-stranded in real life but we just
    keep the top-strand sequence of each fragment as-is.
    """
    if not fragments:
        return ""
    seqs = [f.sequence for f in fragments]
    joined = "".join(seqs)
    return joined


def order_fragments_by_overhangs(
    fragments: List[LinearFragment],
    start_oh: str,
) -> List[LinearFragment]:
    """
    Greedily order fragments so that each right_oh matches the next left_oh.
    start_oh: the overhang that the first fragment should have on its left.
    Returns ordered list, or raises ValueError if ordering fails.
    """
    remaining = list(fragments)
    ordered: List[LinearFragment] = []

    current_oh = start_oh.upper()

    for _ in range(len(fragments)):
        found = None
        for f in remaining:
            if f.left_oh.upper() == current_oh:
                found = f
                break
        if found is None:
            raise ValueError(
                f"Cannot find fragment with left overhang '{current_oh}'. "
                f"Available: {[f.left_oh for f in remaining]}"
            )
        ordered.append(found)
        remaining.remove(found)
        current_oh = found.right_oh.upper()

    return ordered
