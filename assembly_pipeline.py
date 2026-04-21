"""
assembly_pipeline.py – Orchestrate the full Golden Gate assembly pipeline.

Pipeline:
  1. Split gene into 5 Parts × N fragments (via fragment_splitter)
  2. Design primers for every fragment (via primer_designer)
  3. Build PCR products (.dna, linear)
  4. Build Level 0 constructs (.dna, circular) = PCR product in backbone
  5. Build Level 1 constructs (.dna, circular) = all Level 0 inserts for a Part
  6. Build Level 2 construct (.dna, circular)  = all 5 Parts assembled

Entry point: run_pipeline(gene_seq, output_dir)
Returns: pandas DataFrame with the primer table.
All .dna files are written into output_dir.
"""

from __future__ import annotations
import os
import zipfile
from pathlib import Path
from typing import List, Dict, Optional

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font

from fragment_splitter import split_gene, Part, Fragment
from primer_designer import (
    design_primers_for_parts,
    design_fragment_primers,
    build_pcr_product,
    FWD_PREFIX, REV_PREFIX, FWD_PREFIX_A1, REV_PREFIX_E_LAST,
)
from dna_writer import (
    write_dna_file, make_feature, make_primer_feature,
    FLAG_LINEAR_DS, FLAG_CIRCULAR,
)
from backbone_constants import (
    BACKBONE_LEFT_A, BACKBONE_RIGHT_A,
    BACKBONE_LEFT_BCD, BACKBONE_RIGHT_BCD,
    LV1_LEFT_A, LV1_RIGHT_A,
    LV1_LEFT_BCD, LV1_RIGHT_B, LV1_RIGHT_C,
    LEVEL2_PART_ORDER, LEVEL2_LEFT_OHS, LEVEL2_RIGHT_OHS,
)

# Part E uses the same flanks as Part A
LV1_LEFT_E  = LV1_LEFT_A
LV1_RIGHT_E = LV1_RIGHT_A

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# PCR product: the insert begins 13 nt into the product (right after BsmBI cuts
# following the fwd prefix).  The revcomp of the rev prefix (13 nt) follows.
PCR_INSERT_START = 13

# Part label normalisation
_DPRIME = "Dprime"
_PART_LABELS = ["A", "B", "C", _DPRIME, "E"]

# Level 0 backbone sequences per Part family
_LV0_BACKBONE_LEFT = {
    "A": BACKBONE_LEFT_A,
    "B": BACKBONE_LEFT_BCD,
    "C": BACKBONE_LEFT_BCD,
    _DPRIME: BACKBONE_LEFT_BCD,
    "E": BACKBONE_LEFT_BCD,
}
_LV0_BACKBONE_RIGHT = {
    "A": BACKBONE_RIGHT_A,
    "B": BACKBONE_RIGHT_BCD,
    "C": BACKBONE_RIGHT_BCD,
    _DPRIME: BACKBONE_RIGHT_BCD,
    "E": BACKBONE_RIGHT_BCD,
}

# Level 1 left/right flanks per Part
_LV1_LEFT_FLANK = {
    "A": LV1_LEFT_A,
    "B": LV1_LEFT_BCD,
    "C": LV1_LEFT_BCD,
    _DPRIME: LV1_LEFT_BCD,
    "E": LV1_LEFT_E,
}
_LV1_RIGHT_FLANK = {
    "A": LV1_RIGHT_A,
    "B": LV1_RIGHT_B,
    "C": LV1_RIGHT_C,
    _DPRIME: LV1_RIGHT_C,   # same as C
    "E": LV1_RIGHT_E,
}


def revcomp(seq: str) -> str:
    return seq.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def mask_enzyme_sites(seq: str) -> str:
    """
    Silently mutate any BsaI or BsmBI recognition sites in the gene insert
    sequence so they cannot be cut by those enzymes during assembly.

    Each substitution is a single-base change that destroys recognition:
      BsaI  fwd  GGTCTC → GATCTC   (G→A at position 1)
      BsaI  RC   GAGACC → GAGATC   (C→T at position 5, RC of above)
      BsmBI fwd  CGTCTC → CGTATC   (C→A at position 3)
      BsmBI RC   GAGACG → GATACG   (G→T at position 2, RC of above)

    Only applies to the gene-binding / insert sequence — NOT to primer tails or
    backbone sequences, which intentionally carry these sites.
    """
    seq = seq.upper()
    seq = seq.replace("GGTCTC", "GATCTC")   # BsaI fwd
    seq = seq.replace("GAGACC", "GAGATC")   # BsaI rev-complement
    seq = seq.replace("CGTCTC", "CGTATC")   # BsmBI fwd
    seq = seq.replace("GAGACG", "GATACG")   # BsmBI rev-complement  ← one-base change
    return seq


# ---------------------------------------------------------------------------
# Step 3 – Build PCR products
# ---------------------------------------------------------------------------

def build_pcr_dna(
    frag: Fragment,
    is_first_of_part_a: bool = False,
    is_last_of_part_e: bool = False,
) -> str:
    """Return the full PCR product sequence (top strand 5'→3')."""
    return build_pcr_product(frag, is_first_of_part_a, is_last_of_part_e)


def write_pcr_file(
    frag: Fragment,
    path: str,
    is_first_of_part_a: bool = False,
    is_last_of_part_e: bool = False,
) -> None:
    """Write a PCR product .dna file."""
    pcr_seq = build_pcr_dna(frag, is_first_of_part_a, is_last_of_part_e)
    primers_result = design_fragment_primers(frag, is_first_of_part_a, is_last_of_part_e)

    fwd_full = primers_result["fwd"]["full_primer"]
    rev_full = primers_result["rev"]["full_primer"]
    fwd_tm   = primers_result["fwd"]["Tm"]
    rev_tm   = primers_result["rev"]["Tm"]

    # FWD primer covers positions [0 : len(fwd_full)]
    # REV primer covers positions [len(pcr_seq)-len(rev_full) : len(pcr_seq)]
    # Gene insert (just the frag sequence, no tails) is in between
    fwd_len = len(fwd_full)
    rev_start = len(pcr_seq) - len(rev_full)

    features = [
        make_feature(
            name=f"{frag.part_label}.{frag.index} FWD primer",
            start=0,
            end=fwd_len,
            strand=1,
            ftype="primer_bind",
            color="#cc0000",
            note=f"Tm={fwd_tm:.1f}°C",
        ),
        make_feature(
            name=f"gene insert",
            start=fwd_len,
            end=rev_start,
            strand=1,
            ftype="CDS",
            color="#84b0dc",
        ),
        make_feature(
            name=f"{frag.part_label}.{frag.index} REV primer",
            start=rev_start,
            end=len(pcr_seq),
            strand=-1,
            ftype="primer_bind",
            color="#0044cc",
            note=f"Tm={rev_tm:.1f}°C",
        ),
    ]

    write_dna_file(
        path=path,
        sequence=pcr_seq,
        is_circular=False,
        features=features,
        name=f"Part {frag.part_label}.{frag.index} PCR",
    )


# ---------------------------------------------------------------------------
# Step 4 – Build Level 0 constructs
# ---------------------------------------------------------------------------

def build_level0_seq(
    frag: Fragment,
    is_first_of_part_a: bool = False,
    is_last_of_part_e: bool = False,
) -> str:
    """
    Construct the Level 0 circular sequence.
    Level 0 = backbone_left + PCR_product[19:] + backbone_right
    (circularized)
    """
    pcr = build_pcr_dna(frag, is_first_of_part_a, is_last_of_part_e)
    insert = pcr[PCR_INSERT_START:]   # from BsmBI cut site onward (oh5 + insert + oh3-rc + prefix-rc)

    part = frag.part_label
    bb_left  = _LV0_BACKBONE_LEFT[part]
    bb_right = _LV0_BACKBONE_RIGHT[part]

    return bb_left + insert + bb_right


def write_level0_file(
    frag: Fragment,
    path: str,
    is_first_of_part_a: bool = False,
    is_last_of_part_e: bool = False,
) -> None:
    """Write a Level 0 .dna file."""
    seq   = build_level0_seq(frag, is_first_of_part_a, is_last_of_part_e)
    pcr   = build_pcr_dna(frag, is_first_of_part_a, is_last_of_part_e)
    insert = pcr[PCR_INSERT_START:]   # oh5 + [ATATC] + frag.seq + revcomp(oh3) + tail

    bb_left_len  = len(_LV0_BACKBONE_LEFT[frag.part_label])
    bb_right_start = bb_left_len + len(insert)

    primers_result = design_fragment_primers(frag, is_first_of_part_a, is_last_of_part_e)
    fwd_full     = primers_result["fwd"]["full_primer"]
    rev_full     = primers_result["rev"]["full_primer"]
    fwd_full_len = len(fwd_full)
    rev_full_len = len(rev_full)
    fwd_tm       = primers_result["fwd"]["Tm"]
    rev_tm       = primers_result["rev"]["Tm"]

    # In Level 0 the PCR product is embedded starting at bb_left_len - PCR_INSERT_START.
    # The gene sequence (frag.sequence) starts after the oh5 (and ATATC for A1).
    oh5_len   = len(frag.overhang_5)
    ecorv_len = 5 if is_first_of_part_a else 0
    gene_start_in_lv0 = bb_left_len + oh5_len + ecorv_len
    gene_end_in_lv0   = gene_start_in_lv0 + len(frag.sequence)

    features = [
        # Backbone left
        make_feature(
            name="backbone left",
            start=0,
            end=bb_left_len,
            strand=1,
            ftype="misc_feature",
            color="#a6acb3",
        ),
        # 5' overhang region (BsaI-exposed)
        make_feature(
            name=f"oh5 {frag.overhang_5}",
            start=bb_left_len,
            end=bb_left_len + oh5_len,
            strand=1,
            ftype="misc_feature",
            color="#f5a623",
        ),
        # Gene insert
        make_feature(
            name=f"{frag.part_label}.{frag.index} gene insert",
            start=gene_start_in_lv0,
            end=gene_end_in_lv0,
            strand=1,
            ftype="CDS",
            color="#84b0dc",
        ),
        # 3' overhang region
        make_feature(
            name=f"oh3 {frag.overhang_3}",
            start=gene_end_in_lv0,
            end=gene_end_in_lv0 + len(frag.overhang_3),
            strand=1,
            ftype="misc_feature",
            color="#f5a623",
        ),
        # Backbone right
        make_feature(
            name="backbone right",
            start=bb_right_start,
            end=len(seq),
            strand=1,
            ftype="misc_feature",
            color="#a6acb3",
        ),
        # FWD primer binding site (spans backbone-left tail + oh5 + gene-core)
        make_feature(
            name=f"{frag.part_label}.{frag.index} FWD primer",
            start=bb_left_len - PCR_INSERT_START,
            end=bb_left_len - PCR_INSERT_START + fwd_full_len,
            strand=1,
            ftype="primer_bind",
            color="#cc0000",
            note=f"Tm={fwd_tm:.1f}°C  seq={fwd_full}",
        ),
        # REV primer binding site
        make_feature(
            name=f"{frag.part_label}.{frag.index} REV primer",
            start=gene_end_in_lv0,
            end=gene_end_in_lv0 + rev_full_len,
            strand=-1,
            ftype="primer_bind",
            color="#0044cc",
            note=f"Tm={rev_tm:.1f}°C  seq={rev_full}",
        ),
    ]

    write_dna_file(
        path=path,
        sequence=seq,
        is_circular=True,
        features=features,
        name=f"Level 0 Part {frag.part_label}.{frag.index}",
    )


# ---------------------------------------------------------------------------
# Step 5 – Build Level 1 constructs
# ---------------------------------------------------------------------------

def build_level1_seq(part: Part) -> str:
    """
    Build the Level 1 circular sequence for a Part.
    Level 1 = left_flank + assembled_gene_insert + right_flank

    Fragment sequences do NOT include overhangs at either end (they were
    excluded during splitting so each junction overhang is stored only in
    frag.overhang_5 / frag.overhang_3).  We re-insert them here so the
    assembled Level 1 contains the original, unmodified gene sequence.

    Assembly rule:
      gene_insert = oh5(frag0) [+ ATATC if Part A frag0] + seq(frag0)
                  + oh5(frag1) + seq(frag1)
                  + ...
                  + oh5(fragN) + seq(fragN)
                  + oh3(fragN)          ← Part-boundary or CGCT overhang

    This gives:
      Part A: GGAG + ATATC + gene[0..aatg) + AATG
      Part B: AATG + gene[aatg+4..agcc) + AGCC
      ...
      Part E: GCTT + gene[gctt+4..end) + CGCT
    """
    label = part.label
    left_flank  = _LV1_LEFT_FLANK[label]
    right_flank = _LV1_RIGHT_FLANK[label]

    gene_parts = []
    for i, frag in enumerate(part.fragments):
        is_a1 = (label == "A" and i == 0)
        # Always prepend the 5' overhang for this fragment
        gene_parts.append(frag.overhang_5)
        if is_a1:
            # EcoRV context added between GGAG overhang and gene start
            gene_parts.append("ATATC")
        gene_parts.append(frag.sequence)

    # Append the 3' overhang of the last fragment (the Part-boundary overhang)
    if part.fragments:
        gene_parts.append(part.fragments[-1].overhang_3)

    gene_insert = "".join(gene_parts)
    return left_flank + gene_insert + right_flank


def write_level1_file(part: Part, path: str) -> None:
    """Write a Level 1 .dna file."""
    seq = build_level1_seq(part)
    left_flank  = _LV1_LEFT_FLANK[part.label]
    right_flank = _LV1_RIGHT_FLANK[part.label]
    gene_region_start = len(left_flank)
    gene_region_end   = len(seq) - len(right_flank)

    features = [
        # Backbone left
        make_feature(
            name="backbone",
            start=0,
            end=gene_region_start,
            strand=1,
            ftype="misc_feature",
            color="#a6acb3",
        ),
        # Whole Part gene insert
        make_feature(
            name=f"Part {part.label} insert",
            start=gene_region_start,
            end=gene_region_end,
            strand=1,
            ftype="misc_feature",
            color="#c3e5c3",
        ),
        # Backbone right
        make_feature(
            name="backbone",
            start=gene_region_end,
            end=len(seq),
            strand=1,
            ftype="misc_feature",
            color="#a6acb3",
        ),
    ]

    # Walk through the assembled gene insert to annotate each fragment
    # (mirrors the logic in build_level1_seq)
    pos = gene_region_start
    for i, frag in enumerate(part.fragments):
        is_a1 = (part.label == "A" and i == 0)
        oh5_len = len(frag.overhang_5)

        # Overhang feature
        features.append(make_feature(
            name=f"oh5 {frag.overhang_5}",
            start=pos,
            end=pos + oh5_len,
            strand=1,
            ftype="misc_feature",
            color="#f5a623",
        ))
        pos += oh5_len

        if is_a1:
            # EcoRV site (ATATC)
            features.append(make_feature(
                name="EcoRV",
                start=pos,
                end=pos + 5,
                strand=1,
                ftype="misc_feature",
                color="#9b59b6",
            ))
            pos += 5

        # Fragment sequence
        frag_end = pos + len(frag.sequence)
        features.append(make_feature(
            name=f"{part.label}.{frag.index}",
            start=pos,
            end=frag_end,
            strand=1,
            ftype="CDS",
            color="#84b0dc",
        ))
        pos = frag_end

    # Last fragment oh3 (Part-boundary overhang)
    if part.fragments:
        oh3 = part.fragments[-1].overhang_3
        features.append(make_feature(
            name=f"oh3 {oh3}",
            start=pos,
            end=pos + len(oh3),
            strand=1,
            ftype="misc_feature",
            color="#f5a623",
        ))

    write_dna_file(
        path=path,
        sequence=seq,
        is_circular=True,
        features=features,
        name=f"Level 1 Part {part.label}",
    )


# ---------------------------------------------------------------------------
# Step 6 – Build Level 2 construct
# ---------------------------------------------------------------------------

def build_level2_insert(part: Part) -> str:
    """
    Return the gene insert for a Part as exposed by BsmBI digestion of Level 1.
    = Level1[len(left_flank) : len(Level1)-len(right_flank)]
    = oh5 + [ATATC] + frag_seqs + oh3
    """
    lv1_seq     = build_level1_seq(part)
    left_flank  = _LV1_LEFT_FLANK[part.label]
    right_flank = _LV1_RIGHT_FLANK[part.label]
    gene_start  = len(left_flank)
    gene_end    = len(lv1_seq) - len(right_flank)
    return lv1_seq[gene_start:gene_end]


def build_level2_seq(parts: List[Part]) -> str:
    """
    Build the Level 2 sequence by ligating all 5 Part inserts in order A→E.

    In Golden Gate ligation each 4-nt junction overhang appears exactly once
    in the final sequence.  Each Part insert is:
        oh5 + [ATATC if Part A] + gene_region + oh3
    Adjacent Parts share the same 4-nt junction (oh3 of left = oh5 of right),
    so when concatenating digitally we strip the leading oh5 from Parts B–E to
    avoid duplication.

    Final insert:
        Part A oh5(GGAG) + ATATC + gene[0..aatg) + AATG
      + gene[aatg+4..agcc) + AGCC
      + gene[agcc+4..ttcg) + TTCG
      + gene[ttcg+4..gctt) + GCTT
      + gene[gctt+4..cgct) + CGCT
      = GGAG + ATATC + gene[0..cgct+4)
    """
    label_order  = {lab: i for i, lab in enumerate(LEVEL2_PART_ORDER)}
    parts_sorted = sorted(parts, key=lambda p: label_order.get(p.label, 99))

    lv2_left  = _LV1_LEFT_FLANK["A"]
    lv2_right = _LV1_RIGHT_FLANK["E"]

    combined = ""
    for i, part in enumerate(parts_sorted):
        insert = build_level2_insert(part)
        if i > 0:
            # Strip the leading 4-nt oh5: it was already appended as oh3
            # of the previous Part
            insert = insert[4:]
        combined += insert

    return lv2_left + combined + lv2_right


def write_level2_file(parts: List[Part], path: str) -> None:
    """Write the Level 2 .dna file."""
    seq = build_level2_seq(parts)

    label_order = {lab: i for i, lab in enumerate(LEVEL2_PART_ORDER)}
    parts_sorted = sorted(parts, key=lambda p: label_order.get(p.label, 99))

    features = []
    colors = {"A": "#ff9999", "B": "#99ccff", "C": "#99ff99",
              "Dprime": "#ffcc99", "E": "#cc99ff"}

    left_flank_len = len(_LV1_LEFT_FLANK["A"])
    pos = left_flank_len

    for i, part in enumerate(parts_sorted):
        insert = build_level2_insert(part)
        if i > 0:
            insert = insert[4:]   # strip duplicate junction oh5 (same as previous oh3)
        end = pos + len(insert)
        features.append(make_feature(
            name=f"Part {part.label}",
            start=pos,
            end=end,
            strand=1,
            ftype="CDS",
            color=colors.get(part.label, "#a6acb3"),
        ))
        pos = end

    write_dna_file(
        path=path,
        sequence=seq,
        is_circular=True,
        features=features,
        name="Level 2 Full Gene",
    )


# ---------------------------------------------------------------------------
# Ordering Excel files
# ---------------------------------------------------------------------------

def write_ordering_excel(part: Part, gene_name: str, path: str) -> None:
    """
    Write a primer-ordering Excel file for one Part, matching the format used
    in the Athey Lab CYP2D7 ordering sheets exactly:

    Columns (no header row):
      A – Primer name          e.g. splitseq_1_fwd
      B – Gene-specific seq    oh + gene_core (no adapter prefix)
      C – Same as B
      D – Full ordering seq    FWD_PREFIX + B  or  REV_PREFIX + B
      E – Scale                25nm
      F – Purification         STD

    One blank row between each fragment pair.
    Font: Aptos Narrow 12pt (falls back to Calibri on non-Windows).
    Column widths: A=19.83, B=43.33, others default.
    Row height: 16pt.
    """
    wb = Workbook()
    ws = wb.active
    ws.title = "Sheet1"

    # Column widths to match original
    ws.column_dimensions["A"].width = 19.83
    ws.column_dimensions["B"].width = 43.33

    # Default row height
    ws.sheet_format.defaultRowHeight = 16
    ws.sheet_format.customHeight = True

    font = Font(name="Aptos Narrow", size=12)

    n_frags = len(part.fragments)
    current_row = 1

    for i, frag in enumerate(part.fragments):
        is_a1     = (part.label == "A" and i == 0)
        is_e_last = (part.label == "E" and i == n_frags - 1)

        primers = design_fragment_primers(frag, is_a1, is_e_last)
        fwd = primers["fwd"]
        rev = primers["rev"]

        # gene-specific sequences (oh + core, no adapter prefix)
        fwd_gene_seq = fwd["oh"] + fwd["core"]
        rev_gene_seq = rev["oh"] + rev["core"]

        # Full ordering sequences (adapter + gene-specific)
        fwd_full = fwd["full_primer"]
        rev_full = rev["full_primer"]

        frag_label = f"splitseq_{frag.index}"

        for direction, gene_seq, full_seq in [
            ("fwd", fwd_gene_seq, fwd_full),
            ("rev", rev_gene_seq, rev_full),
        ]:
            row_data = [
                f"{frag_label}_{direction}",
                gene_seq.lower(),
                gene_seq.lower(),
                full_seq.lower(),
                "25nm",
                "STD",
            ]
            for col_idx, value in enumerate(row_data, start=1):
                cell = ws.cell(row=current_row, column=col_idx, value=value)
                cell.font = Font(name="Aptos Narrow", size=12)
            current_row += 1

        # Blank row between fragment pairs (but not after the last one)
        if i < n_frags - 1:
            current_row += 1

    wb.save(path)


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def run_pipeline(
    gene_seq: str,
    output_dir: str,
    gene_name: str = "gene",
    target_frag_len: int = 190,
) -> pd.DataFrame:
    """
    Run the full Golden Gate assembly pipeline.

    Parameters
    ----------
    gene_seq    : Full gene/CDS sequence (A/C/G/T, any case)
    output_dir  : Directory to write all .dna files and the primer CSV
    gene_name   : Name used in file names and SnapGene labels
    target_frag_len : Target fragment insert length in bp (default 190)

    Returns
    -------
    pandas DataFrame with the primer table (all fragments, all Parts)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- Step 1: Split gene ----
    parts = split_gene(gene_seq, target_frag_len)

    # ---- Step 1b: Mask internal BsaI / BsmBI sites in insert sequences ----
    for part in parts:
        for frag in part.fragments:
            frag.sequence = mask_enzyme_sites(frag.sequence)

    # ---- Step 2: Design primers ----
    primer_df = design_primers_for_parts(parts)
    primer_df.to_csv(output_dir / f"{gene_name}_primers.csv", index=False)

    # ---- Steps 3 & 4: PCR products + Level 0 ----
    for part in parts:
        n_frags = len(part.fragments)
        part_dir = output_dir / f"Part_{part.label}"
        part_dir.mkdir(exist_ok=True)

        for i, frag in enumerate(part.fragments):
            is_a1    = (part.label == "A" and i == 0)
            is_e_last = (part.label == "E" and i == n_frags - 1)

            frag_dir = part_dir / f"{part.label}.{frag.index}"
            frag_dir.mkdir(exist_ok=True)

            # PCR product
            pcr_path = frag_dir / f"{gene_name}_Part_{part.label}.{frag.index}_PCR.dna"
            write_pcr_file(frag, str(pcr_path), is_a1, is_e_last)

            # Level 0
            lv0_path = frag_dir / f"{gene_name}_Part_{part.label}.{frag.index}.dna"
            write_level0_file(frag, str(lv0_path), is_a1, is_e_last)

    # ---- Step 5: Level 1 constructs + ordering Excel ----
    for part in parts:
        part_dir = output_dir / f"Part_{part.label}"
        lv1_path = part_dir / f"{gene_name}_Part_{part.label}_Level1.dna"
        write_level1_file(part, str(lv1_path))

        xlsx_path = part_dir / f"{gene_name}_Part_{part.label}.xlsx"
        write_ordering_excel(part, gene_name, str(xlsx_path))

    # ---- Step 6: Level 2 construct ----
    lv2_path = output_dir / f"{gene_name}_Level2.dna"
    write_level2_file(parts, str(lv2_path))

    return primer_df


def build_zip(output_dir: str, zip_path: str) -> None:
    """Package all .dna and .csv files into a ZIP for download."""
    output_dir = Path(output_dir)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for file in sorted(output_dir.rglob("*")):
            if file.is_file():
                zf.write(file, file.relative_to(output_dir))
