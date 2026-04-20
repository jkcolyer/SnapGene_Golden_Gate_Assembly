"""
dna_writer.py – Write SnapGene-compatible .dna binary files.

SnapGene .dna format: a sequence of packets.
Each packet: [1-byte type] [4-byte big-endian length] [payload bytes]

Packet types used here:
  0x09 – Cookie / file header (must be first)
  0x00 – DNA sequence + flags
  0x0A – Primers (XML)
  0x0B – Notes (XML)
  0x08 – Features (XML)
  0x07 – Additional properties (XML)
"""

import struct
import datetime
import uuid as _uuid_mod
from typing import Optional, List, Dict

# ---------------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------------
FLAG_LINEAR_DS = 0x02   # linear double-stranded DNA
FLAG_CIRCULAR  = 0x1f   # circular double-stranded DNA (observed in plasmids)

# ---------------------------------------------------------------------------
# Packet assembly helpers
# ---------------------------------------------------------------------------

def _pack_packet(ptype: int, payload: bytes) -> bytes:
    return bytes([ptype]) + struct.pack(">I", len(payload)) + payload


def _cookie_packet() -> bytes:
    """
    Type 0x09 – must be the very first packet.
    Observed payload: b'SnapGene\x00\x01\x00\x0f\x00\x14'
    """
    payload = b"SnapGene\x00\x01\x00\x0f\x00\x14"
    return _pack_packet(0x09, payload)


def _dna_packet(sequence: str, is_circular: bool) -> bytes:
    """
    Type 0x00 – flags byte followed by ASCII sequence bytes.
    """
    flag = FLAG_CIRCULAR if is_circular else FLAG_LINEAR_DS
    payload = bytes([flag]) + sequence.upper().encode("ascii")
    return _pack_packet(0x00, payload)


def _notes_packet(name: str, created_by: str = "Athey lab") -> bytes:
    """
    Type 0x0B – XML notes.  SnapGene shows 'Created with…' here.
    """
    now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    xml = (
        f'<Notes>'
        f'<UUID>{_uuid_mod.uuid4()}</UUID>'
        f'<Type>Natural</Type>'
        f'<CustomMapLabel>{_xml_escape(name)}</CustomMapLabel>'
        f'<UseCustomMapLabel>1</UseCustomMapLabel>'
        f'<CreatedBy>{_xml_escape(created_by)}</CreatedBy>'
        f'<LastModified>{now}</LastModified>'
        f'</Notes>'
    )
    return _pack_packet(0x0B, xml.encode("utf-8"))


def _additional_props_packet(extra_xml: str = "") -> bytes:
    """
    Type 0x07 – AdditionalSequenceProperties.
    For ligation products use phosphorylated ends.
    """
    inner = extra_xml or (
        "<UpstreamModification>FivePrimePhosphorylated</UpstreamModification>"
        "<DownstreamModification>FivePrimePhosphorylated</DownstreamModification>"
    )
    xml = f"<AdditionalSequenceProperties>{inner}</AdditionalSequenceProperties>"
    return _pack_packet(0x07, xml.encode("utf-8"))


# ---------------------------------------------------------------------------
# Feature / annotation helpers
# ---------------------------------------------------------------------------

def _xml_escape(s: str) -> str:
    return (s.replace("&", "&amp;")
             .replace("<", "&lt;")
             .replace(">", "&gt;")
             .replace('"', "&quot;"))


def make_feature(
    name: str,
    start: int,      # 0-based, inclusive
    end: int,        # 0-based, exclusive  →  SnapGene uses 0-based inclusive ranges
    strand: int,     # 1 = forward, -1 = reverse
    ftype: str = "misc_feature",
    color: str = "#a6acb3",
    note: str = "",
) -> Dict:
    """Return a dict representing one feature annotation."""
    return {
        "name": name,
        "start": start,
        "end": end - 1,    # convert to 0-based inclusive for SnapGene
        "strand": strand,
        "type": ftype,
        "color": color,
        "note": note,
    }


def make_primer_feature(
    name: str,
    start: int,      # 0-based inclusive in sequence
    length: int,
    strand: int,     # 1 fwd, -1 rev
    tm: float = 0.0,
    color: str = "#cc0000",
) -> Dict:
    return {
        "name": name,
        "start": start,
        "end": start + length - 1,
        "strand": strand,
        "type": "primer_bind",
        "color": color,
        "tm": tm,
    }


def _features_packet(features: List[Dict]) -> bytes:
    """
    Type 0x08 – Features XML.
    Each feature becomes a <Feature> element.
    """
    if not features:
        return b""
    parts = ['<Features>']
    for f in features:
        strand_val = "1" if f.get("strand", 1) >= 0 else "-1"
        parts.append(
            f'<Feature recentID="0" name="{_xml_escape(f["name"])}" '
            f'directionality="{strand_val}" '
            f'type="{_xml_escape(f.get("type", "misc_feature"))}" '
            f'swappedSegmentNumbering="0" '
            f'allowSegmentOverlaps="0" '
            f'consecutiveTranslationNumbering="0">'
        )
        color = f.get("color", "#a6acb3")
        parts.append(
            f'<Segment range="{f["start"]+1}-{f["end"]+1}" '
            f'color="{color}" type="standard"/>'
        )
        if f.get("note"):
            parts.append(
                f'<Q name="note"><V predef="0" text="{_xml_escape(f["note"])}"/></Q>'
            )
        if f.get("tm"):
            parts.append(
                f'<Q name="Tm"><V predef="0" text="{f["tm"]:.1f}"/></Q>'
            )
        parts.append('</Feature>')
    parts.append('</Features>')
    xml = "".join(parts)
    return _pack_packet(0x08, xml.encode("utf-8"))


def _primers_packet(primers: List[Dict]) -> bytes:
    """
    Type 0x0A – Primers XML.
    primers: list of dicts with keys: name, sequence, strand, start (0-based inclusive), tm
    """
    if not primers:
        return b""
    parts = ['<Primers>']
    for p in primers:
        strand_val = "1" if p.get("strand", 1) >= 0 else "-1"
        tm_str = f'{p.get("tm", 0.0):.1f}'
        parts.append(
            f'<Primer recentID="0" name="{_xml_escape(p["name"])}" '
            f'sequence="{_xml_escape(p["sequence"].upper())}" '
            f'tm="{tm_str}" '
            f'simplified="0">'
        )
        start_1 = p["start"] + 1  # 1-based
        end_1   = p["start"] + len(p["sequence"])
        parts.append(
            f'<BindingSite location="{start_1}-{end_1}" '
            f'boundStrand="{strand_val}" '
            f'meltingTemperature="{tm_str}"/>'
        )
        parts.append('</Primer>')
    parts.append('</Primers>')
    xml = "".join(parts)
    return _pack_packet(0x0A, xml.encode("utf-8"))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def write_dna_file(
    path: str,
    sequence: str,
    is_circular: bool = False,
    features: Optional[List[Dict]] = None,
    primers: Optional[List[Dict]] = None,
    name: str = "construct",
    created_by: str = "Athey lab",
    add_phospho_ends: bool = True,
) -> None:
    """
    Write a SnapGene-compatible .dna file.

    Parameters
    ----------
    path        : output file path (str)
    sequence    : DNA sequence string (A/C/G/T)
    is_circular : True → circular, False → linear
    features    : list of feature dicts from make_feature()
    primers     : list of primer dicts from make_primer_feature()
    name        : construct name shown in SnapGene
    created_by  : creator string for notes packet
    add_phospho_ends : include phosphorylated-ends AdditionalProps packet
    """
    sequence = sequence.upper()
    data = _cookie_packet()
    data += _dna_packet(sequence, is_circular)
    data += _notes_packet(name, created_by)

    feat_list = list(features or [])
    if feat_list:
        data += _features_packet(feat_list)

    primer_list = list(primers or [])
    if primer_list:
        data += _primers_packet(primer_list)

    if add_phospho_ends:
        data += _additional_props_packet()

    with open(path, "wb") as fh:
        fh.write(data)
