"""
app.py – Streamlit UI for the Golden Gate assembly pipeline.

Usage:
    cd "~/Desktop/Athey Lab/gene_assembler"
    streamlit run app.py

The user pastes (or uploads) a gene sequence, clicks Run, and downloads
a ZIP containing:
  - primers CSV
  - PCR product .dna files  (linear, one per fragment)
  - Level 0 .dna files      (circular, fragment in backbone)
  - Level 1 .dna files      (circular, all fragments of each Part assembled)
  - Level 2 .dna file       (circular, all 5 Parts assembled)

A second tab handles nanopore barcode design for each fragment.
"""

import io
import os
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from assembly_pipeline import run_pipeline, build_zip


# ---------------------------------------------------------------------------
# Helpers: fragment table import and merging
# ---------------------------------------------------------------------------

def _merge_parts_lists(base_parts: list, extra_parts: list) -> list:
    """
    Merge two lists of Part objects.
    base_parts take precedence for the same (label, fragment_index).
    extra_parts add any parts/fragments not already present in base_parts.
    """
    from fragment_splitter import Part as _Part

    result: dict = {}   # label → {frag_index: Fragment}
    for part in base_parts:
        result[part.label] = {f.index: f for f in part.fragments}
    for part in extra_parts:
        if part.label not in result:
            result[part.label] = {}
        for frag in part.fragments:
            if frag.index not in result[part.label]:
                result[part.label][frag.index] = frag
            # If already present, keep the base version (it has the sequence)

    part_order = ["A", "B", "C", "Dprime", "E"]
    out = []
    for lbl in sorted(result, key=lambda x: part_order.index(x) if x in part_order else 99):
        frags = sorted(result[lbl].values(), key=lambda f: f.index)
        out.append(_Part(label=lbl, fragments=frags))
    return out


def parse_fragment_table(df: pd.DataFrame) -> list:
    """Convert a user-uploaded fragment table into a list of Part objects."""
    from fragment_splitter import Fragment, Part
    df.columns = [str(c).strip().lower() for c in df.columns]
    required = {"part", "fragment", "oh5", "oh3"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    parts_dict = {}  # label -> list of Fragment
    for _, row in df.iterrows():
        label = str(row["part"]).strip()
        idx   = int(row["fragment"])
        oh5   = str(row["oh5"]).strip().upper()
        oh3   = str(row["oh3"]).strip().upper()
        seq_raw = row.get("sequence", "")
        seq   = "".join(c for c in str(seq_raw).upper() if c in "ACGT") if str(seq_raw).upper() != "NAN" else ""
        frag  = Fragment(name=f"splitseq_{idx}", sequence=seq,
                         overhang_5=oh5, overhang_3=oh3,
                         part_label=label, index=idx)
        parts_dict.setdefault(label, []).append(frag)

    # Sort fragments within each part by index
    part_order = ["A", "B", "C", "Dprime", "D", "E"]
    labels_sorted = sorted(parts_dict.keys(), key=lambda x: part_order.index(x) if x in part_order else 99)
    result = []
    for label in labels_sorted:
        frags = sorted(parts_dict[label], key=lambda f: f.index)
        result.append(Part(label=label, fragments=frags))
    return result


# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="Golden Gate Assembler",
    page_icon="🧬",
    layout="centered",
)

st.title("🧬 Automated Golden Gate Assembly Pipeline")
st.markdown(
    """
Paste or upload your gene/CDS sequence, run the **Golden Gate** tab, then
switch to the **Nanopore Barcodes** tab to generate barcode adapter files.

All `.dna` files open directly in SnapGene.
"""
)

st.divider()

# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------
tab_gg, tab_nb = st.tabs(["⚗️ Golden Gate Design", "🧬 Nanopore Barcodes"])

# ===========================================================================
# Tab 1 – Golden Gate Design  (unchanged logic, now inside a tab)
# ===========================================================================
with tab_gg:
    st.header("Step 1 – Provide your gene sequence")

    col1, col2 = st.columns([2, 1])

    with col1:
        gene_seq_input = st.text_area(
            "Paste gene / CDS sequence (A/C/G/T, any case, spaces/newlines OK)",
            height=160,
            placeholder="ATGGACTACAAGGACCACGACGGAGACTACAAGGACCACGACATCGACTACAAGGATGACG...",
            key="gene_seq_text",
        )

    with col2:
        fasta_file = st.file_uploader(
            "or upload a FASTA file",
            type=["fa", "fasta", "txt"],
            key="fasta_upload",
        )
        if fasta_file is not None:
            raw = fasta_file.getvalue().decode("utf-8")
            seq_lines = [
                line.strip()
                for line in raw.splitlines()
                if line.strip() and not line.startswith(">")
            ]
            gene_seq_input = "".join(seq_lines)
            st.success(f"Loaded FASTA: {len(gene_seq_input)} nt")

    # Clean the input
    gene_seq = "".join(
        c for c in gene_seq_input.upper() if c in "ACGT"
    ) if gene_seq_input else ""

    if gene_seq:
        st.info(f"Sequence length: **{len(gene_seq):,} bp**")

    # ---- Options ----
    st.header("Step 2 – Options")

    col_name, col_len = st.columns(2)
    with col_name:
        gene_name = st.text_input(
            "Gene name (used in file names)",
            value="MyGene",
            help="Alphanumeric + underscores, no spaces.",
        )
        gene_name = "".join(c if c.isalnum() or c == "_" else "_" for c in gene_name) or "MyGene"

    with col_len:
        target_frag_len = st.slider(
            "Target fragment length (bp)",
            min_value=120,
            max_value=270,
            value=190,
            step=10,
            help="NEBridge uses ~200 bp. Shorter = more fragments, longer = fewer.",
        )

    # ---- Run ----
    st.header("Step 3 – Run the pipeline")

    run_button = st.button("▶ Run pipeline", type="primary", use_container_width=True, key="run_gg")

    if run_button:
        if not gene_seq:
            st.error("Please provide a gene sequence.")
        elif len(gene_seq) < 500:
            st.error(
                f"Sequence is only {len(gene_seq)} bp. "
                "The pipeline expects a full gene/CDS (at least 500 bp)."
            )
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                out_dir  = Path(tmpdir) / gene_name
                zip_path = Path(tmpdir) / f"{gene_name}_assembly.zip"

                try:
                    with st.spinner("Splitting gene, designing primers, writing .dna files..."):
                        primer_df = run_pipeline(
                            gene_seq=gene_seq,
                            output_dir=str(out_dir),
                            gene_name=gene_name,
                            target_frag_len=target_frag_len,
                        )
                        build_zip(str(out_dir), str(zip_path))

                    st.success("Pipeline complete!")

                    # Store parts in session state for the barcode tab
                    from fragment_splitter import split_gene
                    st.session_state["gg_parts"]     = split_gene(gene_seq, target_frag_len)
                    st.session_state["gg_gene_name"]  = gene_name
                    st.session_state["gg_gene_seq"]   = gene_seq
                    st.session_state["gg_target_frag"] = target_frag_len

                    # ---- Primer table ----
                    st.subheader("Designed primers")
                    flag_cols = ["gc_clamp_ok", "tm_ok"]
                    display_df = primer_df.copy()
                    for col in flag_cols:
                        if col in display_df.columns:
                            display_df[col] = display_df[col].map({True: "✓", False: "✗"})
                    st.dataframe(display_df, use_container_width=True)

                    bad = primer_df[
                        ~primer_df["gc_clamp_ok"] | ~primer_df["tm_ok"]
                    ] if "gc_clamp_ok" in primer_df.columns else primer_df.iloc[0:0]
                    if len(bad):
                        st.warning(
                            f"{len(bad)} primer(s) have sub-ideal Tm or GC clamp. "
                            "Check the **notes** column for details."
                        )

                    n_parts = primer_df["part"].nunique() if "part" in primer_df.columns else 5
                    n_frags = len(primer_df) // 2
                    st.info(
                        f"**{n_parts} Parts** | **{n_frags} fragments** | "
                        f"**{len(primer_df)} primers** total"
                    )

                    with open(zip_path, "rb") as zf:
                        zip_bytes = zf.read()

                    st.download_button(
                        label=f"⬇ Download {gene_name}_assembly.zip",
                        data=zip_bytes,
                        file_name=f"{gene_name}_assembly.zip",
                        mime="application/zip",
                        use_container_width=True,
                        type="primary",
                    )

                except Exception as exc:
                    st.error(f"Pipeline error: {exc}")
                    import traceback
                    with st.expander("Error details"):
                        st.code(traceback.format_exc())


# ===========================================================================
# Tab 2 – Nanopore Barcodes
# ===========================================================================
with tab_nb:
    st.header("Nanopore Barcode Design")
    st.markdown(
        """
Assigns ONT nanopore barcodes (NB01–NB96) to each fragment's 5' and 3' overhangs
and generates SnapGene `.dna` files for:

- **NBXX_TopStrand.dna** – 40 nt top strand oligo
- **NBXX_BottomStrand.dna** – 50 nt bottom strand oligo
- **XXXX_NBXX.dna** – annealed barcode adapter with sticky ends
- **Ligation.dna** – full ligation product (native adapters + barcodes + insert)

Run the **Golden Gate Design** tab first to load fragment data, or use the
sequence input below if you have not run it yet in this session.
        """
    )

    st.divider()

    # ---- Import from existing gene folder (ZIP) ----
    with st.expander("📁 Import from existing gene folder (ZIP)", expanded=False):
        st.markdown(
            """
Upload a **ZIP of your existing gene folder** (e.g. right-click your CYP2D7 folder
and choose *Compress*).  The tool will scan it for:

- **`XXXX_NBNN.dna` file names** → extracts which barcode is assigned to which overhang
- **`Ligation.dna` files** → reads back the fragment sequence, oh5, and oh3 directly
  from the known fixed structure of the ligation product

This lets you continue from a partially-completed gene without re-running the
pipeline or filling in a manual CSV.  Any parts not yet in the folder (e.g. Parts B–E
for CYP2D7) can be added via the **Import existing fragment data** expander below,
or via the Golden Gate Design tab.
            """
        )

        folder_zip_file = st.file_uploader(
            "Upload gene folder as ZIP",
            type=["zip"],
            key="nb_folder_zip",
        )
        folder_zip_gene_name = st.text_input(
            "Gene name",
            value="CYP2D7",
            key="nb_folder_zip_gene_name",
        )

        if folder_zip_file is not None:
            from barcode_designer import parse_gene_folder_zip
            try:
                with st.spinner("Scanning folder ZIP…"):
                    zip_parts, zip_oh_to_nb, zip_used_pairs = parse_gene_folder_zip(
                        folder_zip_file.getvalue()
                    )

                _clean = "".join(
                    c if c.isalnum() or c == "_" else "_"
                    for c in folder_zip_gene_name
                ) or "ImportedGene"

                st.session_state["zip_imported_parts"] = zip_parts
                st.session_state["zip_oh_to_nb"]       = zip_oh_to_nb
                st.session_state["zip_used_pairs"]     = zip_used_pairs
                st.session_state["gg_gene_name"]       = _clean
                st.session_state["gg_gene_seq"]        = ""

                total_z    = sum(len(p.fragments) for p in zip_parts)
                has_seq_z  = sum(1 for p in zip_parts for f in p.fragments if f.sequence)
                st.success(
                    f"Scanned ZIP: **{len(zip_parts)} Parts**, **{total_z} fragments** "
                    f"({has_seq_z} with sequences from Ligation.dna). "
                    f"Extracted **{len(zip_oh_to_nb)} overhang→barcode** assignments."
                )

                # Preview table
                _rows_z = []
                for _p in zip_parts:
                    for _f in _p.fragments:
                        _rows_z.append({
                            "Part":        _p.label,
                            "Frag #":      _f.index,
                            "oh5":         _f.overhang_5,
                            "NB_oh5":      zip_oh_to_nb.get(_f.overhang_5, "—"),
                            "oh3":         _f.overhang_3,
                            "NB_oh3":      zip_oh_to_nb.get(_f.overhang_3, "—"),
                            "Has sequence": "✅" if _f.sequence else "❌",
                        })
                st.dataframe(pd.DataFrame(_rows_z), use_container_width=True)

            except Exception as _ze:
                st.error(f"Failed to parse folder ZIP: {_ze}")
                import traceback as _tb
                with st.expander("Error details"):
                    st.code(_tb.format_exc())

    # ---- Gene sequence for standalone use ----
    with st.expander("Provide gene sequence (if not already loaded from GG tab)", expanded=False):
        nb_seq_input = st.text_area(
            "Gene / CDS sequence",
            height=120,
            placeholder="ATGGAC...",
            key="nb_gene_seq_text",
        )
        nb_fasta = st.file_uploader(
            "or upload FASTA",
            type=["fa", "fasta", "txt"],
            key="nb_fasta_upload",
        )
        if nb_fasta is not None:
            raw = nb_fasta.getvalue().decode("utf-8")
            seq_lines = [
                line.strip()
                for line in raw.splitlines()
                if line.strip() and not line.startswith(">")
            ]
            nb_seq_input = "".join(seq_lines)
            st.success(f"Loaded FASTA: {len(nb_seq_input)} nt")

        nb_gene_name_input = st.text_input(
            "Gene name",
            value="MyGene",
            key="nb_gene_name_input",
        )
        nb_target_frag = st.slider(
            "Target fragment length (bp)",
            min_value=120, max_value=270, value=190, step=10,
            key="nb_frag_len",
        )
        load_seq_btn = st.button("Load sequence for barcoding", key="nb_load_seq")
        if load_seq_btn:
            nb_seq_clean = "".join(c for c in nb_seq_input.upper() if c in "ACGT")
            if len(nb_seq_clean) < 500:
                st.error("Sequence too short (need ≥ 500 bp).")
            else:
                from fragment_splitter import split_gene
                _clean_name = "".join(
                    c if c.isalnum() or c == "_" else "_"
                    for c in nb_gene_name_input
                ) or "MyGene"
                st.session_state["gg_parts"]      = split_gene(nb_seq_clean, nb_target_frag)
                st.session_state["gg_gene_name"]  = _clean_name
                st.session_state["gg_gene_seq"]   = nb_seq_clean
                st.session_state["gg_target_frag"] = nb_target_frag
                st.success(f"Loaded {len(nb_seq_clean):,} bp — {sum(len(p.fragments) for p in st.session_state['gg_parts'])} fragments across {len(st.session_state['gg_parts'])} Parts.")

    # ---- Import existing fragment data ----
    with st.expander("Import existing fragment data (CSV/Excel)", expanded=False):
        st.markdown(
            "Use this if your fragments were cut at different positions than the automatic "
            "splitter would choose (e.g. an existing wet-lab construct). "
            "The `sequence` column is optional — if omitted, Ligation.dna files will be skipped."
        )

        template_csv = "part,fragment,oh5,oh3,sequence\nA,1,GGAG,AATG,ATGCGATCGATC...\n".encode()
        st.download_button(
            label="⬇ Download template CSV",
            data=template_csv,
            file_name="fragment_template.csv",
            mime="text/csv",
        )

        nb_import_file = st.file_uploader(
            "Upload fragment table (part, fragment, oh5, oh3, sequence)",
            type=["csv", "xlsx", "xls"],
            key="nb_frag_import",
        )

        nb_import_gene_name = st.text_input(
            "Gene name for imported data",
            value="ImportedGene",
            key="nb_import_gene_name",
        )

        if nb_import_file is not None:
            try:
                fname = nb_import_file.name.lower()
                if fname.endswith(".csv"):
                    import_df = pd.read_csv(nb_import_file)
                else:
                    import_df = pd.read_excel(nb_import_file)

                imported_parts = parse_fragment_table(import_df)

                # Merge with any existing ZIP-imported parts
                _zip_base = st.session_state.get("zip_imported_parts", [])
                merged    = _merge_parts_lists(_zip_base, imported_parts)

                st.session_state["gg_parts"]    = merged
                st.session_state["gg_gene_seq"] = ""
                _clean_import_name = "".join(
                    c if c.isalnum() or c == "_" else "_"
                    for c in nb_import_gene_name
                ) or "ImportedGene"
                st.session_state["gg_gene_name"] = _clean_import_name

                n_frags_imported = sum(len(p.fragments) for p in merged)
                n_parts_imported = len(merged)
                _note = " (merged with folder ZIP data)" if _zip_base else ""
                st.success(
                    f"Imported {n_frags_imported} fragments across "
                    f"{n_parts_imported} Parts{_note}."
                )
                st.dataframe(import_df, use_container_width=True)

            except Exception as import_exc:
                st.error(f"Import error: {import_exc}")
                import traceback
                with st.expander("Import error details"):
                    st.code(traceback.format_exc())

    # ---- History Excel upload ----
    st.subheader("Barcode history (optional)")
    history_file = st.file_uploader(
        "Upload history Excel (columns: '4-bp Overhang', 'Barcode label'). "
        "Existing assignments will be reused; leave blank to treat all barcodes as new.",
        type=["xlsx", "xls"],
        key="nb_history_upload",
    )

    # ---- Resolve active parts: merge ZIP data with GG/CSV data ----
    _zip_parts = st.session_state.get("zip_imported_parts", [])
    _gg_parts  = st.session_state.get("gg_parts", [])
    parts = _merge_parts_lists(_zip_parts, _gg_parts) if (_zip_parts or _gg_parts) else None

    if not parts:
        st.info(
            "No fragment data loaded yet. Upload a **gene folder ZIP** above, "
            "run the **Golden Gate Design** tab, or import a fragment CSV."
        )
    else:
        gene_name_nb = st.session_state.get("gg_gene_name", "MyGene")
        total_frags  = sum(len(p.fragments) for p in parts)
        frags_with_seq = sum(1 for p in parts for f in p.fragments if f.sequence)
        _seq_note = (
            f" — {frags_with_seq} have sequences (Ligation.dna will be generated for those)"
            if frags_with_seq < total_frags else ""
        )
        st.success(
            f"Fragment data loaded: **{len(parts)} Parts**, **{total_frags} fragments** "
            f"(gene: {gene_name_nb}){_seq_note}"
        )

        # ---- Run barcode assignment ----
        run_nb_btn = st.button(
            "▶ Assign barcodes & generate files",
            type="primary",
            use_container_width=True,
            key="run_nb",
        )

        if run_nb_btn:
            from barcode_designer import (
                load_history, assign_barcodes,
                generate_all_barcode_files, build_barcode_excel,
                build_barcode_zip, build_combined_zip,
            )

            try:
                # ---- Build history from all sources ----
                # Source 1: folder ZIP (parsed from file names / Ligation.dna)
                oh_to_nb   = dict(st.session_state.get("zip_oh_to_nb",   {}))
                used_pairs = set(st.session_state.get("zip_used_pairs",  set()))

                # Source 2: Excel history upload (merges; Excel takes precedence on conflicts)
                if history_file is not None:
                    with st.spinner("Reading history Excel..."):
                        xl_oh_to_nb, xl_used_pairs = load_history(
                            io.BytesIO(history_file.getvalue())
                        )
                    oh_to_nb.update(xl_oh_to_nb)
                    used_pairs |= xl_used_pairs
                    st.info(
                        f"History: {len(oh_to_nb)} overhang→barcode assignments "
                        f"({len(xl_oh_to_nb)} from Excel"
                        + (f", {len(st.session_state.get('zip_oh_to_nb', {}))} from folder ZIP" if st.session_state.get("zip_oh_to_nb") else "")
                        + f"), {len(used_pairs)} used pairs."
                    )
                elif oh_to_nb:
                    st.info(
                        f"History from folder ZIP: {len(oh_to_nb)} overhang→barcode "
                        f"assignments, {len(used_pairs)} used pairs."
                    )

                with st.spinner("Assigning barcodes and building .dna files..."):
                    # We need a fresh copy of oh_to_nb / used_pairs for display
                    oh_to_nb_copy   = dict(oh_to_nb)
                    used_pairs_copy = set(used_pairs)

                    assignments = assign_barcodes(
                        parts, oh_to_nb_copy, used_pairs_copy
                    )

                    # Generate all .dna file bytes
                    oh_to_nb_gen   = dict(oh_to_nb)
                    used_pairs_gen = set(used_pairs)
                    barcode_files  = generate_all_barcode_files(
                        parts, oh_to_nb_gen, used_pairs_gen
                    )

                    # Build barcode Excel
                    excel_bytes = build_barcode_excel(parts, assignments)

                st.success(
                    f"Barcode assignment complete! "
                    f"{len(barcode_files)} .dna files generated."
                )

                # ---- Assignment table ----
                st.subheader("Barcode assignments")
                rows = []
                for part in parts:
                    for frag in part.fragments:
                        key  = (part.label, frag.index)
                        asgn = assignments.get(key, {})
                        rows.append({
                            "Fragment":   f"{part.label}.{frag.index}",
                            "Part":       part.label,
                            "oh5":        asgn.get("oh5", ""),
                            "NB_oh5":     asgn.get("nb_oh5", ""),
                            "oh5 status": "✅ reused" if asgn.get("oh5_reused") else "🆕 new",
                            "oh3":        asgn.get("oh3", ""),
                            "NB_oh3":     asgn.get("nb_oh3", ""),
                            "oh3 status": "✅ reused" if asgn.get("oh3_reused") else "🆕 new",
                        })
                asgn_df = pd.DataFrame(rows)
                st.dataframe(asgn_df, use_container_width=True)

                n_reused_5 = sum(1 for r in rows if "reused" in r["oh5 status"])
                n_reused_3 = sum(1 for r in rows if "reused" in r["oh3 status"])
                st.caption(
                    f"oh5: {n_reused_5} reused, {len(rows)-n_reused_5} new  |  "
                    f"oh3: {n_reused_3} reused, {len(rows)-n_reused_3} new"
                )

                st.divider()

                # ---- Download buttons ----
                st.subheader("Downloads")

                col_dl1, col_dl2, col_dl3 = st.columns(3)

                with col_dl1:
                    barcode_zip_bytes = build_barcode_zip(barcode_files)
                    st.download_button(
                        label="⬇ Barcode .dna Files",
                        data=barcode_zip_bytes,
                        file_name=f"{gene_name_nb}_barcodes.zip",
                        mime="application/zip",
                        use_container_width=True,
                        type="primary",
                    )

                with col_dl2:
                    st.download_button(
                        label="⬇ Barcode Summary Excel",
                        data=excel_bytes,
                        file_name=f"{gene_name_nb}_barcode_assignments.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True,
                    )

                with col_dl3:
                    # Combined GG + barcode ZIP — re-run GG pipeline to get the
                    # output directory, then merge with barcode files.
                    combined_zip_bytes = None
                    gg_seq = st.session_state.get("gg_gene_seq", "")
                    gg_tfl = st.session_state.get("gg_target_frag", 190)
                    if gg_seq:
                        try:
                            with tempfile.TemporaryDirectory() as tmpdir2:
                                gg_out = Path(tmpdir2) / gene_name_nb
                                run_pipeline(
                                    gene_seq=gg_seq,
                                    output_dir=str(gg_out),
                                    gene_name=gene_name_nb,
                                    target_frag_len=gg_tfl,
                                )
                                combined_zip_bytes = build_combined_zip(
                                    str(gg_out), barcode_files
                                )
                        except Exception:
                            pass  # fallback: GG output unavailable

                    if combined_zip_bytes:
                        st.download_button(
                            label="⬇ All Files Together",
                            data=combined_zip_bytes,
                            file_name=f"{gene_name_nb}_all_files.zip",
                            mime="application/zip",
                            use_container_width=True,
                        )
                    else:
                        st.info(
                            "Run the **Golden Gate Design** tab first (or load a sequence "
                            "above) to enable the combined GG + barcode download."
                        )

            except Exception as exc:
                st.error(f"Barcode design error: {exc}")
                import traceback
                with st.expander("Error details"):
                    st.code(traceback.format_exc())

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption(
    "Built for the Athey Lab · Golden Gate pJUMP assembly system · "
    "All .dna files compatible with SnapGene"
)
