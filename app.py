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
"""

import os
import tempfile
from pathlib import Path

import streamlit as st

from assembly_pipeline import run_pipeline, build_zip

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
Paste or upload your gene/CDS sequence and click **Run pipeline**.
The tool will:

1. **Split** the gene into 5 Parts (A, B, C, D', E) and ~200 bp fragments
2. **Design primers** for each fragment (BsaI + BsmBI tails)
3. **Write `.dna` files** for PCR products, Level 0, Level 1, and Level 2 constructs
4. Package everything into a **ZIP** for download

All `.dna` files open directly in SnapGene.
"""
)

st.divider()

# ---------------------------------------------------------------------------
# Step 1 – Gene sequence input
# ---------------------------------------------------------------------------
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
        # Strip FASTA header lines
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

# ---------------------------------------------------------------------------
# Step 2 – Options
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Step 3 – Run
# ---------------------------------------------------------------------------
st.header("Step 3 – Run the pipeline")

run_button = st.button("▶ Run pipeline", type="primary", use_container_width=True)

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
            out_dir = Path(tmpdir) / gene_name
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

                # ---- Primer table ----
                st.subheader("Designed primers")
                flag_cols = ["gc_clamp_ok", "tm_ok"]
                display_df = primer_df.copy()
                for col in flag_cols:
                    if col in display_df.columns:
                        display_df[col] = display_df[col].map(
                            {True: "✓", False: "✗"}
                        )
                st.dataframe(display_df, use_container_width=True)

                # Flag any sub-ideal primers
                bad = primer_df[
                    ~primer_df["gc_clamp_ok"] | ~primer_df["tm_ok"]
                ] if "gc_clamp_ok" in primer_df.columns else primer_df.iloc[0:0]
                if len(bad):
                    st.warning(
                        f"{len(bad)} primer(s) have sub-ideal Tm or GC clamp. "
                        "Check the **notes** column for details."
                    )

                # ---- Stats ----
                n_parts = primer_df["part"].nunique() if "part" in primer_df.columns else 5
                n_frags = len(primer_df) // 2
                st.info(
                    f"**{n_parts} Parts** | **{n_frags} fragments** | "
                    f"**{len(primer_df)} primers** total"
                )

                # ---- Download ----
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

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption(
    "Built for the Athey Lab · Golden Gate pJUMP assembly system · "
    "All .dna files compatible with SnapGene"
)
