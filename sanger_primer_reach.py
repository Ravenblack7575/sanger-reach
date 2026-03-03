"""
sanger_primer_reach.py

Generates a visual PDF report showing Sanger sequencing primer binding sites
and estimated read coverage on a consensus FASTA sequence.

Usage:
    python sanger_primer_reach.py \
        --primers primers.csv \
        --fasta sequence.fasta \
        --output report.pdf \
        [--sanger-read-length 800] \
        [--max-mismatches 2]

CSV format (with or without header row):
    primer_name, primer_sequence, is_reverse
    Primer_F1,  ATGCGTAAA,        False
    Primer_R1,  TTTAGCCAT,        True

    primer_name     : Name of the primer
    primer_sequence : Primer sequence from 5' to 3'
    is_reverse      : TRUE / FALSE
                      Forward primers → FALSE
                      Reverse primers → TRUE

Dependencies:
    biopython, matplotlib, reportlab
    pip install biopython matplotlib reportlab
"""

import argparse
import csv
import io

import matplotlib
matplotlib.use("Agg")          # non-interactive backend – safe for scripts
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow

from Bio import SeqIO
from Bio.Seq import Seq

from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle,
)
from reportlab.lib import colors


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def load_primers_from_csv(csv_file: str) -> list[tuple[str, str, bool]]:
    """
    Load primers from a CSV file.

    Returns a list of (primer_name, primer_sequence, is_reverse) tuples.
    """
    primers = []
    with open(csv_file, newline="") as fh:
        reader = csv.reader(fh)
        first_row = next(reader)
        # Detect whether the first row is a header or data
        if first_row[2].strip().lower() not in ("true", "false"):
            pass  # It's a header row – skip it
        else:
            name = first_row[0].strip()
            seq  = first_row[1].strip().upper()
            rev  = first_row[2].strip().lower() == "true"
            primers.append((name, seq, rev))

        for row in reader:
            if not row or len(row) < 3:
                continue
            name = row[0].strip()
            seq  = row[1].strip().upper()
            rev  = row[2].strip().lower() == "true"
            primers.append((name, seq, rev))

    return primers


def load_fasta(fasta_file: str) -> tuple[str, str]:
    """
    Read the first record from a FASTA file.

    Returns (record_id, sequence_string).
    """
    record = SeqIO.read(fasta_file, "fasta")
    return record.id, str(record.seq)


# ---------------------------------------------------------------------------
# Primer binding logic (with mismatch tolerance)
# ---------------------------------------------------------------------------

def _count_mismatches(seq_a: str, seq_b: str) -> int:
    """Count mismatches between two equal-length strings."""
    return sum(a != b for a, b in zip(seq_a, seq_b))


def find_primer_binding_sites(
    consensus_seq: str,
    primer_seq: str,
    is_reverse: bool = False,
    max_mismatches: int = 0,
) -> list[tuple[int, int, str, int]]:
    """
    Find all binding sites for a primer on the consensus sequence,
    allowing up to `max_mismatches` mismatches anywhere in the primer.

    For reverse primers, searches for the reverse complement of the primer
    sequence in the forward strand of the consensus.

    Returns a list of (start, end, strand, n_mismatches) tuples,
    sorted by position (ties broken by fewest mismatches).
    """
    consensus = consensus_seq.upper().replace(" ", "").replace("\n", "")
    primer    = primer_seq.upper().replace(" ", "").replace("\n", "")

    if is_reverse:
        search_seq = str(Seq(primer).reverse_complement())
        strand = "reverse"
    else:
        search_seq = primer
        strand = "forward"

    plen = len(search_seq)
    clen = len(consensus)
    binding_sites = []

    for i in range(clen - plen + 1):
        window = consensus[i : i + plen]
        mm = _count_mismatches(window, search_seq)
        if mm <= max_mismatches:
            binding_sites.append((i, i + plen, strand, mm))

    # Sort by position; ties broken by fewest mismatches
    binding_sites.sort(key=lambda x: (x[0], x[3]))
    return binding_sites


# ---------------------------------------------------------------------------
# Figure / plot
# ---------------------------------------------------------------------------

FORWARD_COLOR = "#00CC44"   # green
REVERSE_COLOR = "#0044FF"   # blue


def build_primer_map_figure(
    consensus_seq: str,
    primers: list[tuple[str, str, bool]],
    sequence_id: str = "",
    sanger_read_length: int = 800,
    max_mismatches: int = 0,
    figsize: tuple[float, float] = (14, 6),
) -> plt.Figure:
    """
    Build and return a Matplotlib Figure showing primer binding sites
    and estimated Sanger sequencing reach.
    """
    seq_length = len(consensus_seq)

    fig, ax = plt.subplots(figsize=figsize)

    # ── Consensus sequence backbone ──────────────────────────────────────────
    ax.plot([0, seq_length], [0, 0], "k-", linewidth=2, label="Consensus Sequence")

    tick_interval = max(100, seq_length // 20)
    for tick in range(0, seq_length + 1, tick_interval):
        ax.plot([tick, tick], [-0.1, 0.1], "k-", linewidth=1)
        ax.text(tick, -0.3, str(tick), ha="center", va="top", fontsize=8)

    # ── Primer arrows ────────────────────────────────────────────────────────
    y_positions: list[float] = []

    for primer_name, primer_seq, is_reverse in primers:
        sites = find_primer_binding_sites(
            consensus_seq, primer_seq, is_reverse, max_mismatches
        )
        if not sites:
            continue

        color = REVERSE_COLOR if is_reverse else FORWARD_COLOR

        for start, end, strand, n_mm in sites:
            y_offset = -0.5
            while any(abs(y_offset - y) < 0.6 for y in y_positions):
                y_offset -= 0.6
            y_positions.append(y_offset)

            # Label shows mismatch count if > 0
            label = primer_name if n_mm == 0 else f"{primer_name} ({n_mm}mm)"
            # Fade arrows proportionally for mismatch hits
            alpha_arrow = 0.8 if n_mm == 0 else max(0.35, 0.8 - n_mm * 0.1)
            alpha_line  = 0.6 if n_mm == 0 else max(0.25, 0.6 - n_mm * 0.1)

            if strand == "forward":
                arrow = FancyArrow(
                    start, y_offset, end - start, 0,
                    width=0.3, head_width=0.4,
                    head_length=min(20, (end - start) * 0.2),
                    fc=color, ec="black", linewidth=1.5, alpha=alpha_arrow,
                )
                ax.add_patch(arrow)
                sanger_end = min(end + sanger_read_length, seq_length)
                ax.plot([end, sanger_end], [y_offset, y_offset],
                        "--", color=color, linewidth=2, alpha=alpha_line)
                ax.text(start, y_offset - 0.5, label,
                        fontsize=9, ha="left", va="top", fontweight="bold")
            else:
                arrow = FancyArrow(
                    end, y_offset, start - end, 0,
                    width=0.3, head_width=0.4,
                    head_length=min(20, (end - start) * 0.2),
                    fc=color, ec="black", linewidth=1.5, alpha=alpha_arrow,
                )
                ax.add_patch(arrow)
                sanger_start = max(start - sanger_read_length, 0)
                ax.plot([sanger_start, start], [y_offset, y_offset],
                        "--", color=color, linewidth=2, alpha=alpha_line)
                ax.text(end, y_offset - 0.5, label,
                        fontsize=9, ha="right", va="top", fontweight="bold")

    # ── Axes formatting ──────────────────────────────────────────────────────
    ax.set_xlim(-seq_length * 0.05, seq_length * 1.05)
    max_y = max((abs(y) for y in y_positions), default=0)
    ax.set_ylim(-(max_y + 2), 2)
    ax.set_xlabel("Position (bp)", fontsize=12, fontweight="bold")

    mm_note = f"Max mismatches allowed: {max_mismatches}"
    title_lines = [
        "Primer Binding Sites and Sanger Sequencing Coverage",
        f"Sequence: {sequence_id}   |   Length: {seq_length} bp   |   "
        f"Sanger read length: {sanger_read_length} bp   |   {mm_note}",
    ]
    ax.set_title("\n".join(title_lines), fontsize=13, fontweight="bold")

    ax.set_yticks([])
    for spine in ("left", "right", "top"):
        ax.spines[spine].set_visible(False)

    legend_elements = [
        mpatches.Patch(color=FORWARD_COLOR, label="Forward primer (exact)"),
        mpatches.Patch(color=REVERSE_COLOR, label="Reverse primer (exact)"),
        mpatches.Patch(color=FORWARD_COLOR, alpha=0.4, label="Forward primer (mismatch)"),
        mpatches.Patch(color=REVERSE_COLOR, alpha=0.4, label="Reverse primer (mismatch)"),
        plt.Line2D([0], [0], linestyle="--", color="gray", linewidth=2,
                   label=f"Sanger reach (~{sanger_read_length} bp)"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=9)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Summary table data
# ---------------------------------------------------------------------------

def build_summary_rows(
    consensus_seq: str,
    primers: list[tuple[str, str, bool]],
    sanger_read_length: int,
    max_mismatches: int,
) -> list[list[str]]:
    """
    Return a list of rows (each a list of strings) summarising binding results.
    The first row is the header.
    """
    header = [
        "Primer Name", "Type", "Sequence (5'→3')",
        "Binding site (bp)", "Mismatches", "Seq. coverage (bp)", "Status",
    ]
    rows = [header]

    for primer_name, primer_seq, is_reverse in primers:
        primer_type = "Reverse" if is_reverse else "Forward"
        sites = find_primer_binding_sites(
            consensus_seq, primer_seq, is_reverse, max_mismatches
        )

        if not sites:
            rows.append([
                primer_name, primer_type, primer_seq, "—", "—", "—", "NOT FOUND",
            ])
            continue

        for start, end, strand, n_mm in sites:
            if strand == "forward":
                cov_end  = min(end + sanger_read_length, len(consensus_seq))
                coverage = f"{end}–{cov_end}"
            else:
                cov_start = max(start - sanger_read_length, 0)
                coverage  = f"{cov_start}–{start}"

            status = "Exact match" if n_mm == 0 else f"{n_mm} mismatch{'es' if n_mm > 1 else ''}"
            rows.append([
                primer_name, primer_type, primer_seq,
                f"{start}–{end}", str(n_mm), coverage, status,
            ])

    return rows


# ---------------------------------------------------------------------------
# PDF generation
# ---------------------------------------------------------------------------

def _fig_to_image_flowable(fig: plt.Figure, width_cm: float = 25) -> Image:
    """Render a Matplotlib figure into a ReportLab Image flowable."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    plt.close(fig)

    img = Image(buf)
    aspect = img.imageHeight / img.imageWidth
    img_width  = width_cm * cm
    img_height = img_width * aspect
    img.drawWidth  = img_width
    img.drawHeight = img_height
    return img


def generate_pdf(
    csv_file: str,
    fasta_file: str,
    output_pdf: str,
    sanger_read_length: int = 800,
    max_mismatches: int = 0,
) -> None:
    """
    Main entry point: load inputs, build the figure and summary table,
    then write a PDF report.
    """
    # ── Load data ────────────────────────────────────────────────────────────
    primers     = load_primers_from_csv(csv_file)
    seq_id, seq = load_fasta(fasta_file)

    print(f"Loaded {len(primers)} primers from '{csv_file}'")
    print(f"Loaded sequence '{seq_id}' ({len(seq):,} bp) from '{fasta_file}'")
    print(f"Max mismatches allowed: {max_mismatches}")

    # ── Build figure ─────────────────────────────────────────────────────────
    fig = build_primer_map_figure(
        seq, primers,
        sequence_id=seq_id,
        sanger_read_length=sanger_read_length,
        max_mismatches=max_mismatches,
    )
    map_image = _fig_to_image_flowable(fig, width_cm=24)

    # ── Build summary table ──────────────────────────────────────────────────
    table_data = build_summary_rows(seq, primers, sanger_read_length, max_mismatches)

    col_widths = [3*cm, 2*cm, 5*cm, 3*cm, 2.2*cm, 3*cm, 3*cm]
    tbl = Table(table_data, colWidths=col_widths, repeatRows=1)
    tbl.setStyle(TableStyle([
        ("BACKGROUND",     (0, 0), (-1, 0),  colors.HexColor("#2B4590")),
        ("TEXTCOLOR",      (0, 0), (-1, 0),  colors.white),
        ("FONTNAME",       (0, 0), (-1, 0),  "Helvetica-Bold"),
        ("FONTSIZE",       (0, 0), (-1, 0),  9),
        ("ALIGN",          (0, 0), (-1, -1), "CENTER"),
        ("VALIGN",         (0, 0), (-1, -1), "MIDDLE"),
        ("FONTNAME",       (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE",       (0, 1), (-1, -1), 8),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.white, colors.HexColor("#EEF2FF")]),
        ("GRID",           (0, 0), (-1, -1), 0.5, colors.HexColor("#AAAAAA")),
    ]))

    # Colour Status column per row
    for i, row in enumerate(table_data[1:], start=1):
        status = row[-1]
        if status == "NOT FOUND":
            tbl.setStyle(TableStyle([
                ("TEXTCOLOR", (6, i), (6, i), colors.red),
                ("FONTNAME",  (0, i), (-1, i), "Helvetica-Bold"),
            ]))
        elif status == "Exact match":
            tbl.setStyle(TableStyle([
                ("TEXTCOLOR", (6, i), (6, i), colors.HexColor("#007700")),
            ]))
        else:
            # mismatch hit – amber
            tbl.setStyle(TableStyle([
                ("TEXTCOLOR", (6, i), (6, i), colors.HexColor("#BB6600")),
            ]))

    # ── Assemble PDF ─────────────────────────────────────────────────────────
    styles = getSampleStyleSheet()
    doc    = SimpleDocTemplate(
        output_pdf,
        pagesize=landscape(A4),
        leftMargin=1.5*cm, rightMargin=1.5*cm,
        topMargin=1.5*cm,  bottomMargin=1.5*cm,
    )

    story = [
        Paragraph("Sanger Primer Reach Report", styles["Title"]),
        Spacer(1, 0.3*cm),
        Paragraph(
            f"<b>Input sequence:</b> {seq_id} &nbsp;&nbsp; "
            f"<b>Length:</b> {len(seq):,} bp &nbsp;&nbsp; "
            f"<b>Primers loaded:</b> {len(primers)} &nbsp;&nbsp; "
            f"<b>Sanger read length:</b> {sanger_read_length} bp &nbsp;&nbsp; "
            f"<b>Max mismatches:</b> {max_mismatches}",
            styles["Normal"],
        ),
        Spacer(1, 0.5*cm),
        map_image,
        PageBreak(),
        Paragraph("Primer Binding Summary", styles["Heading1"]),
        Spacer(1, 0.2*cm),
        Paragraph(
            "Status colour key: "
            "<font color='#007700'>Exact match</font> &nbsp;|&nbsp; "
            "<font color='#BB6600'>Mismatch hit</font> &nbsp;|&nbsp; "
            "<font color='red'>NOT FOUND</font>",
            styles["Normal"],
        ),
        Spacer(1, 0.4*cm),
        tbl,
    ]

    doc.build(story)
    print(f"PDF written to: {output_pdf}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a Sanger primer reach PDF report from a CSV primer "
            "file and a FASTA consensus sequence."
        )
    )
    parser.add_argument(
        "--primers", "-p", required=True,
        metavar="PRIMERS_CSV",
        help="Path to the CSV file containing primer information.",
    )
    parser.add_argument(
        "--fasta", "-f", required=True,
        metavar="FASTA_FILE",
        help="Path to the FASTA file containing the consensus sequence.",
    )
    parser.add_argument(
        "--output", "-o", default="sanger_primer_report.pdf",
        metavar="OUTPUT_PDF",
        help="Path for the output PDF report (default: sanger_primer_report.pdf).",
    )
    parser.add_argument(
        "--sanger-read-length", "-s", type=int, default=800,
        metavar="BP",
        help="Expected Sanger sequencing read length in bp (default: 800).",
    )
    parser.add_argument(
        "--max-mismatches", "-m", type=int, default=0,
        choices=range(0, 6),
        metavar="{0-5}",
        help=(
            "Maximum number of mismatches allowed when searching for primer "
            "binding sites (0–5, default: 0 = exact match only)."
        ),
    )
    args = parser.parse_args()

    generate_pdf(
        csv_file=args.primers,
        fasta_file=args.fasta,
        output_pdf=args.output,
        sanger_read_length=args.sanger_read_length,
        max_mismatches=args.max_mismatches,
    )


if __name__ == "__main__":
    main()
