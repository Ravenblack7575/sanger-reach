# Sanger Primer Reach

A command-line tool that maps Sanger sequencing primers onto a consensus sequence and generates a visual PDF report showing binding sites and estimated read coverage.

Originally prototyped as a Jupyter notebook, this script wraps the same logic into a single callable function — no notebook required.

---

## Features

- Finds primer binding sites on a consensus FASTA sequence
- Supports **forward and reverse primers**
- Configurable **mismatch tolerance** (0–5 mismatches)
- Reports **all binding sites** per primer (exact and mismatch hits)
- Generates a **two-page PDF report**:
  - Page 1: visual map with directional arrows and dotted Sanger reach lines
  - Page 2: summary table with binding positions, mismatch counts, and coverage ranges
- Colour-coded output — exact matches, mismatch hits, and not-found primers are visually distinct

---

## Requirements

- Python 3.10+
- [Biopython](https://biopython.org/)
- [Matplotlib](https://matplotlib.org/)
- [ReportLab](https://www.reportlab.com/)

Install with pip:

```bash
pip install biopython matplotlib reportlab
```

Or with Poetry:

```bash
poetry add biopython matplotlib reportlab
```

---

## Input Files

### Primer CSV

A CSV file with one primer per row. A header row is optional and auto-detected.

| Column | Description |
|---|---|
| `primer_name` | Name of the primer |
| `primer_sequence` | Primer sequence written 5' → 3' |
| `is_reverse` | `TRUE` for reverse primers, `FALSE` for forward primers |

**Example:**

```csv
primer_name,primer_sequence,is_reverse
Primer_F1,ATGTCGTCCTTGATCGG,False
Primer_R1,TTTAGCCCCTTATACA,True
Primer_F2,GCACCTTAGGGCGAGTCCTA,False
```

> **Note on reverse primers:** the `primer_sequence` column should contain the primer sequence as ordered from the supplier (5' → 3'). The script automatically searches for the reverse complement on the consensus strand.

### Consensus FASTA

A standard FASTA file containing a single sequence record.

```
>SEQUENCE
ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTCAT...
```

---

## Usage

```bash
python sanger_primer_reach.py \
    --primers primers.csv \
    --fasta sequence.fasta \
    --output report.pdf \
    --sanger-read-length 800 \
    --max-mismatches 0
```

### Arguments

| Argument | Short | Default | Description |
|---|---|---|---|
| `--primers` | `-p` | *(required)* | Path to the primer CSV file |
| `--fasta` | `-f` | *(required)* | Path to the consensus FASTA file |
| `--output` | `-o` | `sanger_primer_report.pdf` | Path for the output PDF |
| `--sanger-read-length` | `-s` | `800` | Expected Sanger read length in bp |
| `--max-mismatches` | `-m` | `0` | Allowed mismatches per primer (0–5) |

---

## Mismatch Tolerance

By default the tool requires an **exact match** (`--max-mismatches 0`). Setting a higher value allows the sliding-window search to tolerate substitutions anywhere in the primer — useful when screening primers against a slightly divergent consensus.

```bash
# Allow up to 2 mismatches
python sanger_primer_reach.py -p primers.csv -f seq.fasta -o report.pdf -m 2
```

**Behaviour with mismatches enabled:**

- All positions where the primer matches within the mismatch threshold are reported.
- Mismatch hits are displayed with **faded arrows** on the map (more mismatches = more transparent).
- Arrow labels include the mismatch count, e.g. `Primer_F1 (2mm)`.
- The summary table has a dedicated **Mismatches** column.
- Status colours: 🟢 Exact match · 🟠 Mismatch hit · 🔴 NOT FOUND

> ⚠️ High mismatch values (4–5) can produce many spurious hits, especially for short primers. Start with 1–2 and increase only if needed.

---

## Output PDF

### Page 1 — Visual Map

- Black horizontal line = consensus sequence (5' → 3', left to right)
- **Green arrows** = forward primers (pointing right →)
- **Blue arrows** = reverse primers (pointing left ←)
- **Dashed lines** extending from each arrow = estimated Sanger read reach
- Faded arrows = mismatch hits; full-opacity arrows = exact matches
- Position tick marks every ~5% of the sequence length

### Page 2 — Summary Table

| Column | Description |
|---|---|
| Primer Name | Name from the CSV |
| Type | Forward or Reverse |
| Sequence (5'→3') | Primer sequence as supplied |
| Binding site (bp) | Start–end positions on the consensus |
| Mismatches | Number of mismatches at this site |
| Seq. coverage (bp) | Estimated range covered by Sanger sequencing |
| Status | Exact match / N mismatch(es) / NOT FOUND |

---



## How It Works

1. **CSV loading** — primers are read and the header row is auto-detected based on whether the third column contains `true`/`false`.
2. **Binding site search** — a sliding window of the primer length scans the consensus. For reverse primers, the reverse complement is searched. At each position, mismatches are counted; positions within the threshold are recorded.
3. **Figure generation** — Matplotlib draws the backbone, arrows, and reach lines. Arrows are stacked vertically to avoid overlap when multiple primers bind nearby.
4. **PDF assembly** — ReportLab composes the figure (as an embedded PNG) and the summary table into a landscape A4 PDF.

---

## Project Structure

```
.
├── sanger_primer_reach.py   # Main script
├── primers.csv              # Your primer file (user-supplied)
├── sequence.fasta           # Your consensus sequence (user-supplied)
└── README.md
```

---

## Limitations

- Only **substitution mismatches** are considered — insertions and deletions are not handled.
- Only the **first FASTA record** in the file is used if the file contains multiple sequences.
- High mismatch settings may produce many off-target hits; review the summary table carefully.

---

## Citation

Peter J. A. Cock, Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, Michiel J. L. de Hoon, **Biopython**: freely available Python tools for computational molecular biology and bioinformatics, Bioinformatics, Volume 25, Issue 11, June 2009, Pages 1422–1423, https://doi.org/10.1093/bioinformatics/btp163

Anthropic. (2025). Claude (Claude Sonnet 4.6) [Large language model]. https://www.anthropic.com

## License

MIT
