<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/" style="text-decoration: none !important;">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/main/oligopool/img/logo.svg"  alt="Oligopool Calculator" width="460" class="center"/>
    </a>
</h1>

<h4><p align="center">Version: 2026.02.01</p></h4>

<p align="center">
  <a href="#features" style="text-decoration: none !important;">✨ Features</a> •
  <a href="#installation" style="text-decoration: none !important;">📦 Installation</a> •
  <a href="#getting-started" style="text-decoration: none !important;">🚀 Getting Started</a> •
  <a href="https://github.com/ayaanhossain/oligopool/blob/master/docs/docs.md" style="text-decoration: none !important;">📚 Docs</a> •
  <a href="https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md" style="text-decoration: none !important;">📋 API</a> •
  <a href="#command-line-interface-cli" style="text-decoration: none !important;">💻 CLI</a> •
  <a href="#citation" style="text-decoration: none !important;">📖 Citation</a> •
  <a href="#license" style="text-decoration: none !important;">⚖️ License</a>
</p>

`Oligopool Calculator` is a Swiss-army knife for [oligo pool libraries](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9300125/): a unified toolkit for high-throughput design, assembly, compression, and analysis of massively parallel assays, designed to integrate seamlessly with Python, the CLI, Jupyter, containers, and AI-assisted workflows.

Design modules generate primers, barcodes, motifs/anchors, and spacers; assembly modules split/pad long constructs; Degenerate Mode compresses similar sequences into IUPAC-degenerate oligos for cost-efficient synthesis (often useful for selection assays); and Analysis Mode packs and counts barcoded reads for activity quantification.

`Oligopool Calculator` has been used to build libraries of tens of thousands of promoters (see [here](https://www.nature.com/articles/s41467-022-32829-5), and [here](https://www.nature.com/articles/s41587-020-0584-2)), ribozymes, and mRNA stability elements (see [here](https://www.nature.com/articles/s41467-024-54059-7)). It has been benchmarked to design pools containing millions of oligos and to process hundreds of millions of sequencing reads per hour on low-cost desktop-grade hardware.

To learn more, please check out [our paper in ACS Synthetic Biology](https://pubs.acs.org/doi/10.1021/acssynbio.4c00661).

<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/" style="text-decoration: none !important;">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/refs/heads/main/oligopool/img/workflow.svg"  alt="Oligopool Calculator Workflow" width="3840" class="center"/>
    </a>
</h1>

**Design and analysis of oligo pool variants using `Oligopool Calculator`.** **(a)** In `Design Mode`, `Oligopool Calculator` generates optimized `barcode`s, `primer`s, `spacer`s, and `motif`s. `Assembly Mode` can `split` longer oligos into shorter `pad`ded fragments for synthesis and assembly. `Degenerate Mode` can `compress` similar variants into IUPAC-degenerate oligos for cost-efficient synthesis (useful for selection-based discovery workflows). **(b)** Once the library is assembled and cloned, barcoded amplicon sequencing data can be processed via `Analysis Mode` for characterization. `Analysis Mode` proceeds by first `index`ing one or more sets of barcodes, `pack`ing the reads, and then producing count matrices either using `acount` (association counting) or `xcount` (combinatorial counting).

<a id="features"></a>
## ✨ Features

- 🧬 **Design mode:** constraint-based design of barcodes, primers, motifs/anchors, and spacers with background screening and utilities (`barcode`, `primer`, `motif`, `spacer`, `background`, `merge`, `revcomp`, `lenstat`, `verify`, `final`).
- 🔧 **Assembly mode:** fragment long oligos into overlapping pieces and add Type IIS primer pads for scarless assembly (`split`, `pad`).
- 🧪 **Degenerate mode:** compress variant libraries with low mutational diversity into IUPAC-degenerate oligos for cost-efficient synthesis (`compress`, `expand`).
- 📈 **Analysis mode:** fast NGS-based activity quantification with read indexing, packing, and barcode/associate counting (`index`, `pack`, `acount`, `xcount`) extensible with callback methods (via Python library).
- 🔁 **Iterative & multiplexed workflows:** `patch_mode` for extending existing pools, cross-set barcode separation, and per-group primer design with cross-compatibility screening.
- ⚡ **Performance:** scalable to very large libraries and high-throughput sequencing datasets, with published benchmarks demonstrating efficient design and analysis on commodity hardware (see paper).
- 🔒 **Rich constraints:** IUPAC sequence constraints, motif exclusion, repeat screening, Hamming-distance barcodes, and primer thermodynamic constraints (including optional paired-primer Tm matching).
- 📊 **DataFrame-centric:** modules operate on CSV/DataFrames and return updated tables plus `stats`; the CLI can emit JSON and supports reproducible stochastic runs (`random_seed`).
- 💻 **CLI + library-first:** full-featured command-line interface with YAML config files, multi-step pipelines (sequential or parallel DAG), **and** a composable Python API for interactive use in scripts and Jupyter notebooks.
- 🤖 **AI-assisted design:** agent-ready documentation for Claude, ChatGPT, and Copilot.


<a id="ai-assisted-design"></a>
## 🤖 AI-Assisted Design

`Oligopool Calculator` is optimized for AI-assisted workflows. Either share the [`docs/agent-skills.md`](https://github.com/ayaanhossain/oligopool/blob/master/docs/agent-skills.md) file with your agent, or share the following raw URL along with a suitable prompt, for direct parsing.
```
https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/agent-skills.md
```
Ensure that your AI/agent explores this document thoroughly. Afterwards, you can chat about the package, your specific design goals, and have the agent plan and execute the design and analysis pipelines.


<a id="installation"></a>
## 📦 Installation

`Oligopool Calculator` is a `Python 3.10+`-exclusive library.

On `Linux`, `macOS`, and `Windows Subsystem for Linux`, you can install `Oligopool Calculator` from [PyPI](https://pypi.org/project/oligopool/), where it is published as the `oligopool` package.
```bash
$ pip install --upgrade oligopool # Installs and/or upgrades oligopool
```
This also installs the command line tools: `oligopool` and `op`.

Or install it directly from GitHub:
```bash
$ pip install git+https://github.com/ayaanhossain/oligopool.git
```
Both approaches should install all dependencies automatically.
> **Note** The GitHub version will always be updated with all recent fixes. The PyPI version should be more stable.

If you are on `Windows` or simply prefer to, `Oligopool Calculator` can also be used via `Docker` (please see [the notes](https://github.com/ayaanhossain/oligopool/blob/master/docs/docker-notes.md)).

Successful installation will look like this.
```python
$ python
>>> import oligopool as op
>>> op.__version__
'2026.02.01'
>>>
```

<a id="getting-started"></a>
## 🚀 Getting Started

`Oligopool Calculator` is carefully designed, easy to use, and stupid fast.

You can import the library and use its various functions either in a script or interactively inside a `Jupyter` environment. Use `help(...)` to read the docs as necessary and follow along.

The [`examples`](https://github.com/ayaanhossain/oligopool/tree/master/examples) directory includes a [design parser](https://github.com/ayaanhossain/oligopool/tree/master/examples/design-assembly-parser), a [library compressor](https://github.com/ayaanhossain/oligopool/tree/master/examples/library-compressor), an [analysis pipeline](https://github.com/ayaanhossain/oligopool/tree/master/examples/analysis-pipeline), and a complete [CLI YAML pipeline](https://github.com/ayaanhossain/oligopool/tree/master/examples/cli-pipeline).

If you want the full end-to-end walkthrough, start with the notebook: [`Oligopool Calculator` in action](https://github.com/ayaanhossain/oligopool/blob/master/examples/OligopoolCalculatorInAction.ipynb).

**Documentation:**
- [User Guide](https://github.com/ayaanhossain/oligopool/blob/master/docs/docs.md) - Comprehensive tutorials, examples, and workflows
- [API Reference](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md) - Complete parameter documentation for all modules
- [AI Agent Guide](https://github.com/ayaanhossain/oligopool/blob/master/docs/agent-skills.md) - Decision trees, best practices, and gotchas for AI-assisted design (Claude, ChatGPT, Copilot)
- [Docker Guide](https://github.com/ayaanhossain/oligopool/blob/master/docs/docker-notes.md) - Run `oligopool` in a container for cross-platform consistency

```python
$ python
>>>
>>> import oligopool as op
>>> help(op)
...
    Automated design and analysis of oligo pool libraries for
    high-throughput genomics and synthetic biology applications.

    Design Mode - build synthesis-ready oligo architectures
        barcode     orthogonal barcodes with Hamming distance guarantees
        primer      Tm-optimized primers with off-target screening
        motif       sequence motifs or anchors
        spacer      neutral fill to reach target length
        background  k-mer database for off-target screening
        merge       collapse columns into single element
        revcomp     reverse complement a column range
        lenstat     length statistics and free-space check
        verify      verify length, motif, and background conflicts
        final       concatenate into synthesis-ready oligos

    Assembly Mode - fragment long oligos for assembly
        split       fragment oligos into overlapping pieces
        pad         Type IIS primer pads for scarless excision

    Degenerate Mode - compress variant libraries for synthesis
        compress    reduce similar variants to IUPAC-degenerate oligos
        expand      expand IUPAC-degenerate oligos into concrete sequences

    Analysis Mode - quantify variants from NGS reads
        index       index barcodes and associated variants
        pack        filter/merge/deduplicate FastQ reads
        acount      association counting (barcode + variant verification)
        xcount      combinatorial counting (single or multiple barcodes)

    Advanced
        vectorDB    LevelDB k-mer storage
        Scry        1-NN barcode classifier

    Usage
        >>> import oligopool as op
        >>> df, stats = op.barcode(input_data='variants.csv', ...)
        >>> help(op.barcode)  # module docs

    Modules return (DataFrame, stats). Chain them iteratively; use patch_mode=True
    to extend pools without overwriting existing designs.

    CLI: `op` | `op COMMAND` | Docs: https://github.com/ayaanhossain/oligopool
...
```

<a id="command-line-interface-cli"></a>
## 💻 Command Line Interface (CLI)

The `oligopool` package installs a CLI with two equivalent entry points: `oligopool` and `op`.

```bash
$ op
$ op cite
$ op manual
$ op manual topics
$ oligopool manual barcode
```

Run `op` with no arguments to see the command list, and run `op COMMAND` to see command-specific options.
```bash
$ op

oligopool v2026.02.01
by ah

Oligopool Calculator is a suite of algorithms for
automated design and analysis of oligo pool libraries.

usage: oligopool COMMAND --argument=<value> ...

COMMANDS Available:

    manual      show module documentation
    cite        show citation information

    pipeline    execute multi-step pipeline from config

    barcode     orthogonal barcodes with cross-set separation
    primer      thermodynamic primers with optional Tm matching
    motif       design or add motifs/anchors
    spacer      neutral spacers to meet length targets

    background  build k-mer background database

    split       break long oligos into overlapping fragments
    pad         add excisable primer pads for scarless assembly

    merge       collapse contiguous columns
    revcomp     reverse-complement a column range

    lenstat     compute length stats and free space
    verify      verify length, motif, and background conflicts
    final       finalize into synthesis-ready oligos

    compress    compress sequences into IUPAC-degenerate oligos
    expand      expand IUPAC oligos to concrete sequences

    index       build barcode/associate index
    pack        preprocess and deduplicate FastQ reads
    acount      association counting (single index)
    xcount      combinatorial counting (multiple indexes)

    complete    print or install shell completion

Run "oligopool COMMAND" to see command-specific options.
```

Install tab-completion to blaze through interactive CLI use (recommended).
```bash
$ op complete --install          # auto-detect shell (restart your shell)
$ op complete --install bash     # or: zsh|fish
```

### YAML Pipelines

Define entire workflows in a single YAML config file and execute with one command:
```bash
$ op pipeline --config mpra_design.yaml
$ op pipeline --config mpra_design.yaml --dry-run  # validate first
```

Pipelines support **sequential** or **parallel DAG** execution-independent steps run concurrently:
```yaml
# mpra_design.yaml
pipeline:
  name: "MPRA Library"
  steps:
    - name: fwd_primer
      command: primer
    - name: rev_primer
      command: primer          # parallel with fwd_primer
    - name: add_barcode
      command: barcode
      after: [fwd_primer]      # waits for fwd_primer
    - name: finalize
      command: final
      after: [add_barcode, rev_primer]

fwd_primer:
  input_data: "variants.csv"
  output_file: "fwd.csv"
  primer_type: forward
  # ...
```

Use `--config` with any command to load parameters from YAML (CLI args override config values). See [`examples/cli-pipeline`](https://github.com/ayaanhossain/oligopool/tree/master/examples/cli-pipeline) for a complete working example.

> **CLI Notes**
> - Commands that write files require an output basename (e.g., `--output-file`, `--index-file`, `--pack-file`, `--count-file`, `--mapping-file`, `--synthesis-file`), unlike library mode where outputs can be returned in-memory.
> - Most `--*-type` parameters accept either integers or descriptive strings (case-insensitive), e.g. `--primer-type forward`, `--barcode-type spectrum`, `--motif-type anchor`, `--pack-type merge`, `--mapping-type sensitive`.
> - For `--primer-sequence-constraint` / `--motif-sequence-constraint`, pass an IUPAC string (`NNNN...`) or a quoted expression like `"'N'*20"` / `'GCC+N*20+CCG'`.
> - `op split` writes separate files per fragment by default (e.g., `out.Split1.oligopool.split.csv`, `out.Split2...`); use `--no-separate-outputs` to write a single combined file.

<a id="citation"></a>
## 📖 Citation

If you use `Oligopool Calculator` in your research publication, please cite our paper.

```
Hossain A, Cetnar DP, LaFleur TL, McLellan JR, Salis HM.
Automated Design of Oligopools and Rapid Analysis of Massively Parallel Barcoded Measurements.
ACS Synth Biol. 2024;13(12):4218-4232. doi:10.1021/acssynbio.4c00661
```

BibTeX:
```bibtex
@article{Hossain2024Oligopool,
  title   = {Automated Design of Oligopools and Rapid Analysis of Massively Parallel Barcoded Measurements},
  author  = {Hossain, Ayaan and Cetnar, Daniel P. and LaFleur, Travis L. and McLellan, James R. and Salis, Howard M.},
  journal = {ACS Synthetic Biology},
  year    = {2024},
  volume  = {13},
  number  = {12},
  pages   = {4218--4232},
  doi     = {10.1021/acssynbio.4c00661}
}
```

You can read the paper online for free at [ACS Synthetic Biology](https://doi.org/10.1021/acssynbio.4c00661).
* PMCID: `PMC11669329`
* PMID: `39641628`

<a id="license"></a>
## ⚖️ License

`Oligopool Calculator` (c) 2026 Ayaan Hossain.

`Oligopool Calculator` is an **open-source software** under the [GPL-3.0](https://opensource.org/license/gpl-3-0) license.

See [LICENSE](https://github.com/ayaanhossain/oligopool/blob/master/LICENSE) file for more details.
