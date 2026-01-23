<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/" style="text-decoration: none !important;">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/main/oligopool/img/logo.svg"  alt="Oligopool Calculator" width="460" class="center"/>
    </a>
</h1>

<h4><p align="center">Version: 2026.01.22</p></h4>

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

`Oligopool Calculator` is a suite of algorithms for automated design and analysis of [oligopool libraries](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9300125/).

It supports scalable design of primers, barcodes, motifs/anchors, and spacers; assembly-aware splitting/padding of long constructs; and rapid packing/counting of barcoded reads for activity quantification.

We have used `Oligopool Calculator` in multiple projects to build libraries of tens of thousands of promoters (see [here](https://www.nature.com/articles/s41467-022-32829-5) and [here](https://www.nature.com/articles/s41587-020-0584-2)), ribozymes, and mRNA stability elements (see [here](https://www.nature.com/articles/s41467-024-54059-7)).

To learn more, please check out [our paper in ACS Synthetic Biology](https://pubs.acs.org/doi/10.1021/acssynbio.4c00661).

`Oligopool Calculator` streamlines the design and analysis of massively parallel barcoded measurements, including iterative workflows where new oligos are continuously added to an existing pool. It has been benchmarked to design millions of compact barcodes and universal primer binding sites, and to process hundreds of millions of sequencing reads per hour on an 8-core desktop.

<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/" style="text-decoration: none !important;">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/refs/heads/main/oligopool/img/workflow.svg"  alt="Oligopool Calculator Workflow" width="3840" class="center"/>
    </a>
</h1>

**Design and analysis of oligopool variants using `Oligopool Calculator`.** **(a)** In `Design Mode`, `Oligopool Calculator` can be used to generate optimized `barcode`s, `primer`s, `spacer`s, `motif`s and `split` longer oligos into shorter `pad`ded fragments for downstream synthesis and assembly. **(b)** Once the library is assembled and cloned, barcoded amplicon sequencing data can be processed via `Analysis Mode` for characterization. `Analysis Mode` proceeds by first `index`ing one or more sets of barcodes, `pack`ing the reads, and then producing count matrices either using `acount` (association counting) or `xcount` (combinatorial counting).

<a id="features"></a>
## ✨ Features

- 🧬 **Design mode:** constraint-based design of primers, barcodes, motifs/anchors, and spacers, with background screening, assembly helpers (`split`, `pad`), and utilities (`merge`, `revcomp`, `lenstat`, `verify`, `final`).
- 🔁 **Iterative & multiplexed workflows:** patch mode for extending existing pools, cross-set barcode separation, and per-group primer design with cross-compatibility screening.
- 📈 **Analysis mode:** fast activity quantification with read indexing, packing, and barcode/associate counting (`index`, `pack`, `acount`, `xcount`) extensible with callback methods (via Python library).
- ⚡ **Performance:** scalable to very large libraries and high-throughput sequencing datasets, with published benchmarks demonstrating efficient design and analysis on commodity hardware (see paper).
- 🔒 **Rich constraints:** IUPAC sequence constraints, motif exclusion, repeat screening, Hamming-distance barcodes, and primer thermodynamic constraints (including optional paired-primer Tm matching).
- 📊 **DataFrame-centric:** modules operate on CSV/DataFrames and return updated tables plus `stats`; the CLI can emit JSON and supports reproducible stochastic runs (`random_seed`).
- 💻 **CLI + library-first:** full-featured command-line pipelines **and** a composable Python API for interactive use in scripts and Jupyter notebooks.


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
Python 3.10.9 | packaged by conda-forge | (main, Feb  2 2023, 20:20:04) [GCC 11.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import oligopool as op
>>> op.__version__
'2026.01.22'
>>>
```

<a id="getting-started"></a>
## 🚀 Getting Started

`Oligopool Calculator` is carefully designed, easy to use, and stupid fast.

You can import the library and use its various functions either in a script or interactively inside a `Jupyter` environment. Use `help(...)` to read the docs as necessary and follow along.

There are examples of a [design parser](https://github.com/ayaanhossain/oligopool/blob/master/examples/design-parser/design_parser.py) and an [analysis pipeline](https://github.com/ayaanhossain/oligopool/blob/master/examples/analysis-pipeline/analysis_pipeline.py) inside the [`examples`](https://github.com/ayaanhossain/oligopool/tree/master/examples) directory.

A notebook demonstrating [`Oligopool Calculator` in action](https://github.com/ayaanhossain/oligopool/blob/master/examples/OligopoolCalculatorInAction.ipynb) is provided there as well.

**Documentation:**
- [User Guide](https://github.com/ayaanhossain/oligopool/blob/master/docs/docs.md) - Comprehensive tutorials, examples, and workflows
- [API Reference](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md) - Complete parameter documentation for all modules
- [AI Agent Guide](https://github.com/ayaanhossain/oligopool/blob/master/docs/agent-skills.md) - Decision trees, best practices, and gotchas for AI-assisted design (Claude, ChatGPT, Copilot)
- [Docker Guide](https://github.com/ayaanhossain/oligopool/blob/master/docs/docker-notes.md) - Run oligopool in a container for cross-platform consistency

```python
$ python
Python 3.12.6 | packaged by conda-forge | (main, Sep 30 2024, 18:08:52) [GCC 13.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>>
>>> import oligopool as op
>>> help(op)
...
    oligopool v2026.01.22
    by ah

    Automated design and analysis of oligopool libraries.

    Modules operate on CSV/DataFrames and return (output_df, stats) or stats.
    Chain them to build libraries iteratively; use Patch Mode (`patch_mode=True`)
    to extend existing pools without overwriting prior designs.

    Design Mode
        barcode     orthogonal barcodes with cross-set separation
        primer      thermodynamic primers with Tm matching
        motif       sequence motifs or constant anchors
        spacer      neutral spacers to meet length targets
        background  k-mer database for off-target screening
        split       break long oligos into overlapping fragments
        pad         add primers + Type IIS sites for assembly
        merge       collapse contiguous columns
        revcomp     reverse complement a column range
        lenstat     length statistics and free-space check
        verify      QC constraints before synthesis
        final       concatenate columns into synthesis-ready oligos

    Analysis Mode
        index       build barcode/associate index
        pack        preprocess and deduplicate FastQ reads
        acount      association counting (single index)
        xcount      combinatorial counting (multiple indexes)

    Advanced
        vectorDB    LevelDB-based k-mer storage
        Scry        1-NN barcode classifier

    Usage
        >>> import oligopool as op
        >>> help(op.barcode)           # module docs
        >>> df, stats = op.barcode(input_data=df, ...)

    CLI: `oligopool` / `op` (run `op` for commands, `op manual barcode` for docs)
    Docs: https://github.com/ayaanhossain/oligopool
...
```

<a id="command-line-interface-cli"></a>
### 💻 Command Line Interface (CLI)

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

oligopool v2026.01.22
by ah

usage: oligopool COMMAND --argument=<value> ...

COMMANDS Available:

    manual      show module documentation
    cite        show citation information

    barcode     orthogonal barcodes with cross-set separation
    primer      thermodynamic primers with Tm matching
    motif       sequence motifs or constant anchors
    spacer      neutral spacers to meet length targets

    background  k-mer database for off-target screening

    split       break long oligos into overlapping fragments
    pad         add primers + Type IIS sites for assembly

    merge       collapse contiguous columns
    revcomp     reverse complement a column range

    lenstat     length statistics and free-space check
    verify      QC constraints before synthesis
    final       concatenate columns into synthesis-ready oligos

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

> **CLI Notes**
> - Commands that write a DataFrame require `--output-file` (unlike in library mode where it is optional).
> - For `--primer-sequence-constraint` / `--motif-sequence-constraint`, pass an IUPAC string (`NNNN...`) or a quoted expression like `"'N'*20"` / `'GCC+N*20+CCG'`.
> - `op split` writes separate files per fragment by default (e.g., `out.Split1.oligopool.split.csv`, `out.Split2...`); use `--no-separate-outputs` for a single combined file.

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
PMCID: `PMC11669329` • PMID: `39641628`

<a id="license"></a>
## ⚖️ License

`Oligopool Calculator` (c) 2026 Ayaan Hossain.

`Oligopool Calculator` is an **open-source software** under the [GPL-3.0](https://opensource.org/license/gpl-3-0) license.

See [LICENSE](https://github.com/ayaanhossain/oligopool/blob/master/LICENSE) file for more details.
