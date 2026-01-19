<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/" style="text-decoration: none !important;">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/main/oligopool/img/logo.svg"  alt="Oligopool Calculator" width="460" class="center"/>
    </a>
</h1>

<h4><p align="center">Version: 2026.01.18</p></h4>

<p align="center">
  <a href="#features" style="text-decoration: none !important;">✨ Features</a> •
  <a href="#installation" style="text-decoration: none !important;">📦 Installation</a> •
  <a href="#getting-started" style="text-decoration: none !important;">🚀 Getting Started</a> •
  <a href="#command-line-interface-cli" style="text-decoration: none !important;">💻 CLI</a> •
  <a href="#citation" style="text-decoration: none !important;">📖 Citation</a> •
  <a href="#license" style="text-decoration: none !important;">⚖️ License</a>
</p>

`Oligopool Calculator` is a suite of algorithms for automated design and analysis of [oligopool libraries](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9300125/).

It supports scalable design of primers, barcodes, motifs/anchors, and spacers; assembly-aware splitting/padding of long constructs; and rapid packing/counting of barcoded reads for activity quantification.

We have used `Oligopool Calculator` in multiple projects to build libraries of tens of thousands of promoters (see [here](https://www.nature.com/articles/s41467-022-32829-5) and [here](https://www.nature.com/articles/s41587-020-0584-2)), ribozymes, and mRNA stability elements (see [here](https://www.nature.com/articles/s41467-024-54059-7)). To learn more, please check out [our paper in ACS Synthetic Biology](https://pubs.acs.org/doi/10.1021/acssynbio.4c00661).

`Oligopool Calculator` streamlines the design and analysis of massively parallel barcoded measurements (including MPRAs), including iterative workflows where new oligos are appended to an existing pool. It has been benchmarked to design millions of compact barcodes and to process hundreds of millions of sequencing reads per hour on a desktop.

<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/" style="text-decoration: none !important;">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/refs/heads/main/oligopool/img/workflow.svg"  alt="Oligopool Calculator Workflow" width="3840" class="center"/>
    </a>
</h1>

**Design and analysis of oligopool variants using `Oligopool Calculator`.** **(a)** In `Design Mode`, `Oligopool Calculator` can be used to generate optimized `barcode`s, `primer`s, `spacer`s, `motif`s and `split` longer oligos into shorter `pad`ded fragments for downstream synthesis and assembly. **(b)** Once the library is assembled and cloned, barcoded amplicon sequencing data can be processed via `Analysis Mode` for characterization. `Analysis Mode` proceeds by first `index`ing one or more sets of barcodes, `pack`ing the reads, and then producing count matrices either using `acount` (association counting) or `xcount` (combinatorial counting).

<a id="features"></a>
## ✨ Features

- 🧬 **Design Mode:** constraint-based design of primers, barcodes, motifs/anchors, and spacers, with optional repeat screening via `background`, plus assembly helpers (`split`, `pad`) and utilities (`merge`, `revcomp`, `lenstat`, `verify`, `final`).
- 🔒 **Rich constraints:** IUPAC sequence constraints, motif exclusion, repeat screening, Hamming-distance barcodes, and primer thermodynamic constraints (including optional paired-primer Tm matching).
- 🧠 **Algorithmic core:** orthogonally symmetric barcode design, adaptive decision trees for primer design, a Scry barcode classifier, and efficient read packing.
- ⚡ **Performance (benchmarked):** design >4 million compact barcodes in ~1.2 h, design universal primer binding sites for 1 million 200-mer oligos in ~15 min, and analyze ~500 million reads/hour on an 8-core desktop (see paper).
- 🧾 **DataFrame + stats:** most modules take a CSV/DataFrame and return either an updated DataFrame + `stats` or a `stats` dictionary; the CLI can emit JSON (`--stats-json` / `--stats-file`) and stochastic modules accept `random_seed`.
- 🔁 **Patch Mode (iterative pool extension):** `patch_mode=True` / `--patch-mode` fills only missing values in an existing output column without overwriting prior designs.
- 🧷 **Cross-set barcodes:** enforce separation from existing barcode columns via `cross_barcode_columns` + `minimum_cross_distance`.
- 🧫 **Multiplexed primer sets:** design primers per group via `oligo_sets` and screen for cross-set compatibility (supports chained primer design via `paired_primer_column`, including per-set pairing).
- ✅ **QC + summaries:** `lenstat` (length/free-space) and `verify` (constraint checks and motif emergence) before ordering/synthesis.
- 🧪 **Analysis Mode (counting readouts):** `index`, `pack`, `acount`, `xcount` for fast packing and barcode/associate counting.
- 💻 **CLI + completion:** `op` / `oligopool` for pipelines, built-in `manual`/`cite`, and tab completion via `op complete --install`.


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

If you are on `Windows` or simply prefer to, `Oligopool Calculator` can also be used via `Docker` (please see [the notes](https://github.com/ayaanhossain/oligopool/blob/master/docker-notes.md)).

**Verifying Installation**

Successful installation will look like this.
```python
$ python
Python 3.10.9 | packaged by conda-forge | (main, Feb  2 2023, 20:20:04) [GCC 11.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import oligopool as op
>>> op.__version__
'2026.01.18'
>>>
```

<a id="getting-started"></a>
## 🚀 Getting Started

`Oligopool Calculator` is carefully designed, easy to use, and stupid fast.

You can import the library and use its various functions either in a script or interactively inside a `Jupyter` environment. Use `help(...)` to read the docs as necessary and follow along.

There are examples of a [design parser](https://github.com/ayaanhossain/oligopool/blob/master/examples/design-parser/design_parser.py) and an [analysis pipeline](https://github.com/ayaanhossain/oligopool/blob/master/examples/analysis-pipeline/analysis_pipeline.py) inside the [`examples`](https://github.com/ayaanhossain/oligopool/tree/master/examples) directory.

A notebook demonstrating [`Oligopool Calculator` in action](https://github.com/ayaanhossain/oligopool/blob/master/examples/OligopoolCalculatorInAction.ipynb) is provided there as well. It shows in-depth use of all major design and analysis methods.

For quick iteration during Design Mode, use `lenstat` to monitor length/free space under an `oligo_length_limit`, and `verify` as a final QC pass before ordering/synthesis.

```python
$ python
Python 3.12.6 | packaged by conda-forge | (main, Sep 30 2024, 18:08:52) [GCC 13.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>>
>>> import oligopool as op
>>> help(op)
...
    oligopool v2026.01.18
    by ah

    Automated design and analysis of oligopool libraries.

    The various modules in Oligopool Calculator can be used
    interactively in a jupyter notebook, or be used to define
    scripts for design and analysis pipelines on the cloud.

    Oligopool Calculator is designed for composable and iterative
    workflows: modules can be chained to add new elements as new
    columns, and element modules support Patch Mode (`patch_mode=True`)
    to fill only missing values when extending an existing pool
    without overwriting prior designs.

    Oligopool Calculator offers two modes of operation
        -   Design Mode for designing oligopool libraries, and
        - Analysis Mode for analyzing oligopool datasets.

    Design Mode workflow

        1. Initialize a pandas DataFrame with core library elements.
            a. The DataFrame must contain a unique 'ID' column serving as primary key.
            b. All other columns in the DataFrame must be DNA sequences.
        2. Define any optional background sequences via the background module.
        3. Add necessary oligopool elements with constraints via element modules.
            a. For multiplexing, design per-set primers via `oligo_sets` in `primer`.
            b. For multiple barcode sets, enforce cross-set separation via
               `cross_barcode_columns` + `minimum_cross_distance` in `barcode`.
            c. For incremental pool extension, use Patch Mode (`patch_mode=True`) in element modules
               to fill only missing values in an existing output column.
        4. Optionally, split long oligos and pad them via assembly modules.
        5. Perform additional maneuvers and finalize library via auxiliary modules.

        Background module available
            - background

        Element modules available
            - primer
            - barcode
            - motif
            - spacer

        Assembly modules available
            - split
            - pad

        Auxiliary modules available
            - merge
            - revcomp
            - lenstat
            - verify
            - final

        Design Mode example sketch

            >>> import pandas as pd
            >>> import oligopool as op
            >>>
            >>> # Read initial library
            >>> init_df = pd.read_csv('initial_library.csv')
            >>>
            >>> # Add oligo elements one by one
            >>> primer_df,  stats = op.primer(input_data=init_df, ...)
            >>> barcode_df, stats = op.barcode(input_data=primer_df, ...)
            ...
            >>> # Check length statistics as needed
            >>> length_stats = op.lenstat(input_data=further_along_df)
            ...
            >>>
            >>> # Split and pad longer oligos if needed
            >>> split_df, stats = op.split(input_data=even_further_along_df, ...)
            >>> first_pad_df,  stats = op.pad(input_data=split_df, ...)
            >>> second_pad_df, stats = op.pad(input_data=split_df, ...)
            ...
            >>>
            >>> # Finalize the library
            >>> final_df, stats = op.final(input_data=ready_to_go_df, ...)
            ...

    Analysis Mode workflow

        1. Index one or more CSVs containing barcode (and associate) data.
        2. Pack all NGS FastQ files, optionally merging them if required.
        3. Use acount for association counting of variants and barcodes.
        4. If multiple barcode combinations are to be counted use xcount.
        5. Combine count DataFrames and perform stats and ML as necessary.

        Indexing module available
            - index

        Packing module available
            - pack

        Counting modules available
            - acount
            - xcount

        Analysis Mode example sketch

            >>> import pandas as pd
            >>> import oligopool as op
            >>>
            >>> # Read annotated library
            >>> bc1_df = pd.read_csv('barcode_1.csv')
            >>> bc2_df = pd.read_csv('barcode_2.csv')
            >>> av1_df = pd.read_csv('associate_1.csv')
            ...
            >>>
            >>> # Index barcodes and any associates
            >>> bc1_index_stats = op.index(barcode_data=bc1_df, ...)
            >>> bc2_index_stats = op.index(barcode_data=bc2_df, ...)
            ...
            >>>
            >>> # Pack experiment FastQ files
            >>> sam1_pack_stats = op.pack(r1_fastq_file='sample_1_R1.fq.gz', ...)
            >>> sam2_pack_stats = op.pack(r1_fastq_file='sample_2_R1.fq.gz', ...)
            ...
            >>>
            >>> # Compute and write barcode combination count matrix
            >>> xcount_df, stats = op.xcount(index_files=['bc1_index', 'bc2_index'], ...)
            ...

    You can learn more about each module using help.
        >>> import oligopool as op
        >>> help(op)
        >>> help(op.primer)
        >>> help(op.barcode)
        ...
        >>> help(op.xcount)

    For advanced uses, the following classes are also available.
        - vectorDB
        - Scry
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

oligopool v2026.01.18
by ah

Oligopool Calculator is a suite of algorithms for
automated design and analysis of oligopool libraries.

usage: oligopool COMMAND --argument=<value> ...

COMMANDS Available:

    manual      show module documentation
    cite        show citation information

    background  build background k-mer database
    barcode     design constrained barcodes
    primer      design constrained primers
    motif       design or add motifs and anchors
    spacer      design or insert spacers

    split       split long oligos into shorter ones
    pad         pad split oligos with primers

    merge       merge oligo elements into one column
    revcomp     reverse complement spanning elements

    lenstat     compute length statistics
    verify      verify constraints and summarize library
    final       finalize library

    index       index barcodes and associates
    pack        preprocess and pack FastQ for counting
    acount      execute association counting
    xcount      execute combinatorial counting

    complete    print or install shell completion

Run "oligopool COMMAND" to see command-specific options.

EST <timestamp>
```

Install tab-completion to blaze through interactive CLI use (recommended).
```bash
$ op complete --install          # auto-detect shell (restart your shell)
$ op complete --install bash     # or: zsh|fish
```

> **CLI Notes**
> - Commands that write a DataFrame require `--output-file` (library mode can return DataFrames in-memory).
> - For `--primer-sequence-constraint` / `--motif-sequence-constraint`, pass an IUPAC string (`NNNN...`) or a quoted expression like `"'N'*20"` / `'GCC+N*20+CCG'`.
> - Run `op COMMAND` to see all options (including `--patch-mode`, `--oligo-sets`, and `--cross-barcode-columns`).

<a id="citation"></a>
## 📖 Citation

If you use `Oligopool Calculator` in your research publication, please cite our paper.

You can also run `op cite` / `oligopool cite`.

```
Hossain A, Cetnar DP, LaFleur TL, McLellan JR, Salis HM.
Automated Design of Oligopools and Rapid Analysis of Massively Parallel Barcoded Measurements.
ACS Synth Biol. 2024;13(12):4218-4232. doi:10.1021/acssynbio.4c00661
```

BibTeX:
```bibtex
@article{hossain2024oligopool,
  title   = {Automated Design of Oligopools and Rapid Analysis of Massively Parallel Barcoded Measurements},
  author  = {Hossain, Ayaan and Cetnar, David P. and LaFleur, Tyler L. and McLellan, James R. and Salis, Howard M.},
  journal = {ACS Synthetic Biology},
  year    = {2024},
  volume  = {13},
  number  = {12},
  pages   = {4218--4232},
  doi     = {10.1021/acssynbio.4c00661},
}
```

You can read the complete article online for free at [ACS Synthetic Biology](https://doi.org/10.1021/acssynbio.4c00661).
PMCID: `PMC11669329` • PMID: `39641628`

<a id="license"></a>
## ⚖️ License

`Oligopool Calculator` (c) 2026 Ayaan Hossain.

`Oligopool Calculator` is an **open-source software** under the [GPL-3.0](https://opensource.org/license/gpl-3-0) license.

See [LICENSE](https://github.com/ayaanhossain/oligopool/blob/master/LICENSE) file for more details.
