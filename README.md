<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/main/oligopool/img/logo.svg"  alt="Oligopool Calculator" width="460" class="center"/>
    </a>
</h1>

<p align="center">
  <a href="#Overview">Overview</a> •
  <a href="#Installation">Installation</a> •
  <a href="#Getting-Started">Getting Started</a> •
  <a href="#License">License</a> •
  <a href="#Citation">Citation</a> •
</p>


## Overview

`Oligopool Calculator` is an end-to-end suite of algorithms for automated design and analysis of oligopool libraries. It enables the scalable design of orthogonal primer sets and unique barcodes, the splitting of long constructs into multiple oligos, and the rapid packing and counting of barcoded reads -- all on a regular 8-core desktop computer. We have demonstrated `Oligopool Calculator` across several projects to build oligopool libraries spanning thousands of promoters, ribozymes, and mRNA stability elements, illustrating the use of a flexible grammar to add multiple barcodes, cut sites, avoid excluded sequences, optimize constraints, and characterized those libraries using highly efficient barcoded read counting functions.

`Oligopool Calculator` facilitates the creative design and application of massively parallel reporter assays by automating solutions to its computational complexities.

<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/refs/heads/main/oligopool/img/workflow.svg"  alt="Oligopool Calculator Workflow" width="3840" class="center"/>
    </a>
</h1>

**Design and analysis of oligopool variants using `Oligopool Calculator`.** **(a)** In `Design Mode`, `Oligopool Calculator` can be used to design optimized `barcode`s, `primer`s, `spacer`s, `motif`s and `split` large oligos into shorter `pad`ded fragments for downstream synthesis and assembly. **(b)** Once the library is assembled and cloned, barcoded amplicon sequencing data can be processed via `Analysis Mode` to generate count matrices for characterizing the variants. `Analysis Mode` proceeds by first `index`ing one or more sets of barcodes, `pack`ing the reads, and then producing the count matrices either using `acount` (association counting) or `xcount` (combinatorial counting).


## Installation

`Oligopool Calculator` is a `Linux`-based `Python3.10+`-exclusive library.

You can install it from [PyPI](https://pypi.org/project/oligopool/), where it is published as the `oligopool` package
```bash
$ pip install oligopool
```
or install it directly from GitHub.
```bash
$ pip install git+https://github.com/ayaanhossain/oligopool.git
```
Both approaches should install all dependencies automatically.

> **Note** This GitHub version will always be the most updated version with all recent fixes. The PyPI version is aimed to be the most stable one.


**Verifying Installation**

Successful installation will look like this.
```python
$ python
Python 3.10.9 | packaged by conda-forge | (main, Feb  2 2023, 20:20:04) [GCC 11.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import oligopool as op
>>> op.__version__
'2024.09.29'
>>>
```

**Uninstalling `Oligopool Calculator`**

You can easily remove the package via
```bash
$ pip uninstall oligopool
```

## Getting Started

`Oligopool Calculator` is very simple and easy to use. Import the library and start using the various functions either in a script or in a `jupyter` `notebook` environment. Just use `help(...)` to read the docs as necessary and follow along!

There are example scripts of a [design parser](https://github.com/ayaanhossain/oligopool/blob/master/examples/design-parser/design_parser.py) and an [analysis pipleine](https://github.com/ayaanhossain/oligopool/blob/master/examples/analysis-pipeline/analysis_pipeline.py) inside the [examples](https://github.com/ayaanhossain/oligopool/tree/master/examples) directory with this repository. A `jupyter` `notebook` demonstrating [`Oligopool Calculator` in action](https://github.com/ayaanhossain/oligopool/blob/master/examples/OligopoolCalculatorInAction.ipynb) is provided there as well.

```python
$ python
Python 3.12.6 | packaged by conda-forge | (main, Sep 30 2024, 18:08:52) [GCC 13.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>>
>>> import oligopool as op
>>> help(op)

...
    oligopool v2024.09.29
    by ah

    Automated design and analysis of oligopool libraries.

    The various modules in Oligopool Calculator can be used
    interactively in a jupyter notebook, or be used to define
    scripts for design and analysis pipelines on the cloud.

    Oligopool Calculator offers two modes of operation
        -   Design Mode for designing oligopool libraries
        - Analysis Mode for analyzing oligopool datasets

    Design Mode workflow

        1. Initialize a pandas DataFrame with core library elements
            a. The DataFrame must contain a unique 'ID' column serving as primary key
            b. All other columns in the DataFrame must be DNA sequences or a - (dash)
        2. Next, define any optional background sequences via the background module
        3. Add necessary oligopool elements with constraints via element modules
        4. Optionally, split long oligos and pad them via assembly modules
        5. Perform length checks as needed and finalize library via auxiliary modules

        Background module available
            - background

        Element modules available
            - primer
            - barcode
            - spacer
            - motif

        Assembly modules available
            - split
            - pad

        Auxiliary modules available
            - lenstat
            - final

        Worfklow example sketch

        >>> import pandas as pd
        >>> import oligopool as op
        >>>
        >>> # Read initial library
        >>> df = pd.read_csv('initial_library.csv')
        >>>
        >>> # Add oligo elements
        >>> barcode_df, stats = op.barcode(input_data=df, ...)
        >>> primer_df,  stats = op.primer(input_data=barcode_df, ...)
        ...
        >>> # Check length statistics
        >>> op.lenstat(input_data=primer_df)
        >>>
        >>> # Finalize the library
        >>> final_df, stats = op.final(input_data=primer_df, ...)
        >>>
        >>> # Split and pad longer oligos if needed
        >>> split_df, stats = op.split(input_data=final_df, ...)
        >>> first_pad_df,  stats = op.pad(input_data=split_df, ...)
        >>> second_pad_df, stats = op.pad(input_data=split_df, ...)
        ...

    Analysis Mode workflow

        1. Index one or more CSVs containing the barcode and anchor information
        2. Pack all NGS FastQ files, optionally merging them if required
        3. If a barcode and its association with the core variant is to be counted use acount
        4. If multiple barcode combinations are to be counted use xcount
        5. Combine count matrices and perform stats and ML as necessary

        Indexing module available
            - index

        Packing module available
            - pack

        Counting modules available
            - acount
            - xcount

        Anamysis Mode worfklow example sketch

        >>> import pandas as pd
        >>> import oligopool as op
        >>>
        >>> # Read marked up library
        >>> df = pd.read_csv('oligo_library.csv')
        >>>
        >>> # Index the barcodes and save the indexes
        >>> bc1_index_stats = op.index(input_data=df, barcode_column='BC1', ...)
        >>> bc2_index_stats = op.index(input_data=df, barcode_column='BC2', ...)
        ...
        >>> # Pack the FastQ files
        >>> sam1_pack_stats = op.pack(r1_file='sample_1_R1.fq.gz', ...)
        >>> sam2_pack_stats = op.pack(r1_file='sample_2_R1.fq.gz', ...)
        ...
        >>> # Compute and write the barcode combination count matrix
        >>> xcount_stats = op.xcount(index_files=['bc1_index', 'bc2_index'],
        ...                          pack_file='sample_1_pack', ...)
        ...
        >>> # Read the count matrix and continue analysis
        >>> count_df = pd.read_csv('sample_1_counts.csv')
        ...

    You can learn more about each module using help.
    >>> import oligopool as op
    >>>
    >>> help(op)
    >>> help(op.primer)
    >>> help(op.barcode)
    ...
    >>> help(op.xcount)
...

```

## License

`Oligpool Calculator` (c) 2024 Ayaan Hossain.

`Oligpool Calculator` is an **open-source software** under [GPLv3](https://opensource.org/license/gpl-3-0) License.

See [LICENSE](https://github.com/ayaanhossain/oligopool/blob/master/LICENSE) file for more details.