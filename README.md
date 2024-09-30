<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/main/oligopool/img/logo.svg"  alt="Oligopool Calculator" width="460" class="center"/>
    </a>
</h1>

<p align="center">
  <a href="#Overview">Overview</a> •
  <a href="#Installation">Installation</a> •
  <a href="#License">License</a> •
  <a href="#Citation">Citation</a> •
  <a href="https://github.com/ayaanhossain/oligopool/blob/master/docs/DOCS.md">Documentation</a>
</p>

## Overview

Massively parallel reporter assays enable the simultaneous construction and characterization of large libraries of designed genetic parts and systems, combining oligopool synthesis with next-generation sequencing. As library sizes grow, it becomes computationally challenging to optimally design oligopool sequences and to rapidly analyze sequencing results.

We present the `Oligopool Calculator`, an end-to-end suite of algorithms, that automates the design and analysis of massively parallel reporter assays with millions of defined variants. New algorithms and data structures enable the scalable design of orthogonal primer sets and unique barcodes, the splitting of long constructs into multiple oligos, and the rapid packing and counting of reads -- all on a regular 8-core desktop computer. We demonstrate the Oligopool Calculator across several projects to build and characterize thousands of promoters, ribozymes, and mRNA stability elements, illustrating how to use a flexible grammar to add multiple barcodes, cut sites, excluded sequences, model-based constraints, read filters, and read processing functions.

The `Oligopool Calculator` facilitates the creative design and application of massively parallel reporter assays by automating solutions to its computational complexities.

<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/refs/heads/main/oligopool/img/workflow.svg"  alt="Oligopool Calculator Workflow" width="3840" class="center"/>
    </a>
</h1>

**Design and analysis of oligopool variants using Oligopool Calculator.** **(a)** In `Design Mode`, `Oligopool Calculator` can be used to design optimized `barcode`s, `primer`s, `spacer`s, `motif`s and `split` large oligos into shorter `pad`ded fragments for downstream synthesis and assembly. **(b)** Once the genetic part variants are assembled and cloned inside a chassis or activated in vitro, barcoded amplicon libraries are sequenced as readout, and the resulting data is processed via `Analysis Mode` to generate count matrices of barcode abundance which serve as a proxy for genetic part activity. `Analysis Mode` proceeds by first `index`ing one or more sets of barcodes, `pack`ing reads to pre-count and collapse the readouts, and then producing the count matrices either using `acount` (association counting) or `xcount` (combinatorial counting).

## Installation

`Oligopool Calculator` is `Linux`-based software, and it is a `Python3.6+`-exclusive library.


**Approach 1: Install from PyPI**

You can install `Oligopool Calculator` from PyPI, where it is published as the `oligopool` package. This is as easy as
```bash
$ pip install --upgrade oligopool --no-cache-dir
```
which will also install all dependencies from PyPI.


**Approach 2: Install from GitHub**

Alternatively, you can install `Oligopool Calculator` from GitHub. To do so, first clone the repository with
```bash
$ git clone https://github.com/ayaanhossain/oligopool.git
```
Once downloaded, navigate into the `oligopool` directory with
```bash
$ cd oligopool
```
and install using `pip`
```bash
$ pip install .
```
which will also install all dependencies. The GitHub version will always be the most updated version.


**Verifying Installation**

If everything went well, `Oligopool Calculator` is now available in your environment under the `oligopool` package. You may verify it like so:
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

You can easily remove the package with
```bash
$ pip uninstall oligopool
```