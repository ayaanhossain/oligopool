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

Oligopool synthesis and next-generation sequencing enable the construction and characterization of large libraries of designed genetic parts and systems. As library sizes grow, it becomes computationally challenging to optimally design large numbers of primer binding sites, barcode sequences, and overlap regions to obtain efficient assemblies and precise measurements.

We present the `Oligopool Calculator`, an *end-to-end suite of algorithms and data structures*, that rapidly designs many thousands of oligonucleotides within an oligopool and rapidly analyzes many billions of barcoded sequencing reads. We introduce several novel concepts that greatly increase the design and analysis throughput, including orthogonally symmetric `barcode` design, adaptive decision trees for `primer` design, a `Scry` barcode classifier, and efficient read packing.

We demonstrated the `Oligopool Calculator`’s capabilities across computational benchmarks and real-data projects, including the design of over four million highly unique and compact barcodes in 1.2 hours, the design of universal primer binding sites for million 200-mer oligos in 15 minutes, and the analysis of about 500 million sequencing reads per hour, all on an 8-core desktop computer.

<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/refs/heads/main/oligopool/img/workflow.svg"  alt="Oligopool Calculator Workflow" width="1080" class="center"/>
    </a>
</h1>

**Design and analysis of oligopool variants using Oligopool Calculator.** **(a)** In `Design Mode`, `Oligopool Calculator` can be used to design optimized `barcode`s, `primer`s, `spacer`s, `motif`s and `split` large oligos into shorter `pad`ded fragments for downstream synthesis and assembly. **(b)** Once the genetic part variants are assembled and cloned inside a chassis or activated in vitro, barcoded amplicon libraries are sequenced as readout, and the resulting data is processed via `Analysis Mode` to generate count matrices of barcode abundance which serve as a proxy for genetic part activity. `Analysis Mode` proceeds by first `index`ing one or more sets of barcodes, `pack`ing reads to pre-count and collapse the readouts, and then producing the count matrices either using `acount` (association counting) or `xcount` (combinatorial counting).