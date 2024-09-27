<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/main/oligopool/img/logo.svg"  alt="Non-Repetitive Parts Calculator" width="460" class="center"/>
    </a>
</h1>

<p align="center">
  <a href="#Overview">Overview</a> •
  <a href="#Installation">Installation</a> •
  <a href="#License">License</a> •
  <a href="#Citation">Citation</a> •
  <a href="https://github.com/ayaanhossain/oligopool/blob/master/docs/DOCS.md">DOCS</a> •
  <a href="#Oligopool-Calculator-in-Action">Oligopool Calculator in Action!</a>
</p>

## Overview

Oligopool synthesis and next-generation sequencing enable the construction and characterization of large libraries of designed genetic parts and systems. As library sizes grow, it becomes computationally challenging to optimally design large numbers of primer binding sites, barcode sequences, and overlap regions to obtain efficient assemblies and precise measurements.

We present the `Oligopool Calculator`, an end-to-end suite of algorithms and data structures, that rapidly designs many thousands of oligonucleotides within an oligopool and rapidly analyzes many billions of barcoded sequencing reads. We introduce several novel concepts that greatly increase the design and analysis throughput, including orthogonally symmetric `barcode` design, adaptive decision trees for `primer` design, a `Scry` barcode classifier, and efficient read packing.

We demonstrate the `Oligopool Calculator`’s capabilities across computational benchmarks and real-data projects, including the design of over four million highly unique and compact barcodes in 1.2 hours, the design of universal primer binding sites for million 200-mer oligos in 15 minutes, and the analysis of about 500 million sequencing reads per hour, all on an 8-core desktop computer.