'''
Lightweight command and module descriptions.

This module exists to keep the package manual (`help(oligopool)`) and the CLI menu
(`op` / `oligopool`) semantically aligned without reusing the exact same phrasing.

It must stay dependency-free (no numpy/pandas) because it is imported during CLI
startup and `import oligopool`.
'''

# Package-facing phrasing (library/manual).
PACKAGE_DESCRIPTIONS = {
    # Design Mode
    'barcode': 'orthogonal barcodes with Hamming distance guarantees',
    'primer': 'Tm-optimized primers with off-target screening',
    'motif': 'sequence motifs or anchors',
    'spacer': 'neutral fill to reach target length',
    'background': 'k-mer database for off-target screening',
    'merge': 'collapse columns into single element',
    'join': 'join two tables on ID with ordered insertion',
    'revcomp': 'reverse complement a column range',
    'final': 'concatenate into synthesis-ready oligos',
    # Assembly Mode
    'split': 'fragment oligos into overlapping pieces',
    'pad': 'Type IIS primer pads for scarless excision',
    # Degenerate Mode
    'compress': 'reduce similar variants to IUPAC-degenerate oligos',
    'expand': 'expand IUPAC-degenerate oligos into concrete sequences',
    # Analysis Mode
    'index': 'index barcodes and associated variants',
    'pack': 'filter/merge/deduplicate FastQ reads',
    'acount': 'association counting (barcode + variant verification)',
    'xcount': 'combinatorial counting (single or multiple barcodes)',
    # QC Mode
    'lenstat': 'length statistics and free-space check',
    'verify': 'verify length, motif, and background conflicts',
    'inspect': 'inspect background/index/pack artifacts',
}


# CLI-facing phrasing (short command help).
CLI_DESCRIPTIONS = {
    # Design Mode
    'barcode': 'orthogonal barcodes with cross-set separation',
    'primer': 'thermodynamic primers with optional Tm matching',
    'motif': 'design or add motifs/anchors',
    'spacer': 'neutral spacers to meet length targets',
    'background': 'build k-mer background database',
    'merge': 'collapse contiguous columns',
    'join': 'join two tables on ID',
    'revcomp': 'reverse-complement a column range',
    'final': 'finalize into synthesis-ready oligos',
    # Assembly Mode
    'split': 'break long oligos into overlapping fragments',
    'pad': 'add excisable primer pads for scarless assembly',
    # Degenerate Mode
    'compress': 'compress sequences into IUPAC-degenerate oligos',
    'expand': 'expand IUPAC oligos to concrete sequences',
    # Analysis Mode
    'index': 'build barcode/associate index',
    'pack': 'preprocess and deduplicate FastQ reads',
    'acount': 'association counting (single index)',
    'xcount': 'combinatorial counting (multiple indexes)',
    # QC Mode
    'lenstat': 'compute length stats and free space',
    'verify': 'detect length, motif, and background conflicts',
    'inspect': 'inspect non-CSV artifacts',
    # Other
    'manual': 'show module documentation',
    'cite': 'show citation information',
    'pipeline': 'execute multi-step pipeline from config',
    'complete': 'print or install shell completion',
}
