import importlib
from typing import Any

__api__ = [
    'barcode',
    'background',
    'primer',
    'motif',
    'spacer',
    'merge',
    'revcomp',
    'lenstat',
    'verify',
    'final',
    'split',
    'pad',
    'index',
    'pack',
    'acount',
    'xcount',
    'vectorDB',
    'Scry',
    '__author__',
    '__version__',
]

# Keep `help(oligopool)` focused on the package-level manual text.
# Users can run `help(oligopool.barcode)` (etc.) for module details.
__all__ = ['__author__', '__version__']

# Setup
__author__ = 'Ayaan Hossain'

__version__ = '2026.01.24'

_LAZY_ATTRS = {
    # Core design functions
    'barcode': ('.barcode', 'barcode'),
    'background': ('.background', 'background'),
    'primer': ('.primer', 'primer'),
    'motif': ('.motif', 'motif'),
    'spacer': ('.spacer', 'spacer'),
    # Auxiliary design functions
    'merge': ('.merge', 'merge'),
    'revcomp': ('.revcomp', 'revcomp'),
    'lenstat': ('.lenstat', 'lenstat'),
    'verify': ('.verify', 'verify'),
    'final': ('.final', 'final'),
    # Assembly design functions
    'split': ('.split', 'split'),
    'pad': ('.pad', 'pad'),
    # Analysis functions
    'index': ('.index', 'index'),
    'pack': ('.pack', 'pack'),
    'acount': ('.acount', 'acount'),
    'xcount': ('.xcount', 'xcount'),
    # Other goodies
    'vectorDB': ('.base.vectordb', 'vectorDB'),
    'Scry': ('.base.scry', 'Scry'),
}


def __getattr__(name: str) -> Any:
    target = _LAZY_ATTRS.get(name)
    if not target:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = target
    module = importlib.import_module(module_name, package=__name__)
    value = getattr(module, attr_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    # Keep `dir(oligopool)` and `help(oligopool)` fast and dependency-light by
    # not advertising lazily-imported attributes (which `pydoc` would then try
    # to resolve via `getattr`, importing heavier scientific dependencies).
    return sorted(globals())

__doc__ = '''
Automated design and analysis of oligopool libraries for high-throughput
functional genomics (MPRAs, CRISPR screens, saturation mutagenesis, etc.).

Design Mode - build synthesis-ready oligo architectures
    barcode     orthogonal barcodes with Hamming distance guarantees
    primer      Tm-optimized primers with off-target screening
    motif       sequence motifs or constant anchors
    spacer      neutral fill to reach target length
    background  k-mer database for off-target screening
    split       fragment long oligos for assembly
    pad         Type IIS primer pads for scarless excision
    merge       collapse columns into single element
    revcomp     reverse complement column range
    lenstat     length statistics and free-space check
    verify      QC before synthesis
    final       concatenate into synthesis-ready oligos

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

CLI: `op --help` | Docs: https://github.com/ayaanhossain/oligopool
'''
