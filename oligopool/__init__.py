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

__version__ = '2026.01.22'

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

__doc__ = f'''
oligopool v{__version__}
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
'''
