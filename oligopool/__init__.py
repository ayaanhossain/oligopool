import importlib
from typing import Any

from .descriptions import PACKAGE_DESCRIPTIONS as _PKG_DESC

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

__version__ = '2026.01.25'

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


def _build_package_doc() -> str:
    return f'''
Automated design and analysis of oligopool libraries for
high-throughput genomics and synthetic biology applications.

Design Mode - build synthesis-ready oligo architectures
    barcode     {_PKG_DESC["barcode"]}
    primer      {_PKG_DESC["primer"]}
    motif       {_PKG_DESC["motif"]}
    spacer      {_PKG_DESC["spacer"]}
    background  {_PKG_DESC["background"]}
    split       {_PKG_DESC["split"]}
    pad         {_PKG_DESC["pad"]}
    merge       {_PKG_DESC["merge"]}
    revcomp     {_PKG_DESC["revcomp"]}
    lenstat     {_PKG_DESC["lenstat"]}
    verify      {_PKG_DESC["verify"]}
    final       {_PKG_DESC["final"]}

Analysis Mode - quantify variants from NGS reads
    index       {_PKG_DESC["index"]}
    pack        {_PKG_DESC["pack"]}
    acount      {_PKG_DESC["acount"]}
    xcount      {_PKG_DESC["xcount"]}

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
'''


__doc__ = _build_package_doc().strip()
