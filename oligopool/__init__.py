import importlib
from typing import Any

__api__ = [
    'background',
    'barcode',
    'primer',
    'motif',
    'spacer',
    'merge',
    'revcomp',
    'lenstat',
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

__version__ = '2026.01.17'

_LAZY_ATTRS = {
    # Core design functions
    'background': ('.background', 'background'),
    'barcode': ('.barcode', 'barcode'),
    'primer': ('.primer', 'primer'),
    'motif': ('.motif', 'motif'),
    'spacer': ('.spacer', 'spacer'),
    # Auxiliary design functions
    'merge': ('.merge', 'merge'),
    'revcomp': ('.revcomp', 'revcomp'),
    'lenstat': ('.lenstat', 'lenstat'),
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

The various modules in Oligopool Calculator can be used
interactively in a jupyter notebook, or be used to define
scripts for design and analysis pipelines on the cloud.

Oligopool Calculator offers two modes of operation
    -   Design Mode for designing oligopool libraries, and
    - Analysis Mode for analyzing oligopool datasets.

Design Mode workflow

    1. Initialize a pandas DataFrame with core library elements.
        a. The DataFrame must contain a unique 'ID' column serving as primary key.
        b. All other columns in the DataFrame must be DNA sequences.
    2. Define any optional background sequences via the background module.
    3. Add necessary oligopool elements with constraints via element modules.
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
'''
