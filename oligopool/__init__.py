# Core Design Functions
from .background import background
from .barcode import barcode
from .primer import primer
from .motif import motif
from .spacer import spacer

# Auxiliary Design Functions
from .lenstat import lenstat
from .final import final

# Assembly Design Functions
from .split import split
from .pad import pad

# Analysis Functions
from .index import index
from .pack import pack
from .acount import acount
from .xcount import xcount

# Other goodies
from .base.vectordb import vectorDB
from .base.scry import Scry


# Setup
__author__ = 'Ayaan Hossain'

__version__ = '2024.09.29'

__doc__ = f'''
oligopool v{__version__}
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
        b. All other columns in the DataFrame must be DNA sequences
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

    Auxiliary modules available
        - lenstat
        - final

    Assembly modules available
        - split
        - pad

    Worfklow example sketch

    >>> import pandas as pd
    >>> import oligopool as op
    >>>
    >>> # Read initial library
    >>> init_df = pd.read_csv('initial_library.csv')
    >>>
    >>> # Add oligo elements one by one
    >>> barcode_df, stats = op.barcode(input_data=init_df, ...)
    >>> primer_df,  stats = op.primer(input_data=barcode_df, ...)
    ...
    >>> # Check length statistics when needed
    >>> op.lenstat(input_data=further_along_df)
    ...
    >>> # Finalize the library
    >>> final_df, stats = op.final(input_data=ready_to_go_df, ...)
    >>>
    >>> # Split and pad longer oligos if needed
    >>> split_df, stats = op.split(input_data=final_df, ...)
    >>> first_pad_df,  stats = op.pad(input_data=split_df, ...)
    >>> second_pad_df, stats = op.pad(input_data=split_df, ...)
    ...

Analysis Mode workflow

    1. Index one or more CSVs containing the barcode information
    2. Pack all NGS FastQ files, optionally merging them if required
    3. Use acount for association counting of variants and barcodes
    4. If multiple barcode combinations are to be counted use xcount
    5. Combine count DataFrames and perform stats and ML as necessary

    Indexing module available
        - index

    Packing module available
        - pack

    Counting modules available
        - acount
        - xcount

    Analysis Mode worfklow example sketch

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
    >>>
    >>> # Pack the FastQ files
    >>> sam1_pack_stats = op.pack(r1_file='sample_1_R1.fq.gz', ...)
    >>> sam2_pack_stats = op.pack(r1_file='sample_2_R1.fq.gz', ...)
    ...
    >>>
    >>> # Compute and write the barcode combination count matrix
    >>> xcount_stats = op.xcount(index_files=['bc1_index', 'bc2_index'],
    ...                          pack_file='sample_1_pack', ...)
    ...
    >>>
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

For advanced use cases, the following classes are also available.
    - vectorDB
    - Scry
'''