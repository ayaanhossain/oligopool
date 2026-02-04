import time as tt

import atexit as ae

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_index as ci


def index(
    barcode_data:str|pd.DataFrame,
    barcode_column:str,
    index_file:str,
    barcode_prefix_column:str|None=None,
    barcode_suffix_column:str|None=None,
    barcode_prefix_gap:int=0,
    barcode_suffix_gap:int=0,
    associate_data:str|pd.DataFrame|None=None,
    associate_column:str|None=None,
    associate_prefix_column:str|None=None,
    associate_suffix_column:str|None=None,
    associate_prefix_gap:int=0,
    associate_suffix_gap:int=0,
    verbose:bool=True) -> dict:
    '''
    Build an index object for mapping barcodes (and optional associates) in NGS reads.
    The resulting `.oligopool.index` is used by `acount`/`xcount` for fast counting via `Scry`.

    Required Parameters:
        - `barcode_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated barcodes.
        - `barcode_column` (`str`): Column name for the barcodes to index.
        - `index_file` (`str`): Filename for output index object.

    Optional Parameters:
        - `barcode_prefix_column` (`str` / `None`): Column for the constant barcode prefix anchor (default: `None`).
        - `barcode_suffix_column` (`str` / `None`): Column for the constant barcode suffix anchor (default: `None`).
        - `barcode_prefix_gap` (`int`): Distance between the prefix and barcode (default: 0).
        - `barcode_suffix_gap` (`int`): Distance between the suffix and barcode (default: 0).
        - `associate_data` (`str` / `pd.DataFrame` / `None`): Path to a CSV file or DataFrame with annotated variants
            (default: `None`).
        - `associate_column` (`str` / `None`): Column name for the associate elements to index (default: `None`).
        - `associate_prefix_column` (`str` / `None`): Column for the constant associate prefix (default: `None`).
        - `associate_suffix_column` (`str` / `None`): Column for the constant associate suffix (default: `None`).
        - `associate_prefix_gap` (`int`): Distance between the prefix and associate (default: 0).
        - `associate_suffix_gap` (`int`): Distance between the suffix and associate (default: 0).
        - `verbose` (`bool`): If `True`, logs progress to stdout (default: `True`).

    Returns:
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `barcode_data` (and `associate_data`, if provided) must contain a unique 'ID' column; all other
            columns used for indexing must be non-empty DNA strings.
        - Provide at least one of `barcode_prefix_column` or `barcode_suffix_column` (and likewise for the
            associate, if used).
        - For error-free indexing, anchor columns should be adjacent to the corresponding barcode/associate column.
        - In reads, anchors must be exactly `barcode_{prefix|suffix}_gap` bases away from the barcode (same for
            associate anchors).
        - For association counting, partial presence of the associate variant can be sufficient; however, its
            prefix/suffix anchors must be adjacent and fully present.
        - Multiple indices may be used simultaneously for combinatorial counting (associate info is ignored).
        - Anchor sequences can be designed with `motif(motif_type=1, ...)` and then referenced via
            `{barcode|associate}_{prefix|suffix}_column`.
    '''

    # Alias Arguments
    barcodedata     = barcode_data
    barcodecol      = barcode_column
    indexfile       = index_file
    barcodeprefix   = barcode_prefix_column
    barcodesuffix   = barcode_suffix_column
    barcodepregap   = barcode_prefix_gap
    barcodepostgap  = barcode_suffix_gap
    associatedata    = associate_data
    associatecol     = associate_column
    associateprefix  = associate_prefix_column
    associatesuffix  = associate_suffix_column
    associatepregap  = associate_prefix_gap
    associatepostgap = associate_suffix_gap
    verbose          = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Indexing Verbiage Print
    liner.send('\n[Oligopool Calculator: Analysis Mode - Index]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass barcodedata Parsing and Validation
    (bardf,
    barcodedata_valid) = vp.get_parsed_indata_info(
        indata=barcodedata,
        indata_field='   Barcode Data   ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)
    input_rows = len(bardf.index) if isinstance(bardf, pd.DataFrame) else 0

    # Sort bardf
    if barcodedata_valid:
        bardf.sort_index(
            inplace=True)

    # Full barcodecol Validation
    (barcodes,
    barcodecol_valid) = vp.get_parsed_column_info(
        col=barcodecol,
        df=bardf,
        col_field='   Barcode Column ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=None,
        adjval=None,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full indexfile Validation
    indexfile_valid = vp.get_outfile_validity(
        outfile=indexfile,
        outfile_suffix='.oligopool.index',
        outfile_field='     Index File   ',
        liner=liner)

    # Adjust indexfile Suffix
    indexfile = ut.get_adjusted_path(
        path=indexfile,
        suffix='.oligopool.index')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full barcodeprefix Validation
    (barcodeprefix,
    barcodeprefix_valid) = vp.get_parsed_column_info(
        col=barcodeprefix,
        df=bardf,
        col_field='   Barcode Prefix ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=barcodecol,
        adjval=+1,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full barcodesuffix Validation
    (barcodesuffix,
    barcodesuffix_valid) = vp.get_parsed_column_info(
        col=barcodesuffix,
        df=bardf,
        col_field='   Barcode Suffix ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=barcodecol,
        adjval=-1,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full barcodepregap Validation
    barcodepregap_valid = vp.get_numeric_validity(
        numeric=barcodepregap,
        numeric_field='   Barcode Pregap ',
        numeric_pre_desc=' Allow Exactly ',
        numeric_post_desc=' bp Gap b/w Prefix and Barcode',
        minval=0,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full barcodepostgap Validation
    barcodepostgap_valid = vp.get_numeric_validity(
        numeric=barcodepostgap,
        numeric_field='   Barcode Postgap',
        numeric_pre_desc=' Allow Exactly ',
        numeric_post_desc=' bp Gap b/w Barcode and Suffix',
        minval=0,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # First Pass associatedata Parsing and Validation
    (assdf,
    associatedata_valid) = vp.get_parsed_associatedata_info(
        associatedata=associatedata,
        associatedata_field=' Associate Data   ',
        required_fields=('ID',),
        bardf=bardf,
        barcodedata_valid=barcodedata_valid,
        liner=liner)

    # Full associatecol Validation
    (associates,
    associatecol_valid) = vp.get_parsed_column_info(
        col=associatecol,
        df=assdf,
        col_field=' Associate Column ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=None,
        adjval=None,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full associateprefix Validation
    (associateprefix,
    associateprefix_valid) = vp.get_parsed_column_info(
        col=associateprefix,
        df=assdf,
        col_field=' Associate Prefix ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=associatecol,
        adjval=+1,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full associatesuffix Validation
    (associatesuffix,
    associatesuffix_valid) = vp.get_parsed_column_info(
        col=associatesuffix,
        df=assdf,
        col_field=' Associate Suffix ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=associatecol,
        adjval=-1,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full associatepregap Validation
    associatepregap_valid = vp.get_numeric_validity(
        numeric=associatepregap,
        numeric_field=' Associate Pregap ',
        numeric_pre_desc=' Allow Exactly ',
        numeric_post_desc=' bp Gap b/w Prefix and Associate',
        minval=0,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full associatepostgap Validation
    associatepostgap_valid = vp.get_numeric_validity(
        numeric=associatepostgap,
        numeric_field=' Associate Postgap',
        numeric_pre_desc=' Allow Exactly ',
        numeric_post_desc=' bp Gap b/w Associate and Suffix',
        minval=0,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # First Pass Validation
    if not all([
        barcodedata_valid,
        barcodecol_valid,
        indexfile_valid,
        barcodeprefix_valid,
        barcodesuffix_valid,
        barcodepregap_valid,
        barcodepostgap_valid,
        associatedata_valid,
        associatecol_valid,
        associateprefix_valid,
        associatesuffix_valid,
        associatepregap_valid,
        associatepostgap_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Setup Warning Dictionary
    warns = {}

    # Parse Barcode Data Feasibility
    liner.send('\n[Step 1: Parsing Barcode Data]\n')

    # Parse bardf
    (parsestatus,
    IDdict,
    barcodedict,
    duplicates,
    barcodecount,
    barcodesuniq,
    barcodelen,
    barcodelenuniq,
    minbarcodelen,
    maxbarcodelen) = ci.get_parsed_barcodes(
        bardf=bardf,
        barcodes=barcodes,
        liner=liner)

    # Barcode Data infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 1,
            'step_name': 'parsing-barcode-data',
            'vars'    : {
                        'duplicates': duplicates,
                     'barcode_count': barcodecount,
                  ' barcodes_unique': barcodesuniq,
                   'min_barcode_len': minbarcodelen,
                   'max_barcode_len': maxbarcodelen,
                'barcode_len_unique': barcodelenuniq},
            'warns'   : warns}
        stats = ut.stamp_stats(
            stats=stats,
            module='index',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return stats

    # Parse Barcode Constant Feasibility
    liner.send('\n[Step 2: Parsing Barcode Constants]\n')

    # Parse barcodeprefix and barcodesuffix
    (parsestatus,
    constantsextracted,
    constantsuniq,
    longconstants,
    barcodeprefix,
    barcodesuffix,
    prefixuniq,
    suffixuniq,
    prefixlen,
    suffixlen) = ci.get_parsed_constants(
        prefixconstant=barcodeprefix,
        suffixconstant=barcodesuffix,
        attachetype=0,
        liner=liner)

    # Barcode Constants infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 2,
            'step_name': 'parsing-barcode-constants',
            'vars'    : {
                'constants_extracted': constantsextracted,
                   'constants_unique': constantsuniq,
                     'long_constants': longconstants,
                     'prefix_unique': prefixuniq,
                     'suffix_unique': suffixuniq,
                        'prefix_len': prefixlen,
                        'suffix_len': suffixlen},
            'warns'   : warns}
        stats = ut.stamp_stats(
            stats=stats,
            module='index',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return stats

    # Do we have associated variants?
    if not associatedata is None:

        # Extract Associate Data
        liner.send('\n[Step 3: Extracting Associate Data]\n')

        # Update Step 3 Warning
        warns[3] = {
            'warn_count': 0,
            'step_name' : 'parsing-associate-data',
            'vars': None}

        # Extract associatedata
        associatedict = ci.get_extracted_associates(
            associates=associates,
            expcount=barcodecount,
            IDdict=IDdict,
            warn=warns[3],
            liner=liner)

        # Parse Associate Constant Feasibility
        liner.send('\n[Step 4: Parsing Associate Constants]\n')

        # Parse associateprefix and associatesuffix
        (parsestatus,
        constantsextracted,
        constantsuniq,
        longconstants,
        associateprefix,
        associatesuffix,
        prefixuniq,
        suffixuniq,
        prefixlen,
        suffixlen) = ci.get_parsed_constants(
            prefixconstant=associateprefix,
            suffixconstant=associatesuffix,
            attachetype=1,
            liner=liner)

        # Associate Constants infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 4,
                'step_name': 'parsing-associate-constants',
                'vars'    : {
                    'constants_extracted': constantsextracted,
                       'constants_unique': constantsuniq,
                         'long_constants': longconstants,
                          'prefix_unique': prefixuniq,
                          'suffix_unique': suffixuniq,
                             'prefix_len': prefixlen,
                             'suffix_len': suffixlen},
                'warns'   : warns}
            stats = ut.stamp_stats(
                stats=stats,
                module='index',
                input_rows=input_rows,
                output_rows=0)

            # Return results
            liner.close()
            return stats

    else:
        associatedict   = {}
        associateprefix = None
        associatesuffix = None

    # Setup Workspace
    (indexfile,
    indexdir) = ut.setup_workspace(
        outfile=indexfile,
        outfile_suffix='.oligopool.index')

    # Schedule indexfile deletion
    ifdeletion = ae.register(
        ut.remove_file,
        indexfile)

    # Prepared Objects Queue
    indexqueue = ut.SafeQueue()

    # Launching Indexing
    liner.send('\n[Step 5: Computing Index]\n')

    # Define Index Stats
    stats = {
        'status'  : True,
        'basis'   : 'solved',
        'step'    : 5,
        'step_name': 'computing-index',
        'vars'    : None,
        'warns'   : warns}

    # Compute Index Objects
    ci.index_engine(
        IDdict=IDdict,
        barcodedict=barcodedict,
        barcodename=barcodecol,
        barcodecount=barcodecount,
        barcodelen=barcodelen,
        barcodeprefix=barcodeprefix,
        barcodesuffix=barcodesuffix,
        barcodepregap=barcodepregap,
        barcodepostgap=barcodepostgap,
        associatedict=associatedict,
        associateprefix=associateprefix,
        associatesuffix=associatesuffix,
        associatepregap=associatepregap,
        associatepostgap=associatepostgap,
        indexdir=indexdir,
        indexqueue=indexqueue,
        liner=liner)

    # Archive All Indexed Objects
    ut.archive(
        objqueue=indexqueue,
        arcfile=indexfile,
        mode='x',
        prodcount=1,
        prodactive=0,
        liner=liner)

    # Indexing Status?
    indexstatus = 'Successful'

    # Indexing Stats
    liner.send('\n[Indexing Statistics]\n')

    liner.send(
        ' Indexing Status : {}\n'.format(
            indexstatus))
    liner.send(
        '     Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Remove Workspace
    ut.remove_directory(
        dirpath=indexdir)

    # Unschedule packfile deletion
    if indexstatus == 'Successful':
        ae.unregister(ifdeletion)

    # Close Liner
    liner.close()

    stats = ut.stamp_stats(
        stats=stats,
        module='index',
        input_rows=input_rows,
        output_rows=0)
    return stats
