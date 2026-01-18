import time as tt

import collections as cx

import textwrap as tw

import numpy as np
import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_lenstat as cl


def verify(
    input_data:str|pd.DataFrame,
    oligo_length_limit:int|None=None,
    excluded_motifs:list|str|pd.DataFrame|None=None,
    verbose:bool=True) -> dict:
    '''
    Verify a library DataFrame by running a lightweight QC pass (architecture summary, length checks,
    and motif checks) without modifying or writing outputs.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with an 'ID' column.

    Optional Parameters:
        - `oligo_length_limit` (`int` / `None`): If provided, checks for length overflow (default: `None`).
        - `excluded_motifs` (`list` / `str` / `pd.DataFrame`): Motifs to track; can be a CSV path or DataFrame
            with 'ID' and 'Exmotif' columns (default: `None`).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A dictionary of verification results (stats only).

    Notes:
        - `verify` is a summary/QC module and does not return a DataFrame.
        - Like `lenstat`, `verify` computes length statistics from sequence columns (ignoring non-sequence
          annotations) and can optionally check against an `oligo_length_limit`.
        - Metadata columns (non-sequence annotations) are tracked and excluded from sequence-only checks.
        - Degenerate/IUPAC sequence columns (e.g., UMIs with 'N') are flagged but not treated as hard errors.
          Most design modules assume concrete DNA (`A`/`T`/`G`/`C`/`-`), so this is intended as a guardrail.
        - Columns with a mix of DNA-like and non-DNA values are flagged as "mixed" (common source of subtle bugs).
        - Excluded-motif check flags motif "emergence" in assembled oligos: for each motif, it reports any
          variants where the motif occurs more often than the minimum occurrence across the library (useful
          for motifs that should appear exactly once, like restriction sites).
    '''

    # Argument Aliasing
    indata     = input_data
    oligolimit = oligo_length_limit
    exmotifs   = excluded_motifs
    verbose    = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Verify Verbiage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Verify]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass input_data Parsing and Validation (no DNA-only constraint)
    (indf,
    data_name,
    indata_valid) = vp.get_parsed_data_info(
        data=indata,
        data_field='    Input Data   ',
        required_fields=None,
        liner=liner)
    input_rows = len(indf.index) if isinstance(indf, pd.DataFrame) else 0

    # Show update on successful parse
    if indata_valid:
        liner.send(
            '    Input Data  : {} w/ {:,} Record(s)\n'.format(
                data_name,
                len(indf.index)))

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full oligolimit Validation (optional)
    oligolimit_valid = True
    if oligolimit is None:
        liner.send('    Oligo Limit : None Specified\n')
    else:
        oligolimit_valid = vp.get_numeric_validity(
            numeric=oligolimit,
            numeric_field='    Oligo Limit ',
            numeric_pre_desc=' At most ',
            numeric_post_desc=' Base Pair(s)',
            minval=4,
            maxval=float('inf'),
            precheck=False,
            liner=liner)

    # Full exmotifs Parsing and Validation (optional)
    (exmotifs,
    exmotifs_valid) = vp.get_parsed_exseqs_info(
        exseqs=exmotifs,
        exseqs_field=' Excluded Motifs',
        exseqs_desc='Unique Motif(s)',
        df_field='Exmotif',
        required=False,
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid,
        exmotifs_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    if not oligolimit is None:
        oligolimit = round(oligolimit)

    # Book-keeping
    stats = None
    warns = {}

    # Step 1: Column Inspection
    liner.send('\n[Step 1: Inspecting Columns]\n')

    def send_wrapped(prefix, text, width=79):
        for line in tw.wrap(
            text,
            width=width,
            initial_indent=prefix,
            subsequent_indent=' ' * len(prefix)):
            liner.send(line + '\n')

    columns = list(indf.columns)
    duplicates = [c for c, n in cx.Counter(columns).items() if n > 1]

    dna_cols = []
    ddna_cols = []
    seq_cols = []
    seq_positions = []

    metadata_cols = []
    metadata_col_types = {}
    non_string_cols = []

    mixed_cols = []
    mixed_col_examples = {}

    non_string_examples = {}
    empty_seq_cols = []
    empty_seq_examples = {}

    for col_idx, col in enumerate(columns):
        series = indf.iloc[:, col_idx]

        # Non-string columns are treated as metadata (annotations).
        is_str = series.map(lambda x: isinstance(x, str))
        if not bool(is_str.all()):
            metadata_cols.append(col)
            non_string_cols.append(col)
            metadata_col_types[f'{col}[{col_idx}]'] = str(series.dtype)
            invalid_mask = ~is_str
            examples = ut.get_row_examples(
                df=indf,
                invalid_mask=invalid_mask,
                id_col='ID',
                limit=5)
            try:
                sample = series[invalid_mask].iloc[0]
                sample_type = type(sample).__name__
                sample_value = repr(sample)
            except Exception:
                sample_type = 'unknown'
                sample_value = '<unavailable>'
            non_string_examples[f'{col}[{col_idx}]'] = {
                'example_type': sample_type,
                'example_value': sample_value,
                'examples': examples,
            }
            continue

        upper = series.str.upper()

        # Determine if the column behaves like a sequence column (DNA/IUPAC) or metadata.
        ddna_mask = upper.map(lambda x: ut.is_DNA(seq=x, dna_alpha=ut.ddna_alpha))
        ddna_fraction = float(ddna_mask.mean()) if len(ddna_mask) else 0.0

        # No DNA-like values: treat as metadata string column.
        if ddna_fraction == 0.0:
            metadata_cols.append(col)
            metadata_col_types[f'{col}[{col_idx}]'] = str(series.dtype)
            continue

        # Mixed DNA-like and non-DNA values: flag (likely a broken sequence column).
        if ddna_fraction < 1.0:
            mixed_cols.append(col)
            invalid_mask = ~ddna_mask
            examples = ut.get_row_examples(
                df=indf,
                invalid_mask=invalid_mask,
                id_col='ID',
                limit=5)
            try:
                sample = upper[invalid_mask].iloc[0]
            except Exception:
                sample = '<unavailable>'
            mixed_col_examples[f'{col}[{col_idx}]'] = {
                'ddna_fraction': ddna_fraction,
                'example_value': str(sample),
                'examples': examples,
            }
            continue

        # Sequence column: all values are DNA/IUPAC (including '-').
        empty_mask = upper.str.len().eq(0)
        if bool(empty_mask.any()):
            empty_seq_cols.append(col)
            empty_seq_examples[f'{col}[{col_idx}]'] = ut.get_row_examples(
                df=indf,
                invalid_mask=empty_mask,
                id_col='ID',
                limit=5)

        seq_cols.append(col)
        seq_positions.append(col_idx)
        indf.iloc[:, col_idx] = upper

        is_dna = upper.map(lambda x: ut.is_DNA(seq=x))
        if bool(is_dna.all()):
            dna_cols.append(col)
        else:
            ddna_cols.append(col)

    # De-duplicate metadata columns (keep order) while preserving duplicates list separately.
    metadata_cols = list(cx.OrderedDict.fromkeys(metadata_cols))

    field_width = 20

    liner.send(
        ' {:>{}}: {:,} Column(s)\n'.format(
            'Sequence Columns',
            field_width,
            len(seq_cols)))
    liner.send(
        ' {:>{}}: {:,} Column(s)\n'.format(
            'DNA Columns',
            field_width,
            len(dna_cols)))
    liner.send(
        ' {:>{}}: {:,} Column(s)\n'.format(
            'Degenerate Columns',
            field_width,
            len(ddna_cols)))
    liner.send(
        ' {:>{}}: {:,} Column(s)\n'.format(
            'Metadata Columns',
            field_width,
            len(metadata_cols)))
    liner.send(
        ' {:>{}}: {:,} Column(s)\n'.format(
            'Mixed Columns',
            field_width,
            len(mixed_cols)))
    liner.send(
        ' {:>{}}: {:,} Column(s)\n'.format(
            'Non-string Columns',
            field_width,
            len(non_string_cols)))
    liner.send(
        ' {:>{}}: {:,} Duplicate(s)\n'.format(
            'Duplicate Names',
            field_width,
            len(duplicates)))
    liner.send(
        ' {:>{}}: {:,} Column(s)\n'.format(
            'Empty Seq. Values',
            field_width,
            len(empty_seq_cols)))

    if not seq_cols:
        liner.send(
            '   - No sequence columns detected [WARNING] (Length/motif checks skipped)\n')

    if seq_cols:
        prefix = '   - Sequence Column(s): '
        send_wrapped(
            prefix=prefix,
            text=', '.join(f"'{c}'" for c in seq_cols))

    if ddna_cols:
        prefix = '   - Degenerate Column(s): '
        send_wrapped(
            prefix=prefix,
            text=', '.join(f"'{c}'" for c in ddna_cols))
        liner.send(
            "   - Degenerate columns contain IUPAC bases [WARNING] (confirm this is intentional)\n")

    if metadata_cols:
        prefix = '   - Metadata Column(s): '
        send_wrapped(
            prefix=prefix,
            text=', '.join(f"'{c}'" for c in metadata_cols))

    if non_string_cols:
        prefix = '   - Non-string Column(s): '
        send_wrapped(
            prefix=prefix,
            text=', '.join(f"'{c}'" for c in non_string_cols))
        for colkey, info in non_string_examples.items():
            example_note = ut.format_row_examples(info.get('examples', []))
            liner.send(
                '   - Column {} is non-string [INFO] (dtype={} example_type={} example_value={}){}\n'.format(
                    colkey,
                    metadata_col_types.get(colkey, '<unknown>'),
                    info.get('example_type', '<unknown>'),
                    info.get('example_value', '<unavailable>'),
                    example_note))

    if mixed_cols:
        prefix = '   - Mixed Column(s): '
        send_wrapped(
            prefix=prefix,
            text=', '.join(f"'{c}'" for c in mixed_cols))
        for colkey, info in mixed_col_examples.items():
            example_note = ut.format_row_examples(info.get('examples', []))
            liner.send(
                '   - Column {} is mixed [WARNING] ({:5.1f} % DNA-like; example_value={!r}){}\n'.format(
                    colkey,
                    info.get('ddna_fraction', 0.0) * 100.,
                    info.get('example_value', '<unavailable>'),
                    example_note))

    if duplicates:
        prefix = '   - Duplicate Name(s): '
        send_wrapped(
            prefix=prefix,
            text=', '.join(f"'{c}'" for c in duplicates))

    if empty_seq_cols:
        prefix = '   - Empty Seq. Column(s): '
        send_wrapped(
            prefix=prefix,
            text=', '.join(f"'{c}'" for c in empty_seq_cols))
        for colkey, examples in empty_seq_examples.items():
            example_note = ut.format_row_examples(examples)
            liner.send(
                '   - Column {} has empty sequence value(s) [WARNING]{}\n'.format(
                    colkey,
                    example_note))

    if duplicates or mixed_cols or empty_seq_cols:
        liner.send(' Verdict: Column Inspection Completed with Violation(s)\n')
    else:
        liner.send(' Verdict: Column Inspection Completed\n')

    # Step 2: Length Stats
    intstats = cx.OrderedDict()
    minspaceavail = None
    maxspaceavail = None
    minoligolen = None
    maxoligolen = None

    if seq_cols:
        liner.send('\n[Step 2: Computing Length Statistics]\n')

        seqdf = indf.iloc[:, seq_positions].copy()

        (intstats,
        minspaceavail,
        maxspaceavail) = cl.lenstat_engine(
            indf=seqdf,
            oligolimit=oligolimit,
            liner=liner)

        if intstats:
            last_key = next(reversed(intstats))
            minoligolen = int(intstats[last_key][3])
            maxoligolen = int(intstats[last_key][4])

        liner.send('\n[Length Statistics]\n')

        statsprint = ut.get_lenstat_statsprint(
            intstats=intstats)

        if verbose and statsprint:
            print('\n{}'.format(statsprint))

        if oligolimit is None:
            liner.send(
                '\n Oligo Length Range: {:,} to {:,} Base Pair(s)\n'.format(
                    minoligolen,
                    maxoligolen))
            liner.send(
                ' Verdict: Oligo Length Range Computed\n')
        else:
            if minspaceavail == maxspaceavail:
                liner.send(
                    '\n Free Space Available: {:,} Base Pair(s)\n'.format(
                        int(minspaceavail)))
            else:
                liner.send(
                    '\n Free Space Available: {:,} to {:,} Base Pair(s)\n'.format(
                        int(minspaceavail),
                        int(maxspaceavail)))
            liner.send(
                ' Verdict: Oligo Length {} Limit\n'.format(
                    ('Within', 'Exceeds')[int(maxoligolen > oligolimit)]))

    # Step 3: Excluded motif emergence check
    motif_stats = None
    motif_ok = True

    if (exmotifs is not None) and seq_cols:
        liner.send('\n[Step 3: Checking Excluded Motifs]\n')

        seqdf = indf.iloc[:, seq_positions]
        if 'CompleteOligo' in seqdf.columns:
            complete = seqdf['CompleteOligo']
            if isinstance(complete, pd.DataFrame):
                complete = complete.iloc[:, 0]
            fullseqs = complete.str.replace('-', '').tolist()
        else:
            fullseqs = list(ut.get_df_concat(df=seqdf))

        motif_stats = {}
        for motif in exmotifs:
            counts = [seq.count(motif) for seq in fullseqs]
            baseline = int(min(counts)) if counts else 0
            maxcount = int(max(counts)) if counts else 0
            viol_mask = np.array([c > baseline for c in counts], dtype=bool)
            viol_rows = int(np.sum(viol_mask))
            examples = ut.get_row_examples(
                df=indf,
                invalid_mask=viol_mask,
                id_col='ID',
                limit=5)
            motif_stats[motif] = {
                'baseline_count': baseline,
                'max_count': maxcount,
                'emergent_rows': viol_rows,
                'examples': examples,
            }

            if viol_rows:
                motif_ok = False

        liner.send(
            ' Motifs Checked: {:,} Unique Motif(s)\n'.format(
                len(exmotifs)))

        prefix = ' Motif Set: '
        send_wrapped(
            prefix=prefix,
            text=', '.join(f"'{m}'" for m in exmotifs))

        liner.send(' Excluded Motif Emergence Summary\n')

        for motif in exmotifs:
            info = motif_stats[motif]
            example_note = ''
            tag = ''
            if info['emergent_rows']:
                tag = ' [WARNING]'
                example_note = ut.format_row_examples(info.get('examples', []))

            liner.send(
                "   - Motif '{}': {} to {} Occurrence(s); Emergent in {:,} Oligo(s){}{}\n".format(
                    motif,
                    info['baseline_count'],
                    info['max_count'],
                    info['emergent_rows'],
                    tag,
                    example_note))

        if motif_ok:
            liner.send(' Verdict: No Emergent Excluded Motif(s) Detected\n')
        else:
            liner.send(' Verdict: Emergent Excluded Motif(s) Detected\n')

    # Determine pass/fail
    length_overflow = False
    if (oligolimit is not None) and (maxoligolen is not None):
        length_overflow = bool(maxoligolen > oligolimit)

    status = True
    if duplicates:
        status = False
    if mixed_cols:
        status = False
    if empty_seq_cols:
        status = False
    if length_overflow:
        status = False
    if not motif_ok:
        status = False

    # Verification Statistics
    liner.send('\n[Verification Statistics]\n')
    field_width = 20

    liner.send(
        ' {:>{}}: {}\n'.format(
            'Verify Status',
            field_width,
            ('Successful', 'Failed')[int(not status)]))
    liner.send(
        ' {:>{}}: {:,} Record(s)\n'.format(
            'Input Records',
            field_width,
            input_rows))
    liner.send(
        ' {:>{}}: {:,} Column(s)\n'.format(
            'Sequence Columns',
            field_width,
            len(seq_cols)))
    liner.send(
        ' {:>{}}: {:,} Column(s)\n'.format(
            'Metadata Columns',
            field_width,
            len(metadata_cols)))
    liner.send(
        ' {:>{}}: {:,} Column(s)\n'.format(
            'Degenerate Columns',
            field_width,
            len(ddna_cols)))
    liner.send(
        ' {:>{}}: {:,} Column(s)\n'.format(
            'Mixed Columns',
            field_width,
            len(mixed_cols)))
    liner.send(
        ' {:>{}}: {:,} Duplicate(s)\n'.format(
            'Duplicate Names',
            field_width,
            len(duplicates)))
    liner.send(
        ' {:>{}}: {:,} Column(s)\n'.format(
            'Empty Seq. Values',
            field_width,
            len(empty_seq_cols)))

    if seq_cols and (minoligolen is not None) and (maxoligolen is not None):
        if minoligolen == maxoligolen:
            liner.send(
                ' {:>{}}: {:,} Base Pair(s)\n'.format(
                    'Oligo Length',
                    field_width,
                    minoligolen))
        else:
            liner.send(
                ' {:>{}}: {:,} to {:,} Base Pair(s)\n'.format(
                    'Oligo Length',
                    field_width,
                    minoligolen,
                    maxoligolen))
        if oligolimit is not None:
            liner.send(
                ' {:>{}}: {:,} Base Pair(s)\n'.format(
                    'Oligo Limit',
                    field_width,
                    oligolimit))
            liner.send(
                ' {:>{}}: {}\n'.format(
                    'Limit Overflow',
                    field_width,
                    ('No', 'Yes')[int(length_overflow)]))

    if exmotifs is not None:
        emergent_motifs = 0
        emergent_rows = 0
        if motif_stats:
            for info in motif_stats.values():
                if info.get('emergent_rows', 0):
                    emergent_motifs += 1
                    emergent_rows += int(info.get('emergent_rows', 0))

        liner.send(
            ' {:>{}}: {:,} Motif(s)\n'.format(
                'Excluded Motifs',
                field_width,
                len(exmotifs)))
        liner.send(
            ' {:>{}}: {:,} Motif(s)\n'.format(
                'Emergent Motifs',
                field_width,
                emergent_motifs))
        if emergent_motifs:
            liner.send(
                ' {:>{}}: {:,} Oligo(s)\n'.format(
                    'Emergent Rows',
                    field_width,
                    emergent_rows))

    # Build Stats Dictionary
    stats = {
        'status'  : status,
        'basis'   : ('solved', 'violations')[int(not status)],
        'step'    : 3 if exmotifs is not None else 2,
        'step_name': 'verifying-input-data',
        'vars'    : {
            'sequence_columns': seq_cols,
            'dna_columns': dna_cols,
            'degenerate_columns': ddna_cols,
            'metadata_columns': metadata_cols,
            'mixed_columns': mixed_cols,
            'non_string_columns': non_string_cols,
            'duplicate_columns': duplicates,
            'empty_sequence_columns': empty_seq_cols,
            'oligo_limit': oligolimit,
            'min_oligo_len': minoligolen,
            'max_oligo_len': maxoligolen,
            'min_space_avail': minspaceavail,
            'max_space_avail': maxspaceavail,
            'length_overflow': length_overflow,
            'len_stat': ut.get_lenstat_dict(
                intstats=intstats),
            'excluded_motifs': exmotifs,
            'excluded_motif_stats': motif_stats,
        },
        'warns'   : warns,
    }

    # Attach examples for debugging (kept in vars for JSON friendliness)
    stats['vars']['metadata_column_types'] = metadata_col_types
    stats['vars']['mixed_column_examples'] = mixed_col_examples
    stats['vars']['non_string_examples'] = non_string_examples
    stats['vars']['empty_sequence_examples'] = empty_seq_examples

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Close Liner
    liner.close()

    # Return Results
    stats = ut.stamp_stats(
        stats=stats,
        module='verify',
        input_rows=input_rows,
        output_rows=0)
    return stats
