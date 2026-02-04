import time as tt

import collections as cx

import textwrap as tw

import numpy as np
import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import vectordb as db
from .base import core_lenstat as cl


def verify(
    input_data:str|pd.DataFrame,
    oligo_length_limit:int|None=None,
    excluded_motifs:list|str|pd.DataFrame|None=None,
    background_directory:str|None=None,
    verbose:bool=True) -> dict:
    '''
    Lightweight QC for an oligo pool DataFrame (column inspection, length stats, and excluded-motif
    emergence), without modifying data or writing outputs.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with an 'ID' column.

    Optional Parameters:
        - `oligo_length_limit` (`int` / `None`): If provided, check for length overflow (default: `None`).
        - `excluded_motifs` (`list` / `str` / `pd.DataFrame` / `None`): Motifs to track; can be a CSV path
            or DataFrame with an 'Exmotif' column (default: `None`).
        - `background_directory` (`str` / `None`): Path to background k-mer database created by
            `background()`. If provided, scans concatenated oligos for any background k-mers and reports
            violations (default: `None`).
        - `verbose` (`bool`): If `True`, logs progress to stdout (default: `True`).

    Returns:
        - A dictionary of verification results (stats only; no DataFrame is returned).

    Notes:
        - Column inspection classifies columns as sequence vs metadata, and flags mixed/degenerate/non-string
            columns (degenerate/IUPAC columns are warnings, not hard errors).
        - Length stats use the same engine as `lenstat` and can be checked against `oligo_length_limit`.
        - Motif emergence means a motif occurs more often than the minimum occurrence across the library
            (useful when a motif should appear exactly once, e.g., restriction sites).
        - Excluded-motif counting is literal substring matching; degenerate/IUPAC bases are not expanded
            as wildcards during motif checks.
        - Motif scans operate on concatenated sequence columns (left-to-right DataFrame order) with `'-'`
            removed. If `CompleteOligo` exists (from `final`), it is used directly.
        - If emergent motifs are detected and multiple sequence columns are present, `verify` attributes
            emergence to column junctions ("edge effects"), unless `CompleteOligo` is provided.
    '''

    # Argument Aliasing
    indata     = input_data
    oligolimit = oligo_length_limit
    exmotifs   = excluded_motifs
    background = background_directory
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
        data_field='      Input Data    ',
        required_fields=None,
        liner=liner)
    input_rows = len(indf.index) if isinstance(indf, pd.DataFrame) else 0

    # Show update on successful parse
    if indata_valid:
        liner.send(
            '      Input Data    : {} w/ {:,} Record(s)\n'.format(
                data_name,
                len(indf.index)))

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full oligolimit Validation (optional)
    oligolimit_valid = True
    if oligolimit is None:
        liner.send('      Oligo Limit   : None Specified\n')
    else:
        oligolimit_valid = vp.get_numeric_validity(
            numeric=oligolimit,
            numeric_field='      Oligo Limit   ',
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
        exseqs_field='   Excluded Motifs  ',
        exseqs_desc='Unique Motif(s)',
        df_field='Exmotif',
        required=False,
        liner=liner)

    # Full background Parsing and Validation
    (background_valid,
    background_type) = vp.get_parsed_background(
        background=background,
        background_field=' Background Database',
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid,
        exmotifs_valid,
        background_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Open Background
    if background_type == 'path':
        background = ut.get_adjusted_path(
            path=background,
            suffix='.oligopool.background')
        background = db.vectorDB(
            path=background,
            maximum_repeat_length=None)

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
    motif_viol_masks = None
    junction_stats = None

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
        motif_viol_masks = {}
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
            motif_viol_masks[motif] = viol_mask

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

        # Step 4: Localize emergent motifs to column junctions (edge effects).
        emergent_motifs = [m for m, info in motif_stats.items() if info.get('emergent_rows', 0)]
        if emergent_motifs:
            liner.send('\n[Step 4: Localizing Emergent Motifs]\n')

            if 'CompleteOligo' in seqdf.columns:
                liner.send(
                    ' Verdict: Junction Attribution Skipped [INFO] (Run verify before `final` to keep separate columns)\n')
            else:
                # Junctions are defined by the order of sequence columns in the DataFrame.
                junction_cols = list(seqdf.columns)

                if len(junction_cols) < 2:
                    liner.send(
                        ' Verdict: Junction Attribution Skipped [INFO] (Need ≥ 2 Sequence Columns)\n')
                else:
                    junction_keys = [
                        '{}|{}'.format(junction_cols[i], junction_cols[i + 1])
                        for i in range(len(junction_cols) - 1)
                    ]

                    junction_stats = {
                        'junction_columns': junction_cols,
                        'motifs': {},
                    }

                    for motif in emergent_motifs:
                        viol_mask = motif_viol_masks.get(motif)
                        viol_idxs = np.where(viol_mask)[0] if viol_mask is not None else np.array([], dtype=int)

                        # Row-level counts: for each junction, number of violating rows with ≥1 junction-spanning hit.
                        counts_by_junction = {k: 0 for k in junction_keys}
                        examples_by_junction = {k: [] for k in junction_keys}
                        any_junction_rows = 0
                        internal_only_rows = 0

                        mot_len = len(motif)
                        baseline = int(motif_stats.get(motif, {}).get('baseline_count', 0))

                        for ridx in viol_idxs:
                            row_vals = seqdf.iloc[ridx, :].tolist()
                            parts = [v.replace('-', '') for v in row_vals]

                            # Boundary positions in the assembled (gap-stripped) oligo.
                            boundaries = []
                            pos = 0
                            for part in parts:
                                pos += len(part)
                                boundaries.append(pos)
                            boundaries = boundaries[:-1]  # N cols -> N-1 junctions

                            full = ''.join(parts)

                            # Find non-overlapping motif occurrences (matches `str.count` semantics used above).
                            occ = []
                            start = 0
                            while True:
                                p = full.find(motif, start)
                                if p < 0:
                                    break
                                occ.append(p)
                                start = p + mot_len

                            crossed = set()
                            # Attribute only emergent occurrences beyond the baseline minimum.
                            for p in occ[baseline:]:
                                end = p + mot_len
                                for j, bpos in enumerate(boundaries):
                                    if p < bpos < end:
                                        crossed.add(j)

                            if crossed:
                                any_junction_rows += 1
                                for j in crossed:
                                    jkey = junction_keys[j]
                                    counts_by_junction[jkey] += 1
                                    if len(examples_by_junction[jkey]) < 5:
                                        examples_by_junction[jkey].append(str(indf.index[ridx]))
                            else:
                                internal_only_rows += 1

                        # Print motif-localization summary
                        liner.send(" Motif '{}': Junction Attribution (Beyond Baseline={})\n".format(
                            motif,
                            baseline))

                        # Only show junctions that contribute at least once (sorted descending).
                        nonzero = [(k, v) for k, v in counts_by_junction.items() if v]
                        if nonzero:
                            nonzero.sort(key=lambda kv: kv[1], reverse=True)
                            for jkey, count in nonzero:
                                left, right = jkey.split('|', 1)
                                example_note = ut.format_row_examples(examples_by_junction.get(jkey, []))
                                liner.send(
                                    "   - Junction '{}' + '{}': {:,} Oligo(s){}\n".format(
                                        left,
                                        right,
                                        count,
                                        example_note))
                        if internal_only_rows:
                            liner.send(
                                '   - Internal (Within One Column): {:,} Oligo(s)\n'.format(
                                    internal_only_rows))

                        liner.send(
                            '   - Junction-Spanning Rows: {:,} Oligo(s)\n'.format(
                                any_junction_rows))

                        junction_stats['motifs'][motif] = {
                            'baseline_count': int(baseline),
                            'any_junction_rows': int(any_junction_rows),
                            'internal_only_rows': int(internal_only_rows),
                            'junction_rows': {k: int(v) for k, v in counts_by_junction.items() if v},
                            'junction_examples': {k: v for k, v in examples_by_junction.items() if v},
                        }

                    liner.send(' Verdict: Motif Junction Attribution Completed\n')

    # Step 4: Background k-mer scan
    background_violations = 0
    background_violation_rows = []
    background_scan_performed = False

    if background_type is not None and seq_cols:
        background_scan_performed = True
        liner.send('\n[Step 4: Checking Background k-mers]\n')

        # Build concatenated oligos
        if 'CompleteOligo' in seq_cols:
            oligos = indf['CompleteOligo'].str.replace('-', '').tolist()
        else:
            seqdf = indf[seq_cols]
            oligos = seqdf.apply(
                lambda row: ''.join(str(v).replace('-', '') for v in row),
                axis=1).tolist()

        # Scan each oligo for background k-mers
        liner.send(' Scanning {:,} Oligo(s) for Background k-mers (k={:,})...\n'.format(
            len(oligos),
            background.K))

        for ridx, oligo in enumerate(oligos):
            if len(oligo) < background.K:
                continue
            # Check all k-mers in oligo
            for i in range(len(oligo) - background.K + 1):
                kmer = oligo[i:i+background.K]
                if kmer in background:
                    background_violations += 1
                    if len(background_violation_rows) < 10:
                        background_violation_rows.append(str(indf.index[ridx]))
                    break  # Only count each oligo once

        if background_violations:
            example_note = ut.format_row_examples(background_violation_rows)
            liner.send(' Found {:,} Oligo(s) with Background k-mer Violation(s){}\n'.format(
                background_violations,
                example_note))
            liner.send(' Verdict: Background Scan Completed with Violation(s)\n')
        else:
            liner.send(' Found 0 Oligo(s) with Background k-mer Violation(s)\n')
            liner.send(' Verdict: Background Scan Completed\n')

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
    if background_violations:
        status = False

    # Verification Statistics
    liner.send('\n[Verification Statistics]\n')
    field_width = 21

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

    if background_scan_performed:
        liner.send(
            ' {:>{}}: {:,}\n'.format(
                'Background k',
                field_width,
                background.K))
        liner.send(
            ' {:>{}}: {:,} Oligo(s)\n'.format(
                'Background Violations',
                field_width,
                background_violations))

    # Determine step number
    step_num = 2  # base: length stats
    if exmotifs is not None:
        step_num = 3
    if background_scan_performed:
        step_num = 4

    # Build Stats Dictionary
    stats = {
        'status'  : status,
        'basis'   : ('solved', 'violations')[int(not status)],
        'step'    : step_num,
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
            'excluded_motif_junction_stats': junction_stats,
            'background_k': background.K if background_scan_performed else None,
            'background_violations': background_violations,
            'background_violation_examples': background_violation_rows,
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

    # Close Background
    if background_type == 'path':
        background.close()

    # Close Liner
    liner.close()

    # Return Results
    stats = ut.stamp_stats(
        stats=stats,
        module='verify',
        input_rows=input_rows,
        output_rows=0)
    return stats
