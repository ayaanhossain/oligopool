import pandas as pd

from . import utils as ut


def verify_engine(
    indf,
    dna_columns,
    constituent_columns,
    oligolimit,
    exmotifs,
    backgrounds,
    liner):
    '''
    Scan an oligo pool for integrity, length,
    excluded motif emergence, and background
    k-mer conflicts.

    :: indf
       type - pd.DataFrame
       desc - a pandas DataFrame storing
              annotated oligo pool variants
    :: dna_columns
       type - list
       desc - list of DNA column names detected
              in the input DataFrame
    :: constituent_columns
       type - list
       desc - list of constituent DNA column names
              (excludes CompleteOligo); may be empty
    :: oligolimit
       type - integer
       desc - maximum oligo length allowed in the
              oligopool, must be 4 or greater
    :: exmotifs
       type - list / None
       desc - list of excluded motif strings (ATGC),
              or None if not specified
    :: backgrounds
       type - list
       desc - list of opened vectorDB background
              databases, or empty list
    :: liner
       type - coroutine
       desc - the liner engine for progress output

    ..Internal use only...
    '''

    # Book-keeping
    num_rows = len(indf)
    plen = ut.get_printlen(value=num_rows)

    # Check if CompleteOligo was detected
    used_complete_oligo = (dna_columns == ['CompleteOligo'])

    # Column list for position attribution:
    # - Scenario A (CompleteOligo + constituents): use constituent columns
    # - Scenario B (CompleteOligo only): use ['CompleteOligo'] as single column
    # - Scenario C (no CompleteOligo): use constituent columns (same as dna_columns)
    attribution_columns = constituent_columns if constituent_columns else ['CompleteOligo']
    # Ordered column list for conflict dicts; may grow to include
    # CompleteOligo when integrity-mismatch rows override boundaries.
    conflict_columns = list(attribution_columns)

    # Report all DNA columns
    all_dna_columns = list(constituent_columns)
    if used_complete_oligo:
        all_dna_columns.append('CompleteOligo')
    liner.send(' DNA Column(s): {} Column(s)\n'.format(len(all_dna_columns)))
    for col in all_dna_columns:
        liner.send('   - {}\n'.format(col))

    # Build concatenated oligos
    oligo_sequences = []
    oligo_lengths = []

    for idx in range(num_rows):
        if (idx + 1) % 1000 == 0 or idx + 1 == num_rows:
            liner.send('     Oligos Constructed: Scanned {:>{},} Variant(s)'.format(idx + 1, plen))
        if used_complete_oligo:
            oligo = str(indf['CompleteOligo'].iloc[idx]).replace('-', '').upper()
        else:
            parts = [str(indf[col].iloc[idx]).replace('-', '').upper() for col in dna_columns]
            oligo = ''.join(parts)
        oligo_sequences.append(oligo)
        oligo_lengths.append(len(oligo))
    liner.send('     Oligos Constructed: Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

    # Consistency check: CompleteOligo vs constituents (Scenario A only)
    complete_oligo_conflicts = [False] * num_rows
    if used_complete_oligo and constituent_columns:
        for idx in range(num_rows):
            if (idx + 1) % 1000 == 0 or idx + 1 == num_rows:
                liner.send(' CompleteOligo Verified: Scanned {:>{},} Variant(s)'.format(idx + 1, plen))
            concat = ''.join(
                str(indf[col].iloc[idx]).replace('-', '').upper()
                for col in constituent_columns)
            if oligo_sequences[idx] != concat:
                complete_oligo_conflicts[idx] = True
        liner.send(' CompleteOligo Verified: Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

    # OligoLength check (all scenarios, when column present)
    oligo_length_conflicts = [False] * num_rows
    oligo_length_series = None
    if 'OligoLength' in indf.columns:
        oligo_length_series = pd.to_numeric(indf['OligoLength'], errors='coerce')
        if oligo_length_series.notna().all():
            for idx in range(num_rows):
                if (idx + 1) % 1000 == 0 or idx + 1 == num_rows:
                    liner.send('   OligoLength Verified: Scanned {:>{},} Variant(s)'.format(idx + 1, plen))
                expected = int(oligo_length_series.iloc[idx])
                actual = oligo_lengths[idx]
                if expected != actual:
                    oligo_length_conflicts[idx] = True
            liner.send('   OligoLength Verified: Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

    # Build per-row integrity conflict data
    integrity_conflicts = []
    integrity_details = []
    for idx in range(num_rows):
        co = complete_oligo_conflicts[idx]
        ol = oligo_length_conflicts[idx]
        if co or ol:
            integrity_conflicts.append(True)
            integrity_details.append({
                'complete_oligo_conflict': co,
                'oligo_length_conflict': ol,
                'expected_oligo_length': int(oligo_length_series.iloc[idx]) if ol else None,
                'actual_oligo_length': oligo_lengths[idx] if ol else None,
            })
        else:
            integrity_conflicts.append(False)
            integrity_details.append(None)

    integrity_conflict_count = sum(integrity_conflicts)
    complete_oligo_conflict_count = sum(complete_oligo_conflicts)
    oligo_length_conflict_count = sum(oligo_length_conflicts)

    min_oligo_length = min(oligo_lengths)
    max_oligo_length = max(oligo_lengths)

    # Precompute per-row column boundaries for attribution
    row_boundaries = []
    for idx in range(num_rows):
        if (idx + 1) % 1000 == 0 or idx + 1 == num_rows:
            liner.send('     Bounds Constructed: Scanned {:>{},} Variant(s)'.format(idx + 1, plen))
        row_boundaries.append(
            ut.get_column_boundaries(idx, attribution_columns, indf))
    liner.send('     Bounds Constructed: Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

    # Override boundaries for rows with CompleteOligo integrity conflicts
    # (constituent boundaries are unreliable when CompleteOligo != concat)
    for idx in range(num_rows):
        if complete_oligo_conflicts[idx]:
            row_boundaries[idx] = [('CompleteOligo', 0, oligo_lengths[idx])]
            if 'CompleteOligo' not in conflict_columns:
                conflict_columns.append('CompleteOligo')

    # Length Conflict Detection
    length_conflicts = []
    length_details = []

    for idx, length in enumerate(oligo_lengths):
        if (idx + 1) % 1000 == 0 or idx + 1 == num_rows:
            liner.send('    Lengths Evaluated  : Scanned {:>{},} Variant(s)'.format(idx + 1, plen))
        if length > oligolimit:
            length_conflicts.append(True)
            # Build column_lengths from boundaries
            column_lengths = {}
            for col_name, start, end in row_boundaries[idx]:
                column_lengths[col_name] = end - start
            length_details.append({
                'oligo_length': length,
                'length_limit': oligolimit,
                'excess_length': length - oligolimit,
                'column_lengths': column_lengths
            })
        else:
            length_conflicts.append(False)
            length_details.append(None)
    liner.send('    Lengths Evaluated  : Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

    # Exmotif Emergence Detection
    exmotif_conflicts = []
    exmotif_details = []
    baselines = {}
    emergence_summary = {}
    excluded_motif_column_conflicts = None

    if exmotifs:
        # Pre-compute baselines
        num_motifs = len(exmotifs)
        mlen = ut.get_printlen(value=num_motifs)
        for midx, motif in enumerate(exmotifs):
            if (midx + 1) % 10 == 0 or midx + 1 == num_motifs:
                liner.send('   Baseline Evaluated  : Scanned {:>{},} Motif(s)'.format(midx + 1, mlen))
            counts = [oligo.count(motif) for oligo in oligo_sequences]
            baselines[motif] = min(counts) if counts else 0
        liner.send('   Baseline Evaluated  : Scanned {:>{},} Motif(s)\n'.format(num_motifs, mlen))

        # Per-row detection
        for row_idx, oligo in enumerate(oligo_sequences):
            if (row_idx + 1) % 1000 == 0 or row_idx + 1 == num_rows:
                liner.send('   Exmotifs Evaluated  : Scanned {:>{},} Variant(s)'.format(row_idx + 1, plen))
            row_details = []
            boundaries = row_boundaries[row_idx]
            for motif in exmotifs:
                count = oligo.count(motif)
                if count > baselines[motif]:
                    # Find all positions
                    positions = []
                    columns = []
                    start = 0
                    while True:
                        pos = oligo.find(motif, start)
                        if pos < 0:
                            break
                        positions.append(pos)
                        columns.append(ut.get_boundary_column(pos, boundaries))
                        # Non-overlapping semantics (match str.count()).
                        start = pos + len(motif)
                    row_details.append({
                        'motif': motif,
                        'occurrences': count,
                        'library_baseline': baselines[motif],
                        'excess_occurrences': count - baselines[motif],
                        'positions': positions,
                        'columns': columns
                    })
            if row_details:
                exmotif_conflicts.append(True)
                exmotif_details.append(row_details)
            else:
                exmotif_conflicts.append(False)
                exmotif_details.append(None)
        liner.send('   Exmotifs Evaluated  : Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

        # Build emergence summary for stats
        for midx, motif in enumerate(exmotifs):
            if (midx + 1) % 10 == 0 or midx + 1 == num_motifs:
                liner.send('  Emergence Evaluated  : Scanned {:>{},} Motif(s)'.format(midx + 1, mlen))
            emergent_rows = sum(
                1 for d in exmotif_details
                if d and any(e['motif'] == motif for e in d))
            max_emergence = max(
                (e['excess_occurrences'] for d in exmotif_details
                 if d for e in d if e['motif'] == motif),
                default=0)
            emergence_summary[motif] = {
                'emergent_row_count': emergent_rows,
                'max_excess_occurrences': max_emergence
            }
        liner.send('  Emergence Evaluated  : Scanned {:>{},} Motif(s)\n'.format(num_motifs, mlen))

        # Compute excluded_motif_column_conflicts: per-column count of oligos with exmotif hits
        excluded_motif_column_conflicts = {col: 0 for col in conflict_columns}
        for row_idx, detail in enumerate(exmotif_details):
            if (row_idx + 1) % 1000 == 0 or row_idx + 1 == num_rows:
                liner.send('     Column Conflicts  : Scanned {:>{},} Variant(s)'.format(row_idx + 1, plen))
            if detail is None:
                continue
            # Collect all columns hit across all motifs in this row
            cols_hit = set()
            for entry in detail:
                cols_hit.update(entry['columns'])
            for col in cols_hit:
                if col in excluded_motif_column_conflicts:
                    excluded_motif_column_conflicts[col] += 1
        liner.send('     Column Conflicts  : Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))
    else:
        exmotif_conflicts = [False] * num_rows
        exmotif_details = [None] * num_rows

    # Background K-mer Conflict Detection
    background_conflicts = []
    background_details = []
    background_column_conflicts = None

    if backgrounds:
        for row_idx, oligo in enumerate(oligo_sequences):
            if (row_idx + 1) % 1000 == 0 or row_idx + 1 == num_rows:
                liner.send(' Background Evaluated  : Scanned {:>{},} Variant(s)'.format(row_idx + 1, plen))
            row_details = {
                'matched_kmers': [],
                'positions': [],
                'columns': [],
                'match_count': 0
            }
            has_conflict = False
            boundaries = row_boundaries[row_idx]

            for bg in backgrounds:
                K = bg.K
                if len(oligo) < K:
                    continue
                for i in range(len(oligo) - K + 1):
                    kmer = oligo[i:i+K]
                    if kmer in bg:
                        has_conflict = True
                        row_details['matched_kmers'].append(kmer)
                        row_details['positions'].append(i)
                        row_details['columns'].append(
                            ut.get_boundary_column(i, boundaries))
                        row_details['match_count'] += 1

            if has_conflict:
                background_conflicts.append(True)
                background_details.append(row_details)
            else:
                background_conflicts.append(False)
                background_details.append(None)
        liner.send(' Background Evaluated  : Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

        # Compute background_column_conflicts
        background_column_conflicts = {col: 0 for col in conflict_columns}
        for row_idx, detail in enumerate(background_details):
            if (row_idx + 1) % 1000 == 0 or row_idx + 1 == num_rows:
                liner.send('     Column Conflicts  : Scanned {:>{},} Variant(s)'.format(row_idx + 1, plen))
            if detail is None:
                continue
            cols_hit = set(detail['columns'])
            for col in cols_hit:
                if col in background_column_conflicts:
                    background_column_conflicts[col] += 1
        liner.send('     Column Conflicts  : Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))
    else:
        background_conflicts = [False] * num_rows
        background_details = [None] * num_rows

    # Build HasAnyConflicts
    has_any_conflicts = [
        l or e or b or i
        for l, e, b, i in zip(length_conflicts, exmotif_conflicts,
                              background_conflicts, integrity_conflicts)
    ]

    # Compute any_column_conflicts: per-column count of rows with any
    # position-level conflict (exmotif or background) in that column
    any_column_conflicts = None
    if excluded_motif_column_conflicts is not None or background_column_conflicts is not None:
        any_column_conflicts = {col: 0 for col in conflict_columns}
        for row_idx in range(num_rows):
            if (row_idx + 1) % 1000 == 0 or row_idx + 1 == num_rows:
                liner.send('   Conflict Aggregation: Scanned {:>{},} Variant(s)'.format(row_idx + 1, plen))
            cols_hit = set()
            ex_detail = exmotif_details[row_idx]
            if ex_detail:
                for entry in ex_detail:
                    cols_hit.update(entry['columns'])
            bg_detail = background_details[row_idx]
            if bg_detail:
                cols_hit.update(bg_detail['columns'])
            for col in cols_hit:
                if col in any_column_conflicts:
                    any_column_conflicts[col] += 1
        liner.send('   Conflict Aggregation: Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

    # Return Results
    return {
        # Per-row lists for DataFrame columns
        'oligo_sequences':        oligo_sequences,
        'oligo_lengths':          oligo_lengths,
        'integrity_conflicts':    integrity_conflicts,
        'length_conflicts':       length_conflicts,
        'exmotif_conflicts':      exmotif_conflicts,
        'background_conflicts':   background_conflicts,
        'has_any_conflicts':      has_any_conflicts,
        'integrity_details':      integrity_details,
        'length_details':         length_details,
        'exmotif_details':        exmotif_details,
        'background_details':     background_details,
        # Scalars for stats
        'used_complete_oligo':            used_complete_oligo,
        'attribution_columns':            attribution_columns,
        'min_oligo_length':               min_oligo_length,
        'max_oligo_length':               max_oligo_length,
        'integrity_conflict_count':       integrity_conflict_count,
        'complete_oligo_conflict_count':  complete_oligo_conflict_count,
        'oligo_length_conflict_count':    oligo_length_conflict_count,
        'baselines':                      baselines,
        'emergence_summary':              emergence_summary,
        # Column conflict dicts
        'excluded_motif_column_conflicts': excluded_motif_column_conflicts,
        'background_column_conflicts':     background_column_conflicts,
        'any_column_conflicts':            any_column_conflicts,
    }
