import time as tt

import atexit as ae

import numpy as np
import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_degenerate as cd

from typing import Tuple


def compress(
    input_data:str|pd.DataFrame,
    mapping_file:str|None=None,
    synthesis_file:str|None=None,
    rollout_simulations:int=100,
    rollout_horizon:int=4,
    random_seed:int|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    '''
    Compress concrete DNA sequences into IUPAC-degenerate oligos for `Degenerate Mode`.

    This reduces synthesis cost for large variant libraries by grouping similar concrete
    sequences into fewer degenerate representations while preserving the full set of
    originals (lossless compression).

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with
            variant sequences. Must contain a unique 'ID' column; all other columns
            are treated as strict-ATGC DNA (A/T/G/C only) and concatenated.

    Optional Parameters:
        - `mapping_file` (`str`): Filename for output mapping DataFrame that links each
            original variant ID to its assigned degenerate oligo ID (default: `None`).
            A `.oligopool.compress.csv` suffix is added if missing.
        - `synthesis_file` (`str`): Filename for output synthesis DataFrame containing
            the degenerate oligos ready for ordering (default: `None`).
            A `.oligopool.compress.csv` suffix is added if missing.
        - `rollout_simulations` (`int`): Number of Monte Carlo simulations per decision
            during compression. Higher values give better compression but take longer
            (default: `100`).
        - `rollout_horizon` (`int`): Number of positions to look ahead during Monte Carlo
            rollouts. Higher values can improve compression quality (default: `4`).
        - `random_seed` (`int` / `None`): Seed for local RNG (default: `None`).
        - `verbose` (`bool`): If `True`, logs progress updates to stdout (default: `True`).

    Returns:
        - `mapping_df` (`pd.DataFrame`): DataFrame mapping original variant IDs to
            degenerate oligo IDs. Columns: `ID`, `Sequence`, `DegenerateID`.
        - `synthesis_df` (`pd.DataFrame`): DataFrame of degenerate oligos for synthesis.
            Columns: `DegenerateID`, `DegenerateSeq`, `Degeneracy`, `OligoLength`.
        - `stats` (`dict`): Statistics dictionary with compression results including
            `compression_ratio`, `min_degeneracy`, `max_degeneracy`, `mean_degeneracy`.

    Notes:
        - Compression is lossless: expanding the degenerate oligos recovers exactly the
          set of unique input sequences (no extras, no missing).
        - Sequences of different lengths are compressed independently by length group.
        - `rollout_simulations` and `rollout_horizon` trade off runtime vs compression ratio.
        - `mapping_df` is your traceability map for mapping sequenced survivors back to
          original `ID`s; `expand` is useful as a sanity check for compression.
    '''

    # Preserve return style when the caller intentionally used ID as index.
    id_from_index = ut.get_id_index_intent(input_data)

    # Argument Aliasing
    indata      = input_data
    mapfile     = mapping_file
    synfile     = synthesis_file
    nsims       = rollout_simulations
    horizon     = rollout_horizon
    random_seed = random_seed
    verbose     = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Compression Verbiage Print
    liner.send('\n[Oligopool Calculator: Degenerate Mode - Compress]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # Parse and validate input data
    (indf,
     dna_cols,
     indata_valid) = vp.get_parsed_compress_data_info(
        data=indata,
        data_field='      Input Data       ',
        liner=liner)

    input_rows = len(indf.index) if isinstance(indf, pd.DataFrame) else 0

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Validate mapping file
    mapfile_valid = vp.get_outdf_validity(
        outdf=mapfile,
        outdf_suffix='.oligopool.compress.csv',
        outdf_field='    Mapping File       ',
        liner=liner)

    # Validate synthesis file
    synfile_valid = vp.get_outdf_validity(
        outdf=synfile,
        outdf_suffix='.oligopool.compress.csv',
        outdf_field='  Synthesis File       ',
        liner=liner)

    # Validate rollout_simulations
    nsims_valid = vp.get_numeric_validity(
        numeric=nsims,
        numeric_field='    Rollout Simulations',
        numeric_pre_desc=' ',
        numeric_post_desc=' Rollout(s) per Decision',
        minval=1,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Validate rollout_horizon
    horizon_valid = vp.get_numeric_validity(
        numeric=horizon,
        numeric_field='    Rollout Horizon    ',
        numeric_pre_desc=' ',
        numeric_post_desc=' Position Lookahead',
        minval=1,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Validate random_seed (do not auto-generate)
    seed_valid = True
    if random_seed is None:
        liner.send('     Random Seed       : None Specified\n')
    elif isinstance(random_seed, (int, np.integer)):
        random_seed = int(random_seed)
        liner.send(
            '     Random Seed       : {:,} ... Assigned\n'.format(
                random_seed))
    else:
        liner.send(
            '     Random Seed       : {} [MUST BE INTEGER OR NONE]\n'.format(
                random_seed))
        seed_valid = False

    # First Pass Validation
    if not all([
        indata_valid,
        mapfile_valid,
        synfile_valid,
        nsims_valid,
        horizon_valid,
        seed_valid]):
        liner.send('\n')
        raise RuntimeError('Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Parameters
    nsims = round(nsims)
    horizon = round(horizon)

    # Setup RNG
    rng = np.random.default_rng(random_seed)

    # Setup stats
    stats = {
        'status': False,
        'basis': 'unsolved',
        'step': 0,
        'step_name': 'initializing',
        'vars': {},
        'warns': {},
    }

    # Schedule file cleanup on error
    if mapfile is not None:
        mapfile = ut.get_adjusted_path(
            path=mapfile,
            suffix='.oligopool.compress.csv')
        ae.register(
            ut.remove_file,
            mapfile)

    if synfile is not None:
        synfile = ut.get_adjusted_path(
            path=synfile,
            suffix='.oligopool.compress.csv')
        ae.register(
            ut.remove_file,
            synfile)

    # Step 1: Parse oligopool
    liner.send('\n[Step 1: Parsing Oligopool]\n')
    stats['step'] = 1
    stats['step_name'] = 'parsing-oligopool'

    # Concatenate DNA columns
    liner.send(' Concatenating {:,} DNA Column(s) ...'.format(len(dna_cols)))
    sequences = ut.get_df_concat(df=indf[dna_cols])
    variant_ids = list(indf.index)
    id_to_sequence = dict(zip(variant_ids, sequences))

    # Count unique sequences
    unique_sequences = set(sequences)
    liner.send(' Unique Sequence(s): {:,} Variant(s)\n'.format(len(unique_sequences)))

    # Group by length
    length_groups = {}
    for seq, vid in zip(sequences, variant_ids):
        seq_len = len(seq)
        if seq_len not in length_groups:
            length_groups[seq_len] = {'sequences': [], 'ids': []}
        length_groups[seq_len]['sequences'].append(seq)
        length_groups[seq_len]['ids'].append(vid)

    # Report length groups
    min_len = min(length_groups.keys())
    max_len = max(length_groups.keys())
    if min_len == max_len:
        liner.send(' Sequence Length(s): {:,} Base Pair(s)\n'.format(min_len))
    else:
        liner.send(' Sequence Length(s): {:,}-{:,} Base Pair(s) ({:,} Length Group(s))\n'.format(
            min_len, max_len, len(length_groups)))

    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Step 2: Compress sequences
    liner.send('\n[Step 2: Compressing Sequences]\n')
    stats['step'] = 2
    stats['step_name'] = 'compressing-sequences'

    all_results = []  # List of (degenerate_seq, covered_ids, degeneracy, length)
    total_input_variants = 0
    total_degenerate_oligos = 0

    for seq_length in sorted(length_groups.keys()):
        group = length_groups[seq_length]
        group_seqs = group['sequences']
        group_ids = group['ids']

        liner.send(' Length Group {:,} bp: {:,} Variant(s)\n'.format(
            seq_length, len(group_seqs)))

        # Compress this length group
        results, covered = cd.get_compressed_group(
            sequences=group_seqs,
            ids=group_ids,
            nsims=nsims,
            horizon=horizon,
            rng=np.random.default_rng(rng.integers(0, 2**31 - 1)),
            liner=liner,
            group_length=seq_length,
        )

        # Collect results with length info
        for (deg_seq, covered_ids, degeneracy) in results:
            all_results.append((deg_seq, covered_ids, degeneracy, seq_length))

        # Compute compression ratio for this group
        unique_in_group = len(set(group_seqs))
        compression_ratio = unique_in_group / len(results) if results else 1.0

        liner.send(' Group Complete: {:,} -> {:,} Degenerate Oligo(s) ({:.1f}x Compression)\n'.format(
            unique_in_group, len(results), compression_ratio))

        total_input_variants += unique_in_group
        total_degenerate_oligos += len(results)

    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Step 3: Build output DataFrames
    liner.send('\n[Step 3: Building Output DataFrames]\n')
    stats['step'] = 3
    stats['step_name'] = 'building-output'

    # Build mapping DataFrame
    mapping_rows = []
    synthesis_rows = []
    degenerate_idx = 0

    for (deg_seq, covered_ids, degeneracy, seq_length) in all_results:
        degenerate_idx += 1
        deg_id = 'D{}'.format(degenerate_idx)

        # Build sequence for each covered variant ID
        for vid in covered_ids:
            orig_seq = id_to_sequence.get(vid)
            mapping_rows.append({
                'ID': vid,
                'Sequence': orig_seq,
                'DegenerateID': deg_id,
            })

        # Build synthesis row
        synthesis_rows.append({
            'DegenerateID': deg_id,
            'DegenerateSeq': deg_seq,
            'Degeneracy': degeneracy,
            'OligoLength': seq_length,
        })

    mapping_df = pd.DataFrame(mapping_rows)
    synthesis_df = pd.DataFrame(synthesis_rows)

    liner.send('   Mapping DataFrame: {:,} Row(s)\n'.format(len(mapping_df)))
    liner.send(' Synthesis DataFrame: {:,} Row(s)\n'.format(len(synthesis_df)))
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Write outputs
    if mapfile is not None:
        ut.write_df_csv(
            df=mapping_df,
            outfile=mapfile,
            sep=',')
    if synfile is not None:
        synthesis_df.to_csv(
            path_or_buf=synfile,
            sep=',',
            index=False)

    # Compute statistics
    compression_ratio = total_input_variants / total_degenerate_oligos if total_degenerate_oligos > 0 else 1.0
    degeneracies = synthesis_df['Degeneracy'].values if len(synthesis_df) > 0 else [1]
    min_degeneracy = int(min(degeneracies))
    max_degeneracy = int(max(degeneracies))
    mean_degeneracy = float(np.mean(degeneracies))

    # Build stats dictionary
    stats['status'] = True
    stats['basis'] = 'solved'
    stats['vars'] = {
        'input_variants': total_input_variants,
        'degenerate_oligos': total_degenerate_oligos,
        'compression_ratio': round(compression_ratio, 2),
        'length_groups': {k: len(set(v['sequences'])) for k, v in length_groups.items()},
        'min_degeneracy': min_degeneracy,
        'max_degeneracy': max_degeneracy,
        'mean_degeneracy': round(mean_degeneracy, 2),
    }

    # Compression Statistics
    liner.send('\n[Compression Statistics]\n')

    liner.send(' Compression Status    : Successful\n')
    liner.send('       Input Variant(s): {:,}\n'.format(total_input_variants))
    liner.send('  Degenerate Oligo(s)  : {:,}\n'.format(total_degenerate_oligos))
    liner.send(' Compression Ratio     : {:.1f}x\n'.format(compression_ratio))
    liner.send('         Min Degeneracy: {:,} Variant(s) per Oligo\n'.format(min_degeneracy))
    liner.send('         Max Degeneracy: {:,} Variant(s) per Oligo\n'.format(max_degeneracy))
    liner.send('        Mean Degeneracy: {:.1f} Variant(s) per Oligo\n'.format(mean_degeneracy))
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Unschedule file deletion
    if (mapfile is not None) or (synfile is not None):
        ae.unregister(ut.remove_file)

    # Close Liner
    liner.close()

    # Stamp stats
    stats = ut.stamp_stats(
        stats=stats,
        module='compress',
        input_rows=input_rows,
        output_rows=len(synthesis_df))
    stats['random_seed'] = random_seed

    mapping_df_return = mapping_df
    if (mapping_df is not None) and id_from_index:
        mapping_df_return = mapping_df.set_index('ID', drop=True)

    return (mapping_df_return, synthesis_df, stats)
