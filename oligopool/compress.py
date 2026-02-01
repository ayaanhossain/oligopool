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
    Compress concrete DNA sequences into IUPAC-degenerate oligos.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligo pool variants.

    Optional Parameters:
        - `mapping_file` (`str`): Filename for output mapping DataFrame (default: `None`).
            A `.oligopool.compress.csv` suffix is added if missing.
        - `synthesis_file` (`str`): Filename for output synthesis DataFrame (default: `None`).
            A `.oligopool.compress.csv` suffix is added if missing.
        - `rollout_simulations` (`int`): Number of Monte Carlo simulations per decision (default: `100`).
        - `rollout_horizon` (`int`): Number of positions to look ahead during rollouts (default: `4`).
        - `random_seed` (`int` / `None`): Seed for local RNG (default: `None`).
        - `verbose` (`bool`): If `True`, logs progress updates to stdout (default: `True`).

    Returns:
        - A tuple of (mapping_df, synthesis_df) DataFrames; saves to files if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column; all other columns must be non-empty strict-ATGC DNA strings.
        - All non-'ID' columns are concatenated (left-to-right) to form the sequence per variant.
        - Compression is lossless: expanding the degenerate oligos recovers exactly the
            set of unique input sequences (no extras, no missing).
        - Sequences of different lengths are compressed independently by length group.
        - `rollout_simulations` and `rollout_horizon` trade off runtime vs compression ratio;
            higher values give better compression but take longer.
        - `mapping_df` links each original variant `ID` to its `DegenerateID` for traceability;
            `synthesis_df` contains the degenerate oligos ready for ordering.
        - Use `expand` to verify that compression output covers exactly the original sequences.
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

    # Parse oligopool sequences
    (sequences,
     variant_ids,
     id_to_sequence,
     length_groups,
     unique_count) = cd.get_parsed_oligopool_sequences(
        indf=indf,
        dna_cols=dna_cols,
        liner=liner)

    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Step 2: Compress sequences
    liner.send('\n[Step 2: Compressing Sequences]\n')
    stats['step'] = 2
    stats['step_name'] = 'compressing-sequences'

    # Compress all length groups
    (all_results,
     total_input_variants,
     total_degenerate_oligos) = cd.get_all_compressed_groups(
        length_groups=length_groups,
        nsims=nsims,
        horizon=horizon,
        rng=rng,
        liner=liner)

    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Step 3: Build output DataFrames
    liner.send('\n[Step 3: Building Output DataFrames]\n')
    stats['step'] = 3
    stats['step_name'] = 'building-output'

    # Build DataFrames
    (mapping_df,
     synthesis_df) = cd.get_compression_dataframes(
        all_results=all_results,
        id_to_sequence=id_to_sequence,
        liner=liner)

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
    compression_ratio = total_input_variants / total_degenerate_oligos \
        if total_degenerate_oligos > 0 else 1.0
    degeneracies = synthesis_df['Degeneracy'].values \
        if len(synthesis_df) > 0 else [1]
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
        'length_groups': {
            k: len(set(v['sequences']))
            for k, v in length_groups.items()},
        'min_degeneracy': min_degeneracy,
        'max_degeneracy': max_degeneracy,
        'mean_degeneracy': round(mean_degeneracy, 2),
    }

    # Compression Statistics
    liner.send('\n[Compression Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'input_variants',
            'degenerate_oligos')))

    liner.send(' Compression Status    : Successful\n')
    liner.send('       Input Variant(s): {:{},d}\n'.format(
        total_input_variants, plen))
    liner.send('  Degenerate Oligo(s)  : {:{},d}\n'.format(
        total_degenerate_oligos, plen))
    liner.send(' Compression Ratio     : {:.1f}x\n'.format(
        compression_ratio))
    liner.send('         Min Degeneracy: {:,} Variant(s) per Oligo\n'.format(
        min_degeneracy))
    liner.send('         Max Degeneracy: {:,} Variant(s) per Oligo\n'.format(
        max_degeneracy))
    liner.send('        Mean Degeneracy: {:.1f} Variant(s) per Oligo\n'.format(
        mean_degeneracy))
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
