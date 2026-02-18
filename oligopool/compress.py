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
        - `mapping_file` (`str` / `None`): Filename for output mapping DataFrame (default: `None`).
            A `.oligopool.compress.mapping.csv` suffix is added if missing.
        - `synthesis_file` (`str` / `None`): Filename for output synthesis DataFrame (default: `None`).
            A `.oligopool.compress.synthesis.csv` suffix is added if missing.
        - `rollout_simulations` (`int`): Number of Monte Carlo simulations per decision (default: `100`).
        - `rollout_horizon` (`int`): Number of positions to look ahead during rollouts (default: `4`).
        - `random_seed` (`int` / `None`): Seed for local RNG (default: `None`).
        - `verbose` (`bool`): If `True`, logs progress to stdout (default: `True`).

    Returns:
        - A pandas DataFrame `mapping_df` mapping each concrete input to a degenerate oligo:
            `ID`, `Sequence`, `DegenerateID`. Saved to `mapping_file` if specified.
        - A pandas DataFrame `synthesis_df` containing degenerate oligos ready for ordering:
            `DegenerateID`, `DegenerateSeq`, `Degeneracy`, `OligoLength`. Saved to `synthesis_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column; all other columns must be non-empty strict-ATGC DNA strings.
        - All non-'ID' columns are concatenated (left-to-right) to form the sequence per variant.
        - Compression is lossless: expanding the degenerate oligos recovers exactly the
            set of unique input sequences (no extras, no missing).
        - Sequences of different lengths are compressed independently by length group.
        - `rollout_simulations` and `rollout_horizon` trade off runtime vs compression ratio;
            higher values give better compression but take longer.
        - Use `expand` to verify that compression output covers exactly the original sequences.
        - Best suited for selection assays where variants are identified by sequencing enrichment
            (rather than barcode-based per-variant quantification).
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

    _MAPPING_SUFFIX = '.oligopool.compress.mapping.csv'
    _SYNTHESIS_SUFFIX = '.oligopool.compress.synthesis.csv'

    # Validate mapping file
    mapfile_valid = vp.get_outdf_validity(
        outdf=mapfile,
        outdf_suffix=_MAPPING_SUFFIX,
        outdf_field='    Mapping File       ',
        liner=liner)

    # Validate synthesis file
    synfile_valid = vp.get_outdf_validity(
        outdf=synfile,
        outdf_suffix=_SYNTHESIS_SUFFIX,
        outdf_field='  Synthesis File       ',
        liner=liner)

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

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
    (random_seed,
    seed_valid) = vp.get_parsed_random_seed_info(
        random_seed=random_seed,
        random_seed_field='     Random Seed       ',
        liner=liner)

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

    if mapfile is not None:
        mapfile = ut.get_adjusted_path(path=mapfile, suffix=_MAPPING_SUFFIX)
    if synfile is not None:
        synfile = ut.get_adjusted_path(path=synfile, suffix=_SYNTHESIS_SUFFIX)

    if (mapfile is not None) and (synfile is not None) and (mapfile == synfile):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s): mapping_file and synthesis_file resolve to the same output path.')

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
        ae.register(
            ut.remove_file,
            mapfile)

    if synfile is not None:
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
        'input_variants'   : total_input_variants,
        'degenerate_oligos': total_degenerate_oligos,
        'compression_ratio': round(compression_ratio, 2),
        'length_groups'    : {k: len(set(v['sequences'])) for k, v in length_groups.items()},
        'min_degeneracy'   : min_degeneracy,
        'max_degeneracy'   : max_degeneracy,
        'mean_degeneracy'  : round(mean_degeneracy, 2),
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
    stats['random_seed'] = random_seed
    stats = ut.stamp_stats(
        stats=stats,
        module='compress',
        input_rows=input_rows,
        output_rows=len(synthesis_df))

    mapping_df_return = mapping_df
    if (mapping_df is not None) and id_from_index:
        mapping_df_return = mapping_df.set_index('ID', drop=True)

    return (mapping_df_return, synthesis_df, stats)
