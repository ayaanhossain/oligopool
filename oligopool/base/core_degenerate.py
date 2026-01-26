import itertools

import numba as nb
import numpy as np

from . import utils as ut


# Bitmask Encoding Constants

BASE_TO_BIT = {'A': 1, 'C': 2, 'G': 4, 'T': 8}
BIT_TO_BASE = {1: 'A', 2: 'C', 4: 'G', 8: 'T'}

# Pre-computed Popcount Lookup Table

_POPCOUNT_TABLE = (0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4)
_NB_POPCOUNT_TABLE = np.array(
    [0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4], dtype=np.int64)

# IUPAC Mapping Functions

def _build_iupac_mappings():
    '''
    Build IUPAC code mappings from utils.ddna_space.
    Internal use only.
    '''

    mask_to_iupac = {}
    iupac_to_mask = {}

    for code, bases in ut.ddna_space.items():
        if code == '-':
            continue
        mask = sum(BASE_TO_BIT[b] for b in bases if b in BASE_TO_BIT)
        if mask > 0:
            mask_to_iupac[mask] = code
            iupac_to_mask[code] = mask

    return mask_to_iupac, iupac_to_mask

MASK_TO_IUPAC, IUPAC_TO_MASK = _build_iupac_mappings()

# IUPAC Ordering (Higher Degeneracy First)

IUPAC_ORDER = ['N', 'B', 'D', 'H', 'V', 'R', 'Y', 'S', 'W', 'K', 'M', 'A', 'C', 'G', 'T']
_NB_IUPAC_ORDER = np.array(
    [IUPAC_TO_MASK[code] for code in IUPAC_ORDER], dtype=np.int64)

# Numba Helper Functions

@nb.njit
def _nb_popcount4(mask):
    '''
    Return popcount for 4-bit mask [0..15].
    Internal use only.
    '''

    return _NB_POPCOUNT_TABLE[mask & 0xF]

@nb.njit
def _nb_degeneracy_product(masks, length):
    '''
    Return total degeneracy = product of
    per-position popcounts. Internal use only.

    :: masks
       type - np.ndarray (int64)
       desc - array of bitmasks
    :: length
       type - int
       desc - number of valid positions in masks
    '''

    total = 1
    for i in range(length):
        total *= _NB_POPCOUNT_TABLE[masks[i] & 0xF]
    return total

@nb.njit
def _nb_compatible_mask_for_position(bases, position, iupac_mask):
    '''
    Return boolean mask indicating which sequences
    are compatible with iupac_mask at position.
    Internal use only.

    :: bases
       type - np.ndarray (uint8, n x L)
       desc - encoded pool bases
    :: position
       type - int
       desc - position index in sequence
    :: iupac_mask
       type - int
       desc - IUPAC bitmask to check
    '''

    n = bases.shape[0]
    result = np.empty(n, dtype=np.bool_)
    for i in range(n):
        base_idx = bases[i, position]
        result[i] = ((iupac_mask >> base_idx) & 1) == 1
    return result

@nb.njit
def _nb_prefix_validity(masks, masks_len, bases, active, prefix_compatible):
    '''
    Validity rule: degeneracy(prefix) <= count
    of active sequences compatible with prefix.
    Internal use only.

    :: masks
       type - np.ndarray (int64)
       desc - array of bitmasks for current prefix
    :: masks_len
       type - int
       desc - number of valid positions in masks
    :: bases
       type - np.ndarray (uint8, n x L)
       desc - encoded pool bases
    :: active
       type - np.ndarray (bool, n)
       desc - mask of active sequences
    :: prefix_compatible
       type - np.ndarray (bool, n)
       desc - mask of sequences compatible with prefix
    '''

    deg = _nb_degeneracy_product(masks, masks_len)
    count = 0
    for i in range(active.shape[0]):
        if active[i] and prefix_compatible[i]:
            count += 1
    return deg <= count

@nb.njit
def _nb_get_observed_mask(bases, active, prefix_compatible, position):
    '''
    Get bitmask of bases observed at position
    among active+compatible sequences.
    Internal use only.
    '''

    observed = 0
    for i in range(bases.shape[0]):
        if active[i] and prefix_compatible[i]:
            observed |= (1 << bases[i, position])
    return observed

@nb.njit
def _nb_get_action_space(bases, active, prefix_compatible, position, iupac_order):
    '''
    Get allowed IUPAC masks at this position.
    Returns array of valid masks and count.
    Internal use only.

    :: bases
       type - np.ndarray (uint8, n x L)
       desc - encoded pool bases
    :: active
       type - np.ndarray (bool, n)
       desc - mask of active sequences
    :: prefix_compatible
       type - np.ndarray (bool, n)
       desc - mask of sequences compatible with prefix
    :: position
       type - int
       desc - position index in sequence
    :: iupac_order
       type - np.ndarray (int64)
       desc - ordered IUPAC masks (high degeneracy first)
    '''

    # Check if any active+compatible
    any_active = False
    for i in range(active.shape[0]):
        if active[i] and prefix_compatible[i]:
            any_active = True
            break

    result = np.empty(16, dtype=np.int64)
    count = 0

    if not any_active:
        return result, count

    observed_mask = _nb_get_observed_mask(
        bases, active, prefix_compatible, position)

    for j in range(iupac_order.shape[0]):
        mask = iupac_order[j]
        if (mask & ~observed_mask) == 0:
            result[count] = mask
            count += 1

    return result, count

# Monte Carlo Rollout Functions

@nb.njit
def _nb_random_completion(
    bases,
    active,
    prefix_masks,
    prefix_len,
    prefix_compatible,
    start_position,
    end_position,
    iupac_order,
    rng_state):
    '''
    Randomly complete from the current prefix
    over a limited horizon, keeping validity.
    Internal use only.

    :: bases
       type - np.ndarray (uint8, n x L)
       desc - encoded pool bases
    :: active
       type - np.ndarray (bool, n)
       desc - mask of active sequences
    :: prefix_masks
       type - np.ndarray (int64)
       desc - array of bitmasks for current prefix
    :: prefix_len
       type - int
       desc - number of valid positions in prefix
    :: prefix_compatible
       type - np.ndarray (bool, n)
       desc - mask of sequences compatible with prefix
    :: start_position
       type - int
       desc - starting position for completion
    :: end_position
       type - int
       desc - ending position (exclusive)
    :: iupac_order
       type - np.ndarray (int64)
       desc - ordered IUPAC masks
    :: rng_state
       type - int
       desc - random seed for this simulation
    '''

    n = bases.shape[0]
    masks = prefix_masks.copy()
    masks_len = prefix_len
    compatible = prefix_compatible.copy()

    # Simple LCG for numba-compatible randomness
    rng = np.uint64(rng_state)

    for position in range(start_position, end_position):
        action_space, action_count = _nb_get_action_space(
            bases, active, compatible, position, iupac_order)

        if action_count == 0:
            break

        # Fisher-Yates shuffle of action space
        for i in range(action_count - 1, 0, -1):
            rng = rng * np.uint64(6364136223846793005) + np.uint64(1442695040888963407)
            j = int(rng % np.uint64(i + 1))
            action_space[i], action_space[j] = action_space[j], action_space[i]

        chosen_mask = -1
        chosen_compatible = np.empty(n, dtype=np.bool_)

        for k in range(action_count):
            candidate_mask = action_space[k]
            position_compat = _nb_compatible_mask_for_position(
                bases, position, candidate_mask)

            next_compatible = np.empty(n, dtype=np.bool_)
            for i in range(n):
                next_compatible[i] = compatible[i] and position_compat[i]

            masks[masks_len] = candidate_mask
            if _nb_prefix_validity(masks, masks_len + 1, bases, active, next_compatible):
                chosen_mask = candidate_mask
                for i in range(n):
                    chosen_compatible[i] = next_compatible[i]
                masks_len += 1
                break

        if chosen_mask == -1:
            break

        compatible = chosen_compatible

    return masks, masks_len

@nb.njit
def _nb_simulate_single_rollout(
    bases,
    active,
    prefix_masks,
    prefix_len,
    prefix_compatible,
    start_position,
    end_position,
    iupac_order,
    rng_state):
    '''
    Run a single Monte Carlo rollout and
    return the degeneracy. Internal use only.
    '''

    completed_masks, completed_len = _nb_random_completion(
        bases, active, prefix_masks, prefix_len, prefix_compatible,
        start_position, end_position, iupac_order, rng_state)
    return _nb_degeneracy_product(completed_masks, completed_len)

@nb.njit(parallel=True)
def _nb_expected_reward_parallel(
    bases,
    active,
    prefix_masks,
    prefix_len,
    prefix_compatible,
    start_position,
    end_position,
    iupac_order,
    nsims,
    base_seed):
    '''
    Monte Carlo rollout estimate with parallel
    simulations. Returns average degeneracy.
    Internal use only.

    :: bases
       type - np.ndarray (uint8, n x L)
       desc - encoded pool bases
    :: active
       type - np.ndarray (bool, n)
       desc - mask of active sequences
    :: prefix_masks
       type - np.ndarray (int64)
       desc - array of bitmasks for current prefix
    :: prefix_len
       type - int
       desc - number of valid positions in prefix
    :: prefix_compatible
       type - np.ndarray (bool, n)
       desc - mask of sequences compatible with prefix
    :: start_position
       type - int
       desc - starting position for rollout
    :: end_position
       type - int
       desc - ending position (exclusive)
    :: iupac_order
       type - np.ndarray (int64)
       desc - ordered IUPAC masks
    :: nsims
       type - int
       desc - number of Monte Carlo simulations
    :: base_seed
       type - int
       desc - base random seed
    '''

    rewards = np.empty(nsims, dtype=np.float64)

    for i in nb.prange(nsims):
        sim_seed = base_seed + i * 2654435761
        rewards[i] = _nb_simulate_single_rollout(
            bases, active, prefix_masks, prefix_len, prefix_compatible,
            start_position, end_position, iupac_order, sim_seed)

    return np.mean(rewards)

# Python Utility Functions

def get_popcount4(mask):
    '''
    Return popcount for 4-bit mask [0..15].
    Internal use only.

    :: mask
       type - integer
       desc - 4-bit bitmask representing bases
    '''

    return _POPCOUNT_TABLE[mask & 0xF]

def get_degeneracy_product(masks):
    '''
    Return total degeneracy = product of
    per-position popcounts. Internal use only.

    :: masks
       type - list / array
       desc - sequence of bitmasks
    '''

    total = 1
    for mask in masks:
        total *= get_popcount4(int(mask))
    return int(total)

def get_validated_concrete_sequence(sequence):
    '''
    Validate that sequence contains only
    A/C/G/T characters. Internal use only.

    :: sequence
       type - string
       desc - DNA sequence to validate
    '''

    return all(ch in BASE_TO_BIT for ch in sequence.strip().upper())

def get_validated_iupac_sequence(sequence):
    '''
    Validate that sequence contains only
    valid IUPAC codes. Internal use only.

    :: sequence
       type - string
       desc - IUPAC sequence to validate
    '''

    return all(ch in IUPAC_TO_MASK for ch in sequence.strip().upper())

def get_iupac_from_masks(masks):
    '''
    Convert bitmasks to IUPAC string.
    Internal use only.

    :: masks
       type - list / array
       desc - sequence of bitmasks
    '''

    return ''.join(MASK_TO_IUPAC[int(mask)] for mask in masks)

def get_masks_from_iupac(sequence):
    '''
    Convert IUPAC string to bitmasks.
    Internal use only.

    :: sequence
       type - string
       desc - IUPAC sequence
    '''

    return [IUPAC_TO_MASK[ch] for ch in sequence.strip().upper()]

# Pool Encoding Functions

def get_encoded_pool(sequences):
    '''
    Encode pool of concrete A/C/G/T sequences
    into uint8 array of base indices 0..3.
    Internal use only.

    :: sequences
       type - list
       desc - list of concrete DNA sequences
              (all same length)

    Returns shape (n, L), dtype uint8,
    values: A->0, C->1, G->2, T->3
    '''

    if not sequences:
        raise ValueError('Empty pool.')

    length = len(sequences[0])
    base_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    bases = np.empty((len(sequences), length), dtype=np.uint8)

    for row_index, seq in enumerate(sequences):
        seq = seq.strip().upper()
        if len(seq) != length:
            raise ValueError(
                'All sequences in a group must have the same length.')
        bases[row_index, :] = [base_to_index[ch] for ch in seq]

    return bases

def get_pool_view(bases, active):
    '''
    Create a pool view dictionary with
    bases and active mask. Internal use only.

    :: bases
       type - np.ndarray
       desc - shape (n, L), uint8 in 0..3
    :: active
       type - np.ndarray
       desc - shape (n,), bool mask of active sequences
    '''

    return {
        'bases': bases,
        'active': active,
        'n_active': int(active.sum()),
        'length': int(bases.shape[1])
    }

# Compatibility and Validity Functions

def get_compatible_mask_for_position(pool, position, iupac_mask):
    '''
    Return boolean mask over all sequences
    indicating compatibility with iupac_mask
    at position. Internal use only.

    :: pool
       type - dict
       desc - pool view dictionary
    :: position
       type - integer
       desc - position index in sequence
    :: iupac_mask
       type - integer
       desc - IUPAC bitmask to check
    '''

    return _nb_compatible_mask_for_position(
        pool['bases'], position, iupac_mask)

def get_prefix_validity(prefix_masks, pool, prefix_compatible):
    '''
    Validity rule: degeneracy(prefix) <= count
    of active sequences compatible with prefix.
    Internal use only.

    :: prefix_masks
       type - list
       desc - list of bitmasks for current prefix
    :: pool
       type - dict
       desc - pool view dictionary
    :: prefix_compatible
       type - np.ndarray
       desc - boolean mask of sequences compatible
    '''

    masks_arr = np.array(prefix_masks, dtype=np.int64)
    return _nb_prefix_validity(
        masks_arr, len(prefix_masks),
        pool['bases'], pool['active'], prefix_compatible)

def get_action_space(pool, position, prefix_compatible):
    '''
    Return allowed IUPAC codes at this position
    (those whose mask is a subset of observed
    bases among active+compatible targets).
    Internal use only.

    :: pool
       type - dict
       desc - pool view dictionary
    :: position
       type - integer
       desc - position index in sequence
    :: prefix_compatible
       type - np.ndarray
       desc - boolean mask of sequences compatible
    '''

    action_arr, count = _nb_get_action_space(
        pool['bases'], pool['active'], prefix_compatible,
        position, _NB_IUPAC_ORDER)
    return list(action_arr[:count])

# Rollout Reward Function

def get_expected_reward(
    rng,
    pool,
    prefix_masks,
    prefix_compatible,
    start_position,
    horizon,
    nsims):
    '''
    Monte Carlo rollout estimate: reward =
    degeneracy(completed_sequence). Uses
    numba-parallel simulations. Internal use only.

    :: rng
       type - np.random.Generator
       desc - random number generator
    :: pool
       type - dict
       desc - pool view dictionary
    :: prefix_masks
       type - list
       desc - list of bitmasks for current prefix
    :: prefix_compatible
       type - np.ndarray
       desc - boolean mask of sequences compatible
    :: start_position
       type - integer
       desc - starting position for rollout
    :: horizon
       type - integer
       desc - number of positions to look ahead
    :: nsims
       type - integer
       desc - number of Monte Carlo simulations
    '''

    if nsims <= 0:
        return float(get_degeneracy_product(prefix_masks))

    length = pool['length']
    end_position = min(length, start_position + horizon)

    if start_position >= end_position:
        return float(get_degeneracy_product(prefix_masks))

    # Prepare arrays for numba
    bases = pool['bases']
    active = pool['active']

    # Pad prefix_masks to max possible length
    max_len = length
    prefix_arr = np.zeros(max_len, dtype=np.int64)
    prefix_len = len(prefix_masks)
    for i, m in enumerate(prefix_masks):
        prefix_arr[i] = m

    # Generate base seed
    base_seed = int(rng.integers(0, 2**31 - 1))

    # Call parallel numba function
    return _nb_expected_reward_parallel(
        bases, active, prefix_arr, prefix_len, prefix_compatible,
        start_position, end_position, _NB_IUPAC_ORDER,
        nsims, base_seed)

# Degenerate Oligo Construction Functions

def get_degenerate_oligo(
    rng,
    pool,
    nsims,
    horizon):
    '''
    Construct one degenerate oligo by greedy
    selection with rollout lookahead.
    Internal use only.

    :: rng
       type - np.random.Generator
       desc - random number generator
    :: pool
       type - dict
       desc - pool view dictionary
    :: nsims
       type - integer
       desc - number of Monte Carlo simulations
    :: horizon
       type - integer
       desc - lookahead positions
    '''

    length = pool['length']
    masks = []
    compatible_mask = np.ones(pool['bases'].shape[0], dtype=bool)

    for position in range(length):
        action_space = get_action_space(pool, position, compatible_mask)
        if not action_space:
            break

        best_score = -1.0
        best_mask = None
        best_next_compatible = None
        best_base_degeneracy = -1

        for candidate_mask in action_space:
            position_compat = get_compatible_mask_for_position(
                pool, position, candidate_mask)
            next_compatible = np.logical_and(compatible_mask, position_compat)

            candidate_prefix = masks + [candidate_mask]
            if not get_prefix_validity(candidate_prefix, pool, next_compatible):
                continue

            score = get_expected_reward(
                rng=rng,
                pool=pool,
                prefix_masks=candidate_prefix,
                prefix_compatible=next_compatible,
                start_position=position + 1,
                horizon=horizon,
                nsims=nsims)

            candidate_deg = get_popcount4(candidate_mask)
            if (score > best_score) or \
               (np.isclose(score, best_score) and candidate_deg > best_base_degeneracy):
                best_score = score
                best_mask = candidate_mask
                best_next_compatible = next_compatible
                best_base_degeneracy = candidate_deg

        if best_mask is None or best_next_compatible is None:
            # Fallback: try low-degeneracy first
            fallback_space = sorted(action_space, key=lambda m: get_popcount4(m))
            extended = False
            for candidate_mask in fallback_space:
                position_compat = get_compatible_mask_for_position(
                    pool, position, candidate_mask)
                next_compatible = np.logical_and(compatible_mask, position_compat)
                candidate_prefix = masks + [candidate_mask]
                if get_prefix_validity(candidate_prefix, pool, next_compatible):
                    best_mask = candidate_mask
                    best_next_compatible = next_compatible
                    extended = True
                    break
            if not extended:
                break

        masks.append(int(best_mask))
        compatible_mask = best_next_compatible

    return masks, compatible_mask

def get_compressed_group(
    sequences,
    ids,
    nsims,
    horizon,
    rng,
    liner=None,
    group_length=None):
    '''
    Compress a concrete pool of a single length
    into a set of IUPAC-degenerate oligos.
    Internal use only.

    :: sequences
       type - list
       desc - list of concrete A/C/G/T sequences
              (all same length)
    :: ids
       type - list
       desc - list of corresponding variant IDs
    :: nsims
       type - integer
       desc - number of Monte Carlo rollouts
    :: horizon
       type - integer
       desc - lookahead positions
    :: rng
       type - np.random.Generator
       desc - random number generator
    :: liner
       type - coroutine / None
       desc - optional liner for progress output
    :: group_length
       type - integer / None
       desc - optional length to display in progress

    Returns tuple of:
        - list of (degenerate_seq, covered_ids, degeneracy)
        - total variants covered
    '''

    # Book-keeping
    unique_map = {}
    for seq, vid in zip(sequences, ids):
        seq_upper = seq.strip().upper()
        if seq_upper not in unique_map:
            unique_map[seq_upper] = vid

    unique_sequences = list(unique_map.keys())
    unique_ids = list(unique_map.values())

    if not unique_sequences:
        return [], 0

    bases = get_encoded_pool(unique_sequences)
    active = np.ones(bases.shape[0], dtype=bool)
    pool = get_pool_view(bases=bases, active=active)

    # Compute print lengths for progress output
    # Max index = len(unique_sequences), max coverage = len(unique_sequences)
    ilen = ut.get_printlen(value=len(unique_sequences))
    clen = ut.get_printlen(value=len(unique_sequences))

    results = []
    total_covered = 0
    degenerate_idx = 0

    while pool['n_active'] > 0:
        degenerate_idx += 1

        masks, _ = get_degenerate_oligo(
            rng=rng,
            pool=pool,
            nsims=nsims,
            horizon=horizon)

        # Check if we got a full-length oligo
        target_length = pool['length']
        if len(masks) < target_length:
            # Fallback: take first remaining sequence literally
            remaining_index = int(np.flatnonzero(pool['active'])[0])
            seq_indices = pool['bases'][remaining_index, :]
            masks = [1 << int(idx) for idx in seq_indices]

        degenerate = get_iupac_from_masks(masks)

        # Compute which active sequences are covered
        covered = pool['active'].copy()
        for position, mask in enumerate(masks):
            covered &= get_compatible_mask_for_position(
                pool, position, int(mask))

        newly_covered = int(covered.sum())
        if newly_covered == 0:
            # Fallback: take first remaining sequence literally
            remaining_index = int(np.flatnonzero(pool['active'])[0])
            seq_indices = pool['bases'][remaining_index, :]
            masks = [1 << int(idx) for idx in seq_indices]
            degenerate = get_iupac_from_masks(masks)
            covered = np.zeros_like(pool['active'])
            covered[remaining_index] = True
            newly_covered = 1

        # Collect covered IDs
        covered_ids = [unique_ids[i] for i in np.flatnonzero(covered)]

        # Mark as inactive
        pool['active'][covered] = False
        pool['n_active'] = int(pool['active'].sum())

        degeneracy = get_degeneracy_product(masks)
        results.append((degenerate, covered_ids, degeneracy))
        total_covered += newly_covered

        # Progress output
        if liner is not None:
            seq_preview = degenerate[:20] + '..' if len(degenerate) > 22 else degenerate
            if group_length:
                liner.send(
                    ' Degenerate {:{},d}: {} Covers {:{},d} Variant(s)\n'.format(
                        degenerate_idx, ilen, seq_preview, newly_covered, clen))
            else:
                liner.send(
                    '|* Degenerate {:{},d}: {} Covers {:{},d} Variant(s)'.format(
                        degenerate_idx, ilen, seq_preview, newly_covered, clen))

    return results, total_covered

# Expansion Functions

def get_expanded_sequences(sequence):
    '''
    Enumerate all concrete A/C/G/T sequences
    from an IUPAC-degenerate string.
    Internal use only.

    :: sequence
       type - string
       desc - IUPAC-degenerate sequence
    '''

    sequence = sequence.strip().upper()
    if not get_validated_iupac_sequence(sequence):
        raise ValueError(
            'Invalid IUPAC sequence: {}'.format(sequence))

    choices_per_position = []
    for ch in sequence:
        mask = IUPAC_TO_MASK[ch]
        bases = []
        for base_bit, base_char in BIT_TO_BASE.items():
            if mask & base_bit:
                bases.append(base_char)
        choices_per_position.append(bases)

    return [''.join(combo) for combo in itertools.product(*choices_per_position)]

def get_expansion_count(sequence):
    '''
    Compute degeneracy (number of concrete
    sequences) without expanding.
    Internal use only.

    :: sequence
       type - string
       desc - IUPAC-degenerate sequence
    '''

    sequence = sequence.strip().upper()
    if not get_validated_iupac_sequence(sequence):
        raise ValueError(
            'Invalid IUPAC sequence: {}'.format(sequence))

    return get_degeneracy_product(get_masks_from_iupac(sequence))

# Orchestration Functions

def get_parsed_oligopool_sequences(
    indf,
    dna_cols,
    liner):
    '''
    Parse oligopool sequences from input DataFrame.
    Concatenate DNA columns and group by length.
    Internal use only.

    :: indf
       type - pd.DataFrame
       desc - input DataFrame with DNA columns
    :: dna_cols
       type - list
       desc - list of DNA column names to concatenate
    :: liner
       type - coroutine
       desc - liner for progress output

    Returns tuple of:
        - sequences: list of concatenated sequences
        - variant_ids: list of variant IDs
        - id_to_sequence: dict mapping ID to sequence
        - length_groups: dict of length -> {sequences, ids}
        - unique_count: number of unique sequences
    '''

    # Concatenate DNA columns
    liner.send(' Concatenating {:,} DNA Column(s) ...'.format(len(dna_cols)))
    sequences = ut.get_df_concat(df=indf[dna_cols])
    variant_ids = list(indf.index)
    id_to_sequence = dict(zip(variant_ids, sequences))

    # Count unique sequences
    unique_sequences = set(sequences)
    unique_count = len(unique_sequences)
    liner.send(' Unique Sequence(s): {:,} Variant(s)\n'.format(unique_count))

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
        liner.send(
            ' Sequence Length(s): {:,}-{:,} Base Pair(s) ({:,} Length Group(s))\n'.format(
                min_len, max_len, len(length_groups)))

    return (sequences, variant_ids, id_to_sequence,
            length_groups, unique_count)

def get_all_compressed_groups(
    length_groups,
    nsims,
    horizon,
    rng,
    liner):
    '''
    Compress all length groups into degenerate oligos.
    Internal use only.

    :: length_groups
       type - dict
       desc - dict of length -> {sequences, ids}
    :: nsims
       type - integer
       desc - number of Monte Carlo rollouts
    :: horizon
       type - integer
       desc - lookahead positions
    :: rng
       type - np.random.Generator
       desc - random number generator
    :: liner
       type - coroutine
       desc - liner for progress output

    Returns tuple of:
        - all_results: list of (deg_seq, covered_ids, degeneracy, length)
        - total_input_variants: count of unique input variants
        - total_degenerate_oligos: count of degenerate oligos
    '''

    all_results = []
    total_input_variants = 0
    total_degenerate_oligos = 0

    # Compute print lengths
    max_group_size = max(
        len(g['sequences']) for g in length_groups.values())
    glen = ut.get_printlen(value=max_group_size)

    for seq_length in sorted(length_groups.keys()):
        group = length_groups[seq_length]
        group_seqs = group['sequences']
        group_ids = group['ids']

        liner.send(' Length Group {:,} bp: {:{},d} Variant(s)\n'.format(
            seq_length, len(group_seqs), glen))

        # Compress this length group
        results, _ = get_compressed_group(
            sequences=group_seqs,
            ids=group_ids,
            nsims=nsims,
            horizon=horizon,
            rng=np.random.default_rng(rng.integers(0, 2**31 - 1)),
            liner=liner,
            group_length=seq_length)

        # Collect results with length info
        for (deg_seq, covered_ids, degeneracy) in results:
            all_results.append((deg_seq, covered_ids, degeneracy, seq_length))

        # Compute compression ratio for this group
        unique_in_group = len(set(group_seqs))
        compression_ratio = unique_in_group / len(results) if results else 1.0

        liner.send(
            ' Group Complete: {:{},d} -> {:{},d} Degenerate Oligo(s) ({:.1f}x)\n'.format(
                unique_in_group, glen, len(results), glen, compression_ratio))

        total_input_variants += unique_in_group
        total_degenerate_oligos += len(results)

    return all_results, total_input_variants, total_degenerate_oligos

def get_compression_dataframes(
    all_results,
    id_to_sequence,
    liner):
    '''
    Build mapping and synthesis DataFrames from
    compression results. Internal use only.

    :: all_results
       type - list
       desc - list of (deg_seq, covered_ids, degeneracy, length)
    :: id_to_sequence
       type - dict
       desc - mapping of variant ID to original sequence
    :: liner
       type - coroutine
       desc - liner for progress output

    Returns tuple of:
        - mapping_df: DataFrame linking variant IDs to degenerate IDs
        - synthesis_df: DataFrame of degenerate oligos for synthesis
    '''

    import pandas as pd

    mapping_rows = []
    synthesis_rows = []
    degenerate_idx = 0

    for (deg_seq, covered_ids, degeneracy, seq_length) in all_results:
        degenerate_idx += 1
        deg_id = 'D{}'.format(degenerate_idx)

        # Build mapping row for each covered variant ID
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

    # Compute print length
    plen = ut.get_printlen(value=max(len(mapping_df), len(synthesis_df)))

    liner.send('   Mapping DataFrame: {:{},d} Row(s)\n'.format(
        len(mapping_df), plen))
    liner.send(' Synthesis DataFrame: {:{},d} Row(s)\n'.format(
        len(synthesis_df), plen))

    return mapping_df, synthesis_df

def get_parsed_degenerate_info(
    indf,
    seqcol,
    liner):
    '''
    Parse degenerate sequences and estimate expansion.
    Internal use only.

    :: indf
       type - pd.DataFrame
       desc - input DataFrame with degenerate sequences
    :: seqcol
       type - string
       desc - column name containing degenerate sequences
    :: liner
       type - coroutine
       desc - liner for progress output

    Returns tuple of:
        - num_input_seqs: count of input degenerate sequences
        - total_estimated: estimated total expansion count
    '''

    num_input_seqs = len(indf)
    liner.send(' Degenerate Oligo(s): {:,}\n'.format(num_input_seqs))

    # Estimate total expansions
    total_estimated = 0
    for seq in indf[seqcol]:
        total_estimated += get_expansion_count(seq)

    liner.send(' Estimated Expansion: {:,} Concrete Sequence(s)\n'.format(
        total_estimated))

    return num_input_seqs, total_estimated

def _expand_single_worker(args):
    '''
    Worker function to expand a single degenerate sequence.
    Module-level for pickling. Internal use only.
    '''

    idx, row_id, deg_seq, base_dict = args
    deg_seq = str(deg_seq)
    degeneracy = get_expansion_count(deg_seq)
    expanded_seqs = get_expanded_sequences(deg_seq)

    rows = []
    for expanded_seq in expanded_seqs:
        row = dict(base_dict)
        row['ID'] = row_id
        row['ExpandedSeq'] = expanded_seq
        row['OligoLength'] = len(expanded_seq)
        rows.append(row)

    return idx, row_id, deg_seq, degeneracy, rows

def get_expanded_dataframe(
    indf,
    seqcol,
    liner):
    '''
    Expand all degenerate sequences in parallel.
    Internal use only.

    :: indf
       type - pd.DataFrame
       desc - input DataFrame with degenerate sequences
    :: seqcol
       type - string
       desc - column name containing degenerate sequences
    :: liner
       type - coroutine
       desc - liner for progress output

    Returns tuple of:
        - expanded_rows: list of dicts for output DataFrame
        - total_expanded: total number of expanded sequences
    '''

    from concurrent.futures import ProcessPoolExecutor
    from concurrent.futures import as_completed

    num_input_seqs = len(indf)

    # Estimate total for print length
    total_estimated = sum(
        get_expansion_count(seq) for seq in indf[seqcol])

    # Compute print lengths
    ilen = ut.get_printlen(value=num_input_seqs)
    dlen = ut.get_printlen(value=total_estimated)

    # Prepare work items
    work_items = []
    for idx, (row_id, deg_seq) in enumerate(indf[seqcol].items(), start=1):
        base_dict = indf.loc[row_id].to_dict()
        work_items.append((idx, row_id, deg_seq, base_dict))

    # Expand sequences
    expanded_rows = []
    total_expanded = 0
    results_by_idx = {}

    # Use parallel expansion for multiple sequences
    if num_input_seqs > 1:
        with ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(_expand_single_worker, item): item[0]
                for item in work_items
            }

            for future in as_completed(futures):
                idx, row_id, deg_seq, degeneracy, rows = future.result()
                results_by_idx[idx] = (row_id, deg_seq, degeneracy, rows)

        # Process results in order
        for idx in sorted(results_by_idx.keys()):
            row_id, deg_seq, degeneracy, rows = results_by_idx[idx]
            expanded_rows.extend(rows)
            total_expanded += len(rows)

            seq_preview = deg_seq[:20] + '..' if len(deg_seq) > 22 else deg_seq
            liner.send(
                ' Expanding {:{},d}: {} (Degeneracy: {:{},d}) -> {:{},d} Sequence(s)\n'.format(
                    idx, ilen, seq_preview, degeneracy, dlen, len(rows), dlen))

    # Single sequence - no parallelization overhead
    else:
        for item in work_items:
            idx, row_id, deg_seq, degeneracy, rows = _expand_single_worker(item)
            expanded_rows.extend(rows)
            total_expanded += len(rows)

            seq_preview = deg_seq[:20] + '..' if len(deg_seq) > 22 else deg_seq
            liner.send(
                ' Expanding {:{},d}: {} (Degeneracy: {:{},d}) -> {:{},d} Sequence(s)\n'.format(
                    idx, ilen, seq_preview, degeneracy, dlen, len(rows), dlen))

    liner.send('\n')
    liner.send(' Expansion Complete: {:,} -> {:,} Concrete Sequence(s)\n'.format(
        num_input_seqs, total_expanded))

    return expanded_rows, total_expanded
