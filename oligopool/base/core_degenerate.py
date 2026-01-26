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
                    ' Degenerate {:>5}: {} Covers {:>5} Variant(s)\n'.format(
                        degenerate_idx, seq_preview, newly_covered))
            else:
                liner.send(
                    '|* Degenerate {:>5}: {} Covers {:>5} Variant(s)'.format(
                        degenerate_idx, seq_preview, newly_covered))

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
