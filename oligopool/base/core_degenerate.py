import itertools

import numpy as np

from . import utils as ut


# Bitmask Encoding Constants

# Bitmask encoding for bases: A=1, C=2, G=4, T=8
BASE_TO_BIT = {'A': 1, 'C': 2, 'G': 4, 'T': 8}
BIT_TO_BASE = {1: 'A', 2: 'C', 4: 'G', 8: 'T'}

# Build MASK_TO_IUPAC and IUPAC_TO_MASK from utils.ddna_space
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

# Deterministic ordering used for tie-breaking and action lists.
# Higher degeneracy first (N > [3-way] > [2-way] > [single])
IUPAC_ORDER = ['N', 'B', 'D', 'H', 'V', 'R', 'Y', 'S', 'W', 'K', 'M', 'A', 'C', 'G', 'T']

# Pre-computed popcount lookup table for 4-bit masks [0..15]
_POPCOUNT_TABLE = (0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4)


# IUPAC Utility Functions

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
    Return total degeneracy = product of per-position popcounts.
    Internal use only.

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
    Validate that sequence contains only A/C/G/T characters.
    Internal use only.

    :: sequence
       type - string
       desc - DNA sequence to validate
    '''

    return all(ch in BASE_TO_BIT for ch in sequence.strip().upper())


def get_validated_iupac_sequence(sequence):
    '''
    Validate that sequence contains only valid IUPAC codes.
    Internal use only.

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


# Pool Representation Functions

def get_encoded_pool(sequences):
    '''
    Encode pool of concrete A/C/G/T sequences into uint8 array of base indices 0..3.
    Internal use only.

    :: sequences
       type - list
       desc - list of concrete DNA sequences (all same length)

    Returns shape (n, L), dtype uint8, values: A->0, C->1, G->2, T->3
    '''

    if not sequences:
        raise ValueError('Empty pool.')

    length = len(sequences[0])
    base_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    bases = np.empty((len(sequences), length), dtype=np.uint8)

    for row_index, seq in enumerate(sequences):
        seq = seq.strip().upper()
        if len(seq) != length:
            raise ValueError('All sequences in a group must have the same length.')
        bases[row_index, :] = [base_to_index[ch] for ch in seq]

    return bases


def get_pool_view(bases, active):
    '''
    Create a pool view dictionary with bases and active mask.
    Internal use only.

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
    For a single position, return boolean mask over all sequences
    indicating compatibility with the given iupac_mask at that position.
    Internal use only.

    :: pool
       type - dict
       desc - pool view dictionary with bases and active arrays
    :: position
       type - integer
       desc - position index in sequence
    :: iupac_mask
       type - integer
       desc - IUPAC bitmask to check compatibility against
    '''

    base_indices = pool['bases'][:, position]
    return ((iupac_mask >> base_indices) & 1).astype(bool)


def get_prefix_validity(prefix_masks, pool, prefix_compatible):
    '''
    Validity rule: degeneracy(prefix) <= number of active sequences compatible with prefix.
    Internal use only.

    :: prefix_masks
       type - list
       desc - list of bitmasks for current prefix
    :: pool
       type - dict
       desc - pool view dictionary
    :: prefix_compatible
       type - np.ndarray
       desc - boolean mask of sequences compatible with prefix
    '''

    deg = get_degeneracy_product(prefix_masks)
    count = int(np.logical_and(pool['active'], prefix_compatible).sum())
    return deg <= count


def get_action_space(pool, position, prefix_compatible):
    '''
    Allowed IUPAC codes at this position are those whose mask is a subset
    of the bases observed among currently active+compatible targets at this position.
    Internal use only.

    :: pool
       type - dict
       desc - pool view dictionary
    :: position
       type - integer
       desc - position index in sequence
    :: prefix_compatible
       type - np.ndarray
       desc - boolean mask of sequences compatible with current prefix
    '''

    active_and_compatible = np.logical_and(pool['active'], prefix_compatible)
    if not active_and_compatible.any():
        return []

    observed_base_indices = pool['bases'][active_and_compatible, position]
    observed_mask = int(np.bitwise_or.reduce(1 << observed_base_indices))

    allowed = []
    for code in IUPAC_ORDER:
        mask = IUPAC_TO_MASK[code]
        if (mask & ~observed_mask) == 0:
            allowed.append(mask)
    return allowed


# Rollout Construction Functions

def get_random_completion(
    rng,
    pool,
    prefix_masks,
    prefix_compatible,
    horizon_positions):
    '''
    Randomly complete from the current prefix over a limited horizon, keeping validity.
    Returns completed masks (may be shorter than full length if blocked).
    Internal use only.

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
       desc - boolean mask of sequences compatible with prefix
    :: horizon_positions
       type - list
       desc - list of position indices to explore
    '''

    masks = list(prefix_masks)
    compatible_mask = prefix_compatible.copy()

    for position in horizon_positions:
        action_space = get_action_space(pool, position, compatible_mask)
        if not action_space:
            break

        shuffled = action_space.copy()
        rng.shuffle(shuffled)

        chosen_compatible = None
        chosen_mask = None

        for candidate_mask in shuffled:
            position_compat = get_compatible_mask_for_position(pool, position, candidate_mask)
            next_compatible = np.logical_and(compatible_mask, position_compat)

            masks.append(candidate_mask)
            if get_prefix_validity(masks, pool, next_compatible):
                chosen_mask = candidate_mask
                chosen_compatible = next_compatible
                break
            masks.pop()

        if chosen_mask is None or chosen_compatible is None:
            break

        compatible_mask = chosen_compatible

    return masks


def get_expected_reward(
    rng,
    pool,
    prefix_masks,
    prefix_compatible,
    start_position,
    horizon,
    nsims):
    '''
    Monte Carlo rollout estimate: reward = degeneracy(completed_sequence) (non-log).
    Internal use only.

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
       desc - boolean mask of sequences compatible with prefix
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
    horizon_positions = list(range(start_position, end_position))

    if not horizon_positions:
        return float(get_degeneracy_product(prefix_masks))

    rewards_sum = 0.0
    for _ in range(nsims):
        sim_rng = np.random.default_rng(rng.integers(0, 2**63 - 1, dtype=np.int64))
        completed_masks = get_random_completion(
            rng=sim_rng,
            pool=pool,
            prefix_masks=prefix_masks,
            prefix_compatible=prefix_compatible,
            horizon_positions=horizon_positions)
        rewards_sum += float(get_degeneracy_product(completed_masks))

    return float(rewards_sum / nsims)


# Compression Functions

def get_degenerate_oligo(
    rng,
    pool,
    nsims,
    horizon):
    '''
    Construct one degenerate oligo by greedy selection with rollout lookahead.
    Returns (masks, compatible_mask_for_full_oligo).
    Internal use only.

    :: rng
       type - np.random.Generator
       desc - random number generator
    :: pool
       type - dict
       desc - pool view dictionary
    :: nsims
       type - integer
       desc - number of Monte Carlo simulations per decision
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
            position_compat = get_compatible_mask_for_position(pool, position, candidate_mask)
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
            if (score > best_score) or (np.isclose(score, best_score) and candidate_deg > best_base_degeneracy):
                best_score = score
                best_mask = candidate_mask
                best_next_compatible = next_compatible
                best_base_degeneracy = candidate_deg

        if best_mask is None or best_next_compatible is None:
            # Fallback: try low-degeneracy first to see if anything is feasible
            fallback_space = sorted(action_space, key=lambda m: get_popcount4(m))
            extended = False
            for candidate_mask in fallback_space:
                position_compat = get_compatible_mask_for_position(pool, position, candidate_mask)
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
    Compress a concrete pool of a single length into a set of IUPAC-degenerate oligos.
    Internal use only.

    :: sequences
       type - list
       desc - list of concrete A/C/G/T sequences (all same length)
    :: ids
       type - list
       desc - list of corresponding variant IDs
    :: nsims
       type - integer
       desc - number of Monte Carlo rollouts per decision
    :: horizon
       type - integer
       desc - lookahead positions
    :: rng
       type - np.random.Generator
       desc - random number generator
    :: liner
       type - coroutine / None
       desc - optional liner coroutine for progress output
    :: group_length
       type - integer / None
       desc - optional length to display in progress

    Returns tuple of:
        - list of (degenerate_seq, covered_ids, degeneracy)
        - total variants covered
    '''

    # Book-keeping
    unique_map = {}  # seq -> id
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

        # Check if we got a full-length degenerate oligo
        # If not (algorithm couldn't extend), fall back to literal sequence
        target_length = pool['length']
        if len(masks) < target_length:
            # Force progress: take first remaining sequence literally
            remaining_index = int(np.flatnonzero(pool['active'])[0])
            seq_indices = pool['bases'][remaining_index, :]
            masks = [1 << int(idx) for idx in seq_indices]

        degenerate = get_iupac_from_masks(masks)

        # Compute which active sequences are covered by this degenerate oligo
        covered = pool['active'].copy()
        for position, mask in enumerate(masks):
            covered &= get_compatible_mask_for_position(pool, position, int(mask))

        newly_covered = int(covered.sum())
        if newly_covered == 0:
            # Force progress: take first remaining sequence literally
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
    Enumerate all concrete A/C/G/T sequences from an IUPAC-degenerate string.
    Internal use only.

    :: sequence
       type - string
       desc - IUPAC-degenerate sequence
    '''

    sequence = sequence.strip().upper()
    if not get_validated_iupac_sequence(sequence):
        raise ValueError('Invalid IUPAC sequence: {}'.format(sequence))

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
    Compute degeneracy (number of concrete sequences) without expanding.
    Internal use only.

    :: sequence
       type - string
       desc - IUPAC-degenerate sequence
    '''

    sequence = sequence.strip().upper()
    if not get_validated_iupac_sequence(sequence):
        raise ValueError('Invalid IUPAC sequence: {}'.format(sequence))

    return get_degeneracy_product(get_masks_from_iupac(sequence))
