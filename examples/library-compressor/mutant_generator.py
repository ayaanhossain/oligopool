'''
Mutant Generator for Degenerate Mode Demonstration.

This module provides functions to generate variant libraries that
demonstrate the compress/expand workflow in oligopool's Degenerate Mode.

Includes:
- Single nucleotide mutants (educational - shows poor compression case)
- Codon variants at a specific position (shows good compression)
- Combinatorial variants at multiple positions (best compression demo)

Author: Ayaan Hossain
'''

import random
import itertools
from typing import Generator, Tuple, List


# Standard DNA bases
DNA_BASES = ('A', 'C', 'G', 'T')

# All 64 codons
CODONS = [''.join(c) for c in itertools.product(DNA_BASES, repeat=3)]


def generate_random_sequence(length: int, seed: int | None = None) -> str:
    '''
    Generate a random DNA sequence of specified length.

    Parameters:
        length: Length of the sequence to generate.
        seed: Random seed for reproducibility (optional).

    Returns:
        A random DNA sequence string.
    '''
    if seed is not None:
        random.seed(seed)
    return ''.join(random.choice(DNA_BASES) for _ in range(length))


def generate_single_mutants(
    sequence: str,
    include_wildtype: bool = True,
) -> Generator[Tuple[str, str, int, str, str], None, None]:
    '''
    Generate all single-nucleotide mutants of a DNA sequence.

    For a sequence of length L, this generates up to 3*L mutants.
    NOTE: These do NOT compress well because mutations are spread
    across the entire sequence.

    Parameters:
        sequence: The wildtype DNA sequence.
        include_wildtype: If True, yield the wildtype as the first entry.

    Yields:
        Tuples of (variant_id, mutant_sequence, position, ref_base, alt_base)
    '''
    sequence = sequence.upper()

    for base in sequence:
        if base not in DNA_BASES:
            raise ValueError(f'Invalid base {base!r}. Use only A, C, G, T.')

    if include_wildtype:
        yield ('WT', sequence, 0, '-', '-')

    seq_list = list(sequence)
    for pos in range(len(sequence)):
        ref_base = sequence[pos]
        for alt_base in DNA_BASES:
            if alt_base != ref_base:
                seq_list[pos] = alt_base
                mutant_seq = ''.join(seq_list)
                seq_list[pos] = ref_base
                variant_id = f'{ref_base}{pos + 1}{alt_base}'
                yield (variant_id, mutant_seq, pos + 1, ref_base, alt_base)


def generate_codon_variants(
    sequence: str,
    codon_position: int,
    include_wildtype: bool = True,
) -> Generator[Tuple[str, str], None, None]:
    '''
    Generate all 64 codon variants at a specific position.

    This is an ideal case for compression - 64 variants that differ
    only at 3 consecutive positions compress to a single degenerate
    oligo with NNN at the codon position.

    Parameters:
        sequence: The wildtype DNA sequence (length must be >= codon_position + 3).
        codon_position: 0-indexed start position of the codon to vary.
        include_wildtype: If True, include wildtype codon.

    Yields:
        Tuples of (variant_id, variant_sequence)

    Example:
        >>> seq = 'ATGATGCCC'  # 9 bp
        >>> for vid, var in generate_codon_variants(seq, 3):
        ...     print(f'{vid}: {var}')
        ATG_AAA: ATGAAACCC
        ATG_AAC: ATGAACCCC
        ...
    '''
    sequence = sequence.upper()

    if len(sequence) < codon_position + 3:
        raise ValueError(
            f'Sequence length {len(sequence)} too short for codon at position {codon_position}'
        )

    wt_codon = sequence[codon_position:codon_position + 3]
    prefix = sequence[:codon_position]
    suffix = sequence[codon_position + 3:]

    for codon in CODONS:
        if not include_wildtype and codon == wt_codon:
            continue
        variant_id = f'{wt_codon}_{codon}'
        variant_seq = prefix + codon + suffix
        yield (variant_id, variant_seq)


def generate_multi_position_variants(
    sequence: str,
    positions: List[int],
    include_wildtype: bool = True,
) -> Generator[Tuple[str, str], None, None]:
    '''
    Generate all combinatorial variants at specified positions.

    Each position can be A, C, G, or T, giving 4^len(positions) variants.
    These compress extremely well - all variants can be represented by
    a single degenerate oligo with 'N' at each varied position.

    Parameters:
        sequence: The wildtype DNA sequence.
        positions: List of 0-indexed positions to vary.
        include_wildtype: If True, include wildtype.

    Yields:
        Tuples of (variant_id, variant_sequence)

    Example:
        >>> seq = 'ATGC'
        >>> list(generate_multi_position_variants(seq, [1, 3]))
        [('pos1A_pos3A', 'AAGA'), ('pos1A_pos3C', 'AAGC'), ...]
    '''
    sequence = sequence.upper()
    seq_list = list(sequence)

    wt_bases = tuple(sequence[p] for p in positions)

    # Generate all combinations
    for bases in itertools.product(DNA_BASES, repeat=len(positions)):
        if not include_wildtype and bases == wt_bases:
            continue

        # Build variant
        for pos, base in zip(positions, bases):
            seq_list[pos] = base

        # Build ID
        variant_id = '_'.join(f'p{p+1}{b}' for p, b in zip(positions, bases))
        variant_seq = ''.join(seq_list)

        yield (variant_id, variant_seq)

        # Restore for next iteration
        for pos, base in zip(positions, wt_bases):
            seq_list[pos] = base


def count_variants(length: int = None, positions: int = None, codons: int = None) -> int:
    '''
    Calculate number of variants for different generation strategies.

    Parameters:
        length: For single mutants, sequence length.
        positions: For multi-position, number of positions.
        codons: For codon variants, number of codons (usually 1).

    Returns:
        Number of variants (excluding wildtype).
    '''
    if length is not None:
        return 3 * length  # Single mutants
    if positions is not None:
        return 4 ** positions  # Multi-position
    if codons is not None:
        return 64 ** codons  # Codon variants
    return 0


def generate_variant_dataframe(
    sequence: str,
    strategy: str = 'codon',
    codon_position: int = 0,
    positions: List[int] = None,
    include_wildtype: bool = True,
):
    '''
    Generate a pandas DataFrame of variants using specified strategy.

    Parameters:
        sequence: The wildtype DNA sequence.
        strategy: One of 'single', 'codon', or 'multi'.
        codon_position: For 'codon' strategy, the start position.
        positions: For 'multi' strategy, list of positions.
        include_wildtype: If True, include wildtype.

    Returns:
        pandas DataFrame with columns: ID, Sequence
    '''
    import pandas as pd

    if strategy == 'single':
        data = [(vid, seq) for vid, seq, *_ in generate_single_mutants(sequence, include_wildtype)]
    elif strategy == 'codon':
        data = list(generate_codon_variants(sequence, codon_position, include_wildtype))
    elif strategy == 'multi':
        if positions is None:
            raise ValueError("positions required for 'multi' strategy")
        data = list(generate_multi_position_variants(sequence, positions, include_wildtype))
    else:
        raise ValueError(f"Unknown strategy: {strategy}")

    return pd.DataFrame(data, columns=['ID', 'Sequence'])


if __name__ == '__main__':
    print('Mutant Generator Demo')
    print('=' * 50)

    # Test sequence
    test_seq = 'ATGCATGCATGC'  # 12 bp

    print(f'\nWildtype: {test_seq}')
    print(f'Length: {len(test_seq)} bp')

    # Codon variants
    print('\n--- Codon Variants (positions 3-5) ---')
    codon_df = generate_variant_dataframe(test_seq, strategy='codon', codon_position=3)
    print(f'Total variants: {len(codon_df)}')
    print(codon_df.head(8).to_string(index=False))

    # Multi-position variants
    print('\n--- Multi-Position Variants (positions 0, 4, 8) ---')
    multi_df = generate_variant_dataframe(test_seq, strategy='multi', positions=[0, 4, 8])
    print(f'Total variants: {len(multi_df)} (4^3 = 64)')
    print(multi_df.head(8).to_string(index=False))

    # Single mutants
    print('\n--- Single Mutants (all positions) ---')
    single_df = generate_variant_dataframe(test_seq, strategy='single')
    print(f'Total variants: {len(single_df)} (1 + 3*12 = 37)')
    print(single_df.head(8).to_string(index=False))
