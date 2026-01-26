'''
Degenerate Mode Demonstration: Compressing Variant Libraries.

This script demonstrates oligopool's Degenerate Mode by:
1. Generating codon variants at specific positions (compresses well)
2. Compressing them into IUPAC-degenerate oligos
3. Verifying lossless roundtrip via expand

Degenerate Mode is ideal for:
- Codon saturation mutagenesis (NNN at a position)
- ML-generated libraries with localized variation
- Any library where variants share most of their sequence

Usage:
    python run_degenerate_demo.py

Author: Ayaan Hossain
'''

import os
import sys

# Add parent to path for oligopool import
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import oligopool as op
from mutant_generator import (
    generate_random_sequence,
    generate_variant_dataframe,
)


def demo_codon_saturation():
    '''Demonstrate compression of codon saturation library.'''

    print('=' * 70)
    print('Demo 1: Codon Saturation Mutagenesis')
    print('=' * 70)

    # Configuration
    SEQUENCE_LENGTH = 150
    CODON_POSITION = 75  # Middle of the sequence
    RANDOM_SEED = 42

    # Step 1: Generate wildtype sequence
    print('\n[Step 1] Generating wildtype sequence')
    print('-' * 40)

    wildtype = generate_random_sequence(SEQUENCE_LENGTH, seed=RANDOM_SEED)
    print(f'Wildtype length: {SEQUENCE_LENGTH} bp')
    print(f'Codon position: {CODON_POSITION}-{CODON_POSITION+2} (0-indexed)')
    print(f'Wildtype codon: {wildtype[CODON_POSITION:CODON_POSITION+3]}')

    # Step 2: Generate all 64 codon variants
    print('\n[Step 2] Generating codon saturation library')
    print('-' * 40)

    df = generate_variant_dataframe(
        wildtype,
        strategy='codon',
        codon_position=CODON_POSITION,
        include_wildtype=True,
    )
    print(f'Generated {len(df)} variants (all 64 codons)')
    print(f'\nDataFrame preview:')
    print(df.head(5).to_string(index=False))
    print('...')

    # Step 3: Compress
    print('\n[Step 3] Compressing with oligopool.compress()')
    print('-' * 40)

    mapping_df, synthesis_df, stats = op.compress(
        input_data=df,
        rollout_simulations=50,
        rollout_horizon=4,
        verbose=True,
        random_seed=RANDOM_SEED,
    )

    # Step 4: Results
    print('\n[Step 4] Compression Results')
    print('-' * 40)

    input_count = stats['vars']['input_variants']
    output_count = stats['vars']['degenerate_oligos']
    ratio = stats['vars']['compression_ratio']

    print(f'Input variants:      {input_count}')
    print(f'Degenerate oligos:   {output_count}')
    print(f'Compression ratio:   {ratio:.1f}x')

    print(f'\nSynthesis DataFrame:')
    print(synthesis_df.to_string(index=False))

    # Verify roundtrip
    print('\n[Step 5] Verifying lossless compression')
    print('-' * 40)

    expand_input = synthesis_df.rename(columns={'DegenerateID': 'ID'})
    expanded_df, _ = op.expand(
        input_data=expand_input,
        sequence_column='DegenerateSeq',
        verbose=False,
    )

    original_seqs = set(df['Sequence'])
    expanded_seqs = set(expanded_df['ExpandedSeq'])

    if original_seqs == expanded_seqs:
        print('VERIFIED: All original sequences recovered exactly!')
    elif original_seqs.issubset(expanded_seqs):
        print(f'All {len(original_seqs)} original sequences recovered.')
        print(f'Degenerate oligo expands to {len(expanded_seqs)} sequences.')

    return ratio


def demo_multi_codon():
    '''Demonstrate compression of multi-codon library.'''

    print('\n')
    print('=' * 70)
    print('Demo 2: Double Codon Saturation (Two Variable Regions)')
    print('=' * 70)

    # Configuration
    SEQUENCE_LENGTH = 150
    RANDOM_SEED = 42

    # Generate wildtype
    print('\n[Step 1] Setup')
    print('-' * 40)

    wildtype = generate_random_sequence(SEQUENCE_LENGTH, seed=RANDOM_SEED)

    # Vary positions at two separate codons (positions 30-32 and 90-92)
    # This creates 64 * 64 = 4096 variants
    positions = [30, 31, 32, 90, 91, 92]  # Two codons

    print(f'Wildtype length: {SEQUENCE_LENGTH} bp')
    print(f'Variable positions: {positions}')
    print(f'Expected variants: 4^6 = {4**6}')

    # Generate variants
    print('\n[Step 2] Generating double codon library')
    print('-' * 40)

    df = generate_variant_dataframe(
        wildtype,
        strategy='multi',
        positions=positions,
        include_wildtype=True,
    )
    print(f'Generated {len(df)} variants')

    # Compress
    print('\n[Step 3] Compressing')
    print('-' * 40)

    mapping_df, synthesis_df, stats = op.compress(
        input_data=df,
        rollout_simulations=30,  # Lower for speed
        rollout_horizon=3,
        verbose=True,
        random_seed=RANDOM_SEED,
    )

    # Results
    print('\n[Step 4] Results')
    print('-' * 40)

    input_count = stats['vars']['input_variants']
    output_count = stats['vars']['degenerate_oligos']
    ratio = stats['vars']['compression_ratio']

    print(f'Input variants:      {input_count:,}')
    print(f'Degenerate oligos:   {output_count}')
    print(f'Compression ratio:   {ratio:.1f}x')

    print(f'\nSynthesis DataFrame:')
    print(synthesis_df.to_string(index=False))

    return ratio


def demo_comparison():
    '''Show compression comparison table.'''

    print('\n')
    print('=' * 70)
    print('Compression Comparison Summary')
    print('=' * 70)

    print('''
Library Type                    Variants    Degenerate    Compression
--------------------------------------------------------------------
Codon saturation (1 codon)          64             1           64x
Double codon (2 codons)          4,096             1        4,096x
Single mutants (150-mer)           451           451            1x

Key insight: Compression works when variants share most of their sequence.
- Codon saturation: All differ at only 3 positions -> compresses to NNN
- Single mutants: Each differs at a unique position -> no compression

Degenerate Mode is ideal for:
- Saturation mutagenesis libraries
- ML-generated variants (often similar)
- Focused evolution at specific residues
''')


def main():
    print('\n')
    print('#' * 70)
    print('#' + ' ' * 20 + 'DEGENERATE MODE DEMO' + ' ' * 28 + '#')
    print('#' * 70)

    # Demo 1: Codon saturation
    ratio1 = demo_codon_saturation()

    # Demo 2: Multi-codon
    ratio2 = demo_multi_codon()

    # Comparison
    demo_comparison()

    print('=' * 70)
    print('DEMO COMPLETE')
    print('=' * 70)
    print(f'''
Results:
- Codon saturation: {ratio1:.0f}x compression
- Double codon: {ratio2:.0f}x compression

These examples demonstrate the power of Degenerate Mode for libraries
where variants share most of their sequence. Order fewer oligos,
save on synthesis costs!
''')


if __name__ == '__main__':
    main()
