# Library Compressor Example

This example demonstrates `oligopool`'s Degenerate Mode for compressing variant libraries into IUPAC-degenerate oligos.

## Files

- `run_degenerate_demo.py` - Main demonstration script
- `mutant_generator.py` - Helper utilities for generating test variants

## Concept

Degenerate Mode compresses similar sequences into fewer IUPAC-degenerate oligos, reducing synthesis costs while maintaining the full variant library.

**Before compression** (64 variants for codon saturation):
```
ATGC...AAA...GCTA  (variant 1)
ATGC...AAC...GCTA  (variant 2)
ATGC...AAG...GCTA  (variant 3)
...
ATGC...TTT...GCTA  (variant 64)
```

**After compression** (1 degenerate oligo):
```
ATGC...NNN...GCTA  (covers all 64 variants)
```

## Usage

```bash
cd examples/library-compressor
python run_degenerate_demo.py
```

## Demo Output

The script runs two demonstrations:

### Demo 1: Codon Saturation
- Generates 64 variants (all codons at one position)
- Compresses to a single NNN degenerate oligo
- **Compression ratio: 64:1**

### Demo 2: Multi-position Variants
- Generates variants with mutations at multiple positions
- Compresses to fewer degenerate oligos
- Verifies lossless roundtrip via `expand()`

## Workflow

```python
import oligopool as op

# Step 1: Compress variants to degenerate oligos
mapping_df, synthesis_df, stats = op.compress(
    input_data=variants_df,
    rollout_simulations=50,
    rollout_horizon=4,
    verbose=True,
)

# Step 2: Expand back to verify (optional)
expanded_df, expand_stats = op.expand(
    input_data=synthesis_df,
    mapping_file=mapping_df,  # Restore original IDs
)
```

## When to Use Degenerate Mode

Ideal for:
- **Codon saturation mutagenesis** (NNN at specific positions)
- **ML-generated libraries** with localized variation
- **Alanine scanning** or similar systematic mutations
- Any library where variants share most of their sequence

Not ideal for:
- Highly diverse libraries (random sequences)
- Libraries where variants differ throughout

## More Information

- [Degenerate Mode documentation](../../docs/docs.md#degenerate-mode)
- [compress API](../../docs/api.md#compress)
- [expand API](../../docs/api.md#expand)
