# Design Assembly Parser Example

This example demonstrates a declarative approach to designing complex oligo architectures using a specification-based parser.

## Files

- `design_assembly_parser.py` - Parser library for declarative design specifications
- `run_design_assembly_parser.py` - Example usage with a promoter library
- `promoters.txt` - Input promoter sequences (~4,500 variants)

## Concept

Instead of calling individual `oligopool` functions sequentially, you define your entire oligo architecture as a specification dictionary:

```python
output = dp.design_parser(
    pool_size=len(promoter_list),
    element_names=['Primer1', 'Cut1', 'Promoter', 'Barcode', 'Primer2', ...],
    elements_spec={
        'Primer1': {
            'element_type': 'primer',
            'oligo_length_limit': 250,
            'primer_sequence_constraint': 'N' * 22,
            ...
        },
        'Cut1': {
            'element_type': 'motif',
            'motif_sequence_constraint': 'NNNGGATCCNNN',
            ...
        },
        'Promoter': {
            'element_type': 'variant',
            'sequences': promoter_list,
        },
        'Barcode': {
            'element_type': 'barcode',
            'barcode_length': 16,
            ...
        },
        ...
    }
)
```

## Usage

```bash
cd examples/design-assembly-parser
python run_design_assembly_parser.py
```

## Architecture

The example designs oligos with this structure:

```
[Primer1]-[Cut1]-[Promoter]-[Barcode]-[Primer2]-[Cut2]-[Primer3]-[Filler]
```

Where:
- `Primer1/2/3` - PCR primers with Tm constraints
- `Cut1/2` - Restriction sites (BamHI, XbaI)
- `Promoter` - Core variant sequences
- `Barcode` - Unique identifiers (16 bp, Hamming ≥ 3)
- `Filler` - Spacer to reach target length

## Element Types

The parser supports these element types:

| Type | Description | Maps to |
|------|-------------|---------|
| `variant` | Core sequences (your variants) | Direct input |
| `primer` | Thermodynamic primers | `op.primer()` |
| `barcode` | Hamming-distance barcodes | `op.barcode()` |
| `motif` | Sequence motifs/anchors | `op.motif()` |
| `spacer` | Neutral fill sequences | `op.spacer()` |

## More Information

- [Design Mode documentation](../../docs/docs.md#design-mode)
- [API Reference](../../docs/api.md)
