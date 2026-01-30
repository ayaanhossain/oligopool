# Analysis Pipeline Example

This example demonstrates NGS-based activity quantification using `oligopool`'s Analysis Mode.

## Files

- `analysis_pipeline.py` - Main analysis script
- `ribozyme_architecture.csv` - Library architecture with barcodes and variants
- `ribozyme_1M_R1.fq.gz` - Simulated paired-end reads (R1)
- `ribozyme_1M_R2.fq.gz` - Simulated paired-end reads (R2)
- `read_simulator.ipynb` - Notebook showing how the test reads were generated

## Workflow

The analysis pipeline follows three steps:

```
index → pack → xcount
```

1. **Index**: Build barcode/associate indexes for mapping
2. **Pack**: Filter, merge, and deduplicate FASTQ reads
3. **Count**: Combinatorial counting across multiple barcode sets

## Usage

```bash
cd examples/analysis-pipeline
python analysis_pipeline.py
```

## Pipeline Details

### Step 1: Index Barcodes

```python
# Index first barcode (BC1)
op.index(
    barcode_data='ribozyme_architecture.csv',
    barcode_column='BC1',
    barcode_prefix_column='OrangeForwardPrimer',
    index_file='BC1',
)

# Index second barcode with associated variants (BC2)
op.index(
    barcode_data='ribozyme_architecture.csv',
    barcode_column='BC2',
    barcode_prefix_column='PinkForwardPrimer',
    barcode_suffix_column='YellowReversePrimer',
    associate_data='ribozyme_architecture.csv',
    associate_column='Variant',
    index_file='BC2',
)
```

### Step 2: Pack Reads

```python
op.pack(
    r1_fastq_file='ribozyme_1M_R1.fq.gz',
    r2_fastq_file='ribozyme_1M_R2.fq.gz',
    r1_read_type=0,  # Forward
    r2_read_type=1,  # Reverse
    minimum_r1_read_quality=30,
    minimum_r2_read_quality=30,
    pack_file='ribozyme_1M',
)
```

### Step 3: Combinatorial Counting

```python
op.xcount(
    pack_file='ribozyme_1M.oligopool.pack',
    index_files=['BC1.oligopool.index', 'BC2.oligopool.index'],
    count_file='ribozyme_counts',
)
```

## Output

The pipeline produces:
- `*.oligopool.index` - Barcode index files
- `*.oligopool.pack` - Packed read file
- `*.oligopool.xcount.csv` - Count matrix with barcode combinations

## More Information

- See [`OligopoolCalculatorInAction.ipynb`](../OligopoolCalculatorInAction.ipynb) for detailed explanations
- [Analysis Mode documentation](../../docs/docs.md#analysis-mode)
