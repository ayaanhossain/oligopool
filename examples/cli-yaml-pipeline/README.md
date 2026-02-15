# CLI Pipeline Examples

This directory contains examples of YAML-based pipeline execution for `oligopool`.

## Files

- `variants.csv` - Sample input data (15 promoter variants)
- `mpra_design_serial.yaml` - Sequential design pipeline example
- `mpra_design_parallel.yaml` - Parallel-branch design example with `join`
- `analysis_single.yaml` - Single-sample analysis DAG example (`index` + `pack` -> `xcount`)
- `analysis_multi.yaml` - Multi-sample analysis DAG (`index once`, `pack/count per sample`)

## Quick Start

```bash
cd examples/cli-yaml-pipeline

# Validate the pipeline config (dry run)
op pipeline --config mpra_design_serial.yaml --dry-run

# Execute the pipeline
op pipeline --config mpra_design_serial.yaml

# Validate parallel-branch design config (dry run)
op pipeline --config mpra_design_parallel.yaml --dry-run

# Validate analysis DAG config
op pipeline --config analysis_single.yaml --dry-run

# Validate multi-sample analysis DAG config
op pipeline --config analysis_multi.yaml --dry-run
```

You can also use the helper runner with explicit design vs analysis targets:

```bash
./run_example.sh design-serial
./run_example.sh design-parallel
./run_example.sh analysis-single
./run_example.sh analysis-multi
./run_example.sh analysis
./run_example.sh analysis-single-run
./run_example.sh analysis-multi-run
./run_example.sh analysis-run
./run_example.sh all-dry
```

## Sequential Pipeline

`mpra_design_serial.yaml` runs four steps in sequence:

```
primer -> barcode -> spacer -> final
```

The sequential example demonstrates basename chaining:
- `output_file: "01_primer"` writes `01_primer.oligopool.primer.csv`
- downstream `input_data: "01_primer"` auto-resolves to that file
- explicit full filenames still work when you want full manual control

Output:
```
Pipeline: MPRA Library Design
Steps: 4

  [1/4] primer ... done
  [2/4] barcode ... done
  [3/4] spacer ... done
  [4/4] final ... done

Pipeline completed successfully.
```

## Parallel-Branch Design Pipeline

`mpra_design_parallel.yaml` shows a rare-but-useful pattern: fan out into independent design branches, then
recombine with `join` into a single annotated table for downstream steps:

```
primer -> { barcode_a, barcode_b } -> join -> spacer -> final
```

## Analysis Pipeline

`analysis_single.yaml` demonstrates a high-value DAG use case for analysis workflows:

```
index_bc1 -------\
index_bc2 --------> count_combinatorial (xcount)
pack_reads -------/
```

Why this is useful:
- `index` and `pack` are independent artifact-build steps and can run in parallel.
- `xcount` runs only after all required artifacts exist.
- Basename chaining stays concise: `index_files: ["a1_bc1", "a1_bc2"]`, `pack_file: "a1_reads"`.

Dry run output:
```
Pipeline: Analysis DAG Demo
Steps: 4 across 2 levels

Dry run: config validation passed.

  Level 1 (parallel):
    index_bc1: index
    index_bc2: index
    pack_reads: pack
  Level 2:
    count_combinatorial: xcount (after: index_bc1, index_bc2, pack_reads)
```

## Multi-Sample Analysis Pipeline

`analysis_multi.yaml` demonstrates the typical production pattern:

```
index_bc1 ----------\
index_bc2 -----------+--> { xcount_sample_a, acount_sample_a }
pack_sample_a -------/

index_bc1 ----------\
index_bc2 -----------+--> { xcount_sample_b, acount_sample_b }
pack_sample_b -------/
```

Why this is useful:
- Build fixed index artifacts once and reuse across many samples.
- Keep sample branches independent so pack/count can scale cleanly.
- Run `xcount` (combinatorial) and `acount` (association verification) in parallel for each sample. (`acount` uses the associate-carrying index only, here `ms_bc2`.)
- Keep each sample output explicit (`ms_s1_xcount`, `ms_s1_acount`, etc) while sharing index basenames.

Dry run output:
```
Pipeline: Multi-Sample Analysis DAG Demo
Steps: 8 across 2 levels

Dry run: config validation passed.

  Level 1 (parallel):
    index_bc1: index
    index_bc2: index
    pack_sample_a: pack
    pack_sample_b: pack
  Level 2 (parallel):
    xcount_sample_a: xcount (after: index_bc1, index_bc2, pack_sample_a)
    acount_sample_a: acount (after: index_bc2, pack_sample_a)
    xcount_sample_b: xcount (after: index_bc1, index_bc2, pack_sample_b)
    acount_sample_b: acount (after: index_bc2, pack_sample_b)
```

## Config Format

### Sequential (simple list)
```yaml
pipeline:
  name: "My Pipeline"
  steps:
    - primer
    - barcode
    - final

primer:
  input_data: "input.csv"
  output_file: "step1"
  # ... parameters

# Include matching 'barcode' and 'final' sections below.
```

### Parallel (DAG with dependencies)
```yaml
pipeline:
  name: "Parallel Pipeline"
  steps:
    - name: index_bc1
      command: index
    - name: index_bc2
      command: index
    - name: pack_reads
      command: pack
    - name: count_combinatorial
      command: xcount
      after: [index_bc1, index_bc2, pack_reads]

index_bc1:
  # ... parameters
pack_reads:
  # ... parameters
```

## Pipeline Facts

- **Output suffixing is append-if-missing**:
  - `output_file: "step1"` -> `step1.oligopool.<module>.csv`
  - `output_file: "step1.csv"` -> `step1.csv.oligopool.<module>.csv`
  - `output_file: "step1.oligopool.<module>.csv"` -> unchanged
- **Basename chaining works for downstream inputs**:
  - if a step writes `output_file: "step1"`, a later step can use `input_data: "step1"`
  - pipeline resolves this to the declared suffixed output path
- **Explicit existing paths always win**:
  - if `input_data` already exists on disk, it is used directly (no alias rewrite)
- **Ambiguous aliases are rejected**:
  - if multiple steps produce the same basename alias, pipeline preflight reports a config error
- **CLI overrides config values**:
  - values passed on the command line take precedence over YAML values
- **`--dry-run` validates and prints execution order only**:
  - sequential: prints step configs
  - DAG: prints dependency levels and parallel groups

## Using --config with Individual Commands

You can also use config files with single commands:

```bash
# Load parameters from config
op barcode --config barcode_config.yaml

# CLI args override config values
op barcode --config barcode_config.yaml --barcode-length 16
```

## More Information

- [User Guide: Config Files](../../docs/docs.md#config-files)
- [CLI Manual](../../docs/docs.md#cli-reference)
- `op manual pipeline` - Pipeline command documentation
