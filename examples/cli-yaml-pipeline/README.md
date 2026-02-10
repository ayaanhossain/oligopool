# CLI Pipeline Examples

This directory contains examples of YAML-based pipeline execution for `oligopool`.

## Files

- `variants.csv` - Sample input data (15 promoter variants)
- `mpra_pipeline.yaml` - Sequential pipeline example
- `parallel_pipeline.yaml` - Parallel DAG pipeline example

## Quick Start

```bash
cd examples/cli-yaml-pipeline

# Validate the pipeline config (dry run)
op pipeline --config mpra_pipeline.yaml --dry-run

# Execute the pipeline
op pipeline --config mpra_pipeline.yaml
```

## Sequential Pipeline

`mpra_pipeline.yaml` runs four steps in sequence:

```
primer -> barcode -> spacer -> final
```

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

## Parallel Pipeline

`parallel_pipeline.yaml` demonstrates DAG execution where independent steps run concurrently:

```
fwd_primer ------> add_spacer --> finalize
                                  ^
barcode_design -------------------|
```

Dry run output:
```
Pipeline: Parallel Design Demo
Steps: 4 across 3 levels

Dry run: config validation passed.

  Level 1 (parallel):
    fwd_primer: primer
    barcode_design: barcode
  Level 2:
    add_spacer: spacer (after: fwd_primer)
  Level 3:
    finalize: final (after: barcode_design, add_spacer)
```

Execution output:
```
Pipeline: Parallel Design Demo
Steps: 4 across 3 levels

  Level 1: fwd_primer, barcode_design (parallel) ... done
  Level 2: add_spacer ... done
  Level 3: finalize ... done

Pipeline completed successfully.
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
  output_file: "step1.csv"
  # ... parameters
```

### Parallel (DAG with dependencies)
```yaml
pipeline:
  name: "Parallel Pipeline"
  steps:
    - name: step_a
      command: primer
    - name: step_b
      command: barcode
      # No 'after' - runs parallel with step_a
    - name: step_c
      command: final
      after: [step_a, step_b]  # Waits for both

step_a:
  # ... parameters
step_b:
  # ... parameters
```

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
