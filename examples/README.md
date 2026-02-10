# Oligopool Calculator Examples

This directory contains working examples demonstrating various `oligopool` workflows.

## Examples

| Directory | Description |
|-----------|-------------|
| [`cli-yaml-pipeline/`](cli-yaml-pipeline/) | YAML-based CLI pipelines (sequential and parallel DAG execution) |
| [`design-assembly-parser/`](design-assembly-parser/) | Declarative design parser for complex oligo architectures |
| [`analysis-pipeline/`](analysis-pipeline/) | NGS analysis with indexing, packing, and combinatorial counting |
| [`library-compressor/`](library-compressor/) | Degenerate Mode: compress variant libraries into IUPAC oligos |

## Notebook

- [`OligopoolCalculatorInAction.ipynb`](OligopoolCalculatorInAction.ipynb) - Interactive notebook demonstrating design and analysis workflows

## Quick Start

```bash
# CLI pipeline example
cd cli-yaml-pipeline
op pipeline --config mpra_pipeline.yaml

# Python examples
cd design-assembly-parser
python run_design_assembly_parser.py

cd analysis-pipeline
python run_analysis_pipeline.py

cd library-compressor
python run_degenerate_demo.py
```

## Documentation

- [User Guide](../docs/docs.md) - Comprehensive tutorials and workflows
- [API Reference](../docs/api.md) - Complete parameter documentation
- [AI Agent Guide](../docs/agent-skills.md) - For AI-assisted design
