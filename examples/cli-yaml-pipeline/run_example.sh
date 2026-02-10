#!/bin/bash
# Run the CLI pipeline examples
# Usage: ./run_example.sh [sequential|parallel|both]

set -e
cd "$(dirname "$0")"

run_sequential() {
    echo "=== Sequential Pipeline ==="
    echo ""
    echo "Dry run:"
    op pipeline --config mpra_pipeline.yaml --dry-run
    echo ""
    echo "Executing:"
    op pipeline --config mpra_pipeline.yaml
    echo ""
}

run_parallel() {
    echo "=== Parallel Pipeline ==="
    echo ""
    echo "Dry run:"
    op pipeline --config parallel_pipeline.yaml --dry-run
    echo ""
    echo "Executing:"
    op pipeline --config parallel_pipeline.yaml
    echo ""
}

clean() {
    echo "Cleaning up output files..."
    rm -f *.oligopool.*.csv
}

case "${1:-both}" in
    sequential)
        clean
        run_sequential
        ;;
    parallel)
        clean
        run_parallel
        ;;
    both)
        clean
        run_sequential
        run_parallel
        ;;
    clean)
        clean
        ;;
    *)
        echo "Usage: $0 [sequential|parallel|both|clean]"
        exit 1
        ;;
esac

echo "Done! Check the output CSV files."
ls -la *.oligopool.*.csv 2>/dev/null || echo "(no output files)"
