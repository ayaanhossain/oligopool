#!/bin/bash
# Run the CLI pipeline examples
# Usage: ./run_example.sh [design-serial|design-parallel|design|analysis-single|analysis-multi|analysis|analysis-single-run|analysis-multi-run|analysis-run|all-dry|clean]
# Legacy aliases still supported: sequential, both

set -e
cd "$(dirname "$0")"

run_sequential() {
    echo "=== Design: Serial Pipeline ==="
    echo ""
    echo "Dry run:"
    op pipeline --config mpra_design_serial.yaml --dry-run
    echo ""
    echo "Executing:"
    op pipeline --config mpra_design_serial.yaml
    echo ""
}

run_parallel_design() {
    echo "=== Design: Parallel Branches + Join ==="
    echo ""
    echo "Dry run:"
    op pipeline --config mpra_design_parallel.yaml --dry-run
    echo ""
    echo "Executing:"
    op pipeline --config mpra_design_parallel.yaml
    echo ""
}

run_analysis_single_dry() {
    echo "=== Analysis: Single-Sample DAG (Dry Run) ==="
    echo ""
    op pipeline --config analysis_single.yaml --dry-run
    echo ""
}

run_analysis_single() {
    echo "=== Analysis: Single-Sample DAG ==="
    echo ""
    op pipeline --config analysis_single.yaml
    echo ""
}

run_analysis_multi_dry() {
    echo "=== Analysis: Multi-Sample DAG (Dry Run) ==="
    echo ""
    op pipeline --config analysis_multi.yaml --dry-run
    echo ""
}

run_analysis_multi() {
    echo "=== Analysis: Multi-Sample DAG ==="
    echo ""
    op pipeline --config analysis_multi.yaml
    echo ""
}

run_all_dry() {
    echo "=== Design + Analysis (Dry Run) ==="
    echo ""
    op pipeline --config mpra_design_serial.yaml --dry-run
    echo ""
    op pipeline --config mpra_design_parallel.yaml --dry-run
    echo ""
    op pipeline --config analysis_single.yaml --dry-run
    echo ""
    op pipeline --config analysis_multi.yaml --dry-run
    echo ""
}

clean() {
    echo "Cleaning up output artifacts..."
    rm -rf -- *.oligopool.*
}

case "${1:-design}" in
    design-serial|sequential)
        clean
        run_sequential
        ;;
    design-parallel)
        clean
        run_parallel_design
        ;;
    design|both)
        clean
        run_sequential
        ;;
    analysis-single)
        run_analysis_single_dry
        ;;
    analysis-multi)
        run_analysis_multi_dry
        ;;
    analysis)
        run_analysis_single_dry
        run_analysis_multi_dry
        ;;
    analysis-single-run)
        clean
        run_analysis_single
        ;;
    analysis-multi-run)
        clean
        run_analysis_multi
        ;;
    analysis-run)
        clean
        run_analysis_single
        run_analysis_multi
        ;;
    all-dry)
        run_all_dry
        ;;
    clean)
        clean
        ;;
    *)
        echo "Usage: $0 [design-serial|design-parallel|design|analysis-single|analysis-multi|analysis|analysis-single-run|analysis-multi-run|analysis-run|all-dry|clean]"
        echo "Legacy aliases: sequential, both"
        exit 1
        ;;
esac

echo "Done! Check the output artifacts."
ls -la *.oligopool.* 2>/dev/null || echo "(no output artifacts)"
