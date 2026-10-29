#!/usr/bin/env python3
"""
run_all.py — StepWise-Q Unified Entry Point
Single command to run benchmarks, models, fairness experiments, figures, and validations.

Usage:
    python run_all.py --all       # Regenerate entire pipeline end-to-end
    python run_all.py --smoke     # Quick smoke test pipeline (< 25s)
    python run_all.py --test      # Run unit test suite
    python run_all.py --validate  # Run validation suite
    python run_all.py --fairness  # Run processor fairness study
"""

from src.runner import main

if __name__ == "__main__":
    main()

