# Notebook Guide

This directory groups the notebooks by intent so the repository root stays focused on the package, documentation, and result artifacts.

## Tutorials

- [tutorials/class_usability.ipynb](tutorials/class_usability.ipynb): package-first usage examples and tomography demonstrations.
- [tutorials/end-to-end.ipynb](tutorials/end-to-end.ipynb): custom end-to-end walkthrough from preprocessing to reconstruction.
- [tutorials/preprocessing_pe_tomography.ipynb](tutorials/preprocessing_pe_tomography.ipynb): step-by-step phase-estimation and tomography construction.
- [tutorials/preprocessing_pe_tomography_reconstruction.ipynb](tutorials/preprocessing_pe_tomography_reconstruction.ipynb): reconstruction-focused continuation of the preprocessing walkthrough.
- [tutorials/resolution_choices_example.ipynb](tutorials/resolution_choices_example.ipynb): how phase-estimation resolution affects the recovered spectrum.

## Benchmarks

- [benchmarks/Benchmark.ipynb](benchmarks/Benchmark.ipynb): spectral reconstruction benchmarks on synthetic matrices.
- [benchmarks/Benchmark-CSV_Primitives.ipynb](benchmarks/Benchmark-CSV_Primitives.ipynb): benchmark pipeline that writes threshold sweeps to [../benchmark/threshold](../benchmark/threshold).
- [benchmarks/real_hw_experiments_primitives.ipynb](benchmarks/real_hw_experiments_primitives.ipynb): primitive-based runtime and hardware experiments.

## Explorations

- [explorations/new_tomography.ipynb](explorations/new_tomography.ipynb): exploratory tomography notebook updated to the current package imports.
- [explorations/initial_tests.ipynb](explorations/initial_tests.ipynb): early exploratory work retained as historical context.
- [explorations/thetas_computation.ipynb](explorations/thetas_computation.ipynb): helper notebook for qRAM angle construction experiments.

## Notes

- Tutorials and benchmarks are the maintained entry points.
- Exploratory notebooks are kept for provenance and idea recovery; they may require more manual inspection than the curated notebooks.
- Each notebook contains a small bootstrap cell so it can import [../QPCA](../QPCA) from its folder without requiring path tweaks.