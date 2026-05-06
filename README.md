# QuPCA

Repository for the paper "Towards an End-to-End Approach for Quantum Principal Component Analysis".

QuPCA collects the code, notebooks, and experimental artifacts used to study an end-to-end implementation of Quantum Principal Component Analysis. The project focuses on the full workflow: matrix encoding, phase estimation, eigenvalue and eigenvector reconstruction through tomography, and downstream benchmarking.

The implementation extends existing state-tomography ideas to recover eigenvectors from the QPCA output state, including the non-trivial case of non-integer eigenvalues. The repository also includes benchmarking notebooks and experiment outputs for synthetic and finance-inspired use cases.

## Quick Start

Install the package in editable mode from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .[notebooks]
```

Minimal example:

```python
from QPCA import QPCA, generate_matrix

input_matrix = generate_matrix(
	matrix_dimension=4,
	replicate_paper=False,
	seed=7,
	eigenvalues_list=[0.65, 0.25, 0.06, 0.04],
)

qpca = QPCA().fit(input_matrix=input_matrix, resolution=8, optimized_qram=True)
eigenvalues, eigenvectors = qpca.eigenvectors_reconstruction(n_shots=10000)
```

## Security Note

- Do not install vulnerable `qiskit-terra` releases up to `0.46.3`.
- The project metadata now targets `qiskit>=1.4.2`, which includes the upstream fix for the unsafe `qiskit.qpy.load()` deserialization issue.
- Normal QuPCA workflows in this repository do not require QPY deserialization. Do not load untrusted QPY files.

Additional details are documented in [SECURITY.md](SECURITY.md).

## Repository Guide

Start here if you want to understand the repository quickly:

1. [QPCA](QPCA) contains the reusable Python package.
2. [notebooks/README.md](notebooks/README.md) is the curated notebook index.
3. [notebooks/tutorials/class_usability.ipynb](notebooks/tutorials/class_usability.ipynb) is the most direct package-oriented notebook.
4. [notebooks/tutorials/end-to-end.ipynb](notebooks/tutorials/end-to-end.ipynb) shows the full custom workflow from preprocessing to reconstruction.
5. [notebooks/benchmarks/Benchmark.ipynb](notebooks/benchmarks/Benchmark.ipynb) and [notebooks/benchmarks/Benchmark-CSV_Primitives.ipynb](notebooks/benchmarks/Benchmark-CSV_Primitives.ipynb) cover benchmark workflows.

### Core Package

- [QPCA/decomposition](QPCA/decomposition) implements the main `QPCA` class.
- [QPCA/quantumUtilities](QPCA/quantumUtilities) contains qRAM, phase estimation, and tomography helpers.
- [QPCA/preprocessingUtilities](QPCA/preprocessingUtilities) provides matrix generation and padding utilities.
- [QPCA/postprocessingUtilities](QPCA/postprocessingUtilities) reconstructs eigenpairs from tomography outputs.
- [QPCA/benchmark](QPCA/benchmark) contains plotting and evaluation logic.

### Notebooks

- [notebooks/tutorials](notebooks/tutorials): curated walkthroughs and package-oriented examples.
- [notebooks/benchmarks](notebooks/benchmarks): benchmark and runtime evaluation workflows.
- [notebooks/explorations](notebooks/explorations): exploratory historical notebooks kept for reference.
- [notebooks/README.md](notebooks/README.md): notebook index, support level, and execution guidance.

All notebooks live under [notebooks](notebooks). They include a small bootstrap cell so imports continue to work when you open them from the new subfolders, but installing the package in editable mode is still the recommended workflow.

### Data And Results

- [benchmark](benchmark) stores benchmark CSV outputs.
- [pay_as_you_go_results](pay_as_you_go_results) stores experiment outputs for synthetic and real-data runs.
- [documents](documents) stores reference papers.

## Documentation

- Sphinx source lives in [docs](docs).
- A built HTML snapshot is currently committed under [docs/_build/html](docs/_build/html).
- The notebook catalog lives in [notebooks/README.md](notebooks/README.md).

## Installation Notes

- The package metadata is defined in [pyproject.toml](pyproject.toml).
- For package development, use `pip install -e .`.
- For notebooks and examples, use `pip install -e .[notebooks]`.
- For docs generation, use `pip install -e .[docs]`.

## Citation

If you use this repository, cite the associated paper referenced above and the original QPCA literature discussed in [documents](documents).
