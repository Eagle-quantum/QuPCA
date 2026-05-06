# Security Notes

## Supported Dependency Baseline

This repository now targets `qiskit>=1.4.2` in [pyproject.toml](pyproject.toml), which avoids the vulnerable `qiskit-terra` releases up to `0.46.3` affected by unsafe QPY deserialization.

## QPY Deserialization

- Normal QuPCA workflows in this repository do not require `qiskit.qpy.load()`.
- Do not deserialize untrusted QPY files.
- If you need to work with serialized circuits, keep Qiskit updated and validate the source of the file before loading it.

## Repository Scope

The package code in [QPCA](QPCA) and the notebooks in the repository focus on circuit construction, tomography, and benchmarking. They do not rely on a QPY import path for normal usage.
