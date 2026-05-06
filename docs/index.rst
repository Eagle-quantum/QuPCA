.. QPCA documentation master file, created by
   sphinx-quickstart on Thu Mar 30 16:07:55 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to QPCA's documentation!
================================
The QPCA package, implemented using the Qiskit SDK, consists of five modules:

* :mod:`~QPCA.decomposition` module: this module contains the implementation of the :class:`~QPCA.QPCA` class, which is the main class for the Quantum PCA algorithm. The :meth:`~QPCA.QPCA.fit` method is used to fit the QPCA model, which involves building the circuit to encode the input matrix and performing phase estimation to encode the eigenvalues in quantum registers. The :meth:`~QPCA.QPCA.eigenvectors_reconstruction` method allows you to reconstruct the original eigenvectors and eigenvalues using the :class:`~QPCA.StateVectorTomography` class. Additionally, the :meth:`~QPCA.QPCA.quantum_input_matrix_reconstruction` method allows you to reconstruct the original input matrix. If the results are not satisfactory, you can consider increasing the number of shots or the number of qubits for the phase estimation resolution.

* :mod:`~QPCA.quantumUtilities` module: this module contains the main quantum routines used by the QPCA algorithm. The :class:`~QPCA.StateVectorTomography` class is particularly important, as it is used by QPCA for eigenvector reconstruction. However, it can also be used as a standalone class to reconstruct the state vector of an arbitrary quantum circuit composed only of real amplitudes.

* :mod:`~QPCA.preprocessingUtilities` module: this module includes the :meth:`~QPCA.generate_matrix` method, which is useful for generating an arbitrary random Hermitian input matrix to be used in the QPCA algorithm.

* :mod:`~QPCA.postprocessingUtilities` module: this module provides the :meth:`~QPCA.general_postprocessing` method, which is used by the QPCA algorithm to reconstruct the eigenvectors and eigenvalues using the information obtained from quantum state tomography.

* :mod:`~QPCA.benchmark` module: this module contains the :class:`~QPCA.Benchmark_Manager` class, which manages all the benchmarking tasks for the QPCA algorithm. It provides methods to benchmark the execution of QPCA and analyze the accuracy of eigenvector and eigenvalue reconstruction. 

Each module plays a specific role in the overall QPCA package, allowing for an end-to-end process from input matrix encoding to eigenvector and eigenvalue reconstruction using quantum routines.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   modules
   Example_modules
   Example_tomography_modules
   Example_benchmark

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
