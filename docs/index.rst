.. Qpca documentation master file, created by
   sphinx-quickstart on Thu Mar 30 16:07:55 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Qpca's documentation!
================================
The Qpca package, implemented using the Qiskit SDK, consists of five modules:

* :mod:`~QPCA.decomposition` module: this module contains the implementation of the :class:`~QPCA.decomposition.QPCA` class, which is the main class for the Quantum PCA algorithm. The :meth:`~QPCA.decomposition.QPCA.fit` method is used to fit the Qpca model, which involves building the circuit to encode the input matrix and performing phase estimation to encode the eigenvalues in quantum registers. The :meth:`~QPCA.decomposition.QPCA.eigenvectors_reconstruction` method allows you to reconstruct the original eigenvectors (and eigenvalues) using the :class:`~QPCA.quantumUtilities.Tomography.StateVectorTomography` class, which implements state vector tomography algorithm. Additionally, the :meth:`~QPCA.decomposition.QPCA.quantum_input_matrix_reconstruction` method allows you to reconstruct the original input matrix. If the results are not satisfactory, you can consider increasing the number of shots or the number of qubits for the phase estimation resolution.

* :mod:`~QPCA.quantumUtilities` module: this module contains the main quantum routines used by the Qpca algorithm. The :class:`~QPCA.quantumUtilities.Tomography.StateVectorTomography` class is particularly important, as it is used by Qpca for eigenvector reconstruction. However, it can also be used as a standalone class to reconstruct the state vector of an arbitrary quantum circuit composed only of real amplitudes.

* :mod:`~QPCA.preprocessingUtilities` module: this module includes the :meth:`~QPCA.preprocessingUtilities.generate_matrix` method, which is useful for generating an arbitrary random Hermitian input matrix to be used in the Qpca algorithm.

* :mod:`~QPCA.postprocessingUtilities` module: this module provides the :meth:`~QPCA.postprocessingUtilities.general_postprocessing` method, which is used by the Qpca algorithm to recontruct the eigenvectors and eigenvalues using the information obtained from quantum state tomography.

* :mod:`~QPCA.benchmark` module: this module contains the :class:`~QPCA.benchmark.Benchmark_Manager` class, which manages all the benchmarking tasks for the Qpca algorithm. It provides methods to benchmark the execution of Qpca and analyze the accuracy of eigenvector and eigenvalue reconstruction. 

Each module plays a specific role in the overall Qpca package, allowing for an end-to-end process from input matrix encoding to eigenvector and eigenvalue reconstruction using quantum routines.

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
