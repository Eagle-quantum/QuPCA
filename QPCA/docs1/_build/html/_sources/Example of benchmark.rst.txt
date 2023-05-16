Example
====

In the :mod:`~QPCA.benchmark` module, you will find 4 methods to benchmark the execution of Qpca algorithm.
Here below there are reported some possible usage of the methods of benchmarking.

First, a seed is set for the reproducibility of the experiments. Then, the resolution is set to 8 qubits and 
the matrix that will be considered is a 4x4 matrix with a custom list of eigenvalues. The number of measures for 
the eigenvectors reconstruction is set to 1000000.

.. code-block:: python

   from QPCA.decomposition.Qpca import QPCA
   import numpy as np
   import matplotlib.pyplot as plt
   from QPCA.preprocessingUtilities.preprocessing_matrix_utilities import generate_matrix

   seed=4747
   resolutions=[8]
   matrix_dimension=4
   eigenvalues_list=[0.65,0.25,0.06,0.04]
   input_matrix=generate_matrix(matrix_dimension=matrix_dimension,replicate_paper=False,seed=seed,eigenvalues_list=eigenvalues_list)
   shots_numbers=[1000000] 

Then, a Qpca object is fitted using the generated input matrix and the chosen resolution. Finally, using :meth:`~QPCA.decomposition.QPCA.eigenvectors_reconstruction`, you are able to reconstruct the eigenvalues/eigenvectors.
The first reference example concerns the accuracy of reconstructing eigenvectors using :meth:`~QPCA.decomposition.QPCA.spectral_benchmarking`

.. code-block:: python

   for resolution in resolutions:
      qpca=QPCA().fit(input_matrix,resolution=resolution)
      for s in shots_numbers:
         reconstructed_eigenvalues,reconstructed_eigenvectors=qpca.eigenvectors_reconstruction(n_shots=s,n_repetitions=1)
         eig_evec_tuple,delta=qpca.spectral_benchmarking(eigenvector_benchmarking=True,sign_benchmarking=False ,eigenvalues_benchmarching=False,print_distances=True,only_first_eigenvectors=False,
                                                        plot_delta=True,distance_type='l2',error_with_sign=True,hide_plot=False,print_error=False)

Eigenvectors benchmark
~~~~~~~~~~~~~~~~~~~~~~

Setting the :obj:`~QPCA.benchmark.eigenvectors_reconstruction.eigenvector_benchmarking` parameter to True, you will obtain a plot like the following one. You can see the eigenvectors that you are able to reconstruct and the values of the reconstructed 
eigenvector (stars) compared to the original ones (circles). In the legend, it is also reported the l2-error distance between the reconstructed and original eigenvector.

.. image:: Images/benchmark1.png

Eigenvalues benchmark
~~~~~~~~~~~~~~~~~~~~~~

You can also benchmark the reconstructed eigenvalues by setting the :obj:`~QPCA.benchmark.eigenvectors_reconstruction.eigenvalues_benchmarching` to True. With this benchmark, a plot showing the reconstructed eigenvalues (red circle) and the 
original ones (black crosses) is shown. If the :obj:`~QPCA.benchmark.eigenvectors_reconstruction.print_error` parameter is set to True, a table showing the absolute error between reconstructed and original eigenvalues it is also reported.

.. image:: Images/benchmark2.png

Eigenvectors reconstruction error benchmark
~~~~~~~~~~~~~~~~~~~~~~

Using the :meth:`~QPCA.decomposition.QPCA.error_benchmarking` method, you can visualize better the trend of the reconstruction error for each eigenvectors as the number of measures and number of 
resolution qubits increase. As before, once the number of measures and resolution qubits are chosen, you can perform the fit and eigenvectors reconstruction procedures. 
Pay attention: it is important to save the results of the benchmark into specific dictionary, as in the code below. This is because the :meth:`~QPCA.decomposition.QPCA.error_benchmarking` function 
expects dictionaries as parameters.

.. code-block:: python
   
   shots_numbers=[100,500,1500,10000,100000,500000,1000000]
   resolutions=[3,5,8]
   resolution_dictionary={}
   resolution_dictionary_shots={}
   for resolution in resolutions:
      error_list=[]
      delta_list=[]
      shots_dict={}
      qpca=QPCA().fit(input_matrix,resolution=resolution)
      for s in shots_numbers:
         
         reconstructed_eigenvalues,reconstructed_eigenvectors=qpca.eigenvectors_reconstruction(n_shots=s,n_repetitions=1)
         eig_evec_tuple,delta=qpca.spectral_benchmarking(eigenvector_benchmarking=True,sign_benchmarking=False ,eigenvalues_benchmarching=False,print_distances=True,only_first_eigenvectors=False,
                                                         plot_delta=True,distance_type='l2',error_with_sign=True,hide_plot=False,print_error=False)
         for e in eig_evec_tuple:
               shots_dict.setdefault(e[0], []).append(s)
         error_list.append(eig_evec_tuple)
         delta_list.append(delta)
      
      resolution_dictionary_shots.update({resolution:shots_dict})
      resolution_dictionary.update({resolution:error_list})

   qpca.error_benchmarking(errors_dict=resolution_dictionary,shots_dict=resolution_dictionary_shots,distance_type='l2')


.. image:: Images/benchmark3.png

With these plots, you can observe the trend of the errors as the number of measures and resolution qubits increases for each eigenvectors.