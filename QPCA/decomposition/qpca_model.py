import math
import numpy as np
from IPython.display import display
from scipy.spatial import distance

from ..benchmark.benchmark import Benchmark_Manager
from ..postprocessingUtilities.postprocessing import general_postprocessing
from ..preprocessingUtilities.preprocessing import check_matrix_dimension
from ..quantumUtilities.phase_estimation_builder import PeCircuitBuilder
from ..quantumUtilities.qram_circuit_builder import QramBuilder
from ..quantumUtilities.state_vector_tomography import StateVectorTomography
from ..warnings_utils.warning_utility import *

__all__ = ["QPCA"]


class QPCA():
    """Quantum Principal component analysis (QPCA).
    Implementation of the QPCA algorithm proposed in "A Low Complexity Quantum Principal Component Analysis Algorithm" paper.

    Parameters
    ----------


    Attributes
    ----------

    input_matrix_trace : NumPy array-like of shape (n_features,)
                        The trace of the input matrix.

    input_matrix : array-like of shape (n_samples, n_features)
                        Input hermitian matrix on which you want to apply QPCA, divided by its trace. Here, `n_samples` represents the number of samples,
                        and `n_features` represents the number of features. In case of non-2^N Hermitian matrix, a zero-padding method is applied to make it a 2^N matrix.

    true_input_matrix : array-like of shape (n_samples, n_features)
                        This is the true input matrix that is given as input.

    resolution : int value
                        The number of qubits used in the phase estimation process to encode the eigenvalues.

    qram_circuit : QuantumCircuit
                        The quantum circuit that encodes the input matrix.

    total_circuit : QuantumCircuit
                        The quantum circuit that performs the encoding of the input matrix and Phase Estimation. The number of qubits will be log(n_samples*n_features)+resolution.

    n_shots : int value
                        Number of measures performed in the tomography process.

    mean_threshold: array-like
                        This array contains the mean between the left and right peaks vertical distance to its neighbouring samples. It is used for the benchmark process to cut out the bad reconstructed eigenvalues.

    reconstructed_eigenvalues : array-like
                        Reconstructed eigenvalues.

    reconstructed_eigenvectors : array-like
                        Reconstructed eigenvectors.

    Notes
    ----------
    It is important to consider that the input matrix is divided by its trace so as to have the eigenvalues between 0 and 1.

    """
    
    def fit(self, input_matrix, resolution, optimized_qram=True, plot_qram=False,plot_pe_circuit=False):
        
        """
        Fit Qpca model. This method generates the encoding matrix circuit and apply the phase estimation operator.

        Parameters
        ----------
        
        input_matrix: array-like of shape (n_samples, n_features)
                    Input hermitian matrix on which you want to apply QPCA, divided by its trace. Here, `n_samples` represents the number of samples,
                    and `n_features` represents the number of features.
                    
        resolution: int value
                    The number of qubits used in the phase estimation process to encode the eigenvalues.
                    
        optimized_qram: bool value, default=True
                        If True, it returns an optimized version of the preprocessing circuit. Otherwise, a custom implementation of a Qram is returned.
                        Unless necessary, it is recommended to keep the optimized version of this circuit.
                    
        plot_qram: bool value, default=False
                    If True, it returns a plot of the Qram circuit that encodes the input matrix.
        
        plot_pe_circuit: bool value, default=False
                    If True, it returns a plot of the circuit composed of Qram and phase estimation operator.
                    

        Returns
        ----------
        self: object
                    Returns the instance itself.
                    
        Notes
        ----------
        """
        
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            
            # 8 qubit resolution because accuracy problems are more common below this value
            
            if resolution<8:
                customWarning.warn(f'You chose {resolution} qubits of resolution. Moreover, since with {resolution} qubits you have an accuracy of {1/2**resolution}, if you know that some eigenvalues are smaller or closer to each other than {1/2**resolution}, please increase the resolution qubits to get better estimates.')
        
        true_input_matrix=input_matrix
        
        # 2^N matrix optimization
        
        input_matrix=check_matrix_dimension(input_matrix)
        
        self.input_matrix_trace=np.trace(input_matrix)
        
        #normalize the input matrix by its trace to obtain eigenvalues between 0 and 1
        
        self.input_matrix=input_matrix/self.input_matrix_trace
        self.true_input_matrix=true_input_matrix/np.trace(true_input_matrix)
        self.resolution=resolution
    
        qc=QramBuilder.generate_qram_circuit(self.input_matrix, optimized_qram=optimized_qram)
        self.qram_circuit=qc
        
        if plot_qram:
            display(qc.draw('mpl'))
        
        pe_circuit=PeCircuitBuilder.generate_PE_circuit(input_matrix=self.input_matrix,resolution=self.resolution,qram_circuit=self.qram_circuit)
        self.total_circuit=pe_circuit
        if plot_pe_circuit:
            display(pe_circuit.draw('mpl'))
        
        return self


    def eigenvectors_reconstruction(self,n_shots=10000,n_repetitions=1,plot_peaks=False, backend=None, eigenvalue_threshold=None, abs_tolerance=None):
        
        """ Method that reconstructs the eigenvalues/eigenvectors once performed Phase Estimation. 

        Parameters
        ----------
        n_shots: int value, default=10000
                Number of measures performed in the tomography process.
                
        n_repetitions: int value, default=1
                Number of times that state vector tomography will be executed. If > 1, the final result will be the average result
                of all the execution of the tomography.
                
        plot_peaks: bool value, defualt=False
                        If True, it returns a plot of the peaks which correspond to the eigenvalues finded by the phase estimation procedure.
        
        backend: Qiskit backend, default value=None.
                    The Qiskit backend used to execute the circuit. If None, the qasm simulator is used by default.
        
        eigenvalue_threshold: float value, default=None
                        It acts as a threshold that cut out the eigenvalues (and the corrseponding eigenvectors) that are smaller than this value.
        
        abs_tolerance: float value, default=None
                        Absolute tolerance parameter used to cut out the eigenvalues estimated badly due to insufficient resolution.
        
        Returns
        ----------
        reconstructed_eigenvalues : array-like
                        Reconstructed eigenvalues.
    
        reconstructed_eigenvectors : array-like
                            Reconstructed eigenvectors.
                    
        Notes
        ----------
        To classically reconstruct the eigenvectors, state vector tomography function is performed (implemented from algorithm 4.1 of "A Quantum Interior Point Method for LPs and SDPs" paper). In this
        way, the statevector of the quantum state is reconstructed and a postprocessing method is executed to get the eigenvectors from the reconstructed statevector.
        """
        
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            
            if n_shots<10000:
                customWarning.warn("You are performing the tomography procedure with less than 10.000 measures. Note that to obtain accurate estimates, it is recommended to carry out at least 10.000 measurements.")
            
        self.n_shots=n_shots
        
        statevector_dictionary=StateVectorTomography.state_vector_tomography(quantum_circuit=self.total_circuit,n_shots=n_shots,n_repetitions=n_repetitions, backend=backend)
        
        eigenvalue_eigenvector_tuple,mean_threshold=general_postprocessing(input_matrix=self.input_matrix,statevector_dictionary=statevector_dictionary,resolution=self.resolution,
                                                                           n_shots=self.n_shots,plot_peaks=plot_peaks,eigenvalue_threshold=eigenvalue_threshold,abs_tolerance=abs_tolerance)
        
        self.mean_threshold=mean_threshold[:len(self.true_input_matrix)]
        
        reconstructed_eigenvalues=np.array([])
        reconstructed_eigenvectors=np.array([])
        
        #reshape the reconstructed eigenvectors. In case of previous padding, remove the unnecessary zero rows/columns
        
        for t in eigenvalue_eigenvector_tuple[:len(self.true_input_matrix)]:
            
            reconstructed_eigenvalues=np.append(reconstructed_eigenvalues,t[0])    
            reconstructed_eigenvectors=np.append(reconstructed_eigenvectors,t[1][:len(self.true_input_matrix)])
        try:
            reconstructed_eigenvectors=reconstructed_eigenvectors.reshape(len(self.true_input_matrix),len(reconstructed_eigenvalues),order='F')
        except:
            raise Exception('QPCA was not able to correctly reconstruct the eigenvectors! Check that you are not considering eigenvalues near to zero. In that case, you can both increase the number of shots or the resolution for the phase estimation.')
        
        self.reconstructed_eigenvalues=reconstructed_eigenvalues
        self.reconstructed_eigenvectors=reconstructed_eigenvectors
        
        return reconstructed_eigenvalues,reconstructed_eigenvectors
    
    def quantum_input_matrix_reconstruction(self):
        
        """ Method to reconstruct the input matrix.

        Parameters
        ----------

        Returns
        ----------
        reconstructed_input_matrix: array-like of shape (n_samples, n_features)
                Reconstructed input matrix.
                    
        Notes
        ----------
        Using the reconstructed eigenvectors and eigenvalues from QPCA, we can reconstruct the original input matrix using the reverse procedure of SVD.
        """
        
        
        k = self.reconstructed_eigenvalues.argsort()[::-1]   
        reconstructed_eigenvalues = self.reconstructed_eigenvalues[k]
        reconstructed_eigenvectors = self.reconstructed_eigenvectors[:,k]
        reconstructed_eigenvalues*=self.input_matrix_trace
        
        #reconstruct the input matrix by multiplying the eigenvectors/eigenvalues matrices 
        
        reconstructed_input_matrix = reconstructed_eigenvectors @ np.diag(reconstructed_eigenvalues) @ reconstructed_eigenvectors.T
        return reconstructed_input_matrix
    
    def transform(self, input_matrix):
        
        """Apply dimensionality reduction to input_matrix.

        Input_matrix is projected on the first principal components previously extracted.

        Parameters
        ----------
        input_matrix : array-like of shape (n_samples, n_features)
            New data, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        Returns
        -------
        input_matrix_transformed : array-like of shape (n_samples, n_components)
            Projection of input_matrix in the first principal components, where `n_samples`
            is the number of samples and `n_components` is the number of the components.
        """
        k = self.reconstructed_eigenvalues.argsort()[::-1]   
        reconstructed_eigenvectors = self.reconstructed_eigenvectors[:,k]
        input_matrix_transformed=np.dot(input_matrix, reconstructed_eigenvectors)
        
        return input_matrix_transformed
    
    def spectral_benchmarking(self, eigenvector_benchmarking=False, eigenvalues_benchmarching=False, sign_benchmarking=False, print_distances=True,
                              only_first_eigenvectors=False, plot_delta=False, distance_type='l2', error_with_sign=False, hide_plot=False, print_error=False):

            """ Method to benchmark the reconstructed eigenvectors/eigenvalues.

            Parameters
            ----------
            eigenvector_benchmarking: bool value, default=True
                    If True, an eigenvectors benchmarking is performed to show the accuracy for the quantum algorithm in estimating the original eigenvectors.

            eigenvalues_benchmarching: bool value, default=False
                    If True, an eigenvalues benchmarking is performed to show the accuracy for the quantum algorithm in estimating the original eigenvalues.
            
            sign_benchmarking: bool value, default=False
                    If True, a table showing correct and wrong signs for each reconstructed eigenvector is returned.

            print_distances: bool value, default=True
                    If True, the distance (defined by the parameter distance_type) between the reconstructed and original eigenvectors is printed in the legend.

            only_first_eigenvectors: bool value, default=True
                    If True, the benchmarking is performed only for the first eigenvector. Otherwise, all the eigenvectors are considered.

            plot_delta: bool value, default=False
                    If True, a plot showing the trend of the tomography error is showed.

            distance_type: string value, default='l2'
                    It defines the distance measure used to benchmark the eigenvectors:

                        -'l2': the l2 distance between original and reconstructed eigenvectors is computed.

                        -'cosine': the cosine distance between original and reconstructed eigenvectors is computed.
            
            error_with_sign: bool value, default=False
                    If True, the eigenvectors' benchmarking is performed considering the reconstructed sign. Otherwise, the benchmark is performed only for the absolute values of the eigenvectors (which means reconstructed eigenvectors with no 
                    reconstructed sign). If True, it has effect only if the eigenvector_benchmarking flag is True.
                    
            hide_plot: bool value, default=False
                    If True, the plot for the eigenvector reconstruction benchmarking is not showed. This is useful to have a cleaner output when executing the eigenvectors reconstruction benchmarking more times (for example for different matrices). 
            
            print_error: bool value, default=False
                    If True, a table showing the absolute error between the original eigenvalues and the reconstructed ones is shown.
                
                        

            Returns
            -------
            results: array-like
                    If eigenvector_benchmarking is True, it contains a list of two elements:
                    
                        - error_list: array-like. 
                            List of distances between the original and reconstructed eigenvectors.
                        - delta: float value.
                            Tomography error value.

            Notes
            -----
            """
            
            bench_manager=Benchmark_Manager(eigenvector_benchmarking=eigenvector_benchmarking, eigenvalues_benchmarching=eigenvalues_benchmarching, sign_benchmarking=sign_benchmarking, print_distances=print_distances,
                                            only_first_eigenvectors=only_first_eigenvectors, plot_delta=plot_delta, distance_type=distance_type, error_with_sign=error_with_sign, hide_plot=hide_plot, print_error=print_error)
            
            results=bench_manager.benchmark(input_matrix=self.true_input_matrix, reconstructed_eigenvalues=self.reconstructed_eigenvalues, 
                                            reconstructed_eigenvectors=self.reconstructed_eigenvectors, mean_threshold=self.mean_threshold, 
                                            n_shots=self.n_shots, resolution=self.resolution)
            
            
            return results
            
    