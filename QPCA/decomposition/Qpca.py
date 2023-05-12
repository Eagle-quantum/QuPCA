import numpy as np
import warnings
import math
from ..quantumUtilities.Tomography import StateVectorTomography
from ..quantumUtilities.qRam_Builder import QramBuilder
from ..quantumUtilities.qPe_Builder import PeCircuitBuilder
from ..postprocessingUtilities.postprocessing_eig_reconstruction import general_postprocessing
from ..preprocessingUtilities.preprocessing_matrix_utilities import next_power_of_2
from ..benchmark.benchmark import eigenvectors_benchmarking,eigenvalues_benchmarking,error_benchmark,sign_reconstruction_benchmarking,distance_function_wrapper
from scipy.spatial import distance
#warnings.filterwarnings("ignore")

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
                        Input hermitian matrix on which you want to apply QPCA divided by its trace, where `n_samples` is the number of samples
                        and `n_features` is the number of features. In case of non-2^N Hermitian matrix, a zero-padding method is applied to make it a 2^N matrix.
    
    true_input_matrix : array-like of shape (n_samples, n_features)
                        This is the true input matrix that is given as input. 
    
    resolution : int value
                        Number of qubits used for the phase estimation process to encode the eigenvalues.
    
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
                    Input hermitian matrix on which you want to apply QPCA divided by its trace, where `n_samples` is the number of samples
                    and `n_features` is the number of features.
                    
        resolution: int value
                    Number of qubits used for the phase estimation process to encode the eigenvalues.
                    
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
        
        true_input_matrix=input_matrix
        matrix_dimension=len(input_matrix)
        
        #check if the matrix dimension is 2^N. If not, pad it with 0
        
        if ((matrix_dimension & (matrix_dimension-1) == 0) and matrix_dimension != 0)==False:
            zeros=np.zeros((matrix_dimension,1))
            zeros_r=np.zeros((1,next_power_of_2(matrix_dimension)))
            for i in range(next_power_of_2(matrix_dimension)-matrix_dimension):
                input_matrix=np.append(input_matrix,zeros,axis=1)
            for i in range(next_power_of_2(matrix_dimension)-matrix_dimension):
                input_matrix=np.append(input_matrix,zeros_r,axis=0)
        
        self.input_matrix_trace=np.trace(input_matrix)
        
        #normalize the input matrix by its trace to obtain eigenvalues between 0 and 1
        
        self.input_matrix=input_matrix/np.trace(input_matrix)
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


    def eigenvectors_reconstruction(self,n_shots=50000,n_repetitions=1,plot_peaks=False):
        
        """ Method that reconstructs the eigenvalues/eigenvectors once performed Phase Estimation. 

        Parameters
        ----------
        n_shots: int value, default=50000.
                Number of measures performed in the tomography process.
                
        n_repetitions: int value, default=1.
                Number of times that state vector tomography will be executed. If > 1, the final result will be the average result
                of all the execution of the tomography.
                
        plot_peaks: bool value, defualt=False
                        If True, it returns a plot of the peaks which correspond to the eigenvalues finded by the phase estimation procedure.
        
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
        
        
        self.n_shots=n_shots
        
        statevector_dictionary=StateVectorTomography.state_vector_tomography(quantum_circuit=self.total_circuit,n_shots=n_shots,n_repetitions=n_repetitions)
        
        eigenvalue_eigenvector_tuple,mean_threshold=general_postprocessing(input_matrix=self.input_matrix,statevector_dictionary=statevector_dictionary,
                                                                           resolution=self.resolution,n_shots=self.n_shots,plot_peaks=plot_peaks)
        
        self.mean_threshold=mean_threshold[:len(self.true_input_matrix)]
        
        reconstructed_eigenvalues=np.array([])
        reconstructed_eigenvectors=np.array([])
        
        #reshape the reconstructed eigenvectors. In case of previous padding, remove the unnecessary zero rows/columns
        
        for t in reconstructed_eigenvalue_eigenvector_tuple[:len(self.true_input_matrix)]:
            
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
    
    
    def spectral_benchmarking(self, eigenvector_benchmarking=True, eigenvalues_benchmarching=False,sign_benchmarking=False,print_distances=True,only_first_eigenvectors=True,plot_delta=False,
                                  distance_type='l2',error_with_sign=False,hide_plot=False,print_error=False):

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
            If eigenvector_benchmarking is True:

                - save_list: array-like. 
                    List of distances between all the original and reconsructed eigenvectors.
                - delta: float value.
                    The tomography error value.

            Notes
            -----
            """

            eigenValues,eigenVectors=np.linalg.eig(self.true_input_matrix)
            idx = eigenValues.argsort()[::-1]   
            original_eigenValues = eigenValues[idx]
            original_eigenVectors = eigenVectors[:,idx]


            if eigenvector_benchmarking:
                error_list, delta=eigenvectors_benchmarking(input_matrix=self.true_input_matrix, original_eigenvalues=original_eigenValues, original_eigenvectors=original_eigenVectors,
                                                            reconstructed_eigenvalues=self.reconstructed_eigenvalues, reconstructed_eigenvectors=self.reconstructed_eigenvectors,
                                                            mean_threshold=self.mean_threshold, n_shots=self.n_shots,print_distances=print_distances, only_first_eigenvectors=only_first_eigenvectors,
                                                            plot_delta=plot_delta,distance_type=distance_type,error_with_sign=error_with_sign,hide_plot=hide_plot)
            if eigenvalues_benchmarching:
                eigenvalues_benchmarking(original_eigenvalues=original_eigenValues,
                                         reconstructed_eigenvalues=self.reconstructed_eigenvalues,
                                         mean_threshold=self.mean_threshold, print_error=print_error)

            if sign_benchmarking:
                sign_reconstruction_benchmarking(input_matrix=self.true_input_matrix, original_eigenvalues=original_eigenValues, original_eigenvectors=original_eigenVectors,
                                                 reconstructed_eigenvalues=self.reconstructed_eigenvalues, reconstructed_eigenvectors=self.reconstructed_eigenvectors,
                                                 mean_threshold=self.mean_threshold,n_shots=self.n_shots)

            if eigenvector_benchmarking:
                
                return error_list, delta
    
            
            
    def error_benchmarking(self,shots_dict=None,errors_dict=None,delta_list=None,plot_delta=False,distance_type='l2'):
        
        """ Method to benchmark the error of the reconstructed eigenvectors.

        Parameters
        ----------
        
        shots_dict: dict-like
                Dictionary that contains the number of shots executed for each reconstructed eigenvalues for a specific resolution value.
                
        errors_dict: dict-like
                Dictionary where a specific value of the resolution used in the phase estimation (key) is associated with a list of tuples where each tuple contains the specific reconstructed eigenvalue with all the error of the related reconstructed eigenvector (for each number of shots).
        
        delta_list: array-like
                List of all the delta (tomography) error for each number of shots executed in the experiments.
                
        plot_delta: bool value, default=False.
                If True, a plot showing the trend of the tomography error is showed.
                
        distance_type: string value, default='l2'
            It defines the distance measure used to benchmark the eigenvectors:

                    -'l2': the l2 distance between original and reconstructed eigenvectors is computed.

                    -'cosine': the cosine distance between original and reconstructed eigenvectors is computed.

        Returns
        -------
                    
        Notes
        -----
        This method should be executed after the execution of eigenvectors_benchmarking() method in the spectral_benchmarking, just to visualize better in specific plots the trends of the reconstruction error for each reconstructed eigenvectors. More precisely, it should be executed after executing eigenvectors_benchmarking method for a different number of shots such that you can visualize better the error trend. You can also use this method to visualize the trends of the reconstruction error for a different number of shots at the increasing of the number of resolution qubits. The important thing to take into consideration is that you need to pass as argument the two dictionary described in the documentation (see the benchmark notebook for a more practical example)
        """

        eigenValues,eigenVectors=np.linalg.eig(self.true_input_matrix)
        idx = eigenValues.argsort()[::-1]   
        original_eigenValues = eigenValues[idx]
        original_eigenVectors = eigenVectors[:,idx]
        dict_original_eigenvalues={}
        
        # dictionary useful to make comparison betweeen the correct eigenvectors
        
        for i in range(len(original_eigenValues)):
            dict_original_eigenvalues.update({original_eigenValues[i]:i})
            
        error_benchmark(shots_dict=shots_dict,error_dict=errors_dict,dict_original_eigenvalues=dict_original_eigenvalues,delta_list=delta_list,plot_delta=plot_delta,label_error=distance_type)