import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library.data_preparation.state_preparation import StatePreparation

from .quantum_utilities import from_binary_tree_to_qcircuit, thetas_computation

__all__ = ["QramBuilder"]

class QramBuilder():
    
    @classmethod
    def generate_qram_circuit(cls,input_matrix, optimized_qram):
        """
        Generate qram circuit.

        Parameters
        ----------
        
        input_matrix: array-like of shape (n_samples, n_features)
                        Input hermitian matrix on which you want to apply QPCA, divided by its trace. Here, `n_samples` represents the number of samples,
                        and `n_features` represents the number of features. 
        
        optimized_qram: bool value
                        If True, it returns an optimized version of the preprocessing circuit. Otherwise, a custom implementation of a Qram is returned.
                        Unless necessary, it is recommended to keep the optimized version of this circuit.

        Returns
        ----------
        qc: QuantumCircuit 
                    Preprocessing circuit.
        
        Notes
        ----------
        This method implements the quantum circuit generation to encode a generic input matrix. For the custom implementation, it is important to note the spatial complexity of the circuit is in the order of log2(n_samples, n_features).
                    
        """
        if optimized_qram:

            flattened_matrix = input_matrix.flatten()
            norm = np.linalg.norm(flattened_matrix)
            state_preparation = StatePreparation(flattened_matrix / norm)
            num_qubits=int(np.ceil(np.log2(len(flattened_matrix))))

            qc = QuantumCircuit(num_qubits)

            qc.append(state_preparation, [i for i in range(num_qubits-1,-1,-1)])
            
        else:
            
            thetas, all_combinations=thetas_computation(input_matrix=input_matrix)
        
            qc=from_binary_tree_to_qcircuit(input_matrix,thetas, all_combinations)      

        return qc

    