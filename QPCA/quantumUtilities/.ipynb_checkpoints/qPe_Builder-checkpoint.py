import numpy as np

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit.algorithms.linear_solvers.matrices.numpy_matrix import NumPyMatrix
from qiskit.circuit.library import PhaseEstimation

class PeCircuitBuilder():
    
    @classmethod
    def generate_PE_circuit(class_,input_matrix,resolution,qram_circuit):
        """
        Generate phase estimation circuit with a number of qubits provided as resolution parameter in the constructor.

        Parameters
        ----------
        
        input_matrix: array-like of shape (n_samples, n_features)
                        Input hermitian matrix on which you want to apply QPCA divided by its trace, where `n_samples` is the number of samples
                        and `n_features` is the number of features.
        
        resolution: int value
                    Number of qubits used for the phase estimation process to encode the eigenvalues.
        
        qram_circuit: QuantumCircuit 
                        The quantum circuit that encodes the input matrix.

        Returns
        ----------
        q_circuit: QuantumCircuit. 
                The quantum circuit that performs the encoding of the input matrix and Phase Estimation.
                    
        Notes
        ----------

        """
        
        u_circuit = NumPyMatrix(input_matrix, evolution_time=2*np.pi)
        pe = PhaseEstimation(resolution, u_circuit, name = "PE")
        tot_qubit = pe.qregs[0].size+qram_circuit.qregs[0].size
        qr_total = QuantumRegister(tot_qubit, 'total')
        q_circuit = QuantumCircuit(qr_total , name='matrix')
        q_circuit.append(qram_circuit.to_gate(), qr_total[pe.qregs[0].size:])
        q_circuit.append(pe.to_gate(), qr_total[0:pe.num_qubits])
        return q_circuit

    