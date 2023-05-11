import numpy as np

from qiskit.circuit.library.standard_gates import RYGate
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit.circuit.library import PhaseEstimation
from qiskit import Aer, transpile
import matplotlib.pyplot as plt
from qiskit.circuit.library.data_preparation.state_preparation import StatePreparation

class StateVectorTomography():

    def __computing_amplitudes(quantum_circuit,q_size,c_size,n_shots,drawing_amplitude_circuit,backend,qubits_to_be_measured):
        
        """ This is the first step of the state vector tomography algorithm described in Algorithm 4.1 in "A Quantum Interior Point Method for LPs and SDPs" paper. It is
            useful to reconstruct each statevector components.

            Parameters
            ----------

            q_size: int value
                        The size of the quantum register under consideration.

            c_size: int value
                        The size of the classical register under consideration.

            qubits_to_be_measured: Union[Qubit, QuantumRegister, int, slice, Sequence[Union[Qubit, int]]]), default=None
                        Qubits to be measured. If None, all the qubits will be measured (like measure_all() instruction).

            Returns
            -------
            probabilities: array-like
                        The reconstructed probabilities statevector (without sign reconstruction).

            Notes
            -----
            This method is used to reconstruct the amplitudes of the values of the statevector that is under consideration. In addition to this function, state vector tomography also includes a sign estimation function.
            """
        
        #initialize a zero vector where we store the estimated probabilities
        
        probabilities=np.zeros(2**c_size)
        quantum_regs_1=QuantumRegister(q_size)
        classical_regs_1 = ClassicalRegister(c_size, 'classical')
        tomography_circuit_1 = QuantumCircuit(quantum_regs_1,classical_regs_1)
        tomography_circuit_1.append(quantum_circuit,quantum_regs_1)
        tomography_circuit_1.measure(quantum_regs_1[qubits_to_be_measured],classical_regs_1)
        
        if drawing_amplitude_circuit:
            display(tomography_circuit_1.draw('mpl'))
        
        job = backend.run(transpile(tomography_circuit_1, backend=backend), shots=n_shots)
        counts = job.result().get_counts()
        
        #compute estimated probabilities as number of observation for the i-th state divided by the total number of shots performed
        
        for i in counts:
            counts[i]/=n_shots
            probabilities[int(i,2)]=counts[i]
        
        return probabilities
    
    def __sign_estimation(quantum_circuit,probabilities,q_size,c_size,n_shots,drawing_sign_circuit,backend,qubits_to_be_measured):
        
        """ This is the second and last step of the state vector tomography algorithm described in Algorithm 4.1 in "A Quantum Interior Point Method for LPs and SDPs" paper. It is
            useful to reconstruct the sign of each statevector's components.

            Parameters
            ----------
            
            probabilities: array-like 
                        The reconstructed probabilities statevector (without sign reconstruction) obtained from computing_amplitudes function.
            
            q_size: int value 
                        The size of the quantum register under consideration.

            c_size: int value 
                        The size of the classical register under consideration.

            qubits_to_be_measured: Union[Qubit, QuantumRegister, int, slice, Sequence[Union[Qubit, int]]]), default=None
                        Qubits to be measured. If None, all the qubits will be measured (like measure_all() instruction).

            Returns
            -------

            statevector_dictionary: dict-like
                        Dictionary where the keys represent the eigenvalues/eigenvectors encoded in the qubits and the values represent the reconstructed statevector's values (with sign).

            Notes
            -----
            This method is used to reconstruct the correct sign of the statevector's values under consideration.
            """
        
        
        qr_total_xi = QuantumRegister(q_size, 'xi')
        n_classical_register=c_size+1
        classical_registers=ClassicalRegister(n_classical_register,'classical')
        qr_control = QuantumRegister(1, 'control_qubit')
        
        #creation of the first controlled operator U from the circuit that contains the matrix encoding and the phase estimation operator
        
        op_U=quantum_circuit.to_gate(label='op_U').control()
        
        #initialize the second controlled operator with the estimated amplitudes computed in the first step of the tomography
        
        op_V = StatePreparation(np.sqrt(probabilities),label='op_V').control()
        
        #implement the sign estimation circuit as an Hadamard test

        tomography_circuit_2 = QuantumCircuit(qr_total_xi,qr_control, classical_registers,name='matrix')
        tomography_circuit_2.h(qr_control)
        tomography_circuit_2.x(qr_control)
        tomography_circuit_2.append(op_U, qr_control[:]+qr_total_xi[:])
        tomography_circuit_2.x(qr_control)
        tomography_circuit_2.append(op_V, qr_control[:]+qr_total_xi[qubits_to_be_measured])
        tomography_circuit_2.h(qr_control)
        tomography_circuit_2.measure(qr_total_xi[qubits_to_be_measured],classical_registers[0:n_classical_register-1])
        tomography_circuit_2.measure(qr_control,classical_registers[n_classical_register-1])
    
        if drawing_sign_circuit:
            display(tomography_circuit_2.draw('mpl'))

        job = backend.run(transpile(tomography_circuit_2, backend=backend), shots=n_shots)
        counts_for_sign = job.result().get_counts()
        tmp=np.zeros(2**c_size)
        
        #check the sign: we consider only the results with control qubit 0
        
        for c in counts_for_sign:
            if c[0]=='0':
                tmp[int(c[1:],2)]=counts_for_sign[c]
        sign_dictionary={}
        sign=0
        for e, (count, prob) in enumerate(zip(tmp, probabilities)):
            if count>0.4*prob*n_shots:
                sign=1
            else:
                sign=-1
            if prob==0:
                sign=1
            sign_dictionary.update({bin(e)[2:].zfill(c_size):sign})

        statevector_dictionary={}
        
        #merge the results of amplitude and sign estimation in a dictionary where the keys are the states and the values are the signed amplitudes
        
        for e,key in enumerate(sign_dictionary):
            statevector_dictionary[key]=sign_dictionary[key]*np.sqrt(probabilities[e])
        return statevector_dictionary
    
    @classmethod
    def state_vector_tomography(class_,quantum_circuit,n_shots,qubits_to_be_measured=None,backend=None,drawing_amplitude_circuit=False,drawing_sign_circuit=False):
        """
        State vector tomography to estimate real vectors.

        Parameters
        ----------

        quantum_circuit: QuantumCircuit 
                    The quantum circuit to be reconstructed. 

        n_shots: int value
                    Number of measures performed in the tomography process.
                    
        qubits_to_be_measured: Union[Qubit, QuantumRegister, int, slice, Sequence[Union[Qubit, int]]]), default=None.
                    Qubits to be measured. If None, all the qubits will be measured (like measure_all() instruction).
        
        backend: Qiskit backend, default value=None.
                    The Qiskit backend used to execute the circuit. If None, the qasm simulator is used by default.
        
        drawing_amplitude_circuit: bool value, default=False.
                    If True, a drawing of the amplitude estimation circuit of the tomography algorithm is displayed. Otherwise, only the reconstructed statevector is returned.
        
        drawing_sign_circuit: bool value, default=False.
                    If True, a drawing of the sign estimation circuit of the tomography algorithm is displayed. Otherwise, only the reconstructed statevector is returned.

        Returns
        -------
        statevector_dictionary: dict-like. 
                    The reconstructed statevector of the input quantum circuit.

        Notes
        -----
        This method reconstruct the real statevector of the input quantum circuit. It is an implementation of the Algorithm 4.1 in "A Quantum Interior Point Method for LPs and SDPs" paper, and it is composed of 
        two parts: amplitudes estimation and sign estimation.
        """
        
        if backend==None:
            backend = Aer.get_backend("qasm_simulator")
            
        q_size=quantum_circuit.qregs[0].size
        if qubits_to_be_measured==None:
            c_size=q_size
            qubits_to_be_measured=list(range(q_size))
        elif isinstance(qubits_to_be_measured,int):
            c_size=1
        else:
            tmp_array=np.array(list(range(q_size)))
            c_size=len(tmp_array[qubits_to_be_measured])

        probabilities=class_.__computing_amplitudes(quantum_circuit,q_size,c_size,n_shots,drawing_amplitude_circuit,backend,qubits_to_be_measured)
        statevector_dictionary=class_.__sign_estimation(quantum_circuit,probabilities,q_size,c_size,n_shots,drawing_sign_circuit,backend,qubits_to_be_measured)

        return statevector_dictionary