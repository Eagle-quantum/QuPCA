import numpy as np
import itertools
import warnings
import math
from numpyarray_to_latex.jupyter import to_jup
from numpyarray_to_latex import to_ltx

def generate_matrix(matrix_dimension,eigenvalues_list=None,replicate_paper=True,seed=None):
        
        """ Hermitian matrix generation.
        
        Parameters
        ----------
        
        matrix_dimension: int value.
                        Dimension of the matrix that you want to generate. Since the matrix need to be squared, this value represent the number of rows and number of columns.
                        
        eigenvalues_list: array-like, default=None.
                        List of eigenvalues that you want to set in your matrix. This is to facilitate the integer eigenvalues setting for the matrix. If None, the matrix will
                        have its own eigenvalues.
                        
        replicate_paper: bool, default=True.
                        If True means that you are going to replicate the experiments with the matrix proposed in the paper (2x2 or 4x4 based on the dimension_matrix value).
                        If False, you are going to generate a new random matrix.
        
        seed: int value, default=None. 
                        Value to pass if you want a reproducible experiment. If None, every time you execute this method, a different matrix is generated.
        
        Returns
        -------
        input_matrix: array-like. The input matrix that you are going to use in the experiments.
        
        Notes
        -----
        With this method you can generate a new random hermitian matrix or the matrices presented in the paper "A Low Complexity Quantum Principal Component Analysis Algorithm" paper.
        """
        
        if replicate_paper == False:
            if seed!=None:
                np.random.seed(seed)
            random_matrix=np.random.rand(matrix_dimension, matrix_dimension) 
            hermitian_matrix=np.dot(random_matrix, random_matrix.T)

            if eigenvalues_list:
                eig, e_v = np.linalg.eig(hermitian_matrix)
                eigenvalues_list=sorted(eigenvalues_list,reverse=True)
                b = np.array(eigenvalues_list)
                input_matrix = e_v @ np.diag(b) @ e_v.T
            else:
                input_matrix=hermitian_matrix
            
        else:
            if eigenvalues_list:
                warnings.warn("Attention! You want to replicate the paper results so the eigenvalues list that you passed will have no effect!")
            if matrix_dimension==2:
                input_matrix = np.array([[1.5, 0.5],[0.5, 1.5]])
            elif matrix_dimension==4:
                input_matrix = np.array([[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 2, 0], [0, 0, 0, 3]])
            else:
                raise Exception('Only matrix 2x2 or 4x4 are allowed to replicate the reference paper')
        
        print(f'Matrix:\n')
        to_jup(input_matrix)
        
        for eigenval, eigenvec in zip(np.linalg.eig(input_matrix)[0][::-1], np.rot90(np.linalg.eig(input_matrix)[1])):
            
            print(f'eigenvalue: {eigenval} - eigenvector: {eigenvec.round(3)}')

        return input_matrix
    
    