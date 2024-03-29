import numpy as np
import itertools
import pandas as pd
import time
from scipy.signal import find_peaks
from ..warnings_utils.warning_utility import *

def general_postprocessing(input_matrix, statevector_dictionary, resolution, n_shots, plot_peaks, eigenvalue_threshold, abs_tolerance):
        
        """ Eigenvectors reconstruction process from the reconstructed statevector.
        
        Parameters
        ----------
        
        input_matrix: array-like
                        Hermitian input matrix.
                        
        statevector_dictionary: dict-like.
                        Dictionary where the keys represent the eigenvalues/eigenvectors encoded in the qubits and the values represent the reconstructed statevector's values.
                        
        resolution: int value
                        Number of qubits used in the phase estimation process to represent the eigenvalues.
                        
        n_shots: int value
                        Number of measures performed in the tomography process.
        
        plot_peaks: bool value
                        If True, it returns a plot of the peaks which correspond to the eigenvalues finded by the phase estimation procedure.
                        
        eigenvalue_threshold: float value, default=None
                        It acts as a threshold that cut out the eigenvalues (and the corrseponding eigenvectors) that are smaller than this value.
        
        abs_tolerance: float value, default=None
                        Absolute tolerance parameter used to cut out the eigenvalues estimated badly due to insufficient resolution.
        
        Returns
        -------
        eigenvalue_eigenvector_tuples: array-like 
                        List of tuples containing as first value the reconstructed eigenvalue and as second value the reconstructed eigenvector (sorted from the highest to the lowest eigenvalues).
        
        mean_threshold: array-like
                        This array contains the mean between the left and right peaks vertical distance to its neighbouring samples. It is useful for the benchmark process to cut out the bad eigenvalues.
        
        Notes
        -----
        """
        
        len_input_matrix=len(input_matrix)
        probabilities=np.array([abs(statevector_dictionary[s_d])**2 for s_d in statevector_dictionary])
        bitstrings=[''.join([''.join(str(j)) for j in i]) for i in list(map(list, itertools.product([0, 1], repeat=resolution+len_input_matrix)))]
        
        #tuple to associate at each state (combination of 0 and 1) the corresponding probabilities
        
        probabilities_bitstringed=list(zip(bitstrings, probabilities))
        
        df=pd.DataFrame(probabilities_bitstringed)
        df.columns=['state','module']
        
        #add a "lambda" column to the dataframe that contains all the qubits where the eigenvalues are encoded (as combination of 0 and 1) in the resolution qubits specified in the phase estimation process
        
        df['lambda']=df['state'].apply(lambda x: x[-resolution:])
        
        #aggregate the module for each lambda and sort them in deascending way to see which are the most probable eigenvalues that phase estimation is able to extract with the given resolution
        
        df_aggregated=df.groupby('lambda').agg({'module':'sum'})
        df_aggregated=df_aggregated.sort_values('module',ascending=False)
        
        tail=df_aggregated.reset_index()
        tail['eigenvalue']=tail['lambda'].apply(lambda x :int(x[::-1],base=2)/(2**resolution))
        if plot_peaks:
            ax=tail[['eigenvalue','module']].sort_values('eigenvalue').set_index('eigenvalue').plot(style='-*',figsize=(15,10))
            
            if eigenvalue_threshold:
                ax.axvline(eigenvalue_threshold,ls='--',c='red',label='eigenvalues threshold')
                ax.legend()
            
            
        binary_eigenvalues, mean_threshold=__peaks_extraction(tail,len_input_matrix,n_shots)
        
        if abs_tolerance==None:
            abs_tolerance=1/n_shots
            
        bad_peaks_mask=np.isclose(mean_threshold,0,atol=abs_tolerance)
        
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            if not any(bad_peaks_mask):
                if abs_tolerance==None:
                    customWarning.warn(f'The default tolerance set is {abs_tolerance}. If some output eigenvalues are not the expected ones, it is recommended to increase the absolute tolerance to cut away the noisy eigenvalues.')
                else:
                    customWarning.warn(f'You set an absolute tolerance of {abs_tolerance}. If some output eigenvalues are not the expected ones, it is recommended to increase the absolute tolerance to cut away the noisy eigenvalues.')
                    
                
        binary_eigenvalues=binary_eigenvalues[~bad_peaks_mask]
        mean_threshold=mean_threshold[~bad_peaks_mask]
        
        df.columns=['state','module','lambda']
        
        #add reconstructed sign to the module column
        
        signs=np.sign(np.array(list(statevector_dictionary.values())))
        signs[np.where(signs==0)]=1
        df['sign']=signs
        df['module']=df['module'].multiply(signs, axis=0)
        df=df.fillna(0)
        
        #for each extracted eigenvalue, reconstruct the corresponding eigenvector using the maximum of maxima procedure
        
        statevector_list=[]
        save_sign=[]
        eigenvalues=[]
        idx_to_remove=[]
        
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
        
            for e, bin_eigenvalue in enumerate(binary_eigenvalues):

                #conversion from binary description of the eigenvalue to integer form, remembering that phase estimation encode the eigenvalue as x/2^resolution
                eigenvalue=int(bin_eigenvalue[::-1],base=2)/(2**resolution)
                eigenvalues.append(eigenvalue)
                module_list=np.array(df.query("state.str.endswith(@bin_eigenvalue)")['module'].values)
                save_sign.append(np.sign(module_list))
                statevector_list.append(np.sqrt(abs(module_list)))

                if eigenvalue_threshold and eigenvalue < eigenvalue_threshold:
                    eigenvalues.pop()
                    save_sign.pop()
                    statevector_list.pop()
                    idx_to_remove.append(e)
            mean_threshold=np.delete(mean_threshold,idx_to_remove)

            for i in range(len(statevector_list)):
                normalization_factor=np.sqrt((1/(sum(statevector_list[i]**2))))
                statevector_list[i]*=normalization_factor
                statevector_list[i]*=save_sign[i]    

            eigenvalue_eigenvector_tuples=[]
            for st, eig in zip(statevector_list,eigenvalues):
                eigenvector=np.zeros(len_input_matrix)
                save_sign=np.sign(st)
                statevector=abs(st)
                max_list=[]
                scaled_statevectors=[]
                for e,i in enumerate(range(0,len(statevector),len_input_matrix)):
                    max_list.append(max(statevector[i:i+len_input_matrix]))
                    scaled_statevectors.append(statevector[i:i+len_input_matrix]/max_list[e])

                idx_max=np.argmax(max_list)
                max_max=max_list[idx_max]
                value=np.sqrt(max_max)
                eigenvector=scaled_statevectors[idx_max]*value*save_sign[len_input_matrix*idx_max:len_input_matrix*idx_max+len_input_matrix]
                eigenvalue_eigenvector_tuples.append((eig,eigenvector))
        
        return eigenvalue_eigenvector_tuples,mean_threshold
    

def __peaks_extraction(df,len_input_matrix,n_shots):
    """ 
    Process to find the correct eigenvalues peaks after the PE procedure.
        
        Parameters
        ----------
        
        df: DataFrame
                        It contains the extracted eigenvalues ordered by module squared.
        
        len_input_matrix: int value
                        Length of the input matrix.
        
        n_shots: int value
                        Number of measures performed in the tomography process.
        
        Returns
        -------
        sorted_peaks: list-like
                        List of peaks in binary format that represent eigenvalues extracted by PE.
                        
        mean_threshold: ndarray-like
                        This array contains the mean between the left and right peaks vertical distance to its neighbouring samples.
        
        
        Notes
        -----
        This algorithm uses the find_peaks SciPy function to extract the right peaks (which correspond to the eigenvalues) of the PE outputs. This is an iterative function that search for the right number of peaks. If an extracted peak is not the right one, it means that the PE resolution needs to be increased to correctly identify the right peaks.
    """
    
    peaks=[]
    offset=1/n_shots
    stop=False
    start=time.time()
    
    with warnings.catch_warnings():
        warnings.simplefilter('always')
        
        while stop==False:
            
            #check if the peaks finder is enter in an infinite loop

            if time.time()-start>30:
                customWarning.warn("The extraction of the eigenvalues is taking a long time. You may have hit a plateau and therefore you may need to restart the execution by increasing the number of resolution qubits and/or the number of measurements performed.")
                warnings.simplefilter('ignore')
                
            #compute the most likely peaks (eigenvalues)

            p_=find_peaks(df.sort_values(['eigenvalue'])['module'],threshold=offset)
            p=p_[0]
            right_thresholds=p_[1]['right_thresholds']
            left_thresholds=p_[1]['left_thresholds']

            #check if the number of peaks are equal (or less) then the expected ones since we don't want a greater number of peaks than expected

            if len(p)>len_input_matrix or len(p)==0:
                offset+=1/n_shots
            else:
                stop=True
                for i in p:
                    el = df.sort_values(['eigenvalue']).iloc[i]
                    peaks.append(el['lambda'])
    
    #mean_threshold helps in showing which are the eigenvalue that we are not able to estimate in a right way. If the phase estimation is not able to extract an eigenvalue in a correct way, this eigenvalue will have the lowest mean_threshold value
    
    mean_threshold=(left_thresholds+right_thresholds)/2
    sorted_peaks=np.array(peaks)[mean_threshold.argsort()[::-1]]
    return sorted_peaks,np.array(sorted(mean_threshold,reverse=True))
    
