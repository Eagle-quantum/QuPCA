import numpy as np
import warnings
import matplotlib.pyplot as plt
from texttable import Texttable
from scipy.spatial import distance 

def eigenvectors_benchmarking(reconstructed_eigenvalue_eigenvector_tuple, original_eigenvectors, original_eigenvalues, mean_threshold, input_matrix, 
                              n_shots, print_distances=True, only_first_eigenvectors=True, plot_delta=False, distance_type='l2', error_with_sign=False, hide_plot=False):

    """ Method to benchmark the quality of the reconstructed eigenvectors.

    Parameters
    ----------
    reconstructed_eigenvalue_eigenvector_tuple: array-like
                Array of tuples where the first element corresponds to an eigenvalue and the second to the corresponding eigenvector.
    
    original_eigenvectors: array-like
                Array representing the original eigenvectors of the input matrix.
    
    original_eigenvalues: array-like
                Array representing the original eigenvalues of the input matrix.
    
    mean_threshold: array-like
                This array contains the mean between the left and right peaks vertical distance to its neighbouring samples. It is useful for the benchmark process to cut out the bad eigenvalues.
                
    input_matrix: array-like of shape (n_samples, n_features)
                Input hermitian matrix on which you want to apply QPCA divided by its trace, where `n_samples` is the number of samples
                and `n_features` is the number of features.
    
    n_shots: int value
                Number of measures performed in the tomography process.
    
    print_distances: bool value, default=True
            If True, the distance (defined by distance_type value) between the original and reconstructed eigenvector is printed in the legend.

    only_first_eigenvectors: bool value, default=True.
            If True, the function returns only the plot relative to the first eigenvalues. Otherwise, all the plot are showed.

    plot_delta: bool value, default=False
            If True, the function also returns the plot that shows how the tomography error decreases as the number of shots increases. 

    distance_type: string value, default='l2'
            It defines the distance measure used to benchmark the eigenvectors:

                    -'l2': the l2 distance between original and reconstructed eigenvectors is computed.

                    -'cosine': the cosine distance between original and reconstructed eigenvectors is computed.
    
    error_with_sign: bool value, default=False
            If True, the benchmarking is performed considering the reconstructed sign of the eigenvectors. Otherwise, the benchmarking is performed between absolute values of eigenvectors.
    
    hide_plot: bool value, default=False
            If True, the plot showing the eigenvector's benchmarking is not showed.
    

    Returns
    -------
    save_list: array-like
            List of tuples where the first element is the eigenvalue and the second is the distance between the corresponding reconstructed eigenvector and the original one.
            
    delta: float value
            The tomography error value.

    Notes
    -----
    The execution of this method shows the distance between original and reconstructed eigenvector's values and allows to visualize the tomography error. In this way, you can check that the reconstruction of the eigenvectors always takes place with an error conforming to the one expressed in the tomography algorithm in the "A Quantum Interior Point Method for LPs and SDPs" paper.
    """
    
    save_list=[]
    #print('rec',reconstructed_eigenvalue_eigenvector_tuple)
    lambdas_num=[eig[0] for eig in reconstructed_eigenvalue_eigenvector_tuple]
    #print(lambdas_num)
    correct_reconstructed_eigenvalues=__remove_usless_peaks(lambdas_num,mean_threshold,original_eigenvalues)
    #print('correct',correct_reconstructed_eigenvalues)
    correct_reconstructed_eigenvectors=[r[1] for c_r_e in correct_reconstructed_eigenvalues for r in reconstructed_eigenvalue_eigenvector_tuple if c_r_e==r[0]]
    #print(correct_reconstructed_eigenvectors)
    original_eigenValues,original_eigenVectors=__reorder_original_eigenvalues_eigenvectors(input_matrix,original_eigenvectors,original_eigenvalues,correct_reconstructed_eigenvalues)
   
    fig, ax = plt.subplots(1,len(correct_reconstructed_eigenvalues),figsize=(30, 10))
    if len(correct_reconstructed_eigenvalues)>1:
        
        for e,chart in enumerate(ax.reshape(-1,order='F')):
            delta=np.sqrt((36*len(original_eigenVectors[:,e%len(input_matrix)])*np.log(len(original_eigenVectors[:,e%len(input_matrix)])))/(n_shots))

            if error_with_sign==True:

                sign_original=np.sign(original_eigenVectors[:,e%len(input_matrix)])
                sign_original[sign_original==0]=1
                sign_reconstructed=np.sign(correct_reconstructed_eigenvectors[e%len(input_matrix)])
                sign_reconstructed[sign_reconstructed==0]=1
                inverse_sign_original=sign_original*-1       
                sign_difference=(sign_original==sign_reconstructed).sum()
                inverse_sign_difference=(inverse_sign_original==sign_reconstructed).sum()

                if sign_difference>=inverse_sign_difference:
                    original_eigenvector=original_eigenVectors[:,e%len(input_matrix)]
                else:
                    original_eigenvector=original_eigenVectors[:,e%len(input_matrix)]*-1
                reconstructed_eigenvector=correct_reconstructed_eigenvectors[e%len(input_matrix)]
            else:
                original_eigenvector=abs(original_eigenVectors[:,e%len(input_matrix)])
                reconstructed_eigenvector=abs(correct_reconstructed_eigenvectors[e%len(input_matrix)])
            if plot_delta:

                for i in range(len(original_eigenVectors[:,(e%len(input_matrix))])):
                    circle=plt.Circle((i+1,original_eigenvector[i]),np.sqrt(7)*delta,color='g',alpha=0.1)
                    chart.add_patch(circle)
                    chart.axis("equal")
                    chart.hlines(original_eigenvector[i],xmin=i+1,xmax=i+1+(np.sqrt(7)*delta))
                    chart.text(i+1+((i+1+(np.sqrt(7)*delta))-(i+1))/2,original_eigenvector[i]+0.01,r'$\sqrt{7}\delta$')
                chart.plot([], [], ' ', label=r'$\delta$='+str(round(delta,4)))   
                chart.plot(list(range(1,len(input_matrix)+1)),reconstructed_eigenvector,marker='*',label='reconstructed',linestyle='None',markersize=12,alpha=0.5,color='r')
                chart.plot(list(range(1,len(input_matrix)+1)),original_eigenvector,marker='o',label='original',linestyle='None',markersize=12,alpha=0.4)
            else:
                chart.plot(list(range(1,len(input_matrix)+1)),reconstructed_eigenvector,marker='*',label='reconstructed',linestyle='None',markersize=12,alpha=0.5,color='r')
                chart.plot(list(range(1,len(input_matrix)+1)),original_eigenvector,marker='o',label='original',linestyle='None',markersize=12,alpha=0.4)


            if print_distances:

                distance=distance_function_wrapper(distance_type,reconstructed_eigenvector,original_eigenvector)
                chart.plot([], [], ' ', label=distance_type+"_error "+str(np.round(distance,4)))

            save_list.append((correct_reconstructed_eigenvalues[e%len(input_matrix)],np.round(distance,4)))
            chart.plot([], [], ' ', label="n_shots "+str(n_shots))
            chart.legend()
            chart.set_ylabel("eigenvector's values")
            chart.set_title('Eigenvectors corresponding to eigenvalues '+str(correct_reconstructed_eigenvalues[e%len(input_matrix)]))
            if only_first_eigenvectors:
                break
    else:
        
        delta=np.sqrt((36*len(original_eigenVectors[:,0])*np.log(len(original_eigenVectors[:,0])))/(n_shots))

        if error_with_sign==True:

            sign_original=np.sign(original_eigenVectors[:,0])
            sign_original[sign_original==0]=1
            sign_reconstructed=np.sign(correct_reconstructed_eigenvectors[0])
            sign_reconstructed[sign_reconstructed==0]=1
            inverse_sign_original=sign_original*-1       
            sign_difference=(sign_original==sign_reconstructed).sum()
            inverse_sign_difference=(inverse_sign_original==sign_reconstructed).sum()

            if sign_difference>=inverse_sign_difference:
                original_eigenvector=original_eigenVectors[:,0]
            else:
                original_eigenvector=original_eigenVectors[:,0]*-1
            reconstructed_eigenvector=correct_reconstructed_eigenvectors[0]
        else:
            original_eigenvector=abs(original_eigenVectors[:,0])
            reconstructed_eigenvector=abs(correct_reconstructed_eigenvectors[0])
        if plot_delta:

            for i in range(len(original_eigenVectors[:,0])):
                circle=plt.Circle((i+1,original_eigenvector[i]),np.sqrt(7)*delta,color='g',alpha=0.1)
                ax.add_patch(circle)
                ax.axis("equal")
                ax.hlines(original_eigenvector[i],xmin=i+1,xmax=i+1+(np.sqrt(7)*delta))
                ax.text(i+1+((i+1+(np.sqrt(7)*delta))-(i+1))/2,original_eigenvector[i]+0.01,r'$\sqrt{7}\delta$')
            ax.plot([], [], ' ', label=r'$\delta$='+str(round(delta,4)))   
            ax.plot(list(range(1,len(input_matrix)+1)),reconstructed_eigenvector,marker='*',label='reconstructed',linestyle='None',markersize=12,alpha=0.5,color='r')
            ax.plot(list(range(1,len(input_matrix)+1)),original_eigenvector,marker='o',label='original',linestyle='None',markersize=12,alpha=0.4)
        else:
            ax.plot(list(range(1,len(input_matrix)+1)),reconstructed_eigenvector,marker='*',label='reconstructed',linestyle='None',markersize=12,alpha=0.5,color='r')
            ax.plot(list(range(1,len(input_matrix)+1)),original_eigenvector,marker='o',label='original',linestyle='None',markersize=12,alpha=0.4)


        if print_distances:

            distance=distance_function_wrapper(distance_type,reconstructed_eigenvector,original_eigenvector)
            ax.plot([], [], ' ', label=distance_type+"_error "+str(np.round(distance,4)))

        save_list.append((correct_reconstructed_eigenvalues[0],np.round(distance,4)))
        ax.plot([], [], ' ', label="n_shots "+str(n_shots))
        ax.legend()
        ax.set_ylabel("eigenvector's values")
        ax.set_title('Eigenvectors corresponding to eigenvalues '+str(correct_reconstructed_eigenvalues[0]))
    if hide_plot:
        plt.close()
    else:
        fig.tight_layout()
        plt.show()

    return save_list,delta

def eigenvalues_benchmarking(reconstructed_eigenvalue_eigenvector_tuple, original_eigenvalues, mean_threshold, print_error):
    """ Method to benchmark the quality of the reconstructed eigenvalues. 

    Parameters
    ----------
    
    reconstructed_eigenvalue_eigenvector_tuple: array-like
                Array of tuples where the first element corresponds to an eigenvalue and the second to the corresponding eigenvector.
    
    original_eigenvalues: array-like
                Array representing the original eigenvalues of the input matrix.
    
    mean_threshold: array-like
                This array contains the mean between the left and right peaks vertical distance to its neighbouring samples. It is useful for the benchmark process to cut out the bad eigenvalues.
    
    print_error: bool value
                If True, a table showing the absolute reconstruction error for each eigenvalue is reported
    
    Returns
    -------

    Notes
    -----
    """    

    fig, ax = plt.subplots(figsize=(10, 10))
    
    lambdas_num=[eig[0] for eig in reconstructed_eigenvalue_eigenvector_tuple]
    correct_reconstructed_eigenvalues=__remove_usless_peaks(lambdas_num,mean_threshold,original_eigenvalues)
    dict_original_eigenvalues = {original_eigenvalues[e]:e+1 for e in range(len(original_eigenvalues))}
    idx_list=[]
    if print_error:
        t = Texttable()
    for eig in correct_reconstructed_eigenvalues:
        x,min_=__find_nearest(original_eigenvalues,eig)
        idx_list.append(dict_original_eigenvalues[x])
        if print_error:
            error=abs(eig-x)
            lista=[['True eigenvalue','Reconstructed eigenvalue' ,'error'], [x,eig, error]]
            t.add_rows(lista)
    ax.plot(idx_list,correct_reconstructed_eigenvalues,marker='o',label='reconstructed',linestyle='None',markersize=25,alpha=0.3,color='r')
    ax.plot(list(range(1,len(original_eigenvalues)+1)),original_eigenvalues,marker='x',label='original',linestyle='None',markersize=20,color='black')
    ax.legend(labelspacing = 3)
    ax.set_ylabel('Eigenvalue')
    ax.set_title('Matching between original and reconstructed eigenvalues')
    plt.show()
    if print_error:
        print(t.draw())
    

            
def sign_reconstruction_benchmarking(input_matrix, original_eigenvectors, reconstructed_eigenvalue_eigenvector_tuple, original_eigenvalues, mean_threshold, n_shots):
    
    """ Method to benchmark the quality of the sign reconstruction for the reconstructed eigenvectors. 

    Parameters
    ----------
    
    input_matrix: array-like of shape (n_samples, n_features)
                Input hermitian matrix on which you want to apply QPCA divided by its trace, where `n_samples` is the number of samples
                and `n_features` is the number of features.
    
    original_eigenvectors: array-like
                Array representing the original eigenvectors of the input matrix.
    
    reconstructed_eigenvalue_eigenvector_tuple:array-like
                Array of tuples where the first element corresponds to an eigenvalue and the second to the corresponding eigenvector.
    
    original_eigenvalues: array-like
                Array representing the original eigenvalues of the input matrix.
    
    mean_threshold: array-like
                This array contains the mean between the left and right peaks vertical distance to its neighbouring samples. It is useful for the benchmark process to cut out the bad eigenvalues.
    
    n_shots: int value
                The number of measures performed for the reconstruction of the eigenvectors.
    

    Returns
    -------

    Notes
    -----
    """ 
    
    lambdas_num=[eig[0] for eig in reconstructed_eigenvalue_eigenvector_tuple]
    
    correct_reconstructed_eigenvalues=__remove_usless_peaks(lambdas_num,mean_threshold,original_eigenvalues)
    
    correct_reconstructed_eigenvectors=[r[1] for c_r_e in correct_reconstructed_eigenvalues for r in reconstructed_eigenvalue_eigenvector_tuple if c_r_e==r[0]]

    original_eigenValues,original_eigenVectors=__reorder_original_eigenvalues_eigenvectors(input_matrix,original_eigenvectors,original_eigenvalues,correct_reconstructed_eigenvalues)
    
    
    t = Texttable()
    for e in range(len(correct_reconstructed_eigenvalues)):
        eigenvalue=original_eigenValues[e]
        o_eigenvector=original_eigenVectors[:,e]
        r_eigenvector=correct_reconstructed_eigenvectors[e]
        sign_original=np.sign(o_eigenvector)
        sign_original[sign_original==0]=1
        sign_reconstructed=np.sign(r_eigenvector)
        sign_reconstructed[sign_reconstructed==0]=1
        inverse_sign_original=np.sign(o_eigenvector)*-1
        inverse_sign_original[inverse_sign_original==0]=1        
        sign_difference=(sign_original==sign_reconstructed).sum()
        inverse_sign_difference=(inverse_sign_original==sign_reconstructed).sum()
        correct_sign=max(sign_difference,inverse_sign_difference)
        wrong_sign=len(o_eigenvector)-correct_sign
        lista=[['Eigenvalue','n_shots' ,'correct_sign','wrong_sign'], [eigenvalue,n_shots,correct_sign,wrong_sign ]]
        t.add_rows(lista)
    print(t.draw())
    

def error_benchmark(shots_dict, error_dict, dict_original_eigenvalues, delta_list=None, plot_delta=False, label_error='l2'):

    """ Method to benchmark the eigenvector's reconstruction error. The execution of this function shows the trend of the error as the number of shots increases.

    Parameters
    ----------
            
    shots_dict: dict-like
            Dictionary that contains as keys the reconstructed eigenvalues and as values the list of shots for which you are able to reconstruct the corresponding eigenvalue.
    
    error_dict: dict-like
            Dictionary that contains as keys the reconstructed eigenvalues and as values the list of reconstruction errors for the corresponding eigenvectors as the number of shots increases.
    
    dict_original_eigenvalues: dict-like
            Dictionary that contains as key the original eigenvalue and as value its index (ordered position).

    delta_list: array-like, default=None
            List that contains all the tomography error computed for each different number of shots. If None, the plot of the tomography error is not showed. 

    plot_delta: bool value, default=False
            If True, the function also returns the plot that shows how the tomography error decreases as the number of shots increases. 

    label_error: string value, default='l2'
            It defines the distance measure used to benchmark the eigenvectors:

                    -'l2': the l2 distance between original and reconstructed eigenvectors is computed.

                    -'cosine': the cosine distance between original and reconstructed eigenvectors is computed.

    Returns
    -------

    Notes
    -----
    """
 
    if plot_delta==False and delta_list:
        warnings.warn("Attention! delta_list that you passed has actually no effect since the flag plot_delta is set to False. Please set it to True if you want to get the delta decreasing plot.")
    
    fig, ax = plt.subplots(1,len(dict_original_eigenvalues),figsize=(30, 10))
    fig.tight_layout()
    original_eigenvalues=list(dict_original_eigenvalues.keys())
    color = {'blue','red','green'}
    for res,c in zip(error_dict,color):

        e_list=[sub for e_ in error_dict[res] for sub in e_]
       
        dict__ = {k: [v for k1, v in e_list if k1 == k] for k, v in e_list}
        
        tmp_dict={}
        tmp_shots_list={}
        for key,value in dict__.items():
            
            x,min_=__find_nearest(original_eigenvalues,key)
            tmp_dict.setdefault(x, [])
            tmp_shots_list.setdefault(x, [])
            tmp_shots_list[x]+=shots_dict[res][key]
            tmp_dict[x]+=value
        for k,j in tmp_dict.items():
            idx=dict_original_eigenvalues[k]
            ax[idx].plot(tmp_shots_list[k],j,'-o',c=c,label='resolution '+str(res))
            ax[idx].set_xticks(tmp_shots_list[k])
            ax[idx].set_xscale('log')
            ax[idx].set_xlabel('n_shots')
            ax[idx].set_ylabel(label_error+'_error')
            ax[idx].set_title(label_error+'_error for eigenvector wrt the eigenvalues {}'.format(k))
            ax[idx].legend()
        
    plt.show()
    
    
def error_benchmark1(shots_dict,error_dict,delta_list=None,plot_delta=False,label_error='l2'):

    """ Method to benchmark the eigenvector's reconstruction error. The execution of this function shows the trend of the error as the number of shots increases.

    Parameters
    ----------
            
    shots_dict: dict-like
            Dictionary that contains as keys the reconstructed eigenvalues and as values the list of shots for which you are able to reconstruct the corresponding eigenvalue.
    
    error_dict: dict-like
            Dictionary that contains as keys the reconstructed eigenvalues and as values the list of reconstruction errors for the corresponding eigenvectors as the number of shots increases.

    delta_list: array-like, default=None
            List that contains all the tomography error computed for each different number of shots. If None, the plot of the tomography error is not showed. 

    plot_delta: bool value, default=False
            If True, the function also returns the plot that shows how the tomography error decreases as the number of shots increases. 

    label_error: string value, default='l2'
            It defines the distance measure used to benchmark the eigenvectors:

                    -'l2': the l2 distance between original and reconstructed eigenvectors is computed.

                    -'cosine': the cosine distance between original and reconstructed eigenvectors is computed.

    Returns
    -------

    Notes
    -----
    """
 
    if plot_delta==False and delta_list:
        warnings.warn("Attention! delta_list that you passed has actually no effect since the flag plot_delta is set to False. Please set it to True if you want to get the delta decreasing plot.")
    
    fig, ax = plt.subplots(1,len(error_dict),figsize=(30, 10))
    fig.tight_layout()

    for e,chart in enumerate(ax.reshape(-1)):

        for res in error_dict:
            
            e_list=[sub for e_ in error_dict[res] for sub in e_]
            dict__ = {k: [v for k1, v in e_list if k1 == k] for k, v in e_list}
            
            try:
                chart.plot(shots_dict[res][list(dict__.keys())[e]],dict__[list(dict__.keys())[e]],'-o')
                chart.set_xticks(shots_dict[res][list(dict__.keys())[e]])
                chart.set_xscale('log')
                chart.set_xlabel('n_shots')
                chart.set_ylabel(label_error+'_error')
                chart.set_title(label_error+'_error for eigenvector wrt the eigenvalues {}'.format(list(dict__.keys())[e]))
            except:
                continue

    fig, ax = plt.subplots(figsize=(10, 10))
    if plot_delta:
        if delta_list==None:
             raise Exception('You need to provide a delta error list!')
        else:
            shots_list = shots_dict[max(shots_dict, key= lambda x: len(set(shots_dict[x])))]
            plt.plot(shots_list,delta_list,'-o')
            plt.xticks(shots_list)
            plt.xscale('log')
            plt.xlabel('n_shots')
            plt.ylabel(r'$\delta$')
            plt.title(r'Tomography error')

    plt.show()

    
def __reorder_original_eigenvalues_eigenvectors(input_matrix, original_eigenVectors, original_eigenValues, lambdas_num):
    
    check_eigenvalues={}
    new_original_eigenvalues=[]
    new_original_eigenvectors=[]
    
    for i in original_eigenValues:
        check_eigenvalues.update({i:0})
    
    for l_n in lambdas_num:
        x,min_=__find_nearest(original_eigenValues,l_n)
        if check_eigenvalues[x]==0:
            check_eigenvalues.update({x:1})
            new_original_eigenvalues.append(x)
        else:
            continue
    
    for d in list(check_eigenvalues.keys()):
        if check_eigenvalues[d]==1:
            check_eigenvalues.pop(d)
           
    if len(check_eigenvalues)>0:
        new_original_eigenvalues+=list(check_eigenvalues.keys())
    
    
    for i in new_original_eigenvalues:
        for j in range(len(original_eigenValues)):
            if i==original_eigenValues[j]:
                new_original_eigenvectors.append(original_eigenVectors[:,j])
    new_original_eigenvectors=np.array(new_original_eigenvectors)
    new_original_eigenvectors=new_original_eigenvectors.reshape(len(input_matrix),len(input_matrix)).T
    return np.array(new_original_eigenvalues),new_original_eigenvectors

def __remove_usless_peaks(lambdas_num, mean_threshold, original_eig):
    
    check_eigenvalues={}
    peaks_to_keep=[]
    peaks_not_to_keep=[]
    
    for i in original_eig:
            check_eigenvalues.update({i:0}) 
    
    if len(np.array(lambdas_num)[np.argwhere(mean_threshold == np.amin(mean_threshold))])>1:

        not_equal_threshold=np.array(lambdas_num)[np.argwhere(mean_threshold != np.amin(mean_threshold))]
        equal_thresholds=np.array(lambdas_num)[np.argwhere(mean_threshold == np.amin(mean_threshold))]
        #print(not_equal_threshold,equal_thresholds)
        dict_for_equal_threshold={}
        
        for n_e_t in not_equal_threshold:
            #print('net',n_e_t)
            x,min_=__find_nearest(original_eig,n_e_t)
            if check_eigenvalues[x]==0:
                check_eigenvalues.update({x:1})
                peaks_to_keep.append(n_e_t)
            else:
                peaks_not_to_keep.append(n_e_t)
        #print('peaks_after_not_equal_threshold',peaks_to_keep)
        for e_t in equal_thresholds:
            x,min_=__find_nearest(original_eig,e_t)
            tuple_=(e_t,min_)
            dict_for_equal_threshold.setdefault(x, []).append(tuple_)
        #print(dict_for_equal_threshold)
        for d_d in dict_for_equal_threshold:
            minimum_to_keep=min(dict_for_equal_threshold[d_d], key = lambda t: t[1])[0]
            not_minimum_to_discard=[e[0] for e in dict_for_equal_threshold[d_d] if e[0]!= minimum_to_keep]
            #print('minimum_to_keep',minimum_to_keep)
            #print('not_minimum_to_keep',not_minimum_to_discard)
            if check_eigenvalues[d_d]==0:   
                check_eigenvalues.update({d_d:1})
                peaks_to_keep.append(minimum_to_keep)
            else:
            
                peaks_not_to_keep.append(minimum_to_keep)
            peaks_not_to_keep+=not_minimum_to_discard
            #print('peaks_to_keep',peaks_to_keep)
            #print('peaks_not_to_keep',peaks_not_to_keep)
    else:
        for n_p in lambdas_num:
            x,min_=__find_nearest(original_eig,n_p)
            if check_eigenvalues[x]==0:
                check_eigenvalues.update({x:1})
                peaks_to_keep.append(n_p)
            else:
                peaks_not_to_keep.append(n_p)
    if len(peaks_not_to_keep)>0:
        idxs=[list(lambdas_num).index(x) for x in peaks_not_to_keep if x in lambdas_num]
        lambdas_num=np.delete(lambdas_num,idxs)
    return sorted(lambdas_num,reverse=True)

def __find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    min_=(np.abs(array - value)).min()
    return array[idx],min_
    
def distance_function_wrapper(distance_type, *args):
    reconstructed_eig= args[0]
    original_eig = args[1]
    if distance_type=='l2':
        return np.linalg.norm(reconstructed_eig-original_eig)
    if distance_type=='cosine':
        return distance.cosine(reconstructed_eig,original_eig)