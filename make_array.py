import numpy as np
def gaussian_smoothing (M,kernel_shape=(3,3),sigma=1):
    '''
    Gaussian smoothing for 16 fraction high resolution repli-seq data
    
    

    Parameters
    ----------
    M : np.ndarray
        Input 2D array, 16 rows, S1-S16.
    kernel_shape : tuple, optional
        Shape of the Gaussian kernel. The default is (3,3).
    sigma : float, optional
        Sigma for Gaussian kernel. The default is 1.

    Returns
    ----------
    newM: np.ndarray
        smoothed repli-seq matrix
        

    '''
    def _gaussian_kernel(shape=kernel_shape,sigma=sigma):
        m,n=[(edge-1)/2 for edge in shape]
        y,x = np.ogrid[-m:m+1,-n:n+1]
        kernel = np.exp(-(x**2 + y**2) / (2 * sigma**2))
        kernel[kernel < np.finfo(kernel.dtype).eps * kernel.max()] = 0
        kernel /= kernel.sum() 
        return kernel
    
    newM = np.zeros_like (M)
    
    _M=np.concatenate((np.array([M[0,:] for i in range((kernel_shape[0]-1)//2)]),M,np.array([M[-1,:] for i in range((kernel_shape[0]-1)//2)])))
    
    _M=np.pad(_M,((0,0),((kernel_shape[0]-1)//2,(kernel_shape[0]-1)//2)),'constant',constant_values=np.nan)
    
    
    
    for i in range((kernel_shape[0] - 1) // 2, _M.shape[0] - (kernel_shape[0] - 1) // 2):
        
        for j in range((kernel_shape[1] - 1) // 2, _M.shape[1] - (kernel_shape[0] - 1) // 2):
            
            box=np.ma.masked_invalid(_M[i-(kernel_shape[0]-1)//2:i+(kernel_shape[0]-1)//2+1,j-(kernel_shape[1]-1)//2:j+(kernel_shape[1]-1)//2+1])
            
            newM[i-(kernel_shape[0]-1)//2,j-(kernel_shape[1]-1)//2] = np.nansum(np.multiply(box,_gaussian_kernel()))
    
    return (newM)


def scale(M):
    
    '''
    This scaling converts Gaussian smoothed values to percentages for each genomic bin across the 16 S-phase fractions (S1-S16).
    
    Parameters
    ----------
    M : np.ndarray
        Gaussian smoothed 2D array, 16 rows, S1-S16.
    Returns
    ----------
    np.ndarray
    
    
    
    '''
        

    col_sums = np.nansum(M, axis=0, keepdims=True)
    
    col_sums[col_sums == 0] = np.nan 
    scaled_M = (M / col_sums) * 100
    
    return np.nan_to_num(scaled_M,nan = 0.0)

