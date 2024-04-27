import numpy as np

def diagonal_interpolation(data):
    """ Diagonal interpolation to impute missing bins in Hi-C data. Currently, tested with armatus only
        
        Parameters
        ----------
        ``data`` : np.ndarray
            A Hi-C matrix with empty lines  
        
        Returns
        -------
        np.ndarray
            An imputed Hi-C matrix
    """
    
    if np.any(np.isnan(data) == True):
        bool_matrix = np.isnan(data, dtype = bool)
        zeroidx = np.where(np.all(bool_matrix == True, axis=1))[0]
    else:
        vmin = np.nanmin(data)
        zeroidx = np.where(np.all(data == vmin, axis=1))[0]
    
    if zeroidx.size!= 0:
        data[zeroidx,:] = np.nan
        data[:,zeroidx] = np.nan
        length = np.shape(data)[0]
        diag = data.diagonal(0)
    
        for i in range(0, length):
            diag = data.diagonal(i)
            nan_idx = np.argwhere(np.isnan(diag))
            if nan_idx.size != 0:
                x = np.arange(len(diag))
                x = np.delete(x, nan_idx)
                y = np.delete(diag, nan_idx)
                if x.size == 0:
                    data[nan_idx, nan_idx+i] = 0
                    data[nan_idx+i, nan_idx] = 0
                else:
                    val = np.interp(nan_idx, x, y)
                    data[nan_idx, nan_idx+i] = val
                    data[nan_idx+i, nan_idx] = data[nan_idx, nan_idx+i]
    return data
