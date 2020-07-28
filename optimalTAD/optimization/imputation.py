import numpy as np

def diagonal_interpolation(data):
    zeroidx = np.where(~data.any(axis=1))[0]
    
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
