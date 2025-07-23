import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parents[3]))
import pandas as pd
import make_array
from featureCall import findFeatures
import numpy as np

df = pd.concat([pd.read_table(f'mESCref_S{i}.bedgraph', header=None, index_col=[0,1,2]) for i in range(1,17)], axis=1)
df = df[ (df.index.get_level_values(2) - df.index.get_level_values(1)) == 50000]
def callIZ(chrom = 'chr1',threshold = 6):
    _df = df.loc[chrom].sort_index(level=2)
    M = np.nan_to_num(_df.values).T
    M = make_array.scale(M)
    M =  make_array.gaussian_smoothing(M)
    M = make_array.scale(M)
    model = findFeatures.cluster_by_birch(M)
    ArgMaxArr = np.array([np.argmax(model.subcluster_centers_, axis=1)[i] * -1 for i in model.labels_[:]])
    features = findFeatures.find_peaks(ArgMaxArr)   
    def get_time_label(S_frac):
        if S_frac <= 2:
            return 'early'
        elif 2 < S_frac <= 5:
            return 'earlymid'
        elif 5 < S_frac <= 8:
            return 'latemid'
        return 'late'
    
    IZcalls = [
        [
            chrom,
            _df.index[f[0]][0],  # start
            _df.index[f[1]][1],  # end
            get_time_label(np.min(np.where(M[:, f[0]] > threshold)[0]))
        ]
        if len(np.where(M[:, f[0]] > threshold)[0]) > 0 
        
        else 
        
        [
            chrom,
            _df.index[f[0]][0],  # start
            _df.index[f[1]][1],  # end
            get_time_label(np.min(np.where(M[:, f[1]+1] > threshold)[0]))
        ] 
        for f in features 
        
    ]
    
    
    return pd.DataFrame(IZcalls)

AutosomeIZs = pd.concat([callIZ(f'chr{i}') for i in range(1, 20)])
AutosomeIZs.to_csv('mESC_F121-9_musAllele_IZ.csv',sep='\t',index=None)
