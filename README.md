# High Resolution Repli-seq
## overview
This directory contains codes used in [Zhao, Sasaki 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-01983-8) for processing high resolution repli-seq data by constructing column-scaled, Gaussian smoothed Repli-Seq data arrays and replication feature extraction (initiation zones, timing transition zones, etc) from these arrays.

Example: mouse embryonic stem cell chr1 map as below:
![alt text](https://github.com/oliviacamel/High-Resolution-RepliSeq/blob/master/doc/fig1.png)
## data
Here, I include high-res repli-seq data of
- human embryonic stem cell lines H1 and H9
- human colon cancer cell line HCT116
- mouse embryonic stem cell line F121-9 (*mus* x *cas* hybrid)
- mouse neural precursor cell differentiated from F121-9 (*mus* x *cas* hybrid)

## analysis outline
### 1. preprocessing and noramlising replication array
  - loading Repli-Seq data from bedgraph files, one bedgraph per S phase fraction

  - Gaussian smoothing and column-wise scaling
```python
import HighResRepliSeq
import pandas as pd
#assuming the current directory contains bedgraph files, one for each S phase fraction, as follows:
#└── musAllele
#    ├── mESCref_S1.bedgraph
#    ├── mESCref_S10.bedgraph
#    ├── mESCref_S11.bedgraph
#    ├── mESCref_S12.bedgraph
#    ├── mESCref_S13.bedgraph
#    ├── mESCref_S14.bedgraph
#    ├── mESCref_S15.bedgraph
#    ├── mESCref_S16.bedgraph
#    ├── mESCref_S2.bedgraph
#    ├── mESCref_S4.bedgraph
#    ├── mESCref_S5.bedgraph
#    ├── mESCref_S6.bedgraph
#    ├── mESCref_S7.bedgraph
#    ├── mESCref_S8.bedgraph
#    ├── mESCref_S9.bedgraph
df = pd.concat([pd.read_table(f'mESCref_S{i}.bedgraph', header=None, index_col=[0,1,2]) for i in range(1,17)], axis=1)
df = df[ (df.index.get_level_values(2) - df.index.get_level_values(1)) == 50000]
processed_df = []
for chrom in [f'chr{i}' for i in range(1,20)]: #we are only looking at autosomes here
  _df = df.loc[chrom].sort_index(level=2)
  M = np.nan_to_num(_df.values).T
  M = HighResRepliSeq.make_array.scale(M)
  M = HighResRepliSeq.make_array.gaussian_smoothing(M)
  M = HighResRepliSeq.make_array.scale(M)
  _df = pd.DataFrame(M.T,index = pd.MultiIndex.from_tuples(
        [(chrom, start, end) for start, end in _df.index],
        names=['chrom', 'start', 'end']), columns=[f'S{i}' for i in range(1, 17)]
    )
  processed_df.append(_df)
df = pd.concat(processed_df, axis=0)
df.to_csv('mESCrefAllele_highresrepliseq_processedArray.csv')
```
### 2. call initiation zones (IZs), termination zones (TZs) and timing transition regions (TTRs) from data
```python
import numpy as np
chrom = 'chr1'
threshold = 6 #define the minimum percentage of replication to assign timing e.g. an IZ needs at least 6%, the threshold defined here, of replication in the first quarter of S phase to be assigned as an 'early IZ'.
#use birch clustering to cluster genomic bins according to their replication profile
_df = df.loc[chrom]
Arr = _df.values.T
model = HighResRepliSeq.findFeatures.cluster_by_birch(Arr)
#make array of clustering result
ArgMaxArr = np.array([np.argmax(model.subcluster_centers_, axis=1)[i] * -1 for i in model.labels_[:]])
#find IZs
features = HighResRepliSeq.findFeatures.find_peaks(ArgMaxArr)  
# find TZs features = HighResRepliSeq.findFeatures.find_valleys(ArgMaxArr)  
# find TTRs features = HighResRepliSeq.findFeatures.find_slopes(ArgMaxArr,direction = 'right') accepted direction:{'right','left'}
calls = pd.DataFrame([
    [
        chrom,
        _df.index[f[0]][0],  # start
        _df.index[f[1]][1],  # end
        HighResRepliSeq.findFeatures.get_time_label(np.min(np.where(Arr[:, f[0]] > threshold)[0]))
    ]
    for f in features
    if len(np.where(Arr[:, f[0]] > threshold)[0]) > 0
])
```

## command line usage examples
