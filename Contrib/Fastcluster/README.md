FastCluster
============


## This is simple workflow for clustering molecules

- Author: iwatobipen
- Date: 201712010
- Current version uses Morgan FP; rad2 is used for clustering

## Requirements

- python3.x
- bayon https://github.com/fujimizu/bayon
- rdkit

## Description

- Users need to install bayon at first, and can find tutorial at following URL. https://github.com/fujimizu/bayon/wiki/Tutorial_English
- Also it is needed to install RDKit for parsing SMILES.
- That's all!

## Basic usage

- input file format is tab delimited text format, "ID" \t "SMILES" \n .....
- $ python fastcluster.py {input; inputfile} {N; number of clusters} { --output; filename of output} {--centroid; filename of centroid information} 
- Example usage
- $ ptyhon fastcluster.py cdk2.smi 5 # clustering 47 compounds to 5 clusters.
```
Fastcluster iwatobipen$ python fastcluster.py cdk2.smi 5

real	0m0.015s
user	0m0.006s
sys	0m0.002s

```


## Output format

- clusterd.tsv is default output format of bayon. List of clusters with similarity points.
```
cluster_1 \t molid1 \t point \t molid2 \t point ... \n
cluster_2 \t molid4 \t point \t molid5 \t point ... \n
....
```

- cluser_parse.tsv is rectangle format of cluster.tsv
```
molid1 \t point \t clusterID1 \n
molid2 \t point \t clusterID2 \n
molid3 \t point \t clusterID3 \n
molid4 \t point \t clusterID4 \n
....
``` 

## Memo

- It will need more cpu time compared with directly using bayon. Because this script converts smiles to fingerprint dataset at first then performs clustering.
