FastCluster
============


## This is simple work flow for clustering molecules

- Author: iwatobipen
- Date: 201712010
- Current version use Morgan FP rad2 is used for clustering

## Requirements

- python3.x
- bayon https://github.com/fujimizu/bayon
- rdkit

## Description

- User need to install bayon at frist, and can found tutorial in following URL. https://github.com/fujimizu/bayon/wiki/Tutorial_English
- Also it is needed to install RDKit for parsing SMILES.
- That's all!

## Basic usage

- input file format is tab delimited text format, "SMILES" \t "ID" \n .....
- $ python fastcluster.py { inputfile } { number of clusters }
- Example usage
- $ ptyhon fast cluster.py cdk2.smi 5 # clustering 47 compounds to 5 clusters.
```
Fastcluster iwatobipen$ python fastcluster.py cdk2.smi 5

real	0m0.015s
user	0m0.006s
sys	0m0.002s
Done!

```


## Output format

- clusterd.tsv is default output format of bayon. List of clusters with similarity points.
```
cluster_1 \t molid1 \t point \t molid2 \t point ... \n
cluster_2 \t molid4 \t point \t molid5 \t point ... \n
....
```

- cluser_parse.tsv is ractangle format of cluster.tsv
```
molid1 \t point \t clusterID1 \n
molid2 \t point \t clusterID2 \n
molid3 \t point \t clusterID3 \n
molid4 \t point \t clusterID4 \n
....
``` 

## Memo

- It will need more cpu time compare with directly using bayon. Because this script convet smiles to fingerprint dataset at first then perform clustering.
