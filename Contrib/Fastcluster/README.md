FastCluster
============


## This is simple work flow for clustering molecules

- Author: iwatobipen
- Date: 201712010

## Requirements

- python3.x
- bayon https://github.com/fujimizu/bayon
- rdkit

## Description

- User need to install bayon at frist, and can found tutorial in following URL. https://github.com/fujimizu/bayon/wiki/Tutorial_English
- Also it is needed to install RDKit for parsing SMILES.
- That's all!

## Basic usage

- $ python fastcluster.py { inputfile } { number of clusters }
- Example usage
+ $ ptyhon fast cluster.py cdk2.smi 5 # clustering 47 compounds to 5 clusters.

