This directory contains scripts for carrying out the calculations and analyis from the publications:

1) A. Vulpetti, U. Hommel, G. Landrum, R. Lewis and C. Dalvit, "Design and NMR-based screening of LEF,
a library of chemical fragments with different Local Environment of Fluorine" J. Am. Chem. Soc. 131 (2009) 12949-12959. https://doi.org/10.1021/ja905207t

2) A. Vulpetti, G. Landrum, S. Ruedisser, P. Erbel and C. Dalvit, "19F NMR Chemical Shift Prediction with 
Fluorine Fingerprint Descriptor" J. of Fluorine Chemistry (2010). https://doi.org/10.1016/j.jfluchem.2009.12.024

The scripts require that the RDKit (www.rdkit.org) be installed and properly configured.

The scripts assume that input db.sdf files only contain molecules having either CF or CF3 moieties. 

Commands to run:

1) Generation of Fluorine Fingerprint (F-FP-L)
------------------------------------------------------------------------------------------
	 python CreateFps.py db.sdf db.layers.pkl > dupes.layers.txt

2) Butina Clustering using Fluorine Fingerprint (F-FP-L)
------------------------------------------------------------------------------------------
	 python ClusterFps.py db.layers.pkl > clusters.layers.txt


3) Cliff Analysis using Fluorine Fingerprint Similarity vs a specified property in propField
-------------------------------------------------------------------------------------------
	python DistancePlot.py db.sdf cliff.txt
	
4) KNN prediction of a property specified in propField using F-FP-L (L= maximum path length in atoms)
-------------------------------------------------------------------------------------------
	python DistancePredict.py --max=L+1 --sim="[0.9,0.8,0.7]" --nbrs=nbrs.txt training.sdf
	test.sdf prediction.txt


In the event you use the scripts for publication please reference the original publications:

1) A. Vulpetti, U. Hommel, G. Landrum, R. Lewis and C. Dalvit, "Design and NMR-based screening of LEF,
a library of chemical fragments with different Local Environment of Fluorine" J. Am. Chem. Soc. 131 (2009) 12949-12959. https://doi.org/10.1021/ja905207t

2) A. Vulpetti, G. Landrum, S. Ruedisser, P. Erbel and C. Dalvit, "19F NMR Chemical Shift Prediction with 
Fluorine Fingerprint Descriptor" J. of Fluorine Chemistry (2010). https://doi.org/10.1016/j.jfluchem.2009.12.024

