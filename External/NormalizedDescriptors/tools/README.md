# Building the normalized descriptor histograms

`make_normalized_histograms.py` maintains the data behind the
NormalizedDescriptors extension, a port of descriptastorus'
`RDKit2DHistogramNormalized` descriptors.

No distribution is fitted. Each descriptor's raw values over a reference set of
molecules are binned with `numpy.histogram`, exactly as descriptastorus'
`data/d_descriptors/make_histdists.py` does, and a value is normalized by
looking it up in the cumulative bins, exactly as descriptastorus does.

## Table

`../data/normalized_descriptor_histograms.txt`. The first 201 histograms are
imported unchanged from descriptastorus' `descriptastorus/descriptors/hists.py`.
The other 17 (SPS, AvgIpc, Phi, BCUT2D_* and the newer counts) are binned from
descriptastorus' `data/chembl_100k.smi`.

    # comment lines start with #
    histogram <name> <nbins>
    <left edge> <cumulative fraction>
    ... (nbins lines, edges ascending)

There are `min(1000, distinct values)` equal-width bins. The cumulative
fraction on each line is the share of values up to and including that bin.

To normalize `v`: take the cumulative fraction of the first bin whose edge is
`>= v`, or `1.0` when `v` is greater than every edge. In python this is
`bins[bisect.bisect(bins, (v,))][1]`. There is no interpolation and no
clipping.

## Commands

    # check the table against a descriptastorus hists.py, lookup for lookup
    python make_normalized_histograms.py check /path/to/descriptastorus/descriptastorus/descriptors/hists.py

    # RDKit descriptors that have no histogram yet
    python make_normalized_histograms.py missing

    # compute them (one d_<name>.gz per descriptor, the same format as
    # descriptastorus' data/d_descriptors) and bin them into the table
    python make_normalized_histograms.py compute \
        --smiles /path/to/descriptastorus/data/chembl_100k.smi \
        --out samples $(python make_normalized_histograms.py missing)
    python make_normalized_histograms.py histogram samples/d_*.gz

## Rebuilding with the current RDKit

The imported histograms were computed with an older RDKit on a different
molecule set. To rebuild every histogram with the installed RDKit on
`chembl_100k.smi`:

    python make_normalized_histograms.py compute \
        --smiles /path/to/descriptastorus/data/chembl_100k.smi \
        --out samples $(python make_normalized_histograms.py names)
    python make_normalized_histograms.py histogram --replace samples/d_*.gz

`compute` uses `rdkit.Chem.Descriptors._descList`. Of the 218 names, only
`RDKit2D_calculated` (descriptastorus' marker column) is not in it, so that
one is kept as imported. AvgIpc is slow on large peptides, and computing it for
100k molecules takes on the order of 40 CPU-minutes.

Requires numpy; `compute` also needs the RDKit python wrappers.
