# GSOC 2022: Integrating xyz2mol Into The RDKit 

### Organization: OpenChemistry
### Student: Sreya Gogineni
### Mentors: Greg Landrum, Joey Storer, Jan H. Jensen

This summer, I worked on integrating 'xyz2mol' into the RDKit, an open source cheminformatics library. '[xyz2mol](https://github.com/jensengroup/xyz2mol)' was originally developed by Professor Jan H. Jensen's research group at the University of Copenhagen, based off of the work published in [this paper](https://onlinelibrary.wiley.com/doi/10.1002/bkcs.10334) (DOI: 10.1002/bkcs.10334). 

The program, given a molecule's charge and the spatial location of each atom, could predict the molecule's most favorable set of bonds.  A user would would pass in the molecule's XYZ file, a file format often used in computational chemistry that delivers each atom's coordinates, and would in return get an RDKit molecule object with predicted bonds in place.

As the original program was written in Python, the nucleus of this project was translation into C++, the language of the RDKit core. 

Integrating xyz2mol into the RDKit required 
- adding an XYZ file parser,
- implementing atomic connectivity determination (knowing which atoms are bonded to each other),
- implementing bond order determination (knowing whether each bond is single, double, or triple), and
- adding Python and Java bindings. 
    
As of the end of the GSOC coding period, the first 3 steps have been completed. The final step, adding bindings to make the features available to RDKit Python and Java users, remains to be finished. 

## The XYZ File Parser

As with other RDKit file parsers (such as the Mol file parser), the XYZ parser constructs an RDKit molecule from the file data. Since the only information an XYZ file contains is the element and location of each atom, the molecule built from the parser contains only atoms and not bonds, as well as a conformer containing the atomic coordinates. The function ```XYZFileToMol()``` calls the file parser.

## Atomic Connectivity Determination

The original xyz2mol offers two methods of predicting connectivity: 'the van der Waals' method and 'Hueckel' method. The former considers atoms' covalent radii to predict bonding, while the Hueckel method uses extended Hueckel theory. 

These two methods were made available through the function ```determineConnectivity()```, which modifies a passed in molecule object in place and adds single bonds wherever a bond is predicted.

## Bond Order Determination

Determining bond order (whether a bond is single, double, or triple) was the largest part of this project. Given a molecule object with bonds corresponding to atomic connectivity, the function, ```determineBondOrdering()``` further modifes the molecule to have a favorable bond ordering. Also added, the function ```determineBonds()``` calls both ```determineConnectivity()``` and ```determineBondOrdering()``` and gives users of the original xyz2mol the ability to use a similar workflow. 

Some interesting tasks while implementing the function included writing an algorithm to calculate the Cartesian product with an arbitrary number of input vectors of arbitrary size and using the Boost graph library.

## Looking Ahead

Through the integration of xyz2mol into the RDKit, its capabilities were made more modular. While the original program did file parsing, connectivity determination, and bond order determination at once, users can now do the three tasks independently of one another, enabling them to potentially swap out atomic connectivity and bond order determination methods or simply read in an XYZ file without using the rest of xyz2mol. 

A lot of progress was made this summer in integrating xyz2mol into the RDKit, but there's yet more work to be done. The first order of business is doing a more comprehensive review of bond order determination for accuracy. This will involve thorough code review and also possibly testing the work with a larger, more diverse set of molecules. And, as mentioned earlier, Python and Java bindings still need to be added. 










