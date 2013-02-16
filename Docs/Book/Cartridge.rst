The RDKit database cartridge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

What is this?
+++++++++++++

This document is a tutorial and reference guide for the RDKit PostgreSQL cartridge.

If you find mistakes, or have suggestions for improvements, please
either fix them yourselves in the source document (the .rst file) or
send them to the mailing list: rdkit-discuss@lists.sourceforge.net 
(you will need to subscribe first)

Tutorial
++++++++

Introduction
************
 
In this example I show how to load a database from the SMILES file of
commercially available compounds that is downloadable from
emolecules.com at URL
http://www.emolecules.com/doc/plus/download-database.php

If you choose to repeat this exact example yourself, please note that
it takes several hours to load the 6 million row database and generate
all fingerprints. The timing information below was collected on a
commodity desktop PC (Dell Studio XPS with a 2.9GHz i7 CPU and 8GB of
RAM) running Ubuntu 12.04 and using PostgreSQL v9.1.4. The database
was installed with default parameters.

To improve performance while loading the database and building the index, 
I changed a couple of postgres configuration settings in `postgresql.conf` ::

  fsync = off				# turns forced synchronization on or off
  synchronous_commit = off		# immediate fsync at commit
  full_page_writes = off			# recover from partial page writes

And to improve search performance, I allowed postgresql to use more memory than the
extremely conservative default settings::

  shared_buffers = 2048MB			# min 128kB
  					# (change requires restart)
  work_mem = 128MB				# min 64kB




Creating the database
*********************

First create the database and install the cartridge::

  ~/RDKit_trunk/Data/emolecules > createdb emolecules
  ~/RDKit_trunk/Data/emolecules > psql -c 'create extension rdkit' emolecules 


Now create and populate a table holding the raw data::

  ~/RDKit_trunk/Data/emolecules > psql -c 'create table raw_data (id SERIAL, smiles text, emol_id integer, parent_id integer)' emolecules
  NOTICE:  CREATE TABLE will create implicit sequence "raw_data_id_seq" for serial column "raw_data.id"
  CREATE TABLE
  ~/RDKit_trunk/Data/emolecules > zcat emolecules-2013-02-01.smi.gz | sed '1d; s/\\/\\\\/g' | psql -c "copy raw_data (smiles,emol_id,parent_id) from stdin with delimiter ' '" emolecules


Create the molecule table, but only for SMILES that the RDKit accepts::

  ~/RDKit_trunk/Data/emolecules > psql emolecules
  psql (9.1.4)
  Type "help" for help.
  emolecules=# select * into mols from (select id,mol_from_smiles(smiles::cstring) m from raw_data) tmp where m is not null;
  WARNING:  could not create molecule from SMILES 'CN(C)C(=[N+](C)C)Cl.F[P-](F)(F)(F)(F)F'
  ... a lot of warnings deleted ...
  SELECT 6008732
  emolecules=# create index molidx on mols using gist(m);
  CREATE INDEX

The last step is only required if you plan to do substructure searches.

Substructure searches
*********************

Example query molecules taken from the `eMolecules home page <http://www.emolecules.com/>`_::

    emolecules=# select count(*) from mols where m@>'c1cccc2c1nncc2' ;
     count 
    -------
       1593
    (1 row)
    
    Time: 3413.018 ms
    emolecules=# select count(*) from mols where m@>'c1ccnc2c1nccn2' ;
     count 
    -------
      3692
    (1 row)

    Time: 760.860 ms
    emolecules=# select count(*) from mols where m@>'c1cncc2n1ccn2' ;
     count 
    -------
      2359
    (1 row)

    Time: 790.864 ms
    emolecules=# select count(*) from mols where m@>'Nc1ncnc(N)n1' ;
     count 
    -------
     14086
    (1 row)

    Time: 2445.430 ms

Notice that the last query is starting to take a while to execute and count all the results. 
This is even more extreme with the next few queries::

    emolecules=# select count(*) from mols where m@>'c1scnn1' ;
     count 
    -------
     108477
    (1 row)
    
    Time: 37925.126 ms
    emolecules=# select count(*) from mols where m@>'c1cccc2c1CNCCN2' ;
     count 
    -------
      2490
    (1 row)
    
    Time: 46126.816 ms
    emolecules=# select count(*) from mols where m@>'c1cccc2c1ncs2' ;
     count 
    -------
     104895
    (1 row)
    
    Time: 77505.272 ms

Given we're searching through 6 million compounds these search times aren't incredibly slow, 
but it would be nice to have them quicker.

One easy way to speed things up, particularly for queries that return a large number of results, is to only 
retrieve a limited number of results::

    emolecules=# select * from mols where m@>'c1cccc2c1ncs2' limit 100 ;
       id    |                                m                                
    ---------+-----------------------------------------------------------------
     5273717 | OC1CC(Nc2nc3ccccc3s2)C1
     5278926 | [I-].CC[n+]1c(/C=C/Nc2ccccc2)sc2ccccc21
     5282075 | COC(=O)c1ccc2nc(Br)sc2c1
     5283354 | CCc1ccc2nc(N(C)CC(=O)O)sc2c1
     5283355 | Cc1ccc2nc(N(C)CC(=O)O)sc2c1
     5283356 | COc1ccc2nc(N(C)CC(=O)O)sc2c1
     5283357 | CCOc1ccc2nc(N(C)CC(=O)O)sc2c1
     ...
     4854425 | NC(=O)c1ccccc1NC(=O)C1CN(c2nc3c(cccc3F)s2)C1
    (100 rows)

    Time: 50.644 ms

SMARTS-based queries
********************

Oxadiazole or thiadiazole::

    emolecules=# select * from mols where m@>'c1[o,s]ncn1'::qmol limit 500;
       id    |                               m                               
    ---------+---------------------------------------------------------------
     5273135 | Cc1nsc(Br)n1
     5284275 | CCCC[Sn](CCCC)(CCCC)c1nc(C)ns1
     5192275 | CCCCCC(CC(=O)OCC)OC(=O)COCc1nc(C)no1
     5188130 | O=c1c2cccnc2ncn1Cc1nc(-c2ccoc2)no1
     5188272 | COCCc1noc(CNCC2CCCN2c2cccnn2)n1
     5188249 | Cc1nc(CN2CCCC(Nc3cc(C)nc4ncnn43)C2)no1
     5188283 | CN(Cc1nc(-c2ccco2)no1)CC1CCCN1c1cccnn1
     5188293 | COCCc1noc(CN(C)CC2CCCN2c2cccnn2)n1
     ...
     5037294 | Cc1noc(COc2cccc([N+](=O)[O-])c2C)n1
    (500 rows)
     
    Time: 313.202 ms

Notice that this is slower than the the pure SMILES query, this is generally true of SMARTS-based queries.

Similarity searches
*******************

Generating fingerprints and indices::

    emolecules=# select id,torsionbv_fp(m) as torsionbv,morganbv_fp(m,2) as mfp2 into fps from mols;
    SELECT 6008732
    Time: 1734537.410 ms
    emolecules=# create index mfp2idx on fps using gist(mfp2);
    CREATE INDEX
    Time: 381025.418 ms
    emolecules=# create index torsionbvidx on fps using gist(torsionbv);
    CREATE INDEX
    Time: 379285.670 ms
    emolecules=# alter table mols add primary key (id);
    alter table fps add foreign key (id) references mols;NOTICE:  ALTER TABLE / ADD PRIMARY KEY will create implicit index "mols_pkey" for table "mols"
    ALTER TABLE
    Time: 50798.813 ms
    emolecules=# alter table fps add foreign key (id) references mols;
    ALTER TABLE
    Time: 39067.348 ms

Basic similarity searching::

    emolecules=# select count(*) from fps where mfp2%morganbv_fp('Cc1ccc2nc(-c3ccc(NC(C4N(C(c5cccs5)=O)CCC4)=O)cc3)sc2c1',2);
     count 
    -------
       513
    (1 row)

    Time: 4044.707 ms

Usually we'd like to find a sorted listed of neighbors along with the accompanying SMILES. 
This SQL function makes that pattern easy::

    emolecules=# create or replace function get_mfp2_neighbors(smiles text)
      returns table(molregno int, m mol, similarity double precision) as
    $$
    select id,m,tanimoto_sml(morganbv_fp($1::mol),mfp2) as similarity
    from fps join mols using (id) 
    where morganbv_fp($1::mol)%mfp2 
    order by morganbv_fp($1::mol)<%>mfp2;
    $$ language sql stable ;
    CREATE FUNCTION
    Time: 600.371 ms
    emolecules=# select * from get_mfp2_neighbors('Cc1ccc2nc(-c3ccc(NC(C4N(C(c5cccs5)=O)CCC4)=O)cc3)sc2c1') limit 10;
     molregno |                             m                             |    similarity     
    ----------+-----------------------------------------------------------+-------------------
      3116265 | Cc1ccc2nc(-c3ccc(NC(=O)[C@@H]4CCCN4C(=O)c4cccs4)cc3)sc2c1 |                 1
      1598902 | Cc1ccc2sc(-c3ccc(NC(=O)C4CCCN4C(=O)c4cccs4)cc3)nc2c1      | 0.888888888888889
      3118194 | O=C(Nc1ccc(-c2nc3ccccc3s2)cc1)[C@@H]1CCCN1C(=O)c1cccs1    |          0.796875
      5695374 | Cc1ccc2nc(NC(=O)C3CCCN3C(=O)c3cccs3)sc2c1                 | 0.777777777777778
      1758570 | Cc1ccc2nc(-c3ccc(NC(=O)C4CCN(C(=O)c5cccs5)CC4)cc3)sc2c1   | 0.772727272727273
      4267350 | Cc1nc2ccc(NC(=O)[C@@H]3CCCN3C(=O)c3cccs3)cc2s1            | 0.738461538461539
      5825487 | Cc1ccc(NC(=O)C2CCCCN2C(=O)c2cccs2)cc1                     | 0.733333333333333
      2682124 | Cc1ccc2nc(-c3ccc(NC(=O)C4CCCN4S(C)(=O)=O)cc3)sc2c1        | 0.701492537313433
      3552075 | Cc1ccc2nc(-c3ccc(NC(=O)C4CCCCN4S(C)(=O)=O)cc3)sc2c1       | 0.686567164179104
      1807011 | CSc1nc2ccc(NC(=O)C3CCCN3C(=O)c3cccs3)cc2s1                | 0.671428571428571
    (10 rows)

    Time: 4156.841 ms
    emolecules=# select * from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1') limit 10;
     molregno |                  m                   |    similarity     
    ----------+--------------------------------------+-------------------
      5283355 | Cc1ccc2nc(N(C)CC(=O)O)sc2c1          |                 1
      5283354 | CCc1ccc2nc(N(C)CC(=O)O)sc2c1         | 0.761904761904762
      5283360 | CN(CC(=O)O)c1nc2ccc(Br)cc2s1         |  0.75609756097561
      5283363 | CN(CC(=O)O)c1nc2ccc(F)cc2s1          | 0.738095238095238
      5283369 | CN(CC(=O)O)c1nc2ccc(Cl)cc2s1         | 0.738095238095238
      5283365 | Cc1cc2nc(N(C)CC(=O)O)sc2cc1C         |             0.725
      5283367 | CN(CC(=O)O)c1nc2ccc(S(C)(=O)=O)cc2s1 | 0.720930232558139
      5283356 | COc1ccc2nc(N(C)CC(=O)O)sc2c1         | 0.704545454545455
      5283362 | CC(C)c1ccc2nc(N(C)CC(=O)O)sc2c1      | 0.704545454545455
      5283358 | CSc1ccc2nc(N(C)CC(=O)O)sc2c1         | 0.704545454545455
    (10 rows)

    Time: 4186.420 ms




Reference Guide
+++++++++++++++

New Types
*********

* `mol` : an rdkit molecule. Can be created from a SMILES via direct type conversion, for example: `'c1ccccc1'::mol` creates a molecule from the SMILES `'c1ccccc1'`
* `qmol` : an rdkit molecule containing query features (i.e. constructed from SMARTS). Can be created from a SMARTS via direct type conversion, for example: `'c1cccc[c,n]1'::qmol` creates a query molecule from the SMARTS `'c1cccc[c,n]1'`
* `sfp` : a sparse count vector fingerprint (`SparseIntVect` in C++ and Python)
* `bfp` : a bit vector fingerprint (`ExplicitBitVect` in C++ and Python)

Parameters
**********

* `rdkit.tanimoto_threshold` : threshold value for the Tanimoto similarity operator. Searches done using Tanimoto similarity will only return results with a similarity of at least this value.
* `rdkit.dice_threshold` : threshold value for the Dice similiarty operator. Searches done using Dice similarity will only return results with a similarity of at least this value.

Operators
*********

Similarity search
-----------------

* `%` : operator used for similarity searches using Tanimoto similarity. Returns whether or not the Tanimoto similarity between two fingerprints (either two `sfp` or two `bfp` values) exceeds `rdkit.tanimoto_threshold`.
* `#` : operator used for similarity searches using Dice similarity. Returns whether or not the Dice similarity between two fingerprints (either two `sfp` or two `bfp` values) exceeds `rdkit.dice_threshold`.
* `<%>` : used for Tanimoto KNN searches (to return ordered lists of neighbors).
* `<#>` : used for Dice KNN searches (to return ordered lists of neighbors).


Substructure and exact structure search
---------------------------------------

* `@>` : substructure search operator. Returns whether or not the `mol` or `qmol` on the right is a substructure of the `mol` on the left.
* `<@` : substructure search operator. Returns whether or not the `mol` or `qmol` on the left is a substructure of the `mol` on the right.
* `@=` : returns whether or not two molecules are the same.

Molecule comparison
-------------------

* `<` : returns whether or not the left `mol` is less than the right `mol`
* `>` : returns whether or not the left `mol` is greater than the right `mol`
* `=` : returns whether or not the left `mol` is equal to the right `mol`
* `<=` : returns whether or not the left `mol` is less than or equal to the right `mol`
* `>=` : returns whether or not the left `mol` is greater than or equal to the right `mol`

*Note* Two molecules are compared by making the following comparisons in order. Later comparisons are only made if the preceding values are equal:

# Number of atoms
# Number of bonds
# Molecular weight
# Number of rings

If all of the above are the same and the second molecule is a substructure of the first, the molecules are declared equal, Otherwise (should not happen) the first molecule is arbitrarily defined to be less than the second.

There are additional operators defined in the cartridge, but these are used for internal purposes.

Functions
*********

Fingerprint Related
-------------------

Generating fingerprints
:::::::::::::::::::::::

* `morgan_fp(mol,int)` : returns an `sfp` which is the count-based Morgan fingerprint for a molecule using connectivity invariants. The second argument provides the radius. This is an ECFP-like fingerprint.
* `morganbv_fp(mol,int)` : returns a `bfp` which is the bit vector Morgan fingerprint for a molecule using connectivity invariants. The second argument provides the radius. This is an ECFP-like fingerprint.
* `featmorgan_fp(mol,int)` : returns an `sfp` which is the count-based Morgan fingerprint for a molecule using chemical-feature invariants. The second argument provides the radius. This is an FCFP-like fingerprint.
* `featmorganbv_fp(mol,int)` : returns a `bfp` which is the bit vector Morgan fingerprint for a molecule using chemical-feature invariants. The second argument provides the radius. This is an FCFP-like fingerprint.
* `rdkit_fp(mol)` : returns a `bfp` which is the RDKit fingerprint for a molecule. This is a daylight-fingerprint using hashed molecular subgraphs.
* `atompair_fp(mol)` : returns an `sfp` which is the count-based atom-pair fingerprint for a molecule.
* `atompairbv_fp(mol)` : returns a `bfp` which is the bit vector atom-pair fingerprint for a molecule.
* `torsion_fp(mol)` : returns an `sfp` which is the count-based topological-torsion fingerprint for a molecule.
* `torsionbv_fp(mol)` : returns a `bfp` which is the bit vector topological-torsion fingerprint for a molecule.
* `layered_fp(mol)` : returns a `bfp` which is the layered fingerprint for a molecule. This is an experimental substructure fingerprint using hashed molecular subgraphs.
* `maccs_fp(mol)` : returns a `bfp` which is the MACCS fingerpring for a molecule (*available from 2013_01 release*).

Working with fingerprints
:::::::::::::::::::::::::

* `tanimoto_sml(fp,fp)` : returns the Tanimoto similarity between two fingerprints of the same type (either two `sfp` or two `bfp` values).
* `dice_sml(fp,fp)` : returns the Dice similarity between two fingerprints of the same type (either two `sfp` or two `bfp` values).
* `size(bfp)` : returns the length of (number of bits in) a `bfp`.
* `add(sfp,sfp)` : returns an `sfp` formed by the element-wise addition of the two `sfp` arguments.
* `subtract(sfp,sfp)` : returns an `sfp` formed by the element-wise subtraction of the two `sfp` arguments.
* `all_values_lt(sfp,int)` : returns a boolean indicating whether or not all elements of the `sfp` argument are less than the `int` argument.
* `all_values_gt(sfp,int)` : returns a boolean indicating whether or not all elements of the `sfp` argument are greater than the `int` argument.

Fingerprint I/O
:::::::::::::::

* `bfp_to_binary_text(bfp)` : returns a bytea with the binary string representation of the fingerprint that can be converted back into an RDKit fingerprint in other software. (*available from Q3 2012 (2012_09) release*)
* `bfp_from_binary_text(bytea)` : constructs a bfp from a binary string representation of the fingerprint. (*available from Q3 2012 (2012_09) release*)


Molecule Related
----------------

Molecule I/O and Validation
::::::::::::::::::::::::::::::::::::

* `is_valid_smiles(smiles)` : returns whether or not a SMILES string produces a valid RDKit molecule.
* `is_valid_ctab(ctab)` : returns whether or not a CTAB (mol block) string produces a valid RDKit molecule.
* `is_valid_smarts(smarts)` : returns whether or not a SMARTS string produces a valid RDKit molecule.
* `is_valid_mol_pkl(bytea)` : returns whether or not a binary string (bytea) can be converted into an RDKit molecule. (*available from Q3 2012 (2012_09) release*)

* `mol_from_smiles(smiles)` : returns a molecule for a SMILES string, NULL if the molecule construction fails.
* `mol_from_smarts(smarts)` : returns a molecule for a SMARTS string, NULL if the molecule construction fails.
* `mol_from_ctab(ctab)` : returns a molecule for a CTAB (mol block) string, NULL if the molecule construction fails.
* `mol_from_pkl(bytea)` : returns a molecule for a binary string (bytea), NULL if the molecule construction fails. (*available from Q3 2012 (2012_09) release*)

* `mol_to_smiles(mol)` : returns the canonical SMILES for a molecule.
* `mol_to_smarts(mol)` : returns SMARTS string for a molecule.
* `mol_to_pkl(mol)` : returns binary string (bytea) for a molecule. (*available from Q3 2012 (2012_09) release*)

Descriptors
:::::::::::

* `mol_amw(mol)` : returns the AMW for a molecule.
* `mol_logp(mol)` : returns the MolLogP for a molecule.
* `mol_tpsa(mol)` : returns the topological polar surface area for a molecule (*available from Q1 2011 (2011_03) release*).
* `mol_fractioncsp3(mol)` : returns the fraction of carbons that are sp3 hybridized (*available from 2013_03 release*).
* `mol_hba(mol)` : returns the number of Lipinski H-bond acceptors (i.e. number of Os and Ns) for a molecule.
* `mol_hbd(mol)` : returns the number of Lipinski H-bond donors (i.e. number of Os and Ns that have at least one H) for a molecule.
* `mol_numatoms(mol)` : returns the total number of atoms in a molecule.
* `mol_numheavyatoms(mol)` : returns the number of heavy atoms in a molecule.
* `mol_numrotatablebonds(mol)` : returns the number of rotatable bonds in a molecule (*available from Q1 2011 (2011_03) release*).
* `mol_numheteroatoms(mol)` : returns the number of heteroatoms in a molecule (*available from Q1 2011 (2011_03) release*).
* `mol_numrings(mol)` : returns the number of rings in a molecule (*available from Q1 2011 (2011_03) release*).
* `mol_numaromaticrings(mol)` : returns the number of aromatic rings in a molecule (*available from 2013_03 release*).
* `mol_numaliphaticrings(mol)` : returns the number of aliphatic (at least one non-aromatic bond) rings in a molecule (*available from 2013_03 release*).
* `mol_numsaturatedrings(mol)` : returns the number of saturated rings in a molecule (*available from 2013_03 release*).
* `mol_numaromaticheterocycles(mol)` : returns the number of aromatic heterocycles in a molecule (*available from 2013_03 release*).
* `mol_numaliphaticheterocycles(mol)` : returns the number of aliphatic (at least one non-aromatic bond) heterocycles in a molecule (*available from 2013_03 release*).
* `mol_numsaturatedheterocycles(mol)` : returns the number of saturated heterocycles in a molecule (*available from 2013_03 release*).
* `mol_numaromaticcarbocycles(mol)` : returns the number of aromatic carbocycles in a molecule (*available from 2013_03 release*).
* `mol_numaliphaticcarbocycles(mol)` : returns the number of aliphatic (at least one non-aromatic bond) carbocycles in a molecule (*available from 2013_03 release*).
* `mol_numsaturatedcarbocycles(mol)` : returns the number of saturated carbocycles in a molecule (*available from 2013_03 release*).
* `mol_inchi(mol)` : returns an InChI for the molecule. (*available from the 2011_06 release, requires that the RDKit be built with InChI support*).
* `mol_inchikey(mol)` : returns an InChI key for the molecule. (*available from the 2011_06 release, requires that the RDKit be built with InChI support*).

Connectivity Descriptors
::::::::::::::::::::::::

* `mol_chi0v(mol)` - `mol_chi4v(mol)` :  returns the ChiXv value for a molecule for X=0-4 (*available from 2012_01 release*).
* `mol_chi0n(mol)` - `mol_chi4n(mol)` :  returns the ChiXn value for a molecule for X=0-4 (*available from 2012_01 release*).
* `mol_kappa1(mol)` - `mol_kappa3(mol)` :  returns the kappaX value for a molecule for X=1-3 (*available from 2012_01 release*).



Other
-----

* `rdkit_version()` : returns a string with the cartridge version number.

There are additional functions defined in the cartridge, but these are used for internal purposes.

License
+++++++

This document is copyright (C) 2013 by Greg Landrum

This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 License.
To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/3.0/ or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.


The intent of this license is similar to that of the RDKit itself.
In simple words: “Do whatever you want with it, but please give us some credit.”
