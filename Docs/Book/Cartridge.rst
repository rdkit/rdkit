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
 


Creating databases
******************

Configuration
-------------

The timing information below was collected on a
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


Creating a database from a file
-------------------------------

In this example I show how to load a database from the SMILES file of
commercially available compounds that is downloadable from
emolecules.com at URL
http://www.emolecules.com/doc/plus/download-database.php

If you choose to repeat this exact example yourself, please note that
it takes several hours to load the 6 million row database and generate
all fingerprints. 

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

Loading ChEMBL
--------------

Start by downloading and installing the postgresql dump from the ChEMBL website 
ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest

Connect to the database, install the cartridge, and create the schema that we'll use::

  chembl_14=# create extension if not exists rdkit;
  chembl_14=# create schema rdk;

Create the molecules and build the substructure search index::

  chembl_14=# select * into rdk.mols from (select molregno,mol_from_ctab(molfile::cstring) m  from compound_structures) tmp where m is not null;
  SELECT 1210823
  chembl_14=# create index molidx on rdk.mols using gist(m);
  CREATE INDEX
  chembl_14=# alter table rdk.mols add primary key (molregno);
  NOTICE:  ALTER TABLE / ADD PRIMARY KEY will create implicit index "mols_pkey" for table "mols"
  ALTER TABLE

Create some fingerprints and build the similarity search index::

  chembl_14=# select molregno,torsionbv_fp(m) as torsionbv,morganbv_fp(m) as mfp2,featmorganbv_fp(m) as ffp2 into rdk.fps from rdk.mols;
  SELECT 1210823
  chembl_14=# create index fps_ttbv_idx on rdk.fps using gist(torsionbv);
  CREATE INDEX
  chembl_14=# create index fps_mfp2_idx on rdk.fps using gist(mfp2);
  CREATE INDEX
  chembl_14=# create index fps_ffp2_idx on rdk.fps using gist(ffp2);
  CREATE INDEX
  chembl_14=# alter table rdk.fps add primary key (molregno);
  NOTICE:  ALTER TABLE / ADD PRIMARY KEY will create implicit index "fps_pkey" for table "fps"
  ALTER TABLE


Substructure searches
*********************

Example query molecules taken from the `eMolecules home page <http://www.emolecules.com/>`_::

  chembl_14=# select count(*) from rdk.mols where m@>'c1cccc2c1nncc2' ;
   count 
  -------
     281
  (1 row)

  Time: 184.043 ms
  chembl_14=# select count(*) from rdk.mols where m@>'c1ccnc2c1nccn2' ;
   count 
  -------
     671
  (1 row)

  Time: 449.998 ms
  chembl_14=# select count(*) from rdk.mols where m@>'c1cncc2n1ccn2' ;
   count 
  -------
     930
  (1 row)

  Time: 568.378 ms
  chembl_14=# select count(*) from rdk.mols where m@>'Nc1ncnc(N)n1' ;
   count 
  -------
    4478
  (1 row)

  Time: 721.758 ms
  chembl_14=# select count(*) from rdk.mols where m@>'c1scnn1' ;
   count 
  -------
   10908
  (1 row)

  Time: 701.036 ms
  chembl_14=# select count(*) from rdk.mols where m@>'c1cccc2c1ncs2' ;
   count 
  -------
   12823
  (1 row)

  Time: 1585.473 ms
  chembl_14=# select count(*) from rdk.mols where m@>'c1cccc2c1CNCCN2' ;
   count 
  -------
    1155
  (1 row)

  Time: 4567.222 ms

Notice that the last two queries are starting to take a while to execute and count all the results. 

Given we're searching through 1.2 million compounds these search times aren't incredibly slow, 
but it would be nice to have them quicker.

One easy way to speed things up, particularly for queries that return a large number of results, is to only 
retrieve a limited number of results::

  chembl_14=# select * from rdk.mols where m@>'c1cccc2c1CNCCN2' limit 100;
   molregno |                                                                                      m                                                                                       
  ----------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    1292129 | Cc1ccc2c(c1)C(=O)N(N(C)C)CC(=O)N2
    1013311 | CCCCC(=O)N1CC(=O)Nc2ccc(F)cc2C1c1ccccc1
    1294754 | COc1cc2c(cc1OCc1ccccc1)NC(=O)[C@@H]1CCCN1C2=O
    1012025 | O=C(c1cc2ccccc2oc1=O)N1CC(=O)Nc2ccc(Br)cc2C1c1ccc(F)cc1
     995226 | CC1Cc2ccccc2N1C(=O)CN1c2ccccc2C(=O)N(C)CC1=O
    1291875 | COC(=O)C1=NN2c3ccccc3CN([C@@H](C)c3ccccc3)C(=O)[C@@H]2[C@H]1c1ccccc1
    ...
    1116370 | COc1ccc(CC(=O)N2CC(=O)Nc3ccc(Br)cc3C2c2ccc(F)cc2)cc1OC
    1114872 | O=C1[C@@H]2[C@H](C(=O)N1Cc1ccccc1)[C@@H]1C(=O)Nc3ccccc3C(=O)N1[C@@H]2c1ccccc1
  Time: 375.747 ms


SMARTS-based queries
--------------------

Oxadiazole or thiadiazole::

    chembl_14=# select * from rdk.mols where m@>'c1[o,s]ncn1'::qmol limit 500;
     molregno |                                                                      m                                                                       
    ----------+----------------------------------------------------------------------------------------------------------------------------------------------
       534296 | Clc1ccccc1CNc1noc(-c2sccc2Br)n1
         1178 | CCCCc1oc2ccccc2c1Cc1cccc(/C(C)=C/Cn2oc(=O)[nH]c2=O)c1
       566382 | COC(=O)CCc1nc(C2CC(c3ccc(O)c(F)c3)=NO2)no1
       499261 | CS/C=C(/C)n1c(=O)onc1C(=O)c1ccc(Br)cc1
       450499 | CS(=O)(=O)c1ccc(Nc2ncnc(N3CCC(c4nc(-c5cccc(C(F)(F)F)c5)no4)CC3)c2[N+](=O)[O-])cc1
       600176 | Cc1nc(-c2c(Cl)cc(Cl)cc2-c2cnc([C@@H](C)NC(=O)N(C)O)c(F)c2)no1
         1213 | CC/C(=C\Cn1oc(=O)[nH]c1=O)c1cccc(OCc2nc(-c3ccc(C(F)(F)F)cc3)oc2C)c1
       659277 | Cn1c(N)c(CCCN)c[n+]1CC1=C(C(=O)O)N2C(=O)[C@@H](NC(=O)/C(=N\OC(C)(C)C(=O)O)c3nsc(N)n3)[C@H]2SC1
         1316 | CCCCCCCC/C(=C\Cn1oc(=O)[nH]c1=O)c1cccc(OCc2nc(-c3ccc(C(F)(F)F)cc3)oc2C)c1
       ...
         1206 | C/C(Cn1oc(=O)[nH]c1=O)=C(/C)c1cccc(OCc2nc(-c3ccc(C(F)(F)F)cc3)oc2C)c1
         1496 | Cc1oc(-c2ccccc2)nc1COc1cccc(C#CC(C)n2oc(=O)[nH]c2=O)c1
    Time: 3365.309 ms


This is slower than the pure SMILES query, this is generally true of SMARTS-based queries.

Using Stereochemistry
---------------------

Note that by default stereochemistry is not taken into account when doing substructure queries::

    chembl_14=# select * from rdk.mols where m@>'NC(=O)[C@@H]1CCCN1C=O' limit 10;
     molregno |                                                                                        m                                                                                         
    ----------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      1295889 | COc1ccc(C[C@@H](C(=O)NCC(N)=O)N(C)C(=O)[C@@H]2CCCN2C(=O)[C@H](CC(C)C)NC(=O)C(C)NC(=O)OCc2ccccc2)cc1
      1293815 | CN1C(=O)C23CC4=CC=CC(O)C4N2C(=O)C1(CO)SS3
      1293919 | CNC(=O)CNC(=O)C(NC(=O)CNC(=O)C1CCCN1C(=O)C(C)NC(=O)C(NC(=O)OC(C)(C)C)C(C)C)C(C)C
      1011887 | COC(=O)C(C)NC(=O)C1CCCN1C(=O)CNC(=O)OCc1ccccc1
      1293021 | CCC(C)C1NC(=O)C(NC(=O)C(CC(C)C)N(C)C(=O)[C@@H]2CC(O)CN2C(=O)[C@H](C)O)C(C)OC(=O)[C@H](Cc2ccc(OC)cc2)N(C)C(=O)[C@@H]2CCCN2C(=O)[C@H](CC(C)C)NC(=O)C(C)C(=O)[C@H](C(C)C)OC(=O)CC1O
      1287353 | CCC(C)C1NC(=O)C(NC(=O)C(CC(C)C)N(C)C(=O)C2CCCN2C(=O)C(C)O)C(C)OC(=O)C(Cc2ccc(OC)cc2)N(C)C(=O)C2CCCN2C(=O)C(CC(C)C)NC(=O)[C@H](C)C(=O)[C@H](C(C)C)OC(=O)CC1O
      1293647 | CCC(C)[C@@H]1NC(=O)[C@@H]2CCCN2C(=O)C(CC(O)CCl)OC(=O)CCNC(=O)[C@H](C)N(C)C(=O)[C@H](C(C)C)N(C)C1=O
      1290320 | C=CCOC(=O)[C@@H]1C[C@@H](OC(C)(C)C)CN1C(=O)[C@@H]1[C@H]2OC(C)(C)O[C@H]2CN1C(=O)OCC1c2ccccc2-c2ccccc21
      1281392 | COC1=CC2C(=O)N(C)[C@@H](C)C(=O)N3NCCC[C@@H]3C(=O)N3[C@@H](C[C@@]4(O)c5ccc(Cl)cc5N[C@@H]34)C(=O)N[C@H](C(C)C)C(=O)N3NCCC[C@@H]3C(=O)N2N=C1
      1014237 | CC(C)COC(=O)N1CC(O)CC1C(=O)Nc1ccc2c(c1)OCO2
    (10 rows)

    Time: 9.447 ms

This can be changed using the `rdkit.do_chiral_sss` configuration variable::

    chembl_14=# set rdkit.do_chiral_sss=true;
    SET
    Time: 0.241 ms
    chembl_14=# select * from rdk.mols where m@>'NC(=O)[C@@H]1CCCN1C=O' limit 10;
     molregno |                                                                                                                                                                                                                                                                                 m                                                                                                                                                                                                                                                                                 
    ----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      1295889 | COc1ccc(C[C@@H](C(=O)NCC(N)=O)N(C)C(=O)[C@@H]2CCCN2C(=O)[C@H](CC(C)C)NC(=O)C(C)NC(=O)OCc2ccccc2)cc1
      1293021 | CCC(C)C1NC(=O)C(NC(=O)C(CC(C)C)N(C)C(=O)[C@@H]2CC(O)CN2C(=O)[C@H](C)O)C(C)OC(=O)[C@H](Cc2ccc(OC)cc2)N(C)C(=O)[C@@H]2CCCN2C(=O)[C@H](CC(C)C)NC(=O)C(C)C(=O)[C@H](C(C)C)OC(=O)CC1O
      1293647 | CCC(C)[C@@H]1NC(=O)[C@@H]2CCCN2C(=O)C(CC(O)CCl)OC(=O)CCNC(=O)[C@H](C)N(C)C(=O)[C@H](C(C)C)N(C)C1=O
      1290320 | C=CCOC(=O)[C@@H]1C[C@@H](OC(C)(C)C)CN1C(=O)[C@@H]1[C@H]2OC(C)(C)O[C@H]2CN1C(=O)OCC1c2ccccc2-c2ccccc21
      1281392 | COC1=CC2C(=O)N(C)[C@@H](C)C(=O)N3NCCC[C@@H]3C(=O)N3[C@@H](C[C@@]4(O)c5ccc(Cl)cc5N[C@@H]34)C(=O)N[C@H](C(C)C)C(=O)N3NCCC[C@@H]3C(=O)N2N=C1
      1007418 | C/C=C\C=C\C(=O)N1CC2(CC(c3cccc(NC(=O)/C=C\C=C/C)c3)=NO2)C[C@H]1C(N)=O
       785530 | C/C=C/C(=O)N1CC2(CC(c3cccc(NC(=O)CC)c3)=NO2)C[C@H]1C(N)=O
      1292152 | CCCCCCCC(=O)N[C@H](C(=O)N[C@H](C(=O)N(C)[C@H](C(=O)N1CCC[C@H]1C(=O)N(C)[C@H](C)C(=O)NCc1ccc(OC)cc1OC)C(C)C)C(C)C)C(C)C
      1281390 | CC(C)[C@@H]1NC(=O)[C@@H]2C[C@@]3(O)c4ccc(Cl)cc4N[C@H]3N2C(=O)[C@H]2CCCNN2C(=O)[C@@H](C)N(C)C(=O)[C@H]2CCCNN2C(=O)[C@@H]2CCCNN2C1=O
      1057962 | CC[C@H](C)[C@@H]1NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCSC)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CCCNC(=N)N)NC(=O)CNC(=O)[C@H](Cc2ccccc2)NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CO)NC(=O)CNC(=O)[C@H](CCC(N)=O)NC(=O)[C@@H](NC(=O)[C@H](CCSC)NC(=O)[C@H](CCCCN)NC(=O)[C@@H]2CCCN2C(=O)[C@@H](N)CO)C(C)C)CSSC[C@@H](C(=O)N[C@@H](CCCCN)C(=O)N[C@H](C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](Cc2cnc[nH]2)C(=O)O)C(C)C)NC(=O)CNC(=O)[C@H](CC(C)C)NC(=O)CNC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)[C@H](CO)NC1=O
    (10 rows)

    Time: 35.383 ms


Similarity searches
*******************

Basic similarity searching::

  chembl_14=# select count(*) from rdk.fps where mfp2%morganbv_fp('Cc1ccc2nc(-c3ccc(NC(C4N(C(c5cccs5)=O)CCC4)=O)cc3)sc2c1');
   count 
  -------
      66
  (1 row)

  Time: 826.886 ms

Usually we'd like to find a sorted listed of neighbors along with the accompanying SMILES. 
This SQL function makes that pattern easy::

  chembl_14=# create or replace function get_mfp2_neighbors(smiles text)
      returns table(molregno integer, m mol, similarity double precision) as
    $$
    select molregno,m,tanimoto_sml(morganbv_fp(mol_from_smiles($1::cstring)),mfp2) as similarity
    from rdk.fps join rdk.mols using (molregno)
    where morganbv_fp(mol_from_smiles($1::cstring))%mfp2
    order by morganbv_fp(mol_from_smiles($1::cstring))<%>mfp2;
    $$ language sql stable ;
  CREATE FUNCTION
  Time: 0.856 ms
  chembl_14=# 
  chembl_14=# select * from get_mfp2_neighbors('Cc1ccc2nc(-c3ccc(NC(C4N(C(c5cccs5)=O)CCC4)=O)cc3)sc2c1') limit 10;
   molregno |                              m                              |    similarity     
  ----------+-------------------------------------------------------------+-------------------
     472512 | Cc1ccc2nc(-c3ccc(NC(=O)C4CCN(C(=O)c5cccs5)CC4)cc3)sc2c1     | 0.772727272727273
     471317 | Cc1ccc2nc(-c3ccc(NC(=O)C4CCCN(S(=O)(=O)c5cccs5)C4)cc3)sc2c1 | 0.657534246575342
     471461 | Cc1ccc2nc(-c3ccc(NC(=O)C4CCN(S(=O)(=O)c5cccs5)CC4)cc3)sc2c1 | 0.647887323943662
     471319 | Cc1ccc2nc(-c3ccc(NC(=O)C4CCN(S(=O)(=O)c5cccs5)C4)cc3)sc2c1  | 0.638888888888889
    1032469 | O=C(Nc1nc2ccc(Cl)cc2s1)[C@@H]1CCCN1C(=O)c1cccs1             | 0.623188405797101
     751668 | COc1ccc2nc(NC(=O)[C@@H]3CCCN3C(=O)c3cccs3)sc2c1             | 0.619718309859155
     471318 | Cc1ccc2nc(-c3ccc(NC(=O)C4CN(S(=O)(=O)c5cccs5)C4)cc3)sc2c1   | 0.611111111111111
     740754 | Cc1ccc(NC(=O)C2CCCN2C(=O)c2cccs2)cc1C                       | 0.606060606060606
     732905 | O=C(Nc1ccc(S(=O)(=O)N2CCCC2)cc1)C1CCCN1C(=O)c1cccs1         | 0.602941176470588
    1087495 | Cc1ccc(NC(=O)C2CCCN2C(=O)c2cccs2)c(C)c1                     | 0.597014925373134
  (10 rows)

  Time: 5453.200 ms
  chembl_14=# select * from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1') limit 10;
   molregno |                           m                           |    similarity     
  ----------+-------------------------------------------------------+-------------------
     412312 | Cc1ccc2nc(N(C)CCN(C)c3nc4ccc(C)cc4s3)sc2c1            | 0.692307692307692
     470082 | CN(CC(=O)O)c1nc2cc([N+](=O)[O-])ccc2s1                | 0.583333333333333
    1040255 | CC(=O)N(CCCN(C)C)c1nc2ccc(C)cc2s1                     | 0.571428571428571
     773946 | Cl.CC(=O)N(CCCN(C)C)c1nc2ccc(C)cc2s1                  | 0.549019607843137
    1044892 | Cc1ccc2nc(N(CCN(C)C)C(=O)c3cc(Cl)sc3Cl)sc2c1          | 0.518518518518518
    1040496 | Cc1ccc2nc(N(CCCN(C)C)C(=O)CCc3ccccc3)sc2c1            | 0.517857142857143
    1049393 | Cc1ccc2nc(N(CCCN(C)C)C(=O)CS(=O)(=O)c3ccccc3)sc2c1    | 0.517857142857143
     441378 | Cc1ccc2nc(NC(=O)CCC(=O)O)sc2c1                        | 0.510204081632653
    1042958 | Cc1ccc2nc(N(CCN(C)C)C(=O)c3ccc4ccccc4c3)sc2c1         | 0.509090909090909
    1047691 | Cc1ccc(S(=O)(=O)CC(=O)N(CCCN(C)C)c2nc3ccc(C)cc3s2)cc1 | 0.509090909090909
  (10 rows)

  Time: 1797.656 ms

Adjusting the similarity cutoff
-------------------------------

By default, the minimum similarity returned with a similarity search is 0.5. This can be adjusted with the `rdkit.tanimoto_threshold` 
(and `rdkit.dice_threshold`) configuration variables::

    chembl_14=# select count(*) from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1');
     count 
    -------
        18
    (1 row)

    Time: 1199.751 ms
    chembl_14=# set rdkit.tanimoto_threshold=0.7;
    SET
    Time: 0.191 ms
    chembl_14=# select count(*) from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1');
     count 
    -------
         0
    (1 row)

    Time: 826.058 ms
    chembl_14=# set rdkit.tanimoto_threshold=0.6;
    SET
    Time: 0.220 ms
    chembl_14=# select count(*) from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1');
     count 
    -------
         1
    (1 row)

    Time: 1092.303 ms
    chembl_14=# set rdkit.tanimoto_threshold=0.5
    chembl_14-# ;
    SET
    Time: 0.257 ms
    chembl_14=# select count(*) from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1');
     count 
    -------
        18
    (1 row)

    Time: 1081.721 ms



Using the MCS code
******************

The most straightforward use of the MCS code is to find the maximum common substructure of a group of molecules::

    chembl_20=# select fmcs(m) from rdk.mols join compound_records using (molregno) where doc_id=3;                                                                                           fmcs                                                                                           
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     [#6]1(-[#7](-[#6](-[#6]2:[#6]:[#6]:[#6](:[#6]:[#6]:2)-[#7]-[#6](-[#6]2:[#6](-[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3):[#6]:[#6]:[#6]:[#6]:2)=[#8])=[#8])-[#6]-[#6]-[#6]):[#6]:[#16]:[#6]:[#6]:1
    (1 row)

    chembl_20=# select fmcs(m) from rdk.mols join compound_records using (molregno) where doc_id=4;
                                      fmcs                                  
    ------------------------------------------------------------------------
     [#6](-[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6])-[#7]-[#6]-[#6](-,:[#6])-,:[#6]
    (1 row)

The same thing can be done with a SMILES column::

    chembl_20=# select fmcs(canonical_smiles) from compound_structures join compound_records using (molregno) where doc_id=4;
                                      fmcs                                  
    ------------------------------------------------------------------------
     [#6](-[#7]-[#6]-[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6])-[#6](-,:[#6])-,:[#6]
    (1 row)

It's also possible to adjust some of the parameters to the FMCS algorithm, though this is somewhat more painful as of this writing (the 2015_03_1 release). 
Here are a couple of examples::

    chembl_20=# select fmcs_smiles(str,'{"Threshold":0.8}') from 
    chembl_20-#   (select string_agg(m::text,' ') as str from rdk.mols 
    chembl_20(#   join compound_records using (molregno) where doc_id=4) as str ;
                                                                               fmcs_smiles                                                                            
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------
     [#6]-[#6]-[#8]-[#6](-[#6](=[#8])-[#7]-[#6](-[#6])-[#6](-,:[#6])-,:[#6])-[#6](-[#8])-[#6](-[#8])-[#6](-[#8]-[#6]-[#6])-[#6]-[#7]-[#6](-[#6])-[#6](-,:[#6])-,:[#6]
    (1 row)

    chembl_20=# select fmcs_smiles(str,'{"AtomCompare":"Any"}') from 
    chembl_20-# (select string_agg(m::text,' ') as str from rdk.mols 
    chembl_20(# join compound_records using (molregno) where doc_id=4) as str ;
                                                                                  fmcs_smiles                                                                               
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     [#6]-,:[#6,#7]-[#8,#6]-[#6,#7](-[#6,#8]-[#7,#6]-,:[#6,#7]-,:[#6,#7]-,:[#7,#6]-,:[#6])-[#6,#7]-[#6]-[#6](-[#8,#6]-[#6])-[#6,#7]-[#7,#6]-[#6]-,:[#6,#8]-,:[#7,#6]-,:[#6]
    (1 row)

*Note* The combination of ``"AtomCompare":"Any"`` and a value of ``"Threshold"`` that is less than 1.0 does a quite generic search and can results in very long search times. 
Using ``"Timeout"`` with this combination is recommended::

    chembl_20=# select fmcs_smiles(str,'{"AtomCompare":"Any","CompleteRingsOnly":true,"Threshold":0.8,"Timeout":60}') from
    chembl_20-#  (select string_agg(m::text,' ') as str from rdk.mols
    chembl_20(#   join compound_records using (molregno) where doc_id=3) as str ;
    WARNING:  findMCS timed out, result is not maximal
                                                                                              fmcs_smiles                                                                                          
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     [#8]=[#6](-[#7]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6](=[#8])-[#7]1-[#6]-[#6]-[#6]-[#6,#7]-[#6]2:[#6]-1:[#6]:[#6]:[#16]:2)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1
    (1 row)

Available parameters and their default values are:

  - MaximizeBonds (true)
  - Threshold (1.0)
  - Timeout (-1, no timeout) 
  - MatchValences (false)
  - MatchChiralTag (false) Applies to atoms
  - RingMatchesRingOnly (false)
  - CompleteRingsOnly (false)
  - MatchStereo (false)  Applies to bonds
  - AtomCompare ("Elements") can be "Elements", "Isotopes", or "Any"
  - BondCompare ("Order") can be "Order", "OrderExact", or "Any"
    
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
* `rdkit.do_chiral_sss` : toggles whether or not stereochemistry is used in substructure matching. (*available from 2013_03 release*).
* `rdkit.sss_fp_size` : the size (in bits) of the fingerprint used for substructure screening.
* `rdkit.morgan_fp_size` : the size (in bits) of morgan fingerprints
* `rdkit.featmorgan_fp_size` : the size (in bits) of featmorgan fingerprints
* `rdkit.layered_fp_size` : the size (in bits) of layered fingerprints
* `rdkit.rdkit_fp_size` : the size (in bits) of RDKit fingerprints
* `rdkit.torsion_fp_size` : the size (in bits) of topological torsion bit vector fingerprints
* `rdkit.atompair_fp_size` : the size (in bits) of atom pair bit vector fingerprints
* `rdkit.avalon_fp_size` : the size (in bits) of avalon fingerprints


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

* `morgan_fp(mol,int default 2)` : returns an `sfp` which is the count-based Morgan fingerprint for a molecule using connectivity invariants. The second argument provides the radius. This is an ECFP-like fingerprint.
* `morganbv_fp(mol,int default 2)` : returns a `bfp` which is the bit vector Morgan fingerprint for a molecule using connectivity invariants. The second argument provides the radius. This is an ECFP-like fingerprint.
* `featmorgan_fp(mol,int default 2)` : returns an `sfp` which is the count-based Morgan fingerprint for a molecule using chemical-feature invariants. The second argument provides the radius. This is an FCFP-like fingerprint.
* `featmorganbv_fp(mol,int default 2)` : returns a `bfp` which is the bit vector Morgan fingerprint for a molecule using chemical-feature invariants. The second argument provides the radius. This is an FCFP-like fingerprint.
* `rdkit_fp(mol)` : returns a `bfp` which is the RDKit fingerprint for a molecule. This is a daylight-fingerprint using hashed molecular subgraphs.
* `atompair_fp(mol)` : returns an `sfp` which is the count-based atom-pair fingerprint for a molecule.
* `atompairbv_fp(mol)` : returns a `bfp` which is the bit vector atom-pair fingerprint for a molecule.
* `torsion_fp(mol)` : returns an `sfp` which is the count-based topological-torsion fingerprint for a molecule.
* `torsionbv_fp(mol)` : returns a `bfp` which is the bit vector topological-torsion fingerprint for a molecule.
* `layered_fp(mol)` : returns a `bfp` which is the layered fingerprint for a molecule. This is an experimental substructure fingerprint using hashed molecular subgraphs.
* `maccs_fp(mol)` : returns a `bfp` which is the MACCS fingerprint for a molecule (*available from 2013_01 release*).

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
:::::::::::::::::::::::::::

* `is_valid_smiles(smiles)` : returns whether or not a SMILES string produces a valid RDKit molecule.
* `is_valid_ctab(ctab)` : returns whether or not a CTAB (mol block) string produces a valid RDKit molecule.
* `is_valid_smarts(smarts)` : returns whether or not a SMARTS string produces a valid RDKit molecule.
* `is_valid_mol_pkl(bytea)` : returns whether or not a binary string (bytea) can be converted into an RDKit molecule. (*available from Q3 2012 (2012_09) release*)

* `mol_from_smiles(smiles)` : returns a molecule for a SMILES string, NULL if the molecule construction fails.
* `mol_from_smarts(smarts)` : returns a molecule for a SMARTS string, NULL if the molecule construction fails.
* `mol_from_ctab(ctab, bool default false)` : returns a molecule for a CTAB (mol block) string, NULL if the molecule construction fails. The optional second argument controls whether or not the molecule's coordinates are saved.
* `mol_from_pkl(bytea)` : returns a molecule for a binary string (bytea), NULL if the molecule construction fails. (*available from Q3 2012 (2012_09) release*)

* `qmol_from_smiles(smiles)` : returns a query molecule for a SMILES string, NULL if the molecule construction fails. Explicit Hs in the SMILES are converted into query features on the attached atom.
* `qmol_from_ctab(ctab, bool default false)` : returns a query molecule for a CTAB (mol block) string, NULL if the molecule construction fails. Explicit Hs in the SMILES are converted into query features on the attached atom. The optional second argument controls whether or not the molecule's coordinates are saved.
* `mol_to_smiles(mol)` : returns the canonical SMILES for a molecule.
* `mol_to_smarts(mol)` : returns SMARTS string for a molecule.
* `mol_to_pkl(mol)` : returns binary string (bytea) for a molecule. (*available from Q3 2012 (2012_09) release*)
* `mol_to_ctab(mol,bool default true)` : returns a CTAB (mol block) string for a molecule. The optional second argument controls whether or not 2D coordinates will be generated for molecules that don't have coordinates. (*available from the 2014_03 release*)


Substructure operations
:::::::::::::::::::::::

* `substruct(mol,mol)` : returns whether or not the second mol is a substructure of the first.
* `substruct_count(mol,mol,bool default true)` : returns the number of substructure matches between the second molecule and the first. The third argument toggles whether or not the matches are uniquified. (*available from 2013_03 release*)


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
* `mol_formula(mol,bool default false, bool default true)` : returns a string with the molecular formula. The second argument controls whether isotope information is included in the formula; the third argument controls whether "D" and "T" are used instead of [2H] and [3H].
(*available from the 2014_03 release*)


Connectivity Descriptors
::::::::::::::::::::::::

* `mol_chi0v(mol)` - `mol_chi4v(mol)` :  returns the ChiXv value for a molecule for X=0-4 (*available from 2012_01 release*).
* `mol_chi0n(mol)` - `mol_chi4n(mol)` :  returns the ChiXn value for a molecule for X=0-4 (*available from 2012_01 release*).
* `mol_kappa1(mol)` - `mol_kappa3(mol)` :  returns the kappaX value for a molecule for X=1-3 (*available from 2012_01 release*).


MCS
:::

* `fmcs(mols)` : an aggregation function that calculates the MCS for a set of molecules
* `fmcs_smiles(text, json default '')` : calculates the MCS for a space-separated set of SMILES. The optional json argument is used to provide parameters to the MCS code.


Other
-----

* `rdkit_version()` : returns a string with the cartridge version number.

There are additional functions defined in the cartridge, but these are used for internal purposes.


Using the Cartridge from Python
+++++++++++++++++++++++++++++++

The recommended adapter for connecting to postgresql is pyscopg2
(https://pypi.python.org/pypi/psycopg2). 

Here's an example of connecting to our local copy of ChEMBL and doing
a basic substructure search:: 

  >>> import psycopg2
  >>> conn = psycopg2.connect(database='chembl_16')
  >>> curs = conn.cursor()
  >>> curs.execute('select * from rdk.mols where m@>%s',('c1cccc2c1nncc2',))
  >>> curs.fetchone()
  (9830, 'CC(C)Sc1ccc(CC2CCN(C3CCN(C(=O)c4cnnc5ccccc54)CC3)CC2)cc1')

That returns a SMILES for each molecule. If you plan to do more work
with the molecules after retrieving them, it is much more efficient to
ask postgresql to give you the molecules in pickled form::

  >>> curs.execute('select molregno,mol_send(m) from rdk.mols where m@>%s',('c1cccc2c1nncc2',))
  >>> row = curs.fetchone()
  >>> row
  (9830, <read-only buffer for 0x...>)

These pickles can then be converted into molecules:: 

  >>> from rdkit import Chem
  >>> m = Chem.Mol(str(row[1]))
  >>> Chem.MolToSmiles(m,True)
  'CC(C)Sc1ccc(CC2CCN(C3CCN(C(=O)c4cnnc5ccccc54)CC3)CC2)cc1'



License
+++++++

This document is copyright (C) 2013 by Greg Landrum

This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 License.
To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/3.0/ or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.


The intent of this license is similar to that of the RDKit itself.
In simple words: “Do whatever you want with it, but please give us some credit.”
