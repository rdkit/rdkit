# The RDKit database cartridge

## What is this?

This document is a tutorial and reference guide for the RDKit PostgreSQL cartridge.

If you find mistakes, or have suggestions for improvements, please either fix them yourselves in the source document (the .md file) or send them to the mailing list: <rdkit-discuss@lists.sourceforge.net> (you will need to subscribe first)

## Tutorial

### Introduction

### Creating databases

#### Configuration

The timing information below was collected on a commodity desktop PC (Dell Studio XPS with a 2.9GHz i7 CPU and 8GB of RAM) running Ubuntu 12.04 and using PostgreSQL v9.1.4. The database was installed with default parameters.

To improve performance while loading the database and building the index, I changed a couple of postgres configuration settings in postgresql.conf :

    synchronous_commit = off      # immediate fsync at commit
    full_page_writes = off            # recover from partial page writes

And to improve search performance, I allowed postgresql to use more memory than the extremely conservative default settings:

    shared_buffers = 2048MB           # min 128kB
                      # (change requires restart)
    work_mem = 128MB              # min 64kB

#### Creating a database from a file

In this example I show how to load a database from the SMILES file of commercially available compounds that is downloadable from emolecules.com at URL <http://downloads.emolecules.com/free/>

If you choose to repeat this exact example yourself, please note that it takes several hours to load the 6 million row database and generate all fingerprints.

First create the database and install the cartridge:

    ~/RDKit_trunk/Data/emolecules > createdb emolecules
    ~/RDKit_trunk/Data/emolecules > psql -c 'create extension rdkit' emolecules

Now create and populate a table holding the raw data:

    ~/RDKit_trunk/Data/emolecules > psql -c 'create table raw_data (id SERIAL, smiles text, emol_id integer, parent_id integer)' emolecules
    NOTICE:  CREATE TABLE will create implicit sequence "raw_data_id_seq" for serial column "raw_data.id"
    CREATE TABLE
    ~/RDKit_trunk/Data/emolecules > zcat emolecules-2013-02-01.smi.gz | sed '1d; s/\\/\\\\/g' | psql -c "copy raw_data (smiles,emol_id,parent_id) from stdin with delimiter ' '" emolecules

Create the molecule table, but only for SMILES that the RDKit accepts:

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

#### Loading ChEMBL

Start by downloading and installing the postgresql dump from the ChEMBL website <ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest>

Connect to the database, install the cartridge, and create the schema that we'll use:

    chembl_25=# create extension if not exists rdkit;
    chembl_25=# create schema rdk;

Create the molecules and build the substructure search index:

    chembl_25=# select * into rdk.mols from (select molregno,mol_from_ctab(molfile::cstring) m  from compound_structures) tmp where m is not null;
    SELECT 1870451
    chembl_25=# create index molidx on rdk.mols using gist(m);
    CREATE INDEX
    chembl_25=# alter table rdk.mols add primary key (molregno);
    ALTER TABLE

Create some fingerprints and build the similarity search index:

    chembl_25=# select molregno,torsionbv_fp(m) as torsionbv,morganbv_fp(m) as mfp2,featmorganbv_fp(m) as ffp2 into rdk.fps from rdk.mols;
    SELECT 1870451
    chembl_25=# create index fps_ttbv_idx on rdk.fps using gist(torsionbv);
    CREATE INDEX
    chembl_25=# create index fps_mfp2_idx on rdk.fps using gist(mfp2);
    CREATE INDEX
    chembl_25=# create index fps_ffp2_idx on rdk.fps using gist(ffp2);
    CREATE INDEX
    chembl_25=# alter table rdk.fps add primary key (molregno);
    ALTER TABLE

Here is a group of the commands used here (and below) in one block so that you can just paste it in at the psql prompt:

    create extension if not exists rdkit;
    create schema rdk;
    select * into rdk.mols from (select molregno,mol_from_ctab(molfile::cstring) m  from compound_structures) tmp where m is not null;
    create index molidx on rdk.mols using gist(m);
    alter table rdk.mols add primary key (molregno);
    select molregno,torsionbv_fp(m) as torsionbv,morganbv_fp(m) as mfp2,featmorganbv_fp(m) as ffp2 into rdk.fps from rdk.mols;
    create index fps_ttbv_idx on rdk.fps using gist(torsionbv);
    create index fps_mfp2_idx on rdk.fps using gist(mfp2);
    create index fps_ffp2_idx on rdk.fps using gist(ffp2);
    alter table rdk.fps add primary key (molregno);
    create or replace function get_mfp2_neighbors(smiles text)
    returns table(molregno bigint, m mol, similarity double precision) as
    $$
    select molregno,m,tanimoto_sml(morganbv_fp(mol_from_smiles($1::cstring)),mfp2) as similarity
    from rdk.fps join rdk.mols using (molregno)
    where morganbv_fp(mol_from_smiles($1::cstring))%mfp2
    order by morganbv_fp(mol_from_smiles($1::cstring))<%>mfp2;
    $$ language sql stable ;

### Substructure searches

Example query molecules taken from the [eMolecules home page](http://www.emolecules.com/):

    chembl_25=# select count(*) from rdk.mols where m@>'c1cccc2c1nncc2' ;
     count
    -------
       461
    (1 row)

    Time: 107.602 ms
    chembl_25=# select count(*) from rdk.mols where m@>'c1ccnc2c1nccn2' ;
     count
    -------
      1124
    (1 row)

    Time: 216.222 ms
    chembl_25=# select count(*) from rdk.mols where m@>'c1cncc2n1ccn2' ;
     count
    -------
      2233
    (1 row)

    Time: 88.266 ms
    chembl_25=# select count(*) from rdk.mols where m@>'Nc1ncnc(N)n1' ;
     count
    -------
      7095
    (1 row)

    Time: 327.855 ms
    chembl_25=# select count(*) from rdk.mols where m@>'c1scnn1' ;
     count
    -------
     16526
    (1 row)

    Time: 568.675 ms
    chembl_25=# select count(*) from rdk.mols where m@>'c1cccc2c1ncs2' ;
     count
    -------
     20745
    (1 row)

    Time: 998.104 ms
    chembl_25=# select count(*) from rdk.mols where m@>'c1cccc2c1CNCCN2' ;
     count
    -------
      1788
    (1 row)

    Time: 1922.273 ms

Notice that the last two queries are starting to take a while to execute and count all the results.

Given we're searching through 1.7 million compounds these search times aren't incredibly slow, but it would be nice to have them quicker.

One easy way to speed things up, particularly for queries that return a large number of results, is to only retrieve a limited number of results:

    chembl_25=# select * from rdk.mols where m@>'c1cccc2c1CNCCN2' limit 100;
     molregno |                                                      m                                                       
    ----------+--------------------------------------------------------------------------------------------------------------
      1671940 | Cc1cccc(C)c1N1C(=O)c2ccccc2NC(=O)C1C(=O)NCc1ccco1
      1318078 | COCN1C(=O)[C@@H]2C[C@@H](O)CN2C(=O)c2ccccc21
      1318783 | O/N=C1/Nc2ccccc2C(=S)N2CSCC12
      1318127 | CC(=O)O[C@H]1C[C@H]2C(=S)Nc3ccccc3C(=S)N2C1
      1308578 | O=C1Nc2cc([N+](=O)[O-])ccc2C(=O)N2CCC[C@@H]12
      1417168 | O=C(NCC(F)(F)F)C1C(=O)Nc2ccccc2C(=O)N1Cc1ccccc1
      ...
       793329 | Cc1ccc2c(c1)C(c1ccccc1)N(C(=O)c1ccc(OC(C)C)cc1)CC(=O)N2
       921215 | O=C1CN(C(=O)c2cc([N+](=O)[O-])ccc2Cl)C(c2ccc(F)cc2)c2cc(F)ccc2N1
       790949 | CCOC(=O)[C@H]1[C@H]2COc3ccc(Cl)cc3[C@@H]2N2C(=O)c3cc(C)ccc3NC(=O)[C@@]12C
       760998 | CC(=O)N1CC(=O)Nc2ccc(Cl)cc2C1c1ccc(F)cc1
    (100 rows)
    
    Time: 97.357 ms

#### SMARTS-based queries

Oxadiazole or thiadiazole:

    chembl_25=# select * from rdk.mols where m@>'c1[o,s]ncn1'::qmol limit 500;
     molregno |                                                 m
    ----------+---------------------------------------------------------------------------------------------------
      1882516 | COc1cccc(CN(C)Cc2nc(C(C)C)no2)c1
      2194441 | Cc1nc([C@](C)(O)C#Cc2ccc3c(c2)-c2nc(C(N)=O)sc2[C@@H](F)CO3)no1
      1881742 | CCOc1ccc(C(F)(F)F)cc1NC(=O)NCc1noc(C)n1
      1949861 | FC(F)(F)c1ccc(-c2nc(-c3ccc4nc[nH]c4c3)no2)cc1
      1949860 | FC(F)(F)c1cccc(-c2nc(-c3ccc4nc[nH]c4c3)no2)c1
      2172627 | O=c1[nH]cc(-c2cc(Cl)ccc2Oc2cc(F)c(S(=O)(=O)Nc3ncns3)cc2F)n2cncc12
      ...
      1848026 | O=C1CCCN1c1cccc(-c2noc([C@H]3CCCCN3C(=O)COc3ccccc3)n2)c1
      1848027 | O=C1CN(c2cccc(-c3noc([C@H]4CCCCN4C(=O)COc4ccccc4)n3)c2)C(=O)N1
      1848036 | CN(C)C(=O)CCC(=O)Nc1cc(F)cc(-c2noc([C@H]3CCCCN3C(=O)COc3ccccc3)n2)c1
      1852688 | CC(Sc1nc(N)cc(N)n1)c1nc(C(C)(C)C)no1
    (500 rows)

    Time: 761.847 ms

This is slower than the pure SMILES query, this is generally true of SMARTS-based queries.

#### Using Stereochemistry

Note that by default stereochemistry is not taken into account when doing substructure queries:

    chembl_25=# select * from rdk.mols where m@>'NC(=O)[C@@H]1CCCN1C=O' limit 10;
     molregno |                                                 m
    ----------+---------------------------------------------------------------------------------------------------
      2213985 | CC[C@H](C)[C@@H]1NC(=O)[C@@H]2CCCN2C(=O)[C@@H]2CCCN2C(=O)[C@H]([C@@H](C)CC)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H]([C@H](C)O)NC(=O)[C@@H]2CSSC[C@H](NC1=O)C(=O)N[C@@H](Cc1cnc[nH]1)C(=O)N[C@H](Cc1ccccc1)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](Cc1c[nH]c3ccccc13)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N2
      1956682 | NC(=O)[C@@H]1CCCN1C(=O)[C@H](Cc1nc(I)[nH]c1I)NC(=O)c1cnccn1
      2212188 | CN1C(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H](Cc2ccc(O)cc2)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc2ccc3ccccc3c2)NC(=O)[C@@H]1CC(=O)O
      2053463 | NCCCC[C@H](NC(=O)[C@H](Cc1ccc(OP(=O)(O)O)cc1)NC(=O)Cc1ccccc1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](Cc1ccccc1)C(N)=O
      2060743 | CCCCCCCCCCCCCCCCNC(=O)CN(CC(=O)NC(C)(C)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CC(N)=O)C(N)=O)C(=O)c1cccnc1
      2060744 | CCCCCCCCCCCCCCCCN(CCCCCCCCCCCCCCCC)CCCCCC(=O)NC(C)(C)C(=O)NC(Cc1ccccc1)C(=O)NC(CC(C)C)C(=O)NC(Cc1ccccc1)C(=O)NC(CCCNC(=N)N)C(=O)N1CCCC1C(=O)NC(CCCNC(=N)N)C(=O)NC(CC(N)=O)C(N)=O
      2077784 | CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCSC)NC1=O
      2077779 | CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC[S+](C)[O-])NC1=O
      2077782 | CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@@H]2CCCN2C(=O)[C@H](CCSC)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC[S+](C)[O-])NC1=O
      2077780 | CC(C)C[C@@H]1NC(=O)[C@H](CC[S+](C)[O-])NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@@H]2CCCN2C(=O)[C@H](CC[S+](C)[O-])NC1=O
    (10 rows)
    
This can be changed using the rdkit.do\_chiral\_sss configuration variable:

    chembl_25=# set rdkit.do_chiral_sss=true;
    SET
    Time: 0.241 ms
    chembl_25=# select * from rdk.mols where m@>'NC(=O)[C@@H]1CCCN1C=O' limit 10;
     molregno |                                                 m
    ----------+---------------------------------------------------------------------------------------------------
      2213985 | CC[C@H](C)[C@@H]1NC(=O)[C@@H]2CCCN2C(=O)[C@@H]2CCCN2C(=O)[C@H]([C@@H](C)CC)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H]([C@H](C)O)NC(=O)[C@@H]2CSSC[C@H](NC1=O)C(=O)N[C@@H](Cc1cnc[nH]1)C(=O)N[C@H](Cc1ccccc1)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](Cc1c[nH]c3ccccc13)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N2
      1956682 | NC(=O)[C@@H]1CCCN1C(=O)[C@H](Cc1nc(I)[nH]c1I)NC(=O)c1cnccn1
      2212188 | CN1C(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H](Cc2ccc(O)cc2)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc2ccc3ccccc3c2)NC(=O)[C@@H]1CC(=O)O
      2053463 | NCCCC[C@H](NC(=O)[C@H](Cc1ccc(OP(=O)(O)O)cc1)NC(=O)Cc1ccccc1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](Cc1ccccc1)C(N)=O
      2060743 | CCCCCCCCCCCCCCCCNC(=O)CN(CC(=O)NC(C)(C)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CC(N)=O)C(N)=O)C(=O)c1cccnc1
      2077784 | CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCSC)NC1=O
      2077779 | CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC[S+](C)[O-])NC1=O
      2077782 | CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@@H]2CCCN2C(=O)[C@H](CCSC)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC[S+](C)[O-])NC1=O
      2077780 | CC(C)C[C@@H]1NC(=O)[C@H](CC[S+](C)[O-])NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@@H]2CCCN2C(=O)[C@H](CC[S+](C)[O-])NC1=O
      2211488 | CC[C@H](C)[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](C(=O)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(=O)N[C@H](CCC(=O)N[C@@H](CCC(=O)N[C@@H](CC(C)C)C(=O)O)Cc1ccccc1)Cc1ccccc1)C(C)C)[C@@H](C)CC
    (10 rows)

    Time: 6.181 ms

#### Tuning queries

It is frequently useful to be able to exert a bit more control over substructure queries without
having to construct complex SMARTS queries. The cartridge function `mol_adjust_query_properties()`
can be used to do just this. Here is an example of the default behavior, using a  
query for 2,6 di-substituted pyridines:

    chembl_25=# select molregno,m from rdk.mols where m@>mol_adjust_query_properties('*c1cccc(NC(=O)*)n1') limit 10;
     molregno |                                                 m
    ----------+---------------------------------------------------------------------------------------------------
      1609520 | Cc1cccc(NC(=O)c2cc(Br)ccc2C(=O)O)n1
      1141456 | CCN(CC)CCCn1cc(NC(=O)Nc2cccc(-c3ccccc3)n2)c2ccccc21
      1431198 | Cc1cccc(NC(=O)c2nc(C)sc2Nc2cccnc2)n1
       734975 | Cc1cccc(NC(=O)CN(C)S(=O)(=O)c2ccc(Cl)cc2)n1
       760426 | Cc1cccc(NC(=O)CCCn2cc([N+](=O)[O-])cn2)n1
       782786 | Cc1cccc(NC(=O)CN2C(=O)NC(C)(c3ccc4ccccc4c3)C2=O)n1
      1478990 | Cc1cccc(NC(=O)Cn2c(=O)sc3cc(C(=O)c4ccccc4)ccc32)n1
      1478787 | Cc1cccc(NC(=O)Cn2c(=O)sc3cc(C(=O)c4ccccc4F)ccc32)n1
      1955608 | C[C@H](N)C(=O)Nc1cccc(N)n1
       773911 | Cc1cccc(NC(=O)c2c(-c3ccccc3)noc2C)n1
    (10 rows)

    Time: 11.895 ms

By default `mol_adjust_query_properties()` makes the following changes to the molecule:

- Converts dummy atoms into "any" queries
- Adds a degree query to every ring atom so that its substitution must match what was provided
- Aromaticity perception is done (if it hasn't been done already)

We can control the behavior by providing an additional JSON argument. Here's an example
where we disable the additional degree queries:

    chembl_25=# select molregno,m from rdk.mols where m@>mol_adjust_query_properties('*c1cccc(NC(=O)*)n1',
    chembl_25(# '{"adjustDegree":false}') limit 10;
     molregno |                                                 m
    ----------+---------------------------------------------------------------------------------------------------
      2146308 | CCn1ncc2cc3nc(c21)NCCOC[C@H](c1ccccc1)NC(=O)N3
      2137309 | CCn1ncc2cc3nc(c21)CCCO[C@@H](O)[C@H](c1ccccc1)NC(=O)N3
      2102593 | CCn1ncc2cc3nc(c21)CCCO[C@@H]([C@@H](C)O)[C@@H](c1ccccc1)NC(=O)N3
      2171613 | CCn1ncc2cc3nc(c21)CCCO[C@@H]([C@H](C)O)[C@@H](c1ccccc1)NC(=O)N3
      2111904 | CCn1ncc2cc3nc(c21)C[C@H](O)COC[C@H](c1cccc(Cl)c1)NC(=O)N3
      2173410 | CCn1ncc2cc3nc(c21)CCCOC[C@H](c1ccccc1)NC(=O)N3
      2189450 | Cn1ncc2cc3nc(c21)CCCOC[C@H](c1ccccc1)NC(=O)N3
      2195752 | CCn1ncc2cc3nc(c21)C[C@H](O)COC[C@H](c1ccccc1)NC(=O)N3
      1609520 | Cc1cccc(NC(=O)c2cc(Br)ccc2C(=O)O)n1
      1141456 | CCN(CC)CCCn1cc(NC(=O)Nc2cccc(-c3ccccc3)n2)c2ccccc21
    (10 rows)

    Time: 10.780 ms

or where we don't add the additional degree queries to ring atoms or dummies (they are only
added to chain atoms):

    chembl_25=# select molregno,m from rdk.mols where m@>mol_adjust_query_properties('*c1cccc(NC(=O)*)n1',
    chembl_25(# '{"adjustDegree":true,"adjustDegreeFlags":"IGNORERINGS|IGNOREDUMMIES"}') limit 10;
     molregno |                                                 m
    ----------+---------------------------------------------------------------------------------------------------
      2146308 | CCn1ncc2cc3nc(c21)NCCOC[C@H](c1ccccc1)NC(=O)N3
      2137309 | CCn1ncc2cc3nc(c21)CCCO[C@@H](O)[C@H](c1ccccc1)NC(=O)N3
      2102593 | CCn1ncc2cc3nc(c21)CCCO[C@@H]([C@@H](C)O)[C@@H](c1ccccc1)NC(=O)N3
      2171613 | CCn1ncc2cc3nc(c21)CCCO[C@@H]([C@H](C)O)[C@@H](c1ccccc1)NC(=O)N3
      2111904 | CCn1ncc2cc3nc(c21)C[C@H](O)COC[C@H](c1cccc(Cl)c1)NC(=O)N3
      2173410 | CCn1ncc2cc3nc(c21)CCCOC[C@H](c1ccccc1)NC(=O)N3
      2189450 | Cn1ncc2cc3nc(c21)CCCOC[C@H](c1ccccc1)NC(=O)N3
      2195752 | CCn1ncc2cc3nc(c21)C[C@H](O)COC[C@H](c1ccccc1)NC(=O)N3
      1609520 | Cc1cccc(NC(=O)c2cc(Br)ccc2C(=O)O)n1
      1141456 | CCN(CC)CCCn1cc(NC(=O)Nc2cccc(-c3ccccc3)n2)c2ccccc21
    (10 rows)

    Time: 12.827 ms

The options available are:

- **adjustDegree** (default: true) : adds a query to match the input atomic degree
- **adjustDegreeFlags** (default: ADJUST_IGNOREDUMMIES | ADJUST_IGNORECHAINS) controls where the degree is adjusted
- **adjustRingCount** (default: false) : adds a query to match the input ring count
- **adjustRingCountFlags** (default: ADJUST_IGNOREDUMMIES | ADJUST_IGNORECHAINS) controls where the ring count is adjusted
- **makeDummiesQueries** (default: true) : convert dummy atoms in the input structure into any-atom queries
- **aromatizeIfPossible** (default: true) : run the aromaticity perception algorithm on the input structure (note: this is largely redundant since molecules built from smiles always have aromaticity perceived)
- **makeBondsGeneric** (default: false) : convert bonds into any-bond queries
- **makeBondsGenericFlags** (default: false) : controls which bonds are made generic
- **makeAtomsGeneric** (default: false) : convert atoms into any-atom queries
- **makeAtomsGenericFlags** (default: false) : controls which atoms are made generic

The various `Flags` arguments mentioned above, which control where particular options are applied,
are constructed by combining operations from the list below with the `|` character.

- **IGNORENONE** : apply the operation to all atoms
- **IGNORERINGS** : do not apply the operation to ring atoms
- **IGNORECHAINS** : do not apply the operation to chain atoms
- **IGNOREDUMMIES** : do not apply the operation to dummy atoms
- **IGNORENONDUMMIES** : do not apply the operation to non-dummy atoms
- **IGNOREALL** : do not apply the operation to any atoms

### Similarity searches

Basic similarity searching:

    chembl_25=# select count(*) from rdk.fps where mfp2%morganbv_fp('Cc1ccc2nc(-c3ccc(NC(C4N(C(c5cccs5)=O)CCC4)=O)cc3)sc2c1');
     count
    -------
        67
    (1 row)

    Time: 177.579 ms

Usually we'd like to find a sorted listed of neighbors along with the accompanying SMILES. This SQL function makes that pattern easy:

    chembl_25=# create or replace function get_mfp2_neighbors(smiles text)
        returns table(molregno bigint, m mol, similarity double precision) as
      $$
      select molregno,m,tanimoto_sml(morganbv_fp(mol_from_smiles($1::cstring)),mfp2) as similarity
      from rdk.fps join rdk.mols using (molregno)
      where morganbv_fp(mol_from_smiles($1::cstring))%mfp2
      order by morganbv_fp(mol_from_smiles($1::cstring))<%>mfp2;
      $$ language sql stable ;
    CREATE FUNCTION
    Time: 0.856 ms
    chembl_25=# select * from get_mfp2_neighbors('Cc1ccc2nc(-c3ccc(NC(C4N(C(c5cccs5)=O)CCC4)=O)cc3)sc2c1') limit 10;
     molregno |                                m                                 |    similarity
    ----------+------------------------------------------------------------------+-------------------
       751668 | COc1ccc2nc(NC(=O)[C@@H]3CCCN3C(=O)c3cccs3)sc2c1                  | 0.619718309859155
       740754 | Cc1ccc(NC(=O)C2CCCN2C(=O)c2cccs2)cc1C                            | 0.606060606060606
       732905 | O=C(Nc1ccc(S(=O)(=O)N2CCCC2)cc1)C1CCCN1C(=O)c1cccs1              | 0.602941176470588
       810850 | Cc1cc(C)n(-c2ccc(NC(=O)C3CCCCN3C(=O)c3cccs3)cc2)n1               | 0.583333333333333
      1224407 | O=C(Nc1cccc(S(=O)(=O)N2CCCC2)c1)C1CCCN1C(=O)c1cccs1              | 0.579710144927536
       779258 | CC1CCN(S(=O)(=O)c2ccc(NC(=O)[C@@H]3CCCN3C(=O)c3cccs3)cc2)CC1     | 0.569444444444444
       472441 | Cc1ccc2nc(-c3ccc(NC(=O)C4CCN(S(=O)(=O)C(C)C)CC4)cc3)sc2c1        | 0.569444444444444
       745651 | Cc1ccc(NC(=O)[C@@H]2CCCN2C(=O)c2cccs2)cc1S(=O)(=O)N1CCCCC1       | 0.567567567567568
       472510 | Cc1ccc2nc(-c3ccc(NC(=O)C4CCN(S(=O)(=O)c5cccc(Cl)c5)CC4)cc3)sc2c1 | 0.565789473684211
      1233426 | Cc1cccc2sc(NC(=O)[C@@H]3CCCN3C(=O)c3cccs3)nc12                   | 0.563380281690141
    (10 rows)

    Time: 28.909 ms
    chembl_25=# select * from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1') limit 10;
     molregno |                                m                         |    similarity
    ----------+----------------------------------------------------------+-------------------
      2138088 | CN(CC(=O)O)c1nc2ccc([N+](=O)[O-])cc2s1                   | 0.673913043478261
      1040255 | CC(=O)N(CCCN(C)C)c1nc2ccc(C)cc2s1                        | 0.571428571428571
       773946 | CC(=O)N(CCCN(C)C)c1nc2ccc(C)cc2s1.Cl                     |              0.56
      1044892 | Cc1ccc2nc(N(CCN(C)C)C(=O)c3cc(Cl)sc3Cl)sc2c1             | 0.518518518518518
       441378 | Cc1ccc2nc(NC(=O)CCC(=O)O)sc2c1                           | 0.510204081632653
      1047691 | Cc1ccc(S(=O)(=O)CC(=O)N(CCCN(C)C)c2nc3ccc(C)cc3s2)cc1    | 0.509090909090909
      1042958 | Cc1ccc2nc(N(CCN(C)C)C(=O)c3ccc4ccccc4c3)sc2c1            | 0.509090909090909
      1015485 | Cc1ccc2nc(N(Cc3cccnc3)C(=O)Cc3ccccc3)sc2c1               |               0.5
       994843 | Cc1ccc(S(=O)(=O)CC(=O)N(CCCN(C)C)c2nc3ccc(C)cc3s2)cc1.Cl |               0.5
       841938 | Cc1ccc2nc(N(CCN(C)C)C(=O)c3ccc4ccccc4c3)sc2c1.Cl         |               0.5
    (10 rows)

    Time: 41.623 ms

#### Adjusting the similarity cutoff

By default, the minimum similarity returned with a similarity search is 0.5. This can be adjusted with the rdkit.tanimoto\_threshold (and rdkit.dice\_threshold) configuration variables:

    chembl_25=# select count(*) from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1');
     count
    -------
        21
    (1 row)

    Time: 181.438 ms
    chembl_25=# set rdkit.tanimoto_threshold=0.7;
    SET
    Time: 0.047 ms
    chembl_25=# select count(*) from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1');
     count
    -------
         0
    (1 row)

    Time: 161.228 ms
    chembl_25=# set rdkit.tanimoto_threshold=0.6;
    SET
    Time: 0.045 ms
    chembl_25=# select count(*) from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1');
     count
    -------
         2
    (1 row)

    Time: 184.275 ms
    chembl_25=# set rdkit.tanimoto_threshold=0.5;
    SET
    Time: 0.055 ms
    chembl_25=# select count(*) from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1');
     count
    -------
        21
    (1 row)

    Time: 181.100 ms

### Using the MCS code

The most straightforward use of the MCS code is to find the maximum common substructure of a group of molecules:

    chembl_25=# select fmcs(m::text) from rdk.mols join compound_records using (molregno) where doc_id=4;
                                      fmcs                                  
    ------------------------------------------------------------------------
     [#6](-[#6]-[#7]-[#6]-[#6](-,:[#6])-,:[#6])-,:[#6]-,:[#6]-,:[#6]-,:[#6]
    (1 row)

    Time: 31.041 ms
    chembl_25=# select fmcs(m::text) from rdk.mols join compound_records using (molregno) where doc_id=5;
                                                                       fmcs                                                                   
    ------------------------------------------------------------------------------------------------------------------------------------------
     [#6]-[#6](=[#8])-[#7]-[#6](-[#6](=[#8])-[#7]1-[#6]-[#6]-[#6]-[#6]-1-[#6](=[#8])-[#7]-[#6](-[#6](=[#8])-[#8])-[#6]-[#6])-[#6](-[#6])-[#6]
    (1 row)

    Time: 705.535 ms

The same thing can be done with a SMILES column:

    chembl_25=# select fmcs(canonical_smiles) from compound_structures join compound_records using (molregno) where doc_id=4;
                                      fmcs                                  
    ------------------------------------------------------------------------
     [#6](-[#7]-[#6]-[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6])-[#6](-,:[#6])-,:[#6]
    (1 row)

    Time: 128.879 ms

It's also possible to adjust some of the parameters to the FMCS algorithm, though this is somewhat more painful as of this writing (the 2017\_03 release cycle). Here are a couple of examples:

    chembl_25=# select fmcs_smiles(str,'{"Threshold":0.8}') from
    chembl_25-#    (select string_agg(m::text,' ') as str from rdk.mols
    chembl_25(#    join compound_records using (molregno) where doc_id=4) as str ;

                                                                               fmcs_smiles                                                                            
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------
     [#6]-[#6]-[#8]-[#6](-[#6](=[#8])-[#7]-[#6](-[#6])-[#6](-,:[#6])-,:[#6])-[#6](-[#8])-[#6](-[#8])-[#6](-[#8]-[#6]-[#6])-[#6]-[#7]-[#6](-[#6])-[#6](-,:[#6])-,:[#6]
    (1 row)

    Time: 9673.949 ms
    chembl_25=#
    chembl_25=# select fmcs_smiles(str,'{"AtomCompare":"Any"}') from
    chembl_25-#    (select string_agg(m::text,' ') as str from rdk.mols
    chembl_25(#    join compound_records using (molregno) where doc_id=4) as str ;
                                                                                  fmcs_smiles                                                                               
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     [#6]-,:[#6,#7]-[#8,#6]-[#6,#7](-[#6,#8]-[#7,#6]-,:[#6,#7]-,:[#6,#7]-,:[#7,#6]-,:[#6])-[#6,#7]-[#6]-[#6](-[#8,#6]-[#6])-[#6,#7]-[#7,#6]-[#6]-,:[#6,#8]-,:[#7,#6]-,:[#6]
    (1 row)

    Time: 304.332 ms

*Note* The combination of `"AtomCompare":"Any"` and a value of `"Threshold"` that is less than 1.0 does a quite generic search and can results in very long search times. Using `"Timeout"` with this combination is recommended:

    chembl_25=# select fmcs_smiles(str,'{"AtomCompare":"Any","CompleteRingsOnly":true,"Threshold":0.8,"Timeout":60}') from
    chembl_25-#    (select string_agg(m::text,' ') as str from rdk.mols
    chembl_25(#    join compound_records using (molregno) where doc_id=3) as str ;

    WARNING:  findMCS timed out, result is not maximal
                                                                                              fmcs_smiles                                                                    

    -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ----------------------
     [#8]=[#6](-[#7]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6](=[#8])-[#7]1-[#6]-[#6]-[#6]-[#6,#7]-[#6]2:[#6]-1:[#6]:[#6]:[#16]:2)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]1:[#6]:
    [#6]:[#6]:[#6]:[#6]:1
    (1 row)

    Time: 60479.753 ms

Available parameters and their default values are:

> -   MaximizeBonds (true)
> -   Threshold (1.0)
> -   Timeout (-1, no timeout)
> -   MatchValences (false)
> -   MatchChiralTag (false) Applies to atoms
> -   RingMatchesRingOnly (false)
> -   CompleteRingsOnly (false)
> -   MatchStereo (false) Applies to bonds
> -   AtomCompare ("Elements") can be "Elements", "Isotopes", or "Any"
> -   BondCompare ("Order") can be "Order", "OrderExact", or "Any"

## Reference Guide

### New Types

-   mol : an rdkit molecule. Can be created from a SMILES via direct type conversion, for example: 'c1ccccc1'::mol creates a molecule from the SMILES 'c1ccccc1'
-   qmol : an rdkit molecule containing query features (i.e. constructed from SMARTS). Can be created from a SMARTS via direct type conversion, for example: 'c1cccc[c,n]1'::qmol creates a query molecule from the SMARTS 'c1cccc[c,n]1'
-   sfp : a sparse count vector fingerprint (SparseIntVect in C++ and Python)
-   bfp : a bit vector fingerprint (ExplicitBitVect in C++ and Python)

### Parameters

-   rdkit.tanimoto\_threshold : threshold value for the Tanimoto similarity operator. Searches done using Tanimoto similarity will only return results with a similarity of at least this value.
-   rdkit.dice\_threshold : threshold value for the Dice similiarty operator. Searches done using Dice similarity will only return results with a similarity of at least this value.
-   rdkit.do\_chiral\_sss : toggles whether or not stereochemistry is used in substructure matching. (*available from 2013\_03 release*).
-   rdkit.sss\_fp\_size : the size (in bits) of the fingerprint used for substructure screening.
-   rdkit.morgan\_fp\_size : the size (in bits) of morgan fingerprints
-   rdkit.featmorgan\_fp\_size : the size (in bits) of featmorgan fingerprints
-   rdkit.layered\_fp\_size : the size (in bits) of layered fingerprints
-   rdkit.rdkit\_fp\_size : the size (in bits) of RDKit fingerprints
-   rdkit.torsion\_fp\_size : the size (in bits) of topological torsion bit vector fingerprints
-   rdkit.atompair\_fp\_size : the size (in bits) of atom pair bit vector fingerprints
-   rdkit.avalon\_fp\_size : the size (in bits) of avalon fingerprints

### Operators

#### Similarity search

-   % : operator used for similarity searches using Tanimoto similarity. Returns whether or not the Tanimoto similarity between two fingerprints (either two sfp or two bfp values) exceeds rdkit.tanimoto\_threshold.
-   \# : operator used for similarity searches using Dice similarity. Returns whether or not the Dice similarity between two fingerprints (either two sfp or two bfp values) exceeds rdkit.dice\_threshold.
-   <%\> : used for Tanimoto KNN searches (to return ordered lists of neighbors).
-   <\#\> : used for Dice KNN searches (to return ordered lists of neighbors).

#### Substructure and exact structure search

-   @\> : substructure search operator. Returns whether or not the mol or qmol on the right is a substructure of the mol on the left.
-   <@ : substructure search operator. Returns whether or not the mol or qmol on the left is a substructure of the mol on the right.
-   @= : returns whether or not two molecules are the same.

#### Molecule comparison

-   < : returns whether or not the left mol is less than the right mol
-   \> : returns whether or not the left mol is greater than the right mol
-   = : returns whether or not the left mol is equal to the right mol
-   <= : returns whether or not the left mol is less than or equal to the right mol
-   \>= : returns whether or not the left mol is greater than or equal to the right mol

*Note* Two molecules are compared by making the following comparisons in order. Later comparisons are only made if the preceding values are equal:

\# Number of atoms \# Number of bonds \# Molecular weight \# Number of rings

If all of the above are the same and the second molecule is a substructure of the first, the molecules are declared equal, Otherwise (should not happen) the first molecule is arbitrarily defined to be less than the second.

There are additional operators defined in the cartridge, but these are used for internal purposes.

### Functions

#### Fingerprint Related

##### Generating fingerprints

-   morgan\_fp(mol,int default 2) : returns an sfp which is the count-based Morgan fingerprint for a molecule using connectivity invariants. The second argument provides the radius. This is an ECFP-like fingerprint.
-   morganbv\_fp(mol,int default 2) : returns a bfp which is the bit vector Morgan fingerprint for a molecule using connectivity invariants. The second argument provides the radius. This is an ECFP-like fingerprint.
-   featmorgan\_fp(mol,int default 2) : returns an sfp which is the count-based Morgan fingerprint for a molecule using chemical-feature invariants. The second argument provides the radius. This is an FCFP-like fingerprint.
-   featmorganbv\_fp(mol,int default 2) : returns a bfp which is the bit vector Morgan fingerprint for a molecule using chemical-feature invariants. The second argument provides the radius. This is an FCFP-like fingerprint.
-   rdkit\_fp(mol) : returns a bfp which is the RDKit fingerprint for a molecule. This is a daylight-fingerprint using hashed molecular subgraphs.
-   atompair\_fp(mol) : returns an sfp which is the count-based atom-pair fingerprint for a molecule.
-   atompairbv\_fp(mol) : returns a bfp which is the bit vector atom-pair fingerprint for a molecule.
-   torsion\_fp(mol) : returns an sfp which is the count-based topological-torsion fingerprint for a molecule.
-   torsionbv\_fp(mol) : returns a bfp which is the bit vector topological-torsion fingerprint for a molecule.
-   layered\_fp(mol) : returns a bfp which is the layered fingerprint for a molecule. This is an experimental substructure fingerprint using hashed molecular subgraphs.
-   maccs\_fp(mol) : returns a bfp which is the MACCS fingerprint for a molecule (*available from 2013\_01 release*).

##### Working with fingerprints

-   tanimoto\_sml(fp,fp) : returns the Tanimoto similarity between two fingerprints of the same type (either two sfp or two bfp values).
-   dice\_sml(fp,fp) : returns the Dice similarity between two fingerprints of the same type (either two sfp or two bfp values).
-   size(bfp) : returns the length of (number of bits in) a bfp.
-   add(sfp,sfp) : returns an sfp formed by the element-wise addition of the two sfp arguments.
-   subtract(sfp,sfp) : returns an sfp formed by the element-wise subtraction of the two sfp arguments.
-   all\_values\_lt(sfp,int) : returns a boolean indicating whether or not all elements of the sfp argument are less than the int argument.
-   all\_values\_gt(sfp,int) : returns a boolean indicating whether or not all elements of the sfp argument are greater than the int argument.

##### Fingerprint I/O

-   bfp\_to\_binary\_text(bfp) : returns a bytea with the binary string representation of the fingerprint that can be converted back into an RDKit fingerprint in other software. (*available from Q3 2012 (2012\_09) release*)
-   bfp\_from\_binary\_text(bytea) : constructs a bfp from a binary string representation of the fingerprint. (*available from Q3 2012 (2012\_09) release*)

#### Molecule Related

##### Molecule I/O and Validation

-   is\_valid\_smiles(smiles) : returns whether or not a SMILES string produces a valid RDKit molecule.
-   is\_valid\_ctab(ctab) : returns whether or not a CTAB (mol block) string produces a valid RDKit molecule.
-   is\_valid\_smarts(smarts) : returns whether or not a SMARTS string produces a valid RDKit molecule.
-   is\_valid\_mol\_pkl(bytea) : returns whether or not a binary string (bytea) can be converted into an RDKit molecule. (*available from Q3 2012 (2012\_09) release*)
-   mol\_from\_smiles(smiles) : returns a molecule for a SMILES string, NULL if the molecule construction fails.
-   mol\_from\_smarts(smarts) : returns a molecule for a SMARTS string, NULL if the molecule construction fails.
-   mol\_from\_ctab(ctab, bool default false) : returns a molecule for a CTAB (mol block) string, NULL if the molecule construction fails. The optional second argument controls whether or not the molecule's coordinates are saved.
-   mol\_from\_pkl(bytea) : returns a molecule for a binary string (bytea), NULL if the molecule construction fails. (*available from Q3 2012 (2012\_09) release*)
-   qmol\_from\_smiles(smiles) : returns a query molecule for a SMILES string, NULL if the molecule construction fails. Explicit Hs in the SMILES are converted into query features on the attached atom.
-   qmol\_from\_ctab(ctab, bool default false) : returns a query molecule for a CTAB (mol block) string, NULL if the molecule construction fails. Explicit Hs in the SMILES are converted into query features on the attached atom. The optional second argument controls whether or not the molecule's coordinates are saved.
-   mol\_to\_smiles(mol) : returns the canonical SMILES for a molecule.
-   mol\_to\_smarts(mol) : returns SMARTS string for a molecule.
-   mol\_to\_pkl(mol) : returns binary string (bytea) for a molecule. (*available from Q3 2012 (2012\_09) release*)
-   mol\_to\_ctab(mol,bool default true) : returns a CTAB (mol block) string for a molecule. The optional second argument controls whether or not 2D coordinates will be generated for molecules that don't have coordinates. (*available from the 2014\_03 release*)
-   mol\_to\_svg(mol,string default '',int default 250, int default 200, string default '') : returns an SVG with a drawing of the molecule. The optional parameters are a string to use as the legend, the width of the image, the height of the image, and a JSON with additional rendering parameters. (*available from the 2016\_09 release*)


##### Substructure operations

-   substruct(mol,mol) : returns whether or not the second mol is a substructure of the first.
-   substruct\_count(mol,mol,bool default true) : returns the number of substructure matches between the second molecule and the first. The third argument toggles whether or not the matches are uniquified. (*available from 2013\_03 release*)
-   mol_adjust_query_properties(mol,string default '') : returns a new molecule with additional query information attached. (*available from the 2016\_09 release*)

##### Descriptors

-   mol\_amw(mol) : returns the AMW for a molecule.
-   mol\_logp(mol) : returns the MolLogP for a molecule.
-   mol\_tpsa(mol) : returns the topological polar surface area for a molecule (*available from Q1 2011 (2011\_03) release*).
-   mol\_fractioncsp3(mol) : returns the fraction of carbons that are sp3 hybridized (*available from 2013\_03 release*).
-   mol\_hba(mol) : returns the number of Lipinski H-bond acceptors (i.e. number of Os and Ns) for a molecule.
-   mol\_hbd(mol) : returns the number of Lipinski H-bond donors (i.e. number of Os and Ns that have at least one H) for a molecule.
-   mol\_numatoms(mol) : returns the total number of atoms in a molecule.
-   mol\_numheavyatoms(mol) : returns the number of heavy atoms in a molecule.
-   mol\_numrotatablebonds(mol) : returns the number of rotatable bonds in a molecule (*available from Q1 2011 (2011\_03) release*).
-   mol\_numheteroatoms(mol) : returns the number of heteroatoms in a molecule (*available from Q1 2011 (2011\_03) release*).
-   mol\_numrings(mol) : returns the number of rings in a molecule (*available from Q1 2011 (2011\_03) release*).
-   mol\_numaromaticrings(mol) : returns the number of aromatic rings in a molecule (*available from 2013\_03 release*).
-   mol\_numaliphaticrings(mol) : returns the number of aliphatic (at least one non-aromatic bond) rings in a molecule (*available from 2013\_03 release*).
-   mol\_numsaturatedrings(mol) : returns the number of saturated rings in a molecule (*available from 2013\_03 release*).
-   mol\_numaromaticheterocycles(mol) : returns the number of aromatic heterocycles in a molecule (*available from 2013\_03 release*).
-   mol\_numaliphaticheterocycles(mol) : returns the number of aliphatic (at least one non-aromatic bond) heterocycles in a molecule (*available from 2013\_03 release*).
-   mol\_numsaturatedheterocycles(mol) : returns the number of saturated heterocycles in a molecule (*available from 2013\_03 release*).
-   mol\_numaromaticcarbocycles(mol) : returns the number of aromatic carbocycles in a molecule (*available from 2013\_03 release*).
-   mol\_numaliphaticcarbocycles(mol) : returns the number of aliphatic (at least one non-aromatic bond) carbocycles in a molecule (*available from 2013\_03 release*).
-   mol\_numsaturatedcarbocycles(mol) : returns the number of saturated carbocycles in a molecule (*available from 2013\_03 release*).
-   mol\_inchi(mol) : returns an InChI for the molecule. (*available from the 2011\_06 release, requires that the RDKit be built with InChI support*).
-   mol\_inchikey(mol) : returns an InChI key for the molecule. (*available from the 2011\_06 release, requires that the RDKit be built with InChI support*).
-   mol\_formula(mol,bool default false, bool default true) : returns a string with the molecular formula. The second argument controls whether isotope information is included in the formula; the third argument controls whether "D" and "T" are used instead of [2H] and [3H]. (*available from the 2014\_03 release*)

##### Connectivity Descriptors

-   mol\_chi0v(mol) - mol\_chi4v(mol) : returns the ChiXv value for a molecule for X=0-4 (*available from 2012\_01 release*).
-   mol\_chi0n(mol) - mol\_chi4n(mol) : returns the ChiXn value for a molecule for X=0-4 (*available from 2012\_01 release*).
-   mol\_kappa1(mol) - mol\_kappa3(mol) : returns the kappaX value for a molecule for X=1-3 (*available from 2012\_01 release*).
-   mol\_numspiroatoms : returns the number of spiro atoms in a molecule (*available from 2015\_09 release*).
-   mol\_numbridgeheadatoms : returns the number of bridgehead atoms in a molecule (*available from 2015\_09 release*).

##### MCS

-   fmcs(mols) : an aggregation function that calculates the MCS for a set of molecules
-   fmcs\_smiles(text, json default '') : calculates the MCS for a space-separated set of SMILES. The optional json argument is used to provide parameters to the MCS code.

#### Other

-   rdkit\_version() : returns a string with the cartridge version number.

There are additional functions defined in the cartridge, but these are used for internal purposes.

## Using the Cartridge from Python

The recommended adapter for connecting to postgresql is pyscopg2 (<https://pypi.python.org/pypi/psycopg2>).

Here's an example of connecting to our local copy of ChEMBL and doing a basic substructure search:

    >>> import psycopg2
    >>> conn = psycopg2.connect(database='chembl_25')
    >>> curs = conn.cursor()
    >>> curs.execute('select * from rdk.mols where m@>%s',('c1cccc2c1nncc2',))
    >>> curs.fetchone()
    (9830, 'CC(C)Sc1ccc(CC2CCN(C3CCN(C(=O)c4cnnc5ccccc54)CC3)CC2)cc1')

That returns a SMILES for each molecule. If you plan to do more work with the molecules after retrieving them, it is much more efficient to ask postgresql to give you the molecules in pickled form:

    >>> curs.execute('select molregno,mol_send(m) from rdk.mols where m@>%s',('c1cccc2c1nncc2',))
    >>> row = curs.fetchone()
    >>> row
    (9830, <memory at 0x...>)

These pickles can then be converted into molecules:

    >>> from rdkit import Chem
    >>> m = Chem.Mol(row[1].tobytes())
    >>> Chem.MolToSmiles(m,True)
    'CC(C)Sc1ccc(CC2CCN(C3CCN(C(=O)c4cnnc5ccccc54)CC3)CC2)cc1'

## License

This document is copyright (C) 2013-2016 by Greg Landrum

This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 License. To view a copy of this license, visit <http://creativecommons.org/licenses/by-sa/4.0/> or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.

The intent of this license is similar to that of the RDKit itself. In simple words: “Do whatever you want with it, but please give us some credit.”
