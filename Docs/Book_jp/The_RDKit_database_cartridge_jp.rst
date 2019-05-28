RDKitデータベースカートリッジ
##############################################
[`The RDKit database cartridge <https://www.rdkit.org/docs/Cartridge.html#the-rdkit-database-cartridge>`__]

このページについて
****************************************************
[`What is this? <https://www.rdkit.org/docs/Cartridge.html#what-is-this>`__]

このページRDKit PostgreSQL cartridgeのチュートリアルとリファレンスガイドです。

間違いを見つけたり、より良い方法を提案する場合は、ソースドキュメント（.mdファイル）をご自身で修正していただくか、
メーリングリスト: rdkit-discuss@lists.sourceforge.net に送ってください。（メーリングリストを利用する場合はまず加入する必要があります。）

チュートリアル
****************************************************
[`Tutorial <https://www.rdkit.org/docs/Cartridge.html#tutorial>`__]

イントロダクション
======================================================
[`Introduction <https://www.rdkit.org/docs/Cartridge.html#introduction>`__]

データベースの作成
======================================================
[`Creating databases <https://www.rdkit.org/docs/Cartridge.html#creating-databases>`__

環境設定
-----------------------------------------------------
[`Configuration <https://www.rdkit.org/docs/Cartridge.html#configuration>`__]

以下のタイミングの情報は、ごく一般的なデスクトップPC（Dell Studio XPS with a 2.9GHz i7 CPU and 8GB of RAM）上で
Ubuntu 12.04を実行し、PostgreSQL v9.1.4.を用いた場合の結果です。データベースはデフォルトのパラメータでインストールされています。

データベースの読み込みとインデックスを作成している間のパフォーマンスを向上させるため、postgresql.confにあるpostgresqlの環境設定を2、3変更しました。

.. code:: python

   synchronous_commit = off      # immediate fsync at commit
   full_page_writes = off            # recover from partial

また、探索の性能を向上させるため、postgresqlが使用可能なメモリを、デフォルトの極めて控えめな設定よりも多くしています。

.. code:: python

   shared_buffers = 2048MB           # min 128kB
                     # (change requires restart)
   work_mem = 128MB              # min 64kB

ファイルからデータベースを作成する
-----------------------------------------------------
[`Creating a database from a file <https://www.rdkit.org/docs/Cartridge.html#creating-a-database-from-a-file>`__]

この例では、emolecules.comのURL http://www.emolecules.com/doc/plus/download-database.php からダウンロード可能な、
購入可能な化合物に関するSMILESファイルからデータベースを読み込む方法をお見せします。

この例をそっくりそのままご自身で繰り返してみようとされる場合、600万行のデータベースの読み込みと、全てのフィンガープリンの生成に数時間要することには気をつけてください。

まず最初にデータベースを作成し、カートリッジをインストールします。

.. code:: sql

   ~/RDKit_trunk/Data/emolecules > createdb emolecules
   ~/RDKit_trunk/Data/emolecules > psql -c 'create extension rdkit' emolecules

次に、生のデータを持つテーブルを作成し、データを格納します:

.. code:: sql

   ~/RDKit_trunk/Data/emolecules > psql -c 'create table raw_data (id SERIAL, smiles text, emol_id integer, parent_id integer)' emolecules
   NOTICE:  CREATE TABLE will create implicit sequence "raw_data_id_seq" for serial column "raw_data.id"
   CREATE TABLE
   ~/RDKit_trunk/Data/emolecules > zcat emolecules-2013-02-01.smi.gz | sed '1d; s/\\/\\\\/g' | psql -c "copy raw_data (smiles,emol_id,parent_id) from stdin with delimiter ' '" emolecules

次に分子のテーブルを作成しますが、RDKitが処理できるSMILESのみのテーブルとなります:

.. code:: sql

   ~/RDKit_trunk/Data/emolecules > psql emolecules
   psql (9.1.4)
   Type "help" for help.
   emolecules=# select * into mols from (select id,mol_from_smiles(smiles::cstring) m from raw_data) tmp where m is not null;
   WARNING:  could not create molecule from SMILES 'CN(C)C(=[N+](C)C)Cl.F[P-](F)(F)(F)(F)F'
   ... a lot of warnings deleted ...
   SELECT 6008732
   emolecules=# create index molidx on mols using gist(m);
   CREATE INDEX

最後のステップは、部分構造検索をするつもりの時のみ必要となります。

ChEMBLの読み込み
-----------------------------------------------------
[`Loading ChEMBL <https://www.rdkit.org/docs/Cartridge.html#loading-chembl>`__]

ChEMBLのウェブサイト ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest から
postsqlgreダンプをダウンロードし、インストールすることから始めます。

データベースに接続、カートリッジをインストールし、そしてこれから使うスキーマを作成します:

.. code:: sql

   chembl_23=# create extension if not exists rdkit;
   chembl_23=# create schema rdk;

分子を作成し、部分構造検索のインデックスを構築します:

.. code:: sql

   chembl_23=# select * into rdk.mols from (select molregno,mol_from_ctab(molfile::cstring) m  from compound_structures) tmp where m is not null;
   SELECT 1727081
   chembl_23=# create index molidx on rdk.mols using gist(m);
   CREATE INDEX
   chembl_23=# alter table rdk.mols add primary key (molregno);
   ALTER TABLE

フィンガープリントをいくつか作成し、類似度検索のインデックスを構築します:

.. code:: sql

   chembl_23=# select molregno,torsionbv_fp(m) as torsionbv,morganbv_fp(m) as mfp2,featmorganbv_fp(m) as ffp2 into rdk.fps from rdk.mols;
   SELECT 1727081
   chembl_23=# create index fps_ttbv_idx on rdk.fps using gist(torsionbv);
   CREATE INDEX
   chembl_23=# create index fps_mfp2_idx on rdk.fps using gist(mfp2);
   CREATE INDEX
   chembl_23=# create index fps_ffp2_idx on rdk.fps using gist(ffp2);
   CREATE INDEX
   chembl_23=# alter table rdk.fps add primary key (molregno);
   ALTER TABLE

psqlプロンプトにペーストするだけでいいように、以下の一つのブロックの中に、ここ（とこれ以降）で使ったコマンドをまとめておきます:

.. code:: sql

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
   returns table(molregno integer, m mol, similarity double precision) as
   $$
   select molregno,m,tanimoto_sml(morganbv_fp(mol_from_smiles($1::cstring)),mfp2) as similarity
   from rdk.fps join rdk.mols using (molregno)
   where morganbv_fp(mol_from_smiles($1::cstring))%mfp2
   order by morganbv_fp(mol_from_smiles($1::cstring))<%>mfp2;
   $$ language sql stable ;

部分構造検索
======================================================
[`Substructure searches <https://www.rdkit.org/docs/Cartridge.html#substructure-searches>`__]

例として用いたクエリの分子は\ `eMolecules home Page <http://www.emolecules.com/>`__\ から取得しました:

.. code:: sql

   chembl_23=# select count(*) from rdk.mols where m@>'c1cccc2c1nncc2' ;
    count
   -------
      447
   (1 row)

   Time: 107.602 ms
   chembl_23=# select count(*) from rdk.mols where m@>'c1ccnc2c1nccn2' ;
    count
   -------
     1013
   (1 row)

   Time: 216.222 ms
   chembl_23=# select count(*) from rdk.mols where m@>'c1cncc2n1ccn2' ;
    count
   -------
     1775
   (1 row)

   Time: 88.266 ms
   chembl_23=# select count(*) from rdk.mols where m@>'Nc1ncnc(N)n1' ;
    count
   -------
     5842
   (1 row)

   Time: 327.855 ms
   chembl_23=# select count(*) from rdk.mols where m@>'c1scnn1' ;
    count
   -------
    15962
   (1 row)

   Time: 568.675 ms
   chembl_23=# select count(*) from rdk.mols where m@>'c1cccc2c1ncs2' ;
    count
   -------
    18986
   (1 row)

   Time: 998.104 ms
   chembl_23=# select count(*) from rdk.mols where m@>'c1cccc2c1CNCCN2' ;
    count
   -------
     1613
   (1 row)

   Time: 1922.273 ms

最後の２つのクエリでは、実行し全ての結果の数を数え始めるまでに時間がかかったことに注意してください。

170万化合物全体から検索していることを考えれば、これらの検索時間が信じられないほど遅いというわけではありませんが、もっと速くなるならそれに越したことはないでしょう。

特に大きな数の結果を返すクエリにおいて、検索を加速する簡単な方法のひとつは、限られた数の結果のみを取得することです:

.. code:: sql

   chembl_23=# select * from rdk.mols where m@>'c1cccc2c1CNCCN2' limit 100;
    molregno |                                                                                             m

   ----------+-----------------------------------------------------------------------------------------------------------------------------------------------------------
   --------------------------------
      908048 | O=C1CN(C(=O)c2ccc(Br)o2)C(c2ccc(F)cc2)c2cc(F)ccc2N1
      931972 | Cl.c1ccc(CC2CNc3ccccc3CN2)cc1
      904450 | CCOC(=O)[C@H]1[C@H]2COc3ccc(Cl)cc3[C@@H]2N2C(=O)c3ccc(Cl)cc3NC(=O)[C@@]12C
      226391 | C/C=C1/CC2C(OC)Nc3cc(OC)c(OC)cc3C(=O)N2C1
      930820 | CN1CC(=O)N(CC(=O)Nc2ccc(N(C)C)cc2)c2ccccc2C1=O
       18576 | CO[C@H]1Nc2c(ccc(C)c2O)C(=O)N2C=C(/C=C/C(N)=O)C[C@@H]12
      249934 | O=C(c1cccc2ccccc12)N1CCN(Cc2cncn2Cc2ccccc2)c2ccccc2C1
      ...
       91020 | CC(C)C[C@H]1C(=O)N2c3ccccc3[C@@](O)(C[C@@H]3NC(=O)c4ccccc4N4C(=O)c5ccccc5NC34)[C@H]2N1C(=O)C(CCCNC(=O)OCc1ccccc1)NC(=O)OC(C)(C)C
       91225 | CC(C)C[C@H]1C(=O)N2c3ccccc3[C@@](O)(C[C@@H]3NC(=O)c4ccccc4N4C(=O)c5ccccc5NC34)[C@H]2N1C(=O)CCC(=O)[O-].[Na+]
      348798 | O=C(O)CN1C(=O)C(c2ccc(Cl)cc2)N(C(C(=O)O)c2ccc(Cl)cc2)C(=O)c2cc(I)ccc21
      348972 | C[C@H](c1ccc(Cl)cc1)N1C(=O)c2cc(I)ccc2N(CCCCC(=O)O)C(=O)[C@@H]1c1ccc(C(F)(F)F)cc1

   ...skipping 23 lines
   Time: 97.357 ms

SMARTSベースのクエリ
-----------------------------------------------------
[`SMARTS-based queries <https://www.rdkit.org/docs/Cartridge.html#smarts-based-queries>`__]

オキサジアゾールあるいはチアジアゾール:

.. code:: sql

   chembl_23=# select * from rdk.mols where m@>'c1[o,s]ncn1'::qmol limit 500;
    molregno |                                                      m
   ----------+--------------------------------------------------------------------------------------------------------------
     1370170 | Fc1cccc(-c2nc(NCC3COc4ccccc4O3)no2)c1F
     1370417 | COc1cc(CN2CCC(Cc3nc(-c4ccc5c(c4)CCO5)no3)C2)ccc1F
     1370526 | Cl.Cn1cc(-c2noc(/C=C3/CCN4CCCC[C@@H]4C3)n2)c2ccccc21
     1379267 | CCC(c1ccccc1)c1noc(CCN(CC)CC)n1
     1404150 | OC[C@H]1O[C@H](c2nc(-c3nc(-c4cccs4)no3)cs2)C[C@@H]1O
     1217463 | CC(C)(C)c1ccc(-c2noc(CCC(=O)N3CCCCC3)n2)cc1
     ...
     1517753 | CC(C)c1noc(N2CCC(CO[C@H]3CC[C@H](c4ccc(S(C)(=O)=O)cc4F)CC3)CC2)n1
     1263024 | COc1cc(Nc2nc3c(s2)CCCC3c2ccccc2)ccc1-c1nc(C)no1
     1264016 | O=C(O)CCc1nc2cc(-c3noc(-c4cc(C(F)(F)F)cc(C(F)(F)F)c4)n3)ccc2[nH]1
     1847733 | Cc1cc(-c2noc([C@H]3CCCCN3C(=O)COc3ccccc3)n2)no1
   (500 rows)

   Time: 761.847 ms

純粋なSMILESのクエリと比較すると遅いですが、これはSMARTSベースのクエリには一般的に当てはまることです。

立体化学の使用方法
-----------------------------------------------------
[`Using Stereochemistry <https://www.rdkit.org/docs/Cartridge.html#using-stereochemistry>`__]

デフォルトでは部分構造クエリを実行する際に立体化学が考慮されないことに気をつけてください:

.. code:: sql

   chembl_23=# select * from rdk.mols where m@>'NC(=O)[C@@H]1CCCN1C=O' limit 10;
    molregno |
                  m

   ----------+-----------------------------------------------------------------------------------------------------------------------------------------------------------
   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
   ---------------------
       87611 | CNCC(=O)N[C@@H](CCCN=C(N)N)C(=O)N1C[C@H]2C[C@H]1C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@H](C(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](Cc1ccccc1)C(=O
   )O)CCSS2
       88372 | CNCCCC[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](CCCCNC)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CO)NC(=O)[C@@H](Cc1ccccc1)NC(=O)[C@@H](Cc1ccccc1)NC(=O)[C@@H](C
   c1ccc2ccccc2c1)NC(C)=O)C(=O)N1CCC[C@@H]1C(=O)N[C@H](C)C(=O)O
       88322 | CC(=O)N[C@H](Cc1ccc2ccccc2c1)C(=O)N[C@H](Cc1ccccc1)C(=O)N[C@H](Cc1ccccc1)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@H](CCCCNC(C)C)C(=O)N[C@@H](Cc1
   ccccc1)C(=O)N[C@@H](CCCCNC(C)C)C(=O)N1CCC[C@@H]1C(=O)N[C@H](C)C(=O)O
       88168 | CC(=O)N[C@H](Cc1ccc2ccccc2c1)C(=O)N[C@H](Cc1ccccc1)C(=O)N[C@H](Cc1ccccc1)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@H](CCCN=C(N)N)C(=O)N[C@@H](Cc1
   ccccc1)C(=O)N[C@@H](CCCCNC1CCCC1)C(=O)N1CCC[C@@H]1C(=O)N[C@H](C)C(=O)O
       88150 | CC(=O)N[C@H](Cc1ccc2ccccc2c1)C(=O)N[C@H](Cc1ccccc1)C(=O)N[C@H](Cc1ccccc1)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@H](CCCN=C(N)N)C(=O)N[C@@H](Cc1
   ccccc1)C(=O)N[C@@H](CCCCNCc1ccc(C)cc1)C(=O)N1CCC[C@@H]1C(=O)N[C@H](C)C(=O)O
       88373 | CC(=O)N[C@H](Cc1ccc2ccccc2c1)C(=O)N[C@H](Cc1ccccc1)C(=O)N[C@H](Cc1ccccc1)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@H](CCCCNC1CCCCC1)C(=O)N[C@@H](
   Cc1ccccc1)C(=O)N[C@@H](CCCCNC1CCCCC1)C(=O)N1CCC[C@@H]1C(=O)N[C@H](C)C(=O)O
       93377 | CC(=O)N[C@@H](Cc1ccc([N+](=O)[O-])cc1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCC/N=C(/N)NS(=O)(=O)c1c(C)c(C)c2c(c1C)CCC(C)(C)O2)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](
   CCC/N=C(/N)NS(=O)(=O)c1c(C)c(C)c2c(c1C)CCC(C)(C)O2)C(=O)N[C@H](C(=O)NCC(=O)N[C@@H](COC(C)(C)C)C(=O)N[C@@H](CCCCNC(=O)c1ccccc1N)C(=O)NCC(=O)O)[C@@H](C)OC(C)(C)C
       94493 | CC(C)C[C@@H]1NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](C(C)C)NC(=O)[C@H](NC(=O)[C@H](CCCCN)NC(=O)[C@@H]2CCCN2C(=O)[C@H](CCC(N)=O)NC
   (=O)CNC(=O)CN)CSSC[C@@H](C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)N[C@@H](CO)C(=O)N[C@H](C(=O)NCC(=O)NCC(N)=O)[C@@H](C)O)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](Cc2ccccc2)NC(=O)CNC
   (=O)[C@@H]2CCCN2C1=O

   ...skipping 1 line
       89559 | CC1(C)SSC(C)(C)[C@@H](C(=O)N[C@@H](Cc2c[nH]cn2)C(=O)N2CCC[C@H]2C(=O)N[C@@H](Cc2ccccc2)C(=O)O)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H]1NC(=O)[C@H](CCCN=C(N)N)N
   C(=O)[C@@H](N)CC(=O)O
   (10 rows)

rdkit.do_chiral_sss構成変数を使うことで設定を変えることができます:

.. code:: sql

   chembl_23=# set rdkit.do_chiral_sss=true;
   SET
   Time: 0.241 ms
   chembl_23=# select * from rdk.mols where m@>'NC(=O)[C@@H]1CCCN1C=O' limit 10;
    molregno |
               m

   ----------+--------------------------------------------------------------------------------------------------------------------------------------------------------------
   -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   ---------------
       87611 | CNCC(=O)N[C@@H](CCCN=C(N)N)C(=O)N1C[C@H]2C[C@H]1C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@H](C(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](Cc1ccccc1)C(=O)O)
   CCSS2
       93377 | CC(=O)N[C@@H](Cc1ccc([N+](=O)[O-])cc1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCC/N=C(/N)NS(=O)(=O)c1c(C)c(C)c2c(c1C)CCC(C)(C)O2)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCC
   /N=C(/N)NS(=O)(=O)c1c(C)c(C)c2c(c1C)CCC(C)(C)O2)C(=O)N[C@H](C(=O)NCC(=O)N[C@@H](COC(C)(C)C)C(=O)N[C@@H](CCCCNC(=O)c1ccccc1N)C(=O)NCC(=O)O)[C@@H](C)OC(C)(C)C
       94493 | CC(C)C[C@@H]1NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](C(C)C)NC(=O)[C@H](NC(=O)[C@H](CCCCN)NC(=O)[C@@H]2CCCN2C(=O)[C@H](CCC(N)=O)NC(=O
   )CNC(=O)CN)CSSC[C@@H](C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)N[C@@H](CO)C(=O)N[C@H](C(=O)NCC(=O)NCC(N)=O)[C@@H](C)O)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](Cc2ccccc2)NC(=O)CNC(=O)[C
   @@H]2CCCN2C1=O
       89558 | NC(N)=NCCC[C@H](NC(=O)[C@@H](N)CC(=O)O)C(=O)N[C@H]1CCSSC[C@@H](C(=O)N[C@@H](Cc2c[nH]cn2)C(=O)N2CCC[C@H]2C(=O)N[C@@H](Cc2ccccc2)C(=O)O)NC(=O)[C@H](Cc2ccc(O)cc
   2)NC1=O
       89559 | CC1(C)SSC(C)(C)[C@@H](C(=O)N[C@@H](Cc2c[nH]cn2)C(=O)N2CCC[C@H]2C(=O)N[C@@H](Cc2ccccc2)C(=O)O)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H]1NC(=O)[C@H](CCCN=C(N)N)NC(=
   O)[C@@H](N)CC(=O)O
      126618 | NC(=O)[C@@H]1CCCN1C(=O)[C@@H]1CCCN1C(=O)[C@@H](O)[C@H](N)Cc1ccccc1
      152339 | O=C(O)CN[C@H](CC1CCCCC1)C(=O)N1CCC[C@H]1C(=O)NCCCc1c[nH]cn1
      152504 | N[C@H](CC1CCCCC1)C(=O)N1[C@H](C(=O)NC/C=C/c2c[nH]cn2)C[C@@H]2CCCC[C@@H]21
      152383 | N[C@H](CC1CCCCC1)C(=O)N1CCC[C@H]1C(=O)NCCCCc1c[nH]cn1
      151837 | N[C@H](CC1CCCCC1)C(=O)N1CCC[C@H]1C(=O)NC/C=C/c1c[nH]cn1
   (10 rows)

   Time: 6.181 ms

クエリの調整
-----------------------------------------------------
[`Tuning queries <https://www.rdkit.org/docs/Cartridge.html#tuning-queries>`__]

複雑なSMARTSクエリを構築する必要無しに、部分構造クエリをもう少しコントロールすることができれば、しばしば役に立つでしょう。
カートリッジの関数\ ``mol_adjust_query_properties()``\ はちょうどこの目的のため使うことができます。
2,6-二置換ピリジンのクエリを使って、デフォルトの場合どのような動きを示すか例をお見せします：

.. code:: sql

   chembl_23=# select molregno,m from rdk.mols where m@>mol_adjust_query_properties('*c1cccc(NC(=O)*)n1') limit 10;
    molregno |                                             m
   ----------+-------------------------------------------------------------------------------------------
     1993749 | Cn1c(Nc2c(Cl)ccc(CNC(=O)C(C)(C)C)c2Cl)nc2cc(C(=O)Nc3cccc(C(F)(F)F)n3)c(N3CCC(F)(F)C3)cc21
     1988455 | Cc1cccc(C(=O)Nc2cccc(Oc3cccnc3)n2)c1
     1870095 | COC(=O)CN(C(=O)C(C)c1c(F)cccc1F)c1cccc(C)n1
     1870023 | CCC(C)C(=O)N(CC(=O)OC)c1cccc(C)n1
     1873944 | Cc1ccc(C(=O)N(C)CC(=O)Nc2cccc(C)n2)cn1
     1873968 | Cc1cccc(NC(=O)CN(C)C(=O)c2ccc(-n3cccc3)nc2)n1
     1882693 | Cc1cccc(NC(=O)CCNCc2c(C)nn(C)c2N(C)C)n1
     1882711 | COc1c(CNCCC(=O)Nc2cccc(C)n2)c(C)nn1C
     1868705 | CCOc1cccc(NC(=O)c2cnc(C)cn2)n1
     1875177 | Cc1cccc(NC(=O)[C@@H]2CCCN2Cc2nc(C)c(C)o2)n1
   (10 rows)

   Time: 11.895 ms

デフォルトでは\ ``mol_adjust_query_properties()``\ は次の変更を分子に施します:

- ダミーアトムを“any”クエリに変換する
- 環の原子全てに級クエリを追加し、置換が与えられているものと合致するようにする
- 芳香族性の認識が行われる（まだ行われていない場合）

追加のJSON引数を与えることで、動きを制御することができます。ここでは追加の級クエリを機能させなくする例を示します：

.. code:: sql

   chembl_23=# select molregno,m from rdk.mols where m@>mol_adjust_query_properties('*c1cccc(NC(=O)*)n1',
   chembl_23(# '{"adjustDegree":false}') limit 10;
    molregno |                                             m
   ----------+-------------------------------------------------------------------------------------------
     1993749 | Cn1c(Nc2c(Cl)ccc(CNC(=O)C(C)(C)C)c2Cl)nc2cc(C(=O)Nc3cccc(C(F)(F)F)n3)c(N3CCC(F)(F)C3)cc21
     1957849 | COc1ccc2ncc(F)c(C[C@H](O)C3CCC(NCc4nc5c(cc4F)OCC(=O)N5)CO3)c2n1
     1959611 | O=C1COc2ccc(CNC3CCN(CCn4c(=O)ccc5ncc(OCc6cccnn6)cc54)CC3)nc2N1
     1988455 | Cc1cccc(C(=O)Nc2cccc(Oc3cccnc3)n2)c1
     1870095 | COC(=O)CN(C(=O)C(C)c1c(F)cccc1F)c1cccc(C)n1
     1870023 | CCC(C)C(=O)N(CC(=O)OC)c1cccc(C)n1
     1873944 | Cc1ccc(C(=O)N(C)CC(=O)Nc2cccc(C)n2)cn1
     1873968 | Cc1cccc(NC(=O)CN(C)C(=O)c2ccc(-n3cccc3)nc2)n1
     1882693 | Cc1cccc(NC(=O)CCNCc2c(C)nn(C)c2N(C)C)n1
     1882711 | COc1c(CNCCC(=O)Nc2cccc(C)n2)c(C)nn1C
   (10 rows)

   Time: 10.780 ms

あるいは、環構造の原子にもダミー（鎖状構造の原子と結合しているもの）にも追加の級クエリを付加しない例を示します:

.. code:: sql

   chembl_23=# select molregno,m from rdk.mols where m@>mol_adjust_query_properties('*c1cccc(NC(=O)*)n1',
   chembl_23(# '{"adjustDegree":true,"adjustDegreeFlags":"IGNORERINGS|IGNOREDUMMIES"}') limit 10;
    molregno |                                             m
   ----------+-------------------------------------------------------------------------------------------
     1993749 | Cn1c(Nc2c(Cl)ccc(CNC(=O)C(C)(C)C)c2Cl)nc2cc(C(=O)Nc3cccc(C(F)(F)F)n3)c(N3CCC(F)(F)C3)cc21
     1957849 | COc1ccc2ncc(F)c(C[C@H](O)C3CCC(NCc4nc5c(cc4F)OCC(=O)N5)CO3)c2n1
     1959611 | O=C1COc2ccc(CNC3CCN(CCn4c(=O)ccc5ncc(OCc6cccnn6)cc54)CC3)nc2N1
     1988455 | Cc1cccc(C(=O)Nc2cccc(Oc3cccnc3)n2)c1
     1873944 | Cc1ccc(C(=O)N(C)CC(=O)Nc2cccc(C)n2)cn1
     1873968 | Cc1cccc(NC(=O)CN(C)C(=O)c2ccc(-n3cccc3)nc2)n1
     1882693 | Cc1cccc(NC(=O)CCNCc2c(C)nn(C)c2N(C)C)n1
     1882711 | COc1c(CNCCC(=O)Nc2cccc(C)n2)c(C)nn1C
     1884388 | Cc1noc(COCC(=O)Nc2ccc(Br)c(C)n2)n1
     1868705 | CCOc1cccc(NC(=O)c2cnc(C)cn2)n1
   (10 rows)

   Time: 12.827 ms

利用可能なオプションは以下の通りです:

-  **adjustDegree** (default: true) :
   入力の原子の級とマッチするようにクエリを付与する
-  **adjustDegreeFlags** (default: ADJUST_IGNOREDUMMIES \|
   ADJUST_IGNORECHAINS) 級が調整されている場所をコントロールする
-  **adjustRingCount** (default: false) :
   入力の環構造の数にマッチするようにクエリを付与する
-  **adjustRingCountFlags** (default: ADJUST_IGNOREDUMMIES \|
   ADJUST_IGNORECHAINS) 環構造の数が調整されている場所をコントロールする
-  **makeDummiesQueries** (default: true) :
   入力構造のダミーアトムをany-atomクエリに変換する
-  **aromatizeIfPossible** (default: true) :
   入力構造に芳香族性を認識するアルゴリズムの処理を行う（注:
   SMILESから構築された分子はいつも芳香族の認識処理が行われているので、これは多くの場合余分な処理です）
-  **makeBondsGeneric** (default: false) :
   結合をany-bondクエリに変換する
-  **makeBondsGenericFlags** (default: false) :
   どの結合を一般的なものにするかコントロールする
-  **makeAtomsGeneric** (default: false) :
   原子をany-atomクエリに変換する
-  **makeAtomsGenericFlags** (default: false) :
   どの原子を一般的なものとするかコントロールする

上で述べた、特定のオプションが適用されている場所をコントロールするさまざまな\ ``Flags``\ 引数は,
下のリストにある操作を\ ``|``\ 記号で連結することで構築されています。

-  **IGNORENONE** : 操作を全ての原子に適用する
-  **IGNORERINGS** : 環構造の原子に操作を適用しない
-  **IGNORECHAINS** : 鎖状構造の原子に操作を適用しない
-  **IGNOREDUMMIES** : ダミー原子に操作を適用しない
-  **IGNORENONDUMMIES** : 非ダミー原子に操作を適用しない
-  **IGNOREALL** : anyアトムに操作を適用しない

類似度検索
======================================================
[`Similarity searches <https://www.rdkit.org/docs/Cartridge.html#similarity-searches>`__]

基本的な類似度検索:

.. code:: sql

   chembl_23=# select count(*) from rdk.fps where mfp2%morganbv_fp('Cc1ccc2nc(-c3ccc(NC(C4N(C(c5cccs5)=O)CCC4)=O)cc3)sc2c1');
    count
   -------
       67
   (1 row)

   Time: 177.579 ms

一般的に我々は、付随するSMILESと一緒に、近傍の分子のソートされた一覧を見つけ出したいと思っています。
そのような場合、このSQL関数をつかうと簡単にできます。

.. code:: sql

   chembl_23=# create or replace function get_mfp2_neighbors(smiles text)
       returns table(molregno integer, m mol, similarity double precision) as
     $$
     select molregno,m,tanimoto_sml(morganbv_fp(mol_from_smiles($1::cstring)),mfp2) as similarity
     from rdk.fps join rdk.mols using (molregno)
     where morganbv_fp(mol_from_smiles($1::cstring))%mfp2
     order by morganbv_fp(mol_from_smiles($1::cstring))<%>mfp2;
     $$ language sql stable ;
   CREATE FUNCTION
   Time: 0.856 ms
   chembl_23=# select * from get_mfp2_neighbors('Cc1ccc2nc(-c3ccc(NC(C4N(C(c5cccs5)=O)CCC4)=O)cc3)sc2c1') limit 10;
    molregno |                             m                              |    similarity
   ----------+------------------------------------------------------------+-------------------
      471319 | Cc1ccc2nc(-c3ccc(NC(=O)C4CCN(S(=O)(=O)c5cccs5)C4)cc3)sc2c1 | 0.638888888888889
     1032469 | O=C(Nc1nc2ccc(Cl)cc2s1)[C@@H]1CCCN1C(=O)c1cccs1            | 0.623188405797101
      751668 | COc1ccc2nc(NC(=O)[C@@H]3CCCN3C(=O)c3cccs3)sc2c1            | 0.619718309859155
      471318 | Cc1ccc2nc(-c3ccc(NC(=O)C4CN(S(=O)(=O)c5cccs5)C4)cc3)sc2c1  | 0.611111111111111
      740754 | Cc1ccc(NC(=O)C2CCCN2C(=O)c2cccs2)cc1C                      | 0.606060606060606
      732905 | O=C(Nc1ccc(S(=O)(=O)N2CCCC2)cc1)C1CCCN1C(=O)c1cccs1        | 0.602941176470588
     1087495 | Cc1ccc(NC(=O)C2CCCN2C(=O)c2cccs2)c(C)c1                    | 0.597014925373134
      471462 | CCS(=O)(=O)N1CCC(C(=O)Nc2ccc(-c3nc4ccc(C)cc4s3)cc2)CC1     | 0.585714285714286
      810850 | Cc1cc(C)n(-c2ccc(NC(=O)C3CCCCN3C(=O)c3cccs3)cc2)n1         | 0.583333333333333
     1224407 | O=C(Nc1cccc(S(=O)(=O)N2CCCC2)c1)C1CCCN1C(=O)c1cccs1        | 0.579710144927536
   (10 rows)

   Time: 28.909 ms
   chembl_23=# select * from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1') limit 10;
    molregno |                           m                           |    similarity
   ----------+-------------------------------------------------------+-------------------
     1044892 | Cc1ccc2nc(N(CCN(C)C)C(=O)c3cc(Cl)sc3Cl)sc2c1          | 0.518518518518518
     1040496 | Cc1ccc2nc(N(CCCN(C)C)C(=O)CCc3ccccc3)sc2c1            | 0.517857142857143
     1049393 | Cc1ccc2nc(N(CCCN(C)C)C(=O)CS(=O)(=O)c3ccccc3)sc2c1    | 0.517857142857143
      441378 | Cc1ccc2nc(NC(=O)CCC(=O)O)sc2c1                        | 0.510204081632653
     1047691 | Cc1ccc(S(=O)(=O)CC(=O)N(CCCN(C)C)c2nc3ccc(C)cc3s2)cc1 | 0.509090909090909
      911501 | Cc1ccc2nc(N(CCN(C)C)C(=O)c3cc(Cl)sc3Cl)sc2c1.Cl       | 0.509090909090909
     1042958 | Cc1ccc2nc(N(CCN(C)C)C(=O)c3ccc4ccccc4c3)sc2c1         | 0.509090909090909
      775269 | Cc1ccc2nc(N(CCCN(C)C)C(=O)CCc3ccccc3)sc2c1.Cl         | 0.508771929824561
     1045663 | Cc1ccc2nc(N(CCCN(C)C)C(=O)COc3ccc(Cl)cc3)sc2c1        |               0.5
     1015485 | Cc1ccc2nc(N(Cc3cccnc3)C(=O)Cc3ccccc3)sc2c1            |               0.5
   (10 rows)

   Time: 41.623 ms

類似度のカットオフの調整
-----------------------------------------------------
[`Adjusting the similarity cutoff <https://www.rdkit.org/docs/Cartridge.html#adjusting-the-similarity-cutoff>`__]

デフォルトでは類似度検索で返される最小の類似度は0.5です。
この値はrdkit.tanimoto_threshold（及びrdkit.dice_threshold）構成変数により調整することができます:

.. code:: sql

   chembl_23=# select count(*) from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1');
    count
   -------
       20
   (1 row)

   Time: 181.438 ms
   chembl_23=# set rdkit.tanimoto_threshold=0.7;
   SET
   Time: 0.047 ms
   chembl_23=# select count(*) from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1');
    count
   -------
        0
   (1 row)

   Time: 161.228 ms
   chembl_23=# set rdkit.tanimoto_threshold=0.6;
   SET
   Time: 0.045 ms
   chembl_23=# select count(*) from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1');
    count
   -------
        1
   (1 row)

   Time: 184.275 ms
   chembl_23=# set rdkit.tanimoto_threshold=0.5;
   SET
   Time: 0.055 ms
   chembl_23=# select count(*) from get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1');
    count
   -------
       20
   (1 row)

   Time: 181.100 ms

MCSコードの使い方
======================================================
[`Using the MCS Code <https://www.rdkit.org/docs/Cartridge.html#using-the-mcs-code>`__]

最も簡単なMCSコードの使い方は、分子のグループから最大共通部分構造を見つけ出すことです:

.. code:: sql

   chembl_23=# select fmcs(m::text) from rdk.mols join compound_records using (molregno) where doc_id=4;
                                     fmcs
   ------------------------------------------------------------------------
    [#6](-[#6]-[#7]-[#6]-[#6](-,:[#6])-,:[#6])-,:[#6]-,:[#6]-,:[#6]-,:[#6]
   (1 row)

   Time: 31.041 ms
   chembl_23=# select fmcs(m::text) from rdk.mols join compound_records using (molregno) where doc_id=5;
                                                                      fmcs
   ------------------------------------------------------------------------------------------------------------------------------------------
    [#6]-[#6](=[#8])-[#7]-[#6](-[#6](=[#8])-[#7]1-[#6]-[#6]-[#6]-[#6]-1-[#6](=[#8])-[#7]-[#6](-[#6](=[#8])-[#8])-[#6]-[#6])-[#6](-[#6])-[#6]
   (1 row)

   Time: 705.535 ms

同じことがSMILESのカラムでできます：

.. code:: sql

   chembl_23=# select fmcs(canonical_smiles) from compound_structures join compound_records using (molregno) where doc_id=4;
                                     fmcs
   ------------------------------------------------------------------------
    [#6](-[#7]-[#6]-[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6])-[#6](-,:[#6])-,:[#6]
   (1 row)

   Time: 128.879 ms

このページを書いている時点（2017_03リリース）では幾分骨の折れる作業ですが、パラメータのいくつかをFMCSアルゴリズムに合わせることも可能です。

.. code:: sql

   chembl_23=# select fmcs_smiles(str,'{"Threshold":0.8}') from
   chembl_23-#    (select string_agg(m::text,' ') as str from rdk.mols
   chembl_23(#    join compound_records using (molregno) where doc_id=4) as str ;

                                                                              fmcs_smiles
   ------------------------------------------------------------------------------------------------------------------------------------------------------------------
    [#6]-[#6]-[#8]-[#6](-[#6](=[#8])-[#7]-[#6](-[#6])-[#6](-,:[#6])-,:[#6])-[#6](-[#8])-[#6](-[#8])-[#6](-[#8]-[#6]-[#6])-[#6]-[#7]-[#6](-[#6])-[#6](-,:[#6])-,:[#6]
   (1 row)

   Time: 9673.949 ms
   chembl_23=#
   chembl_23=# select fmcs_smiles(str,'{"AtomCompare":"Any"}') from
   chembl_23-#    (select string_agg(m::text,' ') as str from rdk.mols
   chembl_23(#    join compound_records using (molregno) where doc_id=4) as str ;
                                                                                 fmcs_smiles
   ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    [#6]-,:[#6,#7]-[#8,#6]-[#6,#7](-[#6,#8]-[#7,#6]-,:[#6,#7]-,:[#6,#7]-,:[#7,#6]-,:[#6])-[#6,#7]-[#6]-[#6](-[#8,#6]-[#6])-[#6,#7]-[#7,#6]-[#6]-,:[#6,#8]-,:[#7,#6]-,:[#6]
   (1 row)

   Time: 304.332 ms

*注* \ ``"AtomCompare":"Any"``\ と1.0よりも小さい\ ``"Threshold"``\ の組み合わせは、極めて一般的な検索となり、とても長い検索時間がかかるという結果になりうる可能性があります。
この組み合わせを使う場合は\ ``"Tmieout"``\ を使うことを推奨します:

.. code:: sql

   chembl_23=# select fmcs_smiles(str,'{"AtomCompare":"Any","CompleteRingsOnly":true,"Threshold":0.8,"Timeout":60}') from
   chembl_23-#    (select string_agg(m::text,' ') as str from rdk.mols
   chembl_23(#    join compound_records using (molregno) where doc_id=3) as str ;

   WARNING:  findMCS timed out, result is not maximal
                                                                                             fmcs_smiles

   -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   ----------------------
    [#8]=[#6](-[#7]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6](=[#8])-[#7]1-[#6]-[#6]-[#6]-[#6,#7]-[#6]2:[#6]-1:[#6]:[#6]:[#16]:2)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]1:[#6]:
   [#6]:[#6]:[#6]:[#6]:1
   (1 row)

   Time: 60479.753 ms

利用可能なパラメータとデフォルトの値は以下の通りです:

-  MaximizeBonds (true)
-  Threshold (1.0)
-  Timeout (-1, no timeout)
-  MatchValences (false)
-  MatchChiralTag (false) :  原子に適用
-  RingMatchesRingOnly (false)
-  CompleteRingsOnly (false)
-  MatchStereo (false) : 結合に適用
-  AtomCompare (“Elements”) : 選択肢は“Elements”、“Isotopes”あるいは“Any”
-  BondCompare (“Order”) : 選択肢は“Order”、“OrderExact”あるいは“Any”

レファレンスガイド
****************************************************
[`Reference Guid <https://www.rdkit.org/docs/Cartridge.html#reference-guide>`__]

新しい型
======================================================
[`New Types <https://www.rdkit.org/docs/Cartridge.html#new-types>`__]

-  mol : rdkitにおける分子（rdkit molecule）。
   SMILESから直接、型を変更することで作ることができます。例えば：\ ``‘c1ccccc1’::mol``\　は、SMILES\ ``‘c1ccccc1’`` で表される分子を作成します。
-  qmol
   :クエリ機能を含むrdkitの分子(すなわち、SMARTSから構築されたもの)。SMARTSから直接、型を変更することで作ることができます。
   例えば:\ ``‘c1cccc[c,n]1’::qmol``\ はSMARTS\ ``‘c1cccc[c,n]1’``\ で表されるクエリ分子を作成します。
-  sfp : スパースカウントベクトルフィンガープリント（sparse count vector fingerprint）（C++とPythonのSparseIntVect）
-  bfp : ビットベクトルフィンガープリント（bit vector fingerprint）（C++とPythonのExplicitBitVect）

パラメータ
======================================================
[`Parameters <https://www.rdkit.org/docs/Cartridge.html#parameters>`__]

-  rdkit.tanimoto_threshold :
   タニモト類似度演算子のための閾値。タニモト類似度を使った検索では、少なくともこの閾値以上の類似度をもつ結果のみ返されます。
-  rdkit.dice_threshold :
   Dice類似度演算子のための閾値。Dice類似度を使った検索では、少なくともこの閾値以上の類似度をもつ結果のみ返されます。
-  rdkit.do_chiral_sss :
   部分構造マッチングで立体化学が使われているかいないかを切り替えます。（2013_03リリースから利用可能）
-  rdkit.sss_fp_size :
   部分構造検索スクリーニングで使われるフィンガープリントのサイズ（ビット数）。
-  rdkit.morgan_fp_size : Morganフィンガープリントのサイズ（ビット数）。
-  rdkit.featmorgan_fp_size :
   FeatMorganフィンガープリントのサイズ（ビット数）。
-  rdkit.layered_fp_size :
   層化（layered）フィンガープリントのサイズ（ビット数）
-  rdkit.rdkit_fp_size : RDKitフィンガープリントのサイズ（ビット数）
-  rdkit.torsion_fp_size :
   トポロジカルトーションビットベクトルフィンガープリントのサイズ（ビット数）
-  rdkit.atompair_fp_size :
   アトムペアビットベクトルフィンガープリントのサイズ（ビット数）
-  rdkit.avalon_fp_size : Avalonフィンガープリントのサイズ（ビット数）

演算子
======================================================
[`Operators <https://www.rdkit.org/docs/Cartridge.html#operators>`__]

類似度検索
-----------------------------------------------------
[`Similarity search <https://www.rdkit.org/docs/Cartridge.html#similarity-search>`__]

-  % :
   タニモト類似度を使った検索で使われる演算子。２つのフィンガープリント（sfp値２つ、あるいはbfp値２つ）の間のタニモト類似度が閾値
   rdkit.tanimoto_threshold を超えているか否かを返します。
-  # :
   Dice類似度を使った検索で使われる演算子。２つのフィンガープリント（sfp値２つ、あるいはbfp値２つ）の間のDice類似度が閾値
   rdkit.dice_threshold を超えているか否かを返します。
-  <%> :
   タニモト k近傍法検索に使われます（近傍の分子の並べ替えられたリストを返すために使われます）。
-  <#> :
   Dice k近傍法検索に使われます（近傍の分子の並べ替えられたリストを返すために使われます）。

部分構造検索と正確な構造の検索
-----------------------------------------------------
[`Substructure and exact strucure search <https://www.rdkit.org/docs/Cartridge.html#substructure-and-exact-structure-search>`__]

-  @> :
   部分構造検索の演算子。演算子の右側におかれたmolあるいはqmolが、左側のmolの部分構造であるか否かを返します。
-  <@ :
   部分構造検索の演算子。演算子の左側におかれたmolあるいはqmolが、右側のmolの部分構造であるか否かを返します。
-  @= : ２つの分子が同一であるか否かを返します。

分子の比較
-----------------------------------------------------
[`Molecule comparison <https://www.rdkit.org/docs/Cartridge.html#molecule-comparison>`__]

-  < : 左側のmolが右側もmolよりも小さいか否かを返します。
-  > : 左側のmolが右側もmolよりも大きいか否かを返します。
-  = : 左側のmolが右側のmolと等しいか否かを返します。
-  <= : 左側のmolが右側のmol以下であるか否かを返します。
-  >= : 左側のmolが右側のmol以上であるか否かを返します。

*注意*
２つの分子は次の比較を順番に行うことで比べられます。より後の順番の比較は、それより前の値が等しい時だけ行われます:

# Number of atoms # Number of bonds # Molecular weight # Number of rings

上記の全てが同じ値で、２つ目の分子が最初の分子の部分構造の場合、分子は等しいと宣言されます。さもなくば（そうなるべきではありませんが）最初の分子は、２番目よりも小さいと勝手に定義されます。

カートリッジには他にも演算子が定義されていますが、これらは内部の目的のためだけに使われます。

関数
======================================================
[`Functions <https://www.rdkit.org/docs/Cartridge.html#functions>`__]

フィンガープリント関連
-----------------------------------------------------
[`Fingerprint Related <https://www.rdkit.org/docs/Cartridge.html#fingerprint-related>`__]

フィンガープリントの生成
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Generating fingerprints <https://www.rdkit.org/docs/Cartridge.html#generating-fingerprints>`__]

-  morgan_fp(mol,int default 2) :
   結合関係の不変量を使って求めた、カウントベースのMorganフィンガープリントであるsfpを返します。２つ目の引数は半径を与えます。ECFPに類似したフィンガープリントです。
-  morganbv_fp(mol,int default 2) :
   結合関係の不変量を使って求めた、ビットベクトルのMorganフィンガープリントであるbfpを返します。２つ目の引数は半径を与えます。ECFPに類似したフィンガープリントです。
-  featmorgan_fp(mol,int default 2) :
   化学的特徴の不変量を使って求めた、カウントベースのMorganフィンガープリントであるsfpを返します。２つ目の引数は半径を与えます。FCFPに類似したフィンガープリントです。
-  featmorganbv_fp(mol,int default 2) :
   化学的特徴の不変量を使って求めた、ビットベクトルのMorganフィンガープリントであるbfpを返します。２つ目の引数は半径を与えます。FCFPに類似したフィンガープリントです。
-  rdkit_fp(mol) :
   RDKitフィンガープリントのbfpを返します。ハッシュ化された分子のサブグラフをつかったDaylightフィンガープリントです。.
-  atompair_fp(mol) :
   カウントベースのアトムペアフィンガープリントであるsfpを返します。
-  atompairbv_fp(mol) :
   ビットベクトルのアトムペアフィンガープリントであるbfpを返します。
-  torsion_fp(mol) :
   カウントベースのトポロジカルトーションフィンガープリントであるsfpを返します。
-  torsionbv_fp(mol) :
   ビットベクトルのトポロジカルトーションフィンガープリントであるbfpを返します。
-  layered_fp(mol) :
   層化（layered）フィンガープリントであるbfpを返します。ハッシュ化された分子のサブグラフを使った実験的な部分構造のフィンガープリントです。
-  maccs_fp(mol) :
   MACCSフィンガープリントであるbfpを返します（\ *2013_01リリースから利用可能*\ )。

フィンガープリントの取り扱い
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Working with fingerpritns <https://www.rdkit.org/docs/Cartridge.html#working-with-fingerprints>`__]

-  tanimoto_sml(fp,fp) :
   同じタイプの２つのフィンガープリント（sfp値2つかbfp値２つ）の間のタニモト類似度を返します。
-  dice_sml(fp,fp) :
   同じタイプの２つのフィンガープリント（sfp値2つかbfp値２つ）の間のDice類似度を返します。
-  size(bfp) : bfp（のビット数）の長さを返します。
-  add(sfp,sfp) :
   ２つのsfp引数の要素ごとの足し算によって作成される１つのsfpを返します。
-  subtract(sfp,sfp) :
   ２つのsfp引数の要素ごとの引き算によって作成される１つのsfpを返します。
-  all_values_lt(sfp,int) :
   sfpの引数の全ての要素が、整数の引数よりも小さいか否かを示すブール値を返します。
-  all_values_gt(sfp,int) :
   sfpの引数の全ての要素が、整数の引数よりも大きいか否かを示すブール値を返します。

フィンガープリントの入出力
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Fingerprint I/O <https://www.rdkit.org/docs/Cartridge.html#fingerprint-i-o>`__]

-  bfp_to_binary_text(bfp) :
   他のソフトウェアでRDKitフィンガープリントに変換し直すことができる、フィンガープリントのバイナリ列表現を保存したbytea型を返します。
   （Q3 2012（2012_09）リリースから利用可能）
-  bfp_from_binary_text(bytea) :
   フィンガープリントのバイナリ列表現からbfpを作成します。

分子関連
-----------------------------------------------------
[`Molecule Related <https://www.rdkit.org/docs/Cartridge.html#molecule-related>`__]

分子の入出力と妥当性の検証
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Molecule I/O and Validation <https://www.rdkit.org/docs/Cartridge.html#molecule-i-o-and-validation>`__]

-  is_valid_smiles(smiles) :
   SMILES文字列が妥当なRDKit分子を生成しているか否かを返します。
-  is_valid_ctab(ctab) :
   CTAB（mol block）文字列が妥当なRDKit分子を生成しているか否かを返します。
-  is_valid_smarts(smarts) :
   SMARTS文字列が妥当なRDKit分子を生成しているか否かを返します。
-  is_valid_mol_pkl(bytea) :
   バイナリ列（bytea型）がRDKit分子に変換しなおすことができるか否かを返します。（\ *Q3 2012（2012_09）リリースから利用可能*\ ）
-  mol_from_smiles(smiles) :
   SMILES文字列に対して分子を返し、分子の構築が失敗した場合はNULLを返します。
-  mol_from_smarts(smarts) :
   SMARTS文字列に対して分子を返し、分子の構築が失敗した場合はNULLを返します。
-  mol_from_ctab(ctab, bool default false) :
   CTAB（mol block）文字列に対して分子を返し、分子の構築が失敗した場合はNULLを返します。オプションの２つ目の引数は分子の座標を保存するか否かをコントロールします。
-  mol_from_pkl(bytea) :
   バイナリ列（bytea型）に対して分子を返し、分子の構築が失敗した場合はNULLを返します。（\ *Q3 2012（2012_09）リリースから利用可能*\ ）
-  qmol_from_smiles(smiles) :
   SMILES文字列にたいしてクエリ分子を返し、分子の構築が失敗した場合NULLを返します。
   SMILESの中で明示的に表されている水素原子は、水素原子が結合している原子のクエリフィーチャーに変換されます。
-  qmol_from_ctab(ctab, bool default false) :
   CTAB（mol block）文字列にたいしてクエリ分子を返し、分子の構築が失敗した場合NULLを返します。
   SMILESの中で明示的に表されている水素原子は、水素原子が結合している原子のクエリフィーチャーに変換されます。
   オプションの２つ目の引数は分子の座標を保存するか否かをコントロールします。
-  mol_to_smiles(mol) : 分子のカノニカルSMILESを返します。
-  mol_to_smarts(mol) : 分子のSMARTS文字列を返します。
-  mol_to_pkl(mol) : 分子のバイナリ列（bytea型）を返します。（\ *Q3 2012（2012_09）リリースから利用可能*\ ）
-  mol_to_ctab(mol,bool default true) :
   分子のCTAB（mol block）文字列を返します。オプションの２つ目の引数は、座標を持っていない分子に対して２D座標を生成するか否かをコントロールします。
-  mol_to_svg(mol,string default ‘’,int default 250, int default 200, string default ‘’) :
   分子の描画を有するSVGを返します。オプションのパラメータは、凡例として使われる文字列、画像の幅、高さ、そして追加のレンダリングについてのパラメータを有するJSONです。（\ *2016_09リリースから利用可能*\ ）

部分構造操作
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Substructure operations <https://www.rdkit.org/docs/Cartridge.html#substructure-operations>`__]

-  substruct(mol,mol) :
   ２つ目のmolが１つ目の部分構造であるか否かを返します。
-  substruct_count(mol,mol,bool default true) :
   ２つ目の分子と１つ目の分子の間でマッチしている部分構造の数を返します。３つ目の引数は部分構造のマッチングを重複を除いたユニークなものとするか否かを切り替えます。（\ *2013_03リリースから利用可能*\ ）
-  mol_adjust_query_properties(mol,string default ‘’) :
   追加のクエリ情報が付加された新しい分子を返します。（2016_09リリースから利用可能）

記述子
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Descriptors <https://www.rdkit.org/docs/Cartridge.html#descriptors>`__]

-  mol_amw(mol) : 分子のAMWを返します。
-  mol_logp(mol) : 分子のMOlLogPを返します。
-  mol_tpsa(mol) : 分子のトポロジカル極性表面積を返します（\ *Q1 2011（2011_03）リリースから利用可能*\ ）。
-  mol_fractioncsp3(mol) :
   sp3混成の炭素の割合を返します（\ *2013_03リリースから利用可能*\ ）。
-  mol_hba(mol) :
   分子のリピンスキーの水素結合アクセプターの数（すなわち、OとNの数）を返します。
-  mol_hbd(mol) :
   分子のリピンスキーの水素結合ドナーの数（すなわち、少なくとも１つHをもつOとNの数）を返します。
-  mol_numatoms(mol) : 分子に含まれる原子の総数を返します。
-  mol_numheavyatoms(mol) : 分子に含まれるヘビーアトムの数を返します。
-  mol_numrotatablebonds(mol) :
   分子に含まれる回転可能な結合の数を返します（\ *Q1 2011（2011_03）リリースから利用可能*\ ）。
-  mol_numheteroatoms(mol) : 分子に含まれるヘテロ原子の数を返します（\ *Q1 2011（2011_03）リリースから利用可能*\ ）。
-  mol_numrings(mol) : 分子に含まれる環構造の数を返します（\ *Q1 2011（2011_03）リリースから利用可能*\ ）。
-  mol_numaromaticrings(mol) :
   分子に含まれる芳香環の数を返します（\ *2013_03リリースから利用可能*\ ）。
-  mol_numaliphaticrings(mol) :
   分子に含まれる脂肪族（少なくとも一つ非芳香族結合をもつ）環構造の数を返します（\ *2013_03リリースから利用可能*\ ）。
-  mol_numsaturatedrings(mol) :
   分子に含まれる飽和環の数を返します（\ *2013_03リリースから利用可能*\ ）。
-  mol_numaromaticheterocycles(mol) :
   分子に含まれる芳香族ヘテロ環の数を返します（\ *2013_03リリースから利用可能*\ ）。
-  mol_numaliphaticheterocycles(mol) :
   分子に含まれる脂肪族（少なくとも一つ非芳香族結合をもつ）環構造の数を返します（\ *2013_03リリースから利用可能*\ ）。
-  mol_numsaturatedheterocycles(mol) :
   分子に含まれる飽和ヘテロ環の数を返します（\ *2013_03リリースから利用可能*\ ）。
-  mol_numaromaticcarbocycles(mol) :
   分子に含まれる芳香族炭素環の数を返します（\ *2013_03リリースから利用可能*\ ）。
-  mol_numaliphaticcarbocycles(mol) :
   分子に含まれる脂肪族（少なくとも一つ非芳香族結合をもつ）炭素環の数を返します（\ *2013_03リリースから利用可能*\ ）。
-  mol_numsaturatedcarbocycles(mol) :
   分子に含まれる飽和炭素環の数を返します（\ *2013_03リリースから利用可能*\ ）。
-  mol_inchi(mol) :
   分子のInChIを返します。（\ *2011_06リリースから利用可能で、RDKitがInChIサポートとともにビルドされている必要があります*\ ）。
-  mol_inchikey(mol) :
   分子のInChI keyを返します。（\ *2011_06リリースから利用可能で、RDKitがInChIサポートとともにビルドされている必要があります*\ ）。
-  mol_formula(mol,bool default false, bool default true) :
   分子式をもつ文字列を返します。２つ目の引数は分子式に同位体の情報を含むかどうかをコントロールします。
   ３つ目の引数は“D”と“T”を、[2H]と[3H]の代わりに使うかどうかをコントロールします。（\ *2014_03リリースから利用可能*\ ）

結合関係記述子
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Connectivity Desctiptors <https://www.rdkit.org/docs/Cartridge.html#connectivity-descriptors>`__]

-  mol_chi0v(mol) - mol_chi4v(mol) :
   x=0-4のChiXvの値を返します（\ *2012_01リリースから利用可能*\ ）
-  mol_chi0n(mol) - mol_chi4n(mol) :
   x=0-4のChiXnの値を返します（\ *2012_01リリースから利用可能*\ ）
-  mol_kappa1(mol) - mol_kappa3(mol) :
   x=1-3のkappaXの値を返します（\ *2012_01リリースから利用可能*\ ）
-  mol_numspiroatoms :
   分子に含まれるスピロ原子の数を返します（\ *2015_09リリースから利用可能*\ ）
-  mol_numbridgeheadatoms :
   分子に含まれる橋頭位の原子の数を返します（\ *2015_09リリースから利用可能*\ ）。

MCS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`MCS <https://www.rdkit.org/docs/Cartridge.html#mcs>`__]

-  fmcs(mols) : 分子のセットについてMCSを計算する集合関数
-  fmcs_smiles(text, json default ‘’) :
   スペースで区切られたSMILESのセットについてMCSを計算します。オプションのjson引数はMCSコードにパラメータを与えるために使われます。

その他
-----------------------------------------------------
[`Other <https://www.rdkit.org/docs/Cartridge.html#other>`__]

-  rdkit_version() :
   カートリッジのバージョン番号をを持つ文字列を返します。

カートリッジには他にも演算子が定義されていますが、これらは内部の目的のためだけに使われます。

Pythonからカートリッジを使う方法
****************************************************
[`Using the Cartridge from Python <https://www.rdkit.org/docs/Cartridge.html#using-the-cartridge-from-python>`__]

postgresqlに接続するためのお勧めのアダブタはpyscopg2です（https://pypi.python.org/pypi/psycopg2 ）。

我々のChEMBLのローカルコピーに接続し、基本的な部分構造検索を行った場合の例を示します：

.. code:: python

   >>> import psycopg2
   >>> conn = psycopg2.connect(database='chembl_16')
   >>> curs = conn.cursor()
   >>> curs.execute('select * from rdk.mols where m@>%s',('c1cccc2c1nncc2',))
   >>> curs.fetchone()
   (9830, 'CC(C)Sc1ccc(CC2CCN(C3CCN(C(=O)c4cnnc5ccccc54)CC3)CC2)cc1')

各分子に対してSMILESを返します。分子を取得した後で、さらに操作を行いたいなら、postgresqlにpickel形式で分子を返すよう要求するのがより効率的です。

.. code:: python

   >>> curs.execute('select molregno,mol_send(m) from rdk.mols where m@>%s',('c1cccc2c1nncc2',))
   >>> row = curs.fetchone()
   >>> row
   (9830, <read-only buffer for 0x...>)

これらのpickelは分子に変換することができます。

.. code:: python

   >>> from rdkit import Chem
   >>> m = Chem.Mol(str(row[1]))
   >>> Chem.MolToSmiles(m,True)
   'CC(C)Sc1ccc(CC2CCN(C3CCN(C(=O)c4cnnc5ccccc54)CC3)CC2)cc1'

ライセンス
**********************
[`License <https://www.rdkit.org/docs/Cartridge.html#license>`__]

この文書の著作権は copyright (C) 2013-2018 by Greg Landrum に所属しています。

この文書はCreative Commons Attribution-ShareAlike 4.0 Licenseのもとでライセンスされています。
このライセンスを見るためにはhttp://creativecommons.org/licenses/by-sa/4.0/ にアクセスするか、
Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.　に手紙を送ってください。

このライセンスの意図はRDKitそのものの意図と似ています。
簡単に言えば“これを使ってなんでもやりたいことをやっていいですが、私たちの功績にも言及してください”
