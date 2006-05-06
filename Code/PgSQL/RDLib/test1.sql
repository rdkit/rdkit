
CREATE FUNCTION rd_hassubstruct(text,text) RETURNS bool
--    AS 'c:/glandrum/RD/trunk/Code/Demos/PgSQL/RDLib/RDLib' LANGUAGE 'c';
    AS '/home2/glandrum/RD/trunk/Code/Demos/PgSQL/RDLib/RDLib' LANGUAGE 'c';
CREATE FUNCTION rd_hassubstructpkl(bytea,bytea) RETURNS bool
--    AS 'c:/glandrum/RD/trunk/Code/Demos/PgSQL/RDLib/RDLib' LANGUAGE 'c';
    AS '/home2/glandrum/RD/trunk/Code/Demos/PgSQL/RDLib/RDLib' LANGUAGE 'c';
CREATE FUNCTION rd_substructfp(text) RETURNS bytea
--    AS 'c:/glandrum/RD/trunk/Code/Demos/PgSQL/RDLib/RDLib' LANGUAGE 'c';
    AS '/home2/glandrum/RD/trunk/Code/Demos/PgSQL/RDLib/RDLib' LANGUAGE 'c';
CREATE FUNCTION rd_fpsize(bytea) RETURNS int
--    AS 'c:/glandrum/RD/trunk/Code/Demos/PgSQL/RDLib/RDLib' LANGUAGE 'c';
    AS '/home2/glandrum/RD/trunk/Code/Demos/PgSQL/RDLib/RDLib' LANGUAGE 'c';
CREATE FUNCTION rd_allprobebitsmatch(bytea,bytea) RETURNS bool
--    AS 'c:/glandrum/RD/trunk/Code/Demos/PgSQL/RDLib/RDLib' LANGUAGE 'c';
    AS '/home2/glandrum/RD/trunk/Code/Demos/PgSQL/RDLib/RDLib' LANGUAGE 'c';

SELECT rd_hassubstruct('cC','c1ccccc1C');
drop table temp;
create table temp (id text,smiles varchar(256));
insert into temp values ('mol-1','CCC');
insert into temp values ('mol-2','Cc1ncccc1');
select * from temp where rd_hassubstruct('Cc1ncccc1',smiles)=1;

select rd_substructfp('c1ccccc1');
select rd_substructfp('');
select rd_fpsize(rd_substructfp(''));
select rd_allprobebitsmatch(rd_substructfp('c1ccccc1'),rd_substructfp('c1ccccc1'));
select rd_allprobebitsmatch(rd_substructfp('c1ccccc1'),rd_substructfp('c1ccccc1C'));
select rd_allprobebitsmatch(rd_substructfp('c1ccccc1C'),rd_substructfp('c1ccccc1'));

drop table testfps;
create table testfps (id text,smiles varchar(256),fp bytea);
insert into testfps values ('mol-1','CC',rd_substructfp('CC'));
insert into testfps values ('mol-2','Cc1ccccc1',rd_substructfp('Cc1ccccc1'));
insert into testfps values ('mol-3','Cc1ccncc1',rd_substructfp('Cc1ccncc1'));
select id,smiles from testfps where rd_allprobebitsmatch(rd_substructfp('c1ccccc1'),fp);


DROP FUNCTION rd_hassubstruct(text,text);
DROP FUNCTION rd_substructfp(text);
DROP FUNCTION rd_fpsize(bytea);
DROP FUNCTION rd_allprobebitsmatch(bytea,bytea);

