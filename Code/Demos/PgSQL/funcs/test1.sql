---------------------------------------------------------------------------
--
-- funcs.sql-
--	  Tutorial on using functions in POSTGRES.
--
--
-- Copyright (c) 1994-5, Regents of the University of California
--
-- $Id: funcs.source,v 1.5 2001/10/26 20:45:33 tgl Exp $
--
---------------------------------------------------------------------------

-----------------------------
-- Creating C Functions
--	in addition to SQL functions, you can also create C functions. 
--	See funcs.c for the definition of the C functions.
-----------------------------

CREATE FUNCTION substruct(text,text) RETURNS int
    AS '/home2/glandrum/RD/trunk/Code/Demos/PgSQL/funcs/funcs' LANGUAGE 'c';
CREATE FUNCTION substructfp(text) RETURNS bytea
    AS '/home2/glandrum/RD/trunk/Code/Demos/PgSQL/funcs/funcs' LANGUAGE 'c';
CREATE FUNCTION fpsize(bytea) RETURNS int
    AS '/home2/glandrum/RD/trunk/Code/Demos/PgSQL/funcs/funcs' LANGUAGE 'c';

SELECT substruct('c1ccccc1C','cC');
drop table temp;
create table temp (id text,smiles varchar(256));
insert into temp values ('mol-1','CCC');
insert into temp values ('mol-2','Cc1ncccc1');
select * from temp where substruct(smiles,'Cc1ncccc1')=1;

-- select * from d1 where substruct(smiles,'Cc1ncccc1')=1;

select substructfp('c1ccccc1');
select substructfp('');
select fpsize(substructfp(''));


DROP FUNCTION substruct(text,text);
DROP FUNCTION substructfp(text);

