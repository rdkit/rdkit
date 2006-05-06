import RDConfig


molFuncCreation="""
DROP FUNCTION rd_libversion();
CREATE FUNCTION rd_libversion() RETURNS text AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
DROP FUNCTION rd_canonsmiles(text);
CREATE FUNCTION rd_canonsmiles(text) RETURNS text AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
DROP FUNCTION rd_molpickle(text);
CREATE FUNCTION rd_molpickle(text) RETURNS bytea AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
DROP FUNCTION rd_hassubstruct(text,text);
CREATE FUNCTION rd_hassubstruct(text,text) RETURNS bool AS '%(base)s/RDLib','rd_hassubstruct_smi' LANGUAGE 'c' immutable;
DROP FUNCTION rd_hassubstruct(bytea,bytea);
CREATE FUNCTION rd_hassubstruct(bytea,bytea) RETURNS bool AS '%(base)s/RDLib','rd_hassubstruct_pkl' LANGUAGE 'c' immutable;
DROP FUNCTION rd_substructcount(text,text);
CREATE FUNCTION rd_substructcount(text,text) RETURNS int AS '%(base)s/RDLib', 'rd_substructcount_smi' LANGUAGE 'c' immutable;
DROP FUNCTION rd_substructcount(bytea,bytea);
CREATE FUNCTION rd_substructcount(bytea,bytea) RETURNS int AS '%(base)s/RDLib','rd_substructcount_pkl' LANGUAGE 'c' immutable;
DROP FUNCTION rd_similarityfp(text);
CREATE FUNCTION rd_similarityfp(text) RETURNS bytea AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
DROP FUNCTION rd_similarityfp_bits(text);
CREATE FUNCTION rd_similarityfp_bits(text) RETURNS bit AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
DROP FUNCTION rd_substructfp(text);
CREATE FUNCTION rd_substructfp(text) RETURNS bytea AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
DROP FUNCTION rd_substructfp_bits(text);
CREATE FUNCTION rd_substructfp_bits(text) RETURNS bit AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
DROP FUNCTION rd_mollogp(text);
CREATE FUNCTION rd_mollogp(text) RETURNS double precision AS '%(base)s/RDLib', 'rd_mollogp_smi' LANGUAGE 'c' immutable;
DROP FUNCTION rd_mollogp(bytea);
CREATE FUNCTION rd_mollogp(bytea) RETURNS double precision AS '%(base)s/RDLib','rd_mollogp_pkl' LANGUAGE 'c' immutable;
DROP FUNCTION rd_amw(text);
CREATE FUNCTION rd_amw(text) RETURNS double precision AS '%(base)s/RDLib', 'rd_amw_smi' LANGUAGE 'c' immutable;
DROP FUNCTION rd_amw(bytea);
CREATE FUNCTION rd_amw(bytea) RETURNS double precision AS '%(base)s/RDLib','rd_amw_pkl' LANGUAGE 'c' immutable;
"""

sigFuncCreation="""
DROP FUNCTION rd_fpsize(bytea);
CREATE FUNCTION rd_fpsize(bytea) RETURNS int AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
DROP FUNCTION rd_allprobebitsmatch(bytea,bytea);
CREATE FUNCTION rd_allprobebitsmatch(bytea,bytea) RETURNS bool AS '%(base)s/RDLib','rd_allprobebitsmatch_pkl' LANGUAGE 'c' immutable;
DROP FUNCTION rd_allprobebitsmatch(bit,bit);
CREATE FUNCTION rd_allprobebitsmatch(bit,bit) RETURNS bool AS '%(base)s/RDLib','rd_allprobebitsmatch_bits' LANGUAGE 'c' immutable;
DROP FUNCTION rd_tanimoto(bytea,bytea);
CREATE FUNCTION rd_tanimoto(bytea,bytea) RETURNS double precision AS '%(base)s/RDLib','rd_tanimoto_pkl' LANGUAGE 'c' immutable;
DROP FUNCTION rd_tanimoto(bit,bit);
CREATE FUNCTION rd_tanimoto(bit,bit) RETURNS double precision AS '%(base)s/RDLib','rd_tanimoto_bits' LANGUAGE 'c' immutable;
DROP FUNCTION rd_dice(bytea,bytea);
CREATE FUNCTION rd_dice(bytea,bytea) RETURNS double precision AS '%(base)s/RDLib','rd_dice_pkl' LANGUAGE 'c' immutable;
DROP FUNCTION rd_dice(bit,bit);
CREATE FUNCTION rd_dice(bit,bit) RETURNS double precision AS '%(base)s/RDLib','rd_dice_bits' LANGUAGE 'c' immutable;
DROP FUNCTION rd_cosine(bytea,bytea);
CREATE FUNCTION rd_cosine(bytea,bytea) RETURNS double precision AS '%(base)s/RDLib','rd_cosine_pkl' LANGUAGE 'c' immutable;
DROP FUNCTION rd_cosine(bit,bit);
CREATE FUNCTION rd_cosine(bit,bit) RETURNS double precision AS '%(base)s/RDLib','rd_cosine_bits' LANGUAGE 'c' immutable;
"""

molAddTrigger="""
CREATE OR REPLACE FUNCTION molAddTrigger()
  RETURNS "trigger" AS
$BODY$
  if TD['event'] != 'INSERT' or TD['when'] != 'BEFORE' or TD['level'] != 'ROW':
    return 'OK'
  id = TD['new']['id']
  smi = TD['new']['smiles']

  if not TD['new']['pickle']:
    if SD.has_key('mol_pickle_plan'):
      plan = SD['mol_pickle_plan']
    else:
      plan = plpy.prepare('select rd_molpickle($1) as pickle',('text',))
      SD['mol_pickle_plan']=plan
    pkl = plpy.execute(plan,(smi,))[0]['pickle']
    if not pkl:
      plpy.error('Could not pickle input smiles. operation aborted')
    TD['new']['pickle']=pkl
      
    res = 'MODIFY'
  else:
    res = 'OK'

  if SD.has_key('mol_insert_plan'):
    plan = SD['mol_insert_plan']
  else:
    plan = plpy.prepare('insert into %(fpTable)s ($1,rd_substructfp_bits($2),rd_similarityfp_bits($2))',
                        ('text','text'))
    SD['mol_insert_plan']=plan
  plpy.execute(plan,(id,smi))

  return res
$BODY$
  LANGUAGE 'plpythonu' VOLATILE;
"""

installTrigger="""
CREATE TRIGGER %(tblName)s_after
  AFTER INSERT OR UPDATE OR DELETE
  ON %(tblName)s
  FOR EACH ROW
  EXECUTE PROCEDURE molAddTrigger();
"""

if __name__=='__main__':
  base = RDConfig.RDBaseDir+'/Code/PgSQL/RDLib'
  print sigFuncCreation%locals()
  print molFuncCreation%locals()
