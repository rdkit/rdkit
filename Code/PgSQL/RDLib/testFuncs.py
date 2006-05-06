# $Id$
#
#  Copyright (C) 2005  greg Landrum and Rational Discovery LLC
#   All Rights Reserved
#
import RDConfig
import unittest,cPickle,os,sys
from Dbase.DbConnection import DbConnect
from Dbase import DbModule
import Chem

def feq(v1,v2,tol=1e-4):
  return abs(v1-v2)<tol

class TestCase(unittest.TestCase):
  def setUp(self):
    self.conn=DbConnect('::RDTests')

  def registerFuncs(self):
    curs = self.conn.GetCursor()
    cmd="""
DROP FUNCTION rd_libversion();
DROP FUNCTION rd_canonsmiles(text);
DROP FUNCTION rd_molpickle(text);
DROP FUNCTION rd_hassubstruct(text,text);
DROP FUNCTION rd_hassubstruct(bytea,bytea);
DROP FUNCTION rd_substructcount(text,text);
DROP FUNCTION rd_substructcount(bytea,bytea);
DROP FUNCTION rd_fpsize(bytea);
DROP FUNCTION rd_allprobebitsmatch(bytea,bytea);
DROP FUNCTION rd_allprobebitsmatch(bit,bit);
DROP FUNCTION rd_tanimoto(bytea,bytea);
DROP FUNCTION rd_tanimoto(bit,bit);
DROP FUNCTION rd_dice(bytea,bytea);
DROP FUNCTION rd_dice(bit,bit);
DROP FUNCTION rd_cosine(bytea,bytea);
DROP FUNCTION rd_cosine(bit,bit);
DROP FUNCTION rd_similarityfp(text);
DROP FUNCTION rd_similarityfp_bits(text);
DROP FUNCTION rd_substructfp(text);
DROP FUNCTION rd_substructfp_bits(text);
DROP FUNCTION rd_mollogp(text);
DROP FUNCTION rd_mollogp(bytea);
DROP FUNCTION rd_molamw(text);
DROP FUNCTION rd_molamw(bytea);
    """
    for line in cmd.split('\n'):
      if line:
        #print line
        try:
          curs.execute(line)
        except:
          pass

    self.conn.Commit()
    base = RDConfig.RDBaseDir+'/Code/PgSQL/RDLib'
    cmd="""
CREATE FUNCTION rd_libversion() RETURNS text AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_canonsmiles(text) RETURNS text AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_molpickle(text) RETURNS bytea AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_hassubstruct(text,text) RETURNS bool AS '%(base)s/RDLib','rd_hassubstruct_smi' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_hassubstruct(bytea,bytea) RETURNS bool AS '%(base)s/RDLib','rd_hassubstruct_pkl' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_substructcount(text,text) RETURNS int AS '%(base)s/RDLib', 'rd_substructcount_smi' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_substructcount(bytea,bytea) RETURNS int AS '%(base)s/RDLib','rd_substructcount_pkl' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_fpsize(bytea) RETURNS int AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_allprobebitsmatch(bytea,bytea) RETURNS bool AS '%(base)s/RDLib','rd_allprobebitsmatch_pkl' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_allprobebitsmatch(bit,bit) RETURNS bool AS '%(base)s/RDLib','rd_allprobebitsmatch_bits' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_tanimoto(bytea,bytea) RETURNS double precision AS '%(base)s/RDLib','rd_tanimoto_pkl' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_dice(bytea,bytea) RETURNS double precision AS '%(base)s/RDLib','rd_dice_pkl' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_cosine(bytea,bytea) RETURNS double precision AS '%(base)s/RDLib','rd_cosine_pkl' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_tanimoto(bit,bit) RETURNS double precision AS '%(base)s/RDLib','rd_tanimoto_bits' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_dice(bit,bit) RETURNS double precision AS '%(base)s/RDLib','rd_dice_bits' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_cosine(bit,bit) RETURNS double precision AS '%(base)s/RDLib','rd_cosine_bits' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_similarityfp(text) RETURNS bytea AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_similarityfp_bits(text) RETURNS bit AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_substructfp(text) RETURNS bytea AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_substructfp_bits(text) RETURNS bit AS '%(base)s/RDLib' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_mollogp(text) RETURNS double precision AS '%(base)s/RDLib', 'rd_mollogp_smi' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_mollogp(bytea) RETURNS double precision AS '%(base)s/RDLib','rd_mollogp_pkl' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_molamw(text) RETURNS double precision AS '%(base)s/RDLib', 'rd_amw_smi' LANGUAGE 'c' immutable;
CREATE FUNCTION rd_molamw(bytea) RETURNS double precision AS '%(base)s/RDLib','rd_amw_pkl' LANGUAGE 'c' immutable;
    """%locals()
    for line in cmd.split('\n'):
      if line:
        try:
          curs.execute(line)
        except:
          #pass
          import traceback
          traceback.print_exc()


    self.conn.Commit()
  def test0(self):
    " registering functions and basic substructs "
    self.registerFuncs()
    curs = self.conn.GetCursor()

    curs.execute("SELECT rd_libversion();")
    res = curs.fetchone()[0]
    rev = '$Rev$'
    self.failUnless(res.split(' ')[0]==rev.split(' ')[0])

    curs.execute("SELECT rd_hassubstruct('cC','c1ccccc1C');")
    res = curs.fetchone()
    self.failUnless(res[0])




  def test1(self):
    " basic fingerprint operations "
    curs = self.conn.GetCursor()
    curs.execute("SELECT rd_fpsize(rd_substructfp('c1ccncc1'))")
    res = curs.fetchone()
    self.failUnless(res[0]==2048)


  def test2(self):
    " bit comparisons "
    curs = self.conn.GetCursor()
    curs.execute("""
    SELECT rd_allprobebitsmatch(rd_substructfp('c1ccccn1'),
                                rd_substructfp('c1cccnc1C'));
    """)
    res = curs.fetchone()
    self.failUnless(res[0])
    curs.execute("""
    SELECT rd_allprobebitsmatch(rd_substructfp('c1cccnc1'),
                                rd_substructfp('c1ccccc1C'));
    """)
    res = curs.fetchone()
    self.failIf(res[0])


  def test3(self):
    " error handling with molecules "
    curs = self.conn.GetCursor()
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("SELECT rd_hassubstruct('QcC','c1ccccc1C');"))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("SELECT rd_hassubstruct('cC','Qc1ccccc1C');"))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("SELECT rd_hassubstruct('nC','c1ccccn1C');"))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("""
    SELECT rd_allprobebitsmatch(rd_substructfp('c1ccccn1C'),
                                rd_substructfp('c1ccccc1C'));
    """))


    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("""
    SELECT rd_allprobebitsmatch(rd_substructfp('c1ccccn1'),
                                rd_substructfp('c1ccccn1C'));
    """))



  def test4(self):
    " Similarity metrics "
    curs = self.conn.GetCursor()
    curs.execute("""
    SELECT rd_tanimoto(rd_substructfp('c1ccncc1'),
                       rd_substructfp('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(feq(res[0],1.0))

    curs.execute("""
    SELECT rd_tanimoto(rd_substructfp('c1ccncc1C'),
                       rd_substructfp('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(res[0]<1.0)

    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("""
    SELECT rd_tanimoto(rd_substructfp('c1ccncn1C'),
                       rd_substructfp('c1ccncc1'));
                       """))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("""
    SELECT rd_tanimoto(rd_substructfp('c1ccncn1'),
                       rd_substructfp('c1ccncn1C'));
                       """))

    curs.execute("""
    SELECT rd_dice(rd_substructfp('c1ccncc1'),
                       rd_substructfp('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(feq(res[0],1.0))

    curs.execute("""
    SELECT rd_dice(rd_substructfp('c1ccncc1C'),
                       rd_substructfp('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(res[0]<1.0)

    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("""
    SELECT rd_dice(rd_substructfp('c1ccncn1C'),
                       rd_substructfp('c1ccncc1'));
                       """))

    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("""
    SELECT rd_dice(rd_substructfp('c1ccncn1'),
                       rd_substructfp('c1ccncn1C'));
                       """))

    curs.execute("""
    SELECT rd_cosine(rd_substructfp('c1ccncc1'),
                       rd_substructfp('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(feq(res[0],1.0))

    curs.execute("""
    SELECT rd_cosine(rd_substructfp('c1ccncc1C'),
                       rd_substructfp('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(res[0]<1.0)

    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("""
    SELECT rd_cosine(rd_substructfp('c1ccncn1C'),
                       rd_substructfp('c1ccncc1'));
                       """))

    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("""
    SELECT rd_cosine(rd_substructfp('c1ccncn1'),
                       rd_substructfp('c1ccncn1C'));
                       """))

  def test5(self):
    " Similarity metrics "
    curs = self.conn.GetCursor()
    curs.execute("""
    SELECT rd_tanimoto(B'11110000',B'11110000');
                       """)
    self.failUnless(feq(curs.fetchone()[0],1.0))
    curs.execute("""
    SELECT rd_tanimoto(B'11110000',B'00001111');
                       """)
    self.failUnless(feq(curs.fetchone()[0],0.0))
    curs.execute("""
    SELECT rd_tanimoto(B'1111',B'1100');
                       """)
    self.failUnless(feq(curs.fetchone()[0],0.5))

    curs.execute("""
    SELECT rd_cosine(B'11110000',B'11110000');
                       """)
    self.failUnless(feq(curs.fetchone()[0],1.0))
    curs.execute("""
    SELECT rd_cosine(B'11110000',B'00001111');
                       """)
    self.failUnless(feq(curs.fetchone()[0],0.0))
    curs.execute("""
    SELECT rd_cosine(B'1111',B'1100');
                       """)
    self.failUnless(feq(curs.fetchone()[0],0.7071))

    curs.execute("""
    SELECT rd_dice(B'11110000',B'11110000');
                       """)
    self.failUnless(feq(curs.fetchone()[0],1.0))
    curs.execute("""
    SELECT rd_dice(B'11110000',B'00001111');
                       """)
    self.failUnless(feq(curs.fetchone()[0],0.0))
    curs.execute("""
    SELECT rd_dice(B'1111',B'1100');
                       """)
    self.failUnless(feq(curs.fetchone()[0],0.66667))

    curs.execute("""
    SELECT rd_tanimoto(B'11110000',B'0000000011110000');
                       """)
    self.failUnless(feq(curs.fetchone()[0],1.0))
    curs.execute("""
    SELECT rd_tanimoto(B'11110000',B'0000000000001111');
                       """)
    self.failUnless(feq(curs.fetchone()[0],0.0))
    curs.execute("""
    SELECT rd_tanimoto(B'11110000',B'00000000110000000');
                       """)
    self.failUnless(feq(curs.fetchone()[0],0.5))

    curs.execute("""
    SELECT rd_tanimoto(B'0000000011110000',B'11110000');
                       """)
    self.failUnless(feq(curs.fetchone()[0],1.0))
    curs.execute("""
    SELECT rd_tanimoto(B'0000000000001111',B'11110000');
                       """)
    self.failUnless(feq(curs.fetchone()[0],0.0))
    curs.execute("""
    SELECT rd_tanimoto(B'00000000110000000',B'11110000');
                       """)
    self.failUnless(feq(curs.fetchone()[0],0.5))

  def test6(self):
    " Similarity metrics using bit strings"
    curs = self.conn.GetCursor()
    curs.execute("""
    SELECT rd_tanimoto(rd_substructfp_bits('c1ccncc1'),
                       rd_substructfp_bits('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(feq(res[0],1.0))

    curs.execute("""
    SELECT rd_tanimoto(rd_substructfp_bits('c1ccncc1C'),
                       rd_substructfp_bits('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(res[0]<1.0)

    curs.execute("""
    SELECT rd_tanimoto(rd_substructfp_bits('c1ccncn1C'),
                       rd_substructfp_bits('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(feq(res[0],0.0))

    curs.execute("""
    SELECT rd_tanimoto(rd_substructfp_bits('c1ccncn1'),
                       rd_substructfp_bits('c1ccncn1C'));
                       """)
    res = curs.fetchone()
    self.failUnless(feq(res[0],0.0))


    curs.execute("""
    SELECT rd_dice(rd_substructfp_bits('c1ccncc1'),
                       rd_substructfp_bits('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(feq(res[0],1.0))

    curs.execute("""
    SELECT rd_dice(rd_substructfp_bits('c1ccncc1C'),
                       rd_substructfp_bits('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(res[0]<1.0)

    curs.execute("""
    SELECT rd_dice(rd_substructfp_bits('c1ccncn1C'),
                       rd_substructfp_bits('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(feq(res[0],0.0))

    curs.execute("""
    SELECT rd_dice(rd_substructfp_bits('c1ccncn1'),
                       rd_substructfp_bits('c1ccncn1C'));
                       """)
    res = curs.fetchone()
    self.failUnless(feq(res[0],0.0))

    curs.execute("""
    SELECT rd_cosine(rd_substructfp_bits('c1ccncc1'),
                       rd_substructfp_bits('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(feq(res[0],1.0))

    curs.execute("""
    SELECT rd_cosine(rd_substructfp_bits('c1ccncc1C'),
                       rd_substructfp_bits('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(res[0]<1.0)

    curs.execute("""
    SELECT rd_cosine(rd_substructfp_bits('c1ccncn1C'),
                       rd_substructfp_bits('c1ccncc1'));
                       """)
    res = curs.fetchone()
    self.failUnless(feq(res[0],0.0))


    curs.execute("""
    SELECT rd_cosine(rd_substructfp_bits('c1ccncn1'),
                       rd_substructfp_bits('c1ccncn1C'));
                       """)
    res = curs.fetchone()
    self.failUnless(feq(res[0],0.0))

  def test7(self):
    " similarity fingerprints "
    curs = self.conn.GetCursor()
    curs.execute("""
    SELECT rd_similarityfp_bits('C1CCC1')
                       """)
    res = curs.fetchone()
    self.failUnless(len(res[0])>=8)

    curs.execute("""
    SELECT rd_tanimoto(rd_similarityfp_bits('C1CCC1'),rd_similarityfp_bits('C1CCC1C'))
                       """)
    res = curs.fetchone()
    self.failUnless(res[0]>=0.0)
    self.failUnless(res[0]<1.0)

    curs.execute("""
    SELECT rd_tanimoto(rd_similarityfp('C1CCC1'),rd_similarityfp('C1CCC1C'))
                       """)
    res = curs.fetchone()
    self.failUnless(res[0]>=0.0)
    self.failUnless(res[0]<1.0)

  def test8(self):
    " Canonical smiles: "
    curs = self.conn.GetCursor()
    curs.execute("SELECT rd_canonsmiles('C1COCCC1');")
    ref = curs.fetchone()[0]
    self.failUnless(ref)
    curs.execute("SELECT rd_canonsmiles('C1CCOCC1');")
    smi = curs.fetchone()[0]
    self.failUnless(ref==smi)
    curs.execute("SELECT rd_canonsmiles('C1CCCCO1');")
    smi = curs.fetchone()[0]
    self.failUnless(ref==smi)

    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("SELECT rd_canonsmiles('Fail');"))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("SELECT rd_canonsmiles('');"))


  def test9(self):
    " substructure counts "
    curs = self.conn.GetCursor()

    curs.execute("SELECT rd_substructcount('O','OCCC(=O)O')")
    self.failUnless(curs.fetchone()[0]==3)
    curs.execute("SELECT rd_substructcount('N','OCCC(=O)O')")
    self.failUnless(curs.fetchone()[0]==0)
    curs.execute("SELECT rd_substructcount('[O,S]','SCCC(=O)O')")
    self.failUnless(curs.fetchone()[0]==3)
    curs.execute("SELECT rd_substructcount('[O,S]','OCCC(=O)O')")
    self.failUnless(curs.fetchone()[0]==3)
    
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("SELECT rd_substructcount('QcC','c1ccccc1C');"))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("SELECT rd_substructcount('QcC','c1ccccc1C');"))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("SELECT rd_substructcount('cC','Qc1ccccc1C');"))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute("SELECT rd_substructcount('nC','c1ccccn1C');"))


    pkl1 = DbModule.binaryHolder(Chem.MolFromSmiles('O').ToBinary())
    pkl2 = DbModule.binaryHolder(Chem.MolFromSmiles('OCCC(=O)O').ToBinary())
    cmd = "SELECT rd_substructcount(cast (%s as bytea),cast (%s as bytea))"
    curs.execute(cmd,(pkl1,pkl2))
    self.failUnless(curs.fetchone()[0]==3)
    pkl1 = DbModule.binaryHolder(Chem.MolFromSmiles('N').ToBinary())
    curs.execute(cmd,(pkl1,pkl2))
    self.failUnless(curs.fetchone()[0]==0)

    pkl1 = DbModule.binaryHolder(Chem.MolFromSmiles('O').ToBinary())
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute(cmd,('',pkl2)))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute(cmd,(pkl1,'')))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute(cmd,('','')))


  def test10(self):
    " substructures with pickles "
    curs = self.conn.GetCursor()

    curs.execute("SELECT rd_hassubstruct(rd_molpickle('O'),rd_molpickle('OCCC(=O)O'))")
    self.failUnless(curs.fetchone()[0])
    curs.execute("SELECT rd_substructcount(rd_molpickle('N'),rd_molpickle('OCCC(=O)O'))")
    self.failIf(curs.fetchone()[0])

  def test11(self):
    " descriptors "
    from Chem import Crippen,Descriptors
    curs = self.conn.GetCursor()

    smi = "c1ncccc1"
    m = Chem.MolFromSmiles(smi)
    pkl= DbModule.binaryHolder(m.ToBinary())
    ref = Crippen.MolLogP(m,addHs=1)
    curs.execute("SELECT rd_mollogp(%s)",(smi,))
    v = curs.fetchone()[0]
    self.failUnlessAlmostEqual(ref,v,4)
    curs.execute("SELECT rd_mollogp(cast (%s as bytea))",(pkl,))
    v = curs.fetchone()[0]
    self.failUnlessAlmostEqual(ref,v,4)
    curs.execute("SELECT rd_mollogp(rd_molpickle(%s))",(smi,))
    v = curs.fetchone()[0]
    self.failUnlessAlmostEqual(ref,v,4)

    ref = Descriptors.MolWt(m)
    curs.execute("SELECT rd_molamw(%s)",(smi,))
    v = curs.fetchone()[0]
    self.failUnlessAlmostEqual(ref,v,4)
    curs.execute("SELECT rd_molamw(cast (%s as bytea))",(pkl,))
    v = curs.fetchone()[0]
    self.failUnlessAlmostEqual(ref,v,4)
    curs.execute("SELECT rd_molamw(rd_molpickle(%s))",(smi,))
    v = curs.fetchone()[0]
    self.failUnlessAlmostEqual(ref,v,4)



    smi = "CCOCC(C)(C)C"
    m = Chem.MolFromSmiles(smi)
    pkl= DbModule.binaryHolder(m.ToBinary())
    ref = Crippen.MolLogP(m,addHs=1)
    curs.execute("SELECT rd_mollogp(%s)",(smi,))
    v = curs.fetchone()[0]
    self.failUnlessAlmostEqual(ref,v,4)
    curs.execute("SELECT rd_mollogp(cast (%s as bytea))",(pkl,))
    v = curs.fetchone()[0]
    self.failUnlessAlmostEqual(ref,v,4)
    curs.execute("SELECT rd_mollogp(rd_molpickle(%s))",(smi,))
    v = curs.fetchone()[0]
    self.failUnlessAlmostEqual(ref,v,4)

    ref = Descriptors.MolWt(m)
    curs.execute("SELECT rd_molamw(%s)",(smi,))
    v = curs.fetchone()[0]
    self.failUnlessAlmostEqual(ref,v,4)
    curs.execute("SELECT rd_molamw(cast (%s as bytea))",(pkl,))
    v = curs.fetchone()[0]
    self.failUnlessAlmostEqual(ref,v,4)
    curs.execute("SELECT rd_molamw(rd_molpickle(%s))",(smi,))
    v = curs.fetchone()[0]
    self.failUnlessAlmostEqual(ref,v,4)


    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute('select rd_mollogp(%s)',('',)))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute('select rd_mollogp(%s)',('QC',)))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute('select rd_mollogp(cast (%s as bytea))',('',)))
    self.failUnlessRaises(DbModule.OperationalError,
                          lambda : curs.execute('select rd_mollogp(cast (%s as bytea))',('randomtext',)))


  def test12(self):
    " using functions: "
    fnDef="""
    create or replace function rd_lipinskicount(smiles text) returns int as
$BODY$
   DECLARE
     res int;
   BEGIN
     res := 0;
     if rd_substructcount('[#7,#8]',smiles) > 10 then
       res := res + 1;
     end if;
     if rd_substructcount('[#7,#8;!H0]',smiles) > 5 then
       res := res + 1;
     end if;
     if rd_molamw(smiles) > 500 then
       res := res + 1;
     end if;
     if rd_mollogp(smiles) > 5 then
       res := res + 1;
     end if;
     RETURN res;
   END
$BODY$
  LANGUAGE plpgsql immutable;
    """
    curs = self.conn.GetCursor()
    curs.execute(fnDef)

    curs.execute("select rd_lipinskicount('Nc1cccc(O)c1')")
    self.failUnless(curs.fetchone()[0]==0)

    curs.execute("select rd_lipinskicount('Cc1ccc(-c2ccc(S(C)(=O)=O)cc2)n1-c1ccc(C)cc1')")
    self.failUnless(curs.fetchone()[0]==1)

    curs.execute("select rd_lipinskicount('Cc1n(-c2ccc(F)cc2)c(-c2ccc(S(C)(=O)=O)cc2)cc1C(Oc1cc(Cl)ccc1)C(F)(F)F')")
    self.failUnless(curs.fetchone()[0]==2)

    curs.execute("select rd_lipinskicount('Cc1cc(Cl)cc(-c2nc(C(F)(F)F)cn2-c2ccc(S(C)(=O)=O)cc2)c1')")
    self.failUnless(curs.fetchone()[0]==1)

if __name__ == '__main__':
  unittest.main()


