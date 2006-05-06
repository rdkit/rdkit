# $Id: MACCSkeys.py 5007 2006-02-22 15:14:41Z glandrum $
#
# Copyright (C) 2001-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" SMARTS definitions for the publically available MACCS keys
and a MACCS fingerprinter

"""
import Chem
import DataStructs
# these are SMARTS patterns corresponding to the MDL MACCS keys
smartsPatts={
  1:('?',0), # ISOTOPE
  #2:('[#103,#104,#105,#106,#107,#106,#109,#110,#111,#112]',0),  # ISOTOPE Not complete
  2:('[#103,#104]',0),  # ISOTOPE Not complete
  3:('[Ge,As,Se,Sn,Sb,Te,Tl,Pb,Bi]',0), # Group IVa,Va,VIa Periods 4-6 (Ge...)  *NOTE* spec wrong
  4:('[Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr]',0), # actinide
  5:('[Sc,Ti,Y,Zr,Hf]',0), # Group IIIB,IVB (Sc...)  *NOTE* spec wrong
  6:('[La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu]',0), # Lanthanide
  7:('[V,Cr,Mn,Nb,Mo,Tc,Ta,W,Re]',0), # Group VB,VIB,VIIB (V...) *NOTE* spec wrong
  8:('[!C;!c;!#1]1~*~*~*~1',0), # QAAA@1
  9:('[Fe,Co,Ni,Ru,Rh,Pd,Os,Ir,Pt]',0), # Group VIII (Fe...)
  10:('[Be,Mg,Ca,Sr,Ba,Ra]',0), # Group IIa (Alkaline earth)
  11:('*1~*~*~*~1',0), # 4M Ring
  12:('[Cu,Zn,Ag,Cd,Au,Hg]',0), # Group IB,IIB (Cu..)
  13:('[O,o]~[N,n](~[C,c])~[C,c]',0), # ON(C)C
  14:('[S,s]-[S,s]',0), # S-S
  15:('[O,o]~[C,c](~[O,o])~[O,o]',0), # OC(O)O
  16:('[!C;!c;!#1]1~*~*~1',0), # QAA@1
  17:('[C,c]#[C,c]',0), #CTC
  18:('[B,Al,Ga,In,Tl]',0), # Group IIIA (B...) *NOTE* spec wrong
  19:('*1~*~*~*~*~*~*~1',0), # 7M Ring
  20:('[Si]',0), #Si
  21:('[C,c]=[C,c](~[!C;!c;!#1])~[!C;!c;!#1]',0), # C=C(Q)Q
  22:('*1~*~*~1',0), # 3M Ring
  23:('[N,n]~[C,c](~[O,o])~[O,o]',0), # NC(O)O
  24:('[N,n]-[O,o]',0), # N-O
  25:('[N,n]~[C,c](~[N,n])~[N,n]',0), # NC(N)N
  26:('[C,c]=;@[C,c](@*)@*',0), # C$=C($A)$A
  27:('[I]',0), # I
  28:('[!C;!c;!#1]~[CH2]~[!C;!c;!#1]',0), # QCH2Q
  29:('P',0),# P
  30:('[C,c]~[!C;!c;!#1](~[C,c])(~[C,c])~*',0), # CQ(C)(C)A
  31:('[!C;!c;!#1]~[F,Cl,Br,I]',0), # QX
  32:('[C,c]~[S,s]~[N,n]',0), # CSN
  33:('[N,n]~[S,s]',0), # NS
  34:('[CH2]=*',0), # CH2=A
  35:('[Li,Na,K,Rb,Cs,Fr]',0), # Group IA (Alkali Metal)
  36:('[$(S@*),$(s@*)]',0), # S Heterocycle
  37:('[N,n]~[C,c](~[O,o])~[N,n]',0), # NC(O)N
  38:('[N,n]~[C,c](~[C,c])~[N,n]',0), # NC(C)N
  39:('[O,o]~[S,s](~[O,o])~[O,o]',0), # OS(O)O
  40:('[S,s]-[O,o]',0), # S-O
  41:('[C,c]#[N,n]',0), # CTN
  42:('F',0), # F
  43:('[!C;!c;!#1;H,H2,H3,H4]~*~[!C;!c;!#1;H,H2,H3,H4]',0), # QHAQH FIX: possibly incomplete
  44:('?',0), # OTHER
  45:('[C,c]=[C,c]~[N,n]',0), # C=CN
  46:('Br',0), # BR
  47:('[S,s]~*~[N,n]',0), # SAN
  48:('[O,o]~[!C;!c;!#1](~[O,o])(~[O,o])~*',0), # OQ(O)O
  49:('[-,-2,-3,-4,+,+2,+3,+4]',0), # CHARGE  FIX: possibly incomplete
  50:('[C,c]=[C,c](~[C,c])~[C,c]',0), # C=C(C)C
  51:('[C,c]~[S,s]~[O,o]',0), # CSO
  52:('[N,n]~[N,n]',0), # NN
  53:('[!C;!c;!#1;H,H2,H3,H4]~*~*~*~[!C;!c;!#1;H,H2,H3,H4]',0), # QHAAAQH FIX: possibly incomplete
  54:('[!C;!c;!#1;H,H2,H3,H4]~*~*~[!C;!c;!#1;H,H2,H3,H4]',0), # QHAAQH FIX: possibly incomplete
  55:('[O,o]~[S,s]~[O,o]',0), #OSO
  56:('[O,o]~[N,n](~[O,o])~[C,c]',0), # ON(O)C
  57:('[$(O@*),$(o@*)]',0), # O Heterocycle
  58:('[!C;!c;!#1]~[S,s]~[!C;!c;!#1]',0), # QSQ
  59:('[S,s]!:*:*',0), # Snot%A%A
  60:('[S,s]=[O,o]',0), # S=O
  61:('*~[S,s](~*)~*',0), # AS(A)A
  62:('*@*!@*@*',0), # A$!A$A
  63:('[N,n]=[O,o]',0), # N=O
  64:('*@*!@[S,s]',0), # A$A!S
  65:('[C,c]:[N,n]',0), # C%N
  66:('[C,c]~[C,c](~[C,c])(~[C,c])~*',0), # CC(C)(C)A
  67:('[!C;!c;!#1]~[S,s]',0), # QS
  68:('[!C;!c;!#1;H,H2,H3,H4]~[!C;!c;!#1;H,H2,H3,H4]',0), # QHQH FIX: possibly incomplete
  69:('[!C;!c;!#1]~[!C;!c;!#1;H,H2,H3,H4]',0), # QQH FIX: possibly incomplete
  70:('[!C;!c;!#1]~[N,n]~[!C;!c;!#1]',0), # QNQ
  71:('[N,n]~[O,o]',0), # NO
  72:('[O,o]~*~*~[O,o]',0), # OAAO
  73:('[S,s]=*',0), # S=A
  74:('[CH3]~*~[CH3]',0), # CH3ACH3
  75:('*!@[N,n]@*',0), # A!N$A
  76:('[C,c]=[C,c](~*)~*',0), # C=C(A)A
  77:('[N,n]~*~[N,n]',0), # NAN
  78:('[C,c]=[N,n]',0), # C=N
  79:('[N,n]~*~*~[N,n]',0), # NAAN
  80:('[N,n]~*~*~*~[N,n]',0), # NAAAN
  81:('[S,s]~*(~*)~*',0), # SA(A)A
  82:('*~[CH2]~[!C;!c;!#1;H,H2,H3,H4]',0), # ACH2QH
  83:('[!C;!c;!#1]1~*~*~*~*~1',0), # QAAAA@1
  84:('[NH2]',0), #NH2
  85:('[C,c]~[N,n](~[C,c])~[C,c]',0), # CN(C)C
  86:('[CH2][!C;!c;!#1][CH2]',0), # CH2QCH2
  87:('[F,Cl,Br,I]!@*@*',0), # X!A$A
  88:('[S,s]',0), # S
  89:('[O,o]~*~*~*~[O,o]',0), # OAAAO
  90:('[!C;!c;!#1;H,H2,H3,H4]~*~*~[CH2]~*',0), # QHAACH2A
  91:('[!C;!c;!#1;H,H2,H3,H4]~*~*~*~[CH2]~*',0), # QHAAACH2A
  92:('[O,o]~[C,c](~[N,n])~[C,c]',0), # OC(N)C
  93:('[!C;!c;!#1]~[CH3]',0), # QCH3
  94:('[!C;!c;!#1]~[N,n]',0), # QN
  95:('[N,n]~*~*~[O,o]',0), # NAAO
  96:('*1~*~*~*~*~1',0), # 5 M ring
  97:('[N,n]~*~*~*~[O,o]',0), # NAAAO
  98:('[!C;!c;!#1]1~*~*~*~*~*~1',0), # QAAAAA@1
  99:('[C,c]=[C,c]',0), # C=C
  100:('*~[CH2]~[N,n]',0), # ACH2N
  101:('[r8,r9,r10,r11,r12]',0), # 8M Ring or larger FIX: This is not exhaustive and it appears that oelib doesn't do this right
  102:('[!C;!c;!#1]~[O,o]',0), # QO
  103:('Cl',0), # CL
  104:('[!C;!c;!#1;H,H2,H3,H4]~*~[CH2]~*',0), # QHACH2A
  105:('[!C;!c;!#1]@*(@*)@*',0), # A$A($A)$A
  106:('[!C;!c;!#1]~*(~[!C;!c;!#1])~[!C;!c;!#1]',0), # QA(Q)Q
  107:('[F,Cl,Br,I]~*(~*)~*',0), # XA(A)A
  108:('[CH3]~*~*~*~[CH2]~*',0), # CH3AAACH2A
  109:('*~[CH2]~[O,o]',0), # ACH2O
  110:('[N,n]~[C,c]~[O,o]',0), # NCO
  111:('[N,n]~*~[CH2]~*',0), # NACH2A
  112:('*~*(~*)(~*)~*',0), # AA(A)(A)A
  113:('[O,o]!:*:*',0), # Onot%A%A
  114:('[CH3]~[CH2]~*',0), # CH3CH2A
  115:('[CH3]~*~[CH2]~*',0), # CH3ACH2A
  116:('[CH3]~*~*~[CH2]~*',0), # CH3AACH2A
  117:('[N,n]~*~[O,o]',0), # NAO
  118:('*~[CH2]~[CH2]~*',1), # ACH2CH2A > 1
  119:('[N,n]=*',0), # N=A
  120:('[!C;!c;R]',1), # Heterocyclic atom > 1
  121:('[$(N@*),$(n@*)]',0), # N Heterocycle
  122:('*~[N,n](~*)~*',0), # AN(A)A
  123:('[O,o]~[C,c]~[O,o]',0), # OCO
  124:('[!C;!c;!#1]~[!C;!c;!#1]',0), # QQ
  125:('?',0), # Aromatic Ring > 1
  126:('*!@[O,o]!@*',0), # A!O!A
  127:('*@*!@[O,o]',1), # A$A!O > 1
  128:('*~[CH2]~*~*~*~[CH2]~*',0), # ACH2AAACH2A
  129:('*~[CH2]~*~*~[CH2]~*',0), # ACH2AACH2A
  130:('[!C;!c;!#1]~[!C;!c;!#1]',1), # QQ > 1 (&...)
  131:('[!C;!c;!#1;H,H2,H3,H4]',1), # QH > 1
  132:('[O,o]~*~[CH2]~*',0), # OACH2A
  133:('*@*!:[N,n]',0), # A$A!N
  134:('[F,Cl,Br,I]',0), # X (HALOGEN)
  135:('[N,n]!:*:*',0), # Nnot%A%A
  136:('[O,o]=*',1), # O=A>1 FIX: maybe not right key
  137:('[!C;!c;R]',0), # Heterocycle
  138:('[!C;!c;!#1]~[CH2]~*',1), # QCH2A>1 (&...)
  139:('[OH,OH2,OH3]',0), # OH
  140:('[O,o]',3), # O > 3
  141:('[CH3]',2), # CH3 > 2
  142:('[N,n]',1), # N > 1
  143:('*@*!@[O,o]',0), # A$A!O
  144:('*!:*:*!:*',0), # Anot%A%Anot%A
  145:('*1~*~*~*~*~*~1',1), # 6M ring > 1
  146:('[O,o]',2), # O > 2
  147:('*~[CH2]~[CH2]~*',0), # ACH2CH2A
  148:('*~[!C;!c;!#1](~*)~*',0), # AQ(A)A
  149:('[CH3]',1), # CH3 > 1
  150:('*!@*@*!@*',0), # A!A$A!A
  151:('[NH,NH2,NH3,NH4]',0), # NH
  152:('[O,o]~[C,c](~[C,c])~[C,c]',0), # OC(C)C
  153:('[!C;!c;!#1]~[CH2]~*',0), # QCH2A
  154:('[C,c]=[O,o]',0), # C=O
  155:('*!@[CH2]!@*',0), # A!CH2!A
  156:('[N,n]~*(~*)~*',0), # NA(A)A
  157:('[C,c]-[O,o]',0), # C-O
  158:('[C,c]-[N,n]',0), # C-N
  159:('[O,o]',1), # O>1
  160:('[CH3]',0), #CH3
  161:('[N,n]',0), # N
  162:('a',0), # Aromatic
  163:('*1~*~*~*~*~*~1',0), # 6M Ring
  164:('[O,o]',0), # O
  165:('[R]',0), # Ring
  166:('?',0), # Fragments  FIX: this should be (*).(*), but that doesn't work properly in oelib
  }

maccsKeys = None

def _InitKeys(keyList,keyDict):
  """ *Internal Use Only*

   generates SMARTS patterns for the keys, run once

  """
  assert len(keyList) == len(keyDict.keys()),'length mismatch'
  for key in keyDict.keys():
    patt,count = keyDict[key]
    if patt != '?':
      try:
        sma = Chem.MolFromSmarts(patt)
      except:
        sma = None
      if not sma:
        print 'SMARTS parser error for key #%d: %s'%(key,patt)
      else:
        keyList[key-1] = sma,count
      
def GenMACCSKeys(mol,**kwargs):
  """ generates the MACCS fingerprint for a molecules

   **Arguments**

     - mol: the molecule to be fingerprinted

     - any extra keyword arguments are ignored
     
   **Returns**

      a _DataStructs.SparseBitVect_ containing the fingerprint.

  >>> m = Chem.MolFromSmiles('CNO')
  >>> bv = GenMACCSKeys(m)
  >>> tuple(bv.GetOnBits())
  (24, 68, 69, 71, 93, 94, 102, 124, 131, 139, 151, 158, 160, 161, 164)
  >>> bv = GenMACCSKeys(Chem.MolFromSmiles('CCC'))
  >>> tuple(bv.GetOnBits())
  (74, 114, 149, 155, 160)

  """
  global maccsKeys
  if maccsKeys is None:
    maccsKeys = [(None,0)]*len(smartsPatts.keys())
    _InitKeys(maccsKeys,smartsPatts)

  res = DataStructs.SparseBitVect(len(maccsKeys)+1)
  for i,(patt,count) in enumerate(maccsKeys):
    if patt is not None:
      matches = mol.GetSubstructMatches(patt)
      if matches:
        if count == 0:
          res[i+1] = 1
        else:
          if len(matches) > count:
            res[i+1] = 1
  return res

FingerprintMol = GenMACCSKeys
  
#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"])

if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)
