// $Id$
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Contribution from Roger Sayle
#include <vector>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolOps.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/no_tracking.hpp>

namespace  {
  struct Patterns {
    RDKit::ROMol *bit_8;
    RDKit::ROMol *bit_11;
    RDKit::ROMol *bit_13;
    RDKit::ROMol *bit_14;
    RDKit::ROMol *bit_15;
    RDKit::ROMol *bit_16;
    RDKit::ROMol *bit_17;
    RDKit::ROMol *bit_19;
    RDKit::ROMol *bit_20;
    RDKit::ROMol *bit_21;
    RDKit::ROMol *bit_22;
    RDKit::ROMol *bit_23;
    RDKit::ROMol *bit_24;
    RDKit::ROMol *bit_25;
    RDKit::ROMol *bit_26;
    RDKit::ROMol *bit_28;
    RDKit::ROMol *bit_30;
    RDKit::ROMol *bit_31;
    RDKit::ROMol *bit_32;
    RDKit::ROMol *bit_33;
    RDKit::ROMol *bit_34;
    RDKit::ROMol *bit_36;
    RDKit::ROMol *bit_37;
    RDKit::ROMol *bit_38;
    RDKit::ROMol *bit_39;
    RDKit::ROMol *bit_40;
    RDKit::ROMol *bit_41;
    RDKit::ROMol *bit_43;
    RDKit::ROMol *bit_44;
    RDKit::ROMol *bit_45;
    RDKit::ROMol *bit_47;
    RDKit::ROMol *bit_48;
    RDKit::ROMol *bit_49;
    RDKit::ROMol *bit_50;
    RDKit::ROMol *bit_51;
    RDKit::ROMol *bit_52;
    RDKit::ROMol *bit_53;
    RDKit::ROMol *bit_54;
    RDKit::ROMol *bit_55;
    RDKit::ROMol *bit_56;
    RDKit::ROMol *bit_57;
    RDKit::ROMol *bit_58;
    RDKit::ROMol *bit_59;
    RDKit::ROMol *bit_60;
    RDKit::ROMol *bit_61;
    RDKit::ROMol *bit_62;
    RDKit::ROMol *bit_63;
    RDKit::ROMol *bit_64;
    RDKit::ROMol *bit_65;
    RDKit::ROMol *bit_66;
    RDKit::ROMol *bit_67;
    RDKit::ROMol *bit_68;
    RDKit::ROMol *bit_69;
    RDKit::ROMol *bit_70;
    RDKit::ROMol *bit_71;
    RDKit::ROMol *bit_72;
    RDKit::ROMol *bit_73;
    RDKit::ROMol *bit_74;
    RDKit::ROMol *bit_75;
    RDKit::ROMol *bit_76;
    RDKit::ROMol *bit_77;
    RDKit::ROMol *bit_78;
    RDKit::ROMol *bit_79;
    RDKit::ROMol *bit_80;
    RDKit::ROMol *bit_81;
    RDKit::ROMol *bit_82;
    RDKit::ROMol *bit_83;
    RDKit::ROMol *bit_84;
    RDKit::ROMol *bit_85;
    RDKit::ROMol *bit_86;
    RDKit::ROMol *bit_87;
    RDKit::ROMol *bit_89;
    RDKit::ROMol *bit_90;
    RDKit::ROMol *bit_91;
    RDKit::ROMol *bit_92;
    RDKit::ROMol *bit_93;
    RDKit::ROMol *bit_94;
    RDKit::ROMol *bit_95;
    RDKit::ROMol *bit_96;
    RDKit::ROMol *bit_97;
    RDKit::ROMol *bit_98;
    RDKit::ROMol *bit_99;
    RDKit::ROMol *bit_100;
    RDKit::ROMol *bit_101;
    RDKit::ROMol *bit_102;
    RDKit::ROMol *bit_104;
    RDKit::ROMol *bit_105;
    RDKit::ROMol *bit_106;
    RDKit::ROMol *bit_107;
    RDKit::ROMol *bit_108;
    RDKit::ROMol *bit_109;
    RDKit::ROMol *bit_110;
    RDKit::ROMol *bit_111;
    RDKit::ROMol *bit_112;
    RDKit::ROMol *bit_113;
    RDKit::ROMol *bit_114;
    RDKit::ROMol *bit_115;
    RDKit::ROMol *bit_116;
    RDKit::ROMol *bit_117;
    RDKit::ROMol *bit_118;
    RDKit::ROMol *bit_119;
    RDKit::ROMol *bit_120;
    RDKit::ROMol *bit_121;
    RDKit::ROMol *bit_122;
    RDKit::ROMol *bit_123;
    RDKit::ROMol *bit_124;
    RDKit::ROMol *bit_126;
    RDKit::ROMol *bit_127;
    RDKit::ROMol *bit_128;
    RDKit::ROMol *bit_129;
    RDKit::ROMol *bit_131;
    RDKit::ROMol *bit_132;
    RDKit::ROMol *bit_133;
    RDKit::ROMol *bit_135;
    RDKit::ROMol *bit_136;
    RDKit::ROMol *bit_137;
    RDKit::ROMol *bit_138;
    RDKit::ROMol *bit_139;
    RDKit::ROMol *bit_140;
    RDKit::ROMol *bit_141;
    RDKit::ROMol *bit_142;
    RDKit::ROMol *bit_144;
    RDKit::ROMol *bit_145;
    RDKit::ROMol *bit_147;
    RDKit::ROMol *bit_148;
    RDKit::ROMol *bit_149;
    RDKit::ROMol *bit_150;
    RDKit::ROMol *bit_151;
    RDKit::ROMol *bit_152;
    RDKit::ROMol *bit_154;
    RDKit::ROMol *bit_155;
    RDKit::ROMol *bit_156;
    RDKit::ROMol *bit_157;
    RDKit::ROMol *bit_158;
    RDKit::ROMol *bit_162;
    RDKit::ROMol *bit_165;
    Patterns() :
      bit_8(RDKit::SmartsToMol("[!#6!#1]1~*~*~*~1")),
      bit_11(RDKit::SmartsToMol("*1~*~*~*~1")),
      bit_13(RDKit::SmartsToMol("[#8]~[#7](~[#6])~[#6]")),
      bit_14(RDKit::SmartsToMol("[#16]-[#16]")),
      bit_15(RDKit::SmartsToMol("[#8]~[#6](~[#8])~[#8]")),
      bit_16(RDKit::SmartsToMol("[!#6!#1]1~*~*~1")),
      bit_17(RDKit::SmartsToMol("[#6]#[#6]")),
      bit_19(RDKit::SmartsToMol("*1~*~*~*~*~*~*~1")),
      bit_20(RDKit::SmartsToMol("[#14]")),
      bit_21(RDKit::SmartsToMol("[#6]=[#6](~[!#6!#1])~[!#6!#1]")),
      bit_22(RDKit::SmartsToMol("*1~*~*~1")),
      bit_23(RDKit::SmartsToMol("[#7]~[#6](~[#8])~[#8]")),
      bit_24(RDKit::SmartsToMol("[#7]-[#8]")),
      bit_25(RDKit::SmartsToMol("[#7]~[#6](~[#7])~[#7]")),
      bit_26(RDKit::SmartsToMol("[#6]=@[#6](@*)@*")),
      bit_28(RDKit::SmartsToMol("[!#6!#1]~[CH2]~[!#6!#1]")),
      bit_30(RDKit::SmartsToMol("[#6]~[!#6!#1](~[#6])(~[#6])~*")),
      bit_31(RDKit::SmartsToMol("[!#6!#1]~[F,Cl,Br,I]")),
      bit_32(RDKit::SmartsToMol("[#6]~[#16]~[#7]")),
      bit_33(RDKit::SmartsToMol("[#7]~[#16]")),
      bit_34(RDKit::SmartsToMol("[CH2]=*")),
      bit_36(RDKit::SmartsToMol("[#16R]")),
      bit_37(RDKit::SmartsToMol("[#7]~[#6](~[#8])~[#7]")),
      bit_38(RDKit::SmartsToMol("[#7]~[#6](~[#6])~[#7]")),
      bit_39(RDKit::SmartsToMol("[#8]~[#16](~[#8])~[#8]")),
      bit_40(RDKit::SmartsToMol("[#16]-[#8]")),
      bit_41(RDKit::SmartsToMol("[#6]#[#7]")),
      bit_43(RDKit::SmartsToMol("[!#6!#1!H0]~*~[!#6!#1!H0]")),
      bit_44(RDKit::SmartsToMol("[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]")),
      bit_45(RDKit::SmartsToMol("[#6]=[#6]~[#7]")),
      bit_47(RDKit::SmartsToMol("[#16]~*~[#7]")),
      bit_48(RDKit::SmartsToMol("[#8]~[!#6!#1](~[#8])~[#8]")),
      bit_49(RDKit::SmartsToMol("[!+0]")),
      bit_50(RDKit::SmartsToMol("[#6]=[#6](~[#6])~[#6]")),
      bit_51(RDKit::SmartsToMol("[#6]~[#16]~[#8]")),
      bit_52(RDKit::SmartsToMol("[#7]~[#7]")),
      bit_53(RDKit::SmartsToMol("[!#6!#1!H0]~*~*~*~[!#6!#1!H0]")),
      bit_54(RDKit::SmartsToMol("[!#6!#1!H0]~*~*~[!#6!#1!H0]")),
      bit_55(RDKit::SmartsToMol("[#8]~[#16]~[#8]")),
      bit_56(RDKit::SmartsToMol("[#8]~[#7](~[#8])~[#6]")),
      bit_57(RDKit::SmartsToMol("[#8R]")),
      bit_58(RDKit::SmartsToMol("[!#6!#1]~[#16]~[!#6!#1]")),
      bit_59(RDKit::SmartsToMol("[#16]!:*:*")),
      bit_60(RDKit::SmartsToMol("[#16]=[#8]")),
      bit_61(RDKit::SmartsToMol("*~[#16](~*)~*")),
      bit_62(RDKit::SmartsToMol("*@*!@*@*")),
      bit_63(RDKit::SmartsToMol("[#7]=[#8]")),
      bit_64(RDKit::SmartsToMol("*@*!@[#16]")),
      bit_65(RDKit::SmartsToMol("c:n")),
      bit_66(RDKit::SmartsToMol("[#6]~[#6](~[#6])(~[#6])~*")),
      bit_67(RDKit::SmartsToMol("[!#6!#1]~[#16]")),
      bit_68(RDKit::SmartsToMol("[!#6!#1!H0]~[!#6!#1!H0]")),
      bit_69(RDKit::SmartsToMol("[!#6!#1]~[!#6!#1!H0]")),
      bit_70(RDKit::SmartsToMol("[!#6!#1]~[#7]~[!#6!#1]")),
      bit_71(RDKit::SmartsToMol("[#7]~[#8]")),
      bit_72(RDKit::SmartsToMol("[#8]~*~*~[#8]")),
      bit_73(RDKit::SmartsToMol("[#16]=*")),
      bit_74(RDKit::SmartsToMol("[CH3]~*~[CH3]")),
      bit_75(RDKit::SmartsToMol("*!@[#7]@*")),
      bit_76(RDKit::SmartsToMol("[#6]=[#6](~*)~*")),
      bit_77(RDKit::SmartsToMol("[#7]~*~[#7]")),
      bit_78(RDKit::SmartsToMol("[#6]=[#7]")),
      bit_79(RDKit::SmartsToMol("[#7]~*~*~[#7]")),
      bit_80(RDKit::SmartsToMol("[#7]~*~*~*~[#7]")),
      bit_81(RDKit::SmartsToMol("[#16]~*(~*)~*")),
      bit_82(RDKit::SmartsToMol("*~[CH2]~[!#6!#1!H0]")),
      bit_83(RDKit::SmartsToMol("[!#6!#1]1~*~*~*~*~1")),
      bit_84(RDKit::SmartsToMol("[NH2]")),
      bit_85(RDKit::SmartsToMol("[#6]~[#7](~[#6])~[#6]")),
      bit_86(RDKit::SmartsToMol("[C;H2,H3][!#6!#1][C;H2,H3]")),
      bit_87(RDKit::SmartsToMol("[F,Cl,Br,I]!@*@*")),
      bit_89(RDKit::SmartsToMol("[#8]~*~*~*~[#8]")),
      bit_90(RDKit::SmartsToMol("[$([!#6!#1!H0]~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[CH2R]1)]")),
      bit_91(RDKit::SmartsToMol("[$([!#6!#1!H0]~*~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~*~[R]1@[R]@[CH2R]1)]")),
      bit_92(RDKit::SmartsToMol("[#8]~[#6](~[#7])~[#6]")),
      bit_93(RDKit::SmartsToMol("[!#6!#1]~[CH3]")),
      bit_94(RDKit::SmartsToMol("[!#6!#1]~[#7]")),
      bit_95(RDKit::SmartsToMol("[#7]~*~*~[#8]")),
      bit_96(RDKit::SmartsToMol("*1~*~*~*~*~1")),
      bit_97(RDKit::SmartsToMol("[#7]~*~*~*~[#8]")),
      bit_98(RDKit::SmartsToMol("[!#6!#1]1~*~*~*~*~*~1")),
      bit_99(RDKit::SmartsToMol("[#6]=[#6]")),
      bit_100(RDKit::SmartsToMol("*~[CH2]~[#7]")),
      bit_101(RDKit::SmartsToMol("[$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1)]")),
      bit_102(RDKit::SmartsToMol("[!#6!#1]~[#8]")),
      bit_104(RDKit::SmartsToMol("[!#6!#1!H0]~*~[CH2]~*")),
      bit_105(RDKit::SmartsToMol("*@*(@*)@*")),
      bit_106(RDKit::SmartsToMol("[!#6!#1]~*(~[!#6!#1])~[!#6!#1]")),
      bit_107(RDKit::SmartsToMol("[F,Cl,Br,I]~*(~*)~*")),
      bit_108(RDKit::SmartsToMol("[CH3]~*~*~*~[CH2]~*")),
      bit_109(RDKit::SmartsToMol("*~[CH2]~[#8]")),
      bit_110(RDKit::SmartsToMol("[#7]~[#6]~[#8]")),
      bit_111(RDKit::SmartsToMol("[#7]~*~[CH2]~*")),
      bit_112(RDKit::SmartsToMol("*~*(~*)(~*)~*")),
      bit_113(RDKit::SmartsToMol("[#8]!:*:*")),
      bit_114(RDKit::SmartsToMol("[CH3]~[CH2]~*")),
      bit_115(RDKit::SmartsToMol("[CH3]~*~[CH2]~*")),
      bit_116(RDKit::SmartsToMol("[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]")),
      bit_117(RDKit::SmartsToMol("[#7]~*~[#8]")),
      bit_118(RDKit::SmartsToMol("[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]")),
      bit_119(RDKit::SmartsToMol("[#7]=*")),
      bit_120(RDKit::SmartsToMol("[!#6R]")),
      bit_121(RDKit::SmartsToMol("[#7R]")),
      bit_122(RDKit::SmartsToMol("*~[#7](~*)~*")),
      bit_123(RDKit::SmartsToMol("[#8]~[#6]~[#8]")),
      bit_124(RDKit::SmartsToMol("[!#6!#1]~[!#6!#1]")),
      bit_126(RDKit::SmartsToMol("*!@[#8]!@*")),
      bit_127(RDKit::SmartsToMol("*@*!@[#8]")),
      bit_128(RDKit::SmartsToMol("[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2R]@[R]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2R]1),$(*~[CH2]~*~[R]1@[R]@[CH2R]1)]")),
      bit_129(RDKit::SmartsToMol("[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[CH2R]1)]")),
      bit_131(RDKit::SmartsToMol("[!#6!#1!H0]")),
      bit_132(RDKit::SmartsToMol("[#8]~*~[CH2]~*")),
      bit_133(RDKit::SmartsToMol("*@*!@[#7]")),
      bit_135(RDKit::SmartsToMol("[#7]!:*:*")),
      bit_136(RDKit::SmartsToMol("[#8]=*")),
      bit_137(RDKit::SmartsToMol("[!C!cR]")),
      bit_138(RDKit::SmartsToMol("[!#6!#1]~[CH2]~*")),
      bit_139(RDKit::SmartsToMol("[O!H0]")),
      bit_140(RDKit::SmartsToMol("[#8]")),
      bit_141(RDKit::SmartsToMol("[CH3]")),
      bit_142(RDKit::SmartsToMol("[#7]")),
      bit_144(RDKit::SmartsToMol("*!:*:*!:*")),
      bit_145(RDKit::SmartsToMol("*1~*~*~*~*~*~1")),
      bit_147(RDKit::SmartsToMol("[$(*~[CH2]~[CH2]~*),$([R]1@[CH2R]@[CH2R]1)]")),
      bit_148(RDKit::SmartsToMol("*~[!#6!#1](~*)~*")),
      bit_149(RDKit::SmartsToMol("[C;H3,H4]")),
      bit_150(RDKit::SmartsToMol("*!@*@*!@*")),
      bit_151(RDKit::SmartsToMol("[#7!H0]")),
      bit_152(RDKit::SmartsToMol("[#8]~[#6](~[#6])~[#6]")),
      bit_154(RDKit::SmartsToMol("[#6]=[#8]")),
      bit_155(RDKit::SmartsToMol("*!@[CH2]!@*")),
      bit_156(RDKit::SmartsToMol("[#7]~*(~*)~*")),
      bit_157(RDKit::SmartsToMol("[#6]-[#8]")),
      bit_158(RDKit::SmartsToMol("[#6]-[#7]")),
      bit_162(RDKit::SmartsToMol("a")),
      bit_165(RDKit::SmartsToMol("[R]")) {}
  };

  boost::flyweight<std::vector<Patterns *>,boost::flyweights::no_tracking> gpats;
  void GenerateFP(const RDKit::ROMol &mol,ExplicitBitVect &fp)
  {
    if(gpats.get().size()==0){
      std::vector<Patterns *> ps;
      ps.push_back(new Patterns());
      gpats = ps;
    }
    const Patterns &pats=*(gpats.get().front());
    PRECONDITION(fp.size()==167,"bad fingerprint");
    fp.clearBits();

    if(!mol.getNumAtoms()) return;

    std::vector<RDKit::MatchVectType> matches;
    RDKit::RWMol::ConstAtomIterator atom;
    RDKit::MatchVectType match;
    unsigned int count;

    for (atom=mol.beginAtoms();atom!=mol.endAtoms();++atom)
      switch ((*atom)->getAtomicNum()) {
      case 3:
      case 11:
      case 19:
      case 37:
      case 55:
      case 87:
        fp.setBit(35);
        break;
      case 4:
      case 12:
      case 20:
      case 38:
      case 56:
      case 88:
        fp.setBit(10);
        break;
      case 5:
      case 13:
      case 31:
      case 49:
      case 81:
        fp.setBit(18);
        break;
      case 9:
        fp.setBit(42);
        fp.setBit(134);
        break;
      case 15:
        fp.setBit(29);
        break;
      case 16:
        fp.setBit(88);
        break;
      case 17:
        fp.setBit(103);
        fp.setBit(134);
        break;
      case 21:
      case 22:
      case 39:
      case 40:
      case 72:
        fp.setBit(5);
        break;
      case 23:
      case 24:
      case 25:
      case 41:
      case 42:
      case 43:
      case 73:
      case 74:
      case 75:
        fp.setBit(7);
        break;
      case 26:
      case 27:
      case 28:
      case 44:
      case 45:
      case 46:
      case 76:
      case 77:
      case 78:
        fp.setBit(9);
        break;
      case 29:
      case 30:
      case 47:
      case 48:
      case 79:
      case 80:
        fp.setBit(12);
        break;
      case 32:
      case 33:
      case 34:
      case 50:
      case 51:
      case 52:
      case 82:
      case 83:
      case 84:
        fp.setBit(3);
        break;
      case 35:
        fp.setBit(46);
        fp.setBit(134);
        break;
      case 53:
        fp.setBit(27);
        fp.setBit(134);
        break;
      case 57:
      case 58:
      case 59:
      case 60:
      case 61:
      case 62:
      case 63:
      case 64:
      case 65:
      case 66:
      case 67:
      case 68:
      case 69:
      case 70:
      case 71:
        fp.setBit(6);
        break;
      case 89:
      case 90:
      case 91:
      case 92:
      case 93:
      case 94:
      case 95:
      case 96:
      case 97:
      case 98:
      case 99:
      case 100:
      case 101:
      case 102:
      case 103:
        fp.setBit(4);
        break;
      case 104:
        fp.setBit(2);
        break;
      }

    if (RDKit::SubstructMatch(mol,*pats.bit_8,match,true))
      fp.setBit(8);
    if (RDKit::SubstructMatch(mol,*pats.bit_11,match,true))
      fp.setBit(11);
    if (RDKit::SubstructMatch(mol,*pats.bit_13,match,true))
      fp.setBit(13);
    if (RDKit::SubstructMatch(mol,*pats.bit_14,match,true))
      fp.setBit(14);
    if (RDKit::SubstructMatch(mol,*pats.bit_15,match,true))
      fp.setBit(15);
    if (RDKit::SubstructMatch(mol,*pats.bit_16,match,true))
      fp.setBit(16);
    if (RDKit::SubstructMatch(mol,*pats.bit_17,match,true))
      fp.setBit(17);
    if (RDKit::SubstructMatch(mol,*pats.bit_19,match,true))
      fp.setBit(19);
    if (RDKit::SubstructMatch(mol,*pats.bit_20,match,true))
      fp.setBit(20);
    if (RDKit::SubstructMatch(mol,*pats.bit_21,match,true))
      fp.setBit(21);
    if (RDKit::SubstructMatch(mol,*pats.bit_22,match,true))
      fp.setBit(22);
    if (RDKit::SubstructMatch(mol,*pats.bit_23,match,true))
      fp.setBit(23);
    if (RDKit::SubstructMatch(mol,*pats.bit_24,match,true))
      fp.setBit(24);
    if (RDKit::SubstructMatch(mol,*pats.bit_25,match,true))
      fp.setBit(25);
    if (RDKit::SubstructMatch(mol,*pats.bit_26,match,true))
      fp.setBit(26);
    if (RDKit::SubstructMatch(mol,*pats.bit_28,match,true))
      fp.setBit(28);
    if (RDKit::SubstructMatch(mol,*pats.bit_30,match,true))
      fp.setBit(30);
    if (RDKit::SubstructMatch(mol,*pats.bit_31,match,true))
      fp.setBit(31);
    if (RDKit::SubstructMatch(mol,*pats.bit_32,match,true))
      fp.setBit(32);
    if (RDKit::SubstructMatch(mol,*pats.bit_33,match,true))
      fp.setBit(33);
    if (RDKit::SubstructMatch(mol,*pats.bit_34,match,true))
      fp.setBit(34);
    if (RDKit::SubstructMatch(mol,*pats.bit_36,match,true))
      fp.setBit(36);
    if (RDKit::SubstructMatch(mol,*pats.bit_37,match,true))
      fp.setBit(37);
    if (RDKit::SubstructMatch(mol,*pats.bit_38,match,true))
      fp.setBit(38);
    if (RDKit::SubstructMatch(mol,*pats.bit_39,match,true))
      fp.setBit(39);
    if (RDKit::SubstructMatch(mol,*pats.bit_40,match,true))
      fp.setBit(40);
    if (RDKit::SubstructMatch(mol,*pats.bit_41,match,true))
      fp.setBit(41);
    if (RDKit::SubstructMatch(mol,*pats.bit_43,match,true))
      fp.setBit(43);
    if (RDKit::SubstructMatch(mol,*pats.bit_44,match,true))
      fp.setBit(44);
    if (RDKit::SubstructMatch(mol,*pats.bit_45,match,true))
      fp.setBit(45);
    if (RDKit::SubstructMatch(mol,*pats.bit_47,match,true))
      fp.setBit(47);
    if (RDKit::SubstructMatch(mol,*pats.bit_48,match,true))
      fp.setBit(48);
    if (RDKit::SubstructMatch(mol,*pats.bit_49,match,true))
      fp.setBit(49);
    if (RDKit::SubstructMatch(mol,*pats.bit_50,match,true))
      fp.setBit(50);
    if (RDKit::SubstructMatch(mol,*pats.bit_51,match,true))
      fp.setBit(51);
    if (RDKit::SubstructMatch(mol,*pats.bit_52,match,true))
      fp.setBit(52);
    if (RDKit::SubstructMatch(mol,*pats.bit_53,match,true))
      fp.setBit(53);
    if (RDKit::SubstructMatch(mol,*pats.bit_54,match,true))
      fp.setBit(54);
    if (RDKit::SubstructMatch(mol,*pats.bit_55,match,true))
      fp.setBit(55);
    if (RDKit::SubstructMatch(mol,*pats.bit_56,match,true))
      fp.setBit(56);
    if (RDKit::SubstructMatch(mol,*pats.bit_57,match,true))
      fp.setBit(57);
    if (RDKit::SubstructMatch(mol,*pats.bit_58,match,true))
      fp.setBit(58);
    if (RDKit::SubstructMatch(mol,*pats.bit_59,match,true))
      fp.setBit(59);
    if (RDKit::SubstructMatch(mol,*pats.bit_60,match,true))
      fp.setBit(60);
    if (RDKit::SubstructMatch(mol,*pats.bit_61,match,true))
      fp.setBit(61);
    if (RDKit::SubstructMatch(mol,*pats.bit_62,match,true))
      fp.setBit(62);
    if (RDKit::SubstructMatch(mol,*pats.bit_63,match,true))
      fp.setBit(63);
    if (RDKit::SubstructMatch(mol,*pats.bit_64,match,true))
      fp.setBit(64);
    if (RDKit::SubstructMatch(mol,*pats.bit_65,match,true))
      fp.setBit(65);
    if (RDKit::SubstructMatch(mol,*pats.bit_66,match,true))
      fp.setBit(66);
    if (RDKit::SubstructMatch(mol,*pats.bit_67,match,true))
      fp.setBit(67);
    if (RDKit::SubstructMatch(mol,*pats.bit_68,match,true))
      fp.setBit(68);
    if (RDKit::SubstructMatch(mol,*pats.bit_69,match,true))
      fp.setBit(69);
    if (RDKit::SubstructMatch(mol,*pats.bit_70,match,true))
      fp.setBit(70);
    if (RDKit::SubstructMatch(mol,*pats.bit_71,match,true))
      fp.setBit(71);
    if (RDKit::SubstructMatch(mol,*pats.bit_72,match,true))
      fp.setBit(72);
    if (RDKit::SubstructMatch(mol,*pats.bit_73,match,true))
      fp.setBit(73);
    if (RDKit::SubstructMatch(mol,*pats.bit_74,match,true))
      fp.setBit(74);
    if (RDKit::SubstructMatch(mol,*pats.bit_75,match,true))
      fp.setBit(75);
    if (RDKit::SubstructMatch(mol,*pats.bit_76,match,true))
      fp.setBit(76);
    if (RDKit::SubstructMatch(mol,*pats.bit_77,match,true))
      fp.setBit(77);
    if (RDKit::SubstructMatch(mol,*pats.bit_78,match,true))
      fp.setBit(78);
    if (RDKit::SubstructMatch(mol,*pats.bit_79,match,true))
      fp.setBit(79);
    if (RDKit::SubstructMatch(mol,*pats.bit_80,match,true))
      fp.setBit(80);
    if (RDKit::SubstructMatch(mol,*pats.bit_81,match,true))
      fp.setBit(81);
    if (RDKit::SubstructMatch(mol,*pats.bit_82,match,true))
      fp.setBit(82);
    if (RDKit::SubstructMatch(mol,*pats.bit_83,match,true))
      fp.setBit(83);
    if (RDKit::SubstructMatch(mol,*pats.bit_84,match,true))
      fp.setBit(84);
    if (RDKit::SubstructMatch(mol,*pats.bit_85,match,true))
      fp.setBit(85);
    if (RDKit::SubstructMatch(mol,*pats.bit_86,match,true))
      fp.setBit(86);
    if (RDKit::SubstructMatch(mol,*pats.bit_87,match,true))
      fp.setBit(87);
    if (RDKit::SubstructMatch(mol,*pats.bit_89,match,true))
      fp.setBit(89);
    if (RDKit::SubstructMatch(mol,*pats.bit_90,match,true))
      fp.setBit(90);
    if (RDKit::SubstructMatch(mol,*pats.bit_91,match,true))
      fp.setBit(91);
    if (RDKit::SubstructMatch(mol,*pats.bit_92,match,true))
      fp.setBit(92);
    if (RDKit::SubstructMatch(mol,*pats.bit_93,match,true))
      fp.setBit(93);
    if (RDKit::SubstructMatch(mol,*pats.bit_94,match,true))
      fp.setBit(94);
    if (RDKit::SubstructMatch(mol,*pats.bit_95,match,true))
      fp.setBit(95);
    if (RDKit::SubstructMatch(mol,*pats.bit_96,match,true))
      fp.setBit(96);
    if (RDKit::SubstructMatch(mol,*pats.bit_97,match,true))
      fp.setBit(97);
    if (RDKit::SubstructMatch(mol,*pats.bit_98,match,true))
      fp.setBit(98);
    if (RDKit::SubstructMatch(mol,*pats.bit_99,match,true))
      fp.setBit(99);
    if (RDKit::SubstructMatch(mol,*pats.bit_100,match,true))
      fp.setBit(100);
    if (RDKit::SubstructMatch(mol,*pats.bit_101,match,true))
      fp.setBit(101);
    if (RDKit::SubstructMatch(mol,*pats.bit_102,match,true))
      fp.setBit(102);
    if (RDKit::SubstructMatch(mol,*pats.bit_104,match,true))
      fp.setBit(104);
    if (RDKit::SubstructMatch(mol,*pats.bit_105,match,true))
      fp.setBit(105);
    if (RDKit::SubstructMatch(mol,*pats.bit_106,match,true))
      fp.setBit(106);
    if (RDKit::SubstructMatch(mol,*pats.bit_107,match,true))
      fp.setBit(107);
    if (RDKit::SubstructMatch(mol,*pats.bit_108,match,true))
      fp.setBit(108);
    if (RDKit::SubstructMatch(mol,*pats.bit_109,match,true))
      fp.setBit(109);
    if (RDKit::SubstructMatch(mol,*pats.bit_110,match,true))
      fp.setBit(110);
    if (RDKit::SubstructMatch(mol,*pats.bit_111,match,true))
      fp.setBit(111);
    if (RDKit::SubstructMatch(mol,*pats.bit_112,match,true))
      fp.setBit(112);
    if (RDKit::SubstructMatch(mol,*pats.bit_113,match,true))
      fp.setBit(113);
    if (RDKit::SubstructMatch(mol,*pats.bit_114,match,true))
      fp.setBit(114);
    if (RDKit::SubstructMatch(mol,*pats.bit_115,match,true))
      fp.setBit(115);
    if (RDKit::SubstructMatch(mol,*pats.bit_116,match,true))
      fp.setBit(116);
    if (RDKit::SubstructMatch(mol,*pats.bit_117,match,true))
      fp.setBit(117);
    if (RDKit::SubstructMatch(mol,*pats.bit_118,matches,true,true) > 1)
      fp.setBit(118);
    if (RDKit::SubstructMatch(mol,*pats.bit_119,match,true))
      fp.setBit(119);
    if (RDKit::SubstructMatch(mol,*pats.bit_120,matches,true,true) > 1)
      fp.setBit(120);
    if (RDKit::SubstructMatch(mol,*pats.bit_121,match,true))
      fp.setBit(121);
    if (RDKit::SubstructMatch(mol,*pats.bit_122,match,true))
      fp.setBit(122);
    if (RDKit::SubstructMatch(mol,*pats.bit_123,match,true))
      fp.setBit(123);
    count = RDKit::SubstructMatch(mol,*pats.bit_124,matches,true,true);
    if (count > 0)
      fp.setBit(124);
    if (count > 1)
      fp.setBit(130);
    if (RDKit::SubstructMatch(mol,*pats.bit_126,match,true))
      fp.setBit(126);
    count = RDKit::SubstructMatch(mol,*pats.bit_127,matches,true,true);
    if (count > 1)
      fp.setBit(127);
    if (count > 0)
      fp.setBit(143);
    if (RDKit::SubstructMatch(mol,*pats.bit_128,match,true))
      fp.setBit(128);
    if (RDKit::SubstructMatch(mol,*pats.bit_129,match,true))
      fp.setBit(129);
    if (RDKit::SubstructMatch(mol,*pats.bit_131,matches,true,true) > 1)
      fp.setBit(131);
    if (RDKit::SubstructMatch(mol,*pats.bit_132,match,true))
      fp.setBit(132);
    if (RDKit::SubstructMatch(mol,*pats.bit_133,match,true))
      fp.setBit(133);
    if (RDKit::SubstructMatch(mol,*pats.bit_135,match,true))
      fp.setBit(135);
    if (RDKit::SubstructMatch(mol,*pats.bit_136,matches,true,true) > 1)
      fp.setBit(136);
    if (RDKit::SubstructMatch(mol,*pats.bit_137,match,true))
      fp.setBit(137);
    count = RDKit::SubstructMatch(mol,*pats.bit_138,matches,true,true);
    if (count > 1)
      fp.setBit(138);
    if (count > 0)
      fp.setBit(153);
    if (RDKit::SubstructMatch(mol,*pats.bit_139,match,true))
      fp.setBit(139);
    count = RDKit::SubstructMatch(mol,*pats.bit_140,matches,true,true);
    if (count > 3)
      fp.setBit(140);
    if (count > 2)
      fp.setBit(146);
    if (count > 1)
      fp.setBit(159);
    if (count > 0)
      fp.setBit(164);
    if (RDKit::SubstructMatch(mol,*pats.bit_141,matches,true,true) > 2)
      fp.setBit(141);
    count = RDKit::SubstructMatch(mol,*pats.bit_142,matches,true,true);
    if (count > 1)
      fp.setBit(142);
    if (count > 0)
      fp.setBit(161);
    if (RDKit::SubstructMatch(mol,*pats.bit_144,match,true))
      fp.setBit(144);
    count = RDKit::SubstructMatch(mol,*pats.bit_145,matches,true,true);
    if (count > 1)
      fp.setBit(145);
    if (count > 0)
      fp.setBit(163);
    if (RDKit::SubstructMatch(mol,*pats.bit_147,match,true))
      fp.setBit(147);
    if (RDKit::SubstructMatch(mol,*pats.bit_148,match,true))
      fp.setBit(148);
    count = RDKit::SubstructMatch(mol,*pats.bit_149,matches,true,true);
    if (count > 1)
      fp.setBit(149);
    if (count > 0)
      fp.setBit(160);
    if (RDKit::SubstructMatch(mol,*pats.bit_150,match,true))
      fp.setBit(150);
    if (RDKit::SubstructMatch(mol,*pats.bit_151,match,true))
      fp.setBit(151);
    if (RDKit::SubstructMatch(mol,*pats.bit_152,match,true))
      fp.setBit(152);
    if (RDKit::SubstructMatch(mol,*pats.bit_154,match,true))
      fp.setBit(154);
    if (RDKit::SubstructMatch(mol,*pats.bit_155,match,true))
      fp.setBit(155);
    if (RDKit::SubstructMatch(mol,*pats.bit_156,match,true))
      fp.setBit(156);
    if (RDKit::SubstructMatch(mol,*pats.bit_157,match,true))
      fp.setBit(157);
    if (RDKit::SubstructMatch(mol,*pats.bit_158,match,true))
      fp.setBit(158);
    if (RDKit::SubstructMatch(mol,*pats.bit_162,match,true))
      fp.setBit(162);
    if (RDKit::SubstructMatch(mol,*pats.bit_165,match,true))
      fp.setBit(165);

    /* BIT 125 */
    RDKit::RingInfo *info = mol.getRingInfo();
    unsigned int ringcount = info->numRings();
    unsigned int nArom = 0;
    for (unsigned int i= 0; i<ringcount; i++) {
      bool isArom = true;
      const std::vector<int> *ring = &info->bondRings()[i];
      std::vector<int>::const_iterator iter;
      for (iter=ring->begin(); iter!=ring->end(); ++iter)
        if (!mol.getBondWithIdx(*iter)->getIsAromatic()) {
          isArom = false;
          break;
        }
      if (isArom) {
        if (nArom) {
          fp.setBit(125);
          break;
        } else nArom++;
      }
    }

    /* BIT 166 */
    std::vector<int> mapping;
    if (RDKit::MolOps::getMolFrags(mol,mapping) > 1)
      fp.setBit(166);
  }
} //end of local anonymous namespace

namespace RDKit {
  namespace MACCSFingerprints {
    ExplicitBitVect *getFingerprintAsBitVect(const ROMol &mol){
      ExplicitBitVect *fp=new ExplicitBitVect(167);
      GenerateFP(mol,*fp);
      return fp;
    }
  }
}
