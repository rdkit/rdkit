// Auto-generated query functions for Abraham features
// These queries are in the EXACT order the trained model expects

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <vector>
#include <memory>
#include <mutex>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Get 241 base feature queries in EXACT model order
// Thread-safe initialization using C++11 static local variable (guaranteed thread-safe)
const std::vector<std::shared_ptr<RWMol>>& GetQueriesAbrahamBaseFeatures() {
    static std::vector<std::shared_ptr<RWMol>> queries;
    static std::once_flag init_flag;
    
    std::call_once(init_flag, []() {
        queries.reserve(241);
        
        // [0] abraham_A_additional_0
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH][CX4][CX4][OH]")));
        // [1] abraham_A_additional_1
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1ccccc1[OH]")));
        // [2] abraham_A_additional_10
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1ccccc1[NX3](=O)=O")));
        // [3] abraham_A_additional_11
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("c1ccccc1[OH]c1ccccc1")));
        // [4] abraham_A_additional_12
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH][CX3]=[CX3]")));
        // [5] abraham_A_additional_13
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[OH][CX4][OH]")));
        // [6] abraham_A_additional_14
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[OH][NX3]")));
        // [7] abraham_A_additional_2
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1ccc([OH])cc1")));
        // [8] abraham_A_additional_3
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3;H2][CX4][OH]")));
        // [9] abraham_A_additional_4
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3;H1][CX4][OH]")));
        // [10] abraham_A_additional_5
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[OH]c1ccccc1")));
        // [11] abraham_A_additional_6
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[OH]c1ccccc1[OH]")));
        // [12] abraham_A_additional_7
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH][CX4][CX3](=O)")));
        // [13] abraham_A_additional_8
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1ccccc1[CX3](=O)")));
        // [14] abraham_A_additional_9
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1ccccc1[F,Cl,Br,I]")));
        // [15] abraham_A_frag_0
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[C][OX2H]")));
        // [16] abraham_A_frag_1
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[c][OX2H]")));
        // [17] abraham_A_frag_10
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=[OX1])[NX3;H1][C]")));
        // [18] abraham_A_frag_11
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=[OX1])[NX3;H1][c]")));
        // [19] abraham_A_frag_12
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[$([SX4](=[OX1])(=[OX1])([!O])[NH,NH2,NH3+]),$([SX4+2]([OX1-])([OX1-])([!O])[NH,NH2,NH3+])]")));
        // [20] abraham_A_frag_13
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3;H1]C(=[OX1])[NX3;H1]")));
        // [21] abraham_A_frag_14
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3;H0]C(=[OX1])[NX3;H1]")));
        // [22] abraham_A_frag_15
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3;H1]C(=[OX1])O")));
        // [23] abraham_A_frag_16
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3;H1]C(=N)[NX3;H0]")));
        // [24] abraham_A_frag_17
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[C]#[CH]")));
        // [25] abraham_A_frag_18
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("P[OH,O-]")));
        // [26] abraham_A_frag_19
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CH][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]")));
        // [27] abraham_A_frag_2
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[C][NX3;H2]")));
        // [28] abraham_A_frag_20
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CH]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]")));
        // [29] abraham_A_frag_21
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([CX3](=O)[OX1H0-,OX2H1])[CX4][CX3](=O)[OX1H0-,OX2H1]")));
        // [30] abraham_A_frag_22
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX3](=O)[OX1H0-,OX2H1]")));
        // [31] abraham_A_frag_23
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[OH]")));
        // [32] abraham_A_frag_24
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][OH]")));
        // [33] abraham_A_frag_25
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[nX3;H1]:n")));
        // [34] abraham_A_frag_26
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[nX3;H1]:c:n")));
        // [35] abraham_A_frag_27
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2;H1]CC[O,N]")));
        // [36] abraham_A_frag_28
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2;H1]C[C,N]=[O,S]")));
        // [37] abraham_A_frag_29
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2;H1]c1ccccc1[O,NX3]")));
        // [38] abraham_A_frag_3
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[c][NX3;H2;!$(NC=O)]")));
        // [39] abraham_A_frag_30
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2;H1]c1ccccc1C=[O,S]")));
        // [40] abraham_A_frag_31
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2;H1]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]")));
        // [41] abraham_A_frag_32
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NH,NH2,NH3+]CC[O,N]")));
        // [42] abraham_A_frag_33
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NH,NH2,NH3+]c1ccccc1[O,N]")));
        // [43] abraham_A_frag_34
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NH,NH2,NH3+]c1ccccc1[C,N]=[O,S]")));
        // [44] abraham_A_frag_35
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2H]c1ccccc1[Cl,Br,I]")));
        // [45] abraham_A_frag_36
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX1]=[C,c]~[C,c]C[OH]")));
        // [46] abraham_A_frag_37
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1cccc2cccnc12")));
        // [47] abraham_A_frag_38
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1")));
        // [48] abraham_A_frag_39
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1")));
        // [49] abraham_A_frag_4
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[C][NX3;H1;!R][C]")));
        // [50] abraham_A_frag_40
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NH,NH2,NH3+]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1")));
        // [51] abraham_A_frag_41
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NH,NH2,NH3+]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1")));
        // [52] abraham_A_frag_42
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)([OX1H0-,OX2H1])c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1")));
        // [53] abraham_A_frag_43
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)([OX1H0-,OX2H1])c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1")));
        // [54] abraham_A_frag_44
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1c([CX4])cccc1[CX4]")));
        // [55] abraham_A_frag_45
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NH,NH2,NH3+]c1c([CX4])cccc1[CX4]")));
        // [56] abraham_A_frag_46
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1c(C[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cccc1")));
        // [57] abraham_A_frag_47
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1cc([CX3](=O)[OX1H0-,OX2H1])ccc1")));
        // [58] abraham_A_frag_48
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1ccc([CX3](=O)[OX1H0-,OX2H1])cc1")));
        // [59] abraham_A_frag_49
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1cc([$([CH](=O)),$(C(=O)C)])ccc1")));
        // [60] abraham_A_frag_5
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[C][NX3;H1;R][C]")));
        // [61] abraham_A_frag_50
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1ccc([$([CH](=O)),$(C(=O)C)])cc1")));
        // [62] abraham_A_frag_6
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[c][NX3;H1;!$(NC=O)][C]")));
        // [63] abraham_A_frag_7
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[c][nX3;H1][c]")));
        // [64] abraham_A_frag_8
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[OX1H0-,OX2H1]")));
        // [65] abraham_A_frag_9
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=[OX1])[NX3;H2]")));
        // [66] abraham_BSEL_frag_0
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H3]")));
        // [67] abraham_BSEL_frag_1
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H2]")));
        // [68] abraham_BSEL_frag_10
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[C][NX3;H1][C]")));
        // [69] abraham_BSEL_frag_11
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[c][NX3;H1]")));
        // [70] abraham_BSEL_frag_12
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[C][NX3;H0](C)[C]")));
        // [71] abraham_BSEL_frag_13
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[c][NX3;H0](C)[C]")));
        // [72] abraham_BSEL_frag_14
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[c][nX3;H0][c]")));
        // [73] abraham_BSEL_frag_15
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("*=[Nv3;!R]")));
        // [74] abraham_BSEL_frag_16
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("*=[Nv3;R]")));
        // [75] abraham_BSEL_frag_17
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[nX2H0,nX3H1+](a)a")));
        // [76] abraham_BSEL_frag_18
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("N#C[A;!#1]")));
        // [77] abraham_BSEL_frag_19
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("N#C[a;!#1]")));
        // [78] abraham_BSEL_frag_2
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H1]")));
        // [79] abraham_BSEL_frag_20
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[$([A;!#1][NX3](=O)=O),$([A;!#1][NX3+](=O)[O-])]")));
        // [80] abraham_BSEL_frag_21
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[$([a;!#1][NX3](=O)=O),$([a;!#1][NX3+](=O)[O-])]")));
        // [81] abraham_BSEL_frag_22
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]")));
        // [82] abraham_BSEL_frag_23
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]")));
        // [83] abraham_BSEL_frag_24
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2;H0;!R]")));
        // [84] abraham_BSEL_frag_25
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2;H0;R]")));
        // [85] abraham_BSEL_frag_26
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[oX2](a)a")));
        // [86] abraham_BSEL_frag_27
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("*=O")));
        // [87] abraham_BSEL_frag_28
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[SX2](*)*")));
        // [88] abraham_BSEL_frag_29
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[sX2](a)a")));
        // [89] abraham_BSEL_frag_3
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H0]")));
        // [90] abraham_BSEL_frag_30
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("*=[SX1]")));
        // [91] abraham_BSEL_frag_31
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[SX3]")));
        // [92] abraham_BSEL_frag_32
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[$([#16X4](=[OX1])(=[OX1])([!#8])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([!#8])[OX2H0])]")));
        // [93] abraham_BSEL_frag_33
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[S,s]")));
        // [94] abraham_BSEL_frag_34
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[P,p]")));
        // [95] abraham_BSEL_frag_35
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("FA")));
        // [96] abraham_BSEL_frag_36
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("Fa")));
        // [97] abraham_BSEL_frag_37
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("Cl")));
        // [98] abraham_BSEL_frag_38
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("Br")));
        // [99] abraham_BSEL_frag_39
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("I")));
        // [100] abraham_BSEL_frag_4
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("*=[CX3H2]")));
        // [101] abraham_BSEL_frag_40
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3;!R](=[OX1])[OX2H0]")));
        // [102] abraham_BSEL_frag_41
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3;R](=[OX1])[OX2H0;R]")));
        // [103] abraham_BSEL_frag_42
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("P(=[OX1])(O)(O)O")));
        // [104] abraham_BSEL_frag_43
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=[OX1])([OX2H0])[OX2H0]")));
        // [105] abraham_BSEL_frag_44
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("nC=[OX1]")));
        // [106] abraham_BSEL_frag_45
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[N;!R]C=[OX1]")));
        // [107] abraham_BSEL_frag_46
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[N;R][C;R]=[OX1]")));
        // [108] abraham_BSEL_frag_47
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]")));
        // [109] abraham_BSEL_frag_48
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("NC(=[OX1])N")));
        // [110] abraham_BSEL_frag_49
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]")));
        // [111] abraham_BSEL_frag_5
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[$(*=[CX3H1]),$([cX3H1](a)a)]")));
        // [112] abraham_BSEL_frag_50
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=[OX1])[NX3][CX3](=[OX1])")));
        // [113] abraham_BSEL_frag_51
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("C1(=[OX1])C=CC(=[OX1])C=C1")));
        // [114] abraham_BSEL_frag_52
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[$([CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])]")));
        // [115] abraham_BSEL_frag_53
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]")));
        // [116] abraham_BSEL_frag_54
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("*1~*2~*(~*3~*(~*~*~*~*3)~*1)~*~*~*1~*2~*~*~*1")));
        // [117] abraham_BSEL_frag_55
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2H]CC[O,N]")));
        // [118] abraham_BSEL_frag_56
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2H]C[C,N]=[O,S]")));
        // [119] abraham_BSEL_frag_57
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2H]c1ccccc1[O,Nv3]")));
        // [120] abraham_BSEL_frag_58
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2H]c1ccccc1C=[O,S]")));
        // [121] abraham_BSEL_frag_59
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2H]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]")));
        // [122] abraham_BSEL_frag_6
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[$(*=[CX3H0]),$([cX3H0](a)(a)A)]")));
        // [123] abraham_BSEL_frag_60
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([OH])[CX4][OH]")));
        // [124] abraham_BSEL_frag_61
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("n:n")));
        // [125] abraham_BSEL_frag_62
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("o:n")));
        // [126] abraham_BSEL_frag_63
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("n:c:n")));
        // [127] abraham_BSEL_frag_64
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("o:c:n")));
        // [128] abraham_BSEL_frag_65
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("n:c:c:n")));
        // [129] abraham_BSEL_frag_66
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[F,Cl,Br,I,N,O,S]-c:c-[F,Cl,Br,I,N,O,S]")));
        // [130] abraham_BSEL_frag_67
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[F,Cl,Br,I,N,O,S]-c:c:c-[F,Cl,Br,I,N,O,S]")));
        // [131] abraham_BSEL_frag_68
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[F,Cl,Br,I,N,O,S]-c:c:c:c-[F,Cl,Br,I,N,O,S]")));
        // [132] abraham_BSEL_frag_69
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("P(=[OX1])N")));
        // [133] abraham_BSEL_frag_7
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("c(a)(a)a")));
        // [134] abraham_BSEL_frag_70
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("Nc:n")));
        // [135] abraham_BSEL_frag_71
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[$(cC[OH]);!$(c[CX3](=O)[OX1H0-,OX2H1])]")));
        // [136] abraham_BSEL_frag_72
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]")));
        // [137] abraham_BSEL_frag_73
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2]-c:c-[OX2]")));
        // [138] abraham_BSEL_frag_8
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("*#C")));
        // [139] abraham_BSEL_frag_9
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[c][NX3;H2]")));
        // [140] abraham_S_additional_0
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[CX4]")));
        // [141] abraham_S_additional_1
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)c1ccccc1")));
        // [142] abraham_S_additional_10
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("C#Nc1ccccc1")));
        // [143] abraham_S_additional_11
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3](=O)=Oc1ccccc1")));
        // [144] abraham_S_additional_12
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[SX4](=O)")));
        // [145] abraham_S_additional_13
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[SX4](=O)(=O)")));
        // [146] abraham_S_additional_14
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[F,Cl,Br,I]")));
        // [147] abraham_S_additional_2
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[CX4][CX4]")));
        // [148] abraham_S_additional_3
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[OX2][CX4]")));
        // [149] abraham_S_additional_4
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[OX2]c1ccccc1")));
        // [150] abraham_S_additional_5
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[NX3]")));
        // [151] abraham_S_additional_6
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[NX3]c1ccccc1")));
        // [152] abraham_S_additional_7
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[CX4][CX3](=O)")));
        // [153] abraham_S_additional_8
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[CX3](=O)")));
        // [154] abraham_S_additional_9
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("c1ccccc1[CX3](=O)c1ccccc1")));
        // [155] abraham_new_[#6](:[#6]:[#6]:[#6]-[#6]):[#6]:[#6]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[#6](:[#6]:[#6]:[#6]-[#6]):[#6]:[#6]")));
        // [156] abraham_new_[#6]:[#6]:[#6]:[#6]:[#6]:[#6]-[#6]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[#6]:[#6]:[#6]:[#6]:[#6]:[#6]-[#6]")));
        // [157] abraham_new_[#6]:[#6]:[#6]:[#6]:[#6]:[#6]-[#8]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[#6]:[#6]:[#6]:[#6]:[#6]:[#6]-[#8]")));
        // [158] abraham_new_[#6]:[#6]:[#6]:[#6]:[#6]:[#6].[#6]-[#6]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[#6]:[#6]:[#6]:[#6]:[#6]:[#6].[#6]-[#6]")));
        // [159] abraham_new_[#6]:[#6]:[#6]:[#6]:[#6]:[#6].[#6]:[#6]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[#6]:[#6]:[#6]:[#6]:[#6]:[#6].[#6]:[#6]")));
        // [160] abraham_new_[#6]:[#6]:[#6]:[#6]:[#6]:[#6]:[#6]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[#6]:[#6]:[#6]:[#6]:[#6]:[#6]:[#6]")));
        // [161] abraham_new_[CX2H0]#[CX2H0]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX2H0]#[CX2H0]")));
        // [162] abraham_new_[CX2H0]#[CX2H]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX2H0]#[CX2H]")));
        // [163] abraham_new_[CX2H0]([CX4H3])#[CX2H]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX2H0]([CX4H3])#[CX2H]")));
        // [164] abraham_new_[CX2H0]([CX4H3])([CX4H3])#[CX2H]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX2H0]([CX4H3])([CX4H3])#[CX2H]")));
        // [165] abraham_new_[CX2H]#[CX2H]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX2H]#[CX2H]")));
        // [166] abraham_new_[CX2H]([CX4H3])#[CX2H]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX2H]([CX4H3])#[CX2H]")));
        // [167] abraham_new_[CX3H0]([CX4H3])([CX4H3])=[CX3H1]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3H0]([CX4H3])([CX4H3])=[CX3H1]")));
        // [168] abraham_new_[CX3H0]([CX4H3])([CX4H3])=[CX3H2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3H0]([CX4H3])([CX4H3])=[CX3H2]")));
        // [169] abraham_new_[CX3H0]([CX4H3])=[CX3H1]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3H0]([CX4H3])=[CX3H1]")));
        // [170] abraham_new_[CX3H0]([CX4H3])=[CX3H2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3H0]([CX4H3])=[CX3H2]")));
        // [171] abraham_new_[CX3H0]=[CX3H0]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3H0]=[CX3H0]")));
        // [172] abraham_new_[CX3H0]=[CX3H1]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3H0]=[CX3H1]")));
        // [173] abraham_new_[CX3H0]=[CX3H2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3H0]=[CX3H2]")));
        // [174] abraham_new_[CX3H1]([CX4H3])=[CX3H1]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3H1]([CX4H3])=[CX3H1]")));
        // [175] abraham_new_[CX3H1]([CX4H3])=[CX3H2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3H1]([CX4H3])=[CX3H2]")));
        // [176] abraham_new_[CX3H1]=[CX3H1]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3H1]=[CX3H1]")));
        // [177] abraham_new_[CX3H1]=[CX3H2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3H1]=[CX3H2]")));
        // [178] abraham_new_[CX3H2]=[CX3H2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3H2]=[CX3H2]")));
        // [179] abraham_new_[CX3](=O)
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)")));
        // [180] abraham_new_[CX3](=O)[NX3][CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[NX3][CX4]")));
        // [181] abraham_new_[CX3](=O)[OX2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3](=O)[OX2]")));
        // [182] abraham_new_[CX3]=[CX3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3]=[CX3]")));
        // [183] abraham_new_[CX3]=[CX3][CX4]([CX4])
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3]=[CX3][CX4]([CX4])")));
        // [184] abraham_new_[CX3]=[CX3][CX4]([CX4])([CX4])
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3]=[CX3][CX4]([CX4])([CX4])")));
        // [185] abraham_new_[CX3]=[CX3][CX4]([CX4])[CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX3]=[CX3][CX4]([CX4])[CX4]")));
        // [186] abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]")));
        // [187] abraham_new_[CX4H0]([CX4H3])([CX4H2])[CX4H2][CX4H3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H0]([CX4H3])([CX4H2])[CX4H2][CX4H3]")));
        // [188] abraham_new_[CX4H0]([CX4H3])([CX4H3])([CX4H3])([CX4H3])
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H0]([CX4H3])([CX4H3])([CX4H3])([CX4H3])")));
        // [189] abraham_new_[CX4H0]([CX4H3])([CX4H3])([CX4H3])[CX4H3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H0]([CX4H3])([CX4H3])([CX4H3])[CX4H3]")));
        // [190] abraham_new_[CX4H0]([CX4H3])([CX4H3])[CX4H2][CX4H3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H0]([CX4H3])([CX4H3])[CX4H2][CX4H3]")));
        // [191] abraham_new_[CX4H0]([CX4H3])[CX2H]#[CX2H]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H0]([CX4H3])[CX2H]#[CX2H]")));
        // [192] abraham_new_[CX4H0]([CX4H3])[CX3H1]=[CX3H2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H0]([CX4H3])[CX3H1]=[CX3H2]")));
        // [193] abraham_new_[CX4H1]([CX4H3])([CX4H2])[CX4H3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H1]([CX4H3])([CX4H2])[CX4H3]")));
        // [194] abraham_new_[CX4H1]([CX4H3])([CX4H3])([CX4H3])
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H1]([CX4H3])([CX4H3])([CX4H3])")));
        // [195] abraham_new_[CX4H1]([CX4H3])([CX4H3])[CX4H3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H1]([CX4H3])([CX4H3])[CX4H3]")));
        // [196] abraham_new_[CX4H1]([CX4H3])[CX2H]#[CX2H]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H1]([CX4H3])[CX2H]#[CX2H]")));
        // [197] abraham_new_[CX4H1]([CX4H3])[CX3H1]=[CX3H2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H1]([CX4H3])[CX3H1]=[CX3H2]")));
        // [198] abraham_new_[CX4H1]([CX4H3])[CX4H2][CX4H2][CX4H3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H1]([CX4H3])[CX4H2][CX4H2][CX4H3]")));
        // [199] abraham_new_[CX4H1]([CX4H3])[CX4H2][CX4H3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H1]([CX4H3])[CX4H2][CX4H3]")));
        // [200] abraham_new_[CX4H2]([CX4H3])([CX4H3])
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H2]([CX4H3])([CX4H3])")));
        // [201] abraham_new_[CX4H2]([CX4H3])[CX2H]#[CX2H]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H2]([CX4H3])[CX2H]#[CX2H]")));
        // [202] abraham_new_[CX4H2]([CX4H3])[CX3H1]=[CX3H2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H2]([CX4H3])[CX3H1]=[CX3H2]")));
        // [203] abraham_new_[CX4H2]([CX4H3])[CX4H2][CX4H2][CX4H3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H2]([CX4H3])[CX4H2][CX4H2][CX4H3]")));
        // [204] abraham_new_[CX4H2]([CX4H3])[CX4H2][CX4H3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H2]([CX4H3])[CX4H2][CX4H3]")));
        // [205] abraham_new_[CX4H2]([CX4H3])[CX4H3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H2]([CX4H3])[CX4H3]")));
        // [206] abraham_new_[CX4H3][CX2H]#[CX2H]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H3][CX2H]#[CX2H]")));
        // [207] abraham_new_[CX4H3][CX3H0]=[CX3H2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H3][CX3H0]=[CX3H2]")));
        // [208] abraham_new_[CX4H3][CX3H1]=[CX3H2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H3][CX3H1]=[CX3H2]")));
        // [209] abraham_new_[CX4H3][CX4H2][CX4H3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H3][CX4H2][CX4H3]")));
        // [210] abraham_new_[CX4H][CX4]([CX4])[CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4H][CX4]([CX4])[CX4]")));
        // [211] abraham_new_[CX4]([CX4])([CX4])([CX4])
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([CX4])([CX4])([CX4])")));
        // [212] abraham_new_[CX4]([CX4])([CX4])[CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([CX4])([CX4])[CX4]")));
        // [213] abraham_new_[CX4]([CX4])[CX3]=[CX3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([CX4])[CX3]=[CX3]")));
        // [214] abraham_new_[CX4]([CX4])[CX3]=[CX3][CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([CX4])[CX3]=[CX3][CX4]")));
        // [215] abraham_new_[CX4]([CX4])[CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([CX4])[CX4]")));
        // [216] abraham_new_[CX4]([CX4])[CX4]([CX4])[CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([CX4])[CX4]([CX4])[CX4]")));
        // [217] abraham_new_[CX4]([CX4])[CX4][CX3](=O)
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([CX4])[CX4][CX3](=O)")));
        // [218] abraham_new_[CX4]([CX4])[CX4][CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([CX4])[CX4][CX4]")));
        // [219] abraham_new_[CX4]([CX4])[CX4][OH]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[CX4]([CX4])[CX4][OH]")));
        // [220] abraham_new_[F,Cl,Br,I]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[F,Cl,Br,I]")));
        // [221] abraham_new_[F,Cl,Br,I][CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[F,Cl,Br,I][CX4]")));
        // [222] abraham_new_[F,Cl,Br,I]c1ccccc1
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[F,Cl,Br,I]c1ccccc1")));
        // [223] abraham_new_[NX3;H0]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3;H0]")));
        // [224] abraham_new_[NX3;H1]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3;H1]")));
        // [225] abraham_new_[NX3;H2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3;H2]")));
        // [226] abraham_new_[NX3][CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3][CX4]")));
        // [227] abraham_new_[NX3]c1ccccc1
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[NX3]c1ccccc1")));
        // [228] abraham_new_[OH][CX3](=O)
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH][CX3](=O)")));
        // [229] abraham_new_[OH][CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH][CX4]")));
        // [230] abraham_new_[OH][CX4][CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH][CX4][CX4]")));
        // [231] abraham_new_[OH]c1ccccc1
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OH]c1ccccc1")));
        // [232] abraham_new_[OX2][CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2][CX4]")));
        // [233] abraham_new_[OX2][CX4][CX4]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2][CX4][CX4]")));
        // [234] abraham_new_[OX2]c1ccccc1
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("[OX2]c1ccccc1")));
        // [235] abraham_new_c1ccccc1
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("c1ccccc1")));
        // [236] abraham_new_c1ccccc1[CX3](=O)
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("c1ccccc1[CX3](=O)")));
        // [237] abraham_new_c1ccccc1[F,Cl,Br,I]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("c1ccccc1[F,Cl,Br,I]")));
        // [238] abraham_new_c1ccccc1[NX3]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("c1ccccc1[NX3]")));
        // [239] abraham_new_c1ccccc1[OH]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("c1ccccc1[OH]")));
        // [240] abraham_new_c1ccccc1[OX2]
        queries.push_back(std::shared_ptr<RWMol>(SmartsToMol("c1ccccc1[OX2]")));
    });
    
    return queries;
}

} // namespace Osmordred
} // namespace Descriptors
} // namespace RDKit
