# Taken from Taylor et al., JCAMD, 26, 451-472 (2012), supporting information.
# 22 Donor SMARTS as 6th in paper is expanded to 2 because there's an
# error in the SMARTS as given - H is used as an explicit atom
# Similarly, the 17th donor in the paper (18th in this file) is incorrect
# SMARTS for much the same reason
[*;+;!H0]     donor1
[NX4;!H0]     donor2
[NX3;!H0]C(=[NX3;!H0])[NX3;!H0]     donor3
[nX3;!H0]c([nX3;!H0])[nX3;!H0]     donor4
[NX3;!H0]C=[NX3;!H0]     donor5
[nX3;!H0]c[nX3;!H0]     donor6
[B,C,N,P,S]=[NX3;!H0]-[B,c,C,n,N,O,F,P,S,Cl,Br,I]     donor7
[B,C,N,P,S]=[NX3;!H0]     donor8
a[nX3;!H0]a     donor9
[O,S]=C[NX3;!H0]     donor10
a-[NX3;!H0]     donor11
[C,N,O,P,S]=[C,N,P,S]-[NX3;!H0]     donor12
[C,N]#C-[NX3;!H0]     donor13
[NX3;!H0]     donor14
[NX2;!H0]     donor15
[OX2;!H0]     donor16
[SX2;!H0]     donor17
[NX3]=[CX3;!H0]     donor18
[nX3]-[CH]     donor19
[nX3][cX3;!H0]     donor20
N=[CX3;!H0]N     donor21
n[cX3;!H0]n     donor22
[CX2]#[CX2;!H0]     donor23
[*-]     acceptor1
[OX1][CX3]=[OX1]     acceptor2
[oX1]~[cX3]~[oX1]     acceptor3
[OX1][PX4](=O)[OX1]     acceptor4
[oX1][pX4](~o)~o     acceptor5
[OX1][PX4](=O)     acceptor6
[oX1]~[pX4]~o     acceptor7
[OX1][SX4](=O)=O     acceptor8
[oX1][sX4](~o)~o     acceptor9
c-[OX1]     acceptor10
[CX4]-[OX1]     acceptor11
[CX3]-[OX1]     acceptor12
[OX2]     acceptor13
[#6]~[NX3](~[#6])[OX1]     acceptor14
# 2nd page of acceptors
[#6][NX4]([#6])([#6])[OX1]     acceptor15
[OX2][SX4]([OX2])(=O)(=O)     acceptor16
O=[SX4](=O)[NX2][B,#6,#7,O,F,P,S,Cl,Br,I]     acceptor17
# extra one, bad H
O=[SX4](=O)[NX2;H]     acceptor18
O=[SX4](=O)[NX3]     acceptor19
[OX1]=S([#6])[#6]     acceptor20
[OX1]S([#6])[#6]     acceptor21
[OX1]     acceptor22
[CX4][NX3;H2]     acceptor23
[NX3][NX3;H2]     acceptor24
[CX4][NX3;H1][CX4]     acceptor25
[CX4][NX3;H1][NX3]     acceptor26
[NX3][NX3;H1][NX3]     acceptor27
[CX4][NX3;H1]O     acceptor28
[CX4][NX3;H1][SX2]     acceptor29
[NX3]([CX4])([CX4])[CX4]     acceptor30
[NX3]([CX4])([CX4])[NX3]     acceptor31
[NX3]([CX4])([NX3])[NX3]     acceptor32
[NX3]([CX4])([CX4])[OX2]     acceptor33
[NX3]([CX4])([CX4])[SX2]     acceptor34
[NX3]([CX4])([CX4])[F,Cl,Br,I]     acceptor35
[nX2]     acceptor36
[B,c,C,n,N,O,F,P,S,Cl,Br,I]=[NX2][B,c,C,n,N,O,F,P,S,Cl,Br,I]     acceptor37
# extra one, bad H
[B,c,C,n,N,O,F,P,S,Cl,Br,I]=[NX2;H]     acceptor38
[B,c,C,n,N,O,F,P,S,Cl,Br,I][NX2][B,c,C,n,N,O,F,P,S,Cl,Br,I]     acceptor39
# extra one, bad H
[B,c,C,n,N,O,F,P,S,Cl,Br,I][NX2;H]     acceptor40
[CX2]#[NX1]     acceptor41
c-[SX1]     acceptor42
[CX4]-[SX1]     acceptor43
[CX3]-[SX1]     acceptor44
[NX3]C(=S)[NX3]     acceptor45
# ones marked IGNORE from the paper
# there's an extra one owing to another erroneous H
[#6,#7]=[#6,#7]-[OR]-[#6,#7]=[#6,#7]     not_acceptor1
[c,n]o[c,n]     not_acceptor2
[#6,#7]=[#6,#7]-[OR]-[c,n]     not_acceptor3
[#6,N]-O[CX3]=[O,S]     not_acceptor4
# extra one
[OH][CX3]=[O,S]     not_acceptor5
c-O-[CH3]     not_acceptor6
C=C-O-[CH3]     not_acceptor7
