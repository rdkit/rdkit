# $Id$
#
# Copyright (C) 2001-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" periodic table data, **obsolete**

  now that the C++ code exposes an interface to the internal PT stuff,
  this data is mostly obsolete

"""
import re
blankExpr = re.compile(r'\ *\t*\ *')
# Num	Symb	RCov	RBO	RVdW	Max Bnd	Mass   nval
periodicData=\
"""
0	X	0.0	0.0	0.0	0	0.000	0  
1	H	0.230	0.330	1.200	1	1.008	1  
2	He	0.930	0.700	1.400	0	4.003	2  
3	Li	0.680	1.230	1.820	1	6.941	1  
4	Be	0.350	0.900	1.700	2	9.012	2  
5	B	0.830	0.820	2.080	3	10.812	3  
6	C	0.680	0.770	1.950	4	12.011	4  
7	N	0.680	0.700	1.850	4	14.007	5  
8	O	0.680	0.660	1.700	2	15.999	6  
9	F	0.640	0.611	1.730	1	18.998	7  
10	Ne	1.120	0.700	1.540	0	20.180	 8 
11	Na	0.970	1.540	2.270	1	22.990	 1 
12	Mg	1.100	1.360	1.730	2	24.305	 2 
13	Al	1.350	1.180	2.050	6	26.982	 3 
14	Si	1.200	0.937	2.100	6	28.086	 4 
15	P	0.750	0.890	2.080	5	30.974	 5 
16	S	1.020	1.040	2.000	6	32.067	 6 
17	Cl	0.990	0.997	1.970	1	35.453	 7 
18	Ar	1.570	1.740	1.880	0	39.948	 8 
19	K	1.330	2.030	2.750	1	39.098	 1 
20	Ca	0.990	1.740	1.973	2	40.078	 2 
21	Sc	1.440	1.440	1.700	6	44.956	 3 
22	Ti	1.470	1.320	1.700	6	47.867	 4 
23	V	1.330	1.220	1.700	6	50.942	 5 
24	Cr	1.350	1.180	1.700	6	51.996	 6 
25	Mn	1.350	1.170	1.700	8	54.938	 7 
26	Fe	1.340	1.170	1.700	6	55.845	 8 
27	Co	1.330	1.160	1.700	6	58.933	 9 
28	Ni	1.500	1.150	1.630	6	58.693	 10
29	Cu	1.520	1.170	1.400	6	63.546	 11
30	Zn	1.450	1.250	1.390	6	65.39	 2 
31	Ga	1.220	1.260	1.870	3	69.723	 3 
32	Ge	1.170	1.188	1.700	4	72.61	 4 
33	As	1.210	1.200	1.850	3	74.922	 5 
34	Se	1.220	1.170	1.900	2	78.96	 6 
35	Br	1.210	1.167	2.100	1	79.904	 7 
36	Kr	1.910	1.910	2.020	0	83.80	 8 
37	Rb	1.470	2.160	1.700	1	85.468	 1 
38	Sr	1.120	1.910	1.700	2	87.62	 2 
39	Y	1.780	1.620	1.700	6	88.906	 3 
40	Zr	1.560	1.450	1.700	6	91.224	 4 
41	Nb	1.480	1.340	1.700	6	92.906	 5 
42	Mo	1.470	1.300	1.700	6	95.94	 6 
43	Tc	1.350	1.270	1.700	6	98.0	 7 
44	Ru	1.400	1.250	1.700	6	101.07	 8 
45	Rh	1.450	1.250	1.700	6	102.906	 9 
46	Pd	1.500	1.280	1.630	6	106.42	 10
47	Ag	1.590	1.340	1.720	6	107.868	 11
48	Cd	1.690	1.480	1.580	6	112.412	 2 
49	In	1.630	1.440	1.930	3	114.818	 3 
50	Sn	1.460	1.385	2.170	4	118.711	 4 
51	Sb	1.460	1.400	2.200	3	121.760	 5 
52	Te	1.470	1.378	2.060	2	127.60	 6 
53	I	1.400	1.387	2.150	1	126.904	 7 
54	Xe	1.980	1.980	2.160	0	131.29	 8 
55	Cs	1.670	2.350	1.700	1	132.905	 1 
56	Ba	1.340	1.980	1.700	2	137.328	 2 
57	La	1.870	1.690	1.700	12	138.906	 3 
58	Ce	1.830	1.830	1.700	6	140.116	 4 
59	Pr	1.820	1.820	1.700	6	140.908	 3 
60	Nd	1.810	1.810	1.700	6	144.24	 4 
61	Pm	1.800	1.800	1.700	6	145.0	 5 
62	Sm	1.800	1.800	1.700	6	150.36	 6 
63	Eu	1.990	1.990	1.700	6	151.964	 7 
64	Gd	1.790	1.790	1.700	6	157.25	 8 
65	Tb	1.760	1.760	1.700	6	158.925	 9 
66	Dy	1.750	1.750	1.700	6	162.50	 10
67	Ho	1.740	1.740	1.700	6	164.930	 11
68	Er	1.730	1.730	1.700	6	167.26	 12
69	Tm	1.720	1.720	1.700	6	168.934	 13
70	Yb	1.940	1.940	1.700	6	173.04	 14
71	Lu	1.720	1.720	1.700	6	174.967	 15
72	Hf	1.570	1.440	1.700	6	178.49	 4 
73	Ta	1.430	1.340	1.700	6	180.948	 5 
74	W	1.370	1.300	1.700	6	183.84	 6 
75	Re	1.350	1.280	1.700	6	186.207	 7 
76	Os	1.370	1.260	1.700	6	190.23	 8 
77	Ir	1.320	1.270	1.700	6	192.217	 9 
78	Pt	1.500	1.300	1.720	6	195.078	 10
79	Au	1.500	1.340	1.660	6	196.967	 11
80	Hg	1.700	1.490	1.550	6	200.59	 2 
81	Tl	1.550	1.480	1.960	3	204.383	 3 
82	Pb	1.540	1.480	2.020	4	207.2	 4 
83	Bi	1.540	1.450	1.700	3	208.980	 5 
84	Po	1.680	1.460	1.700	2	209.0	 6 
85	At	1.700	1.450	1.700	1	210.0	 7 
86	Rn	2.400	2.400	1.700	0	222.0	 8 
87	Fr	2.000	2.000	1.700	1	223.0	 1 
88	Ra	1.900	1.900	1.700	2	226.0	 2 
89	Ac	1.880	1.880	1.700	6	227.0	 3 
90	Th	1.790	1.790	1.700	6	232.038	 4 
91	Pa	1.610	1.610	1.700	6	231.036	 3 
92	U	1.580	1.580	1.860	6	238.029	 4 
93	Np	1.550	1.550	1.700	6	237.0	 5 
94	Pu	1.530	1.530	1.700	6	244.0	 6 
95	Am	1.510	1.070	1.700	6	243.0	 7 
96	Cm	1.500	0.000	1.700	6	247.0	 8 
97	Bk	1.500	0.000	1.700	6	247.0	 9 
98	Cf	1.500	0.000	1.700	6	251.0	 10
99	Es	1.500	0.000	1.700	6	252.0	 11
100	Fm	1.500	0.000	1.700	6	257.0	 12
101	Md	1.500	0.000	1.700	6	258.0	 13
102	No	1.500	0.000	1.700	6	259.0	 14
103	Lr	1.500	0.000	1.700	6	262.0	 15
"""

nameTable = {}
numTable = {}
for line in periodicData.split('\n'):
  splitLine = blankExpr.split(line)
  if len(splitLine)>1:
    nameTable[splitLine[1]] = (int(splitLine[0]),float(splitLine[6]),int(splitLine[7]),\
                               int(splitLine[5]),float(splitLine[2]),float(splitLine[3]),
                               float(splitLine[4]))
    numTable[int(splitLine[0])] = (splitLine[1],float(splitLine[6]),int(splitLine[7]),\
                                   int(splitLine[5]),float(splitLine[2]),float(splitLine[3]),
                                   float(splitLine[4]))

# a list of metals (transition metals, semi-metals, lanthanides and actinides)
metalRanges = ["13","21-32","39-51","57-84","89-103"]
metalNumList = []
for entry in metalRanges:
  t = entry.split('-')
  start = int(t[0])
  if len(t)>1:
    end = int(t[1])
  else:
    end = start
  if start > end:
    start,end = end,start
  metalNumList += range(start,end+1)  
metalNames = map(lambda x:numTable[x][0],metalNumList)

# these are from table 4 of Rev. Comp. Chem. vol 2, 367-422, (1991)
#  the order is [alpha(SP),alpha(SP2),alpha(SP3)]
# where values are not known, None has been inserted
hallKierAlphas = {
  'H':[0.0,0.0,0.0], # removes explicit H's from consideration in the shape
  'C':[-0.22,-0.13,0.0],
  'N':[-0.29,-0.20,-0.04],
  'O':[None,-0.20,-0.04],
  'F':[None,None,-0.07],
  'P':[None,0.30,0.43],
  'S':[None,0.22,0.35],
  'Cl':[None,None,0.29],
  'Br':[None,None,0.48],
  'I':[None,None,0.73]}

