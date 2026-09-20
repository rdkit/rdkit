# InformationContent: which criterion, and why

`CalcInformationContent` computes Basak's neighbourhood-complexity indices —
`IC`, `TIC`, `SIC`, `BIC`, `CIC`, `MIC`, `ZMIC`. Two different equivalence
criteria are available through `ICKeyFlavor`, and they do not agree. This note
records what each one is, and the evidence for the default.

## The primary sources

**Roy, A. B., Basak, S. C., Harriss, D. K. & Magnuson, V. R. (1983).**
"Neighborhood complexities and symmetry of chemical graphs and their biological
applications." In *Mathematical Modelling in Science and Technology* (Avula,
Kalman, Liapis & Rodin, eds.), pp. 745–750. Pergamon Press.
— Defines the indices and gives the worked example used below.

**Basak, S. C. (1987).** "Use of molecular complexity indices in predictive
pharmacology and toxicology: a QSAR approach." *Medical Science Research* **15**,
605–609.
— Review; restates the definitions and the 2-butenol neighbourhoods.

**Majumdar, S., Basak, S. C. et al. (2019).** "Finding Needles in a Haystack:
Determining Key Molecular Descriptors..." *Molecular Informatics*,
`minf201800164`.
— Its supplementary data is the output of **POLLY**, the program Basak's group
wrote and used. 413 molecules with `IC0..IC5`. That file is the oracle this
implementation is validated against; a copy lives in
`test_data/basak_polly_ic_413.csv`.

The papers themselves are not redistributed here — only the citations and the
supplementary data table.

Definitions, all from the 1983 paper, where `A` is the atom count of the
hydrogen-filled graph and `B` the bond count:

```
IC_r  = -sum_i (n_i/A) log2(n_i/A)     Shannon entropy of the class sizes
TIC_r = A * IC_r
SIC_r = IC_r / log2(A)
CIC_r = log2(A) - IC_r
BIC_r = IC_r / log2(B)
```

Only `IC_r` carries information; the other four are algebra. A single wrong
`IC_r` therefore shows up as five reported failures.

## The two criteria, side by side

Both partition the atoms of the hydrogen-filled graph into classes by their
radius-`r` environment, then take the Shannon entropy. They differ in what
"environment" means.

|  | `BASAK` (default) | `MORDRED` |
|---|---|---|
| preprocessing | none | **kekulize**, then add H |
| order 0 | group by atomic number | group by atomic number |
| order ≥ 1 key | flat sorted multiset of **one-step edge descriptors** over the newly reached shell — `(rootNum, rootDeg, bondOrder, neighNum)` | sorted multiset of **complete root-to-leaf path codes** of a BFS tree — every vertex's `(atomicNum, degree)` and every step's bond type |
| radius growth | frontier extended in place, `M`/`SP` carried forward | tree expanded one shell further |
| cost | ~9,700 key integers on the reference set | ~158,400 — about **16×**, growing superlinearly with radius |

`EXTENDED` is a third value: the `BASAK` key with the neighbour degree folded in
as well. It was added "for Mordred parity" and is kept only so that older
osmordred numbers can be reproduced.

### On the radius loop

Both flavours here treat the radius as an **evolution chain**: the structure at
order `r` is grown from the structure at order `r-1`. The mordred package does
not — its `get_code(i, order)` calls `reset(i)` and re-expands from scratch, so
computing orders 0..5 costs `0+1+2+3+4+5 = 15` expansions per atom where 5
suffice. This implementation grows the tree once and samples it after each
expansion. The values are identical; the rebuild is not.

One subtlety: the `visited` set accumulates across expansions and must **not** be
cleared per order. Clearing it lets an already-claimed atom be reached twice.

## Mordred decides atom identity at the wrong level

This is the substantive disagreement, not an implementation detail.

Basak is explicit about what a vertex *is*. From the 1983 paper, p. 746:

> a chemical graph associated with a molecule is denoted by `G(V,E)` where
> `V = {v1, v2 ... vm}` denotes the collection of atoms of the form
> `vi = (ai, di)`, where `ai` stands for an atom's chemical identity of `vi`,
> and `di` (a positive integer) denotes the **valency** of `vi`.

Valency is a *chemical* property — how many bonds the atom forms, counting
multiplicity. Graph degree is a *topological* one — how many neighbours it has.
In a hydrogen-filled graph they coincide for saturated atoms and diverge exactly
at double, triple and aromatic bonds.

Mordred uses the degree. That is a different claim about chemistry: it says an
sp² carbon (degree 3) is a different kind of atom from an sp³ carbon (degree 4),
where Basak says both are tetravalent carbon and belong together.

Basak's own worked example is decided by this. 2-butenol, 1983 Table 1, order 1:

```
    H2  H4  H5  H6
     |   |   |   |
H1 - C1- C2= C3 -C4 - O - H8
     |           |
     H3          H7
```

Every carbon has valency 4 — C2 and C3 included, the double bond counting
twice — so every carbon-bound hydrogen has the same first-order environment and
the seven of them form **one class**. Only H8, bound to a divalent oxygen, is
separate. The paper's partition is `[1,1,1,1,2,7]`, giving `IC1 = 2.0349`.

Reading the same molecule by degree, C1 and C4 have degree 4 while C2 and C3
have degree 3, so the seven hydrogens split 5 + 2. That gives `[1,1,1,1,2,2,5]`
and `IC1 = 2.4997` — which is exactly what mordred returns, and is not what
Basak published.

`BASAK` reproduces the published partition and the published value. `MORDRED`
reproduces mordred.

One honest qualification about *how* `BASAK` gets there. Its key does not encode
the neighbour's degree at all, so every carbon-bound hydrogen keys the same
regardless of which carbon — which lands on Basak's grouping without computing
valency explicitly. It agrees with Basak's answer here; it is not a literal
implementation of the `(element, valency)` vertex.

`ICVertexLabel.VALENCY` makes the literal reading available, and the two are
near-indistinguishable in practice: 87.8% against POLLY versus 88.0% for
`DEGREE`. The level at which the *neighbour* is judged is what matters, and that
is where mordred diverges — not the label on the root.

## Mordred is not POLLY

This is the reason `BASAK` is the default. Measured against POLLY on 411
molecules × `r=0..5` = 2466 values, tolerance 0.0015:

| | IC0 | IC1 | IC2 | IC3 | IC4 | IC5 | overall |
|---|---|---|---|---|---|---|---|
| `BASAK` | 93.7% | **80.3%** | 83.5% | 87.3% | 90.8% | 92.5% | **88.0%** |
| `MORDRED` | 93.7% | **14.4%** | 32.4% | 44.3% | 55.7% | 58.2% | **49.8%** |

Mordred agrees with Basak's own software on about half of all values, and on
one in seven at order 1. Concrete cases — `BASAK` reproduces POLLY exactly while
`MORDRED` does not:

| POLLY | `BASAK` | `MORDRED` | molecule |
|---:|---:|---:|---|
| 2.903 | **2.903** | 3.290 | `NCCc1cn2c(n1)cccc2` |
| 2.724 | **2.724** | 3.193 | `CCCN(CCC)CCc1ccc(c2c1CC(N2)=C)O` |
| 3.296 | **3.296** | 3.624 | `c1(c2c(cc(F)cc2)on1)C1CCN(CCc2c(n3c(...` |
| 2.898 | **2.898** | 3.255 | `c12c(C(c3ccccc3)=NCc3n1c(nn3)C)cc(Cl)cc2` |
| 2.214 | **2.214** | 2.747 | `C1(\c2c(CCc3c1cccc3)cccc2)=C\CCN(C)C` |

(IC1; POLLY values as published.)

The other published worked values — n-propanol `IC2 = 2.855` with partition
`[1,1,1,1,1,2,2,3]`, and the ten aliphatic alcohols of the 1983 Table 2 — are
reproduced by both, so 2-butenol is the only published case that discriminates.
The POLLY table is what makes the difference measurable at scale.

## What `MORDRED` is for

Reproducing the mordred package, not speed. It matches it on 894/894 values and
passes all 344 entries of `test_data/mordred_references/InformationContent.yaml`
— which is what that file describes, now that the default is Basak's criterion.

## Known limits of the default

`BASAK` reaches 88.0%, not 100%. The residual is concentrated where RDKit's
graph breaks a symmetry the molecule actually has:

| | agreement |
|---|---:|
| plain molecules | 92.3% |
| tautomer-ambiguous | 88.2% |
| resonance-asymmetric (nitro, carboxylate, sulfonate) | 33.3% |

RDKit writes nitro as `[N+](=O)[O-]`, so its two chemically equivalent oxygens
take different canonical ranks and different keys; every nitro molecule in the
reference set over-splits. `equalizeDelocalizedBonds` is a partial mitigation
(33.3% → 40.3%), off by default because it does not close the gap and it changes
the descriptor's meaning. Note that every nitro molecule in the reference set is
also aromatic, so the two effects cannot be separated on the data available.
