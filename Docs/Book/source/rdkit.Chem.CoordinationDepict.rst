rdkit.Chem.CoordinationDepict module
====================================

The standard 2D coordinate generators are optimized for organic molecular
graphs. The helpers in this module depict each ligand independently with
CoordGen and arrange the resulting fragments around a single metal centre. The
layout approach was informed by the `metal2d project
<https://github.com/levakrasnovs/metal2d>`_; this module is a self-contained
implementation using existing RDKit APIs.

The coordinate generator modifies the molecule in place and returns the new
conformer id::

    from rdkit import Chem
    from rdkit.Chem import CoordinationDepict

    mol = Chem.MolFromSmiles("[NH3]->[Pt+2](<-[NH3])(<-[Cl-])<-[Cl-]")
    CoordinationDepict.Compute2DCoordinationCoords(mol)

For eta-bound ligands, the drawing helper uses RDKit's haptic-bond
representation and leaves the input molecule unchanged::

    drawing_mol = CoordinationDepict.PrepareCoordinationMolForDrawing(mol)

The initial implementation handles mononuclear complexes. Other molecules
fall back to regular CoordGen coordinate generation. The requested metal-bond
length is an initial value; the layout may increase all metal-ligand distances
uniformly when that is necessary to prevent clashes. Calling the coordinate
generator requires an RDKit build with CoordGen support.

.. automodule:: rdkit.Chem.CoordinationDepict
    :members:
    :undoc-members:
    :show-inheritance:
