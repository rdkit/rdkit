# New CIP labelling

This is a C++ port of https://github.com/SiMolecule/centres, which was
originally written by John Mayfied in Java. The port is originally written
ontop of RDKit, but maintains the abstraction layer present in the orignal
version, so that it should be possible to write adapt it to other frameworks
by writing appropriate Molecule and Labeller classes, and some tests (the ones
included are RDKit-specific).

Currently, this is in a quite un-optimized state, as I wanted to start with
code as close as possible to the original until some testing is available.

### A couple of remarks:
- This does not perceive potential stereochemical configurations in the
molecule: the input already needs to have information on which atoms/bonds
should be labeled. Also, for resolution of tetrahedral chirality, the parity
of the atoms is required, and for double bond stereochemistry, the bond must
already be labeled as cis or trans, and both "anchors" must be known.

- This does not provide CIP rank information on, e.g., neighbors around a
chiral center.

- The input molecule does not require all Hydrogens to be explicit, the code
will take them into account appropriately, although the queries for atomic
and mass numbers of those will be handed down to the molecule class by passing
a nullptr atom to the appropriate methods.

### Some potential points of improvement might be:
- Implement the stereo configuration modes that haven't been ported yet
(atropoisomery, Square planar, extended tetrahedral & cis/trans...). These
haven't been ported yed because RDKit cannot perceive these (or does it?).

- Handle atomic and mass numbers of implicit Hydrogens in the BaseMol class
instead of handing them down to the framework molecule class.

- Improve/remove the templates system, so that it isn't necessary to implement
everything in headers.

- Digraph generation/handling. It might be useful to improve the way nodes,
edges and the digraph itself are handled, so that they don't need to be passed
around as pointers. It might also be useful to port this to use boost::graph.

- Handling of rules & sorting -- also to avoid pointers that might leak.

- Improve handling of vectors of nodes/edges. This was a quick way of porting
over Java lists and arrays, but it might be improved (e.g. by using std::array)
in some situations where the final/expected sizes are known.

- Check for bugs and add (much) more testing, as the current tests barely
scratch the surface.

