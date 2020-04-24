# New CIP labelling

This is a C++ port of https://github.com/SiMolecule/centres, which was
originally written by John Mayfied in Java. The orignal algorithm was
described in:

Hanson, R. M., Musacchio, S., Mayfield, J. W., Vainio, M. J., Yerin, A., Redkin, D. Algorithmic Analysis of Cahn−Ingold−Prelog Rules of Stereochemistry: Proposals for Revised Rules and a Guide for Machine Implementation. J. Chem. Inf. Model. 2018, 58, 1755-1765.


### A couple of remarks:
- This does not perceive potential stereochemical configurations in the
molecule: the input already needs to have information on which atoms/bonds
should be labeled. Also, for resolution of tetrahedral chirality, the parity
of the atoms is required, and for double bond stereochemistry, the bond must
have the appropriate bond directions set on the neighboring bonds.

- This does not provide CIP rank information on, e.g., neighbors around a
chiral center.

- The input molecule does not require all Hydrogens to be explicit, the code
will take them into account appropriately.

- Stereochemical labels will be applied by setting the "_CIPCode" property
on the relevant atoms and bonds.

### Some potential points of improvement might be:
- Implement the stereo configuration modes that haven't been ported yet
(atropoisomery, Square planar, extended tetrahedral & cis/trans...). These
haven't been ported yed because RDKit cannot perceive these.

- Digraph generation/handling. It might be useful to improve the way nodes,
edges and the digraph itself are handled, so that they don't need to be passed
around as pointers. It might also be useful to port this to use boost::graph.

- Handling of rules & sorting -- also to avoid pointers that might leak.

- Improve handling of vectors of nodes/edges. This was a quick way of porting
over Java lists and arrays, but it might be improved (e.g. by using std::array)
in some situations where the final/expected sizes are known.

- Check for bugs and add (much) more testing, as the current tests barely
scratch the surface.

