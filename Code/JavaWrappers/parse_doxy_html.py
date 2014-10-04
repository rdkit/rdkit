#! /usr/bin/python
# parse_doxy_html.py
""""Parse Doxygen-generated html files to get out stuff we want for Javadocs

Most code here works on doxytext:  this is text taken from Doxygen-generated html created by
processing the C++ code.  That html is viewed with Firefox and the appropriate pieces (now starting at
"Detailed Descripted"are just copied and pasted into a text file.  Note that some of the Doxygen-generated
files don't have that section and for now this program can't handle them without some additional annotation
by hand.
"""
from __future__ import print_function

from BeautifulSoup import *
import os
import re

def list_class_files(dir):
    return [name for name in os.listdir(dir)
             if (name.startswith('class_') and (not name.endswith('png')) and name.find('-members') == -1)]
    
def get_detail(fname):
    bs = BeautifulSoup(open(fname).read())
    det = bs.find(text='Detailed Description')
    return [bs,  det]

_example = \
"""int RDKit::Atom::getPerturbationOrder 	( 	INT_LIST  	probe 	 )  	const

returns the perturbation order for a list of integers

This value is associated with chirality.

Parameters:
    	probe 	a list of bond indices. This must be the same length as our number of incoming bonds (our degree).

Returns:
    the number of swaps required to convert the ordering of the probe list to match the order of our incoming bonds: e.g. if our incoming bond order is: [0,1,2,3]

    	getPerturbationOrder([1,0,2,3]) = 1
    	getPerturbationOrder([1,2,3,0]) = 3
    	getPerturbationOrder([1,2,0,3]) = 2
    	

See the class documentation for a more detailed description of our representation of chirality.

Usage

        ... molPtr is a const ROMol & ...
        ... atomPtr is a const Atom * ...
        ROMol::ADJ_ITER nbrIdx,endNbrs;
        boost::tie(nbrIdx,endNbrs) = molPtr.getAtomNeighbors(atomPtr);
        while(nbrIdx!=endNbrs){
          const ATOM_SPTR at=molPtr[*nbrIdx];
          ... do something with the Atom ...
          ++nbrIdx;
        }


Notes:

    * requires an owning molecule

"""

_atom = '''Detailed Description

The class for representing atoms.

Notes:

    * many of the methods of Atom require that the Atom be associated with a molecule (an ROMol).
    * each Atom maintains a Dict of properties:
          o Each property is keyed by name and can store an arbitrary type.
          o Properties can be marked as calculated, in which case they will be cleared when the clearComputedProps() method is called.
          o Because they have no impact upon chemistry, all property operations are const, this allows extra flexibility for clients who need to store extra data on Atom objects.
    * Atom objects are lazy about computing their explicit and implicit valence values. These will not be computed until their values are requested.

Chirality:

The chirality of an Atom is determined by two things:

    * its chiralTag
    * the input order of its bonds (see note below for handling of implicit Hs)

For tetrahedral coordination, the chiralTag tells you what direction you have to rotate to get from bond 2 to bond 3 while looking down bond 1. This is pretty much identical to the SMILES representation of chirality.

NOTE: if an atom has an implicit H, the bond to that H is considered to be at the *end* of the list of other bonds.
Member Enumeration Documentation
enum RDKit::Atom::ChiralType

store type of chirality

Enumerator:
    CHI_UNSPECIFIED 	

    chirality that hasn't been specified
    CHI_TETRAHEDRAL_CW 	

    tetrahedral: clockwise rotation (SMILES @)
    CHI_TETRAHEDRAL_CCW 	

    tetrahedral: counter-clockwise rotation (SMILES @)
    CHI_OTHER 	

    some unrecognized type of chirality

enum RDKit::Atom::HybridizationType

store hybridization

Enumerator:
    UNSPECIFIED 	

    hybridization that hasn't been specified
    OTHER 	

    unrecognized hybridization

Member Function Documentation
int RDKit::Atom::calcExplicitValence 	( 	bool  	strict = true 	 )  	

calculates and returns our explicit valence

Notes:

    * requires an owning molecule

int RDKit::Atom::calcImplicitValence 	( 	bool  	strict = true 	 )  	

calculates and returns our implicit valence

Notes:

    * requires an owning molecule

void RDKit::Atom::clearProp 	( 	const std::string  	key 	 )  	const [inline]

This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
void RDKit::Atom::clearProp 	( 	const char *  	key 	 )  	const [inline]

clears the value of a property

Notes:

    * if no property with name key exists, a KeyErrorException will be thrown.
    * if the property is marked as computed, it will also be removed from our list of computedProperties

Atom * RDKit::Atom::copy 	( 		 )  	const [virtual]

makes a copy of this Atom and returns a pointer to it.

Note: the caller is responsible for deleteing the result

Reimplemented in RDKit::QueryAtom.
unsigned int RDKit::Atom::getDegree 	( 		 )  	const

returns the explicit degree of the Atom (number of bonded neighbors in the graph)

Notes:

    * requires an owning molecule

int RDKit::Atom::getImplicitValence 	( 		 )  	const

returns the implicit valence for this Atom

Notes:

    * requires an owning molecule

unsigned int RDKit::Atom::getNumImplicitHs 	( 		 )  	const

returns the number of implicit Hs this Atom is bound to

Notes:

    * requires an owning molecule

unsigned int RDKit::Atom::getNumRadicalElectrons 	( 		 )  	const [inline]

returns the number of radical electrons for this Atom

Notes:

    * requires an owning molecule

int RDKit::Atom::getPerturbationOrder 	( 	INT_LIST  	probe 	 )  	const

returns the perturbation order for a list of integers

This value is associated with chirality.

Parameters:
    	probe 	a list of bond indices. This must be the same length as our number of incoming bonds (our degree).

Returns:
    the number of swaps required to convert the ordering of the probe list to match the order of our incoming bonds: e.g. if our incoming bond order is: [0,1,2,3]

    	getPerturbationOrder([1,0,2,3]) = 1
    	getPerturbationOrder([1,2,3,0]) = 3
    	getPerturbationOrder([1,2,0,3]) = 2
    	

See the class documentation for a more detailed description of our representation of chirality.

Notes:

    * requires an owning molecule

template<typename T >
void RDKit::Atom::getProp 	( 	const std::string  	key,
		T &  	res	 
	) 			const [inline]

This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
template<typename T >
void RDKit::Atom::getProp 	( 	const char *  	key,
		T &  	res	 
	) 			const [inline]

allows retrieval of a particular property value

Parameters:
    	key 	the name under which the property should be stored. If a property is already stored under this name, it will be replaced.
    	res 	a reference to the storage location for the value.

Notes:

    * if no property with name key exists, a KeyErrorException will be thrown.
    * the boost::lexical_cast machinery is used to attempt type conversions. If this fails, a boost::bad_lexical_cast exception will be thrown.

unsigned int RDKit::Atom::getTotalDegree 	( 		 )  	const

returns the total degree of the Atom (number of bonded neighbors + number of Hs)

Notes:

    * requires an owning molecule

unsigned int RDKit::Atom::getTotalNumHs 	( 	bool  	includeNeighbors = false 	 )  	const

returns the total number of Hs (implicit and explicit) that this Atom is bound to

Notes:

    * requires an owning molecule

bool RDKit::Atom::hasProp 	( 	const std::string  	key 	 )  	const [inline]

This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
bool RDKit::Atom::Match 	( 	const ATOM_SPTR  	what 	 )  	const [virtual]

This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.

Reimplemented in RDKit::QueryAtom.
bool RDKit::Atom::Match 	( 	Atom const *  	what 	 )  	const [virtual]

returns whether or not we match the argument

Notes:

    * for Atom objects, "match" means that atomic numbers are the same.

Reimplemented in RDKit::QueryAtom.
void RDKit::Atom::setIdx 	( 	unsigned int  	index 	 )  	[inline]

sets our index within the ROMol

Notes:

    * this makes no sense if we do not have an owning molecule
    * the index should be < this->getOwningMol()->getNumAtoms()

template<typename T >
void RDKit::Atom::setProp 	( 	const std::string  	key,
		T  	val,
		bool  	computed = false	 
	) 			const [inline]

This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
template<typename T >
void RDKit::Atom::setProp 	( 	const char *  	key,
		T  	val,
		bool  	computed = false	 
	) 			const [inline]

sets a property value

Parameters:
    	key 	the name under which the property should be stored. If a property is already stored under this name, it will be replaced.
    	val 	the value to be stored
    	computed 	(optional) allows the property to be flagged computed.

void RDKit::Atom::updatePropertyCache 	( 	bool  	strict = true 	 )  	

calculates any of our lazy properties

Notes:

    * requires an owning molecule
    * the current lazy properties are implicit and explicit valence

'''
_renote = re.compile('^\w*(Notes?[:]?)(?:.*?$)(.*?)((^\w)|\Z)',  flags=(re.M | re.I | re.DOTALL))
_reparam = re.compile('^\w*(Param(?:eter)?s?[:]?)(?:.*?$)(.*?)((^\w)|\Z)', flags=(re.M | re.I | re.DOTALL))
_rereturn = re.compile('^\w*(Returns[:])(?:.*?$)(.*?)((^\w)|\Z)', flags=(re.M | re.I | re.DOTALL))
_rereturn2 = re.compile('^\w*(Returns)\s+(.*?)((^\w)|\Z)', flags=(re.M | re.I | re.DOTALL))
_reusage = re.compile('^\w*(Usage[:]?)(?:.*?$)(.*?)((^\w)|\Z)', flags=(re.M | re.I | re.DOTALL))

def make_method_doc(doxy_method_text,  class_name):
    for f in (do_usage,  do_note,  do_param,  do_return):
        doxy_method_text = f(doxy_method_text)
    # Create paragraphs
    doxy_method_text = doxy_method_text.replace('\n\n', '\n<p>\n') 
    # But no paragraph markers just before tag
    doxy_method_text = doxy_method_text.replace('<p>@',  '<p>\n@')
    # Get rid of double quotes -- note that this causes an error with initialized string parameters
    doxy_method_text = doxy_method_text.replace('"',  "'")
    # Get lines 
    lines = doxy_method_text.split('\n')
    # Build header -- don't want type there
    header = lines[0][lines[0].find(class_name):]
    start = 1
    while (header.find(')') < 0):
        start += 1
        header += ' ' + lines[start - 1]
       ##  print header
    # Or [] annotation
    header = header.replace('[inline]',  '')
    header = header.replace('[virtual]',  '')
    header = header.replace('[protected]',  '')
    header = header.replace('[explicit]',  '')
    header = '%javamethodmodifiers ' +  header + '"\n/**\n'
    return header +  '\n'.join(lines[start:] ) + '\n*/\npublic";\n'

def make_class_doc(doxy_text,  class_name):
    for f in (do_usage,  do_note):
        doxy_text = f(doxy_text)
    # Create paragraphs
    doxy_text = doxy_text.replace('\n\n', '\n<p>\n') 
    # But no paragraph markers just before tag
    doxy_text = doxy_text.replace('<p>@',  '<p>\n@')
    # Get rid of double quotes
    doxy_text = doxy_text.replace('"',  "'")
    # Get lines
    lines = doxy_text.split('\n')
    # Build header -- don't want type there
    header = '%typemap(javaimports) ' + class_name + ' "\n/** '
    return header +  '\n'.join(lines) + ' */"\n'
    

def do_note(doxy_text):
    m1 = _renote.search(doxy_text)
    if m1 != None:
        repl = m1.group(0)
        if repl[-1] != '\n':
            repl = repl[:-1]
        new_text = '<p>@notes\n'
        for line in m1.group(2).split('\n'):
            line = line.strip()
            if len(line) > 0:
                if line.startswith('*'):
                    line = line[1:].strip()
                new_text = new_text + '<li>' + line + '\n'
        doxy_text = doxy_text.replace(repl, new_text) 
    return doxy_text

def do_param(doxy_text):
    m1 = _reparam.search(doxy_text)
    if m1 != None:
        repl = m1.group(0)
        if repl[-1] != '\n':
            repl = repl[:-1]
        new_text = '<p>@param\n'
        for line in m1.group(2).split('\n'):
            line = line.strip()
            if len(line) > 0:
                new_text = new_text + line + '\n'
        doxy_text = doxy_text.replace(repl, new_text) 
    return doxy_text

def do_return(doxy_text):
    m1 = _rereturn.search(doxy_text)
    if m1 == None:
        m1 = _rereturn2.search(doxy_text)
    if m1 != None:
        repl = m1.group(0)
        if repl[-1] != '\n':
            repl = repl[:-1]
        new_text = '<p>@return\n'
        for line in m1.group(2).split('\n'):
            line = line.strip()
            if len(line) > 0:
                new_text = new_text + line + '\n'
        doxy_text = doxy_text.replace(repl, new_text) 
    return doxy_text

def do_usage(doxy_text):
    m1 = _reusage.search(doxy_text)
    if m1 != None:
        repl = m1.group(0)
        if repl[-1] != '\n':
            repl = repl[:-1]
        new_text = '<p>@example\n<pre><code>\n'
        for line in m1.group(2).split('\n'):
            new_text = new_text + line + '\n'
        ## doxy_text = _reusage.sub(new_text,  doxy_text) + '</code></pre>\n'
        new_text += '</code></pre>\n'
        doxy_text = doxy_text.replace(repl, new_text) 
    return doxy_text

def do_methods(doxy_text,  class_name):
    methods = []
    method_lines = []
    in_method_region = False
    for line in doxy_text.split('\n'):
        if line.find('Function Documentation') >= 0:
            in_method_region = True
        elif in_method_region:
            if line.find(class_name) >= 0 or line.find('Member Data') >= 0:
                if len(method_lines) > 0:
                    methods.append(make_method_doc('\n'.join(method_lines),  class_name))
                if line.find('Member Data') >= 0:
                    in_method_region = False
                else:
                    method_lines = [line]
            else:
                method_lines.append(line)
        
    if len(method_lines) > 0:
        methods.append(make_method_doc('\n'.join(method_lines),  class_name))
    method_lines = [line]
    return methods

def do_class(doxy_text,  class_name):
    in_class_region = False
    class_doc = ''
    for line in doxy_text.split('\n'):
        if line.find('Detailed Description') >= 0:
            in_class_region = True
            class_lines = []
        elif in_class_region:
            if line.strip().endswith('Documentation') :
                if len(class_lines) > 0:
                    class_doc = make_class_doc('\n'.join(class_lines),  class_name)
                    in_class_region = False
            else:
                class_lines.append(line)
                    
    return class_doc

if __name__ == '__main__':
    import sys
    text = open(sys.argv[1]).read()
    class_name = sys.argv[2]
    print(do_class(text,  class_name))
    docs = do_methods(text,  class_name)
    for doc in docs:
        print(doc)
        
        
    
