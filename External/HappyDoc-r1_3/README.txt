HappyDoc Documentation Extraction Tool

  HappyDoc is a tool for extracting documentation from Python source
  code.  It differs from other such applications by the fact that it
  uses the parse tree for a module to derive the information used in
  its output, rather that importing the module directly.  This allows
  the user to generate documentation for modules which need special
  context to be imported.

Download

  Download the latest version of
  "HappyDoc":http://sourceforge.net/projects/happydoc/ from
  "SourceForge":http://sourceforge.net .

Installation

  HappyDoc uses the Distutils package for installation.  Unpack the
  tarball downloaded from SourceForge in a temporary directory.  Then
  run::

    % python ./setup.py install

  to install the package.

Instructions

  HappyDoc uses the Python parser to extract information from
  '__doc__' strings.  It also examines the signatures of functions and
  methods to include function prototypes in the output.  

  To use HappyDoc, simply run it against your source directory.  Use
  the '-h' or '--help' arguments to learn about the arguments and
  options accepted.  See below for more detailed directions about
  configuring your source for documentation purposes.

  Running HappyDoc

    After installation, the HappyDoc command-line application should
    be in your path.  Simply run 'happydoc' with appropriate
    arguments.

  Controlling the Output

    HappyDoc uses two different pluggable classes for controlling
    output.  A **formatter** class is responsible for producing the
    syntax of the actual output (e.g. HTML, XML, SGML, or PDF).  A
    **docset** class is responsible for controlling the formatter and
    managing the logical flow of the information (e.g., writing to
    multiple files or a single file, putting class documentation in a
    different file from the module, etc.).  Formatters and DocSets
    should be implemented so that any two can be combined.  It will
    not always be desirable to do this, but it should be possible.

  Doc-string Format

    *How does an author write documentation so that it will be marked
    up and look fancy?* This is a perennial question for
    "Python":http://www.python.org, and seems to have introduced a
    roadblock into the development of more robust and useful
    documentation tools.  By separating the formatter classes from the
    docset classes, HappyDoc allows a user to create their own
    formatter to interpret comments in any way they see fit.

    The default for the 'HTMLTableFormatter' (the default formatter
    for HappyDoc) is to treat '__doc__' strings as
    "StructuredText":http://www.python.org/sigs/doc-sig/stext.html.
    See also the "StructuredText
    module":HappyDoc/StructuredText.py.html for a description of the
    rules for using StructuredText.

    *Don't like StructuredText?* Write your own formatter that does
    something different and drop it into place.

  Documentation not in Doc-strings

    It is not always desirable to put all documentation in '__doc__'
    strings.  Sometime, notably when working with
    "Zope":http://www.zope.org , special meaning is attached to the
    presence of '__doc__' strings.  For this reason, and to support
    existing code which might not have '__doc__' strings, HappyDoc
    will find and extract documentation in Python comments.  

    Comment documentation can contain all of the same formatting as
    '__doc__' strings.  The preceding comment marker will be stripped
    off and the lines will be assembled and treated as a block of
    StructuredText.

    To use this feature, it is important to place the comments
    **before** the named object which they describe.  In this example:

      #
      # Class documentation goes here
      #
      class ClassWithNoDocStrings:
         "Using __doc__ strings overrides comment documentation."

         def method1(self, params):
             "This method uses a __doc__ string."
             pass

         #
         # Method2 does not use a __doc__ string.
         #
         def method2(self):
             pass

    The output would include the '__doc__' strings for the class and
    for 'method1'.  It would also make it appear that 'method2' had a
    '__doc__' string with the contents '"Method2 does not use a
    __doc__ string."'

  Formatters

    **htmltable** -- The 'HTMLTableFormatter' generates HTML output
    using tables to arrange and align the text.

    **text** -- The 'TextFormatter' generates plain text output with
    little reformatting of the input comments.

  DocSet types

    **file** -- The 'DocSet' sends the output of the formatter to
    multiple files.  Each module and class is documented in its own
    file.  Functions are documented along with the module to which
    they belong.  Methods are documented with the class.

    **single_file** -- The 'SingleFileDocSet' sends the output of the
    formatter to a single output file.

    **stdout** -- The 'StdOutDocSet' sends the output of the formatter
    to stdout.

Using HappyDoc

  Command Line Options

    HappyDoc uses the
    "CommandLineApp":HappyDoc/CommandLineApp.py_CommandLineApp.html
    class to process command line arguments.  To see the syntax help
    for the command line program, use the '-h' or '--help' options.
    The specific options supported are not documented here since they
    might change.

  DocSet and Formatter Parameters

    Many DocSets and Formatters will take parameters.  To pass
    parameters past the command line argument processing of HappyDoc
    and in to the DocSet or Formatter being used, the variable is
    passed as an argument rather than option (no dashes) to HappyDoc.

    To allow DocSets and Formatters to share variable namespace, the
    options passed are prefixed with a value indicating whether the
    variable is for the 'docset_' or 'formatter_'.

    For example::

      % ./happydoc.py -d MySources/doc MySources \
		docset_title='The title' \
		formatter_bgcolor1='#ffccaa'

  Input Types

    HappyDoc accepts 3 basic input types for documentation.  

    1. Any **file name** passed will be treated as a Python source file.
       The file will be parsed (but not imported) and the
       documentation text will be extracted from the resulting parse
       tree.

    2. Any **directory** passed will be interpreted to mean to document
       all of the files in that directory, so HappyDoc will recurse
       into the directory to find files.

    3. A **single, regular, text file** can be passed as the "package
       description file."  This file, defaulting to 'README.txt', will
       be interepreted as StructuredText and included in the place of
       a doc string in the generated 'index.html' file.

Examples

  Two example output documentation sets are available.

  HappyDoc

    Of course HappyDoc is used to produce its own documentation.  The
    most current version is available on the 
    "HappyDoc home page":http://happydoc.sourceforge.net.

  Zope

    A set of 
    "Zope source documentation":http://www.zope.org/Members/hellmann/ZopeSrcDocs 
    based on a recent CVS checkout available on "Zope.org":http://www.zope.org.

Support

  Contact "Doug Hellmann":doughellmann@bigfoot.com with questions.

  Please use the "bug tracker":http://sourceforge.net/bugs/?group_id=9678 
  on the SourceForge project page for HappyDoc to report bugs or
  request new features.

  There are also public forums available on SourceForge for discussing
  plans for future versions of HappyDoc.

Who else is using HappyDoc?

  **Biopython** -- The "Biopython project":http://www.biopython.org
  uses HappyDoc to generate the documentation for their 
  "libraries":http://www.biopython.org/wiki/html/BioPython/BiopythonCode.html

  **Numerical Python** -- "Numerical Python":http://numpy.sourceforge.net 
  adds a fast, compact multidimensional array language facility to Python.

To Do

  Features to be Implemented

    **New Formatters** -- 

      - DocBook

      - PDF

      - MIF

      - Create a 'formatter_structuredtext.py' module to output in
        StructuredText format, rather than just plain-text as
        'formatter_textfile.py' does.

    **New DocSets** -- 

      - DocBook (?)

      - PDF (?)

      - MIF (?)

    **Changes to Output** -- 

      - Add an index file with all objects listed in different ways.

      - Use pprint to format arguments to functions and methods. *Use
        a new class to represent argument default values.  Support
        '__str__' and 'pprint' methods.*

      - Use StructuredTextNG instead of current version to achieve
        better output.

      - Output a list of the inherited methods for a class as links to
        their definition along with the base classes.

      - Recognize references to methods, classes, etc. within
        '__doc__' strings and make them links where possible or
        otherwise highlighted when not.

      - Support variable substitution in the package description file.
        Define some variables, and allow the user to pass in their own.

      - Correctly recognize that two modules can define objects with
        the same name, but those objects are different things.  *Use
        the dotted path specified when something is reference to
        determine the right object.  Other ways?*

      - Recognize 'CommandLineApp' subclasses and generate the usage
        help automatically.  *This will require changing the
        'CommandLineApp' to support generating help without an
        instance, and would require importing the module.*

      - Use separate import lists for standard library and local
        modules.

      - **Fix** singlefile docset so that in HTML mode class names are
        output.  Is this problem part of the formatter?

      - Support for Python 2.0.

    **Distribution** -- 

      - Work on cvstagdiff.py to generate change reports for versions.

      - Add ZEO code to Zope CVS checkout.

    **General** -- 

      - Optimizations to make it run faster and smaller.

	  - Support for Python 2.x

    **Bugs** -- *See the SourceForge "bug
    tracker":http://sourceforge.net/bugs/?group_id=9678 page.*

Version History

  See "CHANGES":CHANGES.html

Known Bugs

  - Classes which inherit from classes with the same name cause too
  much recursion when the base class hierarchy is output.


