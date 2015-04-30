//
//  Copyright (C) 2002-2013 greg landrum, Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_MOLSUPPLIER_H
#define _RD_MOLSUPPLIER_H

#include <RDGeneral/types.h>

#include <string>
#include <list>
#include <vector>
#include <iostream>
#include <GraphMol/ROMol.h>

namespace RDKit {
  std::string strip(const std::string &orig);

  /*! 
  //
  //  Here are a couple of ways one can interact with MolSuppliers:
  //
  //  1) Lazy (ForwardIterator):
  //     while(!supplier.atEnd()){
  //       ROMol *mol = supplier.next();
  //       if(mol){
  //           do something;
  //       }
  //     }
  //  2) Random Access:
  //     for(int i=0;i<supplier.length();i++){
  //       ROMol *mol = supplier[i];
  //       if(mol){
  //           do something;
  //       }
  //     }
  //
  //
  */
  class MolSupplier {
    // this is an abstract base class to supply molecules one at a time
  public:
    MolSupplier() {};
    virtual ~MolSupplier() {};
    virtual void init() = 0;
    virtual void reset() = 0;
    virtual bool atEnd() = 0;
    virtual ROMol *next() = 0;

  private:
    // disable automatic copy constructors and assignment operators
    // for this class and its subclasses.  They will likely be
    // carrying around stream pointers and copying those is a recipe
    // for disaster.
    MolSupplier(const MolSupplier&);
    MolSupplier &operator=(const MolSupplier&);
  protected:
    // stream to read the molecules from:
    std::istream *dp_inStream;
    // do we own dp_inStream?
    bool df_owner; 
  };


  // \brief a supplier from an SD file that only reads forward:
  class ForwardSDMolSupplier : public MolSupplier {
    /*************************************************************************
     * A lazy mol supplier from a SD file. 
     *  - When new molecules are read using "next" their positions in the file are noted. 
     ***********************************************************************************/
  public:
    ForwardSDMolSupplier() { init(); };

    explicit ForwardSDMolSupplier(std::istream *inStream, bool takeOwnership=true,
                                  bool sanitize=true,bool removeHs=true,
                                  bool strictParsing=false);

    virtual ~ForwardSDMolSupplier() {
      if (df_owner && dp_inStream) {
        delete dp_inStream;
        df_owner=false;
        dp_inStream=NULL;
      }
    };

    virtual void init();
    virtual void reset();
    virtual ROMol *next();
    virtual bool atEnd(); 

  protected:
    virtual void checkForEnd();
    ROMol *_next();
    virtual void readMolProps(ROMol *);
    bool df_end; 
    int d_line; // line number we are currently on
    bool df_sanitize,df_removeHs,df_strictParsing;
  };


  // \brief a lazy supplier from an SD file
  class SDMolSupplier : public ForwardSDMolSupplier {
    /*************************************************************************
     * A lazy mol supplier from a SD file. 
     *  - When new molecules are read using "next" their positions in the file are noted. 
     *  - A call to the "length" will automatically parse the entire file and cache all the mol
     *    block positions
     *  - [] operator is used to access a molecule at "idx", calling next following this will result
     *    in the next molecule after "idx"
     ***********************************************************************************/

  public:
    SDMolSupplier() { init(); };

    /*! 
     *   \param fileName - the name of the SD file
     *   \param sanitize - if true sanitize the molecule before returning it
     *   \param removeHs - if true remove Hs from the molecule before returning it
     *                     (triggers sanitization)
     *   \param strictParsing - if not set, the parser is more lax about correctness
     *                          of the contents.
     */
    explicit SDMolSupplier(const std::string &fileName, bool sanitize=true,
                           bool removeHs=true,bool strictParsing=true);
    
    explicit SDMolSupplier(std::istream *inStream, bool takeOwnership=true,
                           bool sanitize=true,bool removeHs=true,bool strictParsing=true);

    
    ~SDMolSupplier() {};
    void init();
    void reset();
    ROMol *next();
    bool atEnd(); 
    void moveTo(unsigned int idx);
    ROMol * operator[](unsigned int idx);
    /*! \brief returns the text block for a particular item
     *  
     *  \param idx - which item to return
     */
    std::string getItemText(unsigned int idx);
    unsigned int length();
    void setData(const std::string &text,bool sanitize=true, bool removeHs=true);

    /*! Resets our internal state and sets the indices of molecules in the stream.
     *  The client should be *very* careful about calling this method, as it's trivial
     *  to end up with a completely useless supplier.
     *
     *   \param locs - the vector of stream positions.
     *
     *  Note that this can be used not only to make reading selected molecules from a
     *  large SD file much faster, but it can also allow subsetting an SD file or
     *  rearranging the order of the molecules.
     */
    void setStreamIndices(const std::vector<std::streampos> &locs);

  private:
    void checkForEnd();
    int d_len; // total number of mol blocks in the file (initialized to -1)
    int d_last; // the molecule we are ready to read
    std::vector<std::streampos> d_molpos;

  };

  //! lazy file parser for Smiles tables
  class SmilesMolSupplier : public MolSupplier {
    /**************************************************************************
     * Lazy file parser for Smiles table file, similar to the lazy SD
     * file parser above
     * - As an when new molecules are read using "next" their
     *    positions in the file are noted.
     *  - A call to the "length" will autamatically parse the entire
     *    file and cache all the mol block positions
     *  - [] operator is used to access a molecule at "idx", calling
     *    next following this will result in the next molecule after
     *    "idx"
     ***************************************************************************/ 
  public:

    /*! 
     *   \param fileName - the name of smiles table file
     *   \param delimiter - delimiting characters between records on a each
     *     line NOTE that this is not a string, the tokenizer looks for
     *     the individual characters in delimiter, not the full string
     *     itself.  So the default delimiter: " \t", means " " or "\t".
     *   \param smilesColumn - column number for the SMILES string (defaults
     *     to the first column)
     *   \param nameColumn - column number for the molecule name (defaults to
     *     the second column) If set to -1 we assume that no name is
     *     available for the molecule and the name is defaulted to the
     *     smiles string
     *   \param titleLine - if true, the first line is assumed to list the
     *     names of properties in order seperated by 'delimiter'. It is
     *     also assume that the 'SMILES' column and the 'name' column
     *     are not specified here if false - no title line is assumed
     *     and the properties are recorded as the "columnX" where "X" is
     *     the column number
     *   \param sanitize - if true sanitize the molecule before returning it
     */
    explicit SmilesMolSupplier(const std::string &fileName, 
			       const std::string &delimiter=" \t",
			       int smilesColumn=0,
			       int nameColumn=1, 
			       bool titleLine=true,		   
                               bool sanitize=true);
    SmilesMolSupplier();
    explicit SmilesMolSupplier(std::istream *inStream, bool takeOwnership=true,
			       const std::string &delimiter=" \t",
			       int smilesColumn=0,
			       int nameColumn=1, 
			       bool titleLine=true,		   
                               bool sanitize=true);                               

    ~SmilesMolSupplier();
    void setData(const std::string &text,
		 const std::string &delimiter=" ",
		 int smilesColumn=0,
		 int nameColumn=1, 
		 bool titleLine=true,		   
		 bool sanitize=true);
    void init();
    void reset();
    ROMol *next();
    bool atEnd();
    void moveTo(unsigned int idx);
    ROMol * operator[](unsigned int idx);
    /*! \brief returns the text block for a particular item
     *  
     *  \param idx - which item to return
     */
    std::string getItemText(unsigned int idx);
    unsigned int length();

  private:
    ROMol *processLine(std::string inLine);
    void processTitleLine();
    std::string nextLine();
    long int skipComments();
    void checkForEnd();
    
    bool df_end; // have we reached the end of the file?
    int d_len; // total number of smiles in the file
    int d_next; // the  molecule we are ready to read
    int d_line; // line number we are currently on
    std::vector<std::streampos> d_molpos; // vector of positions in the file for molecules
    std::vector<int> d_lineNums; 
    std::string d_delim; // the delimiter string
    bool df_sanitize; // sanitize molecules before returning them?
    STR_VECT d_props; // vector of property names
    bool df_title; // do we have a title line?
    int d_smi; // column id for the smile string
    int d_name; // column id for the name
  };

  //! lazy file parser for TDT files
  class TDTMolSupplier : public MolSupplier {
    /**************************************************************************
     * Lazy file parser for TDT files, similar to the lazy SD
     * file parser above
     * - As an when new molecules are read using "next" their
     *    positions in the file are noted.
     *  - A call to the "length" will autamatically parse the entire
     *    file and cache all the mol block positions
     *  - [] operator is used to access a molecule at "idx", calling
     *    next following this will result in the next molecule after
     *    "idx"
     ***************************************************************************/ 
  public:

    /*! 
     *   \param fileName - the name of the TDT file
     *   \param nameRecord - property name for the molecule name.
     *     If empty (the default), the name defaults to be empty
     *   \param confId2D - if >=0 and 2D coordinates are provided, the 2D
     *                   structure (depiction) in the input will be read into the
     *                   corresponding conformer id.
     *   \param confId3D - if >=0 and 3D coordinates are provided, the 3D
     *                   structure (depiction) in the input will be read into the
     *                   corresponding conformer id.
     *   \param sanitize - if true sanitize the molecule before returning it
     */
    explicit TDTMolSupplier(const std::string &fileName, 
			    const std::string &nameRecord="",
			    int confId2D=-1,int confId3D=0,
			    bool sanitize=true);
    explicit TDTMolSupplier(std::istream *inStream, bool takeOwnership=true,
			    const std::string &nameRecord="",
			    int confId2D=-1,int confId3D=0,
			    bool sanitize=true);
    TDTMolSupplier();
    ~TDTMolSupplier();
    void setData(const std::string &text,
		 const std::string &nameRecord="",
		 int confId2D=-1,int confId3D=0,
		 bool sanitize=true);
    void init();
    void reset();
    ROMol *next();
    bool atEnd();
    void moveTo(unsigned int idx);
    ROMol * operator[](unsigned int idx);
    /*! \brief returns the text block for a particular item
     *  
     *  \param idx - which item to return
     */
    std::string getItemText(unsigned int idx);
    unsigned int length();

  private:
    bool advanceToNextRecord();
    void checkForEnd();
    ROMol *parseMol(std::string inLine);

    bool df_end; // have we reached the end of the file?
    int d_len; // total number of mols in the file
    int d_last; // the molecule we are ready to read
    int d_line; // line number we are currently on
    int d_confId2D; // id to use for 2D conformers
    int d_confId3D; // id to use for 3D conformers
    std::vector<std::streampos> d_molpos; // vector of positions in the file for molecules
    bool df_sanitize; // sanitize molecules before returning them?
    std::string d_nameProp; // local storage for the property providing mol names
  };

  //! lazy file parser for PDB files
  class PDBMolSupplier : public MolSupplier {
  public:
    explicit PDBMolSupplier(std::istream *inStream, bool takeOwnership=true,
                            bool sanitize=true, bool removeHs=true,
                            unsigned int flavor=0);
    explicit PDBMolSupplier(const std::string &fname,
                            bool sanitize=true, bool removeHs=true,
                            unsigned int flavor=0);

    virtual ~PDBMolSupplier() {
      if (df_owner && dp_inStream)
        delete dp_inStream;
    };

    virtual void init();
    virtual void reset();
    virtual ROMol *next();
    virtual bool atEnd();

  protected:
    bool df_sanitize,df_removeHs;
    unsigned int d_flavor;
  };
}

#endif
