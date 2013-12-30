//
//  Copyright (C) 2002-2013 Greg Landrum, Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef _RD_MOLWRITERS_H_
#define _RD_MOLWRITERS_H_

#include <RDGeneral/types.h>

#include <string>
#include <iostream>
#include <GraphMol/ROMol.h>

namespace RDKit {

  static int defaultConfId=-1;
  class MolWriter {
  public:
    virtual ~MolWriter() {}
    virtual void write(const ROMol &mol,int confId=defaultConfId) = 0;
    virtual void flush() = 0;
    virtual void close() = 0;
    virtual void setProps(const STR_VECT &propNames)=0;
    virtual unsigned int numMols() const =0;
  };

  //! The SmilesWriter is for writing molecules and properties to
  //! delimited text files.
  class SmilesWriter : public MolWriter {
    /******************************************************************************
     * A Smiles Table writer - this is how it is used
     *  - create a SmilesWriter with a output file name (or a ostream), a delimiter,
     *     and a list of properties that need to be written out
     *  - then a call is made to the write function for each molecule that needs to
     *     be written out
     ******************************************************************************/
  public:
    /*!
      \param fileName       : filename to write to ("-" to write to stdout)
      \param delimiter      : delimiter to use in the text file
      \param nameHeader     : used to label the name column in the output. If this
                              is provided as the empty string, no names will be written.
      \param includeHeader  : toggles inclusion of a header line in the output
      \param isomericSmiles : toggles generation of isomeric SMILES
      \param kekuleSmiles   : toggles the generation of kekule SMILES

     */
    SmilesWriter(std::string fileName, 
		 std::string delimiter=" ",
		 std::string nameHeader="Name",
		 bool includeHeader=true,
                 bool isomericSmiles=false,
                 bool kekuleSmiles=false);
    //! \overload
    SmilesWriter(std::ostream *outStream, 
		 std::string delimiter=" ",
		 std::string nameHeader="Name",
		 bool includeHeader=true,
		 bool takeOwnership=false,
                 bool isomericSmiles=false,
                 bool kekuleSmiles=false);
		 
    ~SmilesWriter();

    //! \brief set a vector of property names that are need to be
    //! written out for each molecule
    void setProps(const STR_VECT &propNames);

    //! \brief write a new molecule to the file
    void write(const ROMol &mol,int confId=defaultConfId);

    //! \brief flush the ostream
    void flush() {
      PRECONDITION(dp_ostream,"no output stream");
      dp_ostream->flush();
    };

    //! \brief close our stream (the writer cannot be used again)
    void close() {
      PRECONDITION(dp_ostream,"no output stream");
      dp_ostream->flush();
      if(df_owner) {
        delete dp_ostream;
        df_owner=false;
      }
      dp_ostream=NULL;
    };

    //! \brief get the number of molecules written so far
    unsigned int numMols() const { return d_molid;} ;

  private:
    // local initialization
    void init(std::string delimiter,std::string nameHeader,
              bool includeHeader,
              bool isomericSmiles,
              bool kekuleSmiles);


    // dumps a header line to the output stream
    void dumpHeader() const;


    std::ostream *dp_ostream;
    bool df_owner;
    bool df_includeHeader; // whether or not to include a title line
    unsigned int d_molid; // the number of the molecules we wrote so far
    std::string d_delim; // delimiter string between various records
    std::string d_nameHeader; // header for the name column in the output file
    STR_VECT d_props; // list of property name that need to be written out
    bool df_isomericSmiles; // whether or not to do isomeric smiles
    bool df_kekuleSmiles; // whether or not to do kekule smiles
  };


  //! The SDWriter is for writing molecules and properties to
  //! SD files 
  class SDWriter : public MolWriter {
    /**************************************************************************************
     * A SD file ( or stream) writer - this is how it is used
     *  - create a SDMolWriter with a output file name (or a ostream),
     *     and a list of properties that need to be written out
     *  - then a call is made to the write function for each molecule that needs to be written out
     **********************************************************************************************/
  public:
    /*!
      \param fileName       : filename to write to ("-" to write to stdout)
     */
    SDWriter(std::string fileName);
    SDWriter(std::ostream *outStream,bool takeOwnership=false);

    ~SDWriter();

    //! \brief set a vector of property names that are need to be
    //! written out for each molecule
    void setProps(const STR_VECT &propNames);

    //! \brief write a new molecule to the file
    void write(const ROMol &mol, int confId=defaultConfId);

    //! \brief flush the ostream
    void flush() { 
      PRECONDITION(dp_ostream,"no output stream");
      dp_ostream->flush();
    } ;

    //! \brief close our stream (the writer cannot be used again)
    void close() {
      PRECONDITION(dp_ostream,"no output stream");
      dp_ostream->flush();
      if(df_owner) {
        delete dp_ostream;
        df_owner=false;
      }
      dp_ostream=NULL;
    };

    //! \brief get the number of molecules written so far
    unsigned int numMols() const { return d_molid; };

    void setForceV3000(bool val) { df_forceV3000=val; };
    bool getForceV3000() const { return df_forceV3000; };    
    
    void setKekulize(bool val) { df_kekulize=val; };
    bool getKekulize() const { return df_kekulize; };    
    
  private:
    void writeProperty(const ROMol &mol, std::string name);

    std::ostream *dp_ostream;
    bool df_owner;
    unsigned int d_molid; // the number of the molecules we wrote so far
    STR_VECT d_props; // list of property name that need to be written out
    bool df_forceV3000; // force writing the mol blocks as V3000
    bool df_kekulize; // toggle kekulization of molecules on writing
  };

  //! The TDTWriter is for writing molecules and properties to
  //! TDT files 
  class TDTWriter : public MolWriter {
    /**************************************************************************************
     * A TDT file ( or stream) writer - this is how it is used
     *  - create a TDTWriter with a output file name (or a ostream),
     *     and a list of properties that need to be written out
     *  - then a call is made to the write function for each molecule that needs to be written out
     **********************************************************************************************/
  public:
    /*!
      \param fileName       : filename to write to ("-" to write to stdout)
     */
    TDTWriter(std::string fileName);
    TDTWriter(std::ostream *outStream,bool takeOwnership=false);

    ~TDTWriter();

    //! \brief set a vector of property names that are need to be
    //! written out for each molecule
    void setProps(const STR_VECT &propNames);

    //! \brief write a new molecule to the file
    void write(const ROMol &mol, int confId=defaultConfId);

    //! \brief flush the ostream
    void flush() { 
      PRECONDITION(dp_ostream,"no output stream");
      dp_ostream->flush();
    };

    //! \brief close our stream (the writer cannot be used again)
    void close() {
      PRECONDITION(dp_ostream,"no output stream");
      dp_ostream->flush();
      if(df_owner) {
        delete dp_ostream;
        df_owner=false;
      }
      dp_ostream=NULL;
    };

    //! \brief get the number of molecules written so far
    unsigned int numMols() const { return d_molid; };

    void setWrite2D(bool state=true) { df_write2D=state; };
    bool getWrite2D() const { return df_write2D; };

    void setWriteNames(bool state=true) { df_writeNames=state; };
    bool getWriteNames() const { return df_writeNames; };

    void setNumDigits(unsigned int numDigits) { d_numDigits=numDigits; };
    unsigned int getNumDigits() const { return d_numDigits;};
    
  private:
    void writeProperty(const ROMol &mol, std::string name);

    std::ostream *dp_ostream;
    bool df_owner;
    unsigned int d_molid; // the number of molecules we wrote so far
    STR_VECT d_props; // list of property name that need to be written out
    bool df_write2D; // write 2D coordinates instead of 3D
    bool df_writeNames; // write a name record for each molecule
    unsigned int d_numDigits; // number of digits to use in our output of coordinates;
  };

  //! The PDBWriter is for writing molecules to Brookhaven Protein
  //! DataBank format files.
  class PDBWriter : public MolWriter {
  public:
    PDBWriter(std::string fileName, unsigned int flavor = 0);
    PDBWriter(std::ostream *outStream, bool takeOwnership=false,
              unsigned int flavor = 0);
    ~PDBWriter();

    //! \brief write a new molecule to the file
    void write(const ROMol &mol, int confId=defaultConfId);

    void setProps(const STR_VECT&) {};

    //! \brief flush the ostream
    void flush() { 
      PRECONDITION(dp_ostream,"no output stream");
      dp_ostream->flush();
    } ;

    //! \brief close our stream (the writer cannot be used again)
    void close() {
      PRECONDITION(dp_ostream,"no output stream");
      dp_ostream->flush();
      if(df_owner) {
        delete dp_ostream;
        df_owner=false;
      }
      dp_ostream=NULL;
    };

    //! \brief get the number of molecules written so far
    unsigned int numMols() const { return d_count;} ;

  private:
    std::ostream *dp_ostream;
    unsigned int d_flavor;
    unsigned int d_count;
    bool df_owner;
  };
 

}

#endif

