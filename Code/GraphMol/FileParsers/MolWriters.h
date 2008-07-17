//
//  Copyright (C) 2002-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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
    virtual void write(ROMol &mol,int confId=defaultConfId) = 0;
    virtual void flush() = 0;
    virtual void setProps(const STR_VECT &propNames)=0;
    virtual unsigned int numMols() const =0;
  };

  class SmilesWriter : public MolWriter {
    /******************************************************************************
     * A Smiles Table writer - this is how it is used
     *  - create a SDMolWriter with a output file name (or a ostream), a delimiter,
     *     and a list of properties that need to be written out
     *  - then a call is made to the write function for each molecule that needs to
     *     be written out
     ******************************************************************************/
  public:
    SmilesWriter(std::string fileName, 
		 std::string delimiter=" ",
		 std::string nameHeader="Name",
		 bool includeHeader=true,
                 bool isomericSmiles=false,
                 bool kekuleSmiles=false);
    SmilesWriter(std::ostream *outStream, 
		 std::string delimiter=" ",
		 std::string nameHeader="Name",
		 bool includeHeader=true,
		 bool takeOwnership=false,
                 bool isomericSmiles=false,
                 bool kekuleSmiles=false);
		 
    ~SmilesWriter();

    // FIX; unfornately we cannot handle this through he constructor for now
    // because of the conversion issues betweem python list of strings and
    // STR_VECT - this is currently handled through a wrapper function and I do not
    // know how to wrap a constructor.
    // set a vector of property names that are need to be
    // written out for each molecule
    void setProps(const STR_VECT &propNames);

    // write a new molecule to the file
    // NOTE: the default value for the confId parameter is set in the
    //   baseClass, do *not* redefine it.
    void write(ROMol &mol,int confId=defaultConfId);

    // flush the ostream
    void flush() {
      PRECONDITION(dp_ostream,"no output stream");
      dp_ostream->flush();
    };

    // get the number of molecules written so far
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


  class SDWriter : public MolWriter {
    /**************************************************************************************
     * A SD file ( or stream) writer - this is how it is used
     *  - create a SDMolWriter with a output file name (or a ostream),
     *     and a list of properties that need to be written out
     *  - then a call is made to the write function for each molecule that needs to be written out
     **********************************************************************************************/
  public:
    SDWriter(std::string fileName);
    SDWriter(std::ostream *outStream,bool takeOwnership=false);

    ~SDWriter();

    // FIX; unfornately we cannot handle this through he constructor for now
    // because of the conversion issues betweem python list of strings and
    // STR_VECT - this is currently handled through a wrapper function and I do not
    // know how to wrap a constructor.
    // set a vector of property names that are need to be
    // written out for each molecule
    void setProps(const STR_VECT &propNames);

    // write a new molecule to the file
    // NOTE: the default value for the confId parameter is set in the
    //   baseClass, do *not* redefine it.
    void write(ROMol &mol, int confId=defaultConfId);

    // flush the ostream
    void flush() { 
      PRECONDITION(dp_ostream,"no output stream");
      dp_ostream->flush();
    } ;

    // get the number of molecules written so far
    unsigned int numMols() const { return d_molid; };

  private:
    void writeProperty(const ROMol &mol, std::string name);

    std::ostream *dp_ostream;
    bool d_owner;
    unsigned int d_molid; // the number of the molecules we wrote so far
    STR_VECT d_props; // list of property name that need to be written out
  };

  class TDTWriter : public MolWriter {
    /**************************************************************************************
     * A TDT file ( or stream) writer - this is how it is used
     *  - create a TDTWriter with a output file name (or a ostream),
     *     and a list of properties that need to be written out
     *  - then a call is made to the write function for each molecule that needs to be written out
     **********************************************************************************************/
  public:
    TDTWriter(std::string fileName);
    TDTWriter(std::ostream *outStream,bool takeOwnership=false);

    ~TDTWriter();

    // FIX; unfornately we cannot handle this through he constructor for now
    // because of the conversion issues betweem python list of strings and
    // STR_VECT - this is currently handled through a wrapper function and I do not
    // know how to wrap a constructor.
    // set a vector of property names that are need to be
    // written out for each molecule
    void setProps(const STR_VECT &propNames);

    // write a new molecule to the file
    // NOTE: the default value for the confId parameter is set in the
    //   baseClass, do *not* redefine it.
    void write(ROMol &mol, int confId=defaultConfId);

    // flush the ostream
    void flush() { 
      PRECONDITION(dp_ostream,"no output stream");
      dp_ostream->flush();
    };

    // get the number of molecules written so far
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
    bool d_owner;
    unsigned int d_molid; // the number of molecules we wrote so far
    STR_VECT d_props; // list of property name that need to be written out
    bool df_write2D; // write 2D coordinates instead of 3D
    bool df_writeNames; // write a name record for each molecule
    unsigned int d_numDigits; // number of digits to use in our output of coordinates;
  };

}

#endif

