//
//  Copyright (C) 2013 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_PDBWRITER_H_
#define _RD_PDBWRITER_H_
#include <GraphMol/FileParsers/MolWriters.h>

namespace RDKit {

  //! The PDBWriter is for writing molecules to Brookhaven Protein
  //! DataBank format files.
  class PDBWriter : public MolWriter {
  public:
    PDBWriter(std::string fileName, unsigned int flavor = 0);
    PDBWriter(std::ostream *outStream, bool takeOwnership=false,
              unsigned int flavor = 0);
    ~PDBWriter();

    //! \brief write a new molecule to the file
    void write(ROMol &mol, int confId=defaultConfId);

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

#endif  // _RD_PDBWRITER_H_

