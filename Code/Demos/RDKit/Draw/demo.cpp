// $Id: sample.cpp 793 2008-08-17 14:33:30Z glandrum $
//
//  Copyright (C) 2009 Greg Landrum
//   All Rights Reserved
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <Demos/RDKit/Draw/MolDrawing.h>

#include <RDGeneral/RDLog.h>
#include <vector>


using namespace RDKit;
void DrawDemo(){
  RWMol *mol=SmilesToMol("Clc1ccc(C(=O)NCc2sccc2)cc1");
  //RWMol *mol=SmilesToMol("c1ncncn1");
  RDKit::MolOps::Kekulize(*mol);

  // generate the 2D coordinates:
  RDDepict::compute2DCoords(*mol);
  std::vector<int> drawing=RDKit::Drawing::DrawMol(*mol);

  std::cout<<"var codes=[";
  std::copy(drawing.begin(),drawing.end(),std::ostream_iterator<int>(std::cout,","));
  std::cout<<"];"<<std::endl;

  delete mol;
}

int
main(int argc, char *argv[])
{
  RDLog::InitLogs();
  DrawDemo();
}
