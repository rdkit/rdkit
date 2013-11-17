rdkit_library(Shape
              alignmentInfo.cpp  
              bestResults.cpp  
              gaussianVolume.cpp  
              shapeAlignment.cpp  
              solutionInfo.cpp
              atomGaussian.cpp   
              moleculeRotation.cpp  
              siMath.cpp
              shape.cpp
              SHARED
              LINK_LIBRARIES  GraphMol RDGeneral RDGeometryLib DataStructs Shape Alignment EigenSolvers FileParsers RDGeneral RDGeometryLib RDBoost)

rdkit_headers(alignmentInfo.h 
              atomGaussian.h 
              bestResults.h 
              config.h 
              gaussianVolume.h 
              moleculeRotation.h 
              shapeAlignment.h 
              siMath.h 
              shape.h
              solutionInfo.h DEST Shape)

add_subdirectory(Wrap)
