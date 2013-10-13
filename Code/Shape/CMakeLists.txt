rdkit_library(Shape
              alignmentInfo.cpp  
              bestResults.cpp  
              gaussianVolume.cpp  
              mainErr.cpp           
              options.cpp           
              printHeader.cpp  
              shapeAlignment.cpp  
              solutionInfo.cpp
              atomGaussian.cpp   
              coordinate.cpp   
              main.cpp            
              moleculeRotation.cpp  
              parseCommandLine.cpp  
              printUsage.cpp   
              siMath.cpp
              shape.cpp
              SHARED
              LINK_LIBRARIES  GraphMol RDGeneral RDGeometryLib DataStructs Shape FileParsers RDGeneral RDGeometryLib RDBoost)

rdkit_headers(alignmentInfo.h 
              atomGaussian.h 
              bestResults.h 
              config.h 
              coordinate.h 
              gaussianVolume.h 
              mainErr.h 
              moleculeRotation.h 
              options.h 
              parseCommandLine.h 
              printHeader.h 
              printUsage.h 
              shapeAlignment.h 
              siMath.h 
              shape.h
              solutionInfo.h DEST Shape)

add_subdirectory(Wrap)