rdkit_library(RDShapeLib
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
              LINK_LIBRARIES RDGeometryLib RDGeneral)

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

