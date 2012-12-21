SOURCES += qtdraw_demo.cc RDKitMolToQPainter.cc

HEADERS += ../../../../GraphMol/MolDrawing/MolDrawing.h

TARGET = qtdraw_demo

CONFIG += qt

INCLUDEPATH += ${RDBASE}/Code

RD_STATIC_LIBS = -lSmilesParse_static -lDepictor_static -lGraphMol_static \
-lRDGeometryLib_static -lRDGeneral_static 

LIBS += -L${RDBASE}/lib $${RD_STATIC_LIBS}

LIBS += -L${BOOST}/lib -lboost_regex -lboost_system
