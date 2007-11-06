
libs:
	$(MAKE) --directory=RDGeneral install;
	$(MAKE) --directory=Geometry install;
	$(MAKE) --directory=DataStructs;
	$(MAKE) --directory=Catalogs;
	$(MAKE) --directory=Query;
	$(MAKE) --directory=Numerics/Optimizer
	$(MAKE) --directory=Numerics/EigenSolvers
	$(MAKE) --directory=Numerics/Alignment;
	$(MAKE) --directory=ForceField
	$(MAKE) --directory=DistGeom
	$(MAKE) --directory=ChemicalFeatures;
	$(MAKE) --directory=GraphMol all;
	$(MAKE) --directory=GraphMol install;
	$(MAKE) --directory=ML/Cluster/Murtagh; 
	$(MAKE) --directory=ML/InfoTheory;
	$(MAKE) --directory=DataManip/MetricMatrixCalc;
	$(MAKE) --directory=SimDivPickers;

cartridge:
	$(MAKE) --directory=PgSQL/RDLib -f Makefile.linux;

wrappers:
	$(MAKE) --directory=RDBoost;
	$(MAKE) --directory=RDBoost wrapper
	$(MAKE) --directory=ML/Cluster/Murtagh wrapper
	$(MAKE) --directory=ML/InfoTheory/Wrap wrapper
	$(MAKE) --directory=ML/Data wrapper
	cd ML/FeatureSelect/Wrap; python setup.py build && python setup.py install --install-lib=$(RDBASE)/Python
	$(MAKE) --directory=DataStructs/Wrap;
	cd Geometry/Wrap; python setup.py build && python setup.py install --install-lib=$(RDBASE)/Python
	$(MAKE) --directory=DataManip/MetricMatrixCalc/Wrap;
	$(MAKE) --directory=SimDivPickers/Wrap;
	cd ForceField/Wrap; python setup.py build && python setup.py install --install-lib=$(RDBASE)/Python
	cd DistGeom/Wrap; python setup.py build && python setup.py install --install-lib=$(RDBASE)/Python
	cd Numerics/Alignment/Wrap; python setup.py build && python setup.py install --install-lib=$(RDBASE)/Python
	cd ChemicalFeatures/Wrap; python setup.py build && python setup.py install --install-lib=$(RDBASE)/Python
	$(MAKE) --directory=GraphMol wrappers;

all:
	$(MAKE) libs;
	$(MAKE) wrappers;
#	$(MAKE) cartridge;


regrs:
	$(MAKE) --directory=RDGeneral regrs;
	$(MAKE) --directory=DataStructs regrs;
	$(MAKE) --directory=Query regrs;
	$(MAKE) --directory=DataManip/MetricMatrixCalc regrs;
	$(MAKE) --directory=Numerics regrs;
	$(MAKE) --directory=Numerics/Optimizer regrs;
	$(MAKE) --directory=Numerics/EigenSolvers regrs;
	$(MAKE) --directory=Numerics/Alignment regrs;
	$(MAKE) --directory=DistGeom regrs;
	$(MAKE) --directory=ForceField regrs;
	$(MAKE) --directory=ChemicalFeatures regrs;
	$(MAKE) --directory=GraphMol allregrs;

docs:
	/usr/bin/doxygen doxygen.config

clean:
	$(MAKE) --directory=RDGeneral clean;
	$(MAKE) --directory=Geometry clean;
	$(MAKE) --directory=RDBoost clean;
	$(MAKE) --directory=DataStructs clean;
	$(MAKE) --directory=Catalogs clean;
	$(MAKE) --directory=Query clean;
	$(MAKE) --directory=DataManip/MetricMatrixCalc clean;
	$(MAKE) --directory=SimDivPickers clean;
	$(MAKE) --directory=GraphMol allclean;
	$(MAKE) --directory=DataStructs/Wrap clean;
	$(MAKE) --directory=DataManip/MetricMatrixCalc/Wrap clean;
	$(MAKE) --directory=SimDivPickers/Wrap clean;
	$(MAKE) --directory=Numerics/Optimizer clean;
	$(MAKE) --directory=Numerics/EigenSolvers clean;
	$(MAKE) --directory=Numerics/Alignment clean;
	$(MAKE) --directory=ForceField clean;
	$(MAKE) --directory=DistGeom clean;
	$(MAKE) --directory=ChemicalFeatures clean;
	$(MAKE) --directory=PgSQL/RDLib -f Makefile.linux clean; 
