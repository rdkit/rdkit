from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps

class RDKitRecipe(ConanFile):
    name = "rdkit"
    version = "0.0.1"
    package_type = "library"

    # Optional metadata
    license = "<Put the package license here>"
    author = "<Put your name here> <And your email here>"
    url = "<Package recipe repository url here, for issues about the package>"
    description = "<Description of RDKit package here>"
    topics = ("<Put some tag here>", "<here>", "<and here>")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": True, "fPIC": True}

    requires = "freetype/2.10.4", "boost/1.83.0"# "eigen/3.3.9"# "boost/1.83.0"

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = "CMakeLists.txt", "*"

    def config_options(self):
        if self.settings.os == "Windows":
            self.options.rm_safe("fPIC")

    def configure(self):
        if self.options.shared:
            self.options.rm_safe("fPIC")

    def layout(self):
        cmake_layout(self)

    def generate(self):
        deps = CMakeDeps(self)
        deps.generate()
        tc = CMakeToolchain(self)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.sourcedirs = ["Code"]
        self.cpp_info.components["RDKitAbbreviations"].libs = ["RDKitAbbreviations"]
        self.cpp_info.components["RDKitAlignment"].libs = ["RDKitAlignment"]
        self.cpp_info.components["RDKitCatalogs"].libs = ["RDKitCatalogs"]
        self.cpp_info.components["RDKitChemicalFeatures"].libs = ["RDKitChemicalFeatures"]
        self.cpp_info.components["RDKitChemReactions"].libs = ["RDKitChemReactions"]
        self.cpp_info.components["RDKitChemTransforms"].libs = ["RDKitChemTransforms"]
        self.cpp_info.components["RDKitCIPLabeler"].libs = ["RDKitCIPLabeler"]
        self.cpp_info.components["RDKitcoordgen"].libs = ["RDKitcoordgen"]
        self.cpp_info.components["RDKitDataStructs"].libs = ["RDKitDataStructs"]
        self.cpp_info.components["RDKitDepictor"].libs = ["RDKitDepictor"]
        self.cpp_info.components["RDKitDeprotect"].libs = ["RDKitDeprotect"]
        self.cpp_info.components["RDKitDescriptors"].libs = ["RDKitDescriptors"]
        self.cpp_info.components["RDKitDistGeometry"].libs = ["RDKitDistGeometry"]
        self.cpp_info.components["RDKitDistGeomHelpers"].libs = ["RDKitDistGeomHelpers"]
        self.cpp_info.components["RDKitEigenSolvers"].libs = ["RDKitEigenSolvers"]
        self.cpp_info.components["RDKitFileParsers"].libs = ["RDKitFileParsers"]
        self.cpp_info.components["RDKitFilterCatalog"].libs = ["RDKitFilterCatalog"]
        self.cpp_info.components["RDKitFingerprints"].libs = ["RDKitFingerprints"]
        self.cpp_info.components["RDKitFMCS"].libs = ["RDKitFMCS"]
        self.cpp_info.components["RDKitForceFieldHelpers"].libs = ["RDKitForceFieldHelpers"]
        self.cpp_info.components["RDKitForceField"].libs = ["RDKitForceField"]
        self.cpp_info.components["RDKitFragCatalog"].libs = ["RDKitFragCatalog"]
        self.cpp_info.components["RDKitga"].libs = ["RDKitga"]
        self.cpp_info.components["RDKitGeneralizedSubstruct"].libs = ["RDKitGeneralizedSubstruct"]
        self.cpp_info.components["RDKitGenericGroups"].libs = ["RDKitGenericGroups"]
        self.cpp_info.components["RDKitGraphMol"].libs = ["RDKitGraphMol"]
        self.cpp_info.components["RDKithc"].libs = ["RDKithc"]
        self.cpp_info.components["RDKitInfoTheory"].libs = ["RDKitInfoTheory"]
        self.cpp_info.components["RDKitmaeparser"].libs = ["RDKitmaeparser"]
        self.cpp_info.components["RDKitMarvinParser"].libs = ["RDKitMarvinParser"]
        self.cpp_info.components["RDKitMMPA"].libs = ["RDKitMMPA"]
        self.cpp_info.components["RDKitMolAlign"].libs = ["RDKitMolAlign"]
        self.cpp_info.components["RDKitMolCatalog"].libs = ["RDKitMolCatalog"]
        self.cpp_info.components["RDKitMolChemicalFeatures"].libs = ["RDKitMolChemicalFeatures"]
        self.cpp_info.components["RDKitMolDraw2D"].libs = ["RDKitMolDraw2D"]
        self.cpp_info.components["RDKitMolEnumerator"].libs = ["RDKitMolEnumerator"]
        self.cpp_info.components["RDKitMolHash"].libs = ["RDKitMolHash"]
        self.cpp_info.components["RDKitMolInterchange"].libs = ["RDKitMolInterchange"]
        self.cpp_info.components["RDKitMolStandardize"].libs = ["RDKitMolStandardize"]
        self.cpp_info.components["RDKitMolTransforms"].libs = ["RDKitMolTransforms"]
        self.cpp_info.components["RDKitO3AAlign"].libs = ["RDKitO3AAlign"]
        self.cpp_info.components["RDKitOptimizer"].libs = ["RDKitOptimizer"]
        self.cpp_info.components["RDKitPartialCharges"].libs = ["RDKitPartialCharges"]
        self.cpp_info.components["RDKitRascalMCES"].libs = ["RDKitRascalMCES"]
        self.cpp_info.components["RDKitRDGeneral"].libs = ["RDKitRDGeneral"]
        self.cpp_info.components["RDKitRDGeometryLib"].libs = ["RDKitRDGeometryLib"]
        self.cpp_info.components["RDKitRDStreams"].libs = ["RDKitRDStreams"]
        self.cpp_info.components["RDKitReducedGraphs"].libs = ["RDKitReducedGraphs"]
        self.cpp_info.components["RDKitRGroupDecomposition"].libs = ["RDKitRGroupDecomposition"]
        self.cpp_info.components["RDKitRingDecomposerLib"].libs = ["RDKitRingDecomposerLib"]
        self.cpp_info.components["RDKitScaffoldNetwork"].libs = ["RDKitScaffoldNetwork"]
        self.cpp_info.components["RDKitShapeHelpers"].libs = ["RDKitShapeHelpers"]
        self.cpp_info.components["RDKitSimDivPickers"].libs = ["RDKitSimDivPickers"]
        self.cpp_info.components["RDKitSLNParse"].libs = ["RDKitSLNParse"]
        self.cpp_info.components["RDKitSmilesParse"].libs = ["RDKitSmilesParse"]
        self.cpp_info.components["RDKitSubgraphs"].libs = ["RDKitSubgraphs"]
        self.cpp_info.components["RDKitSubstructLibrary"].libs = ["RDKitSubstructLibrary"]
        self.cpp_info.components["RDKitSubstructMatch"].libs = ["RDKitSubstructMatch"]
        self.cpp_info.components["RDKitTautomerQuery"].libs = ["RDKitTautomerQuery"]
        self.cpp_info.components["RDKitTrajectory"].libs = ["RDKitTrajectory"]






