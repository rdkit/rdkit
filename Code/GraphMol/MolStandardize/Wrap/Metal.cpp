#include <RDBoost/Wrap.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/MolStandardize/Metal.h>

namespace python = boost::python;

namespace {
RDKit::ROMol *disconnect(RDKit::MolStandardize::MetalDisconnector &self, 
								RDKit::ROMol &mol) {
	RDKit::ROMol* nm = self.disconnect(mol);
	return nm;
}

std::string getMetalNofHelper(RDKit::MolStandardize::MetalDisconnector &self) {
  return RDKit::MolToSmarts(*(self.getMetalNof()));
}

std::string getMetalNonHelper(RDKit::MolStandardize::MetalDisconnector &self) {
  return RDKit::MolToSmarts(*(self.getMetalNon()));
}

void setMetalNonHelper(RDKit::MolStandardize::MetalDisconnector &self, 
								const RDKit::ROMol &mol) {
  self.setMetalNon(mol);
}

void setMetalNofHelper(RDKit::MolStandardize::MetalDisconnector &self, 
								const RDKit::ROMol &mol) {
  self.setMetalNof(mol);
}

} // namespace

BOOST_PYTHON_MODULE(Metal) {
	python::scope().attr("__doc__") = 
			"Module containing functions for molecular standardization";

	std::string docString = "";

	python::class_<RDKit::MolStandardize::MetalDisconnector, boost::noncopyable>(
				"MetalDisconnector", "a class to disconnect metals that are defined as covalently bonded to non-metals",
				python::init<>())
				.add_property("MetalNof", 
												&getMetalNofHelper,
												"RDKit::ROMol containing the metals to disconnect if attached to Nitrogen, Oxygen or Fluorine")
				.add_property("MetalNon", 
												&getMetalNonHelper,
												"RDKit::ROMol containing the metals to disconnect other inorganic elements")
				.def("SetMetalNon", 
												&setMetalNonHelper,
												(python::arg("self"), python::arg("mol")),
												"")
				.def("SetMetalNof", 
												&setMetalNofHelper,
												(python::arg("self"), python::arg("mol")),
												"")
				.def("Disconnect", 
							&disconnect,
	        (python::arg("self"), python::arg("mol")), docString.c_str(),
					python::return_value_policy<python::manage_new_object>());


//  python::def("Disconnect", RDKit::MolStandardize::MetalDisconnector::disconnect,
//              (python::arg("mol")),
//              docString.c_str());
//  python::def("Disconnect", disconnect,
//              (python::arg("mol")), docString.c_str());


}
