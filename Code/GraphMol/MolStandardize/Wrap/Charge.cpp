#include <RDBoost/Wrap.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/Charge.h>

namespace python = boost::python;
using namespace RDKit;

namespace {

std::vector<MolStandardize::ChargeCorrection> defaultChargeCorrections() {
	return MolStandardize::CHARGE_CORRECTIONS;
}

ROMol *reionizeHelper(MolStandardize::Reionizer &self, const ROMol &mol) {
	return self.reionize(mol);
}

} // namespace

BOOST_PYTHON_MODULE(Charge) {
	python::scope().attr("__doc__") = 
					"Module containing functions for charge corrections";
	
	std::string docString = "";

	python::class_<MolStandardize::ChargeCorrection, boost::noncopyable>(
				"ChargeCorrection", python::init<std::string, std::string, int>())
					.def_readwrite("Name", &MolStandardize::ChargeCorrection::Name)
					.def_readwrite("Smarts", &MolStandardize::ChargeCorrection::Smarts)
					.def_readwrite("Charge", &MolStandardize::ChargeCorrection::Charge)
					;

	python::def("CHARGE_CORRECTIONS", defaultChargeCorrections);

	python::class_<MolStandardize::Reionizer, boost::noncopyable>(
				"Reionizer", python::init<>())
					.def("reionize", &reionizeHelper,
							 (python::arg("self"), python::arg("mol")),
							 "",
							 python::return_value_policy<python::manage_new_object>())					;


}

