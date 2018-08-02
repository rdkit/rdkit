#include <RDBoost/Wrap.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/Validate.h>

namespace python = boost::python;
using namespace RDKit;

namespace {

python::list rdkitValidate(MolStandardize::RDKitValidation &self, 
								const ROMol &mol, const bool reportAllFailures) {
	python::list res;
	std::vector<MolStandardize::ValidationErrorInfo> errout = 
					self.validate(mol, reportAllFailures);
  for (auto &query : errout) {
    std::string msg = query.message();
		res.append(msg);
  }
	return res;
}

//MolStandardize::MolVSValidation *getMolVSValidation(python::object validations) {
//	std::vector<MolStandardize::MolVSValidations*> vs;
//
//	std::unique_ptr<std::vector<MolStandardize::MolVSValidations*>> pvect = 
//					pythonObjectToVect<MolStandardize::MolVSValidations*>(validations);
//
//	for (auto v : *pvect) {
//		vs.push_back(new MolStandardize::MolVSValidations(*v));
//	}
//	return new MolStandardize::MolVSValidation(vs);
//}

} // namespace

BOOST_PYTHON_MODULE(Validate) {
	python::scope().attr("__doc__") = 
					"Module containing tools to perform all validations";

	std::string docString = "";

	python::class_<MolStandardize::RDKitValidation, boost::noncopyable>(
									"RDKitValidation", python::init<>() )
					.def("validate", rdkitValidate,
							(python::arg("self"), python::arg("mol"), python::arg("reportAllFailures") = false),
							"")
					;

	python::class_<MolStandardize::MolVSValidations, boost::noncopyable>(
					"MolVSValidations", python::init<>() )
					.def("run", &MolStandardize::MolVSValidations::run, 
							(python::arg("self"), python::arg("mol"), python::arg("reportAllFailures"),
								python::arg("errors")),
							 "")
					;
//	python::class_<MolStandardize::NoAtomValidation, 
//					python::bases<MolStandardize::MolVSValidations>>(
//													"NoAtomValidation", python::init<>() )
//					.def("run", &MolStandardize::NoAtomValidation::run,
//							(python::arg("self"), python::arg("mol"), python::arg("reportAllFailures"),
//								python::arg("errors")),
//							 "")
//					;
					




	python::class_<MolStandardize::MolVSValidation, boost::noncopyable>(
									"MolVSValidation", python::init<>() )
//					.def(python::init<const std::vector<MolStandardize::MolVSValidations*>())
//					.def("__init__", python::make_constructor(&getMolVSValidation))
					.def("Validate", &MolStandardize::MolVSValidation::Validate,
							(python::arg("self"), python::arg("mol"), python::arg("reportAllFailures") = false),
							"")
					;
}
