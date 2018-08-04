#include <RDBoost/Wrap.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/Normalize.h>

namespace python = boost::python;
using namespace RDKit;

namespace {

ROMol* normalizeHelper(MolStandardize::Normalizer &self, const ROMol &mol) {
	return self.normalize(mol);
}

} // namespace

struct normalize_wrapper {
	static void wrap() {
	python::scope().attr("__doc__") = 
					"Module containing tools for normalizing molecules defined by SMARTS patterns";

	std::string docString = "";

	python::class_<MolStandardize::Normalizer, boost::noncopyable>(
					"Normalizer", python::init<>())
					.def(python::init<std::string, unsigned int>())
					.def("normalize", &normalizeHelper,
								(python::arg("self"), python::arg("mol")),
								"",
								python::return_value_policy<python::manage_new_object>())
					;
	}
};

void wrap_normalize() { normalize_wrapper::wrap(); }


