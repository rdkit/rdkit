%{
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/Pipeline.h>

namespace RDKit {
namespace MolStandardize {
    bool operator==(const PipelineLogEntry & rhs, const PipelineLogEntry & lhs) {
        return (rhs.status == lhs.status) && (rhs.detail == lhs.detail);
    }
    bool operator!=(const PipelineLogEntry & rhs, const PipelineLogEntry & lhs) {
        return !(rhs == lhs);
    }
}
}
%}

%include <std_vector.i>
namespace std {
%template(PipelineLog) std::vector<RDKit::MolStandardize::PipelineLogEntry>;
}

%newobject RDKit::MolStandardize::cleanup;
%newobject RDKit::MolStandardize::normalize;
%newobject RDKit::MolStandardize::reionize;
%newobject RDKit::MolStandardize::removeFragments;
%newobject RDKit::MolStandardize::canonicalTautomer;
%newobject RDKit::MolStandardize::tautomerParent;
%newobject RDKit::MolStandardize::fragmentParent;
%newobject RDKit::MolStandardize::stereoParent;
%newobject RDKit::MolStandardize::isotopeParent;
%newobject RDKit::MolStandardize::chargeParent;
%newobject RDKit::MolStandardize::superParent;
%newobject RDKit::MolStandardize::disconnectOrganometallics;

%include <GraphMol/MolStandardize/MolStandardize.h>

%include "enums.swg"
%include <GraphMol/MolStandardize/Pipeline.h>
