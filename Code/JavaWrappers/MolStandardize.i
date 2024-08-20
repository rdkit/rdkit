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

%include <GraphMol/MolStandardize/MolStandardize.h>

#if defined SWIGJAVA
%include "enumtypeunsafe.swg"
#endif
%include <GraphMol/MolStandardize/Pipeline.h>
