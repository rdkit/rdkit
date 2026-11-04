# Converts Data/DefaultMonomers.json into a C++ source file so that the
# built-in MacroMol templates are available without access to the data file.
#
# Usage: cmake -DINPUT=<json file> -DOUTPUT=<cpp file> -P EmbedDefaultMonomers.cmake
#
# The data is written as a sequence of adjacent raw string literals, split at
# line boundaries, because MSVC limits each string literal to 16384 bytes
# (before MSVC 14.51). The file must be ASCII, because MSVC may not read the
# generated source as UTF-8; use \uXXXX escapes for other characters.

set(delimiter "json")
set(maxChunkLength 16000)

file(READ "${INPUT}" contents)

string(REGEX MATCH "[^ -~\t\r\n]" nonAscii "${contents}")
if(NOT nonAscii STREQUAL "")
  message(FATAL_ERROR "${INPUT} contains non-ASCII characters; "
                      "use \\uXXXX escapes instead")
endif()
string(FIND "${contents}" ")${delimiter}\"" delimiterPos)
if(NOT delimiterPos EQUAL -1)
  message(FATAL_ERROR "${INPUT} contains the raw string delimiter "
                      "')${delimiter}\"'")
endif()

file(WRITE "${OUTPUT}"
"// Generated from ${INPUT} by EmbedDefaultMonomers.cmake. Do not edit.
#include <string_view>

namespace RDKit {
namespace detail {
namespace {
constexpr char defaultMonomerDefinitions[] =
")

string(LENGTH "${contents}" remaining)
while(remaining GREATER 0)
  if(remaining GREATER maxChunkLength)
    # end the chunk after the last newline that fits, if there is one
    string(SUBSTRING "${contents}" 0 ${maxChunkLength} chunk)
    string(FIND "${chunk}" "\n" lastNewline REVERSE)
    if(lastNewline GREATER -1)
      math(EXPR chunkLength "${lastNewline} + 1")
      string(SUBSTRING "${contents}" 0 ${chunkLength} chunk)
    endif()
  else()
    set(chunk "${contents}")
  endif()
  string(LENGTH "${chunk}" chunkLength)
  file(APPEND "${OUTPUT}" "R\"${delimiter}(${chunk})${delimiter}\"\n")
  string(SUBSTRING "${contents}" ${chunkLength} -1 contents)
  math(EXPR remaining "${remaining} - ${chunkLength}")
endwhile()

file(APPEND "${OUTPUT}"
"\"\";
}  // namespace

std::string_view getEmbeddedDefaultMonomerDefinitions() {
  return {defaultMonomerDefinitions, sizeof(defaultMonomerDefinitions) - 1};
}
}  // namespace detail
}  // namespace RDKit
")
