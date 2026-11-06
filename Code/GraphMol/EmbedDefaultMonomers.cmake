# Converts Data/DefaultMonomers.json into a C++ source file so that the
# built-in MacroMol templates are available without access to the data file.
#
# Usage: cmake -DINPUT=<json file> -DOUTPUT=<cpp file>
#              [-DSOURCE_ROOT=<dir>] [-DTRACKED_COPY=<file>]
#              -P EmbedDefaultMonomers.cmake
#
# OUTPUT is the file that is compiled. If TRACKED_COPY is given, it is
# overwritten with the same content whenever that content changes; it is
# checked in to git so the generated code can be inspected, following the
# convention used for the flex/bison parsers (e.g. SmilesParse/*.cpp.cmake).
# SOURCE_ROOT, if given, is used to show INPUT as a relative path in the
# generated file, so that the output doesn't depend on the checkout location.
#
# The data is written as a sequence of adjacent raw string literals, split at
# line boundaries, because MSVC limits each string literal to 16384 bytes
# (before MSVC 14.51). The file must be ASCII, because MSVC may not read the
# generated source as UTF-8; use \uXXXX escapes for other characters.

set(delimiter "json")
set(maxChunkLength 16000)

file(READ "${INPUT}" contents)
# carriage returns (e.g. from a Windows checkout) are only JSON whitespace;
# drop them so the output is the same on all platforms
string(REPLACE "\r" "" contents "${contents}")

string(REGEX MATCH "[^ -~\t\n]" nonAscii "${contents}")
if(NOT nonAscii STREQUAL "")
  message(FATAL_ERROR "${INPUT} contains non-ASCII characters; "
                      "use \\uXXXX escapes instead")
endif()
string(FIND "${contents}" ")${delimiter}\"" delimiterPos)
if(NOT delimiterPos EQUAL -1)
  message(FATAL_ERROR "${INPUT} contains the raw string delimiter "
                      "')${delimiter}\"'")
endif()

if(DEFINED SOURCE_ROOT)
  file(RELATIVE_PATH inputName "${SOURCE_ROOT}" "${INPUT}")
else()
  get_filename_component(inputName "${INPUT}" NAME)
endif()

set(generated
"// Generated from ${inputName} by EmbedDefaultMonomers.cmake. Do not edit.
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
  string(APPEND generated "R\"${delimiter}(${chunk})${delimiter}\"\n")
  string(SUBSTRING "${contents}" ${chunkLength} -1 contents)
  math(EXPR remaining "${remaining} - ${chunkLength}")
endwhile()

string(APPEND generated
"\"\";
}  // namespace

std::string_view getEmbeddedDefaultMonomerDefinitions() {
  return {defaultMonomerDefinitions, sizeof(defaultMonomerDefinitions) - 1};
}
}  // namespace detail
}  // namespace RDKit
")

file(WRITE "${OUTPUT}" "${generated}")

if(DEFINED TRACKED_COPY)
  # only rewrite the checked-in copy when its content changes (ignoring line
  # endings), to avoid touching the source tree unnecessarily
  set(trackedContents "")
  if(EXISTS "${TRACKED_COPY}")
    file(READ "${TRACKED_COPY}" trackedContents)
    string(REPLACE "\r" "" trackedContents "${trackedContents}")
  endif()
  if(NOT trackedContents STREQUAL generated)
    file(WRITE "${TRACKED_COPY}" "${generated}")
  endif()
endif()
