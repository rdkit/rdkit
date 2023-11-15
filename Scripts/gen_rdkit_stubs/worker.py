import sys
import os
import pybind11_stubgen
from pybind11_stubgen.parser.mixins.parse import ExtractSignaturesFromPybind11Docstrings
from pybind11_stubgen.printer import Printer
import importlib


"""
worker script
1st param is the temporary directory
2nd param is the module name
"""
if __name__ == "__main__":
    tempdir = sys.argv[1]
    module_name = sys.argv[2]
    outer_dirs = sys.argv[3:]
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    gen_rdkit_stubs = importlib.import_module("gen_rdkit_stubs")

    if hasattr(ExtractSignaturesFromPybind11Docstrings, "parse_function_docstring"):
        parse_function_docstring_orig = ExtractSignaturesFromPybind11Docstrings.parse_function_docstring

        def parse_function_docstring_patched(self, func_name, doc_lines, **kwargs):
            doc_lines = gen_rdkit_stubs.ProcessDocLines.process(module_name, doc_lines)
            return parse_function_docstring_orig(self, func_name, doc_lines, **kwargs)

        ExtractSignaturesFromPybind11Docstrings.parse_function_docstring = parse_function_docstring_patched

    if hasattr(Printer, "print_submodule_import"):
        print_submodule_import_orig = Printer.print_submodule_import

        def print_submodule_import_patched(self, name):
            return [f"from .{name} import *"]

        Printer.print_submodule_import = print_submodule_import_patched

    stored_argv = list(sys.argv)
    try:
        sys.argv = ["",
                    "--root-suffix",
                    "",
                    "--ignore-all-errors",
                    "-o",
                    tempdir,
                    module_name]
        pybind11_stubgen.main()
    except Exception as e:
        if isinstance(e, AssertionError):
            raise
        else:
            print(str(e))
    finally:
        sys.argv = stored_argv
    gen_rdkit_stubs.copy_stubs(os.path.join(tempdir, *module_name.split(".")), outer_dirs)
