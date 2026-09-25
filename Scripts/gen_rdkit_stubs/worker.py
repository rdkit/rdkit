import sys
import os
import argparse
import pybind11_stubgen
from pybind11_stubgen.parser.mixins.parse import ExtractSignaturesFromPybind11Docstrings
from pybind11_stubgen.printer import Printer
import importlib


"""
worker script
1st param is the temporary directory
2nd param is the module name
"""
def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--tempdir",
                        help=f"Temporary directory")
    parser.add_argument("--module-name",
                        help=f"Python module name for which stubs will be generated")
    parser.add_argument("--outer-dirs", nargs="*",
                        help=f"Outer dirs where stubs should be copied")
    parser.add_argument("--keep-incorrect-staticmethods",
                        help=f"Whether incorrectly assigned staticmethods should be kept as such",
                        action="store_true")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    keep_incorrect_staticmethods = args.keep_incorrect_staticmethods or os.environ.get("KEEP_INCORRECT_STATICMETHODS", None)
    if isinstance(keep_incorrect_staticmethods, str):
        try:
            keep_incorrect_staticmethods = bool(int(keep_incorrect_staticmethods))
        except ValueError:
            keep_incorrect_staticmethods = False
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    gen_rdkit_stubs = importlib.import_module("gen_rdkit_stubs")

    if hasattr(ExtractSignaturesFromPybind11Docstrings, "parse_function_docstring"):
        parse_function_docstring_orig = ExtractSignaturesFromPybind11Docstrings.parse_function_docstring

        def parse_function_docstring_patched(self, func_name, doc_lines, **kwargs):
            doc_lines = gen_rdkit_stubs.ProcessDocLines.process(args.module_name, doc_lines, keep_incorrect_staticmethods)
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
                    args.tempdir,
                    args.module_name]
        pybind11_stubgen.main()
    except Exception as e:
        if isinstance(e, AssertionError) or isinstance(e, ImportError):
            raise
        else:
            print(str(e))
    finally:
        sys.argv = stored_argv
    src_path = os.path.join(args.tempdir, *args.module_name.split("."))
    gen_rdkit_stubs.patch_stubs(args.tempdir, src_path)
    gen_rdkit_stubs.copy_stubs(src_path, args.outer_dirs)
