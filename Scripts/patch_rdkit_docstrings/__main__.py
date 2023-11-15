#!/usr/bin/env python

"""
patch_rdkit_docstrings is a runnable Python module.
It will scan RDKit sources searching for docstrings
that lack parameter definitions, or member functions that
do not have an explicit "self" parameter, and will patch
the C++ sources accordingly.

Example usage:

$ cd $RDBASE
$ CLANG_INCLUDE_PATH=/build/llvm-project/libcxx/include \
CLANG_PYTHON_BINDINGS_PATH=/build/llvm-project/clang/bindings/python \
QT_INCLUDE_DIRS=$CONDA_PREFIX/include/qt:$CONDA_PREFIX/include/qt/QtGui \
EIGEN3_INCLUDE_DIR=/usr/prog/Eigen/3.3.9-GCCcore-11.2.0/include \
python -m Scripts.patch_rdkit_docstrings
"""

import sys
import argparse
from . import FixSignatures
from ..gen_rdkit_stubs import purge_rdkit_source_dir_from_sys_path, find_rdkit_include_path, RDKIT_MODULE_NAME


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpp-source-path",
                        help="path where RDKit C++ sources are located (defaults to $RDBASE, or to cwd if RDBASE is not set)",
                        default=FixSignatures.cpp_source_path)
    parser.add_argument("--stubs-path",
                        help=f"path to the {RDKIT_MODULE_NAME}-stubs directory (defaults to ./{RDKIT_MODULE_NAME})",
                        default=FixSignatures.rdkit_stubs_path)
    parser.add_argument("--concurrency",
                        help=f"max number of CPUs to be used (defaults to {FixSignatures.concurrency})",
                        default=FixSignatures.concurrency)
    parser.add_argument("--include-path",
                        help=f"main clang include path (defaults to {FixSignatures.include_path})",
                        default=FixSignatures.include_path)
    parser.add_argument("--python-include-path",
                        help=f"Python clang include path (defaults to {FixSignatures.python_include_path})",
                        default=FixSignatures.python_include_path)
    parser.add_argument("--rdkit-include-path",
                        help=f"RDKit include path (defaults to {FixSignatures.rdkit_include_path})",
                        default=find_rdkit_include_path())
    parser.add_argument("--clang-flags",
                        help=f"flags to be passed to clang (defaults to {FixSignatures.clang_flags})",
                        default=FixSignatures.clang_flags)
    parser.add_argument("--user-clang-flags",
                        help=f"user-defined flags to be passed to clang (defaults to {FixSignatures.user_clang_flags})",
                        default=FixSignatures.clang_flags)
    parser.add_argument("--log-level",
                        help=f"logging level (defaults to {FixSignatures.log_level})",
                        default=FixSignatures.log_level)
    parser.add_argument("--clean", action="store_true", help="force removing all RDKDOCORIG files")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    purge_rdkit_source_dir_from_sys_path()
    if args.rdkit_include_path is None:
        print(f"Failed to find RDKit include path. Please set it through the --rdkit-include-path switch.", file=sys.stderr)
        sys.exit(1)
    fix_signatures = FixSignatures(args)
