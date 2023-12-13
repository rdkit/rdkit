#!/usr/bin/env python

"""
Script to generate Python stubs for RDKit.

This script is invoked as part of the build process
by setting the CMake switch RDK_INSTALL_PYTHON_STUBS=ON.
If you decide to run this script outside the build process,
make sure that the RDKit Python modules for which stubs are
to be generated are the *first* RDKit modules available in
sys.path; otherwise, stubs will not be generated for the
intended RDKit version.

Usage:
./Scripts/gen_rdkit_stubs.py [output_dirs; defaults to $PWD]

Usage example:
$ cd $RDBASE
$ ./Scripts/gen_rdkit_stubs.py
$ cp -R rdkit-stubs $CONDA_PREFIX/lib/python3.*/site-packages

The scripts creates an rdkit-stubs directory in each
directory in output_dirs.
Warnings printed to console can be safely ignored.
"""

import sys
import os
import importlib
import argparse
import multiprocessing
from pathlib import Path
from . import generate_stubs, purge_rdkit_source_dir_from_sys_path, find_rdkit_site_packages_path

def parse_args():
    """Parse command line arguments."""
    default_n_cpus = max(1, multiprocessing.cpu_count() - 2)
    default_output_dirs = [os.getcwd()]
    parser = argparse.ArgumentParser()
    parser.add_argument("--concurrency", type=int,
                        help=f"max number of CPUs to be used (defaults to {default_n_cpus})",
                        default=default_n_cpus)
    parser.add_argument("--verbose",
                        help=f"print non-fatal warnings/errors to stdout (defaults to false)",
                        action="store_true")
    parser.add_argument("--keep-incorrect-staticmethods",
                        help=f"Whether incorrectly assigned staticmethods should be kept as such",
                        action="store_true")
    parser.add_argument("output_dirs", nargs="*",
                        help=f"output directories (defaults to {default_output_dirs[0]})",
                        default=default_output_dirs)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    purge_rdkit_source_dir_from_sys_path()
    site_packages_path = find_rdkit_site_packages_path()
    try:
        importlib.import_module("pybind11_stubgen")
    except ModuleNotFoundError:
        print("Failed to find pybind11_stubgen in PYTHONPATH. "
              "Please pip install pybind11_stubgen (available on PyPI and GitHub).", file=sys.stderr)
        sys.exit(1)
    site_packages_path = Path(site_packages_path)
    generate_stubs(site_packages_path, args)
