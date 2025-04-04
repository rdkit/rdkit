#!/usr/bin/env python

"""
clang_worker script.
Each worker is a daemon listening for input on stdin
and generating output on stdout.
"""

import sys
import os
import importlib

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
patch_rdkit_docstrings = importlib.import_module("patch_rdkit_docstrings")
CppFile = patch_rdkit_docstrings.CppFile
ClangWorkerData = patch_rdkit_docstrings.ClangWorkerData

class MainLoop:
    """Main worker class.

    The loop breaks once there is no more input to read from stdin.
    """
    def __init__(self, clang_worker_data):
        self.clang_worker_data = clang_worker_data

    def run(self):
        while 1:
            line = sys.stdin.readline().strip()
            if not line:
                break
            cpp_class_file_json = line.strip()
            cpp_class_file = CppFile.from_json(cpp_class_file_json)
            res = cpp_class_file.generate_ast(self.clang_worker_data.clang_flags)
            if res:
                cpp_class_file.parse_ast(self.clang_worker_data.arg1_func_byclass_dict)
            sys.stdout.write(cpp_class_file.to_json() + "\n")
            sys.stdout.flush()
        sys.stdout.write("\n")
        sys.stdout.flush()


if __name__ == "__main__":
    clang_worker_data = ClangWorkerData.from_json(sys.argv[1])
    main_loop = MainLoop(clang_worker_data)
    main_loop.run()
