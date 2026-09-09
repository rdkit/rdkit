"""
Test to ensure all generated rdkit-stubs .pyi files are syntactically valid Python.
This prevents regressions where invalid stubs are shipped, breaking mypy for users.
"""
import ast
import pathlib
import sys
import unittest

# 1. Capture the argument BEFORE unittest touches sys.argv
STUBS_DIR = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path("/tmp/stubs_test")

# 2. Truncate sys.argv so unittest doesn't try to interpret the path as a test name
sys.argv = [sys.argv[0]]

class TestStubsSyntax(unittest.TestCase):
    def test_all_pyi_files_parse(self):
        self.assertTrue(STUBS_DIR.exists(), f"Stubs directory not found: {STUBS_DIR}")
        
        pyi_files = list(STUBS_DIR.rglob("*.pyi"))
        self.assertGreater(len(pyi_files), 0, "No .pyi files found to test!")
        
        errors = []
        for pyi_file in pyi_files:
            try:
                ast.parse(pyi_file.read_text(encoding="utf-8"))
            except SyntaxError as e:
                errors.append(f"{pyi_file.relative_to(STUBS_DIR)}:{e.lineno}: {e.msg}")
                
        if errors:
            self.fail("Found syntax errors in generated stubs:\n" + "\n".join(errors))

if __name__ == "__main__":
    unittest.main()